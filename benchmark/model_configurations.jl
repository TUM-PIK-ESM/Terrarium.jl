# Model configuration registry for the Terrarium benchmark suite.
#
# A configuration is one `build_model(::Val{:name}, arch, NF; nlat_half, nz, model_kwargs)` method
# returning `(; model, boundary_conditions, initializers, Δt)`. The architecture and the resolution
# (`nlat_half` horizontally, `nz` vertically) are the ONLY things that vary between benchmark runs;
# everything else is fixed by the configuration so that numbers from different architectures and
# resolutions stay comparable.
#
# This mirrors the registry pattern used by `test/reactant/setup.jl`.
#
# Reactant note: closures that end up inside a traced time step (boundary conditions in particular)
# must capture only `isbits` values — a closure over `NF` captures a `Type` and fails to compile.
# Numeric constants are therefore hoisted into local variables before the closure is formed.

import RingGrids

"""
All configuration names known to the benchmark suite, in the order they appear in tables.
"""
const CONFIGURATIONS = (:land, :land_no_vegetation, :soil_heat)

"""
Horizontal grid for a given `nlat_half`: a full Gaussian grid with `2 * nlat_half` latitude rings.
"""
benchmark_rings(nlat_half::Integer) = RingGrids.FullGaussianGrid(nlat_half)

"""
Number of columns (horizontal grid points) of the benchmark grid at `nlat_half`. The benchmark grids
use an all-land mask, so every grid point of the ring grid is an active column: `8 * nlat_half^2`.
"""
ncolumns(nlat_half::Integer) = RingGrids.get_npoints(benchmark_rings(nlat_half))

"""
Approximate horizontal resolution in degrees latitude at the given `nlat_half`.
"""
resolution_degrees(nlat_half::Integer) = 180 / (2 * nlat_half)

"""
Benchmark grid: `nlat_half` sets the horizontal resolution, `nz` the number of soil layers.

The vertical spacing grows quasi-exponentially from 5 cm at the surface to 1 m at the bottom, as in
`test/coupled_models/land_model_tests.jl`, so that `nz` changes the resolution of a soil column of
roughly fixed depth rather than the depth itself.
"""
function benchmark_grid(arch, ::Type{NF}; nlat_half::Integer, nz::Integer) where {NF}
    spacing = ExponentialSpacing(Δz_min = NF(0.05), Δz_max = NF(1), N = nz)
    return ColumnRingGrid(arch, NF, spacing, benchmark_rings(nlat_half))
end

# Resolve a NamedTuple of per-run model kwargs. Convention (as in SpeedyWeather's benchmark suite):
# a `Type` or `Function` value is called with `NF` — every Terrarium component has an `T(NF)`
# constructor — so a suite can write `model_kwargs = (timestepper = Heun,)`. Everything else is
# passed through unchanged.
_resolve_kwarg_value(v::Type, ::Type{NF}) where {NF} = v(NF)
_resolve_kwarg_value(v::Function, ::Type{NF}) where {NF} = v(NF)
_resolve_kwarg_value(v, ::Type) = v
function resolve_model_kwargs(nt::NamedTuple, ::Type{NF}) where {NF}
    return NamedTuple{keys(nt)}(map(v -> _resolve_kwarg_value(v, NF), values(nt)))
end

"""
Soil texture used by the coupled land configurations: a loam with 20% clay.

This has to be prescribed. The default `SoilTexture` is pure sand (`clay = 0`), for which the SURFEX
field capacity `0.089·(100·clay)^0.35` and wilting point `0.03713·√(100·clay)` both collapse to zero,
so the plant available water `(θ − θ_wp)/(θ_fc − θ_wp)` — and with it the soil moisture limiting
factor β — is not finite, and the Medlyn stomatal conductance asserts during initialization.
"""
benchmark_texture(::Type{NF}) where {NF} = SoilTexture(NF; sand = NF(0.4), clay = NF(0.2), silt = NF(0.4))

"""
Soil process shared by the coupled land configurations: energy, water (Richards equation) and
carbon on a homogeneous loam stratigraphy.
"""
function benchmark_soil(::Type{NF}) where {NF}
    strat = HomogeneousSoilStratigraphy(NF; texture = benchmark_texture(NF))
    hydrology = SoilHydrology(NF, RichardsEq())
    return SoilEnergyWaterCarbon(NF; strat, hydrology)
end

const SECONDS_PER_YEAR = 365.25 * 24 * 3600

"""
Vegetation process shared by the coupled land configurations.

The PALADYN turnover and disturbance rates (`γL`, `γR`, `γS`, `γv_min`) are documented as yr⁻¹ but
are used unconverted in tendencies that are integrated in seconds, which gives the vegetation carbon
pool an e-folding time of ~14 seconds instead of ~14 years. With an explicit time stepper that makes
`carbon_vegetation` oscillate with growing amplitude for any Δt ≳ 25 s until the model aborts; the
unit tests never see it because they take a single time step. Until this is fixed upstream the
benchmark rescales the four rates to s⁻¹, so that the headline configuration measures a converging
simulation. Only parameter *values* change: the code path, and therefore the cost per step, is the
same, and removing this once the rates are converted upstream will not change the timings.
"""
function benchmark_vegetation(::Type{NF}) where {NF}
    per_second(rate) = NF(rate / SECONDS_PER_YEAR)
    carbon_dynamics = PALADYNCarbonDynamics(
        NF;
        γL = per_second(0.3), γR = per_second(0.3), γS = per_second(0.05),
    )
    vegetation_dynamics = PALADYNVegetationDynamics(NF; γv_min = per_second(0.002))
    return VegetationCarbon(NF; carbon_dynamics, vegetation_dynamics)
end

"""
Initial state shared by the coupled land configurations: a warm, variably saturated soil column
under a thin cold snowpack. `carbon_vegetation` is only meaningful with vegetation enabled.
"""
function land_initializers(::Type{NF}; vegetation::Bool) where {NF}
    ## hoisted constants — see the Reactant note at the top of this file
    T_surface = NF(5)
    lapse_rate = NF(0.02)
    saturation_surface = NF(0.8)
    saturation_gradient = NF(0.05)
    saturated = one(NF)
    common = (
        temperature = (x, z) -> T_surface - lapse_rate * z,
        saturation_water_ice = (x, z) -> min(saturated, saturation_surface - saturation_gradient * z),
        snow_water_equivalent = NF(0.2),
        snow_temperature = NF(-2),
    )
    return vegetation ? merge(common, (carbon_vegetation = NF(0.1),)) : common
end

# --- :land — fully coupled land model ---------------------------------------------------
# Vegetation carbon + soil energy/water/carbon with Richards equation + single-layer snow, on top of
# the default surface energy balance, surface hydrology and prescribed atmosphere. This is the
# headline configuration: the cross-architecture overview table is its resolution sweep.
#
# Δt = 600 s: the coupled model raises a DomainError in the turbulent-flux thermodynamics
# (log of a negative saturation) for time steps of roughly half an hour and larger.

function build_model(::Val{:land}, arch, ::Type{NF}; nlat_half::Integer, nz::Integer, model_kwargs = (;)) where {NF}
    grid = benchmark_grid(arch, NF; nlat_half, nz)
    soil = benchmark_soil(NF)
    snow = SingleLayerSnow(NF)
    vegetation = benchmark_vegetation(NF)
    model = LandModel(grid; soil, snow, vegetation, resolve_model_kwargs(model_kwargs, NF)...)
    initializers = land_initializers(NF; vegetation = true)
    return (; model, boundary_conditions = (;), initializers, Δt = NF(600))
end

# --- :land_no_vegetation — same, without the canopy ---------------------------------------
# `vegetation = nothing` also selects `BareGroundEvaporation` and `NoCanopyInterception`, so the
# pair (:land, :land_no_vegetation) isolates the cost of the vegetation and canopy path.

function build_model(::Val{:land_no_vegetation}, arch, ::Type{NF}; nlat_half::Integer, nz::Integer, model_kwargs = (;)) where {NF}
    grid = benchmark_grid(arch, NF; nlat_half, nz)
    soil = benchmark_soil(NF)
    snow = SingleLayerSnow(NF)
    model = LandModel(grid; soil, snow, vegetation = nothing, resolve_model_kwargs(model_kwargs, NF)...)
    initializers = land_initializers(NF; vegetation = false)
    return (; model, boundary_conditions = (;), initializers, Δt = NF(600))
end

# --- :soil_heat — minimal soil heat conduction (Reactant reference) ------------------------
# Mirrors `test/reactant/setup.jl`'s `:soil_heat_global`: the only configuration currently proven to
# compile under Reactant. It is not part of the headline comparison — it exists so that an empty
# Reactant column can be attributed to the land physics rather than to the benchmark harness, which
# is also why it keeps the plain `SoilModel(grid)` defaults (including the default sandy texture)
# rather than the loam of the land configurations.

function build_model(::Val{:soil_heat}, arch, ::Type{NF}; nlat_half::Integer, nz::Integer, model_kwargs = (;)) where {NF}
    grid = benchmark_grid(arch, NF; nlat_half, nz)
    model = SoilModel(grid; resolve_model_kwargs(model_kwargs, NF)...)

    ## smooth surface forcing as a function of column index x and time t; the closure must capture
    ## only isbits values to be a valid kernel argument under Reactant
    amplitude = NF(5)
    ω = NF(2π) / NF(24 * 3600)
    surface_temperature(x, t) = amplitude * sin(ω * t) * cos(x)
    boundary_conditions = PrescribedSurfaceTemperature(:T_ub, surface_temperature)

    T₀ = NF(2)
    lapse_rate = NF(0.05)
    initializers = (temperature = (x, z) -> T₀ * cos(NF(x)) - lapse_rate * z,)
    return (; model, boundary_conditions, initializers, Δt = NF(600))
end

"""
    build_model(name::Symbol, arch, NF; nlat_half, nz, model_kwargs)

Convenience wrapper dispatching on the configuration `name`.
"""
function build_model(name::Symbol, arch, ::Type{NF}; kwargs...) where {NF}
    name in CONFIGURATIONS || error("Unknown configuration `$name`. Known: $(join(CONFIGURATIONS, ", ")).")
    return build_model(Val(name), arch, NF; kwargs...)
end
