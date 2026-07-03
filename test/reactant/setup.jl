# Model registry for the CPU-vs-Reactant correctness suite.
#
# Add a new tested configuration by adding a `build_model(::Val{:name}, arch, NF)` method that
# returns a NamedTuple `(; model, boundary_conditions, initializers, Δt)`. `arch` is the ONLY
# thing that differs between the CPU and Reactant runs — everything else is identical.

import RingGrids

# --- generic helpers --------------------------------------------------------------------

function build_integrator(v::Val, arch, NF)
    cfg = build_model(v, arch, NF)
    return Terrarium.initialize(cfg.model;
        boundary_conditions = cfg.boundary_conditions,
        initializers = cfg.initializers)
end

cpu_dt(v::Val, NF) = build_model(v, CPU(), NF).Δt

# --- :soil_heat_column — minimal single-column SoilModel (heat conduction) --------------

function build_model(::Val{:soil_heat_column}, arch, NF)
    grid = ColumnGrid(arch, NF, ExponentialSpacing(N = 10))
    model = SoilModel(grid)
    # constant surface temperature; default (zero-flux) bottom boundary
    bcs = PrescribedSurfaceTemperature(:T_ub, NF(1))
    # linear initial temperature profile with depth (z ≤ 0)
    inits = (temperature = (x, z) -> NF(-1) - NF(0.05) * z,)
    return (; model, boundary_conditions = bcs, initializers = inits, Δt = NF(600))
end

# --- :soil_heat_global — global ColumnRingGrid SoilModel (heat conduction) --------------
# Mirrors examples/simulations/soil_heat_global.jl, but with a synthetic all-land mask
# (no NetCDF inputs in CI) and a purely functional, smooth surface BC (no array gather).

function build_model(::Val{:soil_heat_global}, arch, NF)
    rings = RingGrids.FullGaussianGrid(4)                 # small: 4 nlat_half
    grid = ColumnRingGrid(arch, NF, ExponentialSpacing(N = 20), rings)  # default all-active mask
    model = SoilModel(grid)

    seconds_per_day = NF(24 * 3600)
    # smooth surface forcing as a function of column index x and time t (no indexed lookup)
    surface_temperature(x, t) = NF(5) * sin(NF(2π) * t / seconds_per_day) * cos(NF(x))
    bcs = PrescribedSurfaceTemperature(:T_ub, surface_temperature)

    # smooth, deterministic initial temperature profile
    inits = (temperature = (x, z) -> NF(2) * cos(NF(x)) - NF(0.05) * z,)
    return (; model, boundary_conditions = bcs, initializers = inits, Δt = NF(600))
end
