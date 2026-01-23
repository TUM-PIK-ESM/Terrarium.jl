"""
    GroundEvaporation{NF, GR} <: AbstractEvapotranspiration

Evaporation scheme for bare ground that calculates the humidity flux as

```math
E = \\beta \\frac{\\Delta q}{r_a}
```
where `Δq` is the vapor pressure deficit in terms of specific humidity, `rₐ` is aerodynamic resistance,
and `β` is an evaporation limiting factor.
"""
struct GroundEvaporation{NF, GR<:AbstractGroundEvaporationResistanceFactor} <: AbstractEvapotranspiration
    "Parameterization for ground resistance to evaporation/sublimation"
    ground_resistance::GR
end

GroundEvaporation(
    ::Type{NF};
    ground_resistance::GR = ConstantEvaporationResistanceFactor(one(NF))
) where {NF, GR} = GroundEvaporation{NF, GR}(; ground_resistance)

@propagate_inbounds surface_humidity_flux(i, j, grid, state, evtr::GroundEvaporation) = state.evaporation_ground[i, j]

# Process methods

variables(::GroundEvaporation) = (
    auxiliary(:evaporation_ground, XY(), units=u"m/s", desc="Ground evaporation contribution to surface humidity flux"),
)

function compute_auxiliary!(state, model, evap::GroundEvaporation)
    grid = get_grid(model)
    atmos = get_atmosphere(model)
    seb = get_surface_energy_balance(model)
    constants = get_constants(model)
    launch!(grid, state, :xy, compute_evaporation_kernel!, evap, seb.skin_temperature, atmos, constants)
end

# Kernels

@kernel function compute_evaporation_kernel!(
    state, grid,
    evap::GroundEvaporation,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    i, j = @index(Global, NTuple)
    Ts = skin_temperature(i, j, grid, state, skinT)
    rₐ = aerodynamic_resistance(i, j, grid, state, aerodynamic_resistance) # aerodynamic resistance
    β = ground_evaporation_resistance_factor(i, j, grid, state, evap.ground_resistance)
    Δq = compute_humidity_vpd(i, j, grid, state, atmos, constants, Ts)
    # Calculate water evaporation flux in m/s (positive upwards)
    state.evaporation_ground[i, j, 1] = β * Δq / rₐ
end

@kernel function compute_evaporation_kernel!(
    state, grid,
    evap::GroundEvaporation,
    soilw::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    i, j = @index(Global, NTuple)
    Ts = skin_temperature(i, j, grid, state, skinT)
    rₐ = aerodynamic_resistance(i, j, grid, state, aerodynamic_resistance) # aerodynamic resistance
    β = ground_evaporation_resistance_factor(i, j, grid, state, evap.ground_resistance, soilw, strat, bgc)
    Δq = compute_humidity_vpd(i, j, grid, state, atmos, constants, Ts)
    # Calculate water evaporation flux in m/s (positive upwards)
    state.evaporation_ground[i, j, 1] = β * Δq / rₐ
end

# Ground resistance to evaporation

@kwdef struct ConstantEvaporationResistanceFactor{NF} <: AbstractGroundEvaporationResistanceFactor
    "Unit interval factor that determines resistance to evaporation; zero corresponds to no evaporation"
    factor::NF = 1.0
end

ConstantEvaporationResistanceFactor(::Type{NF}; kwargs...) where {NF} = ConstantEvaporationResistanceFactor{NF}(; kwargs...)

@inline ground_evaporation_resistance_factor(i, j, grid, state, res::ConstantEvaporationResistanceFactor, args...) = res.factor

"""
    SoilMoistureResistanceFactor <: AbstractGroundEvaporationResistanceFactor

Implements the soil moisture limiting resistance factor of Lee and Pielke (1992),

```math
\\beta = \\begin{cases}
\\frac{1}{4} \\left(1 - \\cos\\left(π \\frac{\\theta_1}{\\theta_{\\text{fc}}} \\right)\\right) & \\theta_1 < \\theta_{\\text{fc}} \\
1 & \\text{otherwise}
\\end{cases}
```
"""
struct SoilMoistureResistanceFactor{NF} <: AbstractGroundEvaporationResistanceFactor end

SoilMoistureResistanceFactor(::Type{NF}) where {NF} = SoilMoistureResistanceFactor{NF}()

# Fallback implementation for interface consistency
ground_evaporation_resistance_factor(i, j, grid, state, ::SoilMoistureResistanceFactor{NF}, args...) where {NF} = one(NF)

@inline function ground_evaporation_resistance_factor(
    i, j, grid, state,
    ::SoilMoistureResistanceFactor{NF},
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
) where {NF}
    fgrid = get_field_grid(grid)
    soil = soil_volume(i, j, fgrid.Nz, grid, state, strat, hydrology, bgc)
    fc = field_capacity(get_hydraulic_properties(hydrology), soil)
    fracs = volumetric_fractions(soil)
    if fracs.water < fc
        β = (1 - cos(π * fracs.water / fc))^2 / 4
    else
        β = NF(1)
    end
    return β
end

# Forcing interface for soil hydrology

"""
    forcing(i, j, k, grid, state, evtr::AbstractEvapotranspiration, ::AbstractSoilHydrology)

Compute and return the evapotranspiration forcing for soil moisture at the given indices `i, j, k`.
The ET forcing is just the `surface_humidity_flux` rescaled by the thickness of layer `k`.
"""
@inline function forcing(i, j, k, grid, state, evtr::AbstractEvapotranspiration, ::AbstractSoilHydrology)
    let Δz = Δzᵃᵃᶜ(i, j, k, grid),
        Qh = surface_humidity_flux(i, j, grid, state, evtr);
        ∂θ∂t = -Qh / Δz # rescale by layer thickness to get water content flux
        return ∂θ∂t * (k == grid.Nz)
    end
end
