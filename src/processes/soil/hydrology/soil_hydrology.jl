"""
Base type for implementations of soil water flow dynamics.
"""
abstract type AbstractSoilWaterFlowOperator <: AbstractOperator end

"""
Represents the simplest case of immobile soil water.
"""
struct NoFlow <: AbstractSoilWaterFlowOperator end

"""
    RichardsEq{NF} <: AbstractSoilWaterFlowOperator

Operator for soil hydrology corresponding to the Richardson-Richards equation for variably saturated
flow in porous media.
"""
@kwdef struct RichardsEq{NF} <: AbstractSoilWaterFlowOperator
    "Closure relation for mapping between water potential (hydraulic head) and saturation"
    saturation_closure = PressureSaturationRelation()
end

get_closure(op::RichardsEq) = op.saturation_closure

"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
struct SoilHydrology{
    NF,
    Operator<:AbstractSoilWaterFlowOperator,
    SoilHydraulicProperties<:AbstractSoilHydraulics{NF}
} <: AbstractSoilHydrology{NF}
    "Soil water flow scheme"
    operator::Operator

    "Soil hydraulic properties parameterization"
    hydraulic_properties::SoilHydraulicProperties
end

SoilHydrology(
    ::Type{NF};
    operator::AbstractSoilWaterFlowOperator = NoFlow(),
    hydraulic_properties::AbstractSoilHydraulics{NF} = SURFEXHydraulics(NF),
) where {NF} = SoilHydrology(operator, hydraulic_properties)

"""
    default_swrc(::AbstractSoilWaterFlowOperator, ::AbstractSoilHydraulics)

Return the default soil water retention curve (SWRC) for the given soil hydrology configuration.
Defaults to `nothing` which represents no use of a pressure-saturation relation.
"""
default_swrc(::AbstractSoilWaterFlowOperator, ::AbstractSoilHydraulics) = nothing

get_hydraulic_properties(hydrology::SoilHydrology) = hydrology.hydraulic_properties

# TODO: This method interface assumes a single water retenction curve for the whole stratigraphy;
# we should ideally relax this assumption for multi-layer stratigraphies
get_soil_water_retention_curve(hydrology::SoilHydrology) = hydrology.swrc

"""
    porosity(idx, state, hydrology::SoilHydrology, strat::AbstractStratigraphy, bgc::AbstractSoilBiogeochemistry)

Return the porosity of the soil volume at `idx` given the current state, hydrology, stratigraphy, and biogeochemistry configurations.
"""
@inline function porosity(idx, state, hydrology::SoilHydrology, strat::AbstractStratigraphy, bgc::AbstractSoilBiogeochemistry)
    props = get_hydraulic_properties(hydrology)
    org = organic_fraction(idx, state, bgc)
    texture = soil_texture(idx, state, strat)
    return (1 - org)*mineral_porosity(props, texture) + org*organic_porosity(idx, state, bgc)
end

# Immobile soil water (NoFlow)

variables(::SoilHydrology{NF,NoFlow}) where {NF} = (
    auxiliary(:saturation_water_ice, XYZ(), domain=UnitInterval(), desc="Saturation level of water and ice in the pore space"),
)

@inline compute_auxiliary!(state, model, hydrology::SoilHydrology) = nothing

@inline compute_tendencies!(state, model, strat::SoilHydrology{NF,NoFlow}) where {NF} = nothing

# TODO: Richardson-Richards equation diffusion/advection

variables(hydrology::SoilHydrology{NF,<:RichardsEq}) where {NF} = (
    prognosic(:matric_potential, XYZ(), get_closure(hydrology.operator)),
    # decalre hydraulic_conductivity as auxiliary variable on grid cell faces
    auxiliary(:hydraulic_conductivity, XYZ(z=Face()), domain=UnitInterval(), units=u"m/s", desc="Hydraulic conductivity at grid cell faces [m/s]"),
)

function compute_auxiliary!(state, model, hydrology::SoilHydrology{NF, <:RichardsEq}) where {NF}
    grid = get_grid(model)
    hydrology = get_soil_hydrology(model)
    strat = get_stratigraphy(model)
    launch!(grid, :xyz, compute_hydraulic_conductivity!, state, grid, hydrology, energy, strat, bgc)
    return nothing
end

@kernel function compute_hydraulic_conductivity!(
    state,
    grid,
    hydrology::SoilHydrology,
    energy::AbstractSoilEnergyBalance,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
)
    idx = @index(Global, NTuple)
    state.hydraulic_conductivity[idx...] = hydraulic_conductivity(idx, state, grid, hydrology, energy, strat, bgc)
end

@inline function hydraulic_conductivity(
    idx, grid, state,
    hydrology::SoilHydrology,
    energy::AbstractSoilEnergyBalance,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
)
    # Get porosity
    por = porosity(idx, state, hydrology, strat, bgc)

    # TODO
end

function compute_tendencies!(state, model, hydrology::SoilHydrology{NF, <:RichardsEq}) where {NF}
    grid = get_grid(model)
    hydrology = get_soil_hydrology(model)
    strat = get_stratigraphy(model)
    launch!(grid, :xyz, compute_saturation_tendency!, state, grid, hydrology, energy, strat, bgc)
    return nothing
end

@kernel function compute_saturation_tendency!(
    state,
    grid,
    hydrology::SoilHydrology,
    energy::AbstractSoilEnergyBalance,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
)
    idx = @index(Global, NTuple)
    state.tendencies.saturation_water_ice[idx...] += saturation_tendency(idx, state, grid, hydrology, energy, strat, bgc)
end

@inline function saturation_tendency(
    idx, state, grid,
    hydrology::SoilHydrology{NF, <:RichardsEq},
    energy::AbstractSoilEnergyBalance,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
) where {NF}
    i, j, k = idx
    # Operators require the underlying Oceananigans grid
    field_grid = get_field_grid(grid)

    # Compute divergence of water fluxes
    dθdt = ∂zᵃᵃᶜ(i, j, k, field_grid, darcy_flux, state, hydrology, energy, strat, bgc)
         + ∂zᵃᵃᶜ(i, j, k, field_grid, gravity_flux, state, hydrology, energy, strat, bgc)

    # Get porosity
    por = porosity(idx, state, hydrology, strat, bgc)

    return -dθdt / por
end

@inline function darcy_flux(i, j, k, grid, state, args...)
    # Get pressure field
    ψ = state.matric_potential
    # Get hydraulic conductivity
    K = state.hydraulic_conductivity
    # Darcy's law: q = -K ∂ψ/∂z
    q = -K * ∂zᵃᵃᶠ(i, j, k, grid, ψ)
    return q
end

# Matric potential <--> saturation closure relation

@kwdef struct PressureSaturationClosure{RetentionCurve<:SWRC} <: AbstractClosureRelation
    "Soil water retention curve(s) from FreezeCurves.jl"
    swrc::RetentionCurve = BrooksCorey()

    PressureSaturationClosure(swrc::SWRC) = new{typeof(swrc)}(ustrip(swrc))
end

getvar(::PressureSaturationClosure) = auxiliary(
    :saturation_water_ice,
    XYZ();
    domain=UnitInterval(),
    desc="Saturation level of water and ice in the pore space",
)
