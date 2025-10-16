"""
Base type for implementations of soil water flow dynamics.
"""
abstract type AbstractSoilWaterFluxOperator <: AbstractOperator end

"""
Represents a hydrology scheme where soil water is immobile.
"""
struct NoFlow <: AbstractSoilWaterFluxOperator end

"""
Base type for evapotranspirative fluxes in soil layers.
"""
abstract type AbstractSoilET{NF} end

"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
struct SoilHydrology{
    NF,
    Operator<:AbstractSoilWaterFluxOperator,
    SoilET<:Union{Nothing, AbstractSoilET{NF}},
    SoilHydraulics<:AbstractSoilHydraulics{NF}
} <: AbstractSoilHydrology{NF}
    "Soil water flux operator"
    operator::Operator

    "Soil evapotranspiration scheme"
    evapotranspiration::SoilET

    "Soil hydraulic properties parameterization"
    hydraulic_properties::SoilHydraulics
end

SoilHydrology(
    ::Type{NF},
    operator::AbstractSoilWaterFluxOperator = NoFlow();
    evapotranspiration::AbstractSoilET = nothing,
    hydraulic_properties::AbstractSoilHydraulics = SoilHydraulicsSURFEX(NF),
) where {NF} = SoilHydrology(operator, evapotranspiration, hydraulic_properties)

"""
    get_swrc(hydrology::SoilHydrology)

Return the soil water retention curve from the `hydraulic_properties` associated with
the given `SoilHydrology` configuration.
"""
get_swrc(hydrology::SoilHydrology) = hydrology.hydraulic_properties.cond_unsat.swrc

"""
    porosity(idx, state, hydrology::SoilHydrology, strat::AbstractStratigraphy, bgc::AbstractSoilBiogeochemistry)

Return the porosity of the soil volume at `idx` given the current state, hydrology, stratigraphy, and biogeochemistry configurations.
"""
@inline function porosity(idx, state, hydrology::SoilHydrology, strat::AbstractStratigraphy, bgc::AbstractSoilBiogeochemistry)
    org = organic_fraction(idx, state, bgc)
    texture = soil_texture(idx, state, strat)
    return (1 - org)*mineral_porosity(hydrology.hydraulic_properties, texture) + org*organic_porosity(idx, state, bgc)
end

# Immobile soil water (NoFlow)

variables(::SoilHydrology{NF, NoFlow}) where {NF} = (
    auxiliary(:saturation_water_ice, XYZ(), domain=UnitInterval(), desc="Saturation level of water and ice in the pore space"),
    auxiliary(:water_table, XY(), units=u"m", desc="Elevation of the water table [m]"),
)

function initialize!(state, model, hydrology::SoilHydrology)
    # Since water content does not change in the NoFlow scheme, we just compute the water table
    # once at initialization time
    grid = get_grid(model)
    launch!(grid, :xy, compute_water_table!, state, grid, hydrology)
end

@inline compute_auxiliary!(state, model, hydrology::SoilHydrology) = nothing

@inline compute_tendencies!(state, model, hydrology::SoilHydrology{NF, NoFlow}) where {NF} = nothing

"""
    compute_water_table!(
        state,
        grid,
        ::SoilHydrology{NF}
    ) where {NF}

Kernel for diagnosing the water table at each grid point given the current soil saturation state.
"""
@kernel function compute_water_table!(
    state,
    grid,
    ::SoilHydrology{NF}
) where {NF}
    i, j = @index(Global, NTuple)
    sat = state.saturation_water_ice
    # get z coordinates of grid cell faces
    zs = znodes(get_field_grid(grid), Center(), Center(), Face())
    # scan z axis starting from the bottom (index 1) to find first non-saturated grid cell
    state.water_table[i, j, 1] = findfirst_z((i, j), <(one(NF)), zs, sat)
end
