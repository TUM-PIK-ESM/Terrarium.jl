"""
Base type for implementations of soil water flow dynamics.
"""
abstract type AbstractVerticalFlow end

"""
Represents a hydrology scheme where soil water is immobile.
"""
struct NoFlow <: AbstractVerticalFlow end

"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
struct SoilHydrology{
    NF,
    VerticalFlow<:AbstractVerticalFlow,
    SoilHydraulics<:AbstractSoilHydraulics{NF},
    VWCForcing<:ForcingType
} <: AbstractSoilHydrology{NF}
    "Soil water vertical flow operator"
    vertflow::VerticalFlow

    "Soil hydraulic properties parameterization"
    hydraulic_properties::SoilHydraulics

    "Forcing for soil moisture (volumetric water content)"
    forcing::VWCForcing
end

function SoilHydrology(
    ::Type{NF},
    vertflow::AbstractVerticalFlow = NoFlow();
    forcing::ForcingType = Forcings(),
    hydraulic_properties::AbstractSoilHydraulics = SoilHydraulicsSURFEX(NF),
) where {NF}
    return SoilHydrology(vertflow, hydraulic_properties, forcing)
end

"""
    get_swrc(hydrology::SoilHydrology)

Return the soil water retention curve from the `hydraulic_properties` associated with
the given `SoilHydrology` configuration.
"""
@inline get_swrc(hydrology::SoilHydrology) = hydrology.hydraulic_properties.cond_unsat.swrc

"""
    get_hydraulic_properties(hydrology::SoilHydrology)

Return the soil hydraulic properties defined by the given `hydrology` process.
"""
@inline get_hydraulic_properties(hydrology::SoilHydrology) = hydrology.hydraulic_properties

"""
    get_closure(::SoilHydrology{NF, NoFlow}) where {NF}

Return the saturation-pressure closure defined by the given `hydrology` process, or `nothing`
if not defined for the given configuration.
"""
@inline get_closure(::SoilHydrology{NF, NoFlow}) where {NF} = nothing

"""
State variables for `SoilHydrology` processes.
"""
variables(hydrology::SoilHydrology{NF}) where {NF} = (
    auxiliary(:saturation_water_ice, XYZ(), domain=UnitInterval(), desc="Saturation level of water and ice in the pore space"),
    auxiliary(:water_table, XY(), water_table, hydrology, units=u"m", desc="Elevation of the water table in meters"),
    input(:liquid_water_fraction, XYZ(), default = 1, domain=UnitInterval(), desc="Fraction of unfrozen water in the pore space"),
)

@inline saturation_water_ice(i, j, k, state, grid, ::AbstractSoilHydrology) = @inbounds state.saturation_water_ice[i, j, k]

@inline liquid_water_fraction(i, j, k, state, grid, ::AbstractSoilHydrology) = @inbounds state.liquid_water_fraction[i, j, k]

"""
Return a computed `Field` that derives the water table from the current `saturation_water_ice` state.
"""
function water_table(hydrology::SoilHydrology{NF}, grid, clock, fields) where {NF}
    is_saturated(i, j, k, grid, op) = @inbounds op.operand[i, j, k] == one(eltype(grid))
    z_faces = znodes(get_field_grid(grid), Center(), Center(), Face())
    # get height of soil column
    hz = abs(z_faces[1] - z_faces[end])
    # get total water + ice content in m³/m² by integrating saturation_water_ice
    sat = ConditionalOperation(fields.saturation_water_ice, condition = is_saturated)
    hw = Field(Integral(sat, dims=3))
    # take difference to get the water table elevation
    return hw - hz
end

# Immobile soil water (NoFlow)

@inline function initialize!(state, model, hydrology::SoilHydrology)
    set!(state.liquid_water_fraction, 1)
    return nothing
end

@inline compute_auxiliary!(state, model, hydrology::SoilHydrology) = nothing

@inline compute_tendencies!(state, model, hydrology::SoilHydrology{NF, NoFlow}) where {NF} = nothing
