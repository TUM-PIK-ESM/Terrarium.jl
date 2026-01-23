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
    VWCForcing<:Union{Nothing, AbstractForcing}
} <: AbstractSoilHydrology{NF}
    "Soil water vertical flow operator"
    vertflow::VerticalFlow

    "Soil hydraulic properties parameterization"
    hydraulic_properties::SoilHydraulics

    "Forcing for soil moisture (volumetric water content)"
    vwc_forcing::VWCForcing
end

function SoilHydrology(
    ::Type{NF},
    vertflow::AbstractVerticalFlow = NoFlow();
    vwc_forcing::Union{Nothing, AbstractForcing} = nothing,
    hydraulic_properties::AbstractSoilHydraulics = SoilHydraulicsSURFEX(NF),
) where {NF}
    return SoilHydrology(vertflow, hydraulic_properties, vwc_forcing)
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
    auxiliary(:water_table, XY(), units=u"m", desc="Elevation of the water table in meters"),
    auxiliary(:hydraulic_conductivity, XYZ(z=Face()), units=u"m/s", desc="Hydraulic conductivity of soil volumes in m/s"),
    input(:liquid_water_fraction, XYZ(), default = 1, domain=UnitInterval(), desc="Fraction of unfrozen water in the pore space"),
)

@propagate_inbounds saturation_water_ice(i, j, k, state, grid, ::AbstractSoilHydrology) = state.saturation_water_ice[i, j, k]

@propagate_inbounds hydraulic_conductivity(i, j, k, state, grid, ::AbstractSoilHydrology) = state.hydraulic_conductivity[i, j, k]

@propagate_inbounds liquid_water_fraction(i, j, k, state, grid, ::AbstractSoilHydrology) = state.liquid_water_fraction[i, j, k]

@propagate_inbounds water_table(i, j, state, grid, ::AbstractSoilHydrology) = state.water_table[i, j, 1]

@inline function compute_water_table!(state, grid, hydrology::AbstractSoilHydrology)
    zs = znodes(get_field_grid(grid), Center(), Center(), Face())
    launch!(grid, :xy, compute_water_table_kernel!, state.water_table, state.saturation_water_ice, zs, hydrology)
end

# Immobile soil water (NoFlow)

@inline function initialize!(state, model, hydrology::SoilHydrology)
    set!(state.liquid_water_fraction, 1)
    return nothing
end

@inline function compute_auxiliary!(state, model, hydrology::SoilHydrology)
    compute_water_table!(state, get_grid(model), hydrology)
end

@inline function compute_tendencies!(state, model, hydrology::SoilHydrology{NF, NoFlow, HP, <:AbstractForcing}) where {NF, HP}
    forcing_kernel = KernelFunctionOperation{Center, Center, Center}(get_grid(model)) do i, j, k, grid
        forcing(i, j, k, state, grid, hydrology.vwc_forcing, hydrology)
    end
    # apply forcing
    set!(state.saturation_water_ice, forcing_kernel)
end

# Kernels

"""
    compute_water_table_kernel!(water_table, sat, z_faces, ::SoilHydrology{NF}) where {NF}

Kernel for diagnosing the water table at each grid point given the current soil saturation profile.
The argument `z_faces` should be the z-coordinates of the grid on the layer faces.
"""
@kernel function compute_water_table_kernel!(water_table, sat, z_faces, ::SoilHydrology{NF}) where {NF}
    i, j = @index(Global, NTuple)
    # scan z axis starting from the bottom (index 1) to find first non-saturated grid cell
    water_table[i, j, 1] = findfirst_z((i, j), <(one(NF)), z_faces, sat)
end
