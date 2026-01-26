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
    SaturationClosure<:AbstractSoilWaterClosure,
    SoilHydraulics<:AbstractSoilHydraulics{NF},
    VWCForcing<:Union{Nothing, AbstractForcing}
} <: AbstractSoilHydrology{NF}
    "Soil water vertical flow operator"
    vertflow::VerticalFlow

    "Closure relation for mapping between saturation water potential (hydraulic head)"
    closure::SaturationClosure

    "Soil hydraulic properties parameterization"
    hydraulic_properties::SoilHydraulics

    "Forcing for soil moisture (volumetric water content)"
    vwc_forcing::VWCForcing
end

function SoilHydrology(
    ::Type{NF},
    vertflow::AbstractVerticalFlow = NoFlow();
    closure::AbstractSoilWaterClosure = SaturationPressureClosure(),
    hydraulic_properties::AbstractSoilHydraulics = SoilHydraulicsSURFEX(NF),
    vwc_forcing::Union{Nothing, AbstractForcing} = nothing,
) where {NF}
    return SoilHydrology(vertflow, closure, hydraulic_properties, vwc_forcing)
end

"""
    get_swrc(hydrology::SoilHydrology)

Return the soil water retention curve from the `hydraulic_properties` associated with
the given `SoilHydrology` configuration.
"""
@inline get_swrc(hydrology::SoilHydrology) = hydrology.hydraulic_properties.swrc

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

@propagate_inbounds saturation_water_ice(i, j, k, grid, state, ::SoilHydrology) = state.saturation_water_ice[i, j, k]

@propagate_inbounds hydraulic_conductivity(i, j, k, grid, state, ::SoilHydrology) = state.hydraulic_conductivity[i, j, k]

@propagate_inbounds liquid_water_fraction(i, j, k, grid, state, ::SoilHydrology) = state.liquid_water_fraction[i, j, k]

@propagate_inbounds water_table(i, j, grid, state, ::SoilHydrology) = state.water_table[i, j]

@inline function compute_water_table!(state, grid, hydrology::SoilHydrology)
    launch!(grid, XY, compute_water_table_kernel!,
            state.water_table, state.saturation_water_ice, hydrology)
end

@inline function adjust_saturation_profile!(state, grid, hydrology::SoilHydrology)
    saturation_water_ice = state.satuation_water_ice
    surface_excess_water = state.surface_excess_water
    fields = (; saturation_water_ice, surface_excess_water)
    launch!(grid, XY, adjust_saturation_profile_kernel!, fields, hydrology)
end

# Immobile soil water (NoFlow)

@inline function initialize!(state, grid, hydrology::SoilHydrology)
    set!(state.liquid_water_fraction, 1)
    return nothing
end

@inline function compute_auxiliary!(state, grid, hydrology::SoilHydrology)
    compute_water_table!(state, get_grid(model), hydrology)
end

@inline function compute_tendencies!(state, grid, hydrology::SoilHydrology{NF, NoFlow, HP, <:AbstractForcing}) where {NF, HP}
    # TODO: do we need to also include ET? does that really make sense for a "no flow" scheme?
    forcing_kernel = KernelFunctionOperation{Center, Center, Center}(get_grid(model)) do i, j, k, grid
        forcing(i, j, k, grid, state, hydrology.vwc_forcing, hydrology)
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
@kernel function compute_water_table_kernel!(water_table, grid, sat, ::SoilHydrology{NF}) where {NF}
    i, j = @index(Global, NTuple)
    zs = znodes(get_field_grid(grid), Center(), Center(), Face())
    # scan z axis starting from the bottom (index 1) to find first non-saturated grid cell
    water_table[i, j, 1] = findfirst_z(i, j, <(one(NF)), zs, sat)
end

"""
    $TYPEDSIGNATURES

Kernel for adjusting saturation profiles to account for oversaturation due to numerical error.
This implementation scans over the saturation profiles at each lateral grid cell and redistributes
excess water upward layer-by-layer until reaching the topmost layer, where any remaining excess
water is added to the `surface_excess_water` pool.
"""
@kernel function adjust_saturation_profile_kernel!(state, grid, ::SoilHydrology{NF}) where {NF}
    i, j = @index(Global, NTuple)
    sat = state.saturation_water_ice
    surface_excess_water = state.surface_excess_water
    field_grid = get_field_grid(grid)
    N = field_grid.Nz
    # First iterate over soil layers from bottom to top
    # TODO: This function might perform badly on GPU....
    # Can we optimize it somehow?
    @inbounds for k in 1:N-1
        if sat[i, j, k] > one(NF)
            # calculate excess saturation
            excess_sat = sat[i, j, k] - one(NF)
            # subtract excess water and add to layer above;
            # note that we need to rescale by the cell thickness to properly conserve mass
            sat[i, j, k] -= excess_sat
            sat[i, j, k+1] += excess_sat * Δzᵃᵃᶜ(i, j, k, field_grid) / Δzᵃᵃᶜ(i, j, k+1, field_grid)
        end
    end
    # then from top to bottom
    @inbounds for k in N:-1:2
        if sat[i, j, k] < zero(NF)
            # calculate saturation deficit
            deficit_sat = -sat[i, j, k]
            # add back saturation deficit and subtract from layer below
            sat[i, j, k] += deficit_sat
            sat[i, j, k-1] -= deficit_sat * Δzᵃᵃᶜ(i, j, k, field_grid) / Δzᵃᵃᶜ(i, j, k-1, field_grid)
        end
    end
    @inbounds if sat[i, j, N] > one(NF)
        # If the uppermost (surface) layer is oversaturated, add to excess water pool
        excess_sat = sat[i, j, N] - one(NF)
        sat[i, j, N] -= excess_sat
        surface_excess_water[i, j, 1] += excess_sat * Δzᵃᵃᶜ(i, j, N, field_grid)
    end
    @inbounds if sat[i, j, 1] < zero(NF)
        # If the uppermost (surface) layer has a deficit, just set to zero.
        # This constitutes a mass balance violation but should not happen under realistic conditions.
        sat[i, j, 1] = zero(NF)
    end
end
