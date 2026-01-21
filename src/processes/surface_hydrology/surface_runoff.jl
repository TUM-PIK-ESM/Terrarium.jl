"""
    $TYPEDEF

Canopy interception and storage implementation following PALADYN (Willeit 2016) considering only liquid water (no snow).

Properties:
$FIELDS
"""
@kwdef struct PALADYNSurfaceRunoff{NF} <: AbstractSurfaceRunoff
    "Damping parameter for saturated fraction of the grid cell"
    f∇::NF = 1.7
end

PALADYNSurfaceRunoff(::Type{NF}; kwargs...) where {NF} = PALADYNSurfaceRunoff{NF}(; kwargs...)

variables(::PALADYNSurfaceRunoff) = (
    auxiliary(:runoff_surface, XY(), units=u"m/s", desc="Surface runoff"),
    auxiliary(:f_sat, XY(), units=u"m/s", desc="Fraction of the grid cell with saturated ground"),
    auxiliary(:infiltration, XY(), units=u"m/s", desc="Infiltration flux"),
)

function compute_auxiliary!(state, model, canopy_hydrology::PALADYNSurfaceRunoff)
    grid = get_grid(model)
    atmos = get_atmosphere(model)
    constants = get_constants(model)
    launch!(state, grid, :xy, compute_auxiliary_kernel!, canopy_hydrology, atmos, constants)
end

# Kernels

@kernel function compute_auxiliary_kernel!(
    state,
    runoff::PALADYNSurfaceRunoff,
    canopy_hydrology::AbstractCanopyHydrology,
    soil_hydrology::AbstractSoilHydrology
)
    
    i, j = @index(Global, NTuple)

    # Get inputs 
    precip_ground = ground_precipitation(i, j, state, grid, canopy_hydrology)

    #TODO: WIP
end
