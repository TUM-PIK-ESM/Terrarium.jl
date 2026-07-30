"""
    $TYPEDEF

Default representation of coupled surface hydrology processes including
canopy rain/snow interception, evapotranspiration, and surface runoff.

Properties:
$FIELDS
"""
struct SurfaceHydrology{
        NF,
        CanopyInterception <: AbstractCanopyInterception,
        Evapotranspiration <: AbstractEvapotranspiration,
        SurfaceRunoff <: AbstractSurfaceRunoff,
    } <: AbstractSurfaceHydrology{NF}
    "Canopy hydrology scheme"
    canopy_interception::CanopyInterception

    "Canopy evapotranspiration scheme"
    evapotranspiration::Evapotranspiration

    "Surface runoff scheme"
    surface_runoff::SurfaceRunoff
end

# No sub-scheme is pinned to `NF` — any of them may hold a promoted (e.g. traced) parameter — so `NF`
# is not derivable from the field types and must be given. `constructorof` rebuilds through this
# constructor so `setproperties` / `Flatten.reconstruct` / `Adapt` retain `NF`.
SurfaceHydrology{NF}(
    canopy_interception::AbstractCanopyInterception,
    evapotranspiration::AbstractEvapotranspiration,
    surface_runoff::AbstractSurfaceRunoff,
) where {NF} = SurfaceHydrology{NF, typeof(canopy_interception), typeof(evapotranspiration), typeof(surface_runoff)}(
    canopy_interception, evapotranspiration, surface_runoff
)

ConstructionBase.constructorof(::Type{<:SurfaceHydrology{NF}}) where {NF} = SurfaceHydrology{NF}

function SurfaceHydrology(
        ::Type{NF};
        canopy_interception::AbstractCanopyInterception = PALADYNCanopyInterception(NF),
        evapotranspiration::AbstractEvapotranspiration = PALADYNCanopyEvapotranspiration(NF),
        surface_runoff::AbstractSurfaceRunoff = DirectSurfaceRunoff(NF)
    ) where {NF}
    return SurfaceHydrology{NF}(canopy_interception, evapotranspiration, surface_runoff)
end

""" $TYPEDSIGNATURES """
function compute_auxiliary!(
        state, grid,
        hydrology::SurfaceHydrology,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        soil::Optional{AbstractSoil} = nothing,
        vegetation::Optional{AbstractVegetation} = nothing,
        args...
    )
    compute_auxiliary!(state, grid, hydrology.canopy_interception, atmos)
    compute_auxiliary!(state, grid, hydrology.evapotranspiration, hydrology.canopy_interception, constants, atmos, soil, vegetation)
    compute_auxiliary!(state, grid, hydrology.surface_runoff, hydrology.canopy_interception, soil)
    return nothing
end

""" $TYPEDSIGNATURES """
function compute_tendencies!(
        state, grid,
        hydrology::SurfaceHydrology,
        args...,
    )
    # Compute tendencies for canopy interception
    compute_tendencies!(state, grid, hydrology.canopy_interception, hydrology.evapotranspiration)
    return nothing
end
