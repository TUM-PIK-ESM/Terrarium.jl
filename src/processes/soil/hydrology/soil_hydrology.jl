"""
Base type for implementations of soil water flow dynamics.
"""
abstract type AbstractVerticalFlow end

"""
Represents a hydrology scheme where soil water is immobile.
"""
struct NoFlow <: AbstractVerticalFlow end

"""
Base type for soil/subsurface runoff schemes.
"""
abstract type AbstractSoilRunoff end

"""
Base type for soil/subsurface evapotranspiration forcing terms.
"""
abstract type AbstractSoilET end

struct NoSoilET <: AbstractSoilET end

"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
struct SoilHydrology{
    NF,
    VerticalFlow<:AbstractVerticalFlow,
    Runoff<:AbstractSoilRunoff,
    SoilET<:AbstractSoilET,
    SoilHydraulics<:AbstractSoilHydraulics{NF},
} <: AbstractSoilHydrology{NF}
    "Soil water vertical flow operator"
    vertflow::VerticalFlow

    "Soil subsurface runoff scheme"
    runoff::Runoff

    "Scheme for distributing ET fluxes in the root zone"
    evapotranspiration::SoilET

    "Soil hydraulic properties parameterization"
    hydraulic_properties::SoilHydraulics
end

function SoilHydrology(
    ::Type{NF},
    vertflow::AbstractVerticalFlow = NoFlow();
    runoff::AbstractSoilRunoff = SoilRunoff(),
    evapotranspiration::AbstractSoilET = NoSoilET(),
    hydraulic_properties::AbstractSoilHydraulics = SoilHydraulicsSURFEX(NF),
) where {NF}
    return SoilHydrology(vertflow, runoff, evapotranspiration, hydraulic_properties)
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
    variables(hydrology.vertflow)...,
    variables(hydrology.runoff)...,
    variables(hydrology.evapotranspiration)...,
)

# Immobile soil water (NoFlow)

variables(::NoFlow) = (
    auxiliary(:saturation_water_ice, XYZ(), domain=UnitInterval(), desc="Saturation level of water and ice in the pore space"),
)

@inline saturation_water_ice(i, j, k, state, hydrology::AbstractSoilHydrology) = @inbounds state.saturation_water_ice[i, j, k]

@inline initialize!(state, model, hydrology::SoilHydrology) = nothing

@inline compute_auxiliary!(state, model, hydrology::SoilHydrology{NF, NoFlow}) where {NF} = nothing

@inline compute_tendencies!(state, model, hydrology::SoilHydrology{NF, NoFlow}) where {NF} = nothing

# Runoff

"""
    $TYPEDEF

Generic scheme for respresenting soil/subsurface runoff as a composition of three components:
excess infiltration, interflow, and drainage.

Not yet implemented.
"""
@kwdef struct SoilRunoff{SR, IF, DR} <: AbstractSoilRunoff
    "Excess infiltration runoff scheme"
    excess_infiltration::SR = nothing
    
    "Subsurface interflow runoff"
    interflow::IF = nothing

    "Groundwater drainage runoff"
    drainage::DR = nothing
end

# Evapotranspiration

"""
    SurfaceEvaporation <: AbstractSoilET

Evapotranspiration scheme for bare soils that allocates the full latent heat flux to evaporation
from the topmost soil layer.
"""
struct SurfaceEvaporation <: AbstractSoilET end

variables(::SurfaceEvaporation) = (
    input(:latent_heat_flux, XY(), units=u"W/m^2", desc="Latent heat flux at the surface [W m⁻²]"),
)

@inline function forcing_ET(i, j, k, grid, state, ::SurfaceEvaporation, constants::PhysicalConstants)
    let Hₗ = state.latent_heat_flux[i, j],
        L = constants.Lsg, # specific latent heat of vaporization
        ρw = constants.ρw, # density of water
        Δz = Δzᵃᵃᶜ(i, j, k, grid);
        q_E = Hₗ / (L * ρw) # ET water flux in m s⁻¹
        ∂θ∂t = -q_E / Δz # rescale by layer thickness to get water content flux
        return ∂θ∂t * (k == grid.Nz) # only nonzero at the surface
    end
end

# Default to zero if evapotranspiration is disabled
@inline forcing_ET(i, j, k, grid, state, ::NoSoilET, args...) = zero(eltype(grid))
