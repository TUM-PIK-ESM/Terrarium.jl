# Compatibility dispatch for vegetation-free configurations
# TODO: This suggests a need for a redesign of the call interface
function compute_auxiliary!(
        state, grid,
        evaporation::AbstractEvapotranspiration,
        interception::NoCanopyInterception,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        soil::Optional{AbstractSoil},
        vegetation::Nothing,
        snow::Optional{AbstractSnow},
        args...
    )
    return compute_auxiliary!(state, grid, evaporation, interception, constants, atmos, soil, snow)
end

"""
    $TYPEDSIGNATURES

Compute an evapotranspiration flux (m/s, positive upwards) as the product of a vapor
conductance `g` (m/s) and a specific humidity difference `Δq` (kg/kg). All evapotranspiration components
— soil/canopy evaporation and transpiration — share the functional form ``Q_h = Δq · g``,
differing only in which conductance and humidity difference are supplied.
"""
@inline humidity_flux(::AbstractEvapotranspiration, Δq, g) = Δq * g

# Forcing interface for soil hydrology

"""
    $TYPEDSIGNATURES

Compute and return the evapotranspiration forcing for soil moisture at the given indices `i, j, k`.
The ET forcing is just the `ground_evapotranspiration_flux` rescaled by the thickness of layer `k`.
"""
@inline function forcing(i, j, k, grid, clock, fields, evapotranspiration::AbstractEvapotranspiration, ::AbstractSoilHydrology, args...)
    let Δz = Δzᵃᵃᶜ(i, j, k, grid)
        E = ground_evapotranspiration_flux(i, j, grid, fields, evapotranspiration) # liquid water flux [m³m⁻²s⁻¹]
        # Rescale by layer thickness and ratio of air to water density to get water content flux
        ∂θ∂t = -E / Δz
        return ∂θ∂t * (k == grid.Nz)
    end
end

# Fallback soil hydrology forcing for ET = nothing
forcing(i, j, k, grid, clock, fields, ::Nothing, ::AbstractSoilHydrology, args...) = zero(eltype(grid))
