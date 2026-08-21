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
— soil/canopy evaporation and transpiration — share the functional form ``E = Δq · g``,
differing only in which conductance and humidity difference are supplied.
"""
@inline compute_evaporation_flux(::AbstractEvapotranspiration, Δq, g) = Δq * g

# Forcing interface for soil hydrology

"""
    forcing(i, j, k, grid, clock, fields, evapotranspiration::AbstractEvapotranspiration, ::AbstractSoilHydrology)

Compute and return the evapotranspiration forcing for soil moisture at the given indices `i, j, k`.
The ET forcing is just the `surface_humidity_flux` rescaled by the thickness of layer `k`.
"""
@inline function forcing(i, j, k, grid, clock, fields, ::AbstractEvapotranspiration, ::AbstractSoilHydrology)
    let Δz = Δzᵃᵃᶜ(i, j, k, grid)
        Q_E = fields.evaporation_soil_water[i, j]
        # Rescale by layer thickness and ratio of air to water density to get water content flux
        ∂θ∂t = -Q_E / Δz
        return ∂θ∂t * (k == grid.Nz)
    end
end
