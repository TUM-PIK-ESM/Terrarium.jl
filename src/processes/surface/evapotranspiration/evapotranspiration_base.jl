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

# Fallback soil hydrology forcing for ET = nothing
forcing(i, j, k, grid, fields, ::Nothing, ::AbstractSoilHydrology, args...) = zero(eltype(grid))
