"""
    $TYPEDEF

Canopy evapotranspiration scheme from PALADYN (Willeit 2016) that includes a canopy
re-evaporation term based on the saturation fraction of canopy water defined by the
canopy hydrology scheme.

```math
E_g &= \\beta \\frac{\\Delta q}{r_a} \\
E_c &= f_{\\text{can}} \\frac{\\Delta q}{r_a} \\
T_c &= \\frac{\\Delta q}{r_a + r_s} \\
```
"""
@kwdef struct PALADYNCanopyEvapotranspiration{NF, GR<:AbstractGroundEvaporationResistanceFactor} <: AbstractEvapotranspiration
    "Drag coefficient for the traansfer of heat and water between the ground and canopy"
    C_can::NF = 0.006

    "Parameterization for ground resistance to evaporation/sublimation"
    ground_resistance::GR = ConstantEvaporationResistance(typeof(C_can))
end

# TODO: The following ET functions all have the same basic functional form and can be generalized to a function
# that takes the humidity gradient/difference and an arbitrary number of resistance terms. We could consider
# making such an abstraction, but we should first consider what the main benefits would be.

"""
    $TYPEDSIGNATURES

Compute potential transpiration from the given humidity gradient, aerodynamic resistance `rₐ` and stomatal conductance `gw_can`.
"""
@inline function compute_transpiration(::PALADYNCanopyEvapotranspiration{NF}, Δq, rₐ, gw_can) where {NF}
    let rₛ = 1 / max(gw_can, sqrt(eps(NF))); # stomatal resistance as reciprocal of conductance
        # Calculate transpriation flux in m/s (positive upwards)
        E_tr = Δq / (rₐ + rₛ)
        return E_tr
    end
end

"""
    $TYPEDSIGNATURES

Compute potential evaporation from the ground below the canopy, following Eq. 5, PALADYN (Willeit 2016);
`Δq` is the humidity gradient, `β` is the ground evaporation resistance factor, `rₐ` is aerodynamic resistance,
and `rₑ` is aerodynamic resistance between the ground and canopy.
"""
@inline function compute_evaporation_ground(::PALADYNCanopyEvapotranspiration, Δq, β, rₐ, rₑ)
    # Calculate ground evaporation flux in m/s (positive upwards)
    E_gnd = β * Δq / (rₐ + rₑ)
    return E_gnd
end

"""
    $TYPEDSIGNATURES

Compute re-evaporation of water intercepted by the canopy from humidity gradient `Δq`, canopy saturation fraction
`f_can`, and aerodynamic resistance `rₐ`.
"""
@inline function compute_evaporation_canopy(::PALADYNCanopyEvapotranspiration, Δq, f_can, rₐ)
    # Calculate canopy re-evaporation flux in m/s (positive upwards)
    E_can = f_can * Δq / rₐ
    return E_can
end

# Process methods

variables(::PALADYNCanopyEvapotranspiration{NF}) where {NF} = (
    auxiliary(:evaporation_canopy, XY(); desc="Canopy re-evaporation contribution to surface humidity flux", units=u"m/s"),
    auxiliary(:evaporation_ground, XY(), units=u"m/s", desc="Ground evaporation contribution to surface humidity flux"),
    auxiliary(:transpiration, XY(), units=u"m/s", desc="Transpiration contribution to surface humidity flux"),
    input(:skin_temperature, XY(); units=u"°C", desc="Skin temperature"),
    input(:ground_temperature, XY(); default=NF(1), units=u"°C", desc="Ground surface temperature"),
    input(:gw_can, XY(); default=NF(1), units=u"m/s", desc="Canopy stomatal conductance")
)

function surface_humidity_flux(i, j, state, grid, evtr::PALADYNCanopyEvapotranspiration)
    @inbounds state.evaporation_canopy[i, j] + state.evaporation_ground[i, j] + state.transpiration[i, j]
end

function compute_auxiliary!(state, model, evtr::PALADYNCanopyEvapotranspiration)
    grid = get_grid(model)
    canopy_hydrology = get_canopy_hydrology(model)
    atmos = get_atmosphere(model)
    constants = get_constants(model)
    launch!(state, grid, :xy, compute_evapotranspiration_kernel!, canopy_hydrology, atmos, constants)
end

function compute_auxiliary!(state, model::AbstractLandModel, evtr::PALADYNCanopyEvapotranspiration)
    grid = get_grid(model)
    canopy_hydrology = get_canopy_hydrology(model)
    soilw = get_soil_hydrology(model)
    strat = get_soil_stratigraphy(model)
    bgc = get_soil_biogeochemistry(model)
    atmos = get_atmosphere(model)
    constants = get_constants(model)
    launch!(state, grid, :xy, compute_evapotranspiration_kernel!, canopy_hydrology, atmos, constants, soilw, strat, bgc)
end

# Kernels

# TODO: Can we somehow refactor this to reduce duplication across dispatches?
@kernel function compute_evapotranspiration_kernel!(
    state, grid,
    evtr::PALADYNCanopyEvapotranspiration{NF},
    canopy_hydrology::AbstractCanopyHydrology,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
) where {NF}
    i, j = @index(Global, NTuple)

    @inbounds let Ts = state.skin_temperature[i, j], # skin temperature (top of canopy)
                  Tg = state.ground_temperature[i, j], # ground temeprature (top snow/soil layer)
                  gw_can = state.gw_can[i, j], # stomatal conductance
                  Δqs = compute_humidity_vpd(i, j, state, grid, atmos, constants, Ts), # humidity gradient between canopy and atmosphere
                  Δqg = compute_humidity_vpd(i, j, state, grid, atmos, constants, Tg), # humidity gradient between ground and canopy
                  rₐ = aerodynamic_resistance(i, j, state, grid, atmos), # aerodynamic resistance
                  rₑ = aerodynamic_resistance(i, j, state, grid, atmos, evtr), # aerodynamic resistance between ground and canopy
                  f_can = saturation_canopy_water(i, j, state, grid, canopy_hydrology);

        # Compute and store ET fluxes
        state.transpiration[i, j, 1] = compute_transpiration(evtr, Δqs, rₐ, gw_can)
        state.evaporation_ground[i, j, 1] = compute_evaporation_ground(evtr, Δqg, NF(1), rₐ, rₑ)
        state.evaporation_canopy[i, j, 1] = compute_evaporation_canopy(evtr, Δqs, f_can, rₐ)
    end
end

@kernel function compute_evapotranspiration_kernel!(
    state, grid,
    evtr::PALADYNCanopyEvapotranspiration,
    canopy_hydrology::AbstractCanopyHydrology,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants,
    soilw::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
)
    i, j = @index(Global, NTuple)

    @inbounds let Ts = state.skin_temperature[i, j], # skin temperature (top of canopy)
                  Tg = state.ground_temperature[i, j], # ground temeprature (top snow/soil layer)
                  gw_can = state.gw_can[i, j], # stomatal conductance
                  Δqs = compute_humidity_vpd(i, j, state, grid, atmos, constants, Ts), # humidity gradient between canopy and atmosphere
                  Δqg = compute_humidity_vpd(i, j, state, grid, atmos, constants, Tg), # humidity gradient between ground and canopy
                  rₐ = aerodynamic_resistance(i, j, state, grid, atmos), # aerodynamic resistance
                  rₑ = aerodynamic_resistance(i, j, state, grid, atmos, evtr), # aerodynamic resistance between ground and canopy
                  f_can = saturation_canopy_water(i, j, state, grid, canopy_hydrology),
                  β = ground_evaporation_resistance_factor(i, j, state, grid, evtr.ground_resistance, soilw, strat, bgc);

        # Compute and store ET fluxes
        state.transpiration[i, j, 1] = compute_transpiration(evtr, Δqs, rₐ, gw_can)
        state.evaporation_ground[i, j, 1] = compute_evaporation_ground(evtr, Δqg, β, rₐ, rₑ)
        state.evaporation_canopy[i, j, 1] = compute_evaporation_canopy(evtr, Δqs, f_can, rₐ)
    end
end

# Kernel functions

"""
    aerodynamic_resistance(i, j, state, atmos::AbstractAtmosphere, evtr::PALADYNCanopyEvapotranspiration)

Compute the aerodynamic resistance between the ground and canopy as a function of LAI and SAI.
"""
@inline function aerodynamic_resistance(i, j, state, grid, atmos::AbstractAtmosphere, evtr::PALADYNCanopyEvapotranspiration)
    @inbounds let LAI = state.LAI[i, j],
                  SAI = state.SAT[i, j],
                  Vₐ = windspeed(i, j, state, grid, atmos),
                  C = ET.C_can; # drag coefficient for the canopy
        rₙ = (1 - exp(-LAI - SAI)) / (C * Vₐ)
        return rₙ
    end
end
