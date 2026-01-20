"""
    $TYPEDEF

Canopy evapotranspiration scheme that includes a canopy re-evaporation term based on the
saturation fraction of canopy water defined by the canopy hydrology scheme.
"""
@kwdef struct CanopyEvapotranspiration{NF, GR<:AbstractEvaporationResistance} <: AbstractEvapotranspiration
    "Drag coefficient for the traansfer of heat and water between the ground and canopy"
    C_can::NF = 0.006

    "Parameterization for ground resistance to evaporation/sublimation"
    ground_resistance::GR = ConstantEvaporationResistance(typeof(C_can))
end

variables(::CanopyEvapotranspiration{NF}) where {NF} = (
    auxiliary(:evaporation_canopy, XY(); desc="Canopy re-evaporation contribution to surface humidity flux", units=u"m/s"),
    auxiliary(:evaporation_ground, XY(), units=u"m/s", desc="Ground evaporation contribution to surface humidity flux"),
    auxiliary(:transpiration, XY(), units=u"m/s", desc="Transpiration contribution to surface humidity flux"),
    input(:gw_can, XY(); default=NF(1), units=u"m/s", desc="Canopy stomatal conductance")
)

function compute_auxiliary!(state, model, ET::CanopyEvapotranspiration)
    grid = get_grid(model)
    canopy_hydrology = get_canopy_hydrology(model)
    soilw = get_soil_hydrology(model)
    strat = get_soil_stratigraphy(model)
    bgc = get_soil_biogeochemistry(model)
    seb = get_surface_energy_balance(model)
    atmos = get_atmosphere(model)
    constants = get_constants(model)
    launch!(state, grid, :xy, compute_evapotranspiration_kernel!, canopy_hydrology, seb.skin_temperature, atmos, soilw, strat, bgc, constants)
end

# Kernels

@kernel function compute_evapotranspiration_kernel!(
    state, grid,
    ET::CanopyEvapotranspiration,
    canopy_hydrology::AbstractCanopyHydrology,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    soilw::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants
)
    i, j = @index(Global, NTuple)

    @inbounds begin
        state.transpiration[i, j, 1] = compute_transpiration(i, j, state, ET, skinT, atmos, constants)
        state.evaporation_ground[i, j, 1] = compute_evaporation_ground(i, j, state, ET, skinT, atmos, soilw, strat, bgc, constants)
        state.evaporation_canopy[i, j, 1] = compute_evaporation_canopy(i, j, state, ET, canopy_hydrology, atmos, constants)
    end
end

# Kernel functions

"""
Compute potential transpiration based on the current skin temperature and .
"""
@inline function compute_transpiration(
    i, j, state,
    ET::CanopyEvapotranspiration,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    @inbounds let Ts = skin_temperature(i, j, state, grid, skinT), # skin temperature (top of canopy)
        Δq = compute_humidity_vpd(i, j, state, atmos, constants, Ts), # humidity gradient
        rₐ = aerodynamic_resistance(i, j, state, grid, atmos), # aerodynamic resistance of air
        rₛ = 1 / state.gw_can[i, j]; # stomatal resistance as reciprocal of conductance
        # Calculate potential transpriation flux in m/s (positive upwards)
        E_tr = Δq / (rₐ + rₛ)
        return E_tr
    end
end

"""
Compute potential evaporation from the ground below the canopy, following Eq. 5, PALADYN (Willeit 2016),
excepting the soil moisture limiting factor, which is applied later.
"""
@inline function compute_evaporation_ground(
    i, j, state,
    ET::CanopyEvapotranspiration,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    soilw::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants
)
    let Ts = ground_temperature(i, j, state, skinT), # temperature of top soil/snow layer
        rₐ = aerodynamic_resistance(i, j, state, grid, atmos), # aerodynamic resistance of air
        rᵥ = aerodynamic_resistance(i, j, state, grid, atmos, ET), # aerodynamic resistance of canopy
        β = ground_evaporation_resistance_factor(i, j, state, grid, ET.ground_resistance, soilw, strat, bgc)
        Δq = compute_humidity_vpd(i, j, state, atmos, constants, Ts); # humidity gradient
        # Calculate potential ground evaporation flux in m/s (positive upwards)
        E_gnd = β * Δq / (rₐ + rᵥ)
        return E_gnd
    end
end

"""
Compute re-evaporation of water intercepted by the canopy.
"""
@inline function compute_evaporation_canopy(
    i, j, state,
    ET::CanopyEvapotranspiration,
    canopy_hydrology::AbstractCanopyHydrology,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    let Ts = ground_temperature(i, j, state, skinT), # temperature of top soil/snow layer
        Δq = compute_humidity_vpd(i, j, state, atmos, constants, Ts), # humidity gradient
        rₐ = aerodynamic_resistance(i, j, state, grid, atmos), # aerodynamic resistance of air
        f_can = saturation_water_canopy(i, j, state, grid, canopy_hydrology);
        # Calculate canopy re-evaporation flux in m/s (positive upwards)
        E_can = f_can * Δq / rₐ
        return E_can
    end
end

"""
    aerodynamic_resistance(i, j, state, atmos::AbstractAtmosphere, ET::CanopyEvapotranspiration)

Compute the aerodynamic resistance between the ground and canopy as a function of LAI and SAI.
"""
@inline function aerodynamic_resistance(i, j, state, grid, atmos::AbstractAtmosphere, ET::CanopyEvapotranspiration)
    @inbounds let LAI = state.LAI[i, j],
                  SAI = state.SAT[i, j],
                  Vₐ = windspeed(i, j, state, atmos),
                  C = ET.C_can; # drag coefficient for the canopy
        rₙ = (1 - exp(-LAI - SAI)) / (C * Vₐ)
        return rₙ
    end
end
