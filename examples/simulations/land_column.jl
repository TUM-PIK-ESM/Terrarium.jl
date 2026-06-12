using Pkg
Pkg.activate("examples/.")
using Terrarium
using CUDA

import CairoMakie as Makie

arch = CPU()
# Define a simple grid with 1 column
grid = ColumnGrid(arch, ExponentialSpacing(Δz_max = 1.0, N = 30))
# Set up Richards model for soil hydrology
swrc = VanGenuchten(α = 2.0, n = 2.0)
hydraulic_properties = ConstantSoilHydraulics(eltype(grid); swrc, unsat_hydraulic_cond = UnsatKVanGenuchten(eltype(grid)))
hydrology = SoilHydrology(eltype(grid), RichardsEq(); hydraulic_properties)
soil = SoilEnergyWaterCarbon(eltype(grid); hydrology)
vegetation = VegetationCarbon(eltype(grid))
# Construct coupled model
land = LandModel(grid; soil, vegetation)
# Variably saturated with water table at roughly 5 m depth
initializers = (
    temperature = 15.0,
    saturation_water_ice = (x, z) -> min(1, 0.5 - 0.1 * z),
    carbon_vegetation = 0.5,
)
integrator = @time initialize(land, ForwardEuler(); initializers);
# plot initial conditions
zs = znodes(integrator.state.temperature)
Makie.plot(interior(integrator.state.temperature)[1,1,:], zs)
Makie.plot(interior(integrator.state.saturation_water_ice)[1,1,:], zs)
# manually set atmospheric inputs to different values
atmospheric_windspeed = 5.0 # m/s
atmospheric_specific_humidity = 1.0e-4 # kg/kg
set!(integrator.state.windspeed, atmospheric_windspeed) # 1 m/s
set!(integrator.state.specific_humidity, atmospheric_specific_humidity) # kg/kg
Δt = 60.0
timestep!(integrator, Δt)
@show integrator.state.latent_heat_flux[1, 1, 1]

#

# Check logic in surface energy balance
# Radiative fluxes
SW_up = integrator.state.surface_shortwave_up[1, 1, 1]
LW_up = integrator.state.surface_longwave_up[1, 1, 1]
SW_down = integrator.state.surface_shortwave_down[1, 1, 1]
LW_down = integrator.state.surface_longwave_down[1, 1, 1]

SW_up ≈ integrator.model.surface_energy_balance.albedo.albedo * SW_down
R_net_positive_down = SW_down + LW_down - SW_up - LW_up
R_net_positive_down > 0
R_net_positive_up = SW_up + LW_up - SW_down - LW_down
R_net_positive_up < 0 
R_net_positive_up ≈ R_n
# Turbulent fluxes
R_n = integrator.state.surface_net_radiation[1, 1, 1]
H = integrator.state.sensible_heat_flux[1, 1, 1]
LE = integrator.state.latent_heat_flux[1, 1, 1]
G = integrator.state.ground_heat_flux[1, 1, 1]

# Compare temperatures
T_air = integrator.state.air_temperature[1, 1, 1]
T_skin = integrator.state.skin_temperature[1, 1, 1]
T_soil_upper = integrator.state.temperature[1, 1, 1]

# Manual flux calculations
ρ_air = 1.225 # kg/m^3 
c_p_air = 1005.0 # J/(kg*K)
r_a = 1 / (integrator.model.atmosphere.aerodynamics.Cₕ[1, 1, 1] * atmospheric_windspeed) # s/m
H_manual = ρ_air * c_p_air * (T_skin - T_air) / r_a
(H_manual - H) < 1 

L = 2.5e6 # J/kg
e_sat = 610.78 * exp(17.27 * T_air / (T_air + 237.3)) # Pa
q_sat = 0.622 * e_sat / 101325 # kg/kg
Δq = q_sat - atmospheric_specific_humidity
LE_manual = ρ_air * L * Δq / r_a
if (LE_manual - LE) < 1
    diff = LE_manual - LE
    println("Latent heat manual - Terrarium: $diff W/m^2")
end
