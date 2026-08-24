using Terrarium
using ClimaParams
using Thermodynamics
using Test


@testset "Equivalence test with Clima" begin
    NF = Float64
    air_temperature = 10.0 # 10 °C
    air_density = 1.225 # kg/m³
    air_pressure = 101_325 # Pa
    air_temperature_K = air_temperature + 273.15 # Convert to Kelvin
    # Clima
    params_clima = Thermodynamics.Parameters.ThermodynamicsParameters(NF)
    # Terrarium
    thermodyn_constants = Terrarium.ThermodynamicConstants(NF)
    # Test all constants are identical
    # Note that for R_d and R_v we use constants with slightly higher precision than ClimaParams
    @test Thermodynamics.Parameters.R_d(params_clima) ≈ round(Thermodynamics.Parameters.R_d(thermodyn_constants), digits = 0)
    @test Thermodynamics.Parameters.R_v(params_clima) ≈ round(Thermodynamics.Parameters.R_v(thermodyn_constants), digits = 1)
    @test Thermodynamics.Parameters.cp_d(params_clima) ≈ Thermodynamics.Parameters.cp_d(thermodyn_constants)
    @test Thermodynamics.Parameters.cp_i(params_clima) ≈ Thermodynamics.Parameters.cp_i(thermodyn_constants)
    @test Thermodynamics.Parameters.cp_l(params_clima) ≈ Thermodynamics.Parameters.cp_l(thermodyn_constants)
    @test Thermodynamics.Parameters.cp_v(params_clima) ≈ Thermodynamics.Parameters.cp_v(thermodyn_constants)
    @test Thermodynamics.Parameters.LH_v0(params_clima) ≈ Thermodynamics.Parameters.LH_v0(thermodyn_constants)
    @test Thermodynamics.Parameters.T_0(params_clima) ≈ Thermodynamics.Parameters.T_0(thermodyn_constants)
    @test Thermodynamics.Parameters.T_freeze(params_clima) ≈ Thermodynamics.Parameters.T_freeze(thermodyn_constants)
    @test Thermodynamics.Parameters.T_triple(params_clima) ≈ Thermodynamics.Parameters.T_triple(thermodyn_constants)
    @test Thermodynamics.Parameters.press_triple(params_clima) ≈ Thermodynamics.Parameters.press_triple(thermodyn_constants)
    # Test that latent heat of sublimation equals the sum of L_f and L_v; note that this differs slightly from Clima in that
    # we derive L_s from L_f and L_v rather than L_f from L_s and L_v.
    @test Thermodynamics.Parameters.LH_s0(thermodyn_constants) ≈ Thermodynamics.Parameters.LH_f0(thermodyn_constants) + Thermodynamics.Parameters.LH_v0(thermodyn_constants)
    # Check that the difference with Clima params is negligible (<1 kJ)
    @test abs(Thermodynamics.Parameters.LH_s0(params_clima) - Thermodynamics.Parameters.LH_s0(thermodyn_constants)) < 1.0e3
    # Test saturation vapor pressure functions
    # For these, we explicitly set the Terrarium value for R_v to the same numeric precision as Clima to check that the resulting values match
    R_v = round(thermodyn_constants.gas_constant_water_vapor, digits = 1)
    e_sat_terrarium = Terrarium.saturation_vapor_pressure(Terrarium.ThermodynamicConstants(NF, gas_constant_water_vapor = R_v), air_temperature)
    q_sat_terrarium = Terrarium.saturation_specific_humidity_vapor(Terrarium.ThermodynamicConstants(NF, gas_constant_water_vapor = R_v), air_temperature, air_density)
    e_sat_clima = Thermodynamics.saturation_vapor_pressure(params_clima, air_temperature_K)
    q_sat_clima = Thermodynamics.q_vap_saturation(params_clima, air_temperature_K, air_density)
    @test e_sat_clima ≈ e_sat_terrarium
    @test q_sat_clima ≈ q_sat_terrarium
end

@testset "Saturated conditions" begin
    air_temperature = 10.0 # 10 °C
    air_density = 1.225 # kg/m³
    air_pressure = 101_325 # Pa
    thermodyn_constants = Terrarium.ThermodynamicConstants()
    # Manual calculation with Tetens formula
    e_sat_tetens = 610.78 * exp(17.27 * air_temperature / (air_temperature + 237.3)) # Pa
    q_sat_tetens = 0.62 * e_sat_tetens / (air_pressure - 0.38 * e_sat_tetens) # kg/kg, Eq. 2.8 Shuttleworth (2012)
    # Using Terrarium function
    e_sat = Terrarium.saturation_vapor_pressure(thermodyn_constants, air_temperature)
    q_sat = Terrarium.saturation_specific_humidity_vapor(thermodyn_constants, air_temperature, air_density)
    # Compare with 5% tolerance
    @test isapprox(e_sat, e_sat_tetens; rtol = 0.05)
    @test isapprox(q_sat, q_sat_tetens; rtol = 0.05)
end
