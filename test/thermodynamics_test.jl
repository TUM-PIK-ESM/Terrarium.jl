using Pkg
Pkg.activate("test/.")
using Terrarium
using ClimaParams
using Thermodynamics
using Test


@testset "Equivalence test with Clima" begin
    FT = Float64
    air_temperature = 10.0 # 10 °C
    air_density = 1.225 # kg/m³
    air_pressure = 101_325 # Pa
    air_temperature_K = air_temperature + 273.15 # Convert to Kelvin
    # Clima
    params_clima = Thermodynamics.Parameters.ThermodynamicsParameters(FT)
    e_sat_clima = Thermodynamics.saturation_vapor_pressure(params_clima, air_temperature_K)
    q_sat_clima = Thermodynamics.q_vap_saturation(params_clima, air_temperature_K, air_density)
    # Terrarium
    thermodyn_constants = Terrarium.ThermodynamicConstants(FT)
    e_sat_terrarium = Terrarium.saturation_vapor_pressure(thermodyn_constants, air_temperature)
    q_sat_terrarium = Terrarium.saturation_specific_humidity_vapor(thermodyn_constants, air_temperature, air_density)
    # Compare
    @test isapprox(e_sat_clima, e_sat_terrarium; rtol = 0.01)
    @test isapprox(q_sat_clima, q_sat_terrarium; rtol = 0.01)
end

@testset "Saturated conditions" begin
    air_temperature = 10.0 # 10 °C
    air_density = 1.225 # kg/m³
    air_pressure = 101_325 # Pa
    thermodyn_constants = Terrarium.ThermodynamicConstants()
    # Manual calcuation
    e_sat_tetens = 610.78 * exp(17.27 * air_temperature / (air_temperature + 237.3)) # Pa
    q_sat_tetens = 0.62 * e_sat_tetens / (air_pressure - 0.38 * e_sat_tetens) # kg/kg
    # Using Terrarium function
    e_sat = Terrarium.saturation_vapor_pressure(thermodyn_constants, air_temperature)
    q_sat = Terrarium.saturation_specific_humidity_vapor(thermodyn_constants, air_temperature, air_density)
    # Compare with 10% tolerance
    @test isapprox(e_sat, e_sat_tetens; rtol = 0.1)
    @test isapprox(q_sat, q_sat_tetens; rtol = 0.1)
end
