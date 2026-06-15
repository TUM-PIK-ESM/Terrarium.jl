using Terrarium
using Test


@testset "Saturated conditions" begin
    air_temperature = 10.0 # 10 °C
    air_density = 1.225 # kg/m³
    air_pressure = 101_325 # Pa
    thermodyn_constants = Terrarium.ThermodynamicConstants()
    # Manual calcuation
    e_sat_tetens = 610.78 * exp(17.27 * air_temperature / (air_temperature + 237.3)) # Pa
    q_sat_tetens = 0.62 * e_sat_tetens / (air_pressure - 0.38 * e_sat_tetens) # kg/kg
    # Using Terrarium function
    e_sat_clima = Terrarium.saturation_vapor_pressure(thermodyn_constants, air_temperature)
    q_sat_clima = Terrarium.saturation_specific_humidity_vapor(thermodyn_constants, air_temperature, air_pressure)
    # Compare with 10% tolerance
    @test isapprox(e_sat_clima, e_sat_tetens; rtol = 0.1)
    @test isapprox(q_sat_clima, q_sat_tetens; rtol = 0.1)
end
