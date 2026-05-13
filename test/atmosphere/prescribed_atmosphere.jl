using Terrarium
using Test
import Thermodynamics as TD


@testset "Vapor pressure deficit" begin
    FT = Float64
    thermodynamic_constants = ThermodynamicConstants(FT)
    temperature_degC = 25.0
    temperature_K = Terrarium.celsius_to_kelvin(thermodynamic_constants, temperature_degC)
    pressure = 101325.0 # Pa
    total_specific_humidity = 0.01 # kg/kg
    vpd_TD = TD.vapor_pressure_deficit(thermodynamic_constants, temperature_K, pressure, total_specific_humidity)
    vpd_terrarium = Terrarium.vapor_pressure_deficit(thermodynamic_constants, temperature_degC, pressure, total_specific_humidity)
    @test vpd_TD ≈ vpd_terrarium
end
