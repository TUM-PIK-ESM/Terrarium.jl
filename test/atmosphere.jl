using Terrarium
using Test

using Oceananigans: Field

@testset "PrescribedAtmosphere" begin
    grid = ColumnGrid(ExponentialSpacing())
    precip = RainSnow()
    solar = TwoBandSolarRadiation()
    tracers = TracerGases(AmbientCO2())
    atm = PrescribedAtmosphere(; grid, precip, solar, T_air, humidity, pressure, windspeed, tracers)
    @test haskey(atm.tracers, :CO2)
    atm_varnames = map(Terrarium.varname, variables(atm))
    @test all(map(∈(atm_varnames), (:T_air, :humidity, :pressure, :windspeed, :rainfall, :snowfall, :surface_shortwave_down, :surface_longwave_down, :CO2)))
end
