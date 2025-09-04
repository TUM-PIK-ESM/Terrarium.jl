using Terrarium
using Test

using Oceananigans: Field

@testset "PrescribedAtmosphere" begin
    grid = ColumnGrid(ExponentialSpacing())
    precip = TwoPhasePrecipitation()
    solar = TwoBandSolarRadiation()
    tracers = TracerGases(AmbientCO2())
    atm = PrescribedAtmosphere(; grid, precip, solar, T_air, humidity, pressure, windspeed, tracers)
    @test haskey(atm.tracers, :CO2)
    atm_varnames = map(Terrarium.varname, variables(atm))
    @test all(map(∈(atm_varnames), (:T_air, :humidity, :pressure, :windspeed, :rainfall, :snowfall, :SwIn, :LwIn, :CO2)))
end
