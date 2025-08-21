using Terrarium
using Test

using Oceananigans: Field

@testset "AtmosphericState" begin
    grid = ColumnGrid(ExponentialSpacing())
    # use nonzero test values to check that everything gets updated correctly
    precip = TwoPhasePrecipitation(rainfall=1.0, snowfall=1.0)
    solar = TwoBandSolarRadiation(shortwave=1.0, longwave=1.0)
    T_air = 1.0
    humidity = 1.0
    pressure = 1.0
    windspeed = 1.0
    tracers = TracerGases(AmbientCO2(1.0))
    atm = AtmosphericState(; grid, precip, solar, T_air, humidity, pressure, windspeed, tracers)
    @test haskey(atm.tracers, :CO2)
    atm_varnames = map(Terrarium.varname, variables(atm))
    @test all(map(∈(atm_varnames), (:T_air, :humidity, :pressure, :windspeed, :rainfall, :snowfall, :SwIn, :LwIn)))
    # check compute_auxiliary!
    state = (
        T_air = Field{Center,Center,Nothing}(get_field_grid(grid)),
        humidity = Field{Center,Center,Nothing}(get_field_grid(grid)),
        pressure = Field{Center,Center,Nothing}(get_field_grid(grid)),
        windspeed = Field{Center,Center,Nothing}(get_field_grid(grid)),
        rainfall = Field{Center,Center,Nothing}(get_field_grid(grid)),
        snowfall = Field{Center,Center,Nothing}(get_field_grid(grid)),
        SwIn = Field{Center,Center,Nothing}(get_field_grid(grid)),
        LwIn = Field{Center,Center,Nothing}(get_field_grid(grid)),
        CO2 = Field{Center,Center,Nothing}(get_field_grid(grid)),
    )
    compute_auxiliary!(state, nothing, atm)
    @test all(state.T_air .≈ 1)
    @test all(state.humidity .≈ 1)
    @test all(state.pressure .≈ 1)
    @test all(state.windspeed .≈ 1)
    @test all(state.rainfall .≈ 1)
    @test all(state.snowfall .≈ 1)
    @test all(state.SwIn .≈ 1)
    @test all(state.LwIn .≈ 1)
    @test all(state.CO2 .≈ 1)
end
