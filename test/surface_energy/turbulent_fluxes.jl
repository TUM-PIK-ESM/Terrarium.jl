using Terrarium
using Test

@testset "Prescribed turbulent fluxes" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N=10))
    turbulent_fluxes = PrescribedTurbulentFluxes()
    surface_energy_balance = SurfaceEnergyBalance(Float64; turbulent_fluxes)
    model = SurfaceEnergyModel(grid, surface_energy_balance)
    model_state = initialize(model)
    state = model_state.state
    @test hasproperty(state.inputs, :sensible_heat_flux)
    @test hasproperty(state.inputs, :latent_heat_flux)
    set!(state.sensible_heat_flux, 10.0)
    set!(state.latent_heat_flux, 5.0)
    @test Terrarium.sensible_heat_flux((1,1), state, turbulent_fluxes) == 10.0
    @test Terrarium.latent_heat_flux((1,1), state, turbulent_fluxes) == 5.0    
end

@testset "Diagnosed turbulent fluxes" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N=10))
    turbulent_fluxes = DiagnosedTurbulentFluxes(Float64)
    surface_energy_balance = SurfaceEnergyBalance(Float64; turbulent_fluxes)
    model = SurfaceEnergyModel(grid, surface_energy_balance)
    model_state = initialize(model)
    state = model_state.state
    @test !hasproperty(state.inputs, :sensible_heat_flux)
    @test !hasproperty(state.inputs, :latent_heat_flux)
    @test hasproperty(state.auxiliary, :sensible_heat_flux)
    @test hasproperty(state.auxiliary, :latent_heat_flux)
    set!(state.skin_temperature, 10.0)
    set!(state.air_temperature, 5.0)
    compute_auxiliary!(state, model, turbulent_fluxes)
    # check that sensible heat fluxes are > 0 (air colder than skin, positive up)
    @test all(state.sensible_heat_flux .> 0)
    set!(state.skin_temperature, 5.0)
    set!(state.air_temperature, 10.0)
    set!(state.specific_humidity, 0.5)
    compute_auxiliary!(state, model, turbulent_fluxes)
    # check that sensible heat fluxes are < 0 (air warmer than skin, negative down)
    @test all(state.sensible_heat_flux .< 0)
end
