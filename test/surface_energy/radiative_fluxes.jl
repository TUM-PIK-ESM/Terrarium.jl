using Terrarium
using Test

@testset "Prescribed radiative fluxes" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N=10))
    radiative_fluxes = PrescribedRadiativeFluxes()
    model = SurfaceEnergyBalanceModel(grid; radiative_fluxes)
    model_state = initialize(model)
    state = model_state.state
    @test hasproperty(state.inputs, :SwOut)
    @test hasproperty(state.inputs, :LwOut)
    set!(state.SwIn, 100.0)
    set!(state.LwIn, 20.0)
    set!(state.SwOut, 50.0)
    set!(state.LwOut, 5.0)
    compute_auxiliary!(state, model, model.radiative_fluxes)
    @test all(state.RadNet .≈ 50.0 - 100.0 + 5.0 - 20.0)
end

@testset "Diagnosed radiative fluxes" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N=10))
    radiative_fluxes = DiagnosedRadiativeFluxes()
    albedo = ConstantAlbedo(albedo=0.5, emissivity=0.9)
    model = SurfaceEnergyBalanceModel(grid; radiative_fluxes, albedo)
    model_state = initialize(model)
    state = model_state.state
    @test !hasproperty(state.inputs, :SwOut)
    @test !hasproperty(state.inputs, :LwOut)
    @test hasproperty(state.auxiliary, :SwOut)
    @test hasproperty(state.auxiliary, :LwOut)
    SwIn = 100.0
    LwIn = 20.0
    set!(state.SwIn, SwIn)
    set!(state.LwIn, LwIn)
    compute_auxiliary!(state, model, model.radiative_fluxes)
    @test all(state.SwOut .≈ 0.5*SwIn)
    @test all(state.LwOut .≈ (1 - 0.9)*LwIn + Terrarium.stefan_boltzmann(model.constants, 273.15, 0.9))
    @test all(state.RadNet .≈ state.SwOut - SwIn + state.LwOut - LwIn)
end
