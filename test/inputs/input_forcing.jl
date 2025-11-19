module ForcingInputTest

using Terrarium
using Terrarium: AbstractLandGrid, prognostic, input
using Test

DEFAULT_NF = Float32

@kwdef struct TestModel{NF, Grid<:AbstractLandGrid{NF}, TS} <: Terrarium.AbstractModel{NF, Grid}
    grid::Grid
    initializer = DefaultInitializer()
end

Terrarium.variables(model::TestModel) = (
    prognostic(:x, XYZ()),
    input(:F, XY()),
)

function Terrarium.compute_tendencies!(state, model::TestModel)
    # set tendency to forcing term
    state.tendencies.x .= state.F
end

function Terrarium.timestep!(state, model::TestModel, euler::ForwardEuler, Δt)
    Terrarium.compute_tendencies!(state, model)
    @. state.x += Δt*state.tendencies.x
end

@testset "Forcing input" begin
    grid = ColumnGrid(CPU(), DEFAULT_NF, ExponentialSpacing())
    model = TestModel(grid)
    F_in = FieldInputSource(eltype(grid), :F)
    integrator = initialize(model, ForwardEuler(eltype(grid)), F_in)
    # check initial values
    @test all(integrator.state.x .≈ 0)
    @test all(integrator.state.F .≈ 1)
    # advance one timestep and check updated values
    timestep!(integrator, 0.1)
    @test all(integrator.state.x .≈ 0.1)
    @test all(integrator.state.F .≈ 1)
end

@testset "Forcing input with time series" begin
    grid = ColumnGrid(CPU(), DEFAULT_NF, ExponentialSpacing())
    model = TestModel(; grid)
    t_F = 0:0.1:1
	F = FieldTimeSeries(grid, XY(), t_F)
	F.data .= ones(size(F));
    F_in = FieldTimeSeriesInputSource(; F)
    sim = initialize(model, ForwardEuler(eltype(grid)), F_in)
    # check initial values
    @test all(sim.state.x .≈ 0)
    @test all(sim.state.F .≈ 1)
    # advance one timestep and check updated values
    timestep!(sim, 0.1)
    @test all(sim.state.x .≈ 0.1)
    @test all(sim.state.F .≈ 1)
end

end