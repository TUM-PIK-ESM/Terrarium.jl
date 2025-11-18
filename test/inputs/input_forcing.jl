module ForcingInputTest

using Terrarium
using Terrarium: AbstractLandGrid, prognostic, input
using Test

DEFAULT_NF = Float32

@kwdef struct TestModel{NF, Grid<:AbstractLandGrid{NF}, TS} <: Terrarium.AbstractModel{NF, Grid, TS}
    grid::Grid
    initializer = DefaultInitializer()
    time_stepping::TS = ForwardEuler{DEFAULT_NF}()
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
    model = TestModel(; grid)
    F_in = FieldInputSource(eltype(grid), :F)
    inputs = InputSources(F_in)
    integrator = initialize(model, ForwardEuler(eltype(grid)); inputs)
    # check initial values
    @test all(integrator.state.x .≈ 0)
    @test all(integrator.state.F .≈ 1)
    # advance one timestep and check updated values
    timestep!(integrator, 0.1)
    @test all(integrator.state.x .≈ 0.1)
    @test all(integrator.state.F .≈ 1)
end

end