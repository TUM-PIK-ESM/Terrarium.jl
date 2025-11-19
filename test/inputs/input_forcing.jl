module ForcingInputTest

using Terrarium
using Terrarium: AbstractLandGrid, prognostic, input
using Test

DEFAULT_NF = Float32

@kwdef struct TestModel{NF, Grid <: AbstractLandGrid{NF}, TS} <: Terrarium.AbstractModel{NF, Grid, TS}
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
    return state.tendencies.x .= state.F
end

function Terrarium.timestep!(state, model::TestModel, euler::ForwardEuler, Δt)
    Terrarium.compute_tendencies!(state, model)
    return @. state.x += Δt * state.tendencies.x
end

@testset "Forcing input" begin
    grid = ColumnGrid(CPU(), DEFAULT_NF, ExponentialSpacing())
    model = TestModel(; grid)
    F = Field(grid, XY())
    set!(F, one(eltype(F)))
    F_in = FieldInputSource(; F)
    sim = initialize(model, F_in)
    # check initial values
    @test all(sim.state.x .≈ 0)
    @test all(sim.state.F .≈ 1)
    # advance one timestep and check updated values
    timestep!(sim, 0.1)
    @test all(sim.state.x .≈ 0.1)
    @test all(sim.state.F .≈ 1)
end

end
