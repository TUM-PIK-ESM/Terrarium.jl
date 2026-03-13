using Terrarium
using Test

# mock a simple model with exponential dynamics (and a constant offset) to test time steppers

@kwdef struct ExpModel{NF, Grid <: Terrarium.AbstractLandGrid{NF}, I} <: Terrarium.AbstractModel{NF, Grid}
    grid::Grid
    initializer::I = DefaultInitializer(eltype(grid))
end

Terrarium.variables(::ExpModel) = (
    Terrarium.prognostic(:u, Terrarium.XY()),
    Terrarium.auxiliary(:v, Terrarium.XY()),
)

# just a constant offset (we could do it differently but this is for testing auxilitary as wel)
function Terrarium.compute_auxiliary!(state, model::ExpModel)
    return state.auxiliary.v .= 0.1
end

# du/dt = u + c = u + 0.1
function Terrarium.compute_tendencies!(state, model::ExpModel)
    return set!(state.tendencies.u, state.prognostic.u + state.auxiliary.v)
end

@testset "ExpModel: Heun and Euler time steppers" begin

    grid = ColumnGrid(CPU(), Float64, UniformSpacing(N = 1))
    model = ExpModel(grid)

    initializers = (u = 0.0, v = 0.1)
    integrator_heun = initialize(model, Heun(); initializers)
    integrator_euler = initialize(model, ForwardEuler(); initializers)

    # test that Heun estimate is more accurate (larger value than Euler here)
    # test that both are what we expect
    timestep!(integrator_heun)
    timestep!(integrator_euler)

    @test integrator_heun.state.u[2] > integrator_euler.state.u[2]

    # Euler: expected value: u = 0.1 * Δt
    dt_euler = default_dt(integrator_euler.timestepper)
    @test integrator_euler.state.u[2] == 0.1 * dt_euler

    # Heun: expected value: u = (0.1Δt + (0.1 * Δt + 0.1) * Δt) / 2
    dt_heun = default_dt(integrator_heun.timestepper)
    @test integrator_heun.state.u[2] == (0.1 * dt_heun + (0.1 * dt_heun + 0.1) * dt_heun) / 2
end

# Use timestep!(state, model, timestepper, Δt) to clip negative values in an super simple example sim
@testset "ExpModel: clip negative values" begin
    grid = ColumnGrid(CPU(), Float64, UniformSpacing(N = 1))
    model = ExpModel(grid)

    Terrarium.timestep!(state, model::ExpModel, timestepper::ForwardEuler, Δt) = begin
        state.u[2] = max(state.u[2], 0.0)
    end

    initializers = (u = -20, v = -5.0)
    integrator = initialize(model, ForwardEuler(); initializers)

    # Test that timestep! clips negative values
    timestep!(integrator)

    @test integrator.state.u[2] >= 0.0
end
