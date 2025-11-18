# mock a simple model with exponential dynamics (and a constant offset) to test time steppers 

struct ExpModel{NF, Grid, I} <: Terrarium.AbstractModel{NF, Grid}
    grid::Grid
    initializer::I
end

ExpModel(grid::Grid, initializer::I) where {Grid, I} = ExpModel{eltype(grid), Grid, I}(grid, initializer)

Terrarium.variables(::ExpModel) = (Terrarium.prognostic(:u, Terrarium.XY()), 
                              Terrarium.auxiliary(:v, Terrarium.XY()))

# just a constant offset (we could do it differently but this is for testing auxilitary as wel)
function Terrarium.compute_auxiliary!(state, model::ExpModel) 
    state.auxiliary.v .= 0.1
end 

# du/dt = u + c = u + 0.1
function Terrarium.compute_tendencies!(state, model::ExpModel) 
    state.tendencies.u .= state.prognostic.u + state.auxiliary.v
end 

@testset "ExpModel: Heun and Euler time steppers" begin 

    grid = ColumnGrid(CPU(), Float64, UniformSpacing(N=1))
    initializer = FieldInitializers(u = 0., v = 0.1)
    model = ExpModel(grid, initializer)

    integrator_heun = initialize(model, Heun())
    integrator_euler = initialize(model, ForwardEuler())

    # test that Heun estimate is more accurate (larger value than Euler here)
    # test that both are what we expect 
    timestep!(integrator_heun)
    timestep!(integrator_euler)

    @test integrator_heun.state.u[2] > integrator_euler.state.u[2]

    # Euler: expected value: u = 0.1*Δt
    dt_euler = default_dt(integrator_euler.timestepper)
    @test integrator_euler.state.u[2] == 0.1*dt_euler

    # Heun: expected value: u = (0.1Δt + (0.1*Δt+0.1)* Δt)/2
    dt_heun = default_dt(integrator_heun.timestepper)
    @test integrator_heun.state.u[2] == (0.1*dt_heun + (0.1*dt_heun+0.1)*dt_heun)/2
end 