# mock a simple model with exponential dynamics (and a constant offset) to test time steppers 

@kwdef struct ExpModel{NF, Grid, I, BC} <: Terrarium.AbstractModel{NF, Grid}
    grid::Grid
    initializer::I = DefaultInitializer()
    boundary_conditions::BC = DefaultBoundaryConditions()
end

ExpModel(grid::Grid, initializer::I, boundary_conditions::BC) where {Grid, I, BC} = ExpModel{eltype(grid), Grid, I, BC}(grid, initializer, boundary_conditions)

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
    model = ExpModel(; grid, initializer)

    modelstate_heun = initialize(model; timestepper=Heun)
    modelstate_euler = initialize(model, timestepper=ForwardEuler)

    # test that Heun estimate is more accurate (larger value than Euler here)
    # test that both are what we expect 
    timestep!(modelstate_heun)
    timestep!(modelstate_euler)

    @test modelstate_heun.state.u[2] > modelstate_euler.state.u[2]

    # Euler: expected value: u = 0.1*Δt
    @test modelstate_euler.state.u[2] == 0.1*default_dt(modelstate_euler.timestepper)

    # Heun: expected value: u = (0.1Δt + (0.1*Δt+0.1)* Δt)/2
    @test modelstate_heun.state.u[2] == (0.1*default_dt(modelstate_heun.timestepper) + (0.1*default_dt(modelstate_heun.timestepper)+0.1)* default_dt(modelstate_heun.timestepper))/2
end 