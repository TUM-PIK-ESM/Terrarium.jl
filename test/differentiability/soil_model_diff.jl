using Terrarium
using Test

using Enzyme
using FreezeCurves
using Oceananigans: Average, Field
using Statistics

grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N=10))
# initial conditions
initializer = FieldInitializers(
    # steady-ish state initial condition for temperature
    temperature = (x,z) -> -1 - 0.02*z,
    # saturated soil
    saturation_water_ice = 1.0,
)
model = SoilModel(; grid, initializer)
sim = initialize(model)
timestep!(sim)

state = sim.state
dstate = make_zero(state)

function dostep!(state, model, timestepper, Δt)
    timestep!(state, model, timestepper, Δt)
    return mean(interior(state.temperature))
    # TODO: Figure out why this is segfaulting in Enzyme
    # Answer: Average operator is not type inferrable, see:
    # https://github.com/CliMA/Oceananigans.jl/issues/4869
    # Tavg = Field(Average(state.temperature, dims=(1, 2, 3)))
    # return Tavg[1,1,1]
end

dostep!(state, model, model.time_stepping, 1.0)

@testset "Soil model: timestep!" begin
    @time Enzyme.autodiff(set_runtime_activity(Reverse), dostep!, Active, Duplicated(state, dstate), Const(model), Const(model.time_stepping), Const(model.time_stepping.Δt))
    @test all(isfinite.(dstate.temperature))
end

@testset "Soil energy: FreeWater freeze-thaw" begin
    # test liquid_water_fraction
    U = -1e7
    sat = 1.0
    por = 0.5
    Lθ = 3.34e8*sat*por
    lwc_grad, = Enzyme.autodiff(Reverse, Terrarium.liquid_water_fraction, Active, Const(FreeWater()), Active(U), Const(Lθ), Const(sat))
    @test lwc_grad[2] ≈ 1 / Lθ
    ## check case where Lθ is zero
    Lθ = 0.0
    lwc_grad, = Enzyme.autodiff(Reverse, Terrarium.liquid_water_fraction, Active, Const(FreeWater()), Active(U), Const(Lθ), Const(sat))
    @test iszero(lwc_grad[2])

    # test energy_to_temperature
    ## case 1: frozen
    sat = 1.0
    por = 0.5
    Lθ = 3.34e8*sat*por
    U = Lθ - 1e7
    C = 2e5
    inv_energy_grad, = Enzyme.autodiff(Reverse, Terrarium.energy_to_temperature, Active, Const(FreeWater()), Active(U), Const(Lθ), Const(C))
    @test inv_energy_grad[2] ≈ 1 / C
    ## case 2: phase change
    sat = 1.0
    por = 0.5
    Lθ = 3.34e8*sat*por
    U = -Lθ/2
    C = 2e5
    inv_energy_grad, = Enzyme.autodiff(Reverse, Terrarium.energy_to_temperature, Active, Const(FreeWater()), Active(U), Const(Lθ), Const(C))
    @test iszero(inv_energy_grad[2])
    ## case 3: thawed
    sat = 1.0
    por = 0.5
    Lθ = 3.34e8*sat*por
    U = Lθ/2
    C = 2e5
    inv_energy_grad, = Enzyme.autodiff(Reverse, Terrarium.energy_to_temperature, Active, Const(FreeWater()), Active(U), Const(Lθ), Const(C))
    @test inv_energy_grad[2] ≈ 1 / C
end
