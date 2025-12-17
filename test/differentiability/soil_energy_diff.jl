using Terrarium
using Test

using Enzyme, FiniteDifferences
using FreezeCurves
using Statistics


function build_soil_energy_model(arch, ::Type{NF}) where {NF}
    grid = ColumnGrid(arch, Float64, ExponentialSpacing(N=10))
    # initial conditions
    initializer = FieldInitializers(
        # steady-ish state initial condition for temperature
        temperature = (x,z) -> -1 - 0.02*z,
        # saturated soil
        saturation_water_ice = 1.0,
    )
    model = SoilModel(; grid, initializer)
    return model
end

function mean_soil_temperature_step!(state, timestepper, model, inputs, Δt)
    timestep!(state, timestepper, model, inputs, Δt)
    return mean(interior(state.temperature))
    # TODO: Figure out why this is segfaulting in Enzyme
    # Answer: Average operator is not type inferrable, see:
    # https://github.com/CliMA/Oceananigans.jl/issues/4869
    # Tavg = Field(Average(state.temperature, dims=(1, 2, 3)))
    # return Tavg[1,1,1]
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

@testset "Soil energy model: timestep!" begin
    model = build_soil_energy_model(CPU(), Float64)
    integrator = initialize(model, ForwardEuler())
    inputs = integrator.inputs
    state = integrator.state
    dstate = make_zero(state)
    stepper = integrator.timestepper
    dstepper = make_zero(stepper)
    @time Enzyme.autodiff(set_runtime_activity(Reverse), mean_soil_temperature_step!, Active, Duplicated(state, dstate), Duplicated(stepper, dstepper), Const(model), Const(integrator.inputs), Const(integrator.timestepper.Δt))
    @test all(isfinite.(dstate.temperature))
end

model = build_soil_energy_model(CPU(), Float64)
dmodel = make_zero(model)
integrator = initialize(model, ForwardEuler())
inputs = integrator.inputs
state = integrator.state
dstate = make_zero(state)
stepper = integrator.timestepper
dstepper = make_zero(stepper)
@time Enzyme.autodiff(
    set_runtime_activity(Reverse),
    mean_soil_temperature_step!,
    Active,
    Duplicated(state, dstate),
    Duplicated(stepper, dstepper),
    Duplicated(model, dmodel),
    Const(integrator.inputs),
    Const(integrator.timestepper.Δt)
)

@inline function thermalconductivity(props::SoilThermalProperties)#, soil::SoilComposition)
    # κs = Terrarium.getproperties(props.cond)
    # fracs = volumetric_fractions(soil)
    return props.cond.mineral*2
    # return sum(fastmap((x, w) -> sqrt(x)*w, κs, fracs))^2
    # apply bulk conductivity weighting
    # return props.cond_bulk(κs, fracs)
end

props = SoilThermalProperties(Float64)
dprops = make_zero(props)
#comp = SoilComposition()
#dcomp = make_zero(comp)
thermalconductivity(props, comp)
Enzyme.autodiff(Reverse, thermalconductivity, Active, Duplicated(props, dprops))#, Duplicated(comp, dcomp))


@kwdef struct Foo{A,B} 
    a::A
    b::B
end  

get_field_arg(foo::Foo) = 2*(foo.a + foo.b)

bar = Foo(a=3.0, b=2.0)
dbar = make_zero(bar)
Enzyme.autodiff(Reverse, get_field_arg, Active, Duplicated(bar, dbar))
