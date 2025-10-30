using Terrarium
using Test

using Enzyme, FiniteDifferences
using FreezeCurves
using Statistics

function build_soil_energy_hydrology_model(arch, ::Type{NF}, flow=RichardsEq(); hydrology_kwargs...) where {NF}
    grid = ColumnGrid(arch, Float64, ExponentialSpacing(N=10))
    # initial conditions
    initializer = FieldInitializers(
        # steady-ish state initial condition for temperature
        temperature = (x,z) -> 1.0,
        # saturated soil
        saturation_water_ice = (x,z) -> min(0.5 - 0.1*z, 1.0),
    )
    hydrology = SoilHydrology(eltype(grid), flow; hydrology_kwargs...)
    model = SoilModel(grid; initializer, hydrology)
    return model
end

@testset "Soil water: Pressure-saturation closure" begin
    por = 0.5
    sat = 0.5
    hydraulic_properties = ConstantHydraulics(Float64, porosity=por)
    model = build_soil_energy_hydrology_model(CPU(), Float64; hydraulic_properties)
    swrc = Terrarium.get_swrc(model.hydrology) # θ(ψₘ)
    swrc_inv = inv(swrc) # ψₘ(θ)
    model_state = initialize(model)
    state = model_state.state
    dstate = make_zero(state)
    set!(state.saturation_water_ice, sat)
    set!(dstate.pressure_head, 1.0) # seed pressure head (output)
    closure = Terrarium.PressureSaturationClosure()
    # first check inverse relation saturation -> pressure
    Enzyme.autodiff(set_runtime_activity(Reverse), Terrarium.invclosure!, Const, Duplicated(state, dstate), Const(model), Const(closure))
    # check that derivatives w.r.t saturation are finite
    @test all(isfinite.(dstate.saturation_water_ice))
    # check that they match analytical derivative (1 / ∂θ∂ψ by inverse function theorem)
    ψm = swrc_inv(sat*por; θsat=por)
    @test all(dstate.saturation_water_ice .≈ por / swrc(FreezeCurves.derivative, ψm; θsat=por))

    # now check pressure -> saturation
    dstate = make_zero(state)
    set!(dstate.saturation_water_ice, 1.0) # seed saturation (output)
    Enzyme.autodiff(set_runtime_activity(Reverse), Terrarium.closure!, Const, Duplicated(state, dstate), Const(model), Const(closure))
    # check that gradients w.r.t saturation are finite
    @test all(isfinite.(dstate.pressure_head))
    # check that saturation is correct
    @test all(state.saturation_water_ice .≈ sat)
    # check derivatives
    @test all(dstate.pressure_head .≈ swrc(FreezeCurves.derivative, ψm; θsat=por) / por)
end

@testset "Soil hydrology: hydraulic_conductivity" begin
    por = 0.5
    cond_unsat = UnsatKVanGenuchten(Float64)
    hydraulic_properties = ConstantHydraulics(Float64; porosity=por, cond_unsat)
    # wrapper function for evaluating hydraulic conductivity
    function eval_hydraulic_cond((por, sat, liq))
        soil = SoilComposition(porosity=por, saturation=sat, liquid=liq, organic=0.0, texture=SoilTexture())
        return Terrarium.hydraulic_conductivity(hydraulic_properties, soil)
    end

    x = (0.5, 0.75, 0.9)
    grads, = Enzyme.autodiff(Reverse, eval_hydraulic_cond, Active, Active(x))
    @test isapprox(grads[1][1], 0, rtol=1e-8, atol=1e-8)
    fd_grads = FiniteDifferences.grad(central_fdm(5, 1), eval_hydraulic_cond, collect(x))
    @test isapprox(collect(grads[1]), fd_grads[1], rtol=1e-8, atol=1e-8)
end

@testset "Soil hydrology: compute_auxiliary! RRE" begin
    hydraulic_properties = SoilHydraulicsSURFEX(Float64)
    model = build_soil_energy_hydrology_model(CPU(), Float64; hydraulic_properties)
    model_state = initialize(model)
    state = model_state.state
    dstate = make_zero(state)
    set!(dstate.hydraulic_conductivity, 1.0) # seed hydraulic cond
    Enzyme.autodiff(set_runtime_activity(Reverse), compute_auxiliary!, Const, Duplicated(state, dstate), Const(model), Const(model.hydrology))
    @test all(isfinite.(dstate.saturation_water_ice))
    @test all(isfinite.(dstate.temperature))
end

@testset "Soil hydrology: compute_auxiliary! RRE" begin
    hydraulic_properties = SoilHydraulicsSURFEX(Float64)
    model = build_soil_energy_hydrology_model(CPU(), Float64; hydraulic_properties)
    model_state = initialize(model)
    state = model_state.state
    # first run compute_auxiliary! for the full model (needed to compute hydraulic conductivities)
    compute_auxiliary!(state, model)
    dstate = make_zero(state)
    set!(dstate.tendencies.saturation_water_ice, 1.0) # seed tendencies
    Enzyme.autodiff(set_runtime_activity(Reverse), compute_tendencies!, Const, Duplicated(state, dstate), Const(model), Const(model.hydrology))
    @test all(isfinite.(dstate.pressure_head))
    @test dstate.pressure_head[1,1,1] > 0 # higher pressure -> weaker gradient -> less outflow
end

@testset "Soil energy/hydrology model: timestep!" begin
    hydraulic_properties = ConstantHydraulics(Float64)
    model = build_soil_energy_hydrology_model(CPU(), Float64; hydraulic_properties)
    model_state = initialize(model)
    state = model_state.state
    dstate = make_zero(state)
    @time Enzyme.autodiff(set_runtime_activity(Reverse), timestep!, Const, Duplicated(state, dstate), Const(model), Const(model.time_stepping))
    @test all(isfinite.(dstate.temperature))
    @test all(isfinite.(dstate.pressure_head))
end
