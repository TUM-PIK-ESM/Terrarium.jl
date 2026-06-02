using Terrarium
using Test

using Terrarium: AbstractLandGrid, AbstractImplicitTimestepper, prognostic, XY

module IMEXTestTypes

    using Terrarium
    using Terrarium: AbstractLandGrid, AbstractImplicitTimestepper, prognostic, XY

    # A mock implicit timestepper used only to verify IMEX routing. Its update is deliberately distinct
    # from forward Euler (u += 2·∂u∂t·Δt) so we can tell which sub-stepper integrated which variable.
    struct MockImplicit{NF} <: AbstractImplicitTimestepper{NF}
        Δt::NF
    end
    MockImplicit(::Type{NF}; Δt = 300.0) where {NF} = MockImplicit{NF}(NF(Δt))

    Terrarium.default_dt(ts::MockImplicit) = ts.Δt
    Terrarium.is_adaptive(::MockImplicit) = false

    function Terrarium.timestep!(integrator, ts::MockImplicit, Δt, names::Tuple)
        Terrarium.update_state!(integrator, compute_tendencies = true)
        state = integrator.state
        for name in names
            state.prognostic[name] .+= 2 .* state.tendencies[name] .* Δt
        end
        return nothing
    end

    # Minimal two-variable model with constant unit tendencies for both prognostic variables.
    @kwdef struct TwoVarModel{NF, Grid <: AbstractLandGrid{NF}, TS <: Terrarium.AbstractTimeStepper} <: Terrarium.AbstractModel{NF, Grid}
        grid::Grid
        initializer = DefaultInitializer(eltype(grid))
        timestepper::TS = ForwardEuler(eltype(grid))
    end

    Terrarium.variables(::TwoVarModel) = (prognostic(:a, XY()), prognostic(:b, XY()))
    Terrarium.compute_auxiliary!(state, ::TwoVarModel) = nothing
    function Terrarium.compute_tendencies!(state, ::TwoVarModel)
        set!(state.tendencies.a, 1.0)
        set!(state.tendencies.b, 1.0)
        return nothing
    end
end

@testset "Single timestepper integrates all prognostic variables" begin
    grid = ColumnGrid(CPU(), Float64, UniformSpacing(N = 1))
    model = IMEXTestTypes.TwoVarModel(grid) # default: single ForwardEuler
    integrator = initialize(model)
    # a single timestepper stores its cache bare (ForwardEuler has none → empty)
    @test keys(integrator.state.timestepper_cache) == ()
    Δt = default_dt(integrator)
    timestep!(integrator)
    # both variables are stepped by forward Euler regardless of their class
    @test all(interior(integrator.state.a) .≈ Δt)
    @test all(interior(integrator.state.b) .≈ Δt)
end

@testset "Single Heun stores a bare cache" begin
    grid = ColumnGrid(CPU(), Float64, UniformSpacing(N = 1))
    model = IMEXTestTypes.TwoVarModel(grid; timestepper = Heun(Float64))
    integrator = initialize(model)
    # the Heun cache is stored directly (not keyed by class)
    @test keys(integrator.state.timestepper_cache) == (:prognostic, :tendencies)
    Δt = default_dt(integrator)
    timestep!(integrator)
    # constant tendency ⇒ Heun matches forward Euler: u += 1·Δt
    @test all(interior(integrator.state.a) .≈ Δt)
    @test all(interior(integrator.state.b) .≈ Δt)
end

@testset "IMEX routes variables to timesteppers by class" begin
    grid = ColumnGrid(CPU(), Float64, UniformSpacing(N = 1))
    timestepper = IMEX(ForwardEuler(Float64), IMEXTestTypes.MockImplicit(Float64))
    model = IMEXTestTypes.TwoVarModel(grid; timestepper)
    # override :b to the implicit class; :a keeps the :explicit default
    integrator = initialize(model; timestepper_classes = (; b = :implicit))
    @test Terrarium.prognostic_names(integrator.state, :explicit) == (:a,)
    @test Terrarium.prognostic_names(integrator.state, :implicit) == (:b,)
    # an IMEX cache keeps the two sub-stepper caches keyed by class
    @test keys(integrator.state.timestepper_cache) == (:explicit, :implicit)
    Δt = default_dt(integrator)
    timestep!(integrator)
    @test all(interior(integrator.state.a) .≈ Δt)       # forward Euler: a += 1·Δt
    @test all(interior(integrator.state.b) .≈ 2 * Δt)   # mock implicit: b += 2·Δt
end

@testset "IMEX with Heun explicit sub-stepper retrieves its cache" begin
    grid = ColumnGrid(CPU(), Float64, UniformSpacing(N = 1))
    timestepper = IMEX(Heun(Float64), IMEXTestTypes.MockImplicit(Float64))
    model = IMEXTestTypes.TwoVarModel(grid; timestepper)
    integrator = initialize(model; timestepper_classes = (; b = :implicit))
    # Heun's (prognostic, tendencies) cache lives in the explicit slot of the IMEX cache
    @test keys(integrator.state.timestepper_cache) == (:explicit, :implicit)
    @test keys(integrator.state.timestepper_cache.explicit) == (:prognostic, :tendencies)
    Δt = default_dt(integrator)
    timestep!(integrator)
    @test all(interior(integrator.state.a) .≈ Δt)       # Heun (constant tendency): a += 1·Δt
    @test all(interior(integrator.state.b) .≈ 2 * Δt)   # mock implicit: b += 2·Δt
end
