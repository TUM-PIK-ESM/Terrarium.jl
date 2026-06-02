using Terrarium
using Test

using Terrarium: AbstractLandGrid, AbstractImplicitTimestepper, EmptyCache, prognostic, XY

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

    # MockImplicit keeps no working state, so it doesn't fetch a cache.
    function Terrarium.timestep!(integrator, ts::MockImplicit, Δt, names::Tuple)
        Terrarium.update_state!(integrator, compute_tendencies = true)
        state = integrator.state
        for name in names
            state.prognostic[name] .+= 2 .* state.tendencies[name] .* Δt
        end
        return nothing
    end

    # Minimal two-variable model with constant unit tendencies. `:a` defaults to the :explicit class and
    # `:b` defaults to :implicit, so a single model exercises both the declared default and IMEX overrides.
    @kwdef struct TwoVarModel{NF, Grid <: AbstractLandGrid{NF}, TS <: Terrarium.AbstractTimeStepper} <: Terrarium.AbstractModel{NF, Grid}
        grid::Grid
        initializer = DefaultInitializer(eltype(grid))
        timestepper::TS = ForwardEuler(eltype(grid))
    end

    Terrarium.variables(::TwoVarModel) = (prognostic(:a, XY()), prognostic(:b, XY(); timestepper = :implicit))
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
    # a single stateless timestepper gets an EmptyCache
    @test integrator.state.timestepper_cache isa EmptyCache
    Δt = default_dt(integrator)
    timestep!(integrator)
    # both variables are stepped by forward Euler regardless of their declared class
    @test all(interior(integrator.state.a) .≈ Δt)
    @test all(interior(integrator.state.b) .≈ Δt)
end

@testset "Single Heun gets a HeunCache" begin
    grid = ColumnGrid(CPU(), Float64, UniformSpacing(N = 1))
    model = IMEXTestTypes.TwoVarModel(grid; timestepper = Heun(Float64))
    integrator = initialize(model)
    @test integrator.state.timestepper_cache isa Terrarium.HeunCache
    Δt = default_dt(integrator)
    timestep!(integrator)
    # constant tendency ⇒ Heun matches forward Euler: u += 1·Δt
    @test all(interior(integrator.state.a) .≈ Δt)
    @test all(interior(integrator.state.b) .≈ Δt)
end

@testset "IMEX routes by declared default class" begin
    grid = ColumnGrid(CPU(), Float64, UniformSpacing(N = 1))
    timestepper = IMEX(ForwardEuler(Float64), IMEXTestTypes.MockImplicit(Float64))
    model = IMEXTestTypes.TwoVarModel(grid; timestepper)
    integrator = initialize(model)
    cache = integrator.state.timestepper_cache
    # the IMEXCache holds the resolved per-variable classes (in prognostic order) as a type parameter
    @test cache isa Terrarium.IMEXCache
    @test Terrarium.timestepper_classes(cache) == (:explicit, :implicit)   # (:a, :b)
    @test cache.explicit isa EmptyCache && cache.implicit isa EmptyCache
    Δt = default_dt(integrator)
    timestep!(integrator)
    @test all(interior(integrator.state.a) .≈ Δt)       # :a default :explicit → forward Euler: a += 1·Δt
    @test all(interior(integrator.state.b) .≈ 2 * Δt)   # :b default :implicit → mock: b += 2·Δt
end

@testset "IMEX timestepper_classes overrides the default; Heun explicit cache" begin
    grid = ColumnGrid(CPU(), Float64, UniformSpacing(N = 1))
    # flip both variables relative to their declared defaults
    timestepper = IMEX(Heun(Float64), IMEXTestTypes.MockImplicit(Float64); timestepper_classes = (; a = :implicit, b = :explicit))
    model = IMEXTestTypes.TwoVarModel(grid; timestepper)
    integrator = initialize(model)
    cache = integrator.state.timestepper_cache
    @test Terrarium.timestepper_classes(cache) == (:implicit, :explicit)   # (:a, :b)
    @test cache.explicit isa Terrarium.HeunCache   # Heun's cache lives in the explicit slot
    Δt = default_dt(integrator)
    timestep!(integrator)
    @test all(interior(integrator.state.a) .≈ 2 * Δt)   # overridden to :implicit → mock: a += 2·Δt
    @test all(interior(integrator.state.b) .≈ Δt)       # overridden to :explicit → Heun: b += 1·Δt
end
