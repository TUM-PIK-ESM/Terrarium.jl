using Terrarium
using Terrarium: Variables, Namespace, initialize!, interior, varname, matches_scope, namespaced_variables
using Test

@testset "inputpath" begin
    @test inputpath(:x) == (:x,)
    @test inputpath(:ns => :x) == (:ns, :x)
    @test inputpath(:a => :b => :c) == (:a, :b, :c)
    @test inputpath((:ns, :x)) == (:ns, :x)
end

@testset "Namespaced input sources" begin
    grid = ColumnGrid(ExponentialSpacing())

    # source with a namespaced name
    X1 = Field(grid, XY())
    src = InputSource(grid, X1; name = :ns1 => :x)
    @test varname(src) == :x
    @test inputpath(src) == (:ns1, :x)
    @test matches_scope(src, (:ns1,))
    @test !matches_scope(src, ())
    @test !matches_scope(src, (:ns2,))

    # variables(src) should declare the input inside the namespace
    vars = variables(src)
    @test length(vars) == 1 && isa(vars[1], Namespace)
    @test varname(vars[1]) == :ns1
    @test vars[1].vars.inputs.x == Terrarium.input(:x, XY())

    # initialize state with a root input :x and two namespaces both declaring :x
    all_vars = Variables(
        Terrarium.input(:x, XY()),
        Terrarium.namespace(:ns1, (Terrarium.input(:x, XY()),)),
        Terrarium.namespace(:ns2, (Terrarium.input(:x, XY()),)),
    )
    state = initialize(all_vars, grid)
    X1 .= 3.0f0
    initialize!(state, InputSources(src))
    ## only ns1.x should be set by the source
    @test all(interior(state.namespaces.ns1.x) .== 3.0f0)
    @test all(interior(state.namespaces.ns2.x) .== 0.0f0)
    @test all(interior(state.inputs.x) .== 0.0f0)

    # a root-level source should not leak into namespaces
    X2 = Field(grid, XY())
    X2 .= 5.0f0
    root_src = InputSource(grid, X2; name = :x)
    initialize!(state, InputSources(root_src))
    @test all(interior(state.inputs.x) .== 5.0f0)
    @test all(interior(state.namespaces.ns1.x) .== 3.0f0)
    @test all(interior(state.namespaces.ns2.x) .== 0.0f0)

    # static field sources do not change values in update_inputs!
    X1 .= 7.0f0
    update_inputs!(state, InputSources(src))
    @test all(interior(state.namespaces.ns1.x) .== 3.0f0)

    # time-varying sources follow the same scoping rules in update_inputs!
    ts = 0.0:1.0:10.0
    S = FieldTimeSeries(grid, XY(), ts)
    for (i, t) in enumerate(ts)
        S[i] .= Float32(t)
    end
    fts_src = InputSource(S; name = :ns2 => :x)
    @test inputpath(fts_src) == (:ns2, :x)
    update_inputs!(state, InputSources(fts_src))
    @test all(interior(state.namespaces.ns2.x) .== 0.0f0)
    Terrarium.tick!(state.clock, 1.0)
    update_inputs!(state, InputSources(fts_src))
    @test all(interior(state.namespaces.ns2.x) .== 1.0f0)
    ## other fields named x are untouched
    @test all(interior(state.namespaces.ns1.x) .== 3.0f0)
    @test all(interior(state.inputs.x) .== 5.0f0)
end

@testset "Namespace merging in Variables" begin
    # same-named namespaces should merge their variables
    vars = Variables(
        Terrarium.namespace(:ns, (Terrarium.input(:a, XY()),)),
        Terrarium.namespace(:ns, (Terrarium.input(:b, XY()),)),
    )
    @test keys(vars.namespaces) == (:ns,)
    @test keys(vars.namespaces.ns.vars.inputs) == (:a, :b)

    # nested namespaces merge recursively
    nested = Variables(
        Terrarium.namespace(:outer, (Terrarium.namespace(:inner, (Terrarium.input(:a, XY()),)),)),
        Terrarium.namespace(:outer, (Terrarium.namespace(:inner, (Terrarium.input(:b, XY()),)),)),
    )
    @test keys(nested.namespaces) == (:outer,)
    @test keys(nested.namespaces.outer.vars.namespaces.inner.vars.inputs) == (:a, :b)
end

@testset "Namespaced inputs with SoilModel" begin
    NF = Float32
    grid = ColumnGrid(CPU(), NF, ExponentialSpacing())
    soil = Terrarium.SoilEnergyWaterCarbon(NF; strat = OARSoilStratigraphy(NF))
    model = SoilModel(grid; soil)

    # prescribe texture for the organic horizon only; note that the fractions
    # must sum to unity in each grid cell
    sand = Field(grid, XY())
    silt = Field(grid, XY())
    clay = Field(grid, XY())
    sand .= 0.25f0
    silt .= 0.5f0
    clay .= 0.25f0
    inputs = InputSources(
        InputSource(grid, sand; name = :organic => :sand),
        InputSource(grid, silt; name = :organic => :silt),
        InputSource(grid, clay; name = :organic => :clay),
    )
    integrator = initialize(model, ForwardEuler(NF); inputs)
    state = integrator.state
    @test Terrarium.namespace_names(state) == (:organic, :surface, :bedrock)
    ## inputs land only in the organic horizon namespace
    @test all(interior(state.namespaces.organic.sand) .== 0.25f0)
    @test all(interior(state.namespaces.organic.silt) .== 0.5f0)
    @test all(interior(state.namespaces.organic.clay) .== 0.25f0)
    ## sibling horizons keep their defaults (sand defaults to one)
    @test all(interior(state.namespaces.surface.sand) .== 1.0f0)
    @test all(interior(state.namespaces.bedrock.sand) .== 1.0f0)

    # simulation runs and inputs are preserved across timesteps
    run!(integrator; steps = 5, Δt = 60.0)
    @test all(interior(state.namespaces.organic.sand) .== 0.25f0)
end
