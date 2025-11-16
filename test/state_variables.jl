module StateVariablesTests

using Terrarium
using Test

import Terrarium: AbstractLandGrid, VarDims, XY, XYZ, prognostic, auxiliary, input, namespace

DEFAULT_NF = Float32

@kwdef struct SubModel{NF, Grid<:AbstractLandGrid{NF}} <: Terrarium.AbstractModel{NF, Grid}
    grid::Grid
    initializer = DefaultInitializer()
end

Terrarium.variables(model::SubModel) = (
    # duplicate naming allowed in new namesapce
    auxiliary(:auxvar2D, XY()),
    # inputs are handled "globally" (i.e. all inputs with a given name refer to the same field)
    input(:forcing, XY()),
)

@kwdef struct TestModel{NF, Grid<:AbstractLandGrid{NF}} <: Terrarium.AbstractModel{NF, Grid}
    grid::Grid
    submodel = SubModel(; grid)
    initializer = DefaultInitializer()
end

struct TestClosure <: Terrarium.AbstractClosureRelation end

Terrarium.closurevar(::TestClosure) = auxiliary(:closurevar, XYZ())

Terrarium.variables(model::TestModel) = (
    prognostic(:progvar3D, XYZ()),
    prognostic(:cprogvar3D, XYZ(), closure=TestClosure()),
    prognostic(:progvar2D, XY()),
    auxiliary(:auxvar3D, XYZ()),
    auxiliary(:auxvar2D, XY()),
    input(:forcing, XY()),
    namespace(:submodel, variables(model.submodel)),
)

@testset "State variable initialization" begin
    grid = ColumnGrid(CPU(), DEFAULT_NF, ExponentialSpacing(N=10))
    model = TestModel(; grid)
    state = initialize(model)
    # Check that all prognostic variables are defined correctly
    @test hasproperty(state, :progvar3D) && isa(state.prognostic.progvar3D, Field{Center,Center,Center})
    @test hasproperty(state, :progvar2D) && isa(state.prognostic.progvar2D, Field{Center,Center,Nothing})
    @test hasproperty(state, :cprogvar3D) && isa(state.prognostic.cprogvar3D, Field{Center,Center,Center})
    # Check that all tendencies are defined correctly
    @test hasproperty(state.tendencies, :progvar3D) && isa(state.tendencies.progvar3D, Field{Center,Center,Center})
    @test hasproperty(state.tendencies, :progvar2D) && isa(state.tendencies.progvar2D, Field{Center,Center,Nothing})
    # Check that tendency for prognostic variable with closure relation is defined only on the clousre variable
    @test hasproperty(state.tendencies, :closurevar) && isa(state.progvar2D, Field{Center,Center,Nothing})
    @test !hasproperty(state.tendencies, :cprogvar3D)
    # Check that all auxiliary variables are defined correctly
    @test hasproperty(state, :auxvar3D) && isa(state.auxvar3D, Field{Center,Center,Center})
    @test hasproperty(state, :auxvar2D) && isa(state.auxvar2D, Field{Center,Center,Nothing})
    # Check that all input variables are defined correctly
    @test hasproperty(state, :forcing) && isa(state.auxvar2D, Field{Center,Center,Nothing})
    # Check submodel namespace
    @test hasproperty(state, :submodel) && isa(state.submodel, StateVariables)
    @test hasproperty(state.submodel, :auxvar2D) && isa(state.submodel.auxvar2D, Field{Center,Center,Nothing})
end

end
