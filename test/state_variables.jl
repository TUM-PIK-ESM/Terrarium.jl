module StateVariablesTests

using Terrarium
using Test

import Terrarium: AbstractLandGrid, VarDims, XY, XYZ, prognostic, auxiliary, input, namespace

DEFAULT_NF = Float32

@kwdef struct SubModel{NF, Grid<:AbstractLandGrid{NF}, TS} <: Terrarium.AbstractModel{NF, Grid, TS}
    grid::Grid
    initializer = DefaultInitializer()
    boundary_conditions = DefaultBoundaryConditions()
    time_stepping::TS = Terrarium.ForwardEuler{DEFAULT_NF}()
end

Terrarium.variables(model::SubModel) = (
    # duplicate naming allowed in new namesapce
    auxiliary(:auxvar2D, XY()),
    # inputs are handled "globally" (i.e. all inputs with a given name refer to the same field)
    input(:forcing, XY()),
)

@kwdef struct TestModel{NF, Grid<:AbstractLandGrid{NF}, TS} <: Terrarium.AbstractModel{NF, Grid, TS}
    grid::Grid
    submodel = SubModel(; grid)
    initializer = DefaultInitializer()
    boundary_conditions = DefaultBoundaryConditions()
    time_stepping::TS = Terrarium.ForwardEuler{DEFAULT_NF}()
end

struct TestClosure <: Terrarium.AbstractClosureRelation end

Terrarium.getvar(::TestClosure) = auxiliary(:closurevar, XYZ())

Terrarium.variables(model::TestModel) = (
    prognostic(:progvar3D, XYZ()),
    prognostic(:cprogvar3D, XYZ(), TestClosure()),
    prognostic(:progvar2D, XY()),
    auxiliary(:auxvar3D, XYZ()),
    auxiliary(:auxvar2D, XY()),
    input(:forcing, XY()),
    namespace(:submodel, variables(model.submodel)),
)

@testset "State variable initialization" begin

    grid = ColumnGrid(CPU(), DEFAULT_NF, ExponentialSpacing(N=10))
    model = TestModel(; grid)
    sim = initialize(model)
    # Check that all prognostic variables are defined correctly
    @test hasproperty(sim.state, :progvar3D) && isa(sim.state.prognostic.progvar3D, Field{Center,Center,Center})
    @test hasproperty(sim.state, :progvar2D) && isa(sim.state.prognostic.progvar2D, Field{Center,Center,Nothing})
    @test hasproperty(sim.state, :cprogvar3D) && isa(sim.state.prognostic.cprogvar3D, Field{Center,Center,Center})
    # Check that all tendencies are defined correctly
    @test hasproperty(sim.state.tendencies, :progvar3D) && isa(sim.state.tendencies.progvar3D, Field{Center,Center,Center})
    @test hasproperty(sim.state.tendencies, :progvar2D) && isa(sim.state.tendencies.progvar2D, Field{Center,Center,Nothing})
    # Check that tendency for prognostic variable with closure relation is defined only on the clousre variable
    @test hasproperty(sim.state.tendencies, :closurevar) && isa(sim.state.progvar2D, Field{Center,Center,Nothing})
    @test !hasproperty(sim.state.tendencies, :cprogvar3D)
    # Check that all auxiliary variables are defined correctly
    @test hasproperty(sim.state, :auxvar3D) && isa(sim.state.auxvar3D, Field{Center,Center,Center})
    @test hasproperty(sim.state, :auxvar2D) && isa(sim.state.auxvar2D, Field{Center,Center,Nothing})
    # Check that all input variables are defined correctly
    @test hasproperty(sim.state, :forcing) && isa(sim.state.auxvar2D, Field{Center,Center,Nothing})
    # Check submodel namespace
    @test hasproperty(sim.state, :submodel) && isa(sim.state.submodel, StateVariables)
    @test hasproperty(sim.state.submodel, :auxvar2D) && isa(sim.state.submodel.auxvar2D, Field{Center,Center,Nothing})
    # Check that input variable is identical across namespaces
    @test sim.state.forcing === sim.state.submodel.forcing
end

end
