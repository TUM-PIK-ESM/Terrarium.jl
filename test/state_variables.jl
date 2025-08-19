using Terrarium
using Test

import Terrarium: VarDims, XY, XYZ, prognostic, auxiliary, namespace
import Oceananigans: Field, Center

DEFAULT_NF = Float32

@testset "State variables for mock type" begin

    @kwdef struct SubModel{NF} <: Terrarium.AbstractModel{NF}
        grid
        initializer = DefaultInitializer()
        boundary_conditions = DefaultBoundaryConditions()
        time_stepping::Terrarium.AbstractTimeStepper{NF} = Terrarium.ForwardEuler{DEFAULT_NF}()
    end

    Terrarium.variables(model::SubModel) = (
        # duplicate naming allowed in new namesapce
        auxiliary(:auxvar2D, XY()),
    )

    @kwdef struct TestModel{NF} <: Terrarium.AbstractModel{NF}
        grid
        submodel = SubModel(; grid)
        initializer = DefaultInitializer()
        boundary_conditions = DefaultBoundaryConditions()
        time_stepping::Terrarium.AbstractTimeStepper{NF} = Terrarium.ForwardEuler{DEFAULT_NF}()
    end

    struct TestClosure <: Terrarium.AbstractClosureRelation end

    Terrarium.getvar(::TestClosure, dims::VarDims) = auxiliary(:closurevar, dims)

    Terrarium.variables(model::TestModel) = (
        prognostic(:progvar3D, XYZ()),
        prognostic(:cprogvar3D, XYZ(), TestClosure()),
        prognostic(:progvar2D, XY()),
        auxiliary(:auxvar3D, XYZ()),
        auxiliary(:auxvar2D, XY()),
        namespace(:submodel),
    )

    grid = ColumnGrid(ExponentialSpacing(N=10))
    model = TestModel{DEFAULT_NF}(; grid)
    sim = initialize(model)
    # Check that all prognostic variables are defined correctly
    @test hasproperty(sim.state, :progvar3D) && isa(sim.state.prognostic.progvar3D, Field{Center,Center,Center})
    @test hasproperty(sim.state, :progvar2D) && isa(sim.state.prognostic.progvar2D, Field{Center,Center,Nothing})
    @test hasproperty(sim.state, :cprogvar3D) && isa(sim.state.prognostic.cprogvar3D, Field{Center,Center,Center})
    # Check that all tendencies are defined correctly
    @test hasproperty(sim.state, :progvar3D_tendency) && isa(sim.state.progvar3D_tendency, Field{Center,Center,Center})
    @test hasproperty(sim.state, :progvar2D_tendency) && isa(sim.state.progvar2D, Field{Center,Center,Nothing})
    # Check that tendency for prognostic variable with closure relation is defined only on the clousre variable
    @test hasproperty(sim.state, :closurevar_tendency) && isa(sim.state.progvar2D, Field{Center,Center,Nothing})
    @test !hasproperty(sim.state, :cprogvar3D_tendency)
    # Check that all auxiliary variables are defined correctly
    @test hasproperty(sim.state, :auxvar3D) && isa(sim.state.auxvar3D, Field{Center,Center,Center})
    @test hasproperty(sim.state, :auxvar2D) && isa(sim.state.auxvar2D, Field{Center,Center,Nothing})
    # Check submodel namespace
    @test hasproperty(sim.state, :submodel) && isa(sim.state.submodel, StateVariables)
    @test hasproperty(sim.state.submodel, :auxvar2D) && isa(sim.state.submodel.auxvar2D, Field{Center,Center,Nothing})

end