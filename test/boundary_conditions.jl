using Terrarium
using Terrarium: prognostic, varname, XY, XYZ
using Oceananigans: CenterField
using Test

@testset "Boundary conditions" begin
    # test default boundary conditions do nothing
    var = prognostic(:test, XYZ())
    grid = ColumnGrid(ExponentialSpacing())
    default_bcs = DefaultBoundaryConditions()
    @test isnothing(Terrarium.get_field_boundary_conditions(default_bcs, grid, var))
    @test variables(default_bcs) == ()
    
    # test prescribed flux
    flux = PrescribedFlux(:flux, 1.0, XY())
    flux_vars = variables(flux)
    @test length(flux_vars) == 1
    @test Terrarium.varname(flux_vars[1]) == :flux
    @test Terrarium.vardims(flux_vars[1]) == XY()
    @test flux.value == 1.0
    fluxfield = Terrarium.create_field(flux_vars[1], DefaultInitializer(), DefaultBoundaryConditions(), grid)
    # check that compute_auxiliary! works
    compute_auxiliary!((flux=fluxfield,), nothing, flux)
    @test all(fluxfield .== 1.0)

    # VerticalBoundaryConditions
    bcs = VerticalBoundaryConditions(top=PrescribedFlux(:fluxtop, 1.0), bottom=PrescribedFlux(:fluxbot, -1.0))
    vars = variables(bcs)
    @test length(vars) == 2
    @test map(varname, vars) == (:fluxtop, :fluxbot)
    fluxtop = Terrarium.create_field(vars[1], DefaultInitializer(), DefaultBoundaryConditions(), grid)
    fluxbot = Terrarium.create_field(vars[2], DefaultInitializer(), DefaultBoundaryConditions(), grid)
    compute_auxiliary!((; fluxtop, fluxbot), nothing, bcs)
    @test all(fluxtop .== 1.0)
    @test all(fluxbot .== -1.0)
    # check that directly specified Field boundary conditions are correctly assigned
    test_bcs = VerticalBoundaryConditions(top=(test=ValueBoundaryCondition(1.0),), bottom=DefaultBoundaryConditions())
    testfield = Terrarium.create_field(var, DefaultInitializer(), test_bcs, grid)
    @test testfield.boundary_conditions.top == ValueBoundaryCondition(1.0)
    test_bcs = VerticalBoundaryConditions(top=DefaultBoundaryConditions(), bottom=(test=ValueBoundaryCondition(1.0),))
    testfield = Terrarium.create_field(var, DefaultInitializer(), test_bcs, grid)
    @test testfield.boundary_conditions.bottom == ValueBoundaryCondition(1.0)
    test_bcs = VerticalBoundaryConditions(top=(test=ValueBoundaryCondition(1.0),), bottom=(test=ValueBoundaryCondition(1.0),))
    testfield = Terrarium.create_field(var, DefaultInitializer(), test_bcs, grid)
    @test testfield.boundary_conditions.top == ValueBoundaryCondition(1.0)
    @test testfield.boundary_conditions.bottom == ValueBoundaryCondition(1.0)
end
