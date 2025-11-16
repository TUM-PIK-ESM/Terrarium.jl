using Terrarium
using Terrarium: prognostic, varname, vardims, XY, XYZ
using Oceananigans: CenterField
using Test

@testset "Boundary conditions" begin
    var = prognostic(:test, XYZ())
    grid = ColumnGrid(ExponentialSpacing())
    # test prescribed flux
    flux = PrescribedFlux(:flux, 1.0, XY())
    flux_vars = variables(flux)
    @test length(flux_vars) == 1
    @test Terrarium.varname(flux_vars[1]) == :flux
    @test Terrarium.vardims(flux_vars[1]) == XY()
    @test flux.value == 1.0
    fluxfield = Field(grid, vardims(flux_vars[1]), DefaultBoundaryConditions())
    # check that compute_auxiliary! works
    compute_auxiliary!((flux=fluxfield,), nothing, flux)
    @test all(fluxfield .== 1.0)

    # ColumnBoundaryConditions
    bcs = ColumnBoundaryConditions(top=PrescribedFlux(:fluxtop, 1.0), bottom=PrescribedFlux(:fluxbot, -1.0))
    vars = variables(bcs)
    @test length(vars) == 2
    @test map(varname, vars) == (:fluxtop, :fluxbot)
    fluxtop = Field(grid, vardims(vars[1]), DefaultBoundaryConditions())
    fluxbot = Field(grid, vardims(vars[2]), DefaultBoundaryConditions())
    compute_auxiliary!((; fluxtop, fluxbot), nothing, bcs)
    @test all(fluxtop .== 1.0)
    @test all(fluxbot .== -1.0)
    # check that directly specified Field boundary conditions are correctly assigned
    test_bcs = ColumnBoundaryConditions(top=(test=ValueBoundaryCondition(1.0),), bottom=DefaultBoundaryConditions())
    testfield = Field(grid, vardims(var), test_bcs)
    @test testfield.boundary_conditions.top == ValueBoundaryCondition(1.0)
    test_bcs = ColumnBoundaryConditions(top=DefaultBoundaryConditions(), bottom=(test=ValueBoundaryCondition(1.0),))
    testfield = Field(grid, vardims(var), test_bcs)
    @test testfield.boundary_conditions.bottom == ValueBoundaryCondition(1.0)
    test_bcs = ColumnBoundaryConditions(top=(test=ValueBoundaryCondition(1.0),), bottom=(test=ValueBoundaryCondition(1.0),))
    testfield = Field(grid, vardims(var), test_bcs)
    @test testfield.boundary_conditions.top == ValueBoundaryCondition(1.0)
    @test testfield.boundary_conditions.bottom == ValueBoundaryCondition(1.0)
end
