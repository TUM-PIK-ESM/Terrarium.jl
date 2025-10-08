using Terrarium
using Test

struct TestClosure
    varname::Symbol
end

Terrarium.getvar(closure::TestClosure) = Terrarium.auxiliary(closure.varname, XYZ())

@testset "Forward Euler" begin
    Δt = 10.0
    euler = ForwardEuler(; Δt)
    @test !is_adaptive(euler)
    @test default_dt(euler) == Δt
    # set up grid and fields
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N=10))
    # here we mock the structure of a `StateVariables` object
    # for a model with prognostic variables at the top level and
    # in a nested namespace.
    state = (
        prognostic = (x = Field(grid, XYZ()), y = Field(grid, XYZ())),
        auxiliary = (z = Field(grid, XYZ()),),
        tendencies = (
            x = Field(grid, XYZ()),
            z = Field(grid, XYZ()),
        ),
        closures = (y = TestClosure(:z),),
        namespaces = (
            inner = (
                prognostic = (x = Field(grid, XYZ()),),
                auxiliary = (;),
                tendencies = (
                    x = Field(grid, XYZ()),
                ),
                namespaces = (;),
                closures = (;),
            ),
        )
    )
    dxdt = 0.1
    dzdt = 0.2
    set!(state.tendencies.x, dxdt)
    set!(state.tendencies.z, dzdt)
    set!(state.namespaces.inner.tendencies.x, dxdt*2)
    Terrarium.explicit_step!(state, grid, ForwardEuler(; Δt), Δt)
    @test all(state.prognostic.x .≈ Δt*dxdt)
    @test all(state.auxiliary.z .≈ Δt*dzdt)
    @test all(state.namespaces.inner.prognostic.x .≈ Δt*dxdt*2)
    # check that y was not changed (inverse closure not evaluated)
    @test all(iszero.(state.prognostic.y))
end