using Terrarium
using Terrarium: ConstantAlbedo, PrescribedAlbedo
using Test

using Oceananigans

@testset "Constant albedo" begin
    albedo = ConstantAlbedo(albedo = 0.4, emissivity = 0.8)
    @test Terrarium.albedo(1, 1, nothing, albedo) == 0.4
    @test Terrarium.emissivity(1, 1, nothing, albedo) == 0.8
end

@testset "PrescribedAlbedo" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N = 10))
    albedo = PrescribedAlbedo()
    state = (
        inputs = (albedo = set!(Field(grid, XY()), 0.4), emissivity = set!(Field(grid, XY()), 0.8))
    )
    @test Terrarium.albedo(1, 1, state, albedo) == 0.4
    @test Terrarium.emissivity(1, 1, state, albedo) == 0.8
end
