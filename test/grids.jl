using Terrarium
using Test

import Oceananigans.Grids: Grids, RectilinearGrid
import Terrarium.RingGrids: FullHEALPixGrid, get_npoints

@testset "Vertical discretizations" begin
    # Uniform spacing
    @test Terrarium.get_spacing(UniformSpacing(Δz=0.1, N=1)) == [0.1]
    @test Terrarium.get_spacing(UniformSpacing(Δz=0.1, N=10)) == repeat([0.1], 10)

    # ExponentialSpacing
    @test Terrarium.get_spacing(ExponentialSpacing(Δz_min=0.1, Δz_max=1.0, N=2)) == [0.1,1.0]
    @test Terrarium.get_spacing(ExponentialSpacing(Δz_min=0.1, Δz_max=1.0, N=3, sig=nothing)) ≈ exp2.(LinRange(log2(0.1), log2(1.0), 3))

    # PrescribedSpacing
    @test Terrarium.get_spacing(PrescribedSpacing(Δz=[0.1,0.2,0.3])) == [0.1,0.2,0.3]
end

@testset "ColumnGrid" begin
    # 2-column grid with 5 linearly spaced points
    num_columns = 2
    grid = ColumnGrid(UniformSpacing(Δz=0.1, N=5), num_columns)
    field_grid = get_field_grid(grid)
    @test isa(field_grid, RectilinearGrid)
    @test field_grid.Nx == num_columns
    @test field_grid.Ny == 1
    @test field_grid.Nz == 5
    @test Grids.z_domain(field_grid) == (-0.5, 0.0)
end

@testset "GlobalRingGrid" begin
    # test with 10-ring HEALPix 
    ring_grid = FullHEALPixGrid(8)
    grid = GlobalRingGrid(UniformSpacing(Δz=0.5, N=10), ring_grid)
    field_grid = get_field_grid(grid)
    @test isa(field_grid, RectilinearGrid)
    @test field_grid.Nx == get_npoints(ring_grid)
    @test field_grid.Ny == 1
    @test field_grid.Nz == 10
    @test Grids.z_domain(field_grid) == (-5.0, 0.0)
end
