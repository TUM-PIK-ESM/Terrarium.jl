using Terrarium
using Test

import Oceananigans.Grids: RectilinearGrid, z_domain
import Terrarium.RingGrids
import Terrarium.RingGrids: FullHEALPixGrid, get_npoints

@testset "Vertical discretizations" begin
    # Uniform spacing
    @test Terrarium.get_spacing(UniformSpacing(Δz = 0.1, N = 1)) == [0.1]
    @test Terrarium.get_spacing(UniformSpacing(Δz = 0.1, N = 10)) == repeat([0.1], 10)

    # ExponentialSpacing
    @test Terrarium.get_spacing(ExponentialSpacing(Δz_min = 0.1, Δz_max = 1.0, N = 2)) == [0.1, 1.0]
    @test Terrarium.get_spacing(ExponentialSpacing(Δz_min = 0.1, Δz_max = 1.0, N = 3, sig = nothing)) ≈ exp2.(LinRange(log2(0.1), log2(1.0), 3))

    # PrescribedSpacing
    @test Terrarium.get_spacing(PrescribedSpacing(Δz = [0.1, 0.2, 0.3])) == [0.1, 0.2, 0.3]
end

@testset "ColumnGrid" begin
    # 2-column grid with 5 linearly spaced points
    num_columns = 2
    grid = ColumnGrid(UniformSpacing(Δz = 0.1, N = 5), num_columns)
    field_grid = get_field_grid(grid)
    @test isa(field_grid, RectilinearGrid)
    @test field_grid.Nx == num_columns
    @test field_grid.Ny == 1
    @test field_grid.Nz == 5
    @test z_domain(field_grid) == (-0.5, 0.0)
end

@testset "ColumnRingGrid" begin
    # test with 10-ring HEALPix
    ring_grid = FullHEALPixGrid(8)
    grid = ColumnRingGrid(UniformSpacing(Δz = 0.5, N = 10), ring_grid)
    field_grid = get_field_grid(grid)
    @test isa(field_grid, RectilinearGrid)
    @test field_grid.Nx == get_npoints(ring_grid)
    @test field_grid.Ny == 1
    @test field_grid.Nz == 10
    @test z_domain(field_grid) == (-5.0, 0.0)

    # Test RingGrids.Field to Oceananigans.Field conversion
    @testset "RingGrids to Oceananigans Field conversion" begin
        # Create a 2D RingGrids field with test data
        ring_grid = FullHEALPixGrid(8)
        grid = ColumnRingGrid(UniformSpacing(Δz = 0.5, N = 10), ring_grid)

        ring_field_2d = rand(ring_grid)

        # Convert to Oceananigans field
        oceananigans_field_2d = Terrarium.Field(ring_field_2d, grid)
        @test isa(oceananigans_field_2d, Terrarium.Field)
        @test size(oceananigans_field_2d) == (get_npoints(ring_grid), 1, 1)

        # Check that values were copied correctly
        ocean_data = Terrarium.interior(oceananigans_field_2d)
        @test all(ocean_data[:, 1, 1] .== ring_field_2d.data[grid.mask.data])

        # Test with 3D field (horizontal + vertical)
        ring_field_3d = rand(ring_grid, 10)

        ocean_field_3d = Terrarium.Field(ring_field_3d, grid)
        @test isa(ocean_field_3d, Terrarium.Field)
        @test size(ocean_field_3d) == (get_npoints(ring_grid), 1, 10)

        # Check that values were copied correctly for each vertical level
        ocean_data_3d = Terrarium.interior(ocean_field_3d)
        for k in 1:10
            @test all(ocean_data_3d[:, 1, k] .== ring_field_3d.data[grid.mask.data, k])
        end

        # Test with masked grid (some points inactive)

        # Create a new ring grid for this test
        mask = rand(Bool, ring_grid)
        # Create a masked ColumnRingGrid using the correct constructor
        masked_grid = ColumnRingGrid(UniformSpacing(Δz = 0.5, N = 10), mask)

        ring_field = rand(ring_grid)

        # Convert to Oceananigans field
        ocean_field_masked = Terrarium.Field(ring_field, masked_grid)
        @test isa(ocean_field_masked, Terrarium.Field)
        @test size(ocean_field_masked) == (sum(mask), 1, 1)

        # Verify only masked (active) points were copied
        ocean_data_masked = Terrarium.interior(ocean_field_masked)

        expected_values = ring_field.data[mask]
        @test all(ocean_data_masked[:, 1, 1] .== expected_values)

        # Test that conversion throws error for RingGrids fields with ndims >= 3
        ring_field_3d_plus = rand(ring_grid, 10, 5)  # horizontal × vertical × extra dimension
        @test_throws ErrorException Terrarium.Field(ring_field_3d_plus, grid)
    end
end
