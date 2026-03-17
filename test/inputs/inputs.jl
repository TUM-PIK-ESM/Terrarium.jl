using Terrarium
using Terrarium: FieldInputSource, FieldTimeSeriesInputSource, Variables, initialize!, interior
using Test
using Unitful

@testset "Input sources" begin
    # Field input source
    grid = ColumnGrid(ExponentialSpacing())
    X1 = Field(grid, XY())
    field_input = InputSource(grid, (; X1))
    @test isa(field_input, FieldInputSource)
    ## check that dimensions were inferred correctly
    @test field_input.dims == XY()
    @test variables(field_input) == (Terrarium.input(:X1, XY()),)
    X2 = Field(grid, XY())
    field_input = InputSource(grid, (; X1, X2))
    @test isa(field_input, FieldInputSource)
    @test field_input.dims == XY()
    @test variables(field_input) == (Terrarium.input(:X1, XY()), Terrarium.input(:X2, XY()))
    ## check state variables are allocated
    state = initialize(Variables(field_input), grid)
    @test hasproperty(state.inputs, :X1)
    @test hasproperty(state.inputs, :X2)

    # FieldTimeSeries input source
    grid = ColumnGrid(ExponentialSpacing())
    ts = 0.0:1.0:10.0
    S1 = FieldTimeSeries(grid, XY(), ts)
    fts_input = InputSource(; S1)
    @test isa(fts_input, FieldTimeSeriesInputSource)
    ## check that dimensions were inferred correctly
    @test fts_input.dims == XY()
    @test fts_input.fts == (; S1)
    S2 = FieldTimeSeries(grid, XY(), ts)
    fts_input = InputSource(; S1, S2)
    @test isa(fts_input, FieldTimeSeriesInputSource)
    @test fts_input.dims == XY()
    @test fts_input.fts == (; S1, S2)
    S3 = FieldTimeSeries(grid, XYZ(), ts)
    ## check that trying to create a FieldTimeSeriesInputSource
    # with inconsistent dimensions fails
    @test_throws AssertionError InputSource(; S1, S2, S3)
    # populate S1 with random data
    S1.data .= randn(size(S1))
    ## check update_inputs!
    fields = (S1 = Field(grid, XY()), S2 = Field(grid, XY()))
    clock = Clock(time = 0)
    update_inputs!(fields, fts_input, clock)
    @test all(fields.S1 .== S1[1])
    @test all(fields.S2 .== S2[1])
    # advance clock and check that inputs are updated
    Terrarium.tick!(clock, 1.0)
    update_inputs!(fields, fts_input, clock)
    @test all(fields.S1 .== S1[2])
    @test all(fields.S2 .== S2[2])
end

@testset "RingGrids.Field InputSource convenience" begin
    import RingGrids

    # Create a ColumnRingGrid
    ring_grid = RingGrids.FullHEALPixGrid(8)
    npoints = RingGrids.get_npoints(ring_grid)
    grid = ColumnRingGrid(UniformSpacing(Δz = 0.5, N = 5), ring_grid)

    # Create RingGrids fields with known data
    ring_field1 = RingGrids.Field(ring_grid)
    ring_field1.data .= collect(1.0:Float32(npoints))

    ring_field2 = RingGrids.Field(ring_grid)
    ring_field2.data .= collect(Float32(npoints):-1.0:1.0)

    # Test single field
    source = InputSource(grid, (; temperature = ring_field1))
    @test isa(source, FieldInputSource)
    @test source.dims == XY()
    @test variables(source) == (Terrarium.input(:temperature, XY()),)

    # Test multiple fields
    source = InputSource(grid, (; temp = ring_field1, pressure = ring_field2))
    @test isa(source, FieldInputSource)
    @test length(variables(source)) == 2

    # Test that data is correctly converted and stored
    state = initialize(Variables(source), grid)
    @test hasproperty(state.inputs, :temp)
    @test hasproperty(state.inputs, :pressure)

    # Initialize and verify data was copied correctly
    initialize!(state.inputs, source)
    expected1 = ring_field1.data[Array(grid.mask)]
    expected2 = ring_field2.data[Array(grid.mask)]
    @test all(interior(state.inputs.temp)[:, 1, 1] .≈ expected1)
    @test all(interior(state.inputs.pressure)[:, 1, 1] .≈ expected2)

    # Test 3D RingGrids field (with vertical dimension)
    ring_field_3d = rand(ring_grid, 5)  # 5 vertical levels
    source_3d = InputSource(grid, (; field3d = ring_field_3d))
    @test isa(source_3d, FieldInputSource)
    @test source_3d.dims == XYZ()
end
