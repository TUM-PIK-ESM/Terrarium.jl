using Terrarium
using Terrarium: FieldInputSource, FieldTimeSeriesInputSource
using Test
using Unitful

@testset "Input sources" begin
    # Field input source
    grid = ColumnGrid(ExponentialSpacing())
    X1 = Field(grid, XY())
    field_input = InputSource(; X1)
    @test isa(field_input, FieldInputSource)
    ## check that dimensions were inferred correctly
    @test field_input.dims == XY()
    @test field_input.fields == (; X1)
    X2 = Field(grid, XY())
    field_input = InputSource(; X1, X2)
    @test isa(field_input, FieldInputSource)
    @test field_input.dims == XY()
    @test field_input.fields == (; X1, X2)
    X3 = Field(grid, XYZ())
    ## check that trying to create a FieldInputSource with inconsistent types or dimensions fails
    @test_throws AssertionError InputSource(; X1, X2, X3)
    @test_throws AssertionError InputSource(; X1, X2, Y = zeros(size(X1)))
    ## check update_inputs!
    set!(X1, 1.0)
    set!(X2, 2.0)
    fields = (X1 = Field(grid, XY()), X2 = Field(grid, XY()))
    clock = Clock(time = 0)
    update_inputs!(fields, field_input, clock)
    @test all(fields.X1 .== X1)
    @test all(fields.X2 .== X2)

    # FieldTimeSeries input source
    grid = ColumnGrid(ExponentialSpacing())
    ts = 0.0:1.0:10.0
    S1 = FieldTimeSeries(grid, XY(), ts)
    fts_input = InputSource(; S1)
    @test isa(field_input, FieldTimeSeriesInputSource)
    ## check that dimensions were inferred correctly
    @test fts_input.dims == XY()
    @test fts_input.fts == (; S1)
    S2 = FieldTimeSeries(grid, XY(), ts)
    fts_input = InputSource(; S1, S2)
    @test isa(field_input, FieldTimeSeriesInputSource)
    @test fts_input.dims == XY()
    @test fts_input.fts == (; S1, S2)
    S3 = FieldTimeSeries(grid, XYZ(), ts)
    ## check that trying to create a FieldInputSource with inconsistent types or dimensions fails
    @test_throws AssertionError InputSource(; S1, S2, S3)
    @test_throws AssertionError InputSource(; X1, S1, S2)
    # populate S1 with random data
    S1.data .= randn(size(S1))
    ## check update_inputs!
    fields = (S1 = Field(grid, XY()), S2 = Field(grid, XY()))
    clock = Clock(time = 0)
    update_inputs!(fields, fts_input, clock)
    @test all(fields.S1 .== S1[1])
    @test all(fields.S2 .== S2[1])
    # advance clock and check that inputs are updated
    Terrarium.tick_time!(clock, 1.0)
    update_inputs!(fields, fts_input, clock)
    @test all(fields.S1 .== S1[2])
    @test all(fields.S2 .== S2[2])
end
