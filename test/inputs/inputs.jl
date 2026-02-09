using Terrarium
using Terrarium: FieldInputSource, FieldTimeSeriesInputSource, Variables
using Test
using Unitful

@testset "Input sources" begin
    # Field input source
    grid = ColumnGrid(ExponentialSpacing())
    field_input = InputSource(eltype(grid), :X1)
    @test isa(field_input, FieldInputSource)
    ## check that dimensions were inferred correctly
    @test field_input.dims == XY()
    @test variables(field_input) == (Terrarium.input(:X1, XY()),)
    field_input = InputSource(eltype(grid), :X1, :X2)
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
