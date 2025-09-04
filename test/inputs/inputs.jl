using Terrarium
using Test
using Unitful

@testset "InputFields" begin
    grid = ColumnGrid(ExponentialSpacing())
    inputs = InputFields(grid)
    @test isempty(get_input_fields(inputs, XY()))
    @test isempty(get_input_fields(inputs, XYZ()))
    test2D_input_field = get_input_field(inputs, :test2D, XY())
    # check that new input field was allocated
    @test haskey(get_input_fields(inputs, XY()), :test2D)
    # check that the 3D (XYZ) inputs are still empty
    @test isempty(get_input_fields(inputs, XYZ()))
    test3D_input_field = get_input_field(inputs, :test3D, XYZ())
    @test haskey(get_input_fields(inputs, XYZ()), :test3D)
    @test get_input_field(inputs, InputVariable(:test2D, XY(), NoUnits, "")) === test2D_input_field
    @test get_input_field(inputs, InputVariable(:test3D, XYZ(), NoUnits, "")) === test3D_input_field
end

@testset "Input sources" begin
    # Field input source
    grid = ColumnGrid(ExponentialSpacing())
    X1 = Field(grid, XY())
    field_input = FieldInputSource(; X1)
    ## check that dimensions were inferred correctly
    @test field_input.dims == XY()
    @test field_input.fields == (; X1)
    X2 = Field(grid, XY())
    field_input = FieldInputSource(; X1, X2)
    @test field_input.dims == XY()
    @test field_input.fields == (; X1, X2)
    X3 = Field(grid, XYZ())
    ## check that trying to create a FieldInputSource with inconsistent types or dimensions fails
    @test_throws AssertionError FieldInputSource(; X1, X2, X3)
    @test_throws AssertionError FieldInputSource(; X1, X2, Y=zeros(size(X1)))
    ## check update_inputs!
    inputs = InputFields(grid)
    clock = Clock(time=0)
    update_inputs!(inputs, field_input, clock)
    @test get_input_field(inputs, :X1, XY()) === X1
    @test get_input_field(inputs, :X2, XY()) === X2

    # FieldTimeSeries input source
    grid = ColumnGrid(ExponentialSpacing())
    ts = 0.0:1.0:10.0
    S1 = FieldTimeSeries(grid, XY(), ts)
    fts_input = FieldTimeSeriesInputSource(; S1)
    ## check that dimensions were inferred correctly
    @test fts_input.dims == XY()
    @test fts_input.fts == (; S1)
    S2 = FieldTimeSeries(grid, XY(), ts)
    fts_input = FieldTimeSeriesInputSource(; S1, S2)
    @test fts_input.dims == XY()
    @test fts_input.fts == (; S1, S2)
    S3 = FieldTimeSeries(grid, XYZ(), ts)
    ## check that trying to create a FieldInputSource with inconsistent types or dimensions fails
    @test_throws AssertionError FieldTimeSeriesInputSource(; S1, S2, S3)
    @test_throws AssertionError FieldTimeSeriesInputSource(; X1, S1, S2)
    # populate S1 with random data
    S1.data .= randn(size(S1))
    ## check update_inputs!
    inputs = InputFields(grid)
    clock = Clock(time=0)
    update_inputs!(inputs, fts_input, clock)
    @test all(get_input_field(inputs, :S1, XY()) .== S1[1])
    @test get_input_field(inputs, :S2, XY()) == S2[1]
    # advance clock and check that inputs are updated
    Terrarium.tick_time!(clock, 1.0)
    update_inputs!(inputs, fts_input, clock)
    @test all(get_input_field(inputs, :S1, XY()) .== S1[2])
    @test get_input_field(inputs, :S2, XY()) == S2[2]
end

@testset "InputProvider" begin
    # Field input source
    grid = ColumnGrid(ExponentialSpacing())
    X1 = Field(grid, XY())
    set!(X1, 1)
    field_input = FieldInputSource(; X1)
    provider = InputProvider(grid, field_input)
    clock = Clock(time=0)
    update_inputs!(provider, clock)
    @test get_input_field(provider.fields, :X1, XY()) === X1
    @test all(get_input_field(provider.fields, :X1, XY()) .== 1)
end
