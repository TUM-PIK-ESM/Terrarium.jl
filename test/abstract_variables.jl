using Terrarium
using Test
using Oceananigans: Center, Face
using Unitful

@testset "varpath" begin
    @test varpath(:x) == (:x,)
    @test varpath(:ns => :x) == (:ns, :x)
    @test varpath(:a => :b => :c) == (:a, :b, :c)
    @test varpath((:ns, :x)) == (:ns, :x)
end

@testset "VarDims types" begin
    # Test XYZ dimensions
    xyz = XYZ()
    @test xyz.x isa Center
    @test xyz.y isa Center
    @test xyz.z isa Center

    xyz_custom = XYZ(x = Face(), y = Center(), z = Face())
    @test xyz_custom.x isa Face
    @test xyz_custom.y isa Center
    @test xyz_custom.z isa Face

    # Test XY dimensions
    xy = XY()
    @test xy.x isa Center
    @test xy.y isa Center

    xy_custom = XY(x = Face(), y = Center())
    @test xy_custom.x isa Face
    @test xy_custom.y isa Center
end

@testset "Variable construction" begin
    # Test basic Variable creation
    var = Variable(:temperature, XYZ(), u"K")
    @test varname(var) == :temperature
    @test vardims(var) isa XYZ
    @test varunits(var) == u"K"

    # Test with default units (NoUnits)
    var_no_units = Variable(:pressure, XY())
    @test varname(var_no_units) == :pressure
    @test varunits(var_no_units) == Terrarium.NoUnits

    # Test convenience constructor
    var_convenience = var(:saturation, XYZ(x = Center(), y = Center(), z = Face()), u"1")
    @test varname(var_convenience) == :saturation
    @test vardims(var_convenience).z isa Face
end

@testset "PrognosticVariable construction" begin
    # Test basic prognostic variable
    prog = prognostic(:temperature, XYZ(); units = u"K")
    @test prog.var.name == :temperature
    @test hasclosure(prog) == false
    @test isa(prog.closure, Nothing)
    @test isa(prog.tendency, AuxiliaryVariable)
    @test varname(prog.tendency.var) == :temperature

    # Test with closure relation
    struct MockClosure <: Terrarium.AbstractClosureRelation end
    prog_with_closure = prognostic(:enthalpy, XYZ(); units = u"J", closure = MockClosure())
    @test hasclosure(prog_with_closure) == true
    @test isa(prog_with_closure.closure, MockClosure)

    # Test with bounds and description
    prog_full = prognostic(
        :moisture,
        XY();
        units = u"kg",
        bounds = (0.0, 1.0),
        desc = "Soil moisture content"
    )
    @test prog_full.desc == "Soil moisture content"
end

@testset "AuxiliaryVariable construction" begin
    # Test basic auxiliary variable without constructor
    aux = auxiliary(:pressure, XYZ(); units = u"Pa")
    @test varname(aux) == :pressure
    @test isa(aux.ctor, Nothing)

    # Test with constructor function
    mock_ctor(params, args...) = params * sum(args)
    aux_with_ctor = auxiliary(:flux, XY(), mock_ctor, 2.0; units = u"m/s")
    @test aux_with_ctor.ctor(3.0, 1.0, 2.0) ≈ 6.0

    # Test with bounds
    aux_bounded = auxiliary(:temperature, XYZ(); units = u"K", bounds = (200.0, 400.0))
    @test aux_bounded.bounds == (200.0, 400.0)
end

@testset "InputVariable construction" begin
    # Test basic input variable without default
    inp = input(:radiation, XYZ(); units = u"W/m^2")
    @test varname(inp) == :radiation
    @test isa(inp.default, Nothing)

    # Test with numeric default
    inp_default_num = input(:precipitation, XY(); default = 0.01, units = u"m/s")
    @test inp_default_num.default ≈ 0.01

    # Test with function default
    f(t) = sin(t)
    inp_default_func = input(:forcing, XYZ(); default = f, units = u"K")
    @test inp_default_func.default(π / 2) ≈ 1.0
end

@testset "tendency helper" begin
    # Test that tendency creates appropriate auxiliary variable
    temp_var = var(:temperature, XYZ(), u"K")
    tend = tendency(temp_var)

    @test varname(tend) == :temperature
    @test vardims(tend) isa XYZ
    # Check units are K/s (upreferred)
    @test string(varunits(tend)) == "K s⁻¹" || string(varunits(tend)) == "K / s"
end

@testset "Namespace construction" begin
    # Test basic namespace creation
    vars = Variables(prognostic(:temp, XYZ(); units = u"K"))
    ns = namespace(:soil, vars)

    @test varname(ns) == :soil
    @test isa(ns.vars.prognostic.temp, PrognosticVariable)

    # Test namespace from tuple of variables
    ns_from_tuple = namespace(
        :atmosphere,
        (prognostic(:humidity, XY(); units = u"kg"), auxiliary(:pressure, XYZ(); units = u"Pa"))
    )
    @test varname(ns_from_tuple) == :atmosphere
    @test hasproperty(ns_from_tuple.vars.prognostic, :humidity)
    @test hasproperty(ns_from_tuple.vars.auxiliary, :pressure)
end

@testset "Variables construction" begin
    # Test basic Variables creation from tuple
    vars = Variables(
        prognostic(:temperature, XYZ(); units = u"K"),
        auxiliary(:pressure, XY(); units = u"Pa"),
        input(:radiation, XYZ(); default = 100.0, units = u"W/m^2")
    )

    @test length(vars.prognostic) == 1
    @test length(vars.auxiliary) == 1
    @test length(vars.inputs) == 1
    @test varname(first(values(vars.prognostic))) == :temperature

    # Test with namespaces
    vars_with_ns = Variables(
        prognostic(:temp, XYZ(); units = u"K"),
        namespace(:soil, (prognostic(:moisture, XY(); units = u"kg"),))
    )

    @test length(vars_with_ns.namespaces) == 1
    @test varname(first(values(vars_with_ns.namespaces))) == :soil

    # Test automatic merging of namespaces with same name
    vars_merged = Variables(
        namespace(:atmosphere, (prognostic(:temp, XYZ(); units = u"K"),)),
        namespace(:atmosphere, (auxiliary(:pressure, XY(); units = u"Pa"),))
    )

    @test length(vars_merged.namespaces) == 1
    ns = first(values(vars_merged.namespaces))
    @test hasproperty(ns.vars.prognostic, :temp)
    @test hasproperty(ns.vars.auxiliary, :pressure)
end

@testset "Variables duplicate detection" begin
    # Test that duplicate variable names raise an error
    @test_throws ErrorException Variables(
        prognostic(:temperature, XYZ(); units = u"K"),
        auxiliary(:temperature, XY(); units = u"Pa")  # Same name!
    )

    # Test that duplicate namespace names are merged (not an error)
    @test_logs (:warn,) Variables(
        namespace(:soil, (prognostic(:temp, XYZ(); units = u"K"),)),
        namespace(:soil, (auxiliary(:pressure, XY(); units = u"Pa"),))
    )
end

@testset "Variable helper functions" begin
    # Test prognostic_variables
    vars = Variables(
        prognostic(:temp, XYZ(); units = u"K"),
        prognostic(:moisture, XY(); units = u"kg"),
        auxiliary(:pressure, XYZ(); units = u"Pa"),
        input(:radiation, XYZ(); default = 100.0)
    )

    prog_vars = prognostic_variables(vars)
    @test length(prog_vars) == 2
    @test varname(first(prog_vars)) in [:temp, :moisture]

    # Test auxiliary_variables
    aux_vars = auxiliary_variables(vars)
    @test length(aux_vars) == 1
    @test varname(first(aux_vars)) == :pressure

    # Test input_variables
    inp_vars = input_variables(vars)
    @test length(inp_vars) == 1
    @test varname(first(inp_vars)) == :radiation
end

@testset "Variable show methods" begin
    # Test that Variable shows correctly
    var = Variable(:temperature, XYZ(), u"K")
    io = IOBuffer()
    show(io, MIME"text/plain"(), var)
    str = String(take!(io))
    @test occursin("temperature", str)
    @test occursin("K", str)

    # Test Variable without units
    var_no_units = Variable(:pressure, XY())
    io = IOBuffer()
    show(io, MIME"text/plain"(), var_no_units)
    str = String(take!(io))
    @test occursin("pressure", str)
end

@testset "Variables show methods" begin
    vars = Variables(
        prognostic(:temp, XYZ(); units = u"K"),
        auxiliary(:pressure, XY(); units = u"Pa")
    )

    io = IOBuffer()
    show(io, vars)
    str = String(take!(io))

    @test occursin("Variables", str)
    @test occursin("Prognostic", str)
    @test occursin("Auxiliary", str)
    @test occursin("temp", str)
    @test occursin("pressure", str)
end

@testset "merge Variables" begin
    vars1 = Variables(
        prognostic(:temp, XYZ(); units = u"K")
    )

    vars2 = Variables(
        auxiliary(:pressure, XY(); units = u"Pa")
    )

    merged = merge(vars1, vars2)

    @test length(merged.prognostic) == 1
    @test length(merged.auxiliary) == 1
    @test hasproperty(merged.prognostic, :temp)
    @test hasproperty(merged.auxiliary, :pressure)
end

@testset "namespaces helper" begin
    nt = (
        soil = (prognostic(:temp, XYZ(); units = u"K"),),
        atmosphere = (auxiliary(:pressure, XY(); units = u"Pa"),),
    )

    ns_tuple = namespaces(nt)

    @test length(ns_tuple) == 2
    names = map(varname, ns_tuple)
    @test :soil in names
    @test :atmosphere in names
end

@testset "with_scope" begin
    var = input(:temperature, XYZ(); units = u"K")

    # Test root scope
    result = with_scope((), var)
    @test length(result) == 1
    @test varname(first(result)) == :temperature

    # Test single namespace
    result = with_scope((:soil,), var)
    @test length(result) == 1
    ns = first(result)
    @test varname(ns) == :soil
    @test varname(first(values(ns.vars.inputs))) == :temperature

    # Test nested namespaces
    result = with_scope((:atmosphere, :boundary, :temp), var)
    @test length(result) == 1
    ns1 = first(result)
    @test varname(ns1) == :atmosphere
    ns2 = first(values(ns1.vars.namespaces))
    @test varname(ns2) == :boundary
    @test varname(first(values(ns2.vars.inputs))) == :temperature
end
