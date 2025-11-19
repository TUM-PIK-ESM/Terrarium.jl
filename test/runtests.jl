using Terrarium
using Test

FLAG_ENZYME_TESTS = "enzyme" in ARGS ? true : false
MAIN_TESTS = !FLAG_ENZYME_TESTS ? true : false

if FLAG_ENZYME_TESTS
    @testset "Enzyme" begin
        include("differentiability/soil_energy_diff.jl")
        include("differentiability/soil_hydrology_diff.jl")
        include("differentiability/vegetation_model_diff.jl")
    end
end

if MAIN_TESTS
    @testset "Grids" begin
        include("grids.jl")
    end

    @testset "State variables" begin
        include("state_variables.jl")
    end
    
    @testset "Timestepping" begin
        include("timestepping/explicit_step.jl")
        include("timestepping/run_simulation.jl")
        include("timestepping/heun.jl")
    end

    @testset "Inputs" begin
        include("inputs/inputs.jl")
        include("inputs/input_forcing.jl")
    end
    
    @testset "Soil model and processes" begin
        include("soil/soil_model_tests.jl")
    end

    @testset "Vegetation processes" begin
        include("vegetation/vegetation_model_tests.jl")
    end

    @testset "Surface energy balance" begin
        include("surface_energy/seb_tests.jl")
    end
end 