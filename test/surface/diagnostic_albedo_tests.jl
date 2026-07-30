using Terrarium
using Test

@testset "DiagnosticAlbedo" begin
    NF = Float64
    grid = ColumnGrid(CPU(), NF, ExponentialSpacing(N = 3))
    albedo = DiagnosticAlbedo(NF)
    snow = SingleLayerSnow(NF)

    @testset "no snow component -> background values" begin
        state = StateVariables(Terrarium.Variables(Terrarium.variables(albedo)...), grid)
        compute_auxiliary!(state, grid, albedo, nothing)
        @test all(state.albedo .≈ albedo.background_albedo)
        @test all(state.emissivity .≈ albedo.background_emissivity)
    end

    @testset "blends background and snow by snow cover fraction" begin
        vars = Terrarium.Variables(Terrarium.variables(snow)..., Terrarium.variables(albedo)...)
        state = StateVariables(vars, grid)

        set!(state.snow_cover_fraction, 0.5)
        compute_auxiliary!(state, grid, albedo, nothing, snow)
        @test all(state.albedo .≈ 0.5 * albedo.background_albedo + 0.5 * snow.albedo.snow_albedo)
        @test all(state.emissivity .≈ 0.5 * albedo.background_emissivity + 0.5 * snow.albedo.snow_emissivity)

        # full snow cover recovers the snow endpoints; bare cover recovers the background
        set!(state.snow_cover_fraction, 1.0)
        compute_auxiliary!(state, grid, albedo, nothing, snow)
        @test all(state.albedo .≈ snow.albedo.snow_albedo)
        @test all(state.emissivity .≈ snow.albedo.snow_emissivity)

        set!(state.snow_cover_fraction, 0.0)
        compute_auxiliary!(state, grid, albedo, nothing, snow)
        @test all(state.albedo .≈ albedo.background_albedo)
    end
end
