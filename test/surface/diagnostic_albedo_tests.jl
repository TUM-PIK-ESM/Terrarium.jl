using Terrarium
using Test

@testset "DiagnosticAlbedo" begin
    NF = Float64
    grid = ColumnGrid(CPU(), NF, ExponentialSpacing(N = 3))
    alb = DiagnosticAlbedo(NF)
    snow = SingleLayerSnow(NF)

    @testset "no snow component -> background values" begin
        state = StateVariables(Terrarium.Variables(Terrarium.variables(alb)...), grid)
        compute_auxiliary!(state, grid, alb, nothing)
        @test all(state.albedo .≈ alb.background_albedo)
        @test all(state.emissivity .≈ alb.background_emissivity)
    end

    @testset "blends background and snow by snow cover fraction" begin
        vars = Terrarium.Variables(Terrarium.variables(snow)..., Terrarium.variables(alb)...)
        state = StateVariables(vars, grid)

        set!(state.snow_cover_fraction, 0.5)
        compute_auxiliary!(state, grid, alb, snow)
        @test all(state.albedo .≈ 0.5 * alb.background_albedo + 0.5 * alb.snow_albedo)
        @test all(state.emissivity .≈ 0.5 * alb.background_emissivity + 0.5 * alb.snow_emissivity)

        # full snow cover recovers the snow endpoints; bare cover recovers the background
        set!(state.snow_cover_fraction, 1.0)
        compute_auxiliary!(state, grid, alb, snow)
        @test all(state.albedo .≈ alb.snow_albedo)
        @test all(state.emissivity .≈ alb.snow_emissivity)

        set!(state.snow_cover_fraction, 0.0)
        compute_auxiliary!(state, grid, alb, snow)
        @test all(state.albedo .≈ alb.background_albedo)
    end
end
