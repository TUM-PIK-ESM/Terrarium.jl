using Terrarium
using Test

using Terrarium: compute_snow_depth, compute_snow_cover_fraction, compute_thermal_conductivity,
    compute_snow_volumetric_energy, snow_density

@testset "Snow properties" begin
    NF = Float64
    snow = SingleLayerSnow(NF)
    constants = PhysicalConstants(NF)
    ρ_w = NF(1000)

    @testset "density" begin
        # the snow process delegates to its density scheme (constant by default)
        @test snow.density isa ConstantSnowDensity
        @test snow_density(snow) == snow_density(snow.density)
        @test snow_density(snow) > 0
        @test snow_density(ConstantSnowDensity(NF; density = NF(300))) == 300
    end

    @testset "snow depth" begin
        # d_snow = W·ρ_w/ρ_snow: zero at W=0, increasing in W, negative W clamped to zero
        ρ_snow = snow_density(snow)
        @test compute_snow_depth(snow, NF(0), ρ_snow, ρ_w) == 0
        @test compute_snow_depth(snow, NF(-1), ρ_snow, ρ_w) == 0
        d1 = compute_snow_depth(snow, NF(0.1), ρ_snow, ρ_w)
        d2 = compute_snow_depth(snow, NF(0.2), ρ_snow, ρ_w)
        @test d1 ≈ 0.1 * ρ_w / snow_density(snow)
        @test d2 > d1 > 0
    end

    @testset "cover fraction" begin
        # f_snow = W/(W+W_ref): in [0,1), zero at W=0, increasing, → 1 as W → ∞
        cover = snow.cover
        @test compute_snow_cover_fraction(cover, NF(0)) == 0
        @test compute_snow_cover_fraction(cover, NF(-1)) == 0
        f_small = compute_snow_cover_fraction(cover, NF(0.005))
        f_large = compute_snow_cover_fraction(cover, NF(0.05))
        @test 0 < f_small < f_large < 1
        @test compute_snow_cover_fraction(cover, NF(1.0e6)) ≈ 1 atol = 1.0e-3
        # at W = W_ref the fraction is exactly one half
        @test compute_snow_cover_fraction(cover, cover.half_coverage) ≈ 0.5
    end

    @testset "thermal conductivity" begin
        # κ = a·(ρ_snow/ρ_w)^b: positive and increasing with bulk density
        ρ_snow = snow_density(snow)
        κ = compute_thermal_conductivity(snow, constants.material, ρ_snow)
        cond = snow.thermal_conductivity
        @test κ > 0
        @test κ ≈ cond.conductivity_coefficient * (ρ_snow / ρ_w)^cond.conductivity_exponent
        denser = SingleLayerSnow(NF; density = ConstantSnowDensity(NF; density = NF(400)))
        @test compute_thermal_conductivity(denser, constants.material, snow_density(denser)) > κ
    end

    @testset "volumetric energy" begin
        d_min = NF(5.0e-3)
        # U_v = Ē/d_snow for d_snow > d_min
        @test compute_snow_volumetric_energy(NF(-1000), NF(0.5), d_min) ≈ -2000 rtol = 1.0e-6
        @test compute_snow_volumetric_energy(NF(0), NF(0.5), d_min) == 0
        # thin snow: depth is floored at d_min so the recovered volumetric energy (hence temperature)
        # stays bounded instead of diverging as d_snow → 0
        @test compute_snow_volumetric_energy(NF(-1000), NF(1.0e-9), d_min) ≈ NF(-1000) / d_min rtol = 1.0e-6
        @test compute_snow_volumetric_energy(NF(0), NF(0), d_min) == 0
    end

    @testset "Float32 type stability" begin
        snow32 = SingleLayerSnow(Float32)
        material32 = PhysicalConstants(Float32).material
        @test compute_snow_depth(snow32, 0.1f0, snow_density(snow32), 1000.0f0) isa Float32
        @test compute_snow_cover_fraction(snow32.cover, 0.02f0) isa Float32
        @test compute_thermal_conductivity(snow32, material32, snow_density(snow32)) isa Float32
    end

    @testset "compute_auxiliary! diagnoses properties" begin
        grid = ColumnGrid(CPU(), NF, ExponentialSpacing(N = 10))
        constants = PhysicalConstants(NF)
        state = StateVariables(snow, grid)
        W = NF(0.1)
        set!(state.snow_water_equivalent, W)
        compute_auxiliary!(state, grid, snow, constants)
        ρ_w = constants.material.density_water
        @test all(state.snow_depth .≈ compute_snow_depth(snow, W, snow_density(snow), ρ_w))
        @test all(state.snow_cover_fraction .≈ compute_snow_cover_fraction(snow.cover, W))
        @test all(isfinite.(state.snow_depth))
    end
end
