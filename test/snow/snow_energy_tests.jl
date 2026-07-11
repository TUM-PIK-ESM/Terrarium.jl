using Terrarium
using Test

@testset "Snow energy closure" begin
    NF = Float64
    grid = ColumnGrid(CPU(), NF, ExponentialSpacing(N = 3))
    snow = SingleLayerSnow(NF)
    constants = PhysicalConstants(NF)
    closure = Terrarium.get_closure(snow)
    state = StateVariables(snow, grid)

    ρ_s = Terrarium.snow_density(snow)
    ρ_w = constants.material.density_water
    L_f = constants.thermodynamics.latent_heat_fusion
    c_i = constants.thermodynamics.specific_heat_capacity_ice
    c_w = constants.thermodynamics.specific_heat_capacity_liquid_water
    W = NF(0.1)
    d_s = W * ρ_w / ρ_s
    Lθ = ρ_s * L_f

    @testset "recovery: frozen (T < 0)" begin
        set!(state.snow_water_equivalent, W)
        T0 = NF(-5)
        U_v = T0 * ρ_s * c_i - Lθ          # frozen: θ_liq = 0
        set!(state.snow_energy, U_v * d_s)
        Terrarium.closure!(state, grid, closure, snow, constants)
        @test all(state.snow_temperature .≈ T0)
        @test all(state.snow_liquid_fraction .≈ 0)
    end

    @testset "phase change: T = 0, partial liquid" begin
        set!(state.snow_water_equivalent, W)
        set!(state.snow_energy, (-Lθ / 2) * d_s)   # half-melted
        Terrarium.closure!(state, grid, closure, snow, constants)
        @test all(state.snow_temperature .≈ 0)
        @test all(state.snow_liquid_fraction .≈ 0.5)
    end

    @testset "temperature clipped at 0 for E > 0" begin
        set!(state.snow_water_equivalent, W)
        # E > 0: pack fully melted at 0°C plus a sensible excess (would give T = 2°C if unclipped)
        set!(state.snow_energy, NF(2) * ρ_s * c_w * d_s)
        Terrarium.closure!(state, grid, closure, snow, constants)
        @test all(state.snow_temperature .≈ 0)     # clipped at 0
        @test all(state.snow_liquid_fraction .≈ 1)
    end

    @testset "round-trip invclosure -> closure (T ≤ 0)" begin
        for T0 in (NF(-8), NF(-0.5), NF(0))
            set!(state.snow_water_equivalent, NF(0.05))
            set!(state.snow_temperature, T0)
            Terrarium.invclosure!(state, grid, closure, snow, constants)
            # scribble over the temperature to ensure the forward closure recomputes it
            set!(state.snow_temperature, NF(-999))
            Terrarium.closure!(state, grid, closure, snow, constants)
            @test all(state.snow_temperature .≈ T0)
        end
    end

    @testset "W -> 0 stays finite (no Inf/NaN)" begin
        set!(state.snow_water_equivalent, NF(0))
        set!(state.snow_energy, NF(0))
        Terrarium.closure!(state, grid, closure, snow, constants)
        @test all(isfinite.(state.snow_temperature))
        @test all(isfinite.(state.snow_liquid_fraction))
    end

    @testset "initialize! sets snow_energy from prescribed temperature" begin
        set!(state.snow_water_equivalent, W)
        set!(state.snow_temperature, NF(-3))
        Terrarium.initialize!(state, grid, snow, constants)
        @test all(isfinite.(state.snow_energy))
        # frozen initial state: energy should be negative
        @test all(state.snow_energy .< 0)
    end
end
