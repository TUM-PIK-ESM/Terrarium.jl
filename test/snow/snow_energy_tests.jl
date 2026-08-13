using Terrarium
using Test

@testset "Snow energy closure" begin
    NF = Float64
    grid = ColumnGrid(CPU(), NF, ExponentialSpacing(N = 3))
    snow = SingleLayerSnow(NF)
    constants = PhysicalConstants(NF)
    closure = Terrarium.get_closure(snow)
    state = StateVariables(snow, grid)

    ρ_snow = Terrarium.snow_density(snow)
    ρ_w = constants.material.density_water
    L_f = constants.thermodynamics.latent_heat_fusion
    c_i = constants.thermodynamics.specific_heat_capacity_ice
    c_w = constants.thermodynamics.specific_heat_capacity_liquid_water
    W = NF(0.1)
    d_snow = W * ρ_w / ρ_snow
    Lθ = ρ_snow * L_f

    @testset "recovery: frozen (T < 0)" begin
        set!(state.snow_water_equivalent, W)
        T0 = NF(-5)
        U_v = T0 * ρ_snow * c_i - Lθ          # frozen: θ_liq = 0
        set!(state.snow_energy, U_v * d_snow)
        Terrarium.closure!(state, grid, closure, snow, constants)
        @test all(state.snow_temperature .≈ T0)
        @test all(state.snow_liquid_fraction .≈ 0)
    end

    @testset "phase change: T = 0, partial liquid" begin
        set!(state.snow_water_equivalent, W)
        set!(state.snow_energy, (-Lθ / 2) * d_snow)   # half-melted
        Terrarium.closure!(state, grid, closure, snow, constants)
        @test all(state.snow_temperature .≈ 0)
        @test all(state.snow_liquid_fraction .≈ 0.5)
    end

    @testset "temperature clipped at 0 for E > 0" begin
        set!(state.snow_water_equivalent, W)
        # E > 0: pack fully melted at 0°C plus a sensible excess (would give T = 2°C if unclipped)
        set!(state.snow_energy, NF(2) * ρ_snow * c_w * d_snow)
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

    @testset "thin snow with residual energy: temperature stays bounded" begin
        # Regression for the unphysically cold snow temperature: the explicit dynamics can leave a
        # vanishing pack (tiny W) holding a non-negligible depth-integrated energy Ē. Recovering the
        # volumetric energy U = Ē/d_snow with only an eps offset then yields an enormous |U| and a
        # temperature thousands of degrees below zero. Flooring the depth at `min_conduction_thickness`
        # bounds the effective heat capacity so the recovered temperature stays physical.
        set!(state.snow_water_equivalent, NF(1.0e-6))   # d_snow ≈ 3e-6 m, far below the 5 mm floor
        set!(state.snow_energy, NF(-1.0e4))             # residual energy inconsistent with the tiny pack
        Terrarium.closure!(state, grid, closure, snow, constants)
        @test all(isfinite.(state.snow_temperature))
        # bounded by ≈ Ē/(C_snow·d_min); without the floor this would be several thousand °C below zero
        @test all(state.snow_temperature .>= -100)
        @test all(state.snow_temperature .<= 0)
    end

    @testset "initialize! sets snow_energy from prescribed temperature" begin
        set!(state.snow_water_equivalent, W)
        set!(state.snow_temperature, NF(-3))
        Terrarium.initialize!(state, grid, snow, constants)
        @test all(isfinite.(state.snow_energy))
        # frozen initial state: energy should be negative
        @test all(state.snow_energy .< 0)
    end

    @testset "basal heat flux: thin-snow depth flooring" begin
        # Q_base = 2·κ·(T_soil − T_snow)/max(d_snow, d_min): the denominator is floored at d_min so the
        # flux stays bounded as the snowpack vanishes (an eps offset would let it diverge and destabilize
        # explicit stepping over frozen soil).
        κ = NF(0.16)
        T_soil = NF(-10)   # frozen soil below
        T_snow = NF(0)
        d_min = snow.min_conduction_thickness
        # thick snow (d_snow ≫ d_min): floor inactive, recovers the plain conduction formula
        d_thick = NF(0.5)
        @test Terrarium.compute_snow_basal_heat_flux(κ, T_soil, T_snow, d_thick, d_min) ≈ 2κ * (T_soil - T_snow) / d_thick
        # vanishing snow: flux is finite and set by the d_min floor, not the near-zero depth
        Q_floor = 2κ * (T_soil - T_snow) / d_min
        for d_thin in (NF(0), NF(1.0e-9), NF(1.0e-6))
            Q = Terrarium.compute_snow_basal_heat_flux(κ, T_soil, T_snow, d_thin, d_min)
            @test isfinite(Q)
            @test Q ≈ Q_floor
        end
    end
end
