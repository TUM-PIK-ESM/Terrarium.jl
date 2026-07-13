using Terrarium
using Test

@testset "Snow model (standalone)" begin
    NF = Float64
    grid = ColumnGrid(CPU(), NF, ExponentialSpacing(N = 3))
    model = SnowModel(grid)
    snow = model.snow
    constants = model.constants

    ρ_s = Terrarium.snow_density(snow)
    ρ_w = constants.material.density_water
    L_f = constants.thermodynamics.latent_heat_fusion
    c_i = constants.thermodynamics.specific_heat_capacity_ice
    W0 = NF(0.1)
    d_s = W0 * ρ_w / ρ_s

    # Build a fresh state with all forcing zeroed; caller sets the prognostic/forcing of interest.
    function fresh_state()
        state = StateVariables(model)
        for f in (:snowfall, :rainfall, :air_temperature, :surface_heat_flux, :basal_heat_flux, :sublimation)
            set!(getproperty(state, f), 0)
        end
        return state
    end

    # Run the auxiliary + closure passes then the tendencies (as a timestep would).
    function step_tendencies!(state)
        compute_auxiliary!(state, model)
        Terrarium.closure!(state, model)
        compute_tendencies!(state, model)
        return state
    end

    @testset "frozen, no forcing: zero tendencies" begin
        state = fresh_state()
        set!(state.snow_water_equivalent, W0)
        set!(state.snow_temperature, NF(-5))
        Terrarium.initialize!(state, model.grid, snow, constants)   # invclosure: T, W -> E
        step_tendencies!(state)
        @test all(state.snow_liquid_fraction .≈ 0)                  # frozen
        @test all(state.tendencies.snow_water_equivalent .≈ 0)      # no melt, no precip/sublimation
        @test all(state.tendencies.snow_energy .≈ 0)                # no fluxes
    end

    @testset "snowfall accumulates SWE and advects cold" begin
        state = fresh_state()
        set!(state.snow_water_equivalent, W0)
        set!(state.snow_temperature, NF(-5))
        Terrarium.initialize!(state, model.grid, snow, constants)
        P_s = NF(1e-6)
        set!(state.snowfall, P_s)
        set!(state.air_temperature, NF(-2))
        step_tendencies!(state)
        @test all(state.tendencies.snow_water_equivalent .≈ P_s)    # dW/dt = snowfall
        # cold snowfall (T < 0) advects negative energy relative to ice at 0°C
        @test all(state.tendencies.snow_energy .≈ ρ_w * P_s * c_i * NF(-2))
    end

    @testset "fully melted pack drains (SWE + energy decrease)" begin
        state = fresh_state()
        set!(state.snow_water_equivalent, W0)
        set!(state.snow_energy, NF(1000))                            # E > 0 -> θ_liq = 1 (fully melted)
        step_tendencies!(state)
        @test all(state.snow_liquid_fraction .≈ 1)
        M_r = snow.hydraulic_properties.saturated_conductivity      # S* = 1 at θ_liq = 1
        @test all(state.tendencies.snow_water_equivalent .≈ -M_r)   # meltwater drains
        @test all(state.tendencies.snow_energy .≈ -ρ_w * L_f * M_r) # latent heat leaves with meltwater
    end

    @testset "no drainage below capillary retention" begin
        state = fresh_state()
        set!(state.snow_water_equivalent, W0)
        Lθ = ρ_s * L_f
        # choose θ_liq just below L_c via a phase-change energy: θ_liq = 1 + U_v/Lθ
        θ_target = snow.hydraulic_properties.capillary_retention / 2
        set!(state.snow_energy, (θ_target - 1) * Lθ * d_s)
        step_tendencies!(state)
        @test all(0 .< state.snow_liquid_fraction .< snow.hydraulic_properties.capillary_retention)
        @test all(state.tendencies.snow_water_equivalent .≈ 0)      # retained, no outflow
    end

    @testset "basal flux warms, surface flux cools (positive upward)" begin
        state = fresh_state()
        set!(state.snow_water_equivalent, W0)
        set!(state.snow_temperature, NF(-5))                        # frozen -> no melt
        Terrarium.initialize!(state, model.grid, snow, constants)
        set!(state.basal_heat_flux, NF(10))                         # soil -> snow: gain
        set!(state.surface_heat_flux, NF(30))                       # snow -> sky: loss
        step_tendencies!(state)
        @test all(state.tendencies.snow_energy .≈ NF(10) - NF(30))  # dE/dt = G_base - G_top
        @test all(state.tendencies.snow_water_equivalent .≈ 0)
    end

    @testset "sublimation removes SWE" begin
        state = fresh_state()
        set!(state.snow_water_equivalent, W0)
        set!(state.snow_temperature, NF(-5))
        Terrarium.initialize!(state, model.grid, snow, constants)
        E_subl = NF(2e-7)
        set!(state.sublimation, E_subl)
        step_tendencies!(state)
        @test all(state.tendencies.snow_water_equivalent .≈ -E_subl)
    end
end
