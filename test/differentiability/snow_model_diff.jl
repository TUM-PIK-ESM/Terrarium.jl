using Terrarium
using Test

using Enzyme

function build_snow_model(arch, ::Type{NF}) where {NF}
    grid = ColumnGrid(arch, NF, ExponentialSpacing(N = 3))
    model = SnowModel(grid)
    return model
end

# Differentiability of the snow-specific physics. These are function-level adjoint checks: the snow energy
# closure reuses the medium-agnostic `FreeWater` maps (`liquid_water_fraction`, `energy_to_temperature`),
# whose adjoints are already covered in `soil_energy_diff.jl`, so here we exercise the snow-specific
# diagnostics and fluxes.

@testset "Snow meltwater outflow: differentiability" begin
    # Darcy outflow M_r = K_sat·S*³; gradient w.r.t. the liquid fraction
    hydraulics = ConstantSnowHydraulics(Float64)
    ## above capillary retention -> nonzero, smooth gradient
    grad, = Enzyme.autodiff(Reverse, Terrarium.compute_meltwater_outflow, Active, Const(hydraulics), Active(0.5))
    @test isfinite(grad[2])
    @test grad[2] > 0
    ## below capillary retention -> no outflow, zero gradient
    grad, = Enzyme.autodiff(Reverse, Terrarium.compute_meltwater_outflow, Active, Const(hydraulics), Active(0.01))
    @test iszero(grad[2])
end

@testset "Snow cover fraction: differentiability" begin
    # f_snow = W/(W + W_ref); gradient w.r.t. SWE is positive for W > 0 and zero for clamped W < 0
    cover = FractionalSnowCover(Float64)
    grad, = Enzyme.autodiff(Reverse, Terrarium.compute_snow_cover_fraction, Active, Const(cover), Active(0.05))
    @test isfinite(grad[2])
    @test grad[2] > 0
    grad, = Enzyme.autodiff(Reverse, Terrarium.compute_snow_cover_fraction, Active, Const(cover), Active(-1.0))
    @test iszero(grad[2])   # clamped: ∂f/∂W = 0 for W < 0
end

@testset "Snow depth: differentiability" begin
    # d_snow = W·ρ_w/ρ_snow; gradient w.r.t. SWE is ρ_w/ρ_snow for W > 0
    snow = SingleLayerSnow(Float64)
    ρ_snow = Terrarium.snow_density(snow)
    ρ_w = 1000.0
    grad, = Enzyme.autodiff(Reverse, Terrarium.compute_snow_depth, Active, Const(snow), Active(0.1), Const(ρ_snow), Const(ρ_w))
    @test grad[2] ≈ ρ_w / ρ_snow
end

@testset "Snow thermal conductivity: differentiability" begin
    # κ = a·(ρ_snow/ρ_w)^b; gradient w.r.t. bulk density is finite and positive
    cond = PowerLawSnowThermalConductivity(Float64)
    material = PhysicalConstants(Float64).material
    grad, = Enzyme.autodiff(Reverse, Terrarium.compute_thermal_conductivity, Active, Const(cond), Const(material), Active(250.0))
    @test isfinite(grad[3])
    @test grad[3] > 0
end

@testset "Snow basal heat flux: differentiability" begin
    # Q_base = 2·κ·(T_soil − T_snow)/max(d_snow, d_min); gradient w.r.t. the soil temperature is 2·κ/d_snow for d_snow > d_min
    κ = 0.3
    d_snow = 0.5
    d_min = 1.0e-3
    grad, = Enzyme.autodiff(Reverse, Terrarium.compute_snow_basal_heat_flux, Active, Const(κ), Active(1.0), Const(-2.0), Const(d_snow), Const(d_min))
    @test grad[2] ≈ 2κ / d_snow
end

@testset "Snow model: timestep!" begin
    model = build_snow_model(CPU(), Float64)
    integrator = initialize(
        model;
        initializers = (
            snow_water_equivalent = 0.3,
            snow_temperature = -5.0,
            surface_heat_flux = 50.0,
            basal_heat_flux = 0.0,
        )
    )
    state = integrator.state
    # perturb into the phase-change band (partial liquid fraction) so the Darcy meltwater outflow
    # nonlinearity is active, then resync the closure so `snow_temperature`/`snow_liquid_fraction` are consistent
    set!(state.snow_energy, state.snow_energy[1, 1, 1] / 4)
    Terrarium.closure!(state, model)
    @test 0 < state.snow_liquid_fraction[1, 1, 1] < 1

    dintegrator = make_zero(integrator)
    dintegrator.state.snow_energy .= 1.0
    dintegrator.state.snow_water_equivalent .= 1.0
    Δt = get_timestepper(integrator.model).Δt
    @time Enzyme.autodiff(set_runtime_activity(Reverse), timestep!, Const, Duplicated(integrator, dintegrator), Const(Δt))
    @test all(isfinite.(dintegrator.state.snow_water_equivalent))
    @test all(isfinite.(dintegrator.state.snow_energy))
    @test all(isfinite.(dintegrator.state.snow_liquid_fraction))
end
