using Terrarium
using Thermodynamics
using Oceananigans
using Test

"""
    test_skin_temperature_solve!(state, grid, seb, model; kwargs...)

Iteratively solve for the implicit skin temperature until convergence by alternating between
computing surface energy fluxes and updating the skin temperature. All input values are set
inside this function from keyword arguments.

# Arguments
- `state`: initialized model state
- `grid`: column grid
- `seb`: [`SurfaceEnergyBalance`](@ref) containing an [`ImplicitSkinTemperature`](@ref)
- `model`: [`SurfaceEnergyModel`](@ref)

# Keyword arguments (all required unless noted)
- `surface_shortwave_down::Float64`: downward shortwave radiation (W/m²)
- `surface_longwave_down::Float64`: downward longwave radiation (W/m²)
- `air_pressure::Float64`: air pressure at the surface (Pa)
- `air_temperature::Float64`: air temperature (°C)
- `relative_humidity::Float64`: relative humidity (0–1)
- `ground_temperature::Float64`: ground temperature (°C)
- `windspeed::Float64`: wind speed at reference height (m/s)
- `κₛ::Float64`: ground thermal conductivity (W/(m·K)) — default 0.5
- `Cₕ::Float64`: aerodynamic heat transfer coefficient — default 2.0e-3
- `max_iterations::Int`: maximum iterations — default 20
- `tolerance::Float64`: convergence tolerance — default `sqrt(eps())`

# Returns
Named tuple with fields: `skin_temperature`, `residual`, `net_radiation`, `latent_heat_flux`,
`sensible_heat_flux`, `ground_heat_flux`.
"""
function test_skin_temperature_solve!(
        state, grid, seb, model;
        surface_shortwave_down,
        surface_longwave_down,
        air_pressure,
        air_temperature,
        relative_humidity,
        ground_temperature,
        windspeed
    )
    set!(state.surface_shortwave_down, surface_shortwave_down)
    set!(state.surface_longwave_down, surface_longwave_down)
    set!(state.air_pressure, air_pressure)
    set!(state.air_temperature, air_temperature)
    thermodynamics_constants = model.constants.thermodynamics
    air_density = Thermodynamics.air_density(thermodynamics_constants, air_temperature + 273.15, air_pressure)
    specific_humidity_saturation = Terrarium.saturation_specific_humidity_vapor(thermodynamics_constants, air_temperature, air_density)
    set!(state.specific_humidity, relative_humidity * specific_humidity_saturation)
    set!(state.ground_temperature, ground_temperature)
    set!(state.windspeed, windspeed)
    # Seed the implicit solve with a physical initial guess (the ground temperature),
    # mirroring the SEB initializer used in the coupled model.
    set!(state.skin_temperature, ground_temperature)

    @time compute_auxiliary!(state, model)
    @test all(isfinite.(state.skin_temperature))

    # Check if ground heat flux converged to gradient flux
    field_grid = get_field_grid(grid)
    Δz = Oceananigans.Δzᵃᵃᶜ(1, 1, field_grid.Nz, field_grid)
    κₛ = model.surface_energy_balance.skin_temperature.κₛ
    Ts = state.skin_temperature[1, 1]
    Tg = state.ground_temperature[1, 1]
    G = state.ground_heat_flux[1, 1]
    G_gradient = -κₛ * (Ts - Tg) / (Δz / 2)
    Ts_implicit = Tg - G * Δz / (2 * κₛ)
    residual = Ts - Ts_implicit

    @test G_gradient ≈ G
    @test abs(residual) < sqrt(eps())

    @info "skin temperature: $Ts residual: $residual R: $(state.surface_net_radiation[1, 1]) Hₛ: $(state.sensible_heat_flux[1, 1]) Hₗ: $(state.latent_heat_flux[1, 1]) G: $G"

    return (
        skin_temperature = state.skin_temperature,
        residual = residual,
        surface_net_radiation = state.surface_net_radiation,
        latent_heat_flux = state.latent_heat_flux,
        sensible_heat_flux = state.sensible_heat_flux,
        ground_heat_flux = state.ground_heat_flux,
    )
end

@testset "Prescribed skin temperature" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N = 10))
    skin_temperature = PrescribedSkinTemperature(eltype(grid))
    seb = SurfaceEnergyBalance(Float64; skin_temperature)
    model = SurfaceEnergyModel(grid, surface_energy_balance = seb)
    state = StateVariables(model)
    @test hasproperty(state.inputs, :skin_temperature)
    set!(state.skin_temperature, 1.0)
    compute_auxiliary!(state, grid, skin_temperature)
    @test all(state.skin_temperature .≈ 1.0)
end

@testset "Implicit skin temperature" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N = 10))
    solver = Terrarium.default_skin_temperature_solver(eltype(grid))
    skin_temperature = ImplicitSkinTemperature(eltype(grid); κₛ = 0.5, solver)
    aerodynamics = Terrarium.ConstantAerodynamics(Cₕ = 2.0e-3)
    atmosphere = Terrarium.PrescribedAtmosphere(Float64; aerodynamics = aerodynamics)
    seb = SurfaceEnergyBalance(Float64; skin_temperature)
    model = SurfaceEnergyModel(grid, surface_energy_balance = seb, atmosphere = atmosphere)

    # Sunny and dry in Bergen Norway (Figure 5.11, Shuttleworth 2012)
    state = StateVariables(model)
    results = test_skin_temperature_solve!(
        state, grid, seb, model;
        surface_shortwave_down = 600.0,
        surface_longwave_down = 300.0,
        air_pressure = 101_325,
        air_temperature = 10.0,
        relative_humidity = 0.75,
        ground_temperature = 13.0,
        windspeed = 5.0,
    )
    @test all(results.sensible_heat_flux .> 0)
    @test all(results.latent_heat_flux .> 0)
    @test all(results.ground_heat_flux .< 0)
    @test all(results.surface_net_radiation .< 0)

    # Sunny and humid
    state = StateVariables(model)
    results = test_skin_temperature_solve!(
        state, grid, seb, model;
        surface_shortwave_down = 600.0,
        surface_longwave_down = 300.0,
        air_pressure = 101_325,
        air_temperature = 10.0,
        relative_humidity = 0.97,
        ground_temperature = 13.0,
        windspeed = 5.0,
    )
    @test all(results.sensible_heat_flux .> 0)
    @test all(results.latent_heat_flux .> 0)
    @test all(results.ground_heat_flux .< 0)
    @test all(results.surface_net_radiation .< 0)

    # Cloudy and dry
    state = StateVariables(model)
    results = test_skin_temperature_solve!(
        state, grid, seb, model;
        surface_shortwave_down = 150.0,
        surface_longwave_down = 200.0,
        air_pressure = 101_325,
        air_temperature = 5.0,
        relative_humidity = 0.6,
        ground_temperature = 10.0,
        windspeed = 10.0
    )
    @test all(results.sensible_heat_flux .< 0)
    @test all(results.latent_heat_flux .> 0)
    @test all(results.ground_heat_flux .> 0)
    @test all(results.surface_net_radiation .> 0)

    # Cloudy and humid
    state = StateVariables(model)
    results = test_skin_temperature_solve!(
        state, grid, seb, model;
        surface_shortwave_down = 150.0,
        surface_longwave_down = 200.0,
        air_pressure = 101_325,
        air_temperature = 5.0,
        relative_humidity = 0.9,
        ground_temperature = 10.0,
        windspeed = 10.0
    )
    @test all(results.sensible_heat_flux .> 0)
    @test all(results.latent_heat_flux .> 0)
    @test all(results.ground_heat_flux .> 0)
    @test all(results.surface_net_radiation .> 0)
end
