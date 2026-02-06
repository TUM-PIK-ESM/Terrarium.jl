using Terrarium
using Test

using Oceananigans.BoundaryConditions: BoundaryCondition, Flux

@testset "VegetationSoilModel" begin
    grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max=1.0, N=50))
    swrc = VanGenuchten(α=2.0, n=2.0)
    hydraulic_properties = ConstantSoilHydraulics(eltype(grid); swrc, unsat_hydraulic_cond=UnsatKVanGenuchten(eltype(grid)))
    hydrology = SoilHydrology(eltype(grid), RichardsEq(); hydraulic_properties)
    soil = SoilEnergyWaterCarbon(eltype(grid); hydrology)
    # Variably saturated with water table
    initializer = FieldInitializers(
        temperature = (x, z) -> 1.0 - 0.02 * z,
        saturation_water_ice = (x, z) -> min(1, 0.8 - 0.05*z)
    )
    vegetation = VegetationCarbon(eltype(grid))
    vegsoil = VegetationSoilModel(grid; soil, vegetation)
    integrator = initialize(vegsoil, ForwardEuler())
    # Check that infiltration is correctly coupled to soil hydrology
    set!(integrator.state.infiltration, 1e-8)
    sat_top_bc = integrator.state.saturation_water_ice.boundary_conditions.top
    @test isa(sat_top_bc, BoundaryCondition{<:Flux})
    @test all(Field(sat_top_bc.condition) .≈ -1e-8) # note the negative sign
    # Check that ground heat flux is coupled to soil energy
    energy_top_bc = integrator.state.internal_energy.boundary_conditions.top
    @test isa(energy_top_bc, BoundaryCondition{<:Flux})
    @test energy_top_bc.condition == integrator.state.ground_heat_flux
    # Advance one timestep
    timestep!(integrator, 60.0)
    @test all(isfinite.(integrator.state.saturation_water_ice))
    @test all(isfinite.(integrator.state.internal_energy))
    # TODO: also check ET and veg processes once they are working...
end