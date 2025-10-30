using Terrarium
using Test

using FreezeCurves
using Oceananigans

import Terrarium: SurfaceEvaporation, hydraulic_conductivity

@testset "Hydraulic properties (constant)" begin
    # For prescribed hyraulic properties, just check that the returned values
    # match what was set.
    hydraulic_props = ConstantHydraulics(
        Float64;
        cond_sat = 1e-6,
        porosity = 0.3,
        field_capacity = 0.1,
        wilting_point = 0.02,
    )
    @test saturated_hydraulic_conductivity(hydraulic_props) == hydraulic_props.cond_sat
    @test mineral_porosity(hydraulic_props) == hydraulic_props.porosity
    @test field_capacity(hydraulic_props) == hydraulic_props.field_capacity
    @test wilting_point(hydraulic_props) == hydraulic_props.wilting_point
end

@testset "Hydraulic properties (SURFEX)" begin
    # mineral porosity
    # TODO: should the hydraulic properties struct constructors also enforce parameter bounds?
    hydraulic_props = SoilHydraulicsSURFEX(Float64)
    por0 = mineral_porosity(hydraulic_props, SoilTexture(sand=0.0, silt=0.7, clay=0.3))
    @test por0 ≈ hydraulic_props.porosity
    for sand in 0.1:0.1:1.0
        silt = (1 - sand)*0.7
        clay = (1 - sand)*0.3
        por = mineral_porosity(hydraulic_props, SoilTexture(; sand, silt, clay))
        # test that increasing sand content decreases the mineral porosity;
        # we could also reproduce the calculation here, but that seems a bit redundant
        # and of course depends on the choice of parameters.
        @test 0 < por < por0
    end

    # check that wilting point is equal to zero when there is no clay
    wp0 = wilting_point(hydraulic_props, SoilTexture(sand=0.5, silt=0.5, clay=0.0))
    @test iszero(wp0)
    for clay in 0.1:0.1:1.0
        sand = (1-clay)*0.7
        silt = (1-clay)*0.3
        wp = wilting_point(hydraulic_props, SoilTexture(; sand, silt, clay))
        @test 0 < wp < 1
    end

    # check that field capacity is equal to zero when there is no clay
    fc0 = field_capacity(hydraulic_props, SoilTexture(sand=0.5, silt=0.5, clay=0.0))
    @test iszero(fc0)
    for clay in 0.1:0.1:1.0
        sand = (1-clay)*0.7
        silt = (1-clay)*0.3
        fc = field_capacity(hydraulic_props, SoilTexture(; sand, silt, clay))
        @test 0 < fc < 1
    end
end

@testset "Unsaturated hydraulic conductivity (linear)" begin
    hydraulics = ConstantHydraulics(Float64; cond_unsat=UnsatKLinear(Float32))

    # saturated case
    soil = SoilComposition()
    K = hydraulic_conductivity(hydraulics, soil)
    @test K ≈ hydraulics.cond_sat

    # unsaturated
    soil = SoilComposition(saturation=0.5);
    K = hydraulic_conductivity(hydraulics, soil)
    @test 0 < K < hydraulics.cond_sat

    # dry
    soil = SoilComposition(saturation=0.0);
    K = hydraulic_conductivity(hydraulics, soil)
    @test iszero(K)

    # frozen
    soil = SoilComposition(liquid=0.0);
    K = hydraulic_conductivity(hydraulics, soil)
    @test iszero(K)
end

@testset "Unsaturated hydraulic conductivity (Van Genuchten)" begin
    hydraulics = ConstantHydraulics(Float64; cond_unsat=UnsatKVanGenuchten(Float64))

    # saturated case
    soil = SoilComposition()
    K = hydraulic_conductivity(hydraulics, soil)
    @test K ≈ hydraulics.cond_sat

    # unsaturated
    soil = SoilComposition(saturation=0.5);
    K = hydraulic_conductivity(hydraulics, soil)
    @test 0 < K < hydraulics.cond_sat

    # dry
    soil = SoilComposition(saturation=0.0);
    K = hydraulic_conductivity(hydraulics, soil)
    @test iszero(K)

    # frozen
    soil = SoilComposition(liquid=0.0);
    K = hydraulic_conductivity(hydraulics, soil)
    @test iszero(K)
end

@testset "SoilHydrology: Richardson-Richards' equation" begin
    grid = ColumnGrid(UniformSpacing(Δz=0.1, N=100))
    swrc = VanGenuchten(α=2.0, n=2.0)
    hydraulic_properties = ConstantHydraulics(Float64; cond_unsat=UnsatKVanGenuchten(Float64; swrc))
    hydrology = SoilHydrology(eltype(grid), RichardsEq(); hydraulic_properties)

    # Fully saturated, steady state
    initializer = FieldInitializers(
        saturation_water_ice = (x, z) -> 1.0
    )
    model = SoilModel(; grid, hydrology, initializer)
    state = initialize(model)
    # check that initial water table depth is correctly calculated from initial condition
    @test all(iszero.(state.state.water_table))
    # also check that pressure head is zero everywhere
    @test all(iszero.(state.state.pressure_head))
    compute_auxiliary!(state.state, model)
    # check that all hydraulic conductivities are finite and equal to K_sat
    @test all(isfinite.(state.state.hydraulic_conductivity))
    @test all(state.state.hydraulic_conductivity .≈ hydraulic_properties.cond_sat)
    compute_tendencies!(state.state, model)
    # check that all tendencies are zero
    @test all(iszero.(state.state.tendencies.saturation_water_ice))
    # check timestep!
    timestep!(state)
    @test all(state.state.saturation_water_ice .≈ 1)

    # Variably saturated with water table
    initializer = FieldInitializers(
        saturation_water_ice = (x, z) -> min(1, 0.5 - 0.1*z)
    )
    model = SoilModel(; grid, hydrology, initializer)
    state = initialize(model)
    water_table = state.state.water_table
    hydraulic_cond = state.state.hydraulic_conductivity
    saturation = state.state.saturation_water_ice
    # check that initial water table depth is correctly calculated from initial condition
    @test all(water_table .≈ -5.0)
    # also check that pressure head is negative
    @test all(state.state.pressure_head .< 0)
    compute_auxiliary!(state.state, model)
    # check that all hydraulic conductivities are finite and positive
    @test all(isfinite.(hydraulic_cond))
    @test all(hydraulic_cond .> 0)
    compute_tendencies!(state.state, model)
    # check that all tendencies are finite
    @test all(isfinite.(state.state.tendencies.saturation_water_ice))
    # do timestep! and compute total water mass
    ∫sat₀ = Field(Integral(saturation, dims=3))
    compute!(∫sat₀)
    timestep!(state)
    ∫sat₁ = Field(Integral(saturation, dims=3))
    compute!(∫sat₁)
    # check saturation levels are all finite and valid
    @test all(isfinite.(saturation))
    @test all(0 .<= saturation .<= 1)
    # check mass conservation
    @test ∫sat₀[1,1,1] ≈ ∫sat₁[1,1,1]
    # run for one simulation hour and check that mass is still conserved
    run!(state, period=Hour(1), Δt=60.0)
    ∫sat₂ = Field(Integral(saturation, dims=3))
    compute!(∫sat₂)
    @test all(isfinite.(saturation))
    @test all(0 .<= saturation .<= 1)
    @test ∫sat₀[1,1,1] ≈ ∫sat₁[1,1,1] ≈ ∫sat₂[1,1,1]
end

@testset "Soil ET" begin
    Nz = 10
    grid = ColumnGrid(UniformSpacing(Δz=0.1, N=Nz))
    swrc = VanGenuchten(α=2.0, n=2.0)
    hydraulic_properties = ConstantHydraulics(Float64; cond_unsat=UnsatKVanGenuchten(Float64; swrc))
    evapotranspiration = SurfaceEvaporation()
    hydrology = SoilHydrology(eltype(grid), RichardsEq(); hydraulic_properties, evapotranspiration)
    # Variably saturated with water table
    initializer = FieldInitializers(
        temperature = 10.0, # positive soil temperature
        saturation_water_ice = 1.0 # fully saturated
    )
    model = SoilModel(; grid, hydrology, initializer)
    state = initialize(model)
    # check that forcing_ET is zero when no latent heat flux is supplied
    @test iszero(Terrarium.forcing_ET(1, 1, Nz, grid.grid, state.state, evapotranspiration, model.constants))
    # negative latent heat flux → condensation
    set!(state.state.latent_heat_flux, -100.0)
    ET_flux = Terrarium.forcing_ET(1, 1, Nz, grid.grid, state.state, evapotranspiration, model.constants)
    @test ET_flux > 0
    # positive latent heat flux → evaporation
    set!(state.state.latent_heat_flux, 100.0)
    ET_flux = Terrarium.forcing_ET(1, 1, Nz, grid.grid, state.state, evapotranspiration, model.constants)
    @test ET_flux < 0
    # check tendency calculation
    dθdt = Terrarium.volumetric_water_content_tendency((1, 1, Nz), grid, state.state, hydrology, model.constants)
    @test dθdt == ET_flux
    # take one timestep and check that water was evaporated
    dt = 60.0
    timestep!(state, dt)
    @test state.state.saturation_water_ice[1, 1, Nz] .≈ 1 + ET_flux*dt / hydraulic_properties.porosity
end
