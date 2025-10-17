"""
    RichardsEq{PS} <: AbstractSoilWaterFluxOperator

Operator for soil hydrology corresponding to the Richardson-Richards equation for variably saturated
flow in porous media.
"""
@kwdef struct RichardsEq{PS} <: AbstractSoilWaterFluxOperator
    "Closure relation for mapping between water potential (hydraulic head) and saturation"
    saturation_closure::PS = PressureSaturationClosure()
end

get_closure(op::RichardsEq) = op.saturation_closure

variables(hydrology::SoilHydrology{NF, <:RichardsEq}) where {NF} = (
    prognostic(:pressure_head, XYZ(), get_closure(hydrology.operator)),
    auxiliary(:hydraulic_conductivity, XYZ(), units=u"m/s", desc="Hydraulic conductivity of soil volumes [m/s]"),
    auxiliary(:water_table, XY(), units=u"m", desc="Elevation of the water table [m]"),
)

function initialize!(state, model, ::SoilHydrology{NF, <:RichardsEq}) where {NF}
    invclosure!(state, model, state.closures.pressure_head)
    return nothing
end

function compute_auxiliary!(state, model, hydrology::SoilHydrology{NF, <:RichardsEq}) where {NF}
    grid = get_grid(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)
    launch!(grid, :xyz, compute_hydraulic_conductivity!, state, grid, hydrology, strat, bgc)
    return nothing
end

@kernel function compute_hydraulic_conductivity!(
    state,
    grid,
    hydrology::SoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
)
    i, j, k = @index(Global, NTuple)
    soil = soil_composition((i, j, k), state, strat, hydrology, bgc)
    state.hydraulic_conductivity[i, j, k] = hydraulic_conductivity(hydrology.hydraulic_properties, soil)
    # Manually set hydraulic conductivity at the boundaries to the conductivity of the cell centers at the boundaries;
    # note that this is an assumption that could potentially be relaxed in the future.
    # TODO: maybe this should be done with boundary conditions?
    field_grid = get_field_grid(grid) 
    if k == 1
        state.hydraulic_conductivity[i, j, 0] = state.hydraulic_conductivity[i, j, 1]
    elseif k == field_grid.Nz
        state.hydraulic_conductivity[i, j, field_grid.Nz+1] = state.hydraulic_conductivity[i, j, field_grid.Nz]
    end
end

function compute_tendencies!(state, model, hydrology::SoilHydrology{NF, <:RichardsEq}) where {NF}
    grid = get_grid(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)
    launch!(grid, :xyz, compute_saturation_tendency!, state, grid, hydrology, strat, bgc)
    return nothing
end

@kernel function compute_saturation_tendency!(
    state,
    grid,
    hydrology::SoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
)
    i, j, k = @index(Global, NTuple)
    state.tendencies.saturation_water_ice[i, j, k] += saturation_tendency((i, j, k), grid, state, hydrology, strat, bgc)
end

@inline function saturation_tendency(
    idx, grid, state,
    hydrology::SoilHydrology{NF, <:RichardsEq},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
) where {NF}
    i, j, k = idx
    # Operators require the underlying Oceananigans grid
    field_grid = get_field_grid(grid)

    # Compute divergence of water fluxes
    # ∂θ∂t = ∇⋅K(θ)(∇ψ + 1) = ∇⋅K(θ)∇ψ
    ∂θ∂t = ∂zᵃᵃᶜ(i, j, k, field_grid, darcy_flux, state)

    # Get porosity
    por = porosity(idx, state, hydrology, strat, bgc)

    # Rescale volumetric flux to saturation flux
    return -∂θ∂t / por
end

@inline function darcy_flux(i, j, k, grid, state)
    # Get pressure head
    ψ = state.pressure_head
    # Get hydraulic conductivity
    K = state.hydraulic_conductivity
    # Darcy's law: q = -K ∂ψ/∂z
    # TODO: also account for effect of temperature on matric potential
    ∇ψ = ∂zᵃᵃᶠ(i, j, k, grid, ψ)
    ## Take minimum of hydraulic conductivities in the direction of flow
    Kₖ = (∇ψ < 0)*min(K[i, j, k-1], K[i, j, k]) +
         (∇ψ >= 0)*min(K[i, j, k], K[i, j, k+1])
    ## Note that elevation and hydrostatic pressure are assumed to be accounted
    ## for in the computation of ψ, so we do need any additional terms.
    q = -Kₖ*∇ψ
    return q
end

# Matric potential <--> saturation closure relation

@kwdef struct PressureSaturationClosure <: AbstractClosureRelation end

closurevar(::PressureSaturationClosure) = auxiliary(
    :saturation_water_ice,
    XYZ();
    domain=UnitInterval(),
    desc="Saturation level of water and ice in the pore space",
)

function closure!(state, model::AbstractSoilModel, ::PressureSaturationClosure)
    grid = get_grid(model)
    hydrology = get_soil_hydrology(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)
    # determine saturation from pressure
    launch!(grid, :xyz, pressure_to_saturation!, state, hydrology, strat, bgc)
    # diagnose water table elevation
    launch!(grid, :xyz, compute_water_table!, state, grid, hydrology)
    return nothing
end

function invclosure!(state, model::AbstractSoilModel, ::PressureSaturationClosure)
    grid = get_grid(model)
    hydrology = get_soil_hydrology(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)
    # diagnose water table elevation from current saturation state
    launch!(grid, :xyz, compute_water_table!, state, grid, hydrology)
    # determine pressure head from saturation
    launch!(grid, :xyz, saturation_to_pressure!, state, hydrology, strat, bgc)
    return nothing
end

@kernel function pressure_to_saturation!(
    state,
    hydrology::SoilHydrology{NF, <:RichardsEq},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
) where {NF}
    idx = @index(Global, NTuple)
    pressure_to_saturation!(idx, state, hydrology, strat, bgc)
end

@inline function pressure_to_saturation!(
    idx, state,
    hydrology::SoilHydrology{NF, <:RichardsEq},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
) where {NF}
    i, j, k = idx
    ψ = state.pressure_head[i, j, k] # assumed given
    swrc = get_swrc(hydrology)
    por = porosity(idx, state, hydrology, strat, bgc)
    vol_water_ice_content = swrc(ψ; θsat=por)
    state.saturation_water_ice[i, j, k] = vol_water_ice_content / por
    return nothing
end

@kernel function saturation_to_pressure!(
    state,
    hydrology::SoilHydrology{NF, <:RichardsEq},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
) where {NF}
    idx = @index(Global, NTuple)
    saturation_to_pressure!(idx, state, hydrology, strat, bgc)
end

@inline function saturation_to_pressure!(
    idx, state,
    hydrology::SoilHydrology{NF, <:RichardsEq},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
) where {NF}
    i, j, k = idx
    sat = state.saturation_water_ice[i, j, k] # assumed given
    # TODO: To avoid invalid saturation values, we will just clip it to the
    # valid range here for now... but this of course results in mass balance violations.
    # sat = state.saturation_water_ice[i, j, k] = clamp(sat, zero(sat), one(sat))
    # get inverse of SWRC
    inv_swrc = inv(get_swrc(hydrology))
    por = porosity(idx, state, hydrology, strat, bgc)
    # compute matric pressure head
    ψm = inv_swrc(sat*por; θsat=por)
    # compute elevation pressure head
    zs = znodes(state.saturation_water_ice)
    ψz = zs[k]
    # compute hydrostatic pressure head assuming impermeable lower boundary
    z₀ = state.water_table[i, j, 1]
    ψh = max(0, z₀ - zs[k])
    # compute total pressure head as sum of ψh + ψm + ψz
    # note that ψh and ψz will cancel out in the saturated zone
    state.pressure_head[i, j, k] = ψh + ψm + ψz
    return nothing
end
