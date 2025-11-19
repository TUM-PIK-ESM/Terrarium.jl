struct SurfaceEnergyModel{
        NF,
        GridType <: AbstractLandGrid{NF},
        SEB <: AbstractSurfaceEnergyBalance,
        Atmosphere <: AbstractAtmosphere,
        Initializer <: AbstractInitializer,
    } <: AbstractSurfaceEnergyModel{NF, GridType}
    "Spatial grid"
    grid::GridType

    "Atmospheric inputs"
    atmosphere::Atmosphere

    "Surface energy balance scheme"
    surface_energy_balance::SEB

    "Physical constants"
    constants::PhysicalConstants

    "State variable initializer"
    initializer::Initializer
end

function SurfaceEnergyModel(
        grid::AbstractLandGrid{NF},
        surface_energy_balance::AbstractSurfaceEnergyBalance = SurfaceEnergyBalance(NF);
        atmosphere::AbstractAtmosphere = PrescribedAtmosphere(NF),
        constants::PhysicalConstants = PhysicalConstants(NF),
        initializer::AbstractInitializer = DefaultInitializer()
    ) where {NF}
    return SurfaceEnergyModel(grid, atmosphere, surface_energy_balance, constants, initializer)
end

variables(model::SurfaceEnergyModel) = tuplejoin(variables(model.atmosphere), variables(model.surface_energy_balance))

function compute_auxiliary!(state, model::SurfaceEnergyModel)
    compute_auxiliary!(state, model, model.atmosphere)
    return compute_auxiliary!(state, model, model.surface_energy_balance)
end

function compute_tendencies!(state, ::SurfaceEnergyModel)
    compute_tendencies!(state, model, model.atmosphere)
    return compute_tendencies!(state, model, model.surface_energy_balance)
end
