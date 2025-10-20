struct SurfaceEnergyBalanceModel{
    NF,
    GridType<:AbstractLandGrid{NF},
    Albedo<:AbstractAlbedo,
    RadiativeFluxes<:AbstractRadiativeFluxes,
    TurbulentFluxes<:AbstractTurbulentFluxes,
    Initializer<:AbstractInitializer,
    TimeStepper<:AbstractTimeStepper,
} <: AbstractLandModel{NF, GridType, TimeStepper}
    "Spatial grid"
    grid::GridType

    "Scheme for parameterizing surface albedo"
    albedo::Albedo

    "Scheme for determining the net radiation budget"
    radiative_fluxes::RadiativeFluxes

    "Scheme for computing turbulent (sensible and latent) heat fluxes"
    turbulent_fluxes::TurbulentFluxes

    "Scheme for determining skin temperature and ground heat flux"
    skin_temperature::SkinTemperature

    "State variable initializer"
    initializer::Initializer
    
    "Time stepping scheme"
    time_stepping::TimeStepper
end


