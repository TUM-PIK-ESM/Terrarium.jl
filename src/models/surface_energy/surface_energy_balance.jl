struct SurfaceEnergyBalanceModel{
    NF,
    GridType<:AbstractLandGrid,
    RadiativeBudget<:AbstractRadiativeBudget,
    TurbulentFluxes<:AbstractTurbulentFluxes,
    Initializer<:AbstractInitializer,
    TimeStepper<:AbstractTimeStepper,
} <: AbstractLandModel{NF, GridType, TimeStepper}
    "Spatial grid"
    grid::GridType

    "Scheme for determining the net radiation budget"
    radiative_budget::RadiativeBudget

    "Scheme for computing turbulent (sensible and latent) heat fluxes"
    turbulent_fluxes::TurbulentFluxes

    "Scheme for determining skin temperature and ground heat flux"
    skin_temperature::SkinTemperature

    "State variable initializer"
    initializer::Initializer
    
    "Time stepping scheme"
    time_stepping::TimeStepper
end