# Initialization

function initialize!(state, model::SoilModel)
    # run model initializer
    initialize!(state, model, model.initializer)
    # launch kernel for generic initailization routine
    grid = get_grid(model)
    launch!(
        grid,
        :xyz,
        initialize_kernel!,
        state,
        model.energy,
        model.hydrology,
        model.strat,
        model.biogeochem,
        model.constants
    )
end

@kernel function initialize_kernel!(
    state,
    energy::AbstractSoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants,
)
    idx = @index(Global, NTuple)
    # TODO: need a more comprehensive initialization scheme for all soil model components
    # Note that this assumes temperature has already been iniitialized!
    fc = get_freezecurve(hydrology)
    temperature_to_energy!(idx, state, fc, energy, hydrology, strat, bgc, constants)
end
