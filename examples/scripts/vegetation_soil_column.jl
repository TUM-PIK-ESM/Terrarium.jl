using Terrarium
using CUDA

# Define a simple grid with 1 column
grid = ColumnGrid(GPU(), ExponentialSpacing(Δz_max=1.0, N=30))
# Set up Richards model for soil hydrology
swrc = VanGenuchten(α=2.0, n=2.0)
hydraulic_properties = ConstantHydraulics(eltype(grid), cond_unsat=UnsatKVanGenuchten(eltype(grid); swrc))
hydrology = SoilHydrology(eltype(grid), RichardsEq(); hydraulic_properties)
# Variably saturated with water table at roughly 5 m depth
initializer = FieldInitializers(
    saturation_water_ice = (x, z) -> min(1, 0.5 - 0.1*z)
)
# Construct submodels
soil = SoilModel(grid; hydrology, initializer)
vegetation = VegetationModel(grid)
# Construct coupled model
coupled_model = VegetationSoilModel(grid; soil, vegetation)
# TODO: this is currently slow
integrator = @time initialize(coupled_model, ForwardEuler());
Δt = 60.0
timestep!(integrator, Δt) # TODO: currently produces NaNs for some variables
