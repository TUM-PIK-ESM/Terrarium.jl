using Terrarium
using CUDA

arch = CPU()
# Define a simple grid with 1 column
grid = ColumnGrid(arch, ExponentialSpacing(Δz_max=1.0, N=30))
# Set up Richards model for soil hydrology
swrc = VanGenuchten(α=2.0, n=2.0)
hydraulic_properties = ConstantSoilHydraulics(eltype(grid); swrc, unsat_hydraulic_cond=UnsatKVanGenuchten(eltype(grid)))
hydrology = SoilHydrology(eltype(grid), RichardsEq(); hydraulic_properties)
soil = SoilEnergyWaterCarbon(eltype(grid); hydrology)
vegetation = VegetationCarbon(eltype(grid))
# Variably saturated with water table at roughly 5 m depth
initializer = FieldInitializers(
    saturation_water_ice = (x, z) -> min(1, 0.5 - 0.1*z),
    carbon_vegetation = 0.1,
)
# Construct coupled model
vegsoil = VegetationSoilModel(grid; soil, vegetation, initializer)
# TODO: this is currently slow
integrator = @time initialize(vegsoil, ForwardEuler());
Δt = 60.0
@time timestep!(integrator, Δt)
