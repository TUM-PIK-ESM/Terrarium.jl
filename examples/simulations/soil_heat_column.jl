# # [Soil heat conduction in a 1D vertical column](@id soil_heat_column)
# This example shows how to set up a simple model of nonlinear heat conduction
# in a single vertical soil column (similar to the example shown in the README).
# The [`SoilEnergyBalance`](@ref) process in Terrarium solves the nonlinear form
# of the heat equation with phase change. This allows for the simulation of
# freeze/thaw dynamics in both seasonally and perennially frozen soils.

using Terrarium
## For plotting
import CairoMakie as Makie

# We start by creating a single column grid with 10 exponentially spaced soil layers:
grid = ColumnGrid(CPU(), Float32, ExponentialSpacing(N = 10))

# Next we specify an initializer suitable for [`SoilModel`](@ref). We will use
# a quasi-steady-state initialization for soil temperature (linear with depth)
# and a fully saturated state for soil water/ice content.
initializer = SoilInitializer(
    eltype(grid),
    energy = QuasiThermalSteadyState(eltype(grid), T₀ = -1.0),
    hydrology = ConstantSaturation(eltype(grid), sat = 1.0)
)

# We're now ready to create our model:
model = SoilModel(grid; initializer)

# Boundary conditions are imposed directly on the corresponding `Field`s during
# initialization. We here set a constant surface temperature of 1°C, while the lower
# boundary is left with its default zero flux condition. Note that the name assigned
# to the boundary condition `:T_ub` is arbitrary; you can set it to whatever you want.
boundary_conditions = PrescribedSurfaceTemperature(:T_ub, 1.0)

# Next we choose a timestepper (here just [`ForwardEuler`](@ref)) and initialize the model:
timestepper = ForwardEuler(eltype(grid))
integrator = initialize(model, timestepper; boundary_conditions)

# We can now try taking a single one timestep and `@time` it; note that the first evaluation
# will be slower due to compilation.
timestep!(integrator)
@time timestep!(integrator)

# We an also run the simulation forward for a set period of time:
@time run!(integrator, period = Day(3))

# Now let's extract the relevant state variables for inspection. The function
# [`interior`](@extref Oceananigans.Fields.interior) comes from Oceananigans
# and is used to extract the interior values of the spatial field on grid cell centers,
# excluding the halo regions used for representing boundary conditions.
T = interior(integrator.state.temperature)[1, 1, :]
f = interior(integrator.state.liquid_water_fraction)[1, 1, :]

# Finally, we plot the temperature and liquid fraction profiles. With
# [`znodes`](@extref Oceananigans.Grids.znodes), the positions of the interior
# nodes in the $z$-direction can be extracted.
zs = znodes(integrator.state.temperature)
let fig = Makie.Figure()
    ax1 = Makie.Axis(fig[1, 1], ylabel = "Depth / m", xlabel = "Temperature / °C")
    ax2 = Makie.Axis(fig[1, 2], xlabel = "Liquid fraction")

    Makie.scatterlines!(ax1, T, zs)
    Makie.scatterlines!(ax2, f, zs)
    fig
end
