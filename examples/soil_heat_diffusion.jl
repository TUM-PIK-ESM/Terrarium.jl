# # Soil heat diffusion with periodic surface forcing
#
# This example demonstrates how to simulate 1D soil heat conduction
# in a single vertical column using Terrarium.jl.
#
# We verify the numerical solution against the analytical solution
# of the heat equation with sinusoidal surface forcing:
#
#     ∂T/∂t = α ∂²T/∂z²
#
# where the analytical solution is
#
#     T(z,t) = T₀ + A exp(-z√(π/(αP))) sin(2πt/P − z√(π/(αP)))
#
# This example validates:
#
# - exponential damping with depth
# - phase lag with depth
# - agreement between analytic and numerical solution

using Terrarium
using CairoMakie

# Physical parameters

T₀ = 2.0               # mean surface temperature
A = 1.0                # forcing amplitude
P = 24 * 3600          # forcing period (seconds)
k = 2.0                # thermal conductivity
c = 1.0e6              # volumetric heat capacity
α = k / c              # thermal diffusivity
w = 2*pi/P              # angular frequency of surface forcing
d = sqrt(2*α/w)         # thermal damping depth

# Known Analytical solution

function heat_conduction_solution(T₀, A, P, α)
    w = 2*pi/P              # angular frequency of surface forcing
    d = sqrt(2*α/w)         # thermal damping depth
    T(z, t) = T₀ +
              A *
              exp(-z / d) *
              sin(w * t - (z/d))
    return T
end

T_sol = heat_conduction_solution(T₀, A, P, α)

# Model setup

grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(Δz_min = 0.05, Δz_max = 100.0, N = 100,))

# Configure soil properties for pure heat diffusion
# set carbon content to zero so the soil has only a mineral constituent
biogeochem = ConstantSoilCarbonDensity(ρ_soc = 0.0)
# set porosity to zero to remove influence of pore space;
soil_porosity = ConstantSoilPorosity(mineral_porosity = 0.0)
strat = HomogeneousStratigraphy(Float64; porosity = soil_porosity)
# set thermal properties
thermal_properties = SoilThermalProperties(
    eltype(grid);
    conductivities = SoilThermalConductivities(mineral = k),
    heat_capacities = SoilHeatCapacities(mineral = c),
)
energy = SoilEnergyBalance(eltype(grid); thermal_properties)
soil = SoilEnergyWaterCarbon(eltype(grid); energy, strat, biogeochem)
model = SoilModel(grid; soil)

# Apply periodic surface temperature forcing

upper_bc(z, t) = T₀ + A * sin(2π * t / P)
bcs = PrescribedSurfaceTemperature(:Tsurf, upper_bc)

# Temp IC
initializers = (temperature = (x, z) -> T_sol(-z, 0.0),)
integrator = initialize(model, ForwardEuler(); initializers, boundary_conditions=bcs)

# Run simulation

Δt = 60.0
Ts_buffer = [deepcopy(integrator.state.temperature)]
times = [0.0]

# run for 1 hour, saving every time step
while current_time(integrator) < 2P
    timestep!(integrator, Δt)
    push!( Ts_buffer, deepcopy(integrator.state.temperature))
    push!(times,current_time(integrator))
end

# Extract numerical solution
z = znodes(integrator.state.temperature)

T_numeric = reduce(hcat, Ts_buffer)[1, :, :]

# Compute analytical solution on same grid

T_target = T_sol.(reshape(-z, 1, :), reshape(times, :, 1))

# Compute error

relative_error = abs.((T_numeric .- T_target) ./ T_target)

max_error = maximum(relative_error)

println("Maximum relative error = ", max_error)


# Verify solver accuracy

@assert max_error < 0.1


# Plot comparison at final timestep

fig = Figure()
empty!(ax)
ax = Axis(fig[1, 1], xlabel = "Temperature", ylabel = "Depth",)
ylims!(ax, -5d, 0.0)
lines!(ax, T_numeric[end, :], z, label = "Numerical")
lines!(ax, T_target[end, :], z, linestyle = :dash, label = "Analytical")
axislegend(ax)
fig