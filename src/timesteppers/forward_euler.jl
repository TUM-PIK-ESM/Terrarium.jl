"""
    $TYPEDEF

Simple forward Euler time stepping scheme.
"""
@kwdef struct ForwardEuler{NF} <: AbstractTimeStepper{NF}
    "Initial timestep size in seconds"
    Δt::NF = 300.0
end

default_dt(euler::ForwardEuler) = euler.Δt

is_adaptive(euler::ForwardEuler) = false

is_initialized(euler::ForwardEuler) = true

function timestep!(driver::ModelDriver, timestepper::ForwardEuler, Δt)
    # Euler step
    explicit_step!(driver.state, get_grid(driver.model), timestepper, Δt)
    # Apply closure relations
    closure!(driver.state, driver.model)
    # Update clock
    tick!(driver.state.clock, Δt)
end
