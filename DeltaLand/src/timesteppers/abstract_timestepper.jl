# AbstractTimeStepper

abstract type AbstractTimeStepper{NF} end

"""
    get_dt(timestepper::AbstractTimeStepper)

Get the current timestep size for the time stepper.
"""
function get_dt end

"""
    is_adaptive(timestepper::AbstractTimeStepper)

Returns `true` if the given time stepper is adaptive, false otherwise.
"""
function is_adaptive end
