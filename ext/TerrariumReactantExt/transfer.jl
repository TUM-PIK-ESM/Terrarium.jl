# Device initialization.
#
# A model whose grid lives on `ReactantState` is initialized on the device: the state
# `Field`s are allocated directly on the device grid and `initialize!` populates them in place.
# This works because the eager KernelAbstractions launches invoked by `initialize!` (halo fills,
# closure/auxiliary kernels) run on the Reactant backend, and Oceananigans' `set_to_function!`
# for a device field internally detours through the CPU. Only `run!`/`timestep!` are traced and
# compiled (see integrator.jl).
#
# The sole device-specific ingredient is the clock: a traced `ConcreteRNumber`-backed
# `Clock(field_grid)` so that time advances *inside* the compiled step. With that clock in hand
# the rest is the generic, architecture-agnostic `Terrarium.initialize`, reached with `@invoke`
# (a direct call would re-dispatch to this same, more-specific method and recurse).

function Terrarium.initialize(
        model::ReactantModel{NF},
        params = nothing;
        clock = nothing,
        inputs = Terrarium.InputSources(NF),
        boundary_conditions = (;),
        initializers = (;),
        fields = (;),
    ) where {NF}
    # Traced device clock (unless the caller supplied one).
    device_clock = isnothing(clock) ?
        Oceananigans.TimeSteppers.Clock(get_field_grid(get_grid(model))) : clock
    return @invoke Terrarium.initialize(
        model::AbstractModel, params;
        clock = device_clock, inputs, boundary_conditions, initializers, fields
    )
end
