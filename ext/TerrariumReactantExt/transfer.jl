# Build-on-CPU, transfer-to-device initialization.
#
# Eager KernelAbstractions launches (e.g. the halo fills and closure kernels invoked by
# `initialize!`) cannot run on the Reactant backend outside a compiled region. So when the model's
# grid lives on `ReactantState`, we build a CPU twin, initialize it eagerly on the CPU, then move
# the initialized state to the device. Kernel launches then only occur inside the compiled
# `timestep!` (see integrator.jl), where they trace normally.

# Rebuild a model with a new grid, re-inferring type parameters. `setproperties` keeps the original
# concrete parameters and would try to `convert` the device grid to the CPU grid type.
function rebuild_model(model::AbstractModel, newgrid)
    ModelWrapper = Base.typename(typeof(model)).wrapper
    return ModelWrapper(; getproperties(model)..., grid = newgrid)
end

function Terrarium.initialize(
        model::ReactantModel{NF},
        params = nothing;
        clock = nothing,
        inputs = Terrarium.InputSources(NF),
        boundary_conditions = (;),
        initializers = (;),
        fields = (;),
    ) where {NF}
    device_grid = get_grid(model)
    # 1. CPU twin: same processes/params, CPU grid.
    cpu_model = rebuild_model(model, on_architecture(CPU(), device_grid))
    cpu_clock = isnothing(clock) ? Terrarium.Clock(time = zero(NF)) : on_architecture(CPU(), clock)
    cpu_integrator = Terrarium.initialize(cpu_model, params;
        clock = cpu_clock, inputs, boundary_conditions, initializers, fields)
    # 2. Transfer the fully-initialized CPU integrator to the device.
    return to_device(cpu_integrator, model, boundary_conditions)
end

# Assemble a device integrator: fresh device state (allocated on the device grid, no eager
# initializer kernels) with the CPU state's data copied in, and a traced device clock.
function to_device(cpu_integrator, device_model::ReactantModel{NF}, boundary_conditions) where {NF}
    device_grid = get_grid(device_model)
    device_field_grid = get_field_grid(device_grid)
    device_clock = Oceananigans.TimeSteppers.Clock(device_field_grid)   # traced ConcreteRNumber clock
    inputs = cpu_integrator.inputs
    device_state = StateVariables(device_model;
        clock = device_clock,
        input_variables = Terrarium.variables(inputs),
        boundary_conditions,
    )
    copy_state_data!(device_state, cpu_integrator.state)
    return ModelIntegrator(device_clock, device_model, inputs, device_state, cpu_integrator.initializers)
end

# Copy field data (and clock) from a source state onto a destination state, pairing by name and
# recursing into namespaces. Only mutable `Field`s are copied directly; views (e.g. the
# ground-temperature view into `temperature`) track their parent and update automatically.
function copy_state_data!(dst, src)
    for group in (:prognostic, :tendencies, :auxiliary, :inputs)
        _copy_group!(getfield(dst, group), getfield(src, group))
    end
    for nsname in namespace_names(dst)
        copy_state_data!(getproperty(getfield(dst, :namespaces), nsname),
                         getproperty(getfield(src, :namespaces), nsname))
    end
    _copy_clock!(getfield(dst, :clock), getfield(src, :clock))
    return dst
end

function _copy_group!(dst_nt::NamedTuple, src_nt::NamedTuple)
    for name in keys(dst_nt)
        dst = dst_nt[name]
        dst isa Oceananigans.Fields.AbstractField || continue
        copyto!(interior(dst), Array(interior(src_nt[name])))
    end
    return dst_nt
end

_scalar(x::Number) = x
_scalar(x) = Reactant.to_number(x)

function _copy_clock!(dst_clock, src_clock)
    dst_clock.time = convert(typeof(_scalar(dst_clock.time)), _scalar(src_clock.time))
    dst_clock.iteration = convert(typeof(_scalar(dst_clock.iteration)), _scalar(src_clock.iteration))
    return dst_clock
end
