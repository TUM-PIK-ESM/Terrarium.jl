# Return true if debug mode is enabled, false otherwise
@inline debug_mode() = DEBUG[]

"""
    $SIGNATURES

Check whether the given `field` has any `NaN` or `Inf` values and raise an error if `NaN`s are detected.
"""
checkfinite!(field::AbstractField, name = nothing) = any(!isfinite, parent(field)) && error("Found NaN/Inf values in Field $name: $field")
function checkfinite!(nt::NamedTuple)
    for key in keys(nt)
        checkfinite!(nt[key], key)
    end
    return nothing
end

"""
    $SIGNATURES

Provides a "hook" for handling debug calls from relevant callsites. Default implementations for
`Field` and `NamedTuple` (assumed to be of `Field`s) simply forward to [`checkfinite!`](@ref).
"""
@inline debughook!(args...) = nothing
@inline debughook!(field::AbstractField) = checkfinite!(field)
@inline debughook!(nt::NamedTuple) = checkfinite!(nt)

"""
    $SIGNATURES

Utility method that forwards `args` to `debughook!` *if and only if debug mode is enabled*. Debug mode is set by
the global variable `DEBUG` which can be toggled by the user facing API [`debug!`](@ref).
"""
@inline function debugsite!(args...)
    if debug_mode()
        debughook!(args...)
    end
    return nothing
end

function kernel_operation(func, state, grid, args...; location = (Center, Center, Nothing), with_clock = false)
    fields = get_fields(state, args...)
    fgrid = get_field_grid(grid)
    if with_clock
        return KernelFunctionOperation{location...}(func, fgrid, state.clock, fields, args...)
    else
        return KernelFunctionOperation{location...}(func, fgrid, fields, args...)
    end
end

kernel_operation2D(func, state, grid, args...; X = Center, Y = Center, with_clock = false) = kernel_operation(func, state, grid, args...; location = (X, Y, Nothing), with_clock)
kernel_operation3D(func, state, grid, args...; X = Center, Y = Center, Z = Center, with_clock = false) = kernel_operation(func, state, grid, args...; location = (X, Y, Z), with_clock)
