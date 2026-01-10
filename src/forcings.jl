
const OceananigansForcing{P} = Union{DiscreteForcing{P}, ContinuousForcing{LX, LY, LZ, P}} where {LX, LY, LZ}

"""
    $TYPEDEF

Container type that wraps a `Dict` of forcings, which may be `AbstractProcess`es (`AbstractProcess`es that implement
the [`forcing`](@ref) method) or `Oceananigans` forcing functions, i.e. `ContinuousForcing` or `DiscreteForcing`.
"""
struct Forcings
    forcings::Dict{Symbol, Union{AbstractProcess, OceananigansForcing}}
end

const ForcingType = Union{AbstractProcess, OceananigansForcing, Forcings}

function Forcings(; kwargs...)
    return Forcings(Dict(kwargs...))
end

"""
    forcing(i, j, k, state, grid, forcing::AbstractProcess, target::AbstractProcess, args...)

Compute a forcing tendency contribution of `proc` to the `target` process at volume `i, j, k`
on `grid` given the current `state`. Typically `proc` represents a source/sink that acts on
a prognostic state variable defined by `target`. Default implementation returns `zero(eltype(grid))`.
"""
@inline function forcing(i, j, k, state, grid, forcing::AbstractProcess, target::AbstractProcess, args...)
    return zero(eltype(grid))
end

"""
    forcing(i, j, k, state, grid, target::AbstractProcess, args...)

Convenience dispatch of `forcing` that assumes `target` to be a process with a property `forcing`
corresponding to an `AbstractForcing` or Oceananigans forcing function.
"""
@inline function forcing(i, j, k, state, grid, target::AbstractProcess, args...)
    return forcing(i, j, k, state, grid, target.forcing, target, args...)
end

"""
    forcing(i, j, k, state, grid, fs::Forcings, target::AbstractProcess, args...)

Invoke and sum the forcing tendencies for all components defined in `forcings`.
"""
@inline function forcing(i, j, k, state, grid, fs::Forcings, target::AbstractProcess, args...)
    return sum(map(f -> forcing(i, j, k, state, grid, f, target, args...), values(fs.forcings)), init=zero(eltype(grid)))
end

"""
    forcing(i, j, k, state, grid, forcing::OceananigansForcing, target::AbstractProcess, args...)

Return the value computed by the given `Oceananigans` forcing function, which should be an instance of
either `DiscreteForcing` or `ContinuousForcing`. Note that `target` and additional `args` are only included
for interface consistency and are not passed through to `forcing`.
"""
@inline function forcing(i, j, k, state, grid, forcing::OceananigansForcing, target::AbstractProcess, args...)
    return forcing(i, j, k, get_field_grid(grid), state.clock, state)
end
