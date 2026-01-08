"""
    $TYPEDEF

Base type for `AbstractProcess`es that provide forcing tendencies to other `AbstractProcess`es.
"""
abstract type AbstractForcing <: AbstractProcess end

const OceananigansForcing{P} = Union{DiscreteForcing{P}, ContinuousForcing{LX, LY, LZ, P}} where {LX, LY, LZ}
const AnyForcing = Union{AbstractForcing, OceananigansForcing}

"""
    $TYPEDEF

Container type that wraps a `Dict` of forcings, which may be `AbstractForcing`s (`AbstractProcess`es that implement
the [`forcing`](@ref) method) or `Oceananigans` forcing functions, i.e. `ContinuousForcing` or `DiscreteForcing`.
"""
struct Forcings <: AbstractForcing
    forcings::Dict{Symbol, AnyForcing}
end

function Forcings(; kwargs...)
    return Forcings(Dict(kwargs...))
end

"""
    forcing(i, j, k, state, grid, forcing::AbstractForcing, target::AbstractProcess, args...)

Compute a forcing tendency contribution of `proc` to the `target` process at volume `i, j, k`
on `grid` given the current `state`. Typically `proc` represents a source/sink that acts on
a prognostic state variable defined by `target`. 
"""
function forcing end

"""
"""
@inline function forcing(i, j, k, state, grid, target::AbstractProcess, args...)
    return forcing(i, j, k, state, grid, f, target.forcings, target, args...)
end

"""
    forcing(i, j, k, state, grid, forcings::Forcings, target::AbstractProcess, args...)

Invoke and sum the forcing tendencies for all components defined in `forcings`.
"""
@inline function forcing(i, j, k, state, grid, forcings::Forcings, target::AbstractProcess, args...)
    return sum(map(f -> forcing(i, j, k, state, grid, f, target, args...), values(forcings.forcings)))
end

"""
    forcing(i, j, k, state, grid, forcing::OceananigansForcing, target::AbstractProcess, args...)

Return the value computed by the given `Oceananigans` forcing function, which should be an instance of
either `DiscreteForcing` or `ContinuousForcing`. Note that `target` and additional `args` are only included
for interface consistency and are not passed to `forcing`.
"""
@inline function forcing(i, j, k, state, grid, forcing::OceananigansForcing, target::AbstractProcess, args...)
    return forcing(i, j, k, get_field_grid(grid), state.clock, state)
end
