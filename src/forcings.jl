const OceananigansForcing{P} = Union{DiscreteForcing{P}, ContinuousForcing{LX, LY, LZ, P}} where {LX, LY, LZ}
const AnyForcing = Union{AbstractProcess, OceananigansForcing}

"""
    $TYPEDEF

Container type that wraps a `Dict` of forcing types, which may be `AbstractProcess`es that implement
the [`forcing`](@ref) method or `Oceananigans` forcing types, i.e. `ContinuousForcing` or `DiscreteForcing`.
"""
struct Forcings
    forcings::Dict{Symbol, AnyForcing}
end

function Forcings(; kwargs...)
    return Forcings(Dict(kwargs...))
end

"""
    forcing(i, j, k, state, grid, proc::AbstractProcess, target::AbstractProcess, args...)

Compute a forcing tendency contribution of `proc` to the `target` process at volume `i, j, k`
on `grid` given the current `state`. Typically `proc` represents a source/sink that acts on
a prognostic state variable defined by `target`. 
"""
function forcing end

"""
    forcing(i, j, k, state, grid, forcings::Forcings, target::AbstractProcess, args...)

Invoke and sum the forcing tendencies for all components defined in `forcings`.
"""
@inline function forcing(i, j, k, state, grid, forcings::Forcings, target::AbstractProcess, args...)
    return sum(map(f -> forcing(i, j, k, state, grid, f, target, args...), values(forcings.forcings)))
end

"""
    forcing(i, j, k, state, grid, forcing::OceananigansForcing, ::AbstractProcess, args...)

Return the value computed by the given `Oceananigans` forcing function, which should be an instance of
either `DiscreteForcing` or `ContinuousForcing`.
"""
@inline function forcing(i, j, k, state, grid, forcing::OceananigansForcing, ::AbstractProcess, args...)
    return forcing(i, j, k, get_field_grid(grid), state.clock, state)
end
