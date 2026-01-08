"""
    $TYPEDEF

Static vegetation root distribution implementation in PALADYN (Willeit 2016)
based on the scheme proposed by Zeng (2001). The root distribution is modeled as:

```math
f()
```

Properties:
$TYPEDFIELDS
"""
@kwdef struct StaticExponentialRootDistribution{NF} <: AbstractRootDistribution
    d1::NF = 7.0 # TODO PFT-specific parameter (here needleleaf tree)
    d2::NF = 2.0 # TODO PFT-specific parameter (here needleleaf tree)
end

StaticExponentialRootDistribution(::Type{NF}; kwargs...) where {NF} = StaticExponentialRootDistribution{NF}(; kwargs...)

variables(roots::StaticExponentialRootDistribution) = (
    auxiliary(:root_fraction, XY(), root_fraction, roots), # Static root fraction defined as function
)

"""
Returns a `FunctionField` that lazily computes the static root distribution on a 1D column grid.
"""
function root_fraction(roots::StaticExponentialRootDistribution{NF}, grid::AbstractColumnGrid, clock, fields) where {NF}
    fgrid = get_field_grid(grid)
    Δz = zspacings(fgrid, Center(), Center(), Center())
    ∂r∂z = FunctionField(fgrid, parameters=roots) do x, z, p
        0.5 * (p.d1 * exp(p.d1 * z) + p.d2 * exp(p.d2 * z))
    end
    r = ∂r∂z * Δz
    r_norm = r / sum(r, dims=3)
    return r_norm
end
