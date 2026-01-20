"""
    $TYPEDEF

Static vegetation root distribution implementation in PALADYN (Willeit 2016)
based on the scheme proposed by Zeng (2001). The PDF of the root distribution
is modeled as

```math
\\frac{\\partial R}{\\partial z} = \\frac{1}{2} \\left( d_1 \\exp(d_1 z) + d_2 \\exp(d_2 z) \\right)
```
which is then integrated over the soil column and normalized to sum to unity. Note that
this is effectively the average of two exponential distributions with rates `d1` and `d2`.

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
    # define pdf of root distribution as a continuous function of depth
    ∂R∂z = FunctionField(fgrid, parameters=roots) do x, z, p
        0.5 * (p.d1 * exp(p.d1 * z) + p.d2 * exp(p.d2 * z))
    end
    # scale by the layer thicknesses
    Δz = zspacings(fgrid, Center(), Center(), Center())
    R = ∂R∂z * Δz
    # and normalize
    R_norm = R / sum(r, dims=3)
    return R_norm
end
