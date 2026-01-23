"""
    $TYPEDEF

Static vegetation root distribution implementation in PALADYN (Willeit 2016)
based on the scheme proposed by Zeng (2001). The PDF of the root distribution
is modeled as

```math
\\frac{\\partial R}{\\partial z} = \\frac{1}{2} \\left( a \\exp(a z) + b \\exp(b z) \\right)
```
which is then integrated over the soil column and normalized to sum to unity. Note that
this is effectively the average of two exponential distributions with rates `a` and `b`, both
with units m⁻¹. The resulting CDF of this distribution determines the root distribution.

Properties:
$FIELDS
"""
@kwdef struct StaticExponentialRootDistribution{NF} <: AbstractRootDistribution
    "First empirical rate parameter for root distribution"
    a::NF = 7.0 # TODO PFT-specific parameter (here needleleaf tree)
    
    "Second empirical rate parameter for root distribution"
    b::NF = 2.0 # TODO PFT-specific parameter (here needleleaf tree)
end

StaticExponentialRootDistribution(::Type{NF}; kwargs...) where {NF} = StaticExponentialRootDistribution{NF}(; kwargs...)

variables(roots::StaticExponentialRootDistribution) = (
    auxiliary(:root_fraction, XYZ(), root_fraction, roots), # Static root fraction defined as function
)

"""
Returns a `FunctionField` that lazily computes the static root distribution on a 1D column grid.
"""
function root_fraction(roots::StaticExponentialRootDistribution{NF}, grid::AbstractColumnGrid, clock, fields) where {NF}
    fgrid = get_field_grid(grid)
    # define pdf of root distribution as a continuous function of depth
    ∂R∂z = FunctionField{Center, Center, Center}(fgrid, parameters=roots) do x, z, p
        0.5 * (p.a * exp(p.a * z) + p.b * exp(p.b * z))
    end
    # scale by the layer thicknesses
    Δz = zspacings(fgrid, Center(), Center(), Center())
    R = ∂R∂z * Δz
    # and normalize
    R_norm = R / sum(R, dims=3)
    return R_norm
end
