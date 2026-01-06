"""
    $TYPEDEF

Vegetation root distributions implementation in PALADYN (Willeit 2016)
based on the scheme proposed by Zeng (2001).

Properties:
$TYPEDFIELDS
"""
@kwdef struct PALDAYNRootDistribution{NF} <: AbstractRootDistribution
    d1::NF = 7.0 # TODO PFT-specific parameter (here needleleaf tree)
    d2::NF = 2.0 # TODO PFT-specific parameter (here needleleaf tree)
end

PALDAYNRootDistribution(::Type{NF}; kwargs...) where {NF} = PALDAYNRootDistribution{NF}(; kwargs...)

variables(rd::PALDAYNRootDistribution) = (
    auxiliary(:root_fraction, XY(), root_fraction, rd), # Static root fraction defined as function
)

function root_fraction(rd::PALDAYNRootDistribution{NF}, grid, clock, fields) where {NF}
    return FunctionField(grid, parameters=rd) do x, y, z, p
        NF(1) - NF(0.5) * (exp(-p.d1 * z) + exp(-p.d2 * z))
    end
end
