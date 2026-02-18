"""
    $TYPEDEF

Base type for vertical discretizations.
"""
abstract type AbstractVerticalSpacing{NF} end

"""
    $SIGNATURES

Return the number of vertical layers defined by this discretization.
"""
num_layers(spacing::AbstractVerticalSpacing) = spacing.N

"""
    $SIGNATURES

Return a `Vector` of vertical layer thicknesses according to the given discretization.
"""
get_spacing(spacing::AbstractVerticalSpacing) = spacing.(1:num_layers(spacing))

"""
    $TYPEDEF

Uniform vertical discretization with `N` layers of size `Δz`.

Properties:
$TYPEDFIELDS
"""
Base.@kwdef struct UniformSpacing{NF} <: AbstractVerticalSpacing{NF}
    Δz::NF = 0.1
    N::Int = 100
end

(spacing::UniformSpacing)(i::Int) = spacing.Δz

"""
    $TYPEDEF

Variably-spaced vertical discretization with `N` layers increasing quasi-exponentially in thickness from
`Δz_min` at the top (surface) to `Δz_max` at the bottom. The integer property `sig` determines to what
significant digit each layer thickness should be rounded.

Properties:
$TYPEDFIELDS
"""
Base.@kwdef struct ExponentialSpacing{NF, ST <: Union{Nothing, Integer}} <: AbstractVerticalSpacing{NF}
    "Minimum layer thickness at the surface"
    Δz_min::NF = 0.05

    "Maximum layer thickness at the bottom"
    Δz_max::NF = 100.0

    "Number of layers"
    N::Int = 50

    "Number of significant digits for rounding or `nothing`"
    sig::ST = 3

    function ExponentialSpacing(Δz_min::NF, Δz_max::NF, N, sig) where {NF}
        @assert N > 1 "number of grid points for exponential spacing must be > 1"
        return new{NF, typeof(sig)}(Δz_min, Δz_max, N, sig)
    end
end

function (spacing::ExponentialSpacing)(i::Int)
    @assert 0 < i <= num_layers(spacing) "index $i out of range"
    logΔz₀ = log2(spacing.Δz_min)
    logΔzₙ = log2(spacing.Δz_max)
    logΔzᵢ = logΔz₀ + (i - 1) * (logΔzₙ - logΔz₀) / (num_layers(spacing) - 1)
    if isnothing(spacing.sig)
        return exp2(logΔzᵢ)
    else
        return round(exp2(logΔzᵢ); sigdigits = spacing.sig)
    end
end

"""
    $TYPEDEF

Vertical discretization with prescribed thicknesses for each layer. The number of layers
is equal to the length of the given vector.

Properties:
$TYPEDFIELDS
"""
@kwdef struct PrescribedSpacing{NF} <: AbstractVerticalSpacing{NF}
    Δz::Vector{NF}
end

(spacing::PrescribedSpacing)(i::Int) = spacing.Δz[i]

num_layers(spacing::PrescribedSpacing) = length(spacing.Δz)
