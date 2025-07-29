num_digits(x::Real) = Int(abs(ceil(log10(x))))

abstract type AbstractVerticalSpacing end

"""
    $SIGNATURES

Return the number of vertical layers defined by this discretization.
"""
get_npoints(spacing::AbstractVerticalSpacing) = spacing.N

"""
    $SIGNATURES

Return a `Vector` of vertical layer thicknesses according to the given discretization.
"""
get_spacing(spacing::AbstractVerticalSpacing) = spacing.(1:get_npoints(spacing))

Base.@kwdef struct UniformSpacing{NF} <: AbstractVerticalSpacing
    Δz::NF = 0.5
    N::Int = 100
end

(spacing::UniformSpacing)(i::Int) = spacing.Δz

Base.@kwdef struct ExponentialSpacing{NF} <: AbstractVerticalSpacing
    Δz_min::NF = 0.1
    Δz_max::NF = 500.0
    sig::Int = max(num_digits(Δz_min), num_digits(Δz_max))
    N::Int = 100
end

function (spacing::ExponentialSpacing)(i::Int)
    @assert 0 < i <= get_npoints(spacing) "index $i out of range"
    Δz₀ = log2(spacing.Δz_min)
    Δzₙ = log2(spacing.Δz_max)
    Δzᵢ = Δz₀ + (Δzₙ - Δz₀) / get_npoints(spacing)*i
    return round(exp2(Δzᵢ); sigdigits=spacing.sig)
end

struct ManualSpacing{NF} <: AbstractVerticalSpacing
    Δz::Vector{NF}
end

(spacing::ManualSpacing)(i::Int) = spacing.Δz[i]

get_npoints(spacing::ManualSpacing) = length(spacing.Δz)

