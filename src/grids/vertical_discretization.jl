abstract type AbstractVerticalSpacing{NF} end

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

Base.@kwdef struct UniformSpacing{NF} <: AbstractVerticalSpacing{NF}
    Δz::NF = 0.5
    N::Int = 50
end

(spacing::UniformSpacing)(i::Int) = spacing.Δz

Base.@kwdef struct ExponentialSpacing{NF,ST<:Union{Nothing,Integer}} <: AbstractVerticalSpacing{NF}
    Δz_min::NF = 0.1
    Δz_max::NF = 100.0
    sig::ST = 3
    N::Int = 50

    function ExponentialSpacing(Δz_min::NF, Δz_max::NF, sig, N) where {NF}
        @assert N > 1 "number of grid points for exponential spacing must be > 1"
        return new{NF,typeof(sig)}(Δz_min, Δz_max, sig, N)
    end
end

function (spacing::ExponentialSpacing)(i::Int)
    @assert 0 < i <= get_npoints(spacing) "index $i out of range"
    logΔz₀ = log2(spacing.Δz_min)
    logΔzₙ = log2(spacing.Δz_max)
    logΔzᵢ = logΔz₀ + (i-1)*(logΔzₙ - logΔz₀) / (get_npoints(spacing)-1)
    if isnothing(spacing.sig)
        return exp2(logΔzᵢ)
    else
        return round(exp2(logΔzᵢ); sigdigits=spacing.sig)
    end
end

@kwdef struct PrescribedSpacing{NF} <: AbstractVerticalSpacing{NF}
    Δz::Vector{NF}
end

(spacing::PrescribedSpacing)(i::Int) = spacing.Δz[i]

get_npoints(spacing::PrescribedSpacing) = length(spacing.Δz)

