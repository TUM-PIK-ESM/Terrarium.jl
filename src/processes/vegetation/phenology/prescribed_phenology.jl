"""
    $TYPEDEF

Prescribed vegetation phenology where `leaf_area_index` is treated a (possibly time-varying) input variable.

Properties:
$TYPEDFIELDS
"""
@parameterized @kwdef struct PrescribedPhenology{NF} <: AbstractPhenology{NF}
    "Seasonal maximum LAI used to determine the phenology factor"
    @param max_leaf_area_index::NF = 5.0 (bounds = Positive,)
end

PrescribedPhenology(::Type{NF}; kwargs...) where {NF} = PrescribedPhenology{NF}(; kwargs...)

variables(phenol::PrescribedPhenology) = (
    auxiliary(:phenology_factor, XY(), kernel(phenology_factor, phenol)),
    input(:leaf_area_index, XY()), # Leaf Area Index [m²/m²]
)

function phenology_factor(i, j, grid, fields, phenology::PrescribedPhenology{NF}) where {NF}
    LAI = fields.leaf_area_index[i, j]
    LAI_max = phenology.max_leaf_area_index
    ϕ = clamp(LAI / LAI_max, NF(0), NF(1))
    return ϕ
end

@inline compute_auxiliary!(state, grid, phenology::PrescribedPhenology, args...) = nothing

@inline compute_tendencies!(state, grid, phenology::PrescribedPhenology, args...) = nothing
