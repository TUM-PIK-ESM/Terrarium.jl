"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
@kwdef struct PrescribedAlbedo <: AbstractAlbedo end

variables(::PrescribedAlbedo) = (
    input(:albedo, XY(), domain=UnitInterval(), desc="Surface albedo, i.e. ratio of outgoing to incoming shortwave radiation [-]"),
    input(:emissivity, XY(), domain=UnitInterval(), desc="Surface emissivity, i.e. efficiency of longwave emission [-]"),
)

@inline albedo(i, j, state, ::PrescribedAlbedo) = @inbounds state.albedo[i, j]

@inline emissivity(i, j, state, ::PrescribedAlbedo) = @inbounds state.emissivity[i, j]

"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
@kwdef struct ConstantAlbedo{NF} <: AbstractAlbedo
    "Surface albedo, i.e. ratio of outgoing to incoming shortwave radiation [-]"
    albedo::NF = 0.3

    "Surface emissivity, i.e. fraction of thermal radiation emitted from the surface [-]"
    emissivity::NF = 0.97
end

ConstantAlbedo(::Type{NF}; kwargs...) where {NF} = ConstantAlbedo{NF}(; kwargs...)

@inline albedo(i, j, state, albedo::ConstantAlbedo) = albedo.albedo

@inline emissivity(i, j, state, albedo::ConstantAlbedo) = albedo.emissivity
