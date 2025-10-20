"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
@kwdef struct PrescribedAlbedo{NF} <: AbstractAlbedo
    "Surface emissivity, i.e. fraction of thermal radiation emitted from the surface [-]"
    emissivity::NF = 0.97
end

variables(::PrescribedAlbedo) = (
    input(:albedo, XY(), desc="Surface albedo, i.e. ratio of outgoing to incoming shortwave radiation [-]"),
)

@inline albedo(idx, state, ::PrescribedAlbedo) = state.albedo[idx...]

@inline emissivity(idx, state, albedo::PrescribedAlbedo) = albedo.emissivity

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

@inline albedo(idx, state, ::ConstantAlbedo) = albedo.albedo

@inline emissivity(idx, state, albedo::ConstantAlbedo) = albedo.emissivity
