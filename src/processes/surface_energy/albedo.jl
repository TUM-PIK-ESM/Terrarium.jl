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

@inline albedo(idx, state, ::PrescribedAlbedo) = state.albedo[idx...]

@inline emissivity(idx, state, ::PrescribedAlbedo) = state.emissivity[idx...]

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

@inline albedo(idx, state, albedo::ConstantAlbedo) = albedo.albedo

@inline emissivity(idx, state, albedo::ConstantAlbedo) = albedo.emissivity
