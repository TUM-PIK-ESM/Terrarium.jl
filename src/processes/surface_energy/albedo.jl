"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
@kwdef struct PrescribedAlbedo{NF} <: AbstractAlbedo{NF} end

PrescribedAlbedo(::Type{NF}) where {NF} = PrescribedAlbedo{NF}()

variables(::PrescribedAlbedo) = (
    input(:albedo, XY(), domain = UnitInterval(), desc = "Surface albedo, i.e. ratio of outgoing to incoming shortwave radiation [-]"),
    input(:emissivity, XY(), domain = UnitInterval(), desc = "Surface emissivity, i.e. efficiency of longwave emission [-]"),
)

"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
@kwdef struct ConstantAlbedo{NF} <: AbstractAlbedo{NF}
    "Surface albedo, i.e. ratio of outgoing to incoming shortwave radiation [-]"
    albedo::NF = 0.3

    "Surface emissivity, i.e. fraction of thermal radiation emitted from the surface [-]"
    emissivity::NF = 0.97
end

ConstantAlbedo(::Type{NF}; kwargs...) where {NF} = ConstantAlbedo{NF}(; kwargs...)

@inline albedo(i, j, grid, fields, albedo::ConstantAlbedo) = albedo.albedo

@inline emissivity(i, j, grid, fields, albedo::ConstantAlbedo) = albedo.emissivity
