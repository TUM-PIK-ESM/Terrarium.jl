# Interface

"""
Base type for unsaturated hydraulic conductivity parameterizations.
"""
abstract type AbstractUnsatK end

"""
Base type for soil hydraulic properties and parameterization schemes.
"""
abstract type AbstractSoilHydraulics{UK<:AbstractUnsatK} end

"""
    $SIGNATURES

Compute hydraulic conductivity at saturation.
"""
function saturated_hydraulic_conductivity end

"""
    $SIGNATURES

Compute (variably saturated) hydraulic conductivity based on the given hydraulic properties,
soil water retention curve (SWRC), and volumetric fractions.
"""
function hydraulic_conductivity end

"""
    $SIGNATURES

Compute the natural porosity of the mineral soil constitutents, i.e. excluding organic material.
"""
function mineral_porosity end

"""
    $SIGNATURES

Compute the empirical wilting point of the soil.
"""
function wilting_point end

"""
    $SIGNATURES

Compute the empirical field capacity of the soil.
"""
function field_capacity end

# Hydraulic parameterizations

"""
    ConstantHydraulics{NF} <: AbstractSoilHydraulics

Represents a simple case where soil hydraulic properties are given as constant values.
This is mostly provided just for testing, although it may be useful in certain cases where direct
measurements of hydraulic properites are available.

Properties:
$TYPEDFIELDS
"""
@kwdef struct ConstantHydraulics{NF, UK} <: AbstractSoilHydraulics{UK}
    "Hydraulic conductivity at saturation [m/s]"
    cond_sat::NF = 1e-5

    "Unsaturated hydraulic conductivity formulation; defaults to `cond_sat`"
    cond_unsat::UK = cond_sat

    "Prescribed soil porosity [-]"
    porosity::NF = 0.49

    "Prescribed field capacity [-]"
    field_capacity::NF = 0.1

    "Prescribed wilting point [-]"
    wilting_point::NF = 0.05
end

ConstantHydraulics(::Type{NF}; kwargs...) where {NF} = ConstantHydraulics{NF}(; kwargs...)

@inline saturated_hydraulic_conductivity(hydraulics::ConstantHydraulics, args...) = hydraulics.cond_sat

@inline mineral_porosity(hydraulics::ConstantHydraulics, args...) = hydraulics.porosity

@inline wilting_point(hydraulics::ConstantHydraulics, args...) = hydraulics.wilting_point

@inline field_capacity(hydraulics::ConstantHydraulics, args...) = hydraulics.field_capacity

"""
    $TYPEDEF

SURFEX parameterization of mineral soil porosity (Masson et al. 2013).

Properties:
$TYPEDFIELDS
"""
@kwdef struct SoilHydraulicsSURFEX{NF, UK} <: AbstractSoilHydraulics{UK}
    "Hydraulic conductivity at saturation [m/s]"
    cond_sat::NF = 1e-5

    "Unsaturated hydraulic conductivity formulation; defaults to `cond_sat`"
    cond_unsat::UK = cond_sat

    "Base porosity of soil without any sand [-]"
    porosity::NF = 0.49

    "Linear coeficient of porosity adjustment due to sand content [-]"
    porosity_sand_coef::NF = -1.1e-3
    
    "Linear coeficient of wilting point adjustment due to clay content [-]"
    wilting_point_coef::NF = 37.13e-3
    
    "Linear coeficient of field capacity adjustment due to clay content [-]"
    field_capacity_coef::NF = 89.0e-3
    
    "Exponent of field capacity adjustment due to clay content [-]"
    field_capacity_exp::NF = 0.35
end

SoilHydraulicsSURFEX(::Type{NF}; kwargs...) where {NF} = SoilHydraulicsSURFEX{NF}(; kwargs...)

# TODO: this is not quite correct, SURFEX uses a hydraulic conductivity function that decreases exponentially with depth
@inline saturated_hydraulic_conductivity(hydraulics::SoilHydraulicsSURFEX, args...) = hydraulics.cond_sat

@inline function mineral_porosity(hydraulics::SoilHydraulicsSURFEX, texture::SoilTexture)
    p₀ = hydraulics.porosity
    β_s = hydraulics.porosity_sand_coef
    por = p₀ + β_s*texture.sand*100
    return por
end

@inline function wilting_point(hydraulics::SoilHydraulicsSURFEX, texture::SoilTexture)
    β_w = hydraulics.wilting_point_coef
    wp = β_w*sqrt(texture.clay*100)
    return wp
end

@inline function field_capacity(hydraulics::SoilHydraulicsSURFEX, texture::SoilTexture)
    η = hydraulics.field_capacity_exp
    β_c = hydraulics.field_capacity_coef
    fc = β_c*(texture.clay*100)^η
    return fc
end

# Unsaturated hydraulic conductivity schemes

"""
    $TYPEDEF

Simple formulation of hydraulic conductivity as a linear function of the liquid water saturated fraction,
i.e. `soil.water / (soil.water + soil.ice + soil.air)`.
"""
struct UnsatKLinear{RetentionCurve<:SWRC} <: AbstractUnsatK
    "Soil water retention curve"
    swrc::RetentionCurve = ustrip(FreezeCurves.BrooksCorey())
end

function hydraulic_conductivity(hydraulics::AbstractSoilHydraulics{<:UnsatKLinear}, soil::SoilComposition)
    let θw = soil.water, # unfrozen water content
        θsat = soil.water + soil.ice + soil.air, # water + ice content at saturation (porosity)
        K_sat = saturated_hydraulic_conductivity(hydraulics, soil.texture);
        K = K_sat * θw / θsat
        return K
    end
end

"""
    $TYPEDEF

Formulation of hydraulic conductivity as a function of saturated hydraulic conductivity `K_sat` and
volumetric fractions, assumed to include those of water, ice, and air.

See van Genuchten (1980) and Westermann et al. (2023).
"""
@kwdef struct UnsatKVanGenuchten{NF, RetentionCurve<:FreezeCurves.VanGenuchten} <: AbstractUnsatK
    "Exponential scaling factor for ice impedance"
    Ω::NF = 7

    "Van Genuchten soil water retention curve"
    swrc::RetentionCurve = ustrip(FreezeCurves.VanGenuchten())
end

function hydraulic_conductivity(
    hydraulics::AbstractSoilHydraulics{<:UnsatKVanGenuchten},
    soil::SoilComposition,
)
    let n = swrc.n, # van Genuchten parameter `n`
        θw = soil.water, # unfrozen water content
        θwi = soil.water + soil.ice, # total water + ice content
        θsat = θwi + soil.air, # water + ice content at saturation (porosity)
        Ω = hydraulics.cond_unsat.Ω, # scaling parameter for ice impedance
        I_ice = 10^(-Ω*(1 - θw/θtot)), # ice impedance factor
        K_sat = saturated_hydraulic_conductivity(hydraulics, soil.texture);
        # We use `complex` types here to permit illegal state values which may occur when using adaptive time steppers.
        K = abs(K_sat*I_ice*sqrt(complex(θw/θsat))*(1 - complex(1 - complex(θw/θsat)^(n/(n+1)))^((n-1)/n))^2)
        return K
    end
end
