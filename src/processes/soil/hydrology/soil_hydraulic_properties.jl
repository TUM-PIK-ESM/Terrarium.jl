# Interface

"""
Base type for unsaturated hydraulic conductivity parameterizations.
"""
abstract type AbstractUnsatK end

"""
    get_swrc(::AbstractUnsatK)

Return the soil water retention curve associated with the given unsaturated hydraulic conductivity scheme.
"""
function get_swrc end

"""
Base type for soil hydraulic properties and parameterization schemes.
"""
abstract type AbstractSoilHydraulics{NF, UnsatK} end

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

# Hydraulics parameterizations

"""
    ConstantHydraulics{NF} <: AbstractSoilHydraulics

Represents a simple case where soil hydraulic properties are given as constant values.
This is mostly provided just for testing, although it may be useful in certain cases where direct
measurements of hydraulic properites are available.

Properties:
$TYPEDFIELDS
"""
@kwdef struct ConstantHydraulics{NF, UnsatK} <: AbstractSoilHydraulics{NF, UnsatK}
    "Unsaturated hydraulic conductivity formulation; defaults to `cond_sat`"
    cond_unsat::UnsatK

    "Hydraulic conductivity at saturation [m/s]"
    cond_sat::NF = 1e-5

    "Prescribed soil porosity [-]"
    porosity::NF = 0.49

    "Prescribed field capacity [-]"
    field_capacity::NF = 0.25

    "Prescribed wilting point [-]"
    wilting_point::NF = 0.05
end

ConstantHydraulics(::Type{NF}; cond_unsat=UnsatKLinear(NF), kwargs...) where {NF} = ConstantHydraulics{NF, typeof(cond_unsat)}(; cond_unsat, kwargs...)

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
@kwdef struct SoilHydraulicsSURFEX{NF, UnsatK} <: AbstractSoilHydraulics{NF, UnsatK}
    "Unsaturated hydraulic conductivity formulation; defaults to `cond_sat`"
    cond_unsat::UnsatK

    "Hydraulic conductivity at saturation [m/s]"
    cond_sat::NF = 1e-5

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

SoilHydraulicsSURFEX(::Type{NF}; cond_unsat=UnsatKLinear(NF), kwargs...) where {NF} = SoilHydraulicsSURFEX{NF, typeof(cond_unsat)}(; cond_unsat, kwargs...)

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
    swrc::RetentionCurve
end

UnsatKLinear(::Type{NF}; swrc=FreezeCurves.BrooksCorey()) where {NF} = UnsatKLinear(adapt(NF, ustrip(swrc)))

function hydraulic_conductivity(
    hydraulics::AbstractSoilHydraulics{NF, <:UnsatKLinear},
    soil::SoilVolume
) where {NF}
    let fracs = volumetric_fractions(soil),
        θw = fracs.water, # unfrozen water content
        θsat = fracs.water + fracs.ice + fracs.air, # water + ice content at saturation (porosity)
        K_sat = saturated_hydraulic_conductivity(hydraulics, soil.solid.texture);
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
struct UnsatKVanGenuchten{NF, RetentionCurve<:FreezeCurves.VanGenuchten} <: AbstractUnsatK
    "Exponential scaling factor for ice impedance"
    impedance::NF

    "Van Genuchten soil water retention curve"
    swrc::RetentionCurve
end

UnsatKVanGenuchten(::Type{NF}; impedance::NF = NF(7), swrc=FreezeCurves.VanGenuchten()) where {NF} = UnsatKVanGenuchten(NF(impedance), adapt(NF, ustrip(swrc)))

function hydraulic_conductivity(
    hydraulics::AbstractSoilHydraulics{NF, <:UnsatKVanGenuchten},
    soil::SoilVolume,
) where {NF}
    let n = hydraulics.cond_unsat.swrc.n, # van Genuchten parameter `n`
        fracs = volumetric_fractions(soil),
        θw = fracs.water, # unfrozen water content
        θwi = fracs.water + fracs.ice, # total water + ice content
        θsat = porosity(soil), # porosity
        f = liquid_fraction(soil),
        Ω = hydraulics.cond_unsat.impedance, # scaling parameter for ice impedance
        I_ice = 10^(-Ω*(1 - f)), # ice impedance factor
        K_sat = saturated_hydraulic_conductivity(hydraulics, soil.solid.texture);
        # We use `complex` types here to permit illegal state values which may occur when using adaptive time steppers.
        K = abs(K_sat*I_ice*sqrt(complex(θw/θsat))*(1 - complex(1 - complex(θw/θsat)^(n/(n+1)))^((n-1)/n))^2)
        return K
    end
end
