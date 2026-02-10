# Interface

"""
    $TYPEDEF

Base type for unsaturated hydraulic conductivity parameterizations.
"""
abstract type AbstractUnsatK{NF} end

"""
    get_swrc(::AbstractUnsatK)

Return the soil water retention curve associated with the given unsaturated hydraulic conductivity scheme.
"""
function get_swrc end

"""
    $TYPEDEF

Base type for soil hydraulic properties and parameterization schemes.
"""
abstract type AbstractSoilHydraulics{NF, RC <: SWRC, UnsatK <: AbstractUnsatK} end

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
    $TYPEDEF

Represents a simple case where soil hydraulic properties are given as constant values.
This is mostly provided just for testing, although it may be useful in certain cases where direct
measurements of hydraulic properites are available.

Properties:
$TYPEDFIELDS
"""
@kwdef struct ConstantSoilHydraulics{NF, RC, UnsatK <: AbstractUnsatK{NF}} <: AbstractSoilHydraulics{NF, RC, UnsatK}
    "Soil water retention curve"
    swrc::RC

    "Unsaturated hydraulic conductivity formulation; defaults to `sat_hydraulic_cond`"
    unsat_hydraulic_cond::UnsatK

    "Hydraulic conductivity at saturation [m/s]"
    sat_hydraulic_cond::NF = 1.0e-5

    "Prescribed field capacity [-]"
    field_capacity::NF = 0.25

    "Prescribed wilting point [-]"
    wilting_point::NF = 0.05
end

function ConstantSoilHydraulics(
        ::Type{NF};
        swrc = BrooksCorey(),
        unsat_hydraulic_cond = UnsatKLinear(NF),
        kwargs...
    ) where {NF}
    swrc = adapt(NumberFormatAdaptor{NF}(), ustrip(swrc))
    return ConstantSoilHydraulics{NF, typeof(swrc), typeof(unsat_hydraulic_cond)}(; swrc, unsat_hydraulic_cond, kwargs...)
end

@inline saturated_hydraulic_conductivity(hydraulics::ConstantSoilHydraulics, args...) = hydraulics.sat_hydraulic_cond

@inline wilting_point(hydraulics::ConstantSoilHydraulics, args...) = hydraulics.wilting_point

@inline field_capacity(hydraulics::ConstantSoilHydraulics, args...) = hydraulics.field_capacity

"""
    $TYPEDEF

Soil hydraulics parameterization that includes the SURFEX (Masson et al. 2013) formulation of field capacity
and wilting point as a function of soil texture.

Properties:
$TYPEDFIELDS
"""
@kwdef struct SoilHydraulicsSURFEX{NF, RC, UnsatK <: AbstractUnsatK{NF}} <: AbstractSoilHydraulics{NF, RC, UnsatK}
    "Soil water retention curve"
    swrc::RC

    "Unsaturated hydraulic conductivity formulation; defaults to `sat_hydraulic_cond`"
    unsat_hydraulic_cond::UnsatK

    "Hydraulic conductivity at saturation [m/s]"
    sat_hydraulic_cond::NF = 1.0e-5

    "Linear coeficient of wilting point adjustment due to clay content [-]"
    wilting_point_coef::NF = 37.13e-3

    "Linear coeficient of field capacity adjustment due to clay content [-]"
    field_capacity_coef::NF = 89.0e-3

    "Exponent of field capacity adjustment due to clay content [-]"
    field_capacity_exp::NF = 0.35
end

function SoilHydraulicsSURFEX(
        ::Type{NF};
        swrc = BrooksCorey(),
        unsat_hydraulic_cond = UnsatKLinear(NF),
        kwargs...
    ) where {NF}
    swrc = adapt(NumberFormatAdaptor{NF}(), ustrip(swrc))
    return SoilHydraulicsSURFEX{NF, typeof(swrc), typeof(unsat_hydraulic_cond)}(; swrc, unsat_hydraulic_cond, kwargs...)
end

# TODO: this is not quite correct, SURFEX uses a hydraulic conductivity function that decreases exponentially with depth
@inline saturated_hydraulic_conductivity(hydraulics::SoilHydraulicsSURFEX, args...) = hydraulics.sat_hydraulic_cond

@inline function wilting_point(hydraulics::SoilHydraulicsSURFEX, texture::SoilTexture)
    β_w = hydraulics.wilting_point_coef
    wp = β_w * sqrt(texture.clay * 100)
    return wp
end

@inline function field_capacity(hydraulics::SoilHydraulicsSURFEX, texture::SoilTexture)
    η = hydraulics.field_capacity_exp
    β_c = hydraulics.field_capacity_coef
    fc = β_c * (texture.clay * 100)^η
    return fc
end

# Unsaturated hydraulic conductivity schemes

"""
    $TYPEDEF

Simple formulation of hydraulic conductivity as a linear function of the liquid water saturated fraction,
i.e. `soil.water / (soil.water + soil.ice + soil.air)`.
"""
struct UnsatKLinear{NF} <: AbstractUnsatK{NF} end

UnsatKLinear(::Type{NF}) where {NF} = UnsatKLinear{NF}()

function hydraulic_conductivity(
        hydraulics::AbstractSoilHydraulics{NF, RC, UnsatKLinear{NF}},
        soil::SoilVolume
    ) where {NF, RC}
    let fracs = volumetric_fractions(soil),
            θw = fracs.water, # unfrozen water content
            θsat = fracs.water + fracs.ice + fracs.air, # water + ice content at saturation (porosity)
            K_sat = saturated_hydraulic_conductivity(hydraulics, soil.solid.texture)
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
struct UnsatKVanGenuchten{NF} <: AbstractUnsatK{NF}
    "Exponential scaling factor for ice impedance"
    impedance::NF
end

UnsatKVanGenuchten(::Type{NF}; impedance::NF = NF(7)) where {NF} = UnsatKVanGenuchten(NF(impedance))

function hydraulic_conductivity(
        hydraulics::AbstractSoilHydraulics{NF, <:VanGenuchten, UnsatKVanGenuchten{NF}},
        soil::SoilVolume
    ) where {NF}
    # TODO: The SWRC parameters will need to also be spatially varying at some point
    let n = hydraulics.swrc.n, # van Genuchten parameter `n`
            fracs = volumetric_fractions(soil),
            θw = fracs.water, # unfrozen water content
            θwi = fracs.water + fracs.ice, # total water + ice content
            θsat = porosity(soil), # porosity
            f = liquid_fraction(soil),
            Ω = hydraulics.unsat_hydraulic_cond.impedance, # scaling parameter for ice impedance
            I_ice = 10^(-Ω * (1 - f)), # ice impedance factor
            K_sat = saturated_hydraulic_conductivity(hydraulics, soil.solid.texture)
        # We use `complex` types here to permit illegal state values which may occur when using adaptive time steppers.
        K = abs(K_sat * I_ice * sqrt(complex(θw / θsat)) * (1 - complex(1 - complex(θw / θsat)^(n / (n + 1)))^((n - 1) / n))^2)
        return K
    end
end
