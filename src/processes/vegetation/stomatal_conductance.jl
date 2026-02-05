"""
    $TYPEDEF
Stomatal conductance implementation from PALADYN (Willeit 2016) following the optimal stomatal conductance model
(Medlyn et al. 2011).

Authors: Maha Badri and Matteo Willeit

Properties:
$TYPEDFIELDS
"""
@kwdef struct MedlynStomatalConductance{NF} <: AbstractStomatalConductance{NF}
    "Parameter in optimal stomatal conductance formulation representing the quasi-linear
    relationship between conductance and net assimilation, Lin et al. 2015 [-], PFT specific"
    g₁::NF = 2.3 # TODO: value for Needleleaf tree PFT

    "Minimum stomatal condutance parameter [mm s⁻¹]"
    g_min::NF = 0.5
end

MedlynStomatalConductance(::Type{NF}; kwargs...) where {NF} = MedlynStomatalConductance(; kwargs...)

variables(::MedlynStomatalConductance) = (
    auxiliary(:gw_can, XY(), units=u"m/s"), # Canopy conducatance for water vapor
    auxiliary(:λc, XY()), # Ratio of leaf-internal and air CO2 concentration [-]
)

@inline @propagate_inbounds stomatal_conductance(i, j, grid, fields, ::MedlynStomatalConductance) = fields.gw_can[i, j]

@inline function compute_gw_can(
    stomcond::MedlynStomatalConductance{NF},
    photo::LUEPhotosynthesis{NF},
    vpd, An, co2, LAI, β
) where NF
    let g_min = stomcond.g_min / 1000, # convert mm/s to m/s
        g₁ = stomcond.g₁,
        k_ext = photo.k_ext;
        g₀ = g_min * (1 - exp(-k_ext * LAI)) * β
        g_can = g₀ + (1 + g₁ / sqrt(vpd)) * An / co2 * NF(1e6)
        return g_can
    end
end

"""
    $SIGNATURES

Computes the ratio of leaf-internal and air CO2 concentration `λc`, 
derived from the optimal stomatal conductance model (Medlyn et al. 2011),
Eq. 71, PALADYN (Willeit 2016).
"""
@inline function compute_λc(stomcond::MedlynStomatalConductance{NF}, vpd) where NF
    λc = NF(1.0) - NF(1.6) / (NF(1.0) + stomcond.g₁ / sqrt(vpd * NF(1.0e-3)))
    return λc
end

# Process methods

function compute_auxiliary!(
    state, grid,
    stomcond::MedlynStomatalConductance,
    photo::LUEPhotosynthesis,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants,
    args...
)
    out = auxiliary_fields(state, stomcond)
    fields = get_fields(state, stomcond, photo, atmos, constants; except = out)
    launch!(grid, XY, compute_stomatal_conductance_kernel!, out, fields, stomcond, photo, atmos, constants)
end

# Kernel functions

@propagate_inbounds function compute_stomatal_conductance(i, j, grid, fields, stomcond::MedlynStomatalConductance{NF}, photo::LUEPhotosynthesis{NF}, atmos::AbstractAtmosphere, constants::PhysicalConstants, args...) where NF
    # Get inputs
    An = fields.An[i, j]
    CO2 = fields.CO2[i, j]
    LAI = fields.LAI[i, j]
    β = fields.SMLF[i, j]

    # Compute vpd [Pa]
    vpd = compute_vpd(i, j, grid, fields, atmos, constants)

    # Compute conductance gw_can and internal CO2 ratio λc
    gw_can = compute_gw_can(stomcond, photo, vpd, An, CO2, LAI, β)
    λc = compute_λc(stomcond, vpd)

    return gw_can, λc
end

@propagate_inbounds function compute_stomatal_conductance!(out, i, j, grid, fields, stomcond::MedlynStomatalConductance{NF}, photo::LUEPhotosynthesis{NF}, atmos::AbstractAtmosphere, constants::PhysicalConstants, args...) where NF
    gw_can, λc = compute_stomatal_conductance(i, j, grid, fields, stomcond, photo, atmos, constants, args...)
    out.gw_can[i, j, 1] = gw_can
    out.λc[i, j, 1] = λc
    return out
end

# Kernels

@kernel function compute_stomatal_conductance_kernel!(out, grid, fields, stomcond::AbstractStomatalConductance, args...)
    i, j = @index(Global, NTuple)
    compute_stomatal_conductance!(out, i, j, grid, fields, stomcond, args...)
end
