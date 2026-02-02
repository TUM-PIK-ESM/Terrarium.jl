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

@inline @propagate_inbounds stomatal_conductance(i, j, grid, state, ::MedlynStomatalConductance) = state.gw_can[i, j]

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

function compute_auxiliary!(state, model, stomcond::MedlynStomatalConductance)
    grid = get_grid(model)
    atmos = get_atmosphere(model)
    constants = get_constants(model)
    photo = get_photosynthesis(model)
    launch!(grid, XY, compute_auxiliary_kernel!, state, stomcond, photo, atmos, constants)
end

@kernel function compute_auxiliary_kernel!(
    state, grid,
    stomcond::MedlynStomatalConductance{NF},
    photo::LUEPhotosynthesis{NF},
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
) where NF
    i, j = @index(Global, NTuple)

    @inbounds let An = state.An[i, j],
                  CO2 = state.CO2[i, j],
                  LAI = state.LAI[i, j],
                  β = state.SMLF[i, j];
                
        # Compute vpd [Pa]
        vpd = compute_vpd(i, j, grid, state, atmos, constants)

        # Compute conducatance
        gw_can = compute_gw_can(stomcond, photo, vpd, An, CO2, LAI, β)

        # Compute λc
        λc = compute_λc(stomcond, vpd)

        # Store result
        state.gw_can[i, j, 1] = gw_can
        state.λc[i, j, 1] = λc
    end
end
