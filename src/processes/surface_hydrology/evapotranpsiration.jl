"""
    $TYPEDEF

Evapotranspiration scheme from PALADYN (Willeit 2016) that includes a canopy
evaporation term based on the saturation fraction of canopy water defined by the
canopy hydrology scheme.

```math
E_g &= \\beta \\frac{\\Delta q}{r_a} \\
E_c &= f_{\\text{can}} \\frac{\\Delta q}{r_a} \\
T_c &= \\frac{\\Delta q}{r_a + r_s} \\
```

Properties:
$FIELDS
"""
struct PALADYNCanopyEvapotranspiration{NF, GR <: AbstractGroundEvaporationResistanceFactor} <: AbstractEvapotranspiration{NF}
    "Drag coefficient for the traansfer of heat and water between the ground and canopy"
    C_can::NF

    "Parameterization for ground resistance to evaporation/sublimation"
    ground_resistance::GR
end

function PALADYNCanopyEvapotranspiration(
        ::Type{NF};
        C_can::NF = 0.006,
        ground_resistance = ConstantEvaporationResistanceFactor(typeof(C_can))
    ) where {NF}
    return PALADYNCanopyEvapotranspiration{NF, typeof(ground_resistance)}(C_can, ground_resistance)
end

# TODO: The following ET functions all have the same basic functional form and can be generalized to a function
# that takes the humidity gradient/difference and an arbitrary number of resistance terms. We could consider
# making such an abstraction, but we should first consider what exactly the benefits would be since it also
# obfuscates that actual calculations, and the equations are quite simple. Perhaps one benefit would be reduced
# unit testing overhead?

"""
    $TYPEDSIGNATURES

Compute potential transpiration from the given humidity gradient, aerodynamic resistance `rₐ` and stomatal conductance `gw_can`.
"""
@inline function compute_transpiration(::PALADYNCanopyEvapotranspiration{NF}, Δq, rₐ, gw_can) where {NF}
    let rₛ = 1 / max(gw_can, sqrt(eps(NF)))  # stomatal resistance as reciprocal of conductance
        # Calculate transpriation flux in m/s (positive upwards)
        E_trp = Δq / (rₐ + rₛ)
        return E_trp
    end
end

"""
    $TYPEDSIGNATURES

Compute potential evaporation from the ground below the canopy, following Eq. 5, PALADYN (Willeit 2016);
`Δq` is the humidity gradient, `β` is the ground evaporation resistance factor, `rₐ` is aerodynamic resistance,
and `rₑ` is aerodynamic resistance between the ground and canopy.
"""
@inline function compute_evaporation_ground(::PALADYNCanopyEvapotranspiration, Δq, β, rₐ, rₑ)
    # Calculate ground evaporation flux in m/s (positive upwards)
    E_gnd = β * Δq / (rₐ + rₑ)
    return E_gnd
end

"""
    $TYPEDSIGNATURES

Compute evaporation of water intercepted by the canopy from humidity gradient `Δq`, canopy saturation fraction
`f_can`, and aerodynamic resistance `rₐ`.
"""
@inline function compute_evaporation_canopy(::PALADYNCanopyEvapotranspiration, Δq, f_can, rₐ)
    # Calculate canopy evaporation flux in m/s (positive upwards)
    E_can = f_can * Δq / rₐ
    return E_can
end

# Process methods

variables(::PALADYNCanopyEvapotranspiration{NF}) where {NF} = (
    auxiliary(:evaporation_canopy, XY(); desc = "Canopy evaporation contribution to surface humidity flux", units = u"m/s"),
    auxiliary(:evaporation_ground, XY(), units = u"m/s", desc = "Ground evaporation contribution to surface humidity flux"),
    auxiliary(:transpiration, XY(), units = u"m/s", desc = "Transpiration contribution to surface humidity flux"),
    input(:skin_temperature, XY(); units = u"°C", desc = "Skin temperature"),
    input(:ground_temperature, XY(); default = NF(1), units = u"°C", desc = "Ground surface temperature"),
    input(:canopy_water_conductance, XY(); default = NF(1), units = u"m/s", desc = "Canopy stomatal conductance"), # consider direct coupling in the future
)

@propagate_inbounds function surface_humidity_flux(i, j, grid, fields, ::PALADYNCanopyEvapotranspiration)
    E_gnd = fields.evaporation_ground[i, j]
    E_can = fields.evaporation_canopy[i, j]
    T_can = fields.transpiration[i, j]
    return E_gnd + E_can + T_can
end

function compute_auxiliary!(
        state, grid,
        evapotranspiration::PALADYNCanopyEvapotranspiration,
        canopy_interception::AbstractCanopyInterception,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants,
        soil::Optional{AbstractSoil} = nothing
    )
    out = auxiliary_fields(state, evapotranspiration)
    fields = get_fields(state, evapotranspiration, canopy_interception, atmos, soil; except = out)
    launch!(grid, XY, compute_auxiliary_kernel!, out, fields, evapotranspiration, canopy_interception, atmos, constants, soil)
    return nothing
end

# Kernel functions

"""
    $TYPEDEF

Compute `transpiration`, `evaporation_ground`, and `evaporation_canopy` fluxes on `grid`
for the given scheme `evtr` and process dependencies.
"""
@propagate_inbounds function compute_evapotranspiration!(
        out, i, j, grid, fields,
        evapotranspiration::PALADYNCanopyEvapotranspiration,
        canopy_interception::AbstractCanopyInterception,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants,
        soil::Optional{AbstractSoil} = nothing
    )
    # Get inputs
    Ts = fields.skin_temperature[i, j] # skin temperature (top of canopy)
    Tg = fields.ground_temperature[i, j] # ground temperature (top snow/soil layer)
    gw_can = fields.canopy_water_conductance[i, j] # stomatal conductance

    # Compute VPD and resistance terms
    Δqs = compute_humidity_vpd(i, j, grid, fields, atmos, constants, Ts) # humidity gradient between canopy and atmosphere
    Δqg = compute_humidity_vpd(i, j, grid, fields, atmos, constants, Tg) # humidity gradient between ground and canopy
    rₐ = aerodynamic_resistance(i, j, grid, fields, atmos) # aerodynamic resistance
    rₑ = aerodynamic_resistance(i, j, grid, fields, atmos, evapotranspiration) # aerodynamic resistance between ground and canopy
    f_can = saturation_canopy_water(i, j, grid, fields, canopy_interception)
    β = ground_evaporation_resistance_factor(i, j, grid, fields, evapotranspiration.ground_resistance, soil)

    # Compute and store ET fluxes
    out.transpiration[i, j, 1] = compute_transpiration(evapotranspiration, Δqs, rₐ, gw_can)
    out.evaporation_ground[i, j, 1] = compute_evaporation_ground(evapotranspiration, Δqg, β, rₐ, rₑ)
    out.evaporation_canopy[i, j, 1] = compute_evaporation_canopy(evapotranspiration, Δqs, f_can, rₐ)
    return out
end

"""
    $TYPEDSIGNATURES

Compute the aerodynamic resistance between the ground and canopy as a function of LAI and SAI.
"""
@inline function aerodynamic_resistance(i, j, grid, fields, atmos::AbstractAtmosphere, evapotranspiration::PALADYNCanopyEvapotranspiration)
    @inbounds let LAI = fields.leaf_area_index[i, j],
            SAI = fields.SAI[i, j],
            Vₐ = windspeed(i, j, grid, fields, atmos),
            C = evtr.C_can  # drag coefficient for the canopy
        rₙ = (1 - exp(-LAI - SAI)) / (C * Vₐ)
        return rₙ
    end
end

# Ground resistance to evaporation

"""
    $TYPEDEF

Represents a spatiotemporally constant, prescribed ground evaporation resistance factor.
"""
@kwdef struct ConstantEvaporationResistanceFactor{NF} <: AbstractGroundEvaporationResistanceFactor
    "Unit interval factor that determines resistance to evaporation; zero corresponds to no evaporation"
    factor::NF = 1.0
end

ConstantEvaporationResistanceFactor(::Type{NF}; kwargs...) where {NF} = ConstantEvaporationResistanceFactor{NF}(; kwargs...)

@inline ground_evaporation_resistance_factor(i, j, grid, fields, res::ConstantEvaporationResistanceFactor, args...) = res.factor

"""
    $TYPEDEF

Implements the soil moisture limiting resistance factor of Lee and Pielke (1992),

```math
\\beta =
\\frac{1}{4} \\left[1 - \\cos\\left(π \\theta_1/\\theta_{\\text{fc}} \\right)\\right] \\quad \\text{for } \\theta_1 < \\theta_{\\text{fc}}
```
otherwise ``\\beta=1``.
"""
struct SoilMoistureResistanceFactor{NF} <: AbstractGroundEvaporationResistanceFactor end

SoilMoistureResistanceFactor(::Type{NF}) where {NF} = SoilMoistureResistanceFactor{NF}()

# Fallback implementation for interface consistency
ground_evaporation_resistance_factor(i, j, grid, fields, ::SoilMoistureResistanceFactor{NF}, args...) where {NF} = one(NF)

@inline function ground_evaporation_resistance_factor(
        i, j, grid, fields,
        ::SoilMoistureResistanceFactor{NF},
        soil::AbstractSoil
    ) where {NF}
    fgrid = get_field_grid(grid)
    strat = get_stratigraphy(soil)
    hydrology = get_hydrology(soil)
    bgc = get_biogeochemistry(soil)
    soil = soil_volume(i, j, fgrid.Nz, grid, fields, strat, hydrology, bgc)
    fc = field_capacity(get_hydraulic_properties(hydrology), soil)
    fracs = volumetric_fractions(soil)
    if fracs.water < fc
        β = (1 - cos(π * fracs.water / fc))^2 / 4
    else
        β = NF(1)
    end
    return β
end

# Forcing interface for soil hydrology

"""
    forcing(i, j, k, grid, clock, fields, evtr::AbstractEvapotranspiration, ::AbstractSoilHydrology)

Compute and return the evapotranspiration forcing for soil moisture at the given indices `i, j, k`.
The ET forcing is just the `surface_humidity_flux` rescaled by the thickness of layer `k`.
"""
@inline function forcing(i, j, k, grid, clock, fields, evapotranspiration::AbstractEvapotranspiration, ::AbstractSoilHydrology)
    let Δz = Δzᵃᵃᶜ(i, j, k, grid),
            Qh = surface_humidity_flux(i, j, grid, fields, evapotranspiration)
        ∂θ∂t = -Qh / Δz # rescale by layer thickness to get water content flux
        return ∂θ∂t * (k == grid.Nz)
    end
end

# Kernels

@kernel inbounds = true function compute_auxiliary_kernel!(out, grid, fields, evapotranspiration::AbstractEvapotranspiration, args...)
    i, j = @index(Global, NTuple)
    compute_evapotranspiration!(out, i, j, grid, fields, evapotranspiration, args...)
end
