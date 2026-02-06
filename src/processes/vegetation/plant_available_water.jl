"""
    $TYPEDEF

Implementation of vegetation water availability (a.k.a "plant available water")
that computes the wilting fraction

```math
W_i = \\min(\\frac{\\theta_{\\text{w},i} - \\theta_{\\text{wp},i}}{\\theta_{\\text{fc},i} - \\theta_{\\text{wp},i}}}, 1)
```

where ``\\theta_{\\text{w},i}`` is the volumetric water content of the i'th soil layer, ``\\theta_{\\text{fc},i}``
is the "field capacity", and ``\\theta_{\\text{wp},i}`` is the "wilting point". The water availability

Properties:
$TYPEDFIELDS
"""
@kwdef struct FieldCapacityLimitedPAW{NF} <: AbstractPlantAvailableWater{NF} end

FieldCapacityLimitedPAW(::Type{NF} = Float32) where {NF} = FieldCapacityLimitedPAW{NF}()

variables(paw::FieldCapacityLimitedPAW{NF}) where {NF} = (
    auxiliary(:plant_available_water, XYZ(), desc="Fraction of soil water available for plant root water uptake"),
    auxiliary(:SMLF, XY(), soil_moisture_limiting_factor, paw), # soil moisture limiting factor
    input(:root_fraction, XYZ(), desc="Fraction of roots in each soil layer")
)

"""
    $TYPEDSIGNATURES

Field constructor for the soil moisture limiting factor. Returns a derived `Field` that calculates
the integral of `W(z) * r(z)` where `W` is the water availability coefficient and `r` is the root fraction.
"""
function soil_moisture_limiting_factor(::FieldCapacityLimitedPAW, grid, clock, fields)
    Δz = zspacings(get_field_grid(grid), Center(), Center(), Center())
    PAW = Integral(fields.plant_available_water * fields.root_fraction / Δz, dims=3)
    return Field(PAW)
end

# Process methods

function compute_auxiliary!(
    state, grid,
    paw::FieldCapacityLimitedPAW,
    soil::AbstractSoil,
    args...
)
    out = auxiliary_fields(state, paw)
    fields = get_fields(state, paw, soil; except = out)
    launch!(grid, XYZ, compute_paw_kernel!, out, fields, paw, soil)
    # compute the derived soil moisture limiting factor field
    compute!(state.SMLF)
end

# Kernel functions

@propagate_inbounds function compute_paw(
    i, j, k, grid, fields,
    paw::FieldCapacityLimitedPAW{NF},
    soil::AbstractSoil
) where {NF}
    # Unpack soil processes
    strat = get_stratigraphy(soil)
    hydrology = get_hydrology(soil)
    bgc = get_biogeochemistry(soil)
    # Compute soil composition and hydraulic properties
    vol = soil_volume(i, j, k, grid, fields, strat, hydrology, bgc)
    θfc = field_capacity(hydrology.hydraulic_properties, vol.solid.texture)
    θwp = wilting_point(hydrology.hydraulic_properties, vol.solid.texture)
    # Compute liquid water content
    fracs = volumetric_fractions(vol)
    θw = fracs.water
    # Compute PAW
    return max(min(NF(1), (θw - θwp) / (θfc - θwp)), NF(0))
end

@propagate_inbounds function compute_paw!(
    out, i, j, k, grid, fields,
    paw::FieldCapacityLimitedPAW{NF},
    soil::AbstractSoil,
    args...
) where {NF}
    PAW = compute_paw(i, j, k, grid, fields, paw, soil)
    out.plant_available_water[i, j, k] = PAW
    return out
end

# Kernels

@kernel inbounds=true function compute_paw_kernel!(out, grid, fields, paw::AbstractPlantAvailableWater, args...)
    i, j, k = @index(Global, NTuple)
    compute_paw!(out, i, j, k, grid, fields, paw, args...)
end
