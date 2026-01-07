"""
    $TYPEDEF

Implementation of vegetation water availability (a.k.a "plant available water")
that computes the wilting fraction

```math
W_i = \\min(\\frac{\\theta_{\\text{w},i} - \\theta_{\\text{fc},i}}{\\theta_{\\text{fc},i} - \\theta_{\\text{wp},i}}}, 1)
```

where ``\\theta_{\\text{w},i}`` is the volumetric water content of the i'th soil layer, ``\\theta_{\\text{fc},i}``
is the "field capacity", and ``\\theta_{\\text{wp},i}`` is the "wilting point". The water availability

Properties:
$TYPEDFIELDS
"""
@kwdef struct FieldCapacityLimitedPAW{NF} <: AbstractPlantAvailableWater end

FieldCapacityLimitedPAW(::Type{NF} = Float32) where {NF} = FieldCapacityLimitedPAW{NF}()

variables(paw::FieldCapacityLimitedPAW) = (
    auxiliary(:plant_available_water, XYZ(), desc="Fraction of soil water available for plant root water uptake"),
    auxiliary(:SMLF, XY(), soil_moisture_limiting_factor, paw), # soil moisture limiting factor
    input(:root_fraction, XYZ())
)

"""
Field constructor for the soil moisture limiting factor. Returns an `Integral` of
`W(z) * r(z)` where `W` is the water availability coefficient and `r` is the root fraction.
"""
function soil_moisture_limiting_factor(paw::FieldCapacityLimitedPAW, grid, clock, fields)
    return Integral(fields.plant_available_water * fields.root_fraction)
end

function compute_auxiliary!(state, model, paw::FieldCapacityLimitedPAW)
    grid = get_grid(model)
    hydrology = get_soil_hydrology(model)
    strat = get_straigraphy(model)
    energy = get_soil_energy_balance(model)
    bgc = get_soil_biogeochemistry(model)
    launch!(state, grid, :xyz, compute_paw_kernel!, paw, hydrology, strat, energy, bgc)
end

@kernel function compute_paw_kernel!(
    state, grid,
    paw::FieldCapacityLimitedPAW,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    energy::AbstractSoilEnergyBalance,
    bgc::AbstractSoilBiogeochemistry
)
    i, j, k = @index(Global, NTuple)

    @inbounds let soil = soil_volume(i, j, k, state, grid, strat, energy, hydrology, bgc),
                  θfc = field_capacity(hydrology.hydraulic_properties, soil.solid.texture),
                  θwp = wilting_point(hydrology.hydraulic_properties, soil.solid.texture),
                  vol = volumetric_fractions(soil[i, j, k]),
                  θw = vol.water;
        # compute PAW
        state.plant_available_water[i, j, k] = (θw - θfc) / (θfc - θwp)
    end
end
