@parameterized @kwdef struct PrescribedVegetation{NF, Phenology, StomatalConductance, RootDistribution, PAW} <: AbstractVegetation{NF}
    "Phenology scheme"
    @component phenology::Phenology

    "Stomatal conductance"
    @component stomatal_conductance::StomatalConductance

    "Vegetation root distribution"
    @component root_distribution::RootDistribution

    "Plant available water"
    @component plant_available_water::PAW

    "Plant physical traits"
    @component traits::PlantTraits{NF}
end

function PrescribedVegetation(
        ::Type{NF};
        phenology = PrescribedPhenology(NF),
        stomatal_conductance = MedlynStomatalConductance(NF),
        root_distribution = StaticExponentialRootDistribution(NF),
        plant_available_water = FieldCapacityLimitedPAW(NF),
        traits = PlantTraits(NF)
    ) where {NF}
    return PrescribedVegetation(; phenology, stomatal_conductance, root_distribution, plant_available_water, traits)
end

"""
    $TYPEDSIGNATURES

Compute auxiliary variables for all vegetation component processes based on the given
atmospheric inputs defined by `atmos` and (optionally) `soil` state. If `soil = nothing`,
stress factors due to soil temperature and moisture availability will be ignored.
"""
function compute_auxiliary!(
        state, grid,
        veg::PrescribedVegetation,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        soil::Optional{AbstractSoil} = nothing,
        args...
    )
    # Compute auxiliary variables for each component
    # Roots: need soil state and computes root_fraction
    compute_auxiliary!(state, grid, veg.root_distribution, soil)

    # PAW: needs soil saturation profile and computes soil_moisture_limiting_factor
    compute_auxiliary!(state, grid, veg.plant_available_water, soil)

    # Phenology: needs LAI(t) and air temperature(t) as inputs and computes phen(t)
    compute_auxiliary!(state, grid, veg.phenology, atmos)

    # Stomatal conductance: needs atm. inputs(t) and computes λc(t)
    compute_auxiliary!(state, grid, veg.stomatal_conductance, veg.traits, constants, atmos)
    return nothing
end

# No tendencies for prescribed vegetation
@inline compute_tendencies!(state, grid, veg::PrescribedVegetation, args...) = nothing
