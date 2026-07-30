"""
    $TYPEDEF

Represents a generic coupling of vegetation carbon processes.
"""
@kwdef struct VegetationCarbon{
        NF,
        Photosynthesis <: AbstractPhotosynthesis,
        StomatalConductance <: AbstractStomatalConductance,
        AutotrophicRespiration <: AbstractAutotrophicRespiration,
        Phenology <: AbstractPhenology,
        CarbonDynamics <: AbstractVegetationCarbonDynamics,
        VegetationDynamics <: Optional{AbstractVegetationDynamics},
        RootDistribution <: Optional{AbstractRootDistribution},
        PAW <: Optional{AbstractPlantAvailableWater},
    } <: AbstractVegetation{NF}
    "Photosynthesis scheme"
    photosynthesis::Photosynthesis # not prognostic

    "Stomatal conductance scheme"
    stomatal_conductance::StomatalConductance # not prognostic

    "Autotrophic respiration scheme"
    autotrophic_respiration::AutotrophicRespiration # not prognostic

    "Phenology scheme"
    phenology::Phenology # not prognostic

    "Vegetation carbon pool dynamics"
    carbon_dynamics::CarbonDynamics # prognostic

    "Vegetation population density or coverage fraction dynamics"
    vegetation_dynamics::VegetationDynamics # prognostic

    "Vegetation root distribution"
    root_distribution::RootDistribution

    "Plant available water"
    plant_available_water::PAW
end

# No component process is pinned to `NF` — any of them may hold a promoted (e.g. traced) parameter —
# so `NF` is not derivable from the field types and must be given. `constructorof` rebuilds through
# this constructor so `setproperties` / `Flatten.reconstruct` / `Adapt` retain `NF`.
VegetationCarbon{NF}(
    photosynthesis::AbstractPhotosynthesis,
    stomatal_conductance::AbstractStomatalConductance,
    autotrophic_respiration::AbstractAutotrophicRespiration,
    phenology::AbstractPhenology,
    carbon_dynamics::AbstractVegetationCarbonDynamics,
    vegetation_dynamics::Optional{AbstractVegetationDynamics},
    root_distribution::Optional{AbstractRootDistribution},
    plant_available_water::Optional{AbstractPlantAvailableWater},
) where {NF} = VegetationCarbon{
    NF,
    typeof(photosynthesis),
    typeof(stomatal_conductance),
    typeof(autotrophic_respiration),
    typeof(phenology),
    typeof(carbon_dynamics),
    typeof(vegetation_dynamics),
    typeof(root_distribution),
    typeof(plant_available_water),
}(
    photosynthesis,
    stomatal_conductance,
    autotrophic_respiration,
    phenology,
    carbon_dynamics,
    vegetation_dynamics,
    root_distribution,
    plant_available_water,
)

ConstructionBase.constructorof(::Type{<:VegetationCarbon{NF}}) where {NF} = VegetationCarbon{NF}

function VegetationCarbon(
        ::Type{NF};
        photosynthesis = LUEPhotosynthesis(NF),
        stomatal_conductance = MedlynStomatalConductance(NF),
        autotrophic_respiration = PALADYNAutotrophicRespiration(NF),
        phenology = PALADYNPhenology(NF),
        carbon_dynamics = PALADYNCarbonDynamics(NF),
        vegetation_dynamics = PALADYNVegetationDynamics(NF),
        root_distribution = StaticExponentialRootDistribution(NF),
        plant_available_water = FieldCapacityLimitedPAW(NF)
    ) where {NF}
    return VegetationCarbon{NF}(
        photosynthesis,
        stomatal_conductance,
        autotrophic_respiration,
        phenology,
        carbon_dynamics,
        vegetation_dynamics,
        root_distribution,
        plant_available_water
    )
end

"""
    $TYPEDSIGNATURES

Compute auxiliary variables for all vegetation component processes based on the given
atmospheric inputs defined by `atmos` and (optionally) `soil` state. If `soil = nothing`,
stress factors due to soil temperature and moisture availability will be ignored.
"""
function compute_auxiliary!(
        state, grid,
        veg::VegetationCarbon,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        soil::Optional{AbstractSoil} = nothing,
        args...
    )
    # Compute auxiliary variables for each component
    # PAW: needs soil saturation profile and computes soil_moisture_limiting_factor
    compute_auxiliary!(state, grid, veg.plant_available_water, soil)

    # Veg. carbon dynamics: needs C_veg(t-1) and computes LAI_b(t-1)
    compute_auxiliary!(state, grid, veg.carbon_dynamics)

    # Phenology: needs LAI_b(t-1) and computes LAI(t-1) and phen(t-1)
    compute_auxiliary!(state, grid, veg.phenology)

    # Stomatal conductance: needs atm. inputs(t) and computes λc(t)
    # TODO: Note the (implicit) circular dependency between photosynthesis and stomatal conductance;
    # can this be refactored?
    compute_auxiliary!(state, grid, veg.stomatal_conductance, veg.photosynthesis, constants, atmos)

    # Photosynthesis: needs atm. inputs(t), λc(t), LAI(t-1), and computes Rd(t) and GPP(t)
    compute_auxiliary!(state, grid, veg.photosynthesis, veg.stomatal_conductance, constants, atmos)

    # Autotrophic respiration: needs atm. inputs(t), GPP(t), Rd(t), C_veg(t-1), phen(t-1) and computes Ra(t) and NPP(t)
    compute_auxiliary!(state, grid, veg.autotrophic_respiration, veg.carbon_dynamics, atmos)

    # Note: vegetation_dynamics compute_auxiliary! does nothing for now
    compute_auxiliary!(state, grid, veg.vegetation_dynamics)
    return nothing
end

"""
    $TYPEDSIGNATURES

Compute tendencies for carbon and vegetation dynamics.
"""
function compute_tendencies!(state, grid, veg::VegetationCarbon, args...)
    # Needs NPP(t), C_veg(t-1), LAI_b(t-1) and computes tendency for C_veg
    compute_tendencies!(state, grid, veg.carbon_dynamics)

    # Needs NPP(t), C_veg(t-1), LAI_b(t-1), ν(t-1) and computes tendency for ν
    compute_tendencies!(state, grid, veg.vegetation_dynamics, veg.carbon_dynamics)

    return nothing
end
