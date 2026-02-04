"""
    $TYPEDEF

Model for natural (unmanaged) vegetation processes for a single plant functional type (PFT).
Multiple PFTs can be later handled with a `TiledVegetationModel` type that composes multiple
`VegetationModel`s with different parameters for each PFT.

Properties:
$TYPEDFIELDS
"""
@kwdef struct VegetationModel{
    NF,
    Vegetation<:AbstractVegetation{NF},
    Atmosphere<:AbstractAtmosphere{NF},
    GridType<:AbstractLandGrid{NF},
    Constants<:PhysicalConstants{NF},
    Initializer<:AbstractInitializer,
} <: AbstractVegetationModel{NF, GridType}
    "Spatial grid type"
    grid::GridType

    "Atmospheric input configuration"
    atmosphere::Atmosphere = PrescribedAtmosphere(eltype(grid))

    "Vegetation processes"
    vegetation::Vegetation = VegetationCarbon(eltype(grid))

    "Physical constants"
    constants::Constants = PhysicalConstants(eltype(grid))

    "State variable initializer"
    initializer::Initializer = DefaultInitializer()
end

function compute_auxiliary!(state, model::VegetationModel)
    # Unpack vegetation model
    (; grid, atmosphere, vegetation, constants) = model
    # Compute auxiliary variables for coupled vegetation processes
    compute_auxiliary!(state, grid, vegetation, atmosphere, constants)
end

function compute_tendencies!(state, model::VegetationModel)
    # Unpack vegetation model
    (; grid, vegetation, constants) = model
    # Compute auxiliary variables for coupled vegetation processes
    compute_tendencies!(state, grid, vegetation, constants)
end

