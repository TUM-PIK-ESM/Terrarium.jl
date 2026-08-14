# Generic methods for all AbstractVegetation implementations

# TODO: Remove once dedicated vegetation surface parameterizations are added
# Also may be a good reason to rename VegetationCarbonCycle to NaturalVegetation or similar
@inline compute_albedo(i, j, grid, fields, veg::AbstractVegetation{NF}) where {NF} = veg.traits.albedo

@inline compute_emissivity(i, j, grid, fields, veg::AbstractVegetation{NF}) where {NF} = veg.traits.emissivity

# TODO: will need to change once PFTs are added
@inline vegetation_area_fraction(i, j, grid, fields, ::Union{Nothing, AbstractVegetation}) = zero(eltype(grid))
