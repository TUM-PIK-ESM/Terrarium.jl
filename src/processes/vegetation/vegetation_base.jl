# Generic methods for all AbstractVegetation implementations

# TODO: Remove once dedicated vegetation surface parameterizations are added
# Also may be a good reason to rename VegetationCarbonCycle to NaturalVegetation or similar
@propagate_inbounds compute_albedo(i, j, grid, fields, veg::AbstractVegetation{NF}) where {NF} = veg.traits.albedo

# TODO: will need to change once PFTs are added
@propagate_inbounds vegetation_area_fraction(i, j, grid, fields, ::Nothing) = zero(eltype(grid))
@propagate_inbounds function vegetation_area_fraction(i, j, grid, fields, veg::AbstractVegetation)
    if isnothing(veg.vegetation_dynamics)
        LAI = fields.leaf_area_index[i, j, 1]
        LAI_max = maximum_leaf_area_index(i, j, grid, fields, veg.traits)
        f_veg = LAI / LAI_max
        return f_veg
    else
        return vegetation_area_fraction(i, j, grid, fields, veg.vegetation_dynamics)
    end
end
