"""
    $TYPEDEF

Plant-specific parameters for a single functional type. This is a temporary, partial construct
that will be soon be replaced will full support for PFTs and other trait parameterizations.
"""
@parameterized @kwdef struct PlantTraits{NF}
    "Specific leaf area ([kattgeTRYGlobalDatabase2011](@cite)). PFT specific."
    @param SLA::NF = 10.0 (units = u"m^2/kg", bounds = Positive) # Value for Needleleaf tree PFT

    "Allometric coefficient, modified from [coxDescriptionTRIFFIDDynamic2001](@cite) to account for bwl=1. PFT specific."
    @param awl::NF = 2.0 (units = u"kg/m^2", bounds = Positive) # Value for Needleleaf tree PFT

    "Minimum Leaf Area Index modified from [clarkJointUKLand2011](@cite). PFT specific."
    @param LAI_min::NF = 1.0 (bounds = Positive,) # Value for Needleleaf tree PFT

    "Maximum Leaf Area Index modified from [clarkJointUKLand2011](@cite). PFT specific."
    @param LAI_max::NF = 6.0 (bounds = Positive,) # Value for Needleleaf tree PFT

    "Extinction coefficient for radiation through vegetation"
    @param k_ext::NF = 0.5 (bounds = Positive,)

    "Snow-free visble canopy albedo for shortwave radiation"
    @param albedo::NF = 0.02 (bounds = UnitInterval,)
end

PlantTraits(::Type{NF}) where {NF} = PlantTraits{NF}()

@inline specific_leaf_area(i, j, grid, fields, traits::PlantTraits) = traits.SLA

@inline allometric_coefficient(i, j, grid, fields, traits::PlantTraits) = traits.awl

@inline minimum_leaf_area_index(i, j, grid, fields, traits::PlantTraits) = traits.LAI_min

@inline maximum_leaf_area_index(i, j, grid, fields, traits::PlantTraits) = traits.LAI_max
