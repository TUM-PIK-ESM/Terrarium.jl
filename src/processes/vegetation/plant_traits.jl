"""
    $TYPEDEF

Plant-specific parameters for a single functional type. This is a temporary, partial construct
that will be soon be replaced will full support for PFTs and other trait parameterizations.
"""
@parameterized @kwdef struct PlantTraits{NF}
    "Specific leaf area ([kattgeTRYGlobalDatabase2011](@cite)). PFT specific."
    @param specific_leaf_area::NF = 10.0 (units = u"m^2/kg", bounds = Positive) # Value for Needleleaf tree PFT

    "Allometric coefficient, modified from [coxDescriptionTRIFFIDDynamic2001](@cite) to account for bwl=1. PFT specific."
    @param awl::NF = 2.0 (units = u"kg/m^2", bounds = Positive) # Value for Needleleaf tree PFT

    "Minimum Leaf Area Index modified from [clarkJointUKLand2011](@cite). PFT specific."
    @param minimum_leaf_area_index::NF = 1.0 (bounds = Positive,) # Value for Needleleaf tree PFT

    "Maximum Leaf Area Index modified from [clarkJointUKLand2011](@cite). PFT specific."
    @param maximum_leaf_area_index::NF = 6.0 (bounds = Positive,) # Value for Needleleaf tree PFT

    "Extinction coefficient for radiation through vegetation"
    @param extinction_coefficient::NF = 0.5 (bounds = Positive,)

    "Snow-free canopy albedo for broadband shortwave radiation"
    @param albedo::NF = 0.15 (bounds = UnitInterval,)

    "Snow-free canopy emissivity for broadband longwave radiation (default value from PALADYN: 0.96)"
    @param emissivity::NF = 0.96 (bounds = UnitInterval,)
end

PlantTraits(::Type{NF}; kwargs...) where {NF} = PlantTraits{NF}(; kwargs...)

@inline specific_leaf_area(i, j, grid, fields, traits::PlantTraits) = traits.specific_leaf_area

@inline allometric_coefficient(i, j, grid, fields, traits::PlantTraits) = traits.awl

@inline minimum_leaf_area_index(i, j, grid, fields, traits::PlantTraits) = traits.minimum_leaf_area_index

@inline maximum_leaf_area_index(i, j, grid, fields, traits::PlantTraits) = traits.maximum_leaf_area_index

@inline extinction_coefficient(i, j, grid, fields, traits::PlantTraits) = traits.extinction_coefficient

@inline albedo(i, j, grid, fields, traits::PlantTraits) = traits.albedo
