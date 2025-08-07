# TODO maybe change the name later, if the PALADYN phenology approach has a more specific name

@kwdef struct PALADYNPhenology{NF} <: AbstractPhenology
    # TODO add phenology parameters

end

variables(phenol::PALADYNPhenology) = (
    auxiliary(:LAI_b, XY()), # Balanced Leaf Area Index 
    auxiliary(:LAI, XY()), # Leaf Area Index 
)

@inline function compute_auxiliary!(idx, state, model::AbstractVegetationModel, phenol::PALADYNPhenology{NF}) where NF
    i, j = idx

    # TODO add phenology implementation from PALADYN
    # Compute f_deciduous, a factor for smooth transition between evergreen and deciduous
    # For now, set f_deciduous to 0.0 (evergreen PFT)
    f_deciduous = NF(0.0)

    # Compute phen 
    # For now, set phen to 1.0 (full leaf-out, evergreen PFT)
    phen = NF(1.0)

    # Compute LAI
    state.LAI[i, j] = (f_deciduous * phen + (NF(1.0) - f_deciduous)) * state.LAI_b[i, j]
end
   