@kwdef struct PhenologyModel{NF} <: AbstractPhenology
    # TODO add phenology parameters

end

variables(phenol::PhenologyModel) = (
    auxiliary(:LAI, XY()), # Leaf Area Index (LAI)
)

function compute_auxiliary!(idx, state, model, phenol::PhenologyModel)
    i, j = idx

    # TODO add phenology implementation from Paladyn
    # Compute f_deciduous, a factor for smooth transition between evergreen and deciduous
    # For now, set f_deciduous to 0.0 (evergreen PFT)
    f_deciduous = 0.0

    # Compute phen 
    # For now, set phen to 1.0 (full leaf-out, evergreen PFT)
    phen = 1.0

    # Compute LAI
    state.LAI[i, j] = (f_deciduous * phen + (1 - f_deciduous)) * state.LAI_b[i, j]
end
   