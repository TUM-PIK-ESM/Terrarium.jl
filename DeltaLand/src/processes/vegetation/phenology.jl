# Note: maybe change the name later, if the PALADYN phenology approach has a more specific name
"""
    $TYPEDEF

Vegetation phenology implementation from PALADYN (Willeit 2016).

Authors: Maha Badri and Matteo Willeit

Properties:
$TYPEDFIELDS
"""
@kwdef struct PALADYNPhenology{NF} <: AbstractPhenology
    # TODO add phenology parameters
end

variables(phenol::PALADYNPhenology) = (
    auxiliary(:LAI_b, XY()), # Balanced Leaf Area Index 
    auxiliary(:phen, XY()), # Phenology factor
    auxiliary(:LAI, XY()), # Leaf Area Index 
)

function compute_auxiliary!(state, model, phenol::PALADYNPhenology)
    grid = get_grid(model)
    launch!(grid, :xy, compute_auxiliary_kernel!, state, phenol)
end

@kernel function compute_auxiliary_kernel!(state, phenol::PALADYNPhenology)
    i, j = @index(Global, NTuple)

    # TODO add phenology implementation from PALADYN
    # Compute f_deciduous, a factor for smooth transition between evergreen and deciduous
    # For now, set f_deciduous to 0.0 (evergreen PFT)
    f_deciduous = zero(NF)

    # Compute phen 
    # For now, set phen to 1.0 (full leaf-out, evergreen PFT)
    phen = NF(1.0)

    # Compute LAI
    state.LAI[i, j] = (f_deciduous * phen + (NF(1.0) - f_deciduous)) * state.LAI_b[i, j]
end
   