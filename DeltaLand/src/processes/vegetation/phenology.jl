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

"""
    $SIGNATURES

Computes `f_deciduous`, a factor for smooth transition between evergreen and deciduous.
"""
@inline function compute_f_deciduous(phenol::PALADYNPhenology{NF}) where NF
    # TODO add phenology implementation from PALADYN
    # For now, set f_deciduous to 0.0 (evergreen PFT)
    f_deciduous = zero(NF)

    return f_deciduous
end

"""
    $SIGNATURES

Computes `phen`, the phenology factor.
"""
@inline function compute_phen(phenol::PALADYNPhenology{NF}) where NF
    # TODO add phenology implementation from PALADYN
    # For now, set phen to 1.0 (full leaf-out, evergreen PFT)
    phen = NF(1.0)

    return phen
end


function compute_auxiliary!(state, model, phenol::PALADYNPhenology)
    grid = get_grid(model)
    launch!(grid, :xy, compute_auxiliary_kernel!, state, phenol)
end

@kernel function compute_auxiliary_kernel!(state, phenol::PALADYNPhenology)
    i, j = @index(Global, NTuple)

    # Compute f_deciduous
    f_deciduous = compute_f_deciduous(phenol)

    # Compute phen 
    phen = compute_phen(phenol)

    # Compute LAI
    state.LAI[i, j] = (f_deciduous * phen + (NF(1.0) - f_deciduous)) * state.LAI_b[i, j]
end
   