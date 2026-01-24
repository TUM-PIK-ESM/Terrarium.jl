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

PALADYNPhenology(::Type{NF}; kwargs...) where {NF} = PALADYNPhenology{NF}(; kwargs...)

variables(::PALADYNPhenology) = (
    auxiliary(:phen, XY()), # Phenology factor [-]
    auxiliary(:LAI, XY()), # Leaf Area Index [m²/m²]
    input(:LAI_b, XY()), # Balanced Leaf Area Index [m²/m²]
)

"""
    $SIGNATURES

Computes `f_deciduous`, a factor for smooth transition between evergreen and deciduous [-].
"""
@inline function compute_f_deciduous(phenol::PALADYNPhenology{NF}) where NF
    # TODO add phenology implementation from PALADYN
    # For now, set f_deciduous to 0.0 (evergreen PFT)
    f_deciduous = zero(NF)

    return f_deciduous
end

"""
    $SIGNATURES

Computes `phen`, the phenology factor [-].
"""
@inline function compute_phen(phenol::PALADYNPhenology{NF}) where NF
    # TODO add phenology implementation from PALADYN
    # For now, set phen to 1.0 (full leaf-out, evergreen PFT)
    phen = NF(1.0)

    return phen
end

"""
    $SIGNATURES

Computes `LAI`, based on the balanced Leaf Area Index `LAI_b`:
"""
@inline function compute_LAI(phenol::PALADYNPhenology{NF}, LAI_b::NF) where NF
 # Compute f_deciduous
    f_deciduous = compute_f_deciduous(phenol)

    # Compute phen 
    phen = compute_phen(phenol)

    # Compute LAI
    LAI = (f_deciduous * phen + (NF(1.0) - f_deciduous)) * LAI_b

    return LAI
end

function compute_auxiliary!(state, model, phenol::PALADYNPhenology)
    grid = get_grid(model)
    launch!(grid, XY, compute_auxiliary_kernel!, state, phenol)
end

@kernel function compute_auxiliary_kernel!(state, grid, phenol::PALADYNPhenology)
    i, j = @index(Global, NTuple)

    # Get input
    LAI_b = state.LAI_b[i, j]

    # Compute phen
    phen = compute_phen(phenol)

    # Compute LAI
    LAI = compute_LAI(phenol, LAI_b)

    # Store results
    state.phen[i, j] = phen
    state.LAI[i, j] = LAI
end
   