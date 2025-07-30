abstract type AbstractStateVariables end

# temporary solution for holding state variables
@kwdef struct StateVariables{
    prognames, tendnames, auxnames, bcnames,
    ProgVars, TendVars, AuxVars, BCVars
} <: AbstractStateVariables
    prognostic::NamedTuple{prognames, ProgVars} = (;)
    tendencies::NamedTuple{tendnames, TendVars} = (;)
    auxiliary::NamedTuple{auxnames, AuxVars} = (;)
    boundary_vars::NamedTuple{bcnames, BCVars} = (;)
end

Base.propertynames(
    vars::StateVariables{prognames, tendnames, auxnames, bcnames}
) where {prognames, tendnames, auxnames, bcnames} = (
    prognames...,
    tendnames...,
    auxnames...,
    bcnames...,
    fieldnames(typeof(vars))...,
)

function Base.getproperty(
    vars::StateVariables{prognames, tendnames, auxnames, bcnames},
    name::Symbol
) where {prognames, tendnames, auxnames, bcnames}
    # forward getproperty calls to variable groups
    if name ∈ prognames
        return getproperty(getfield(vars, :prognostic), name)
    elseif name ∈ tendnames
        return getproperty(getfield(vars, :tendencies), name)
    elseif name ∈ auxnames
        return getproperty(getfield(vars, :auxiliary), name)
    elseif name ∈ bcnames
        return getproperty(getfield(vars, :boundary_vars), name)
    else
        return getfield(vars, name)
    end
end

function reset_tendencies!(state::StateVariables)
    for var in getfield(state, :tendencies)
        set!(var, zero(eltype(var)))
    end
end
