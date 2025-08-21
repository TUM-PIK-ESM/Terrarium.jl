"""
Alias for `PrescribedFlux` with name `Q_g` representing the 
"""
GroundHeatFlux(value) = PrescribedFlux(:Q_g, value, XY())

"""
Alias for `PrescribedFlux` with name `Q_geo` representing the geothermal heat flux
at the bottom of a deep soil column.
"""
GeothermalHeatFlux(value) = PrescribedFlux(:Q_geo, value, XY())

"""
Alias for `PrescribedFlux` with name `Q_inf` representing topsoil infiltration.
"""
InfiltrationFlux(value) = PrescribedFlux(:Q_inf, value, XY())

"""
Alias for `PrescribedFlux` with name `Q_out` representing water drainage at the bottom of
a soil column. When this flux is set to zero, it corresponds to an impermeable boundary.
"""
FreeDrainage(value) = PrescribedFlux(:Q_out, value, XY())
ImpermeableBoundary(::Type{NF}) where {NF<:AbstractFloat} = PrescribedFlux(:Q_out, zero(NF), XY())

# Convenience constructors for common soil upper and lower boundary conditions
SoilBoundaryConditions(
    grid::AbstractLandGrid;
    top = SoilUpperBoundaryConditions(grid),
    bottom = SoilLowerBoundaryConditions(grid),
) = VerticalBoundaryConditions(; top, bottom)

SoilUpperBoundaryConditions(
    grid::AbstractLandGrid;
    energy = GroundHeatFlux(zero(eltype(grid))),
    hydrology = InfiltrationFlux(zero(eltype(grid)))
) = SoilBoundaryCondition(energy, hydrology)

SoilLowerBoundaryConditions(
    grid::AbstractLandGrid;
    energy = GeothermalHeatFlux(zero(eltype(grid))),
    hydrology = ImpermeableBoundary(eltype(grid)),
) = SoilBoundaryCondition(energy, hydrology)

"""
    $TYPEDEF

Boundary condition type for soil models that represents the conditions applied at one
of the boundaries of the soil domain.
"""
struct SoilBoundaryCondition{EnergyBC, WaterBC} <: AbstractBoundaryConditions
    "Boundary conditions for the soil energy balance"
    energy::EnergyBC

    "Boundary conditions for the soil water balance"
    hydrology::WaterBC
end

variables(bc::SoilBoundaryCondition) = tuplejoin(variables(bc.hydrology), variables(bc.energy))

function compute_auxiliary!(state, model, bc::SoilBoundaryCondition)
    compute_auxiliary!(state, model, bc.hydrology)
    compute_auxiliary!(state, model, bc.energy)
end

function get_field_boundary_conditions(bcs::SoilBoundaryCondition, grid::AbstractLandGrid)
    water_bc = get_field_boundary_conditions(bcs.hydrology, grid)
    energy_bc = get_field_boundary_conditions(bcs.energy, grid)
    return merge(water_bc, energy_bc)
end
