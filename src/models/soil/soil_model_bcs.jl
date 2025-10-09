"""
    $TYPEDEF

Boundary condition type for soil models that provides boundary conditions for each
of the relevant internal processes, i.e. energy, hydrology, ...
"""
struct SoilBC{EnergyBC, WaterBC} <: AbstractBoundaryConditions
    "Boundary condition for the soil energy balance"
    energy::EnergyBC

    "Boundary condition for the soil water balance"
    hydrology::WaterBC
end

variables(bc::SoilBC) = tuplejoin(variables(bc.hydrology), variables(bc.energy))

function compute_auxiliary!(state, model, bc::SoilBC)
    compute_auxiliary!(state, model, bc.hydrology)
    compute_auxiliary!(state, model, bc.energy)
end

function get_field_boundary_conditions(bcs::SoilBC, grid::AbstractLandGrid)
    water_bc = get_field_boundary_conditions(bcs.hydrology, grid)
    energy_bc = get_field_boundary_conditions(bcs.energy, grid)
    return merge(water_bc, energy_bc)
end

"""
Alias for `PrescribedFlux` with name `Q_g` representing the net ground heat flux at the soil surface.
"""
GroundHeatFlux(init) = PrescribedFlux(:U, Input(:Q_g, init, units=u"W/m^2"))

"""
Alias for `PrescribedFlux` with name `Q_geo` representing the geothermal heat flux at the bottom
boundary of the soil column.
"""
GeothermalHeatFlux(init) = PrescribedFlux(:U, Input(:Q_geo, init, units=u"W/m^2"))

"""
Alias for `PrescribedFlux` with name `Q_inf` representing liquid water infiltration at the soil surface.
"""
InfiltrationFlux(init) = PrescribedFlux(:saturation_water_ice, Input(:Q_inf, init, units=u"m/s"))

"""
Alias for `NoFlux` representing a zero-flux boundary condition for water flow (prognostic variable `saturation_water_ice`).
"""
ImpermeableBoundary() = NoFlux(:saturation_water_ice)

"""
Alias for `PrescribedGradient` representing a Neumann-type zero pressure gradient at the bottom of the soil
column, thereby allowing free drainage of water.
"""
FreeDrainage(::Type{NF}) where {NF} = PrescribedGradient(:water_potential, zero(NF))

"""
Alias for `ColumnBoundaryConditions` with defaults suitable for `SoilModel`s.
"""
SoilBoundaryConditions(::Type{NF}; top=default_soil_upperbc(NF), bottom=default_soil_lowerbc(NF)) where {NF} = ColumnBoundaryConditions(; top, bottom)

default_soil_upperbc(::Type{NF}, energy=GroundHeatFlux(zero(NF)), hydrology=InfiltrationFlux(zero(NF))) where {NF} = SoilBC(energy, hydrology)
default_soil_lowerbc(::Type{NF}, energy=GeothermalHeatFlux(zero(NF)), hydrology=ImpermeableBoundary()) where {NF} = SoilBC(energy, hydrology)
