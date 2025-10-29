"""
    $TYPEDEF

Boundary condition type for soil models that provides boundary conditions for each
of the relevant internal processes, i.e. energy, hydrology, ...
"""
@kwdef struct SoilBC{EnergyBC, WaterBC} <: AbstractBoundaryConditions
    "Boundary condition for the soil energy balance"
    energy::EnergyBC = GroundHeatFlux()

    "Boundary condition for the soil water balance"
    hydrology::WaterBC = ImpermeableBoundary()
end

variables(bc::SoilBC) = tuplejoin(variables(bc.hydrology), variables(bc.energy))

function compute_auxiliary!(state, model, bc::SoilBC)
    compute_auxiliary!(state, model, bc.hydrology)
    compute_auxiliary!(state, model, bc.energy)
end

function compute_tendencies!(state, model, bc::SoilBC)
    compute_tendencies!(state, model, bc.hydrology)
    compute_tendencies!(state, model, bc.energy)
end

function get_field_boundary_conditions(bcs::SoilBC, grid::AbstractLandGrid)
    water_bc = get_field_boundary_conditions(bcs.hydrology, grid)
    energy_bc = get_field_boundary_conditions(bcs.energy, grid)
    return merge(water_bc, energy_bc)
end

# Energy BCs

"""
Alias for `PrescribedFlux` on `internal_energy` with name `ground_heat_flux` representing the net ground heat flux at the soil surface.
"""
GroundHeatFlux(init=nothing) = PrescribedFlux(:internal_energy, Input(:ground_heat_flux, init, units=u"W/m^2"))

"""
Alias for `PrescribedFlux` on `internal_energy` with name `geothermal_heat_flux` representing the geothermal heat flux at the bottom
boundary of the soil column.
"""
GeothermalHeatFlux(init=nothing) = PrescribedFlux(:internal_energy, Input(:geothermal_heat_flux, init, units=u"W/m^2"))

"""
Alias for `PrescribedValue` on `temperature` with the given name.
"""
PrescribedTemperature(name::Symbol, init=nothing) = PrescribedValue(:temperature, Input(name, init, units=u"°C"))
PrescribedTemperature(condition) = PrescribedValue(:temperature, condition)

# Hydrology BCs

"""
Alias for `PrescribedFlux` with name `infiltration` representing liquid water infiltration at the soil surface.
"""
InfiltrationFlux(init=nothing) = PrescribedFlux(:saturation_water_ice, Input(:infiltration, init, units=u"m/s"))

"""
Alias for `NoFlux` representing a zero-flux boundary condition for water flow (prognostic variable `saturation_water_ice`).
"""
ImpermeableBoundary() = NoFlux(:saturation_water_ice)

"""
Alias for `PrescribedGradient` representing a Neumann-type zero pressure gradient at the bottom of the soil
column, thereby allowing free drainage of water.
"""
FreeDrainage(::Type{NF}) where {NF} = PrescribedGradient(:pressure_head, zero(NF))

# SoilBoundaryConditions constructor

"""
Creates `ColumnBoundaryConditions` with defaults suitable for the given soil process configurations.
"""
SoilBoundaryConditions(
    ::Type{NF},
    energy::SoilEnergyBalance = SoilEnergyBalance(NF),
    hydrology::SoilHydrology = SoilHydrology(NF);
    top=SoilBC(energy = default_upperbc(energy), hydrology = default_upperbc(hydrology)),
    bottom=SoilBC(energy = default_lowerbc(energy), hydrology = default_lowerbc(hydrology))
) where {NF} = ColumnBoundaryConditions(; top, bottom)

# Default boundary conditions for various soil process configurations

default_upperbc(::SoilEnergyBalance) = GroundHeatFlux()
default_upperbc(::SoilHydrology{NF, NoFlow}) where {NF} = nothing
default_upperbc(::SoilHydrology) = InfiltrationFlux()

default_lowerbc(::SoilEnergyBalance) = GeothermalHeatFlux()
default_lowerbc(::SoilHydrology{NF, NoFlow}) where {NF} = nothing
default_lowerbc(::SoilHydrology) = ImpermeableBoundary()
