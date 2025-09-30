@kwdef struct PrescribedSoil{NF, Grid<:AbstractLandGrid{NF}} <: AbstractBoundaryConditions
    "Spatial discretization"
    grid::Grid

    "Prescribed soi layer thickness in meters [m]"
    thickness::NF = 1.0

    "Prescribed soil layer porosity (soil moisture at saturation)"
    porosity::NF = 0.5
end

variables(::PrescribedSoil) = (
    input(:soil_temperature, XY(), units=u"°C", desc="Near-surface average soil temperature [°C]"),
    input(:soil_moisture, XY(), desc="Root zone soil moisture as volumetric fraction [m^3/m^3]"),
    input(:soil_field_capacity, XY(), desc="Root zone field capacity as volumetric fraction [m^3/m^3]"),
    input(:soil_wilting_point, XY(), desc="Root zone wilting point as volumetric fraction [m^3/m^3]"),
)

# Convenience constructors for common soil upper and lower boundary conditions
VegetationBoundaryConditions(
    grid::AbstractLandGrid;
    top = PrescribedAtmosphere(; grid),
    bottom = PrescribedSoil(; grid),
) = ColumnBoundaryConditions(; top, bottom)
