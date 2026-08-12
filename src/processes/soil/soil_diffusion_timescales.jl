#####
##### Soil heat diffusion timescale
#####

# Fallback
cell_diffusion_timescale(state, grid, ::AbstractSoilThermodynamics{NF}, soil, args...) where {NF} = convert(NF, Inf)

"""
    $TYPEDSIGNATURES

Return the minimum thermal diffusion timescale over the soil column for the `SoilThermodynamics`
energy process. The per-cell timescale is the reciprocal of the diagonal decay rate of the discrete
heat-conduction operator, ``τ = C·Δzᶜ / (κᶠ_k/Δzᶠ_k + κᶠ_{k+1}/Δzᶠ_{k+1})``, where `C` is the bulk
volumetric heat capacity (J m⁻³ K⁻¹), `Δzᶜ` the cell thickness, `κᶠ` the thermal conductivity
interpolated to a face (as used by the conduction operator), and `Δzᶠ` the center-to-center spacing at
a face. This is the Gershgorin stability bound of the explicit operator: it reduces to `Δz²·C/(2κ)` on
a uniform grid but is much stricter at thin surface layers of a graded grid, where the center-to-center
spacing `Δzᶠ` — not the cell thickness — sets the limit.

The sensible heat capacity is used (not the apparent capacity, which includes the latent heat of
phase change); since the apparent capacity is never smaller, the reported timescale is conservative
near the freezing front.
"""
function cell_diffusion_timescale(
        state, grid,
        energy::SoilThermodynamics,
        soil::AbstractSoil,
        args...
    )
    strat = get_stratigraphy(soil)
    hydrology = get_hydrology(soil)
    bgc = get_biogeochemistry(soil)
    field_grid = get_field_grid(grid)
    fields = get_fields(state, energy, hydrology, strat, bgc)
    τ = KernelFunctionOperation{Center, Center, Center}(
        compute_thermal_diffusion_timescale, field_grid, fields, energy, hydrology, strat, bgc
    )
    return minimum(τ)
end

"""
    $TYPEDSIGNATURES

Kernel function returning the thermal diffusion timescale ``τ_k = C_k Δz_k^C / (κ_k^F/Δz_k^F +
κ_{k+1}^F/Δz_{k+1}^F)`` at cell `i, j, k`. The heat capacity `C` is computed at the cell centre from
the local soil composition; the face conductivities `κ^F` are interpolated from the adjacent cell
centres exactly as in the heat conduction operator, and `Δz^F` is the center-to-center spacing at the
face. Boundary faces (below cell `1` and above cell `Nz`) carry the surface/bottom boundary flux
rather than an interior conduction term and are excluded from the sum.
"""
@inline @propagate_inbounds function compute_thermal_diffusion_timescale(
        i, j, k, grid, fields,
        energy::SoilThermodynamics,
        hydrology::AbstractSoilHydrology,
        strat::AbstractStratigraphy,
        bgc::AbstractSoilBiogeochemistry
    )
    soil = soil_volume(i, j, k, grid, fields, strat, hydrology, bgc)
    C = compute_heat_capacity(energy.thermal_properties, soil)
    Δz = Δzᵃᵃᶜ(i, j, k, grid)
    Nz = grid.Nz
    # Conductance of the lower face k (between cells k-1 and k) and the upper face k+1 (between cells
    # k and k+1): κ interpolated to the face, divided by the center-to-center spacing at the face.
    # These are the same face conductivities and spacings the explicit conduction operator uses.
    κ_lower = ℑzᵃᵃᶠ(i, j, k, grid, compute_thermal_conductivity, fields, energy, hydrology, strat, bgc)
    κ_upper = ℑzᵃᵃᶠ(i, j, k + 1, grid, compute_thermal_conductivity, fields, energy, hydrology, strat, bgc)
    conductance_lower = ifelse(k > 1, κ_lower / Δzᵃᵃᶠ(i, j, k, grid), zero(Δz))
    conductance_upper = ifelse(k < Nz, κ_upper / Δzᵃᵃᶠ(i, j, k + 1, grid), zero(Δz))
    # diagonal decay rate of the discrete operator; τ is its reciprocal
    inverse_timescale = (conductance_lower + conductance_upper) / (C * Δz)
    return ifelse(inverse_timescale > zero(inverse_timescale), inv(inverse_timescale), convert(eltype(Δz), Inf))
end

"""
    $TYPEDSIGNATURES

Return the minimum diffusive timescale over the coupled soil processes as the minimum of the
individual energy (heat conduction) and hydrology (Richards) timescales.
"""
function cell_diffusion_timescale(state, grid, soil::SoilEnergyWaterCarbon, constants::PhysicalConstants)
    τ_energy = cell_diffusion_timescale(state, grid, get_energy_balance(soil), soil, constants)
    τ_water = cell_diffusion_timescale(state, grid, get_hydrology(soil), soil, constants)
    return min(τ_energy, τ_water)
end

#####
##### Soil hydrology (Richards) timescale
#####

# Fallback
cell_diffusion_timescale(state, grid, ::AbstractSoilHydrology{NF}, args...) where {NF} = NF(Inf)

"""
    $TYPEDSIGNATURES

Return the minimum hydraulic diffusion timescale ``τ = Δz² (∂θ/∂ψ) / K`` over the soil column for a
Richards-equation `SoilHydrology`, where `K` is the (variably saturated) hydraulic conductivity
(m s⁻¹) and `∂θ/∂ψ` the specific moisture capacity (m⁻¹). This is the diffusive stability limit of
the explicit Richards operator, whose effective moisture diffusivity is ``D = K ∂ψ/∂θ = K / (∂θ/∂ψ)``.

As the soil saturates, ``∂θ/∂ψ → 0`` and the timescale tends to zero: explicit Richards is stiff
near saturation, which is precisely the regime motivating implicit timestepping. The
`TimeStepWizard`'s `min_change` bound keeps the step from collapsing to zero in a single adjustment.
"""
function cell_diffusion_timescale(
        state, grid,
        hydrology::SoilHydrology{NF, RichardsEq},
        soil::AbstractSoil,
        args...
    ) where {NF}
    strat = get_stratigraphy(soil)
    bgc = get_biogeochemistry(soil)
    field_grid = get_field_grid(grid)
    fields = get_fields(state, hydrology, strat, bgc)
    τ = KernelFunctionOperation{Center, Center, Center}(
        compute_hydraulic_diffusion_timescale, field_grid, fields, hydrology, strat, bgc
    )
    return minimum(τ)
end

"""
    $TYPEDSIGNATURES

Kernel function returning the hydraulic diffusion timescale ``Δz² (∂θ/∂ψ) / K`` at cell `i, j, k`.
The hydraulic conductivity `K` is evaluated at the cell centre from the local soil composition; the
specific moisture capacity `∂θ/∂ψ` is the analytic derivative of the soil-water retention curve
(from FreezeCurves) evaluated at the matric potential `ψₘ`, reconstructed from the stored total
`pressure_head` by removing its hydrostatic and elevation components. Cells with zero hydraulic
conductivity (fully dry or frozen) impose no restriction and return `Inf`.
"""
@inline @propagate_inbounds function compute_hydraulic_diffusion_timescale(
        i, j, k, grid, fields,
        hydrology::SoilHydrology{NF, RichardsEq},
        strat::AbstractStratigraphy,
        bgc::AbstractSoilBiogeochemistry
    ) where {NF}
    soil = soil_volume(i, j, k, grid, fields, strat, hydrology, bgc)
    # hydraulic conductivity at the cell centre
    K = hydraulic_conductivity(get_hydraulic_properties(hydrology), soil)
    # reconstruct the matric potential ψₘ = ψ - ψ_hydrostatic - ψ_elevation from the total head
    ψ = fields.pressure_head[i, j, k]
    z = znode(i, j, k, grid, Center(), Center(), Center())
    z_ref = znode(i, j, grid.Nz + 1, grid, Center(), Center(), Face())
    ψz = z - z_ref
    z₀ = fields.water_table[i, j, 1]
    ψh = max(zero(z), z₀ - z)
    ψm = ψ - ψh - ψz
    # specific moisture capacity ∂θ/∂ψ (m⁻¹) from the analytic SWRC derivative
    swrc = get_swrc(hydrology)
    por = porosity(i, j, k, grid, fields, strat, bgc)
    ∂θ∂ψ = convert(NF, swrc(FreezeCurves.derivative, ψm; θsat = por))
    # moisture diffusivity D = K / (∂θ/∂ψ) ⟹ τ = Δz² / D = Δz² (∂θ/∂ψ) / K
    Δz = Δzᵃᵃᶜ(i, j, k, grid)
    return ifelse(K > zero(K), Δz^2 * abs(∂θ∂ψ) / K, convert(NF, Inf))
end
