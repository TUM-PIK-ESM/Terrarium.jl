# Simple single-layer snow scheme

> Status: **planned**. Initial design

Date of initial draft: 2026-07-01

Base revision: 462c6422b22829b88e491457894264419704d691

## Originating prompt

> Please review the existing Terrarium processes and prepare a plan for the implementation of a
> simple, single-layer snow scheme. The scheme should be loosely based on the Utah energy balance
> model (see attachment) and should consist of a single layer with homogeneous density and thermal
> properties. The surface temperature should be handled by the surface energy balance whose
> interface will need to be modified to optionally except an AbstractSnow type. Write your initial
> plan into the `scratch/design` folder. Include this prompt text at the top. Feel free to ask any
> clarifying questions that arise.

## Revision log

Follow-up prompt (revision 2), verbatim:

> Make sure to append my prompts to the design document to keep track of this process. Regarding the
> open question, I agree that this is a problem and would update the plan to make the depth-integrated
> energy content the prognostic variable while exposing the volumetric energy content to the closure.
> This should be done via a method and without allocating an auxiliary variable.
>
> A couple of additional problems need to be addressed, however, before starting the implementation:
> 1. Currently, `ground_heat_flux` is directly wired to the soil temperature state variable and is
>    also used in the SEB. This becomes problematic here because our snow surface heat flux would be
>    distinct from the conductive heat flux at the ground-snow interface. Investigate this ambiguity
>    and suggest a workaround or design change to the current BC handling.
> 2. The existing soil energy and hydrology code is set up and named to be for soil. However, much of
>    the physics is transferrable as noted in the design plan. We should consider an initial
>    refactoring where the existing "soil" module is renamed to "ground" with soil-specific
>    parameterizations where appropriate. This should be done as a separate, preliminary work item
>    (and associated pull request) to avoid too many changes in one go.
> Please provide your feedback regarding these points and update the plan accordingly.

Resolutions adopted in this revision: (i) the prognostic snow energy is the **depth-integrated**
(column) energy content (J/m²); the volumetric energy seen by the closure is computed by a **method**
with no auxiliary field (see "Physical formulation"). (ii) The `ground_heat_flux` dual-role ambiguity
is addressed by splitting it into `surface_heat_flux` and `ground_heat_flux` (see "Boundary-condition
disambiguation"). (iii) The soil→ground refactor is adopted as a scoped, abstraction-first preliminary
work item (see "Preliminary refactoring work items").

## Resolved design decisions

The following forks were resolved with the model author before drafting:

1. **Thermal state representation — depth-integrated energy prognostic, volumetric closure
   (revised).** The prognostic snow energy is the **depth-integrated (column) energy content** `E`
   (J/m²), which is the conserved quantity and keeps the energy tendency free of the moving-boundary
   `1/d_s` factor that the original volumetric prognostic introduced. The **volumetric** energy
   `U_v = E/d_s` is exposed to the existing soil enthalpy closure (`FreeWater`,
   [`soil_energy_closures.jl`](../../src/processes/soil/energy/soil_energy_closures.jl)) via a
   **method**, with no auxiliary field allocated for `U_v`. Snow layer depth `d_s` is a diagnostic
   derived from SWE and a homogeneous bulk density. This realizes the original intent (reuse the
   volumetric phase-change closure) while fixing the conservation/regularization concern raised in
   review.
2. **Scope — process first, then couple.** Deliver an `AbstractSnow` process plus a minimal
   standalone `SnowModel` that can be unit- and Enzyme-tested in isolation, *then* wire it into
   `LandModel` as an optional component (default `NoSnow`).
3. **Meltwater — continuous Darcy-type outflow.** Liquid fraction is diagnosed from energy content;
   meltwater drains to the soil surface as a continuous outflow flux (UEB eqns 23–24), preserving
   continuous-time, differentiable dynamics per `AGENTS.md`.
4. **Surface transition — fractional snow cover.** A diagnosed sub-grid snow-covered area fraction
   `f_snow ∈ [0,1]` blends snow-surface and bare-soil/soil contributions in the surface energy
   balance and in the soil boundary fluxes.

## Background: the Utah Energy Balance (UEB) model

Tarboton, Chowdhury & Jackson (1994) describe a lumped, depth-averaged, single-layer snowpack with
two prognostic state variables: water equivalence `W` (m) and energy content `U` (kJ/m²) relative to
ice at 0 °C. The governing ODEs (their eqns 1–2) are

```
dU/dt = Q_sn + Q_li + Q_p + Q_g − Q_le + Q_h + Q_e − Q_m
dW/dt = P_r + P_s − M_r − E
```

with net shortwave `Q_sn`, incoming/outgoing longwave `Q_li`/`Q_le`, advected precipitation heat
`Q_p`, ground heat flux `Q_g`, sensible/latent fluxes `Q_h`/`Q_e`, advected meltwater heat `Q_m`,
rain/snow inputs `P_r`/`P_s`, meltwater outflow `M_r`, and sublimation `E`. Depth-averaged snow
temperature follows from `(U, W)` by phase partitioning (eqns 6–8); the **snow surface temperature**
`T_s` is solved separately from an equilibrium surface energy balance (eqns 20–22) because snow is a
strong insulator and `T_s ≠ T_avg`. Meltwater outflow uses a Darcy relation on the liquid fraction
in excess of capillary retention (eqns 23–24).

### Adaptation to Terrarium

- UEB's predictor–corrector time integration is **not** carried over. Terrarium supplies the
  integrator; the snow scheme contributes only well-formed ODE tendencies via `compute_tendencies!`.
- UEB's effective soil interaction depth `D_e` (which lumps a 0.4 m soil slab into `U`) is **dropped**.
  Terrarium already resolves an explicit multi-layer soil column, so the snow energy is snow-only and
  snow↔soil exchange is an explicit conductive flux at the snow base. This is a deliberate departure.
- UEB's equilibrium surface-temperature parameterization is realized by Terrarium's existing
  **skin-temperature solve** in the surface energy balance, which already solves the same equilibrium
  balance via a Newton iteration ([`skin_temperature.jl`](../../src/processes/surface_energy/skin_temperature.jl),
  [`surface_energy_balance.jl`](../../src/processes/surface_energy/surface_energy_balance.jl)). The
  snow scheme supplies the *sub-surface* side of that balance (temperature, conductance, layer
  thickness) instead of the soil top layer.
- Precipitation phase partitioning is already available upstream: `RainSnow`
  ([`prescribed_atmosphere.jl`](../../src/processes/atmosphere/prescribed_atmosphere.jl)) provides
  `rainfall` and `snowfall` inputs. The snow scheme consumes these directly rather than re-deriving
  the temperature-threshold partition (UEB eqn 12).

## How the relevant existing machinery works

These are the integration points the scheme attaches to (verified against current source):

- **Skin temperature / SEB.** `SurfaceEnergyBalance{NF, SkinTemperature, TurbulentFluxes,
  RadiativeFluxes, Albedo}`. The `ImplicitSkinTemperature` solve closes
  `R_net(T_s) + H_s(T_s) + H_l(T_s) − G(T_s, T_g) = 0` (positive-upward convention). The conduction
  term uses `compute_skin_temperature(skinT, Tg, G, Δz) = Tg − G·Δz/(2κₛ)`, where `Tg` is the
  `ground_temperature` input and `Δz` is the **topmost soil cell** thickness
  ([`skin_temperature.jl:159-172`](../../src/processes/surface_energy/skin_temperature.jl#L159-L172)).
  An in-code TODO at line 85 already anticipates pulling `κₛ`/thickness from snow/canopy.
- **Ground temperature source.** `SoilEnergyBalance` registers `ground_temperature` as a `view` of the
  uppermost soil layer, with the TODO "Revisit this if/when we extend the vertical layers to include
  snow and canopy" ([`soil_energy.jl:51-57`](../../src/processes/soil/energy/soil_energy.jl#L51-L57)).
  The snow scheme is the resolution of that TODO.
- **Soil top boundary fluxes.** `initialize(model::LandModel)` wires the `ground_heat_flux` auxiliary
  as a `FluxBoundaryCondition` on soil `internal_energy`, and `infiltration` (negated) as a flux on
  `saturation_water_ice` ([`land_model.jl:55-65`](../../src/models/coupled/land_model.jl#L55-L65),
  [`soil_model_bcs.jl`](../../src/models/soil/soil_model_bcs.jl)). The snow scheme reuses these two
  fields unchanged — it only changes *how their values are computed* when snow is present.
- **Coupling order.** `LandModel.compute_auxiliary!` runs atmos → soil → vegetation →
  surface_hydrology → SEB; `compute_tendencies!` runs surface_hydrology → soil → vegetation
  ([`land_model.jl:80-96`](../../src/models/coupled/land_model.jl#L80-L96)).
- **Enthalpy closure to reuse.** `FreeWater` provides `liquid_water_fraction(fc, U, Lθ, sat)` and
  `energy_to_temperature(fc, U, Lθ, C)` operating on volumetric energy `U`, latent-heat content `Lθ`,
  and heat capacity `C` ([`soil_energy_closures.jl:135-163`](../../src/processes/soil/energy/soil_energy_closures.jl#L135-L163)).
  The snow closure is the same enthalpy map with snow-specific `Lθ` and `C`.
- **Constants.** `PhysicalConstants` exposes `latent_heat_fusion`, `latent_heat_sublimation`,
  `specific_heat_capacity_ice`, `specific_heat_capacity_liquid_water`, `density_water`, `density_ice`,
  and the Stefan–Boltzmann constant ([`physical_constants.jl`](../../src/processes/physical_constants.jl)).
- **Process template.** `BareGroundEvaporation`
  ([`bare_ground_evaporation.jl`](../../src/processes/surface_hydrology/evapotranspiration/bare_ground_evaporation.jl))
  and the degree-day example ([`examples/extending/simple_snow_ddm.jl`](../../examples/extending/simple_snow_ddm.jl))
  are the structural patterns for a new `XY()`-fielded process and its kernels.

## Physical formulation

All snow fields are `XY()` (single layer, one value per column, indexed `[i, j, 1]`). Bulk snow
density `ρ_s` is a homogeneous, constant parameter (the "simple" assumption); homogeneous thermal
properties follow from it.

### State variables (prognostic)

- `snow_energy` `E` (J/m²) — depth-integrated (column) energy content relative to ice at 0 °C. This
  is the conserved quantity; surface and basal fluxes add to it directly.
- `snow_water_equivalent` `W` (m) — total water substance (ice + retained liquid), per UEB `W`.

### Diagnostics (auxiliary / closure)

- `snow_depth` `d_s = W·ρ_w/ρ_s` (m).
- `snow_temperature` `T_snow` (°C) and `snow_liquid_fraction` `θ_liq ∈ [0,1]` — from the enthalpy
  closure.
- **Volumetric energy `U_v = E/d_s` (J/m³) is *not* an auxiliary field.** It is returned by a method
  (e.g. `volumetric_snow_energy(i, j, grid, fields, snow)`) evaluated inside the closure kernel, using
  `safediv` to remain finite as `W → 0` (the indeterminate limit is masked downstream by `f_snow → 0`).
- `snow_cover_fraction` `f_snow ∈ [0,1]` — smooth function of `W`, e.g. `f_snow = W/(W + W₀)` or
  `tanh(W/W_ref)`; differentiable, → 0 as `W → 0`.
- `snow_thermal_conductivity` `κ_snow` (W/m/K) — from `ρ_s` (e.g. Sturm/Yen power law `κ = a·ρ_s^b`),
  constant for fixed `ρ_s`.
- `ground_heat_flux` (reused), `infiltration` (reused) — see soil coupling below.

### Enthalpy closure (reusing `FreeWater`)

Treat the bulk snow as a water-substance medium of volumetric mass `ρ_s` (kg/m³), fraction `θ_liq`
liquid and `1−θ_liq` ice. Then

```
U_v  = E / d_s                        # volumetric energy via method (not stored)
Lθ_v = ρ_s · L_f                      # volumetric latent heat for complete melt (J/m³)
C_v  = ρ_s · (θ_liq·c_w + (1−θ_liq)·c_i)   # volumetric heat capacity (J/m³/K)
θ_liq   = liquid_water_fraction(FreeWater(), U_v, Lθ_v, 1)
T_snow  = energy_to_temperature(FreeWater(), U_v, Lθ_v, C_v)
```

This is the identical enthalpy map used for soil with `sat·por → 1` (snow is "fully saturated" with
water substance in the lumped sense). It is implemented as a `closure!`/`invclosure!` pair so that
`T_snow`/`θ_liq` are recovered from `U_v` exactly as soil temperature is recovered from soil energy.
The only division by `d_s` lives here, inside the closure (guarded with `safediv`), rather than in the
prognostic tendency. `invclosure!` (initialization from a prescribed `T_snow`) sets
`E = (T·C_v − Lθ_v·(1−θ_liq))·d_s`, so `W` must be initialized before `snow_energy`.

### Mass balance (continuous)

```
dW/dt = f_snow·0 + P_s + R_on_snow − M_r − E_subl
```

- `P_s` snowfall (mass source; can initiate a pack from `W = 0`).
- `R_on_snow = f_snow·rainfall` — rainfall intercepted by the snow-covered fraction (adds to `W`,
  advects sensible/latent heat into `U`). Rainfall on the bare fraction bypasses the pack (→ soil).
- `M_r` meltwater outflow (Darcy, below).
- `E_subl` sublimation/evaporation from the snow surface, supplied by the latent-heat flux of the SEB
  over the snow-covered fraction (couples to the existing latent-flux machinery).

### Meltwater outflow (UEB eqns 23–24, continuous)

```
S* = clamp((θ_liq/(1−θ_liq) − L_c) / (ρ_w/ρ_s − ρ_w/ρ_i − L_c), 0, 1)
M_r = K_sat · S*^3
```

with capillary retention `L_c` and saturated snow conductivity `K_sat` as parameters. `S*` is the
relative saturation in excess of capillary retention; the cubic Darcy form gives continuous outflow
that persists while liquid is present (even under a negative energy balance), as in UEB. Output is
added to soil infiltration (below). `clamp` is expressed with `min`/`max` (kernel-safe).

### Energy balance (continuous, depth-integrated)

With `E` (J/m²) prognostic, the column energy tendency is a direct flux sum — no `1/d_s` factor and no
moving-boundary term:

```
dE/dt = Q_top − Q_base + Q_precip − Q_melt
```

- `Q_top` — conductive flux delivered to the snow surface, equal to the SEB surface heat flux over the
  snow-covered fraction (sign set so energy leaving the snow top reduces `E`); see "Boundary-condition
  disambiguation".
- `Q_base` — snow→soil basal conduction (below); equals the soil-top `ground_heat_flux`.
- `Q_precip` — advected heat of `P_s` (at `min(T_air,0)`) and `R_on_snow` (latent + sensible relative
  to 0 °C ice reference), UEB eqn 13.
- `Q_melt = ρ_w·L_f·M_r` — heat advected out with meltwater at 0 °C (UEB eqn 25).

The depth dependence now appears only inside the enthalpy closure (`U_v = E/d_s`), where it is guarded
and masked by `f_snow → 0` in the thin-snow limit.

## Surface-energy-balance interface change

The SEB must "optionally accept an `AbstractSnow`". Plan:

1. **Thread an optional `snow` argument** through the SEB entry points, mirroring the existing
   `hydrology::Optional{AbstractSurfaceHydrology} = nothing` pattern:
   - `compute_auxiliary!(state, grid, seb, constants, atmos, hydrology, snow=nothing)`
   - `solve_surface_energy_balance!`, `compute_surface_energy_fluxes!`, and the fused kernels gain a
     trailing `snow::Optional{AbstractSnow} = nothing`.
2. **Snow-aware sub-surface conduction.** Replace the direct reads in
   `compute_skin_temperature(i, j, grid, fields, skinT)` with accessor functions dispatched on the
   presence of snow, returning the fraction-weighted conduction target and conductance:
   ```
   Tg_eff = f_snow·T_snow + (1−f_snow)·T_soil_top
   κ_eff  = f_snow·κ_snow + (1−f_snow)·κₛ
   Δz_eff = f_snow·d_s     + (1−f_snow)·Δz_soil_top
   ```
   The default method (`snow === nothing`) returns today's soil-only values, so non-snow behavior is
   byte-for-byte unchanged. This also resolves the `κₛ`/thickness TODO at `skin_temperature.jl:85`.
   The SEB writes its closure flux into `surface_heat_flux` (renamed from `ground_heat_flux` in
   preliminary PR-B); the snow energy balance reads it as `Q_top`.
3. **Snow-aware albedo/emissivity.** Add a `SnowAlbedo <: AbstractAlbedo` (or make the SEB blend),
   returning `f_snow`-weighted snow vs underlying albedo/emissivity. A first pass uses constant snow
   albedo/emissivity (UEB Table 2: `A_bg = 0.25`, snow `ε_s = 0.99`); snow-age decay (UEB / BATS) is
   deferred.
4. **Single composite skin temperature.** Fractional cover is realized as a *composite blended
   surface* (one skin temperature per column, fluxes weighted by `f_snow`) rather than separate
   energy-balance tiles. This is the pragmatic single-layer realization; independent snow/bare tiles
   are noted as future work.

## Boundary-condition disambiguation (surface vs ground heat flux)

### Investigation of the current ambiguity

`ground_heat_flux` is a single `XY` auxiliary `Field` that currently serves **two distinct roles**:

1. *SEB-internal.* The skin-temperature solve writes it as `G = R_net + H_s + H_l`
   ([`skin_temperature.jl:98-101,179-191`](../../src/processes/surface_energy/skin_temperature.jl#L98-L101))
   and reads it back to derive `Ts = Tg − G·Δz/(2κ)` ([`skin_temperature.jl:159-172`](../../src/processes/surface_energy/skin_temperature.jl#L159-L172)).
   This is the *conductive flux from the skin into the medium immediately below it*.
2. *Soil coupling.* `initialize(model::LandModel)` wraps the same field as a `FluxBoundaryCondition`
   on soil `internal_energy` (top), via `SoilHeatFlux(...)`
   ([`land_model.jl:55-65`](../../src/models/coupled/land_model.jl#L55-L65),
   [`soil_model_bcs.jl:6`](../../src/models/soil/soil_model_bcs.jl#L6)); the soil energy tendency
   consumes it through `compute_z_bcs!` ([`boundary_conditions.jl:36-38`](../../src/boundary_conditions.jl#L36-L38)).

Without snow these two roles describe the **same interface** (the skin sits directly on the soil top),
so one field is consistent. **With snow they are two different fluxes:** the SEB closes the balance
against the *snow surface*, so `G` is the flux into the *top* of the snowpack, whereas the soil sees
the conduction across the *snow base*, `Q_base(T_snow, T_soil_top)`. The two differ by snowpack energy
storage `dE/dt` plus advected meltwater heat. Feeding the SEB's `G` into the soil BC would inject the
surface flux directly into the soil, bypass the snowpack, and break energy conservation.

### Recommended design change: split the field by semantics

Decouple the two roles into two named `XY` fluxes with unambiguous meaning:

- **`surface_heat_flux`** — the SEB skin-temperature closure flux: conduction from the skin into the
  medium directly beneath it (snow top, or ground top when no snow). Owned and written by the SEB; it
  is the `G` of the skin-temperature equation and the top boundary of the uppermost sub-surface medium.
- **`ground_heat_flux`** — the conductive flux into the top of the soil/ground column; feeds the soil
  energy BC (wiring otherwise unchanged). This name now aligns with the soil→ground refactor below
  ("flux into the ground").

Coupling rule (dispatched on the uppermost medium):

- **No snow** (`NoSnow`): `ground_heat_flux ≡ surface_heat_flux`.
- **Snow present**: `ground_heat_flux = f_snow·Q_base + (1−f_snow)·surface_heat_flux`, where `Q_base`
  is the snow→soil conduction (Fourier between `T_snow` and `T_soil_top` across series snow+soil
  conductances). `surface_heat_flux` is the top boundary of the **snow** energy balance (`Q_top`).

`Q_base` depends only on `T_snow` and `T_soil_top` (both available after the soil and snow auxiliary
passes, before the SEB), so there is no circular dependency with the skin-temperature solve.

### Two implementation strategies for selecting the soil-BC source

- **(A) Static wiring by model type (recommended).** The snow component type is known at
  `initialize(model::LandModel)` time. Wire the soil-top energy BC to `surface_heat_flux` when the
  component is `NoSnow`, and to `ground_heat_flux` (written by the snow scheme) otherwise. No
  dual-writer field, no runtime branch, and the no-snow path is byte-for-byte identical to today (only
  the field is renamed).
- **(B) Single writer via an explicit coupling step.** Always wire the soil BC to `ground_heat_flux`
  and add `compute_ground_heat_flux_coupling!(state, grid, snow, seb, ...)` after both the snow and SEB
  auxiliary passes, dispatching on the `snow` type to copy `surface_heat_flux` (NoSnow) or the blended
  `Q_base` (snow). More uniform, at the cost of one extra pass.

Strategy (A) is preferred. Either way, the `surface_heat_flux` rename is a small, **snow-independent**
change that is best landed as a preliminary PR (see below), as it touches the SEB, the
`SoilHeatFlux`/`PrescribedSurfaceTemperature` aliases, and the `LandModel`/`SurfaceEnergyModel` wiring.

### Water-flux analogue (`infiltration`)

`infiltration` (feeds soil `saturation_water_ice` top BC) has the same dual-role structure and is
handled identically: with snow, `infiltration = f_snow·M_r + (1−f_snow)·(rainfall + throughfall)` plus
the existing surface-hydrology contributions, so meltwater and rain-through reach the soil as a
continuous flux. The energy flux is the headline ambiguity; the water flux is resolved by the same
split-and-blend pattern and should be addressed in the same preliminary PR for consistency.

## Preliminary refactoring work items

Two structural changes should land as their own PRs **before** the snow scheme, to keep each change
reviewable and to avoid bundling unrelated churn (cf. `AGENTS.md` "scope creep in PRs").

### PR-A: soil/substrate abstraction (soil → ground, scoped)

> Detailed design and revised recommendation: [`pr_a_ground_abstraction.md`](pr_a_ground_abstraction.md).
> That document supersedes the "additive vs full rename" recommendation below — the author has opted for
> the directory rename `soil/` → `soil/` with an `AbstractGround…` abstraction layer.

I agree with the direction: much of the soil **energy** machinery is medium-agnostic and the snow
scheme reuses it. However, a blanket directory/type rename of everything under `processes/soil` would
be both too wide and partly incorrect, because several subsystems are intrinsically soil-specific.
Recommended scoping:

- **Generalize (move under a `ground`/`substrate` abstraction):**
  - the 1D heat-conduction operator `ExplicitTwoPhaseHeatConduction` (Fourier divergence);
  - the `FreeWater` enthalpy closure and the scalar maps `energy_to_temperature` /
    `liquid_water_fraction` (water/ice phase change — reused verbatim by snow);
  - an `AbstractGround` supertype and the shared conductive-coupling interface (`ground_temperature`
    accessor and the snow-aware conduction accessors the SEB dispatches on).
- **Keep soil-specific (as the soil *parameterization* of a ground medium):** stratigraphy
  (`SoilTexture`, `SoilHorizon`, porosity, `soil_volume`), Richards hydrology and hydraulic closures,
  the thermal *constituent* properties (mineral/organic/quartz conductivities, texture-dependent λ —
  cf. [`texture_dependent_thermal_conductivity.md`](texture_dependent_thermal_conductivity.md)), and
  soil biogeochemistry.

Two ways to realize this, in increasing churn:

- **Additive (recommended first step):** introduce `AbstractGround{NF}` with `AbstractSoil <:
  AbstractGround`, and lift only the generic energy operator + closure into a shared location, leaving
  public types (`SoilModel`, `SoilEnergyWaterCarbon`, `SoilEnergyBalance`) and the `soil/` tree in
  place. This unlocks snow reuse with minimal disruption to imports (ExplicitImports tests), docs,
  examples, and the public API.
- **Full rename:** rename the `soil` module/tree to `ground`. Higher cost (wide import/doc/export/API
  churn) for largely cosmetic benefit beyond the additive step; defer unless there is a concrete need.

Recommendation: do the **additive** abstraction extraction now; treat a full cosmetic rename as
optional later work. This still satisfies the goal (snow reuses the ground energy/closure machinery)
without a sweeping rename.

### PR-B: surface/ground heat-flux split

> Detailed design: [`pr_b_surface_ground_heat_flux_split.md`](pr_b_surface_ground_heat_flux_split.md).
> Note the refinement there: the `infiltration` blend is deferred to the snow PR (no dual-role exists
> for the water flux until snow is present), so PR-B is scoped to the energy-flux split only.

Implement the `surface_heat_flux` / `ground_heat_flux` disambiguation (and the `infiltration`
analogue) described above, using static wiring strategy (A). Snow-independent; verifiable by a
no-snow regression (current results unchanged).

Suggested order: **PR-A then PR-B** (PR-B's `ground_heat_flux` naming and the shared conductive-coupling
interface read more naturally on top of the ground abstraction), though the two are largely
independent and could be reordered. The snow PR depends on both.

## Module structure (new and modified files)

```
src/processes/snow/
├── abstract_types.jl     # AbstractSnow{NF}, NoSnow, interface accessors
│                         #   (snow_temperature, snow_depth, snow_cover_fraction,
│                         #    snow_thermal_conductivity, ...)
├── snow_properties.jl    # ρ_s, depth d_s(W), κ_snow(ρ_s), f_snow(W), volumetric_snow_energy method,
│                         #   albedo helpers
├── snow_energy.jl        # enthalpy closure (reuse FreeWater via U_v = E/d_s), column energy tendency
├── snow_mass.jl          # SWE mass balance, Darcy meltwater outflow, sublimation coupling
└── snow.jl               # SingleLayerSnow concrete process + variables/initialize!/
                          #   compute_auxiliary!/compute_tendencies!/closure!

src/models/snow/
└── snow_model.jl         # minimal standalone SnowModel for unit + AD testing

modified (snow PR — builds on preliminary PR-A and PR-B):
src/processes/surface_energy/skin_temperature.jl     # snow-aware conduction accessors
src/processes/surface_energy/surface_energy_balance.jl # optional snow arg threading
src/processes/surface_energy/albedo.jl               # SnowAlbedo (or blend)
src/models/coupled/land_model.jl                     # @component snow (default NoSnow), ordering,
                                                     #   pass snow to SEB, blend ground_heat_flux +
                                                     #   infiltration
src/Terrarium.jl                                     # includes + exports

(The surface_heat_flux/ground_heat_flux split and the AbstractGround abstraction land earlier in the
preliminary PRs; see "Preliminary refactoring work items".)
```

Naming: concrete type `SingleLayerSnow` (UEB-inspired); `NoSnow` is the inert default analogous to
`NoCanopyInterception`. `AbstractSnow{NF} <: AbstractCoupledProcesses{NF}`.

## Implementation phases

**Preliminary (separate PRs, before snow):**

- **PR-A — soil/substrate abstraction.** Additive `AbstractGround` supertype + extraction of the
  generic heat-conduction operator and `FreeWater` closure into a shared location; soil-specific
  parameterizations stay put. No behavior change.
- **PR-B — surface/ground heat-flux split.** Rename the SEB closure flux to `surface_heat_flux`,
  reserve `ground_heat_flux` for the soil-top BC, static wiring strategy (A); same for `infiltration`.
  Verified by a no-snow regression.

**Snow scheme:**

1. **Process skeleton + properties.** `AbstractSnow`/`NoSnow`/`SingleLayerSnow`, `variables`,
   `snow_properties.jl` (depth, density, conductivity, cover fraction, `volumetric_snow_energy` method).
   No coupling yet.
2. **Energy closure.** `closure!`/`invclosure!` reusing `FreeWater` via `U_v = E/d_s` (guarded);
   verify `T_snow`/`θ_liq` recovery and round-trip `invclosure!`→`closure!`.
3. **Tendencies.** Mass balance + Darcy outflow + depth-integrated energy tendency `dE/dt`, driven by
   prescribed top/base fluxes (standalone `SnowModel`).
4. **SEB interface change.** Thread optional `snow`; snow-aware conduction accessors; `SnowAlbedo`;
   confirm non-snow path unchanged.
5. **LandModel integration.** Add `@component snow = NoSnow(...)`; insert snow auxiliary (closure +
   `Q_base`) before SEB and snow tendencies after surface hydrology; blend `ground_heat_flux` and
   `infiltration`.
6. **Docs.** A `docs/` page following the AGENTS.md process-page template (Overview / Implementations /
   Methods / Kernel functions), with a "scheme under development" warning.

## Testing

- **Unit:** `f_snow`, `d_s`, `κ_snow` monotonicity and limits; closure recovers `T_snow`/`θ_liq` and
  round-trips with `invclosure!`; `M_r` zero below capillary retention and increasing in `θ_liq`.
- **Conservation:** integrate the standalone `SnowModel` with zero surface/base flux and check column
  energy `E = U_v·d_s` and water `W` conservation under accumulation, melt, and complete ablation;
  verify no energy/water created as `W → 0`.
- **Coupling regression:** `LandModel` with `NoSnow` reproduces current results bit-for-bit (the
  default conduction accessors return soil-only values).
- **Differentiability:** Enzyme adjoints through `compute_tendencies!` and the closure, following
  `test/differentiability`; the `1/d_s` regularization must keep gradients finite as `W → 0`.
- **GPU/type stability:** `SingleLayerSnow` kernels are `Float32`-clean and allocation-free; no literal
  `0.0`/`1.0`; `ifelse`/`min`/`max` instead of branches.

## Open questions / deferred work

1. **(Resolved)** ~~`1/d_s` regularization near `W → 0`.~~ Resolved in revision 2: the prognostic is the
   depth-integrated energy `E` (J/m²), so the tendency carries no `1/d_s` factor. The single remaining
   division is inside the closure (`U_v = E/d_s`), guarded with `safediv` and masked by `f_snow → 0`.
2. **Composite surface vs separate tiles.** The single blended skin temperature is an approximation of
   fractional cover. Separate snow/bare energy-balance tiles (two skin temperatures) are more rigorous
   but heavier; deferred.
3. **Snow albedo dynamics.** Constant snow albedo first; age-/temperature-dependent decay (UEB/BATS,
   Dickinson et al. 1993) and the shallow-snow albedo interpolation (UEB `r = (1−z/h)e^{−z/2h}`)
   deferred.
4. **Snow density evolution.** `ρ_s` is constant ("simple" scheme). Compaction/densification (and the
   resulting time-varying `κ_snow`, `d_s`) is a clear future extension and would relax the homogeneity
   assumption.
5. **Aerodynamic roughness over snow.** UEB reduces turbulent transfer with a snow roughness; the first
   pass keeps the existing aerodynamic resistance. Worth revisiting if turbulent fluxes over snow are
   biased.
6. **Stability corrections.** Like UEB, neutral turbulent transfer is assumed; Richardson/Monin–Obukhov
   corrections are out of scope.

## References

- Tarboton, D. G., Chowdhury, T. G., & Jackson, T. H. (1994). *A Spatially Distributed Energy Balance
  Snowmelt Model.* Utah Water Research Laboratory Working Paper WP-94-HWR-DGT/003. (Attached.)
- Male, D. H., & Gray, D. M. (1981). Snowcover ablation and runoff. In *Handbook of Snow*. (UEB melt
  outflow, eqn 9.45.)
- Dickinson, R. E., Henderson-Sellers, A., & Kennedy, P. J. (1993). BATS v1e. NCAR/TN-387+STR. (Snow
  albedo; deferred.)
