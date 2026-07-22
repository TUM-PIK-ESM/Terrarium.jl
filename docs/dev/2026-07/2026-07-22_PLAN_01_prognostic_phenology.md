# Prognostic phenology: prescribed LAI + GDD-based PALADYN phenology factor

> Status: **completed** (implementation + verification). Part A (`PrescribedPhenology` derives the
> phenology factor) and Part B (prognostic GDD reformulation of `PALADYNPhenology`) are both
> implemented and verified; documentation pages and Enzyme adjoint tests remain as follow-up. Part A
> required finishing the previously-unexercised `kernel(...)` lazy-auxiliary API. Covers (a) the
> `PrescribedPhenology` scheme's derivation of the phenology factor from prescribed LAI, and (b) a
> continuous-time, prognostic reformulation of the PALADYN GDD phenology that removes the need to
> store previous timesteps' inputs.

Date of initial draft: 2026-07-22

Base revision: b53747d7fb9bd12c213d3a144a3fc064cf519dfd

## Revision log

> 2026-07-22 — Initial draft.
>
> 2026-07-22 — Part A implemented and verified. `PrescribedPhenology` derives `phenology_factor`
> via the `kernel(...)` lazy-auxiliary API. Getting this to work end-to-end required finishing that
> API, which `PrescribedPhenology` is the first consumer of:
> - `src/abstract_variables.jl`: added an `auxiliary(::Variable, ::KernelFunction, ::Nothing)`
>   dispatch (the convenience constructor previously only accepted `ctor::Function` or `::Nothing`),
>   and relaxed the `AuxiliaryVariable` `ctor` field type bound from `Union{Nothing, Function}` to
>   unconstrained (a `KernelFunction` is a callable struct, not a `Function`).
> - `src/utils/kernel_utils.jl`: the `KernelFunction` call operators now (i) unwrap the Terrarium
>   grid via `get_field_grid` before constructing the Oceananigans `KernelFunctionOperation`,
>   (ii) pass location *types* (`map(typeof, location(...))`) as required by `KernelFunctionOperation`
>   indexing, and (iii) invert the `clock` kwarg so `clock = false` (the default) yields the
>   autonomous variant that does not pass the clock to the kernel function.
> - `src/processes/vegetation/phenology/prescribed_phenology.jl`: `variables(phenol::...)` binds the
>   process instance and passes it through `kernel(phenology_factor, phenol)` so the kernel function
>   receives the parameters. Verified: for prescribed `LAI = 2.0`, `max_leaf_area_index = 4.0`, the
>   computed `phenology_factor` is `0.5`. Full vegetation, state-variable, and abstract-variable test
>   suites pass (the one `parameter_handling` failure is pre-existing and unrelated).
>
> 2026-07-22 — Part B implemented. `PALADYNPhenology` rewritten as a prognostic GDD scheme:
> - New parameters `T_gdd_base`, `gdd_crit`, `T_senescence_range`, `gdd_relaxation_time`,
>   `f_deciduous` (placeholder defaults; PALADYN Table 5 values flagged TODO).
> - New prognostic variable `growing_degree_days` (K⋅day) with tendency
>   `d(gdd)/dt = max(T_air − T_gdd_base, 0)/spd − 𝟙[T_air<T_gdd_base]·gdd/(gdd_relaxation_time·spd)`,
>   where `spd = seconds_per_day` (accumulation in degree-days while the timestepper integrates in
>   seconds). Diagnostic `ϕ = f_deciduous·min(gdd/gdd_crit, senescence_ramp) + (1−f_deciduous)`,
>   `LAI = ϕ·LAI_b`.
> - Air temperature is obtained through the `air_temperature(i, j, grid, fields, atmos)` accessor, so
>   `atmosphere` is now threaded into the phenology `compute_auxiliary!`/`compute_tendencies!`. This
>   required adding `atmosphere` to `VegetationCarbon`'s `compute_tendencies!` and to the
>   `VegetationModel` / `LandModel` tendency couplings (`compute_tendencies!(..., constants,
>   atmosphere)`); `PrescribedPhenology`'s no-op methods absorb the extra argument via `args...`.
> - Added the `sitchEvaluationEcosystemDynamics2003` bibliography entry (cited by the new docstrings).
> - Rewrote `test/vegetation/phenology_tests.jl`: 26 tests across green-up/mature/senescence/dormant
>   regimes, the evergreen limit, the GDD tendency (warm accumulation, cold relaxation, zero at base),
>   and in-model `PALADYNPhenology`/`PrescribedPhenology`. Full vegetation suite green. Regime values
>   verified against hand calculations (e.g. gdd=100/gdd_crit=200 ⇒ ϕ=0.5; T_air=2 °C ⇒ ϕ=0.4).

## Problem description

The autotrophic respiration scheme `PALADYNAutotrophicRespiration` requires a `phenology_factor`
input `ϕ ∈ [0, 1]` (leaf-out fraction), consumed in the root-maintenance-respiration term
(`src/processes/vegetation/respiration/autotrophic_respiration.jl`, `compute_Rm`). This factor must
be produced by whichever `AbstractPhenology` implementation is coupled into the vegetation model.

Two phenology implementations exist:

- `PrescribedPhenology` — treats `leaf_area_index` as a (possibly time-varying) input. Prior to this
  work it produced *no* `phenology_factor`, so pairing it with `PALADYNAutotrophicRespiration` left
  that input unsatisfied.
- `PALADYNPhenology` — currently a stub. `compute_phenology_factor` returns a hardcoded `1.0` and
  `compute_f_deciduous` returns `0.0` (evergreen). The real PALADYN scheme is GDD-based and thus
  history-dependent; a naive implementation would require storing previous timesteps' temperature
  inputs, which conflicts with the differentiable, continuous-time design of Terrarium.

This plan specifies how each scheme derives `ϕ` in a manner consistent with the project rules:
continuous-time dynamics only, no discrete/per-day updates, no stored input history, GPU- and
Enzyme-compatible.

## Background

### PALADYN phenology (Willeit & Ganopolski, 2016; after Sitch et al., 2003)

Phenology type is set by the coldest-month temperature: if it falls below a PFT-specific threshold
`T_phen_cmon`, the PFT is cold-deciduous, otherwise evergreen. Needleleaf trees are always evergreen.
For deciduous PFTs:

- `LAI = ϕ · LAI_b`  (balanced/maximum LAI `LAI_b`).                                        (eq. 82)
- `ϕ = gdd / gddcrit`, where `gdd` accumulates degree-days above a base temperature `T_gdd_base`
  at a PFT-specific critical value `gddcrit`.                                                (eq. 83)
- After `ϕ` reaches its maximum of 1 it is held constant until air temperature drops below
  `T_gdd_base`. Senescence then proceeds, with all leaves lost at `T_gdd_base − 5 °C`.
- Raingreen (moisture-driven) phenology is not represented.

The scheme has three history-dependent elements: (1) the GDD accumulation, (2) the implicit annual
GDD reset used by LPJ/Sitch, and (3) the coldest-month-temperature classification. The senescence
ramp is an instantaneous function of current air temperature and needs no history.

### Design constraint

Terrarium prohibits discrete-time updates (e.g. "reset once per year", "sum once per day"). All
memory of the past must live in **prognostic state variables** integrated by the timestepper, and
all dynamics must be well-posed ODEs. This is what lets us discard stored input history.

## Summary of changes

### Part A — `PrescribedPhenology` derives the phenology factor

Derive `ϕ` from the prescribed LAI by inverting eq. 82 for the deciduous case, `ϕ = LAI / LAI_b`,
using a configured reference (annual-maximum) LAI:

```
ϕ = clamp(leaf_area_index / max_leaf_area_index, 0, 1)
```

Implementation:

- Add a `max_leaf_area_index` parameter (documented explicitly as the *reference / annual-maximum*
  LAI, not a generic constant — the seasonal amplitude of `ϕ` is entirely set by how well it matches
  the true annual maximum of the prescribed series).
- Emit `phenology_factor` as an `auxiliary` variable so `PrescribedPhenology` becomes a drop-in
  producer of `{leaf_area_index, phenology_factor}`, matching `PALADYNPhenology`'s outputs.
- The derivation is a pure algebraic map of a single cell's input, differentiable except for the
  `clamp` kinks (acceptable; subdifferentiable).

This part is essentially already drafted in `src/processes/vegetation/phenology/prescribed_phenology.jl`
via the `kernel(phenology_factor)` auxiliary-variable constructor. See "Testing and verification"
for a wiring check that must be completed before this is considered done.

### Part B — Prognostic GDD reformulation of `PALADYNPhenology`

Introduce a prognostic growing-degree-day state variable and derive `ϕ` diagnostically.

**Prognostic GDD (growth + continuous reset):**

```
d(gdd)/dt = max(T_air − T_gdd_base, 0) − (gdd / τ) · 𝟙[T_air < T_gdd_base]
```

with `t` in days so `gdd` carries units of degree-days (matching `gddcrit`). The first term is the
standard heat sum. The second term replaces LPJ's discrete annual reset with a continuous cold-season
leak (timescale `τ`, order weeks): `gdd` drains toward zero over winter so each spring restarts the
linear ramp from ≈0, keeping the variable bounded and event-free.

**Diagnostic phenology factor (unifies growth cap and senescence ramp):**

```
ϕ_deciduous = min( gdd / gddcrit,
                   clamp((T_air − (T_gdd_base − 5)) / 5, 0, 1) )
```

- Growth (`T_air ≥ T_gdd_base`): temperature term saturates at 1 ⇒ `ϕ = gdd/gddcrit` (eq. 83).
- Mature (`gdd/gddcrit ≥ 1`, warm): `ϕ = 1`, held.
- Senescence (`T_air < T_gdd_base`): temperature term falls below 1 and draws `ϕ` linearly to 0 at
  `T_gdd_base − 5 °C`.

**Deciduous/evergreen classification:** the coldest-month-temperature test is an annual statistic
that does not reduce to a clean ODE. Represent `f_deciduous` as a prescribed PFT/site property (or an
infrequently-updated input), consistent with the current `compute_f_deciduous` stub and with the fact
that needleleaf PFTs are always evergreen. Combine as in the existing `compute_LAI`:

```
ϕ_effective = f_deciduous · ϕ_deciduous + (1 − f_deciduous) · 1
LAI = ϕ_effective · LAI_b        (evergreen: ϕ_effective → 1, LAI = LAI_b)
```

**New parameters on `PALADYNPhenology`** (replacing the `# TODO add phenology parameters` stub),
each PFT-specific in PALADYN:

- `T_gdd_base` — base temperature for degree-day accumulation and senescence onset (°C).
- `gddcrit` — critical degree-day sum for full leaf-out (K·day).
- `τ` — cold-season GDD relaxation timescale (day).
- `f_deciduous` — prescribed deciduous fraction (or promote to an input).

**State/variable changes:** declare `gdd` as `prognostic(:growing_degree_days, XY())`, add a
`compute_tendencies!(state, grid, phenol::PALADYNPhenology, ...)` dispatch and a
`compute_gdd_tendency` kernel function (following the `PALADYNVegetationDynamics` /
`compute_ν_tendency` pattern in `src/processes/vegetation/dynamics/vegetation_dynamics.jl`), and keep
`phenology_factor` and `leaf_area_index` as auxiliaries computed in `compute_auxiliary!`.

## Testing and verification

- **Wiring check for the `kernel(...)` auxiliary-variable API (Part A).** DONE — see the revision
  log. `PrescribedPhenology` was the first user of the `kernel(func)` constructor in a
  `variables(...)` declaration, so the full path had to be repaired (auxiliary dispatch, ctor field
  bound, grid unwrapping, location types, `clock` kwarg inversion). Verified end-to-end: the lazy
  `KernelFunctionOperation`-backed `phenology_factor` field computes `LAI / max_leaf_area_index`
  (clamped) when read inside the downstream autotrophic-respiration kernel.
- **Unit / regime tests (Part B).** Verify `ϕ` reproduces the three PALADYN regimes: linear ramp
  during warm-up, held at 1 when mature, linear decline to 0 between `T_gdd_base` and
  `T_gdd_base − 5`. Verify `gdd` drains over a synthetic cold period and re-ramps the following warm
  period (no cross-year ratcheting).
- **Conservation / bounds.** Assert `0 ≤ ϕ ≤ 1` and `gdd ≥ 0` across a seasonal forcing cycle.
- **Type stability & allocations.** `@code_warntype` / allocation checks on the new kernel functions.
- **Differentiability.** Enzyme adjoint test through `compute_tendencies!` for `gdd` and through the
  diagnostic `ϕ`, mirroring `test/differentiability`. If the kinks at `T_gdd_base` degrade the
  adjoint, substitute softplus for `max(·,0)` and a smooth logistic for the cold-season indicator.
- **Evergreen limit.** With `f_deciduous = 0`, confirm `ϕ_effective = 1` and `LAI = LAI_b`
  independent of temperature.

## Documentation changes

- Add / update the phenology doc page following the standard Overview / Implementations / Methods /
  Kernel functions structure, including non-canonical docstrings for `PrescribedPhenology` and
  `PALADYNPhenology` and their `compute_auxiliary!` / `compute_tendencies!` dispatches.
- Document `max_leaf_area_index` as the reference (annual-maximum) LAI, and note the evergreen caveat
  for `PrescribedPhenology` (below).
- Cite eqs. 82–83 of Willeit & Ganopolski (2016) and Sitch et al. (2003) in the PALADYN section.
- Use `jldoctest` blocks and unicode math per project rules.

## Known limitations

- **`PrescribedPhenology` assumes deciduous-like behavior.** `ϕ = LAI/LAI_max` conflates a
  *structural* LAI ratio with a *phenological* leaf-out state. For an evergreen PFT (where PALADYN
  sets `ϕ = 1` regardless of LAI), a dip in prescribed LAI would spuriously suppress root
  respiration. Document this; if evergreen support is needed, allow `phenology_factor` to be a second
  prescribed input or add an evergreen flag.
- **Spring-frost clipping (Part B).** The diagnostic `min(...)` has no phase memory, so a hard frost
  during spring green-up would transiently clip `ϕ` where strict LPJ (a state machine) would not.
  Judged physically acceptable (frost damage) and avoids carrying a discrete phase flag, but it is a
  deviation from the reference.
- **Coldest-month classification is not dynamic.** Treating `f_deciduous` as prescribed omits online
  biogeographic reclassification. Acceptable given needleleaf-always-evergreen and climate-timescale
  variation; revisit if dynamic vegetation requires it.
- **`gdd_relaxation_time` is a modeling choice**, not from the PALADYN paper; it emulates the annual
  reset and needs tuning so a full winter drains `gdd` but short warm spells do not.
- **`growing_degree_days` is unbounded in a permanently-warm deciduous climate** (no cold season to
  trigger the relaxation). The `min(gdd/gdd_crit, …)` cap means this has no effect on `ϕ`, but the
  raw accumulator grows without bound over very long runs; acceptable for now.
- **Parameter defaults are placeholders.** `T_gdd_base`, `gdd_crit`, `T_senescence_range` and the
  deciduous classification are PFT-specific (PALADYN Table 5) and must be calibrated per PFT.
- **Bare-`VegetationModel` timestepping is blocked by an unrelated, pre-existing assertion**
  (`@assert isfinite(LAI)` in `medlyn_stomatal_conductance.jl`) that fires when `carbon_vegetation`
  is zero-initialized, making `balanced_leaf_area_index` non-finite. This is independent of the
  phenology changes (GDD-tendency correctness was verified directly) but should be revisited — it is
  also a throw path in a kernel, which the kernel rules disallow.

## Future work / remaining TODOs observed nearby

These are pre-existing TODOs in the vegetation code, surfaced here for tracking; they are **out of
scope** for this plan unless noted:

- `paladyn_phenology.jl` — replace the stubbed `compute_f_deciduous` (returns 0) and
  `compute_phenology_factor` (returns 1) with the implementation above; add the phenology parameters.
- `respiration/autotrophic_respiration.jl`:
  - `cn_sapwood`, `cn_root` — physical meaning and units unconfirmed (`# TODO check ... + add unit`).
  - `compute_f_temp` hardcodes the constants `308.56`, `56.02`, `46.02`; per the rule against
    hardcoded literals these should move into the process struct or `PhysicalConstants`.
  - The `T_soil > 7 °C` hard bound is inherited from CLIMBER-X/PALADYN and is not further justified;
    flagged as a candidate for improvement / data-driven replacement.
  - `compute_resp10` returns a placeholder constant `0.066`; needs a real implementation, meaning,
    and units.
- `dynamics/vegetation_dynamics.jl` & `dynamics/carbon_dynamics.jl` — several parameters are annual
  and marked to be converted to daily; the disturbance rate `γv` is a placeholder minimum awaiting
  the soil-moisture-dependent PALADYN formulation; leaf/stem/root carbon splitting is not yet
  implemented.
