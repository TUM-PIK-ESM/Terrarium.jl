# Reactant support for `VegetationModel` with `PrescribedVegetation` (static LAI)

> Status: **in progress**. Exploration and blocker identification complete and reproduced in a
> Terrarium-free-ish MWE; no fixes implemented yet, and the design decisions below need sign-off. Scope
> is deliberately limited to a *static, constant* prescribed LAI — the time-varying LAI climatology is
> known to be blocked by input interpolation and is explicitly out of scope for this step.

Date of initial draft: 2026-08-14

Base revision: 99c748b79711e295e99a5a6370d386853652cf06

## Originating prompt

> We continue with the Reactant adjustment, but we move to a new model. Also start a new markdown
> document for this. We want to make the `VegetationModel` with `PrescribedVegetation` Reactant
> compatible. However, in a first step we want to do this with just a static constant input and not
> with a temporally evolving prescribed LAI climatology (this will be blocked by the interpolation not
> Reactant compatible)

## Revision log

- 2026-08-14: initial draft, written during exploration and before the blocker list was confirmed by
  running the model.
- 2026-08-14: blockers confirmed by running. The suspected data-dependent-branch problem in the
  photosynthesis kernel could **not** be reached — initialization fails first, so the traced step has
  never executed. Two initialization-path blockers are confirmed instead, both rooted in lazy
  `AbstractOperation`s, and reproduced in
  `2026-08-14-reactant_abstract_operation_mwe.jl`. The `:vegetation_column` test configuration has been
  added to `test/reactant/setup.jl` (not yet registered in `runtests.jl`, since it does not pass).
- 2026-08-14: blocker #2 narrowed to a pure Oceananigans reproducer and handed off for an upstream fix.
  The initial characterisation here ("`compute!` on an operation-backed `Field` fails under Reactant")
  was too broad: it only fails when the grid's **vertical coordinate is array-valued**, which Terrarium's
  `ExponentialSpacing` column grids always are. See
  `2026-08-14-OCEANANIGANS_computed_field_stretched_grid_issue.md` and
  `2026-08-14-oceananigans_computed_field_mwe.jl`.

## Problem description

`test/reactant/` covers standalone `SoilModel` and `SnowModel` configurations and the coupled
`LandModel` (soil + snow, no vegetation; see `2026-08-04-PLAN_reactant_coupled_land_model.md`). The
vegetation side has never been traced — that plan lists "Vegetation is still unsupported" as a known
limitation, pointing at a root-fraction `exp(::TracedRArray)` blocker.

The goal here is the *narrowest useful* vegetation configuration: a standalone
`VegetationModel(grid; vegetation = PrescribedVegetation(NF))` whose leaf area index is a **constant**
input field, driven by the default `PrescribedAtmosphere`.

### Explicitly out of scope: the LAI climatology

`examples/simulations/vegetation_global.jl` drives the same model with an annually cycling ERA5-Land
LAI climatology:

```julia
lai_input = InputSource(grid, lai_highveg; source_grid = ..., name = :leaf_area_index, cycle = true)
integrator = initialize(model; inputs = lai_input)
```

That path performs temporal interpolation of the input source on every `update_inputs!`, which does not
trace under Reactant. Making it work is a separate problem (and a separate document). This step
prescribes LAI as a constant instead, via the `initializers` mechanism already used for constant inputs
elsewhere in the suite (e.g. `air_temperature = NF(-2)` in `:land_soil_snow`).

## Background

### What `PrescribedVegetation` actually runs

`PrescribedVegetation` (`src/processes/vegetation/prescribed_vegetation.jl`) composes six components:

| Component | Default | Role |
| --- | --- | --- |
| `phenology` | `PrescribedPhenology` | LAI is an *input*; derives `phenology_factor` |
| `photosynthesis` | `LUEPhotosynthesis` | `net_assimilation`, `leaf_respiration`, `gross_primary_production` |
| `stomatal_conductance` | `MedlynStomatalConductance` | `canopy_water_conductance` |
| `root_distribution` | `StaticExponentialRootDistribution` | `root_fraction` (field constructor) |
| `plant_available_water` | `FieldCapacityLimitedPAW` | `plant_available_water`, `soil_moisture_limiting_factor` |
| `traits` | `PlantTraits` | scalar PFT parameters |

`compute_tendencies!` is a **no-op** — nothing is prognostic. Everything happens in
`compute_auxiliary!`, which under `VegetationModel` is called *without* a soil component:

```julia
function compute_auxiliary!(state, model::VegetationModel)
    (; grid, atmosphere, vegetation, constants) = model
    compute_auxiliary!(state, grid, vegetation, constants, atmosphere)   # note: no `soil`
    return nothing
end
```

so `soil === nothing` in `compute_auxiliary!(state, grid, veg::PrescribedVegetation, …)`. That matters:
both soil-dependent components short-circuit to no-ops during stepping

- `compute_auxiliary!(state, grid, ::StaticExponentialRootDistribution, args...) = nothing`
- `compute_auxiliary!(state, grid, ::AbstractPlantAvailableWater, soil::Nothing, args...) = nothing`

leaving only two kernels in the traced step: `LUEPhotosynthesis` and `MedlynStomatalConductance`
(`phenology`'s `compute_auxiliary!` is also a no-op; `phenology_factor` is a `kernel(...)`-backed
auxiliary).

This splits the problem cleanly:

- **Eager (initialization) path**: `root_fraction` and `soil_moisture_limiting_factor` field
  constructors, and `initialize!` for PAW. These run once, on the device, outside tracing.
- **Traced (stepping) path**: the photosynthesis and stomatal-conductance kernels.

### Confirmed blockers

Reproduced by building `:vegetation_column` on a `ReactantState` grid. **Initialization fails, so the
traced stepping path has never run** — the blockers below are all on the eager initialization path, and
further blockers (see "Not yet reachable") are expected once these are cleared.

All cases are reproduced by `docs/dev/2026-08/2026-08-14-reactant_abstract_operation_mwe.jl`:

```
A. sum(a * b, dims=3)                        [Fields only]                        OK
B. sum(a * b / Δz, dims=3)                   [+ Δz KernelFunctionOperation]       DimensionMismatch
C. Field(Integral(a * b, dims=3)) |> compute!                                     DimensionMismatch
D. Field(Integral(a * b / Δz, dims=3)) |> compute!  [= soil_moisture_limiting_factor] DimensionMismatch
E. Field(Δz) |> compute!                     [Δz alone]                           InvalidIRError
F. Field(a * b) |> compute!                  [elementwise, no reduction]          InvalidIRError
G. root_fraction(grid, …)                    [FunctionField + sum]                MethodError: exp(::TracedRArray)
```

#### 1. Reactant materializes lazy operations by pushing *whole arrays* through scalar `getindex`

This is the root cause of B, C, D and G. `root_fraction`'s `sum(R, dims = 3)` reaches Reactant's
`overloaded_mapreduce`, which calls `TracedUtils.materialize_traced_array` on the lazy
`BinaryOperation`. That evaluates the operation's *scalar* `getindex` chain with whole traced arrays
instead of scalar indices, so the `FunctionField` closure receives `z::TracedRArray{Float32, 3}`:

```
MethodError: no method matching exp(::Reactant.TracedRArray{Float32, 3})
 [1] root_density        @ src/processes/vegetation/hydraulics/root_distribution.jl:39
 [2] (::var"#root_fraction##0#1")(x::StepRangeLen{…}, z::Reactant.TracedRArray{Float32, 3}, params::…)
 [3] call_func           @ Oceananigans/src/Fields/function_field.jl:52
 [4] getindex            @ Oceananigans/src/Fields/function_field.jl:56
 [8] materialize_traced_array(x::BinaryOperation{…})   @ Reactant.TracedUtils
 [9] overloaded_mapreduce(…)                            @ Reactant.TracedRArrayOverrides
```

`exp` on a `TracedRArray` is *correctly* undefined: `exp(::Matrix)` in Julia means the matrix
exponential, not the elementwise one, so Reactant declines to guess. The bug is upstream of that — the
scalar closure should never have been handed an array.

The same mechanism breaks `Δz`: a `KernelFunctionOperation` whose `getindex` reads grid internals, which
under materialization receives `Base.OneTo` ranges for `i, j, k` and fails with `DimensionMismatch`
inside `getspacing_3d`. Case A is the control: with every operand a materialized `Field`, plain
elementwise math *does* survive materialization, so the reduction machinery itself is not the problem.

Case C shows `Integral` (a `Scan`) fails even without `Δz`, so `soil_moisture_limiting_factor` is
blocked on two independent counts.

#### 2. `compute!` on an operation-backed `Field` does not compile on a stretched grid

Cases E and F. `compute!` launches `gpu__compute!` on the Reactant KA backend and fails:

```
InvalidIRError: compiling MethodInstance for Oceananigans.AbstractOperations.gpu__compute!(…)
Reason: unsupported call to an unknown function (call to jl_f_throw_methoderror)
```

This fails for a *plain elementwise* `Field(a * b)` — no reduction, no `Δz`, no `FunctionField` — so it
is independent of blocker #1 and will not be fixed by fixing it.

Narrowed down in a pure Oceananigans MWE (`2026-08-14-oceananigans_computed_field_mwe.jl`): the trigger
is an **array-valued (stretched) vertical coordinate**, not the operation itself. The same code works on
the CPU, works on `ReactantState` with a uniform `z = (-1, 0)`, and works when the arithmetic is done by
broadcasting instead of through an `AbstractOperation`. In the compiled kernel signature the grid's
float parameter has become `CuTracedRNumber{Float32, 1}` and the operation's eltype
`Reactant.TracedRNumber{Float32}`, while the destination array is plain `Float32` — hence the
unresolvable store. Full write-up, including suggested starting points for a fix, in
`2026-08-14-OCEANANIGANS_computed_field_stretched_grid_issue.md`.

Terrarium calls `compute!` on an operation-backed field in exactly one place — `FieldCapacityLimitedPAW`
(`plant_available_water.jl:61,67`, in both `compute_auxiliary!` and `initialize!`) — which is why the
soil and snow models never hit this and the vegetation model does immediately. Terrarium's column grids
use `ExponentialSpacing`, so they are always on the failing side of the trigger.

**This one is being taken upstream** (2026-08-14). The local kernel-based workaround in "Proposed
approach" below removes Terrarium's dependence on it either way.

### Not yet reachable: data-dependent `if`/`else` in the photosynthesis kernel

`compute_respiration_assimilation` (`lue_photosynthesis.jl`) branches on traced values, with different
code paths and different variables assigned per branch:

```julia
if swdown > zero(NF) && T_air > NF(-3.0)
    ...
    if LAI > zero(NF)
        ...  # ~8 derived quantities, then An = Ag - Rd
    else
        An = zero(NF); Rd = zero(NF)
    end
else
    An = zero(NF); Rd = zero(NF)
end
```

and `compute_temperature_stress` does the same:

```julia
if photo.T_CO2_low < T_air < photo.T_CO2_high
    low  = NF(1) / (NF(1) + exp(k1 * (k2 - T_air)))
    high = NF(1) - NF(0.01) * exp(k3 * (T_air - photo.T_photos_high))
    T_stress = low * high
else
    T_stress = zero(NF)
end
```

This is the AGENTS.md rule "use `ifelse` — never short-circuiting `if`/`else` in kernels". The
short-circuiting `&&` and the nested branches are the concern; whether they actually fail to raise or
merely produce poor IR needs to be measured, since (unlike the `LandModel` masked-store case) there are
no *stores* in the branches, only value selection.

Note the rewrite is not purely mechanical: `ifelse` evaluates **both** arms, and the "no leaves" arm
exists precisely because the lit-and-leafy arm divides by quantities that are degenerate when
`LAI == 0` or `swdown == 0` (e.g. `Vc_max` divides by `pres_i - Γ_star`; `compute_λc` and
`compute_stomatal_conductance` take `sqrt(vpd)`). Guard values must be clamped so both arms are total,
as was done for the van Genuchten conductivity in blocker #2 of the `LandModel` plan.

#### 2. `root_fraction`: `FunctionField` + lazy reduction, evaluated on device

`root_fraction` (`root_distribution.jl`) is a *field constructor* — the declared auxiliary is the lazy
object it returns:

```julia
function root_fraction(grid::AbstractColumnGrid, clock, fields, rootdist::StaticExponentialRootDistribution{NF}) where {NF}
    fgrid = get_field_grid(grid)
    ∂R∂z = FunctionField{Center, Center, Center}(fgrid, parameters = rootdist) do x, z, params
        root_density(params, z)      # NF(0.5) * (a * exp(a * z) + b * exp(b * z))
    end
    Δz = zspacings(fgrid, Center(), Center(), Center())
    R = ∂R∂z * Δz
    R_norm = R / sum(R, dims = 3)    # lazy AbstractOperation + reduction
    return R_norm
end
```

This is the `exp(::TracedRArray)` blocker referenced in the `LandModel` plan (case G above).

Mitigating factor: with `soil === nothing` this is only touched at initialization, which is eager, and
the distribution is *static* — it depends only on the grid and two scalar parameters, never on state.

#### The other affected site: `soil_moisture_limiting_factor`

```julia
function soil_moisture_limiting_factor(grid, clock, fields, ::FieldCapacityLimitedPAW)
    Δz = zspacings(get_field_grid(grid), Center(), Center(), Center())
    β = Integral(fields.plant_available_water * fields.root_fraction / Δz, dims = 3)
    return Field(β)
end
```

Hit by blocker #1 twice (`Integral`/`Scan`, and `Δz`) and by blocker #2 (`compute!` is called on it from
both `initialize!` and `compute_auxiliary!`). With no soil the inputs are trivial — PAW is `set!` to 1 —
but the field is still constructed and computed.

## Proposed approach

Both blockers are on the initialization path and both stem from expressing these two quantities as lazy
`AbstractOperation`s. The natural fix for Terrarium is to stop doing that and compute them with kernels
into plain `Field`s, which is also what AGENTS.md prescribes ("Never loop over grid points outside
kernels — use `launch!`") and is cheaper besides: `root_fraction` is *static*, yet is currently
re-evaluated lazily on every access.

**Option A — materialize both into `Field`s computed by kernels.** `root_fraction` becomes a plain
`XYZ` auxiliary filled by an `XY` kernel that walks the column (two passes over `k`: accumulate the
normalizing sum, then write the normalized fractions; `Nz` is a compile-time constant, so the loop
unrolls). `soil_moisture_limiting_factor` likewise becomes an `XY` kernel that accumulates
`Σ_k PAW[i,j,k] * root_fraction[i,j,k]` directly, dropping the `Integral`/`Δz` operation entirely.
This sidesteps blockers #1 and #2 without needing anything upstream, and changes no numerics *if* the
same midpoint-rule discretization is kept.

**Option B — additionally use the analytic CDF for the root distribution.** The current code
approximates the layer integral by the midpoint rule (`∂R∂z * Δz`) and then normalizes. The
distribution has a closed-form CDF — for `∂R/∂z = ½(a·e^{az} + b·e^{bz})`, `R(z) = ½(e^{az} + e^{bz})` —
so the fraction in a layer is exactly `R(z_top) − R(z_bottom)`, normalized over the truncated column.
That is both more accurate and simpler, but it **changes results** for existing users, so it should be a
deliberate, separately-reviewed decision rather than a side effect of Reactant work.

Recommendation: **Option A**, keeping the current discretization so this stays a pure
Reactant-compatibility change, with Option B proposed separately if the accuracy improvement is wanted.

Whether blockers #1 and #2 should *also* be reported upstream is worth deciding independently — #2 in
particular (`compute!` on a plain elementwise operation-backed `Field` failing to compile) looks like a
general Oceananigans/Reactant gap that will affect other users, not something specific to Terrarium.

## Summary of changes

Not yet implemented — pending sign-off on the approach above. Anticipated:

| File | Change |
| --- | --- |
| `src/processes/vegetation/hydraulics/root_distribution.jl` | Compute `root_fraction` into a `Field` via a kernel instead of returning a lazy operation |
| `src/processes/vegetation/hydraulics/plant_available_water.jl` | Compute `soil_moisture_limiting_factor` via a kernel; drop `Integral`/`compute!` |
| `src/processes/vegetation/photosynthesis/lue_photosynthesis.jl` | Branchless (`ifelse`) rewrite with totalized arms, *if* the traced step proves to need it |
| `test/reactant/setup.jl` | New `:vegetation_column` configuration (constant LAI, default atmosphere) — **added** |
| `test/reactant/runtests.jl` | Register `:vegetation_column` in the testset |

## Testing and verification

`:vegetation_column` follows the existing registry pattern: a single `ColumnGrid` column on an
`ExponentialSpacing` grid, `PrescribedVegetation` with all-default components, constant
`leaf_area_index` and `air_temperature`, compared CPU-vs-Reactant across every field after N steps by
the shared `test_model` harness.

Because `compute_tendencies!` is a no-op, *nothing evolves*: every field is a pure function of the
(constant) inputs, so the comparison is effectively a check that the auxiliary kernels compute the same
values on both backends. This is a weaker test than `:land_soil_snow` and should be recognised as such —
it verifies compilation and pointwise numerics, not integration. It becomes a real dynamical test only
once a time-varying LAI input is supported (out of scope here) or the model is coupled into `LandModel`.

## Documentation changes

None expected for the model itself. If `root_fraction` or `soil_moisture_limiting_factor` change from
lazy to materialized, their docstrings must be updated to match (both currently document themselves as
returning a lazily-computed field).

## Known limitations

- **Time-varying LAI remains unsupported** — the `InputSource` interpolation path does not trace. This
  is the motivating restriction for this step, not an incidental one.
- **The test is static** — see the caveat under Testing and verification.
- **The traced stepping path is still unverified.** Initialization fails before any kernel is traced, so
  the photosynthesis and stomatal-conductance kernels have never been compiled under Reactant. Further
  blockers there are likely (see "Not yet reachable") and cannot be enumerated until the initialization
  blockers are fixed.

## Future work

- Reactant support for interpolated/time-varying `InputSource`s, which would let the ERA5-Land LAI
  climatology in `examples/simulations/vegetation_global.jl` run under Reactant.
- Extend to `VegetationCarbonCycle` (prognostic carbon pool, `PaladynPhenology`, autotrophic
  respiration, vegetation dynamics), which has real tendencies and would give a meaningful dynamical
  Reactant test.
- Wire vegetation into the coupled `LandModel` Reactant configuration once both halves work
  standalone.
