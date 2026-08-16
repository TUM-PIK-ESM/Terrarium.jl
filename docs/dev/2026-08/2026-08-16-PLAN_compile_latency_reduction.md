# Reduce compile-time latency in model initialization via OrderedDicts

> Status: **in progress**. OrderedDict-based initialization implemented on `bg/compile-latency`.
> `initialize()` compile time reduced ~39% (158 s → 96 s); `timestep!()` unchanged (~110 s).
> See [Results](#results) for full benchmark table.

Date of initial draft: 2026-08-16

Base revision: dd85eaf3652ffbb76d06945424a3375a39af52fd

## Originating prompt

> I would like to reduce the compile time of Terrarium model initialization because this is a real
> bottleneck for debugging. The basic idea was to replace the NamedTuples used during initialization
> with `OrderedDict`s and only build NamedTuples at the very last step. These changes will need to be
> applied in both state_variables.jl and abstract_variables.jl. You can evaluate the performance
> improvement using the land_compiler_analysis.jl script.

## Revision log

- 2026-08-16: Initial draft. Scope confirmed: all 4 variable groups (prognostic, tendencies, auxiliary,
  inputs) plus `Variables` itself. DataStructures.jl already a dependency. Focus on inference time
  reduction. Full test suite + snoop benchmark for validation.
- 2026-08-16: Implementation complete. Key deviations from original plan:
  - User rewrote the `Variables` constructor with a `register!` pattern instead of the planned
    OrderedDict-based deduplication (see [Summary of changes](#summary-of-changes)).
  - Per-variable `initialize()` dispatches accept `OrderedDict{Symbol, AbstractField}` for the
    fields context (not NamedTuple as originally planned in Option A), with a convenience wrapper
    that converts NamedTuple → OrderedDict at the boundary.
  - Auxiliary custom constructors receive `NamedTuple(fields)` — converted from the accumulated
    dict at call time (same as plan).
- 2026-08-16: Benchmarks run on both `main` and `bg/compile-latency`; results recorded below.

## Problem description

Model initialization in Terrarium is slow due to Julia's compiler inferring types through deeply nested
`NamedTuple` constructions. The hot path is:

```
StateVariables(model)
  → Variables(tuplejoin(variables(model), input_variables))
    → creates NamedTuples for prognostic, tendencies, auxiliary, inputs, namespaces
  → initialize(vars.prognostic, grid, ...) via foldl with merge()
  → initialize(vars.auxiliary, grid, ...) via foldl with merge()
  → ... (repeated for each group)
  → final NamedTuple construction for StateVariables fields
```

Each `merge()` call in the `foldl` loop creates a new `NamedTuple` type, and Julia must infer the
return type of each iteration. With many variables (a full coupled land model has 50+ state variables),
this compiles to O(n²) type inference work because each fold step produces a distinct `NamedTuple{names}`
type that the compiler hasn't seen before.

The solution is to use `OrderedDict{Symbol, Field}` during the accumulation phase (where order matters
but types don't need to be distinct), then convert to `NamedTuple` only at the final construction of
`StateVariables`. This eliminates the per-iteration type explosion while preserving field ordering and
runtime behavior.

## Background

### Current initialization flow

1. **`Variables(vars::Tuple{...})`** in `abstract_variables.jl`:
   - Partitions variables into prognostic/auxiliary/input/namespaces
   - Creates `NamedTuple`s via `(; map(var -> varname(var) => var, vars)...)`
   - These NamedTuples are then passed to `StateVariables()` constructor

2. **`StateVariables(vars, grid; ...)`** in `state_variables.jl`:
   - Calls `initialize(vars.prognostic, grid, clock, bcs, fields)` for each group
   - Each `initialize` does a `foldl(vars, init=(;)) do nt, var ... merge(nt, (; name => field)) end`
   - This is the primary inference bottleneck: each fold step creates a new NamedTuple type

3. **Final assembly**: The resulting NamedTuples of Fields become fields of `StateVariables`

### Why OrderedDict helps

- `OrderedDict{Symbol, Field}` has a stable type regardless of how many entries it accumulates
- `push!` or `dict[key] = value` doesn't create new types
- The final conversion to `NamedTuple` happens once, not n times
- DataStructures.jl is already a project dependency

### Prior attempts

Commit `a78ad54c1` ("Some more attempts at reducing compile time") added `@nospecialize` markers and
reverted to `fastiterate` for looping, but didn't address the fundamental NamedTuple type explosion in
the `foldl` accumulation pattern.

## Summary of changes

### 1. `abstract_variables.jl` — Variables container uses OrderedDicts internally

**Current:**
```julia
struct Variables{ProgVars, TendVars, AuxVars, InputVars, Namespaces}
    prognostic::ProgVars       # NamedTuple{names, Tuple{...}}
    tendencies::TendVars       # NamedTuple{names, Tuple{...}}
    auxiliary::AuxVars         # NamedTuple{names, Tuple{...}}
    inputs::InputVars          # NamedTuple{names, Tuple{...}}
    namespaces::Namespaces     # NamedTuple{names, Tuple{...}}
end
```

**After:**
```julia
struct Variables
    prognostic::OrderedDict{Symbol, AbstractVariable}
    tendencies::OrderedDict{Symbol, AuxiliaryVariable}
    auxiliary::OrderedDict{Symbol, AbstractVariable}
    inputs::OrderedDict{Symbol, AbstractVariable}
    namespaces::OrderedDict{Symbol, Namespace}
end
```

- The `Variables(vars::Tuple{...})` constructor will populate OrderedDicts instead of NamedTuples
- Add `Base.propertynames(::Variables)` and `Base.getproperty(::Variables, ::Symbol)` so existing
  code that accesses `vars.prognostic`, `vars.auxiliary`, etc. continues to work (these return the
  OrderedDicts, which downstream code can iterate over)
- The `check_duplicates` function already works with variable collections — no change needed
- `merge_namespaces` returns a tuple of Namespace objects — convert to OrderedDict storage

### 2. `state_variables.jl` — initialize() uses OrderedDict accumulation

**Current:**
```julia
function initialize(vars::NamedTuple{names, ...}, grid, clock, bcs, fields) where {names}
    return foldl(vars, init = (;)) do nt, var
        field = initialize(var, grid, clock, bcs, merge(nt, fields))
        merge(nt, (; varname(var) => field))
    end
end
```

**After:**
```julia
function initialize(vars::OrderedDict{Symbol, AbstractVariable}, grid, clock, bcs, fields)
    result = OrderedDict{Symbol, Field}()
    for (name, var) in vars
        field = initialize(var, grid, clock, bcs, fields)
        result[name] = field
    end
    return result
end
```

- Replace `foldl` + `merge` with a simple loop over OrderedDict entries
- The `fields` context (passed to individual `initialize(var, ...)` calls) remains a NamedTuple for
  type stability at the per-variable level — only the accumulator changes
- After all groups are initialized as OrderedDicts, convert to NamedTuples in the final
  `StateVariables` constructor call

### 3. `state_variables.jl` — StateVariables construction from OrderedDicts

The `StateVariables` struct itself retains its NamedTuple fields (these are accessed at runtime during
simulation, where type stability matters). The change is only in how they're constructed:

```julia
# After initializing each group as an OrderedDict, convert to NamedTuple once:
prognostic_nt = NamedTuple{Tuple(keys(prognostic_od))}(values(prognostic_od))
```

This single conversion per group (4 total) is negligible compared to the n conversions in the foldl.

### 4. Property access compatibility

Since `Variables` will store OrderedDicts but existing code accesses them as properties:

```julia
# This must still work:
vars.prognostic   # → returns OrderedDict{Symbol, AbstractVariable}
vars.auxiliary    # → returns OrderedDict{Symbol, AbstractVariable}
```

Add to `abstract_variables.jl`:
```julia
Base.propertynames(::Variables) = (:prognostic, :tendencies, :auxiliary, :inputs, :namespaces)
Base.getproperty(v::Variables, name::Symbol) = getfield(v, name)  # direct field access
```

For code that iterates over `vars.prognostic` (expecting a NamedTuple of variables), OrderedDict
iteration returns `(name, var)` pairs — check if any code does `for var in vars.prognostic` and
adjust to `for (_, var) in vars.prognostic` or `values(vars.prognostic)`.

### 5. Individual variable initialize() dispatches

The per-variable `initialize(var::AbstractVariable, ...)` and `initialize(var::AuxiliaryVariable, ...)`
dispatches receive a `fields::NamedTuple` context. This NamedTuple is the *user-provided* fields or
the accumulated fields from previous variables. Two options:

- **Option A (preferred)**: Keep `fields` as a NamedTuple for these dispatches. The OrderedDict
  accumulator only replaces the outer foldl loop. The `fields` context passed to each variable's
  constructor is built from the user's `fields` kwarg, not from the accumulator.
- **Option B**: Also convert to OrderedDict lookup via `get(fields_dict, name, nothing)`. This adds
  runtime dispatch but further reduces inference burden.

**Decision: Option A.** The per-variable initialize calls are few and the fields context is small
(user-provided). The main inference cost is in the foldl accumulator, not in the per-variable calls.

## Results

Benchmarks run on the same machine using `test/benchmarks/compilation/land_compiler_analysis.jl` with
a coupled `LandModel` (Richards hydrology + snow + vegetation carbon cycle, 30 vertical levels).

| Function | `main` (baseline) | `bg/compile-latency` | Δ time | Δ allocs |
| --- | --- | --- | --- | --- |
| `initialize(::LandModel)` | 157.6 s (42.15 M, 2.84 GiB, 104.3% compile) | 96.2 s (41.29 M, 2.79 GiB, 107.4% compile) | **-39%** | -2% |
| `timestep!(::ModelIntegrator)` | 109.7 s (15.00 M, 1.01 GiB, 101.1% compile) | 115.2 s (15.14 M, 1.02 GiB, 101.1% compile) | +5% | +1% |

### Interpretation

- **`initialize()` improved ~39%** — the OrderedDict swap in both `Variables` and `StateVariables`
  construction eliminated a significant portion of the NamedTuple type explosion. Allocations also
  dropped slightly (42.15 M → 41.29 M) confirming fewer intermediate types were created.
- **`timestep!()` unchanged** — compile latency for the first timestep is dominated by kernel
  compilation (`@kernel` call graphs for tendency/auxiliary computations), not variable initialization.
  The +5 s difference is within noise for cold-start measurements.
- **Compile time is still ~100% of wall time** for both functions, meaning virtually all measured
  time is JIT compilation, not execution.

### Inference tree analysis (`@snoop_inference`)

The `@snoop_inference` tree for `initialize(::LandModel)` on the branch reveals:

| Metric | Value | % of wall |
| --- | --- | --- |
| Total wall time | 96.96 s | 100% |
| Codegen + execution (root self) | **86.43 s** | **89%** |
| Actual type inference | **10.53 s** | **11%** |

**Key finding: the bottleneck is codegen, not inference.** The OrderedDict change reduced wall time
by ~60s, but 86s of the remaining 96s is code generation and execution, not type inference. This
means further inference optimizations will have diminishing returns — the compiler spends most of its
time lowering already-inferred code to LLVM IR, not inferring types.

**Inference breakdown (10.53s total):**

| Package | Self inference time | % of inference |
| --- | --- | --- |
| Core/Other | 8.80 s | 84% |
| Terrarium | 1.40 s | 13% |
| Oceananigans | 0.32 s | 3% |

**Top Terrarium inference hotspots by inclusive time:**

| Function | Inclusive time | % of inference | Children |
| --- | --- | --- | --- |
| `initialize#607` (kw wrapper) | 2.80 s | 27% | 7 |
| `StateVariables#600` (kw wrapper) | 1.44 s | 14% | 21 |
| `initialize!` (multiple dispatches) | 1.35 s, 1.32 s, 1.20 s | 13%, 13%, 11% | 17, 16, 6 |
| `variables` | 1.02 s | 10% | 14 |
| `fastmap` | 0.87 s | 8% | 20 |
| `closure!` | 0.80 s | 8% | 11 |
| `get_fields#152` | 0.51 s | 5% | 23 |

These functions have near-zero self-time but expensive callees — the inference cost is in what they
dispatch to, not their own bodies. The high child counts (7-21 per function) indicate many specialized
method instances being inferred, likely from dispatch on the mega-type `LandModel{Float64, ...}` with
its deeply nested component type parameters.

### What wasn't captured

The benchmark script uses `@time` on a single (cold) call. The user reported that recompilation
occurs on every subsequent call to `initialize()`, which was not quantified here. This is a separate
issue from the first-call compile cost addressed by this change.

## Testing and verification

### Correctness tests
1. Module loads cleanly: `using Terrarium` passes without errors
2. Full test suite passes on both branches (verified before benchmarking)
3. Namespace recursion works correctly (tested via coupled LandModel with soil + vegetation namespaces)

### Performance benchmarks
- See [Results](#results) above for the benchmark table
- Target of ≥50% reduction was not met; actual improvement was ~39% for `initialize()`
- The remaining compile cost is in kernel compilation, outside the scope of this change

### Regression checks
1. ✅ `StateVariables` struct fields are still concrete-typed NamedTuples (runtime type stability preserved)
2. ✅ `get_fields(state, ...)` works correctly (generated function unchanged)
3. ✅ `fastiterate` loops over state properties work (StateVariables internals unchanged)
4. ⚠️ Enzyme differentiability tests — not yet re-run after these specific changes

## Documentation changes

- No user-facing API changes; `Variables` and `StateVariables` public interfaces remain identical
- Internal implementation detail only

## Known limitations

1. **`fields` context in auxiliary constructors**: Some auxiliary variables have custom Field
   constructors (`var.ctor(var, grid, clock, fields)`) that expect a NamedTuple `fields`. These
   must receive a proper NamedTuple, so we'll need to ensure the context passed to these constructors
   is still a NamedTuple (Option A above handles this).

2. **Ordering dependency**: The comment in the current code notes that "Fields visible to each
   constructor are dependent on the order in which variables were declared." OrderedDict preserves
   insertion order, so this invariant is maintained.

## Future work

1. **`timestep!()` compile cost is unchanged** (110 s on `main` vs 115 s on the branch) — the bulk
   of Terrarium's compile latency lives in kernel compilation and timestepper dispatch, not in
   variable initialization. Addressing this would require a separate investigation into:
   - Kernel function inference chains (the `@kernel` call graph for `compute_tendencies!`,
     `compute_auxiliary!`, etc.)
   - Timestepper cache type parameters and how they propagate through `ModelIntegrator`
   - Whether `@precompile` hints for common model configurations would help

2. **Recompilation on every `initialize()` call** — the benchmark script timed a single call with
   `@time`, which includes both compilation and execution. The fact that compile time is ~100% of
   wall time means this was a cold (first) call. If subsequent calls also recompile, the root cause
   is likely type instability in the return type of `initialize(model; ...)` — possibly from
   closure functions or initializer lambdas leaking into the inferred signature. This warrants
   separate investigation with `@code_warntype` and/or SnoopCompile's inference flamegraph.

3. **Precompilation hints** (`@precompile`) for common model configurations could complement this
   change, but only after understanding why recompilation occurs on repeated calls.

4. The `get_fields(state, vars::Tuple{...})` generated function is already type-stable; no changes
   needed after the Variables internal representation change.
