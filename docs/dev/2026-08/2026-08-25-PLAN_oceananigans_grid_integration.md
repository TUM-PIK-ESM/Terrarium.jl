# Oceananigans AbstractGrid Integration for Terrarium

> Status: **planned**. Draft integration plan for making `AbstractLandGrid` a subtype of Oceananigans `AbstractGrid` with support for multi-domain vertical discretizations (Ground, Snow, Canopy).

Date of initial draft: 2026-08-25

Base revision: 4f3841955af84ac6cb4a6f88b04538fc0dd7d658

## Originating prompt

> Please review the Oceananigans `AbstractGrid` interface and grid implementations. Then review the current Terrarium `AbstractLandGrid`s and draft a plan for how to integrate the two interfaces. The general idea is that `AbstractLandGrid` should be a valid subtype of `AbstractGrid`, but implementations of `AbstractLandGrid` should allow for separate vertical discretizations in three domains: Ground, Snow, and Canopy. The existing column grids should still be based on `RectilinearGrid`, but there will need to be a new more generic `LandGrid` that wraps an underlying Oceananigans grid and creates three instances, one for each domain. Review the instructions in AGENTS.md for drafting plans and stop to ask any clarifying questions that are necessary.
>
> **Clarifications**:
> - Domain boundaries: **Fixed (pre-computed offsets)** for GPU compatibility
> - Horizontal grid: `LandGrid` wraps **any Oceananigans grid** (not just RectilinearGrid); `ColumnGrid` becomes `ColumnLandGrid` as a sibling implementation
> - Backward compatibility: **None** — breaking changes are acceptable
> - Priority: **Ground domain first**, then Snow and Canopy incrementally
> - GPU testing: Available for validation
> - Timeline: Paper submission in 1-2 months

## Revision log

> 2026-08-25: Initial draft created with clarifications from user feedback. Key decisions: fixed domain boundaries, no backward compatibility, Ground-first priority, LandGrid as generic Oceananigans grid wrapper.

## Problem description

Currently, Terrarium's `AbstractLandGrid` interface is loosely coupled to Oceananigans:
- `AbstractLandGrid{NF, Arch}` is a standalone abstract type (not a subtype of `Oceananigans.AbstractGrid`)
- Current implementations (`ColumnGrid`, `ColumnRingGrid`) wrap a single `Oceananigans.RectilinearGrid`
- The TODO comment in `grids.jl` explicitly acknowledges this is a prototype: "These grid types should be replaced with proper implementations of Oceananigans `AbstractGrid` at some point"
- All three domains (Ground, Snow, Canopy) currently share the same vertical discretization via a single underlying grid

The desired architecture requires:
1. `AbstractLandGrid` to be a proper subtype of `Oceananigans.AbstractGrid`
2. Support for **three separate vertical discretizations** (Ground, Snow, Canopy) within a single land grid
3. **No backward compatibility** — breaking changes are acceptable for this refactoring
4. Seamless integration with Oceananigans' field operations, kernel launching, and node/spacing APIs
5. **Fixed domain boundaries** (pre-computed offsets) for GPU/Reactant compatibility
6. `LandGrid` as a generic wrapper for any Oceananigans grid type; `ColumnLandGrid` as sibling to existing column patterns

## Background

### Oceananigans AbstractGrid interface

Key characteristics of Oceananigans grids:
- **Type parameters**: `AbstractGrid{FT, TX, TY, TZ, Arch}` where topology `{TX, TY, TZ}` and architecture `Arch`
- **Required fields/methods**: `Nx, Ny, Nz, Hx, Hy, Hz, architecture`, plus topology-dependent coordinate arrays
- **Core API**: `size()`, `topology()`, `nodes()`, `halo_size()`, `architecture()`
- **Location system**: Fields defined at `Center`, `Face`, or `Nothing` in each dimension
- **Underlying grid concept**: `AbstractUnderlyingGrid` for primary grids, with curvilinear extensions

### Current Terrarium AbstractLandGrid

```julia
abstract type AbstractLandGrid{NF, Arch} end

# Current wrapper pattern
get_field_grid(grid::AbstractLandGrid) :: Oceananigans.AbstractGrid
Base.size(grid::AbstractLandGrid) = size(get_field_grid(grid))
Architectures.architecture(grid::AbstractLandGrid) = architecture(get_field_grid(grid))
```

Current implementations:
- `ColumnGrid`: Wraps a single `RectilinearGrid` with 1D horizontal (column index) + vertical
- `ColumnRingGrid`: Wraps a `RectilinearGrid` with lateral discretization from `RingGrids.AbstractGrid`

### Multi-domain requirements

Land models need three vertically-stacked domains:
1. **Ground**: Soil profile, typically 50-100+ layers, exponentially spaced
2. **Snow**: Seasonal snowpack, 1-20+ layers, variable thickness
3. **Canopy**: Vegetation layers, 1-10 layers, typically near-surface

Each domain has:
- Independent vertical discretization (different `Nz`, different spacing)
- Shared horizontal discretization (same columns/ring grid)
- Coupled boundary conditions at interfaces (ground-snow, snow-canopy, canopy-atmosphere)

## Summary of changes

### 1. Redefine AbstractLandGrid hierarchy

```julia
# Make AbstractLandGrid a subtype of AbstractGrid with topology parameters
abstract type AbstractLandGrid{NF, TX, TY, Arch} <: Oceananigans.AbstractGrid{NF, TX, TY, Flat, Arch} end
```

### 2. Create LandGrid as generic Oceananigans grid wrapper

`LandGrid` wraps **any** Oceananigans grid (not just `RectilinearGrid`), enabling multi-domain vertical discretizations:

```julia
"""
    LandGrid{NF, TX, TY, Arch, UnderlyingGrid} <: AbstractLandGrid{NF, TX, TY, Arch}

Wraps any Oceananigans grid and creates three domain-specific vertical discretizations
for Ground, Snow, and Canopy with independent vertical layering.

Domain boundaries are **fixed** (pre-computed offsets) for GPU/Reactant compatibility.
"""
struct LandGrid{NF, TX, TY, Arch, UnderlyingGrid<:Oceananigans.AbstractGrid} <: AbstractLandGrid{NF, TX, TY, Arch}
    "Underlying Oceananigans grid for horizontal discretization"
    horizontal_grid::UnderlyingGrid
    
    "Ground domain vertical discretization"
    ground_vert::AbstractVerticalSpacing{NF}
    
    "Snow domain vertical discretization (may be empty/nothing)"
    snow_vert::Union{Nothing, AbstractVerticalSpacing{NF}}
    
    "Canopy domain vertical discretization (may be empty/nothing)"
    canopy_vert::Union{Nothing, AbstractVerticalSpacing{NF}}
    
    "Pre-computed domain offsets for indexing (fixed boundaries)"
    ground_offset::Int
    snow_offset::Int  
    canopy_offset::Int
    
    "Cached domain-specific grids (lazy initialization)"
    ground_grid::Union{Nothing, DomainGrid}
    snow_grid::Union{Nothing, DomainGrid}
    canopy_grid::Union{Nothing, DomainGrid}
end
```

### 3. Define DomainGrid type for domain-specific views

```julia
"""
    DomainGrid{NF, TX, TY, Arch} <: AbstractLandGrid{NF, TX, TY, Arch}

Represents a single domain (Ground/Snow/Canopy) within a LandGrid,
providing a view with domain-specific vertical extent and discretization.
"""
struct DomainGrid{NF, TX, TY, Arch, ParentGrid<:LandGrid} <: AbstractLandGrid{NF, TX, TY, Arch}
    parent::ParentGrid
    domain::Symbol  # :ground, :snow, :canopy
    vertical_start::Int
    vertical_end::Int
    vertical_spacing::AbstractVerticalSpacing{NF}
end
```

### 4. Transform ColumnGrid to ColumnLandGrid

`ColumnLandGrid` is a **sibling** of `LandGrid`, not a wrapper around it:

```julia
"""
    ColumnLandGrid{NF, Arch} <: AbstractLandGrid{NF, Periodic, Flat, Arch}

Represents a set of laterally independent vertical columns with multi-domain
vertical discretizations. This is the 1D column grid implementation.
"""
struct ColumnLandGrid{NF, Arch, UnderlyingGrid<:Oceananigans.AbstractGrid} <: AbstractLandGrid{NF, Periodic, Flat, Arch}
    "Underlying Oceananigans rectilinear grid"
    grid::UnderlyingGrid
    
    "Vertical discretization (Ground domain)"
    vertical_spacing::AbstractVerticalSpacing{NF}
    
    "Domain this grid represents (default: :ground)"
    domain::Symbol
end
```

Existing `ColumnGrid` and `ColumnRingGrid` will be replaced by `ColumnLandGrid`.
No backward compatibility layer — all user code must update to new API.

### 5. Implement AbstractGrid interface for AbstractLandGrid

Required methods (type-stable, allocation-free where possible):
- `Base.size(grid::AbstractLandGrid)` → returns total `(Nx, Ny, Nz_total)`
- `Base.size(grid::AbstractLandGrid, dim)` → size in dimension `dim`
- `topology(grid::AbstractLandGrid)` → `(TX, TY, Flat)`
- `architecture(grid::AbstractLandGrid)` → `Arch`
- `halo_size(grid::AbstractLandGrid)` → `(Hx, Hy, Hz_total)`
- `nodes(grid::AbstractLandGrid, location)` → coordinate arrays
- `eltype(grid::AbstractLandGrid)` → `NF`
- `isrectilinear(grid::AbstractLandGrid)` → `isrectilinear(grid.horizontal_grid)`

Domain-specific dispatch (fixed offsets, pre-computed):
- `size(grid::LandGrid, domain::Symbol)` → domain-specific `(Nx, Ny, Nz_domain)`
- `nodes(grid::LandGrid, domain::Symbol, location)` → domain-specific coordinates
- `get_domain_grid(grid::LandGrid, ::Val{:ground})` → cached `DomainGrid`

### 6. Update Field construction and kernel launching

Modify `grids/grid_utils.jl`:
- `Field(grid::AbstractLandGrid, dims, ...)` → default to Ground domain
- `Field(grid::LandGrid, domain::Symbol, dims, ...)` → domain-specific field
- `launch!(grid::AbstractLandGrid, workspec, ...)` → launch on full grid
- `launch!(grid::LandGrid, domain::Symbol, workspec, ...)` → launch on domain slice

### 7. Domain-aware state variables

- State variables tagged with domain: `PrognosticVariable(:soil_temperature, domain=:ground, ...)`
- Initializers dispatch on domain: `initialize!(state, grid, ::Val{:ground})`
- Default domain is `:ground` for single-domain grids

## Testing and verification

### Unit tests

1. **Grid construction tests** (`test/grids/land_grid_construction.jl`):
   - Create `LandGrid` with three independent vertical discretizations
   - Verify `size()`, `topology()`, `architecture()` methods
   - Test domain-specific size queries
   - Verify offset calculations

2. **AbstractGrid interface tests** (`test/grids/abstract_grid_interface.jl`):
   - Verify `AbstractLandGrid` is subtype of `AbstractGrid`
   - Test `nodes()` for each domain
   - Test `halo_size()` computation
   - Verify `eltype()` propagation

3. **Field construction tests** (`test/grids/field_construction.jl`):
   - Create fields on `LandGrid` (default domain)
   - Create fields on specific domains: `Field(grid, :snow, ...)`
   - Verify field grid associations
   - Test boundary condition application per domain

4. **Kernel launching tests** (`test/grids/kernel_launching.jl`):
   - Launch kernels on full grid
   - Launch kernels on specific domains
   - Verify index mapping (global vs domain-local)
   - Test GPU architecture compatibility

### 4. Backward compatibility tests

**No backward compatibility layer** — all existing code must be updated to use the new API.
This is a breaking change that simplifies the refactoring:
- Old `ColumnGrid` → replaced by `ColumnLandGrid`
- Old `ColumnRingGrid` → replaced by `ColumnLandGrid` with ring grid support
- All examples and tests must be updated

### Integration tests

1. **Multi-domain simulation** (`test/models/multi_domain_land_model.jl`):
   - Run `LandModel` with all three domains active
   - Verify coupling at domain interfaces
   - Test energy/mass conservation across domains

2. **Snow dynamics test** (`test/models/snow_with_ground.jl`):
   - Ground + Snow domains
   - Snow accumulation/melt
   - Ground heat flux through snow

3. **Canopy-atmosphere coupling** (`test/models/canopy_tests.jl`):
   - Canopy + Ground domains
   - Vegetation transpiration
   - Canopy energy balance

### Differentiability tests

Ensure Reactant/Enzyme compatibility:
- `test/reactant/land_grid_ad.jl`: Differentiate through multi-domain initialization
- Verify no throw paths in domain indexing kernels
- Test gradient computation across domain boundaries

## Documentation changes

### API documentation (`docs/src/grids.md`)

- New section: "Multi-domain Land Grids"
- Explain `LandGrid` constructor and domain concepts
- Code examples:
  ```julia
  # Create multi-domain land grid
  ground_vert = ExponentialSpacing(0.05, 100.0, 50)
  snow_vert = UniformSpacing(0.1, 20)
  canopy_vert = UniformSpacing(1.0, 3)
  
  grid = LandGrid(CPU(), Float32, ground_vert, snow_vert, canopy_vert; num_columns=10)
  
  # Query domain properties
  size(grid, :ground)  # (10, 1, 50)
  size(grid, :snow)    # (10, 1, 20)
  size(grid, :canopy)  # (10, 1, 3)
  
  # Create domain-specific field
  soil_temp = Field(grid, :ground, (Center(), Center(), Center()))
  ```

### Model documentation updates

- Update `LandModel` docstring to explain domain keyword arguments
- Add examples of multi-domain simulations
- Document domain-specific boundary conditions

### Migration guide

**Breaking changes documentation**:
- "Migrating from ColumnGrid to ColumnLandGrid"
- Code transformation examples (old → new)
- API differences and rationale

## Known limitations

1. **Performance overhead**: Domain offset calculations add indirection; may impact GPU performance if not inlined properly
2. **Memory usage**: Cached `DomainGrid` instances increase memory footprint (mitigated by lazy initialization)
3. **Complexity**: Multi-domain indexing is more error-prone than single-domain
4. **Breaking changes**: All existing user code must be updated — no backward compatibility
5. **RingGrid integration**: `ColumnLandGrid` multi-domain support requires careful handling of ring grid masks per domain

## Future work

### Phase 2: Dynamic domain activation
- Snow domain appears/disappears based on snow water equivalent
- Canopy domain seasonal activation
- Dynamic grid resizing (challenging for GPU)

### Phase 3: Lateral heterogeneity
- Different vertical discretizations per column
- Terrain-following coordinates for ground
- Adaptive vertical refinement

### Phase 4: Oceananigans upstream contributions
- Propose `MultiDomainGrid` abstraction to Oceananigans
- Contribute land-specific grid optimizations
- Share RingGrid integration patterns

### Phase 5: 2D/3D horizontal grids
- Extend beyond column grids to full 2D horizontal discretizations
- Integrate with Oceananigans `LatitudeLongitudeGrid`
- Support for unstructured horizontal grids

## Implementation steps (phased approach)

**Timeline: Paper submission in 1-2 months**

### Phase 1: Core infrastructure + Ground domain (2 weeks)
1. Define `LandGrid` and `DomainGrid` types with Ground-only support initially
2. Implement `AbstractGrid` interface methods for `AbstractLandGrid`
3. Add unit tests for grid construction and queries (CPU + GPU)
4. Ensure CPU kernel launching works on Ground domain
5. **Milestone**: Paper can use single-domain (Ground) implementation

### Phase 2: Snow domain integration (1 week)
1. Add Snow domain vertical discretization to `LandGrid`
2. Implement domain-aware field construction for Snow
3. Test snow-ground coupling at interface
4. GPU compatibility testing for multi-domain kernels

### Phase 3: Canopy domain + full integration (1 week)
1. Add Canopy domain vertical discretization
2. Full three-domain simulation tests
3. Performance benchmarking (CPU + GPU)
4. Reactant/Enzyme differentiability tests

### Phase 4: Cleanup and documentation (1 week)
1. Complete API documentation with multi-domain examples
2. Update all existing examples to new API
3. Code review and refactoring
4. Final integration tests before paper submission

## Clarifying questions (resolved)

| Question | Answer |
|----------|--------|
| Domain coupling strategy | **Fixed boundaries** (pre-computed offsets) for GPU compatibility |
| Horizontal discretization | `LandGrid` wraps **any Oceananigans grid**; `ColumnLandGrid` is sibling implementation |
| Field location semantics | TBD — can fields exist at domain interfaces? |
| Backward compatibility scope | **None** — breaking changes acceptable |
| Testing infrastructure | GPU resources available for testing |
| Reactant compatibility | Must avoid throw paths; test AD early |
| Priority domains | **Ground first**, then Snow, then Canopy incrementally |
| Timeline constraints | Paper submission in 1-2 months |

### Remaining questions

1. **Field location semantics**: Can fields exist at domain interfaces (e.g., snow-ground boundary)? Should we support `Face` locations that span domains?
2. **Snow depth variation**: With fixed boundaries, how do we handle seasonal snow accumulation/melt? Pre-allocate max snow layers and use masking?
3. **Horizontal grid types**: Which Oceananigans grid types should `LandGrid` initially support? All of them, or just `RectilinearGrid` + `LatitudeLongitudeGrid`?

## Dependencies and prerequisites

- Oceananigans.jl: Current stable release (verify `AbstractGrid` API stability)
- RingGrids.jl: For `ColumnRingGrid` integration
- KernelAbstractions.jl: For GPU kernel launching
- Documenter.jl: For documentation updates
- TestEnv.jl: For isolated test environments

## Risk assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Breaking existing user code | High | Medium | Acceptable — no backward compat, update all tests/examples in Phase 1 |
| GPU performance degradation | Medium | Medium | Profile early, inline offset calculations, benchmark vs single-domain |
| Reactant compilation failures | Medium | High | Avoid throw paths, test AD early, hoist domain logic out of kernels |
| Complexity overwhelm | High | Medium | Phased approach (Ground→Snow→Canopy), rigorous code review, keep core types simple |
| Oceananigans API changes | Low | Medium | Pin Oceananigans version, track upstream changes |
| Timeline pressure (paper) | Medium | High | Ground-only implementation sufficient for paper; multi-domain is future work |

---

*This plan document should be reviewed and signed off before implementation begins. Revise based on feedback from maintainers and potential users.*
