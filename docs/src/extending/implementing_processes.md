# Implementing processes

```@meta
CurrentModule = Terrarium
```

This page is a step-by-step guide for implementing a new concrete process type in Terrarium.
We assume familiarity with the [Basic concepts](@ref) and [Core interfaces](@ref) pages, both of which cover the `AbstractProcess` interface and variable system.
For a concrete, self-contained worked example that follows the workflow below step by step, see the [Implementing a process: 1D linear heat diffusion](@ref) example.

## Control flow

Every process in Terrarium is implemented across three levels of abstraction:
I. **Top-level interface methods** for `AbstractProcess`, i.e. `variables`, `initialize!`, `compute_auxiliary!`, and `compute_tendencies!` (see doc section on the [`AbstractProcess` interface](@ref "The AbstractProcess interface"))
II. **Kernel functions**; this includes both the `@kernel` entry point and inlined functions with signatures of the form `compute_*(i,j[,k], grid, fields, ::Process, args...)`
III. **Process methods** defined directly by the subtype of `AbstractProcess`. These are typically scalar-valued functions corresponding to individual terms or expressions in the mathematical formulation of the process physics.

The flow of execution is I → II  → III: top-level methods are invoked by the enclosing `AbstractModel`, these methods then `launch!` their corresponding `@kernel`s which in turn call the inner kernel functions.

## Development workflow

We recommend to start implementing new processes using the following workflow:
1. Start by drafting an initial version of the process `struct` subtyping `AbstractProcess`; it's good to start simple by directly putting all relevant parameters into this `struct`. You can later consider grouping these parameters into logical parameterization types.
2. Define `variables(::Process)` for the new process type. Identify which variables are `prognostic`, which are `auxiliary`, and which might be assumed to be `input`s, provided either by external forcing data or by another process.
3. Write down unit tests and initial implementations of the the scalar-valued **process methods** (**level III**). These should correspond to the key equations or terms in the mathematical formulation of the physics. Note that these methods are invoked on **each elementary volume** and thus cannot be used to define equations involving spatial gradients or fluxes.
4. Implement the `compute_*(i, j[, k], grid, fields, ...)` kernel functions with the corresponding `@kernel` entry points `compute_auxiliary_kernel!` and `compute_tendencies_kernel!` (**level II**). Functions which might be independently useful outside of the process interface (e.g. [`adjust_saturation_profile!`](@ref) in `SoilHydrology`) can also define their own `@kernel` entry points.
5. Finally, write implement the top-level process methods for `initialize!`, `compute_auxiliary!` and `compute_tendencies!` that `launch!` their respective `@kernel`s (**level I**).

See the example [Linear heat conduction](@ref) for a detailed walkthrough of this workflow.
