---
title: 'Terrarium.jl: Flexible, GPU-accelerated, and fully auto-differentiable land modeling in Julia'
tags:
  - Julia
  - land
  - land surface
  - climate
  - numerical
  - differentiable
  - GPU-accelerated
authors:
  - name: Brian Robert Groenke
    orcid: 0000-0003-2570-9342
    corresponding: true
    equal-contrib: true
    affiliation: 1
  - name: Maximilian Gelbrecht
    orcid: 0000-0002-0729-6671
    equal-contrib: true
    affiliation: "1, 2"
  - name: Maha Badri
    affiliation: "1, 2"
  - name: Olivier Bonte
    orcid: 0000-0003-1806-7572
    affiliation: 3
affiliations:
 - name: Potsdam Institute for Climate Impact Research (PIK), Germany
   index: 1
   ror: "03e8s1d88"
 - name: Technical University of Munich, Germany
   index: 2
   ror: "02kkvpp62"
 - name: Ghent University, Belgium 
   index: 3
   ror: "00cv9y106"
date: July 2026
bibliography: paper.bib
---

# Summary

The land surface is a crucial component of Earth's climate system. Terrestrial ecosystems, as well as the solid earth that hosts them, exchange heat, water, and carbon with the atmosphere on timescales ranging from seconds to millenia, driving key feedbacks that significantly affect both local and global climate conditions. Land surface and terrestrial ecosystem models aim to simulate the physical processes governing the land's response to changes in both weather and climate. They are therefore crticially important components of so-called *Earth System Models* (ESMs) which aim to simulate the natural variability of the Earth system, as well as its long-term response to changes in both natural and anthropogenic forcing.

Terrarium.jl is a software package in the Julia programming language [@bezansonJuliaFreshApproach2017] that aims to enable fast, flexible, and user-friendly land surface and terrestrial ecosystem modeling both standalone and coupled to other Earth system component models. It enables the user to quickly and flexibly configure and run simulations of terrestrial heat, water, and carbon dynamics at both global and regional scales. The numerical implementation is based on Oceananigans.jl [@ramadhanOceananigansjlFastFriendly2020; @wagnerHighlevelHighresolutionOcean2025], a widely used software package for "ocean-flavored" fluid dynamics that provides a rich set of numerical primitives for GPU-accelerated and automatic differentiation (AD) compatible finite volume computations on both rectangular and spherical grids. Terrarium.jl uses these primitives to solve the relevant governing equations for nonlinear transport in porous media (i.e. soil and rock), as well as to facilitate the efficient computation of point-scale processes like the vegetation carbon cycle, surface hydrology, and the exchange of water and energy between the land surface and the atmosphere.

Terrarium.jl is being developed as part of the NumericalEarth project, which is an international consortium of researchers and software developers dedicated to building modern and accessible tools for multi-scale, multi-architecture, and fully differentiable Earth system modeling, centered around Oceananigans as a numerical core.

# Statement of need

Most existing ESM land components, such as those involved in the international Coupled Model Intercomparison Project (CMIP), are implemented in languages such as Fortran and C++; these languages, despite their advantages, incur a significant engineering burden and limit the accessibility of the model code for many scientists and researchers. Furthermore, Fortran/C++ make it far more difficult for these models to take advantage of state-of-the-art tools for high-level GPU programming and AD which have served as the foundation of recent advances in machine learning and data-driven modeling [@gelbrechtDifferentiableProgrammingEarth2023]. As such, there is a general need for a new generation of ESM components that are built using modern numerical computing tools in higher level languages such as Julia or Python, which can greatly improve accessibility and productivity.

# State of the field

Terrarium.jl represents a significant step towards filling this need, though it should be noted that it is not the first. There have been recent efforts towards implementing modern land surface models in Python/JAX, such as DifferLand [@fangDifferentiableLandModel2026] and JAX-CanVeg [@jiangJAXCanVegDifferentiableLand2025], both of which enable GPU-accelerated and AD-compatible land modeling in Python. Similarly, recent work by the Climate Modeling Alliance (CliMA) on ClimaLand.jl [@deckClimaLandLandSurface2026] has demonstrated the promise of GPU-accelerated land modeling in Julia as part of the Clima projects efforts to build a new, CMIP-class ESM capable of being automatically calibrated using global observational datasets. In contrast, Terrarium.jl enables not only GPU-accelerated, global land simulations but also full comaptibility with both forward- and reverse-mode AD via Enzyme.jl and Reactant.jl [@mosesInsteadRewritingForeign2020; mosesDJ4EarthDifferentiablePerformancePortable2026]. In addition, Terrarium.jl is designed to allow for fast and efficient coupling with other global or regional, GPU-accelerated ESM components such as the SpeedyWeather.jl [@klowerSpeedyWeatherjlReinventingAtmospheric2024] and Breeze.jl (REF) atmosphere models.

# Software design

Terrarium.jl is designed to be highly modular, providing a framework for constructing a multitude of different land model configurations rather than implementing a single, monolithic model design. This allows Terrarium to serve as a common set of numerical tools and physical process implementations which can serve as the basis for a wide range of land and ecosystem models, spanning both highly detailed local-scale simulations with prescribed boundary conditions and intermediate-complexity global-scale simulations as part of a coupled Earth system model. Terrarium models are composed of four basic components: one or more `Process`es (subtyping `AbstractProcess`) implementing specific physical relationship and governing equations, a `grid` that describes the discretization of the underlying spatial domain, a `timestepper` that determines how the prognosic variables of the model are advanced in time, and an `initializer` that specifies how the initial state of the prognostic variables is determined. Each `Process` is required to implement the following interface:
- A method `variables(::Process)` that returns a tuple of state variables required by the `Process`. Each variable must be declared as `prognostic`, `auxiliary`, or `input`,
- A method `compute_auxiliary!(state, grid, ::Process, args...)` which computes all `auxiliary` variables based on the current prognostic state,
- A method `compute_tendencies!(state, grid, ::Process, args...)` which computes the tendencies of the `prognostic` variables at the current time step.
This interface helps to enforce a mathematically consistent model design whereby the state is fully specified by the prognostic variables and all other "auxiliary" variables are derived from them. Each `Process` type is free to specify additional arguments (here represented by `args...`) corresponding to other `Process` or helper types on which they are dependent. Terrarium therefore relies heavily on Julia's system for multiple dispatch to choose the correct method implementation at runtime based on the types of the given arguments.

While the modular design of Terrarium.jl is highly convenient for mixing and matching various process implementations and parameterizations, it also comes with some important trade-offs. Terrarium also prioritizes support for multi-architecture parallelization via KernelAbstractions.jl and automatic differentiation via Enzyme.jl and Reactant.jl. Maximizing modularity would necessitate each `Process` launching its own GPU kernel, which can incur significant overhead in models that involve coupling together several different highly interdependent physical processes. To balance modularity with efficient parallelization, Terrarium `Process`es are typically implemented at three different levels of abstraction: a top-level interface (`compute_auxiliary!` and `compute_tendencies!`) which can be used to evaluate and test each `Process` type standalone, a grid-aware kernel function interface (e.g. `compute_energy_tendency(i, j, k, grid, fields, ::Process, args...)`) called at the given grid indices `i, j, k` within a device-agnostic compute kernel, and primitive functions (e.g. `physical_quantity(::Process, args...)`) that implement scalar equations independent of any particular spatial grid. This pattern allows for most of the actual physics code (defined at the kernel level) to be reused across different kernels. Terrarium makes use of this in so-called `CoupledProcess` types which define fused kernels that bundle multiple `Process`es into a single component. This allows develoeprs to balance modularity with GPU efficiency at the cost of some additional boilerplate incurred by the need to split functions across the different abstraction levels.

# Research impact statement

WIP

Terrarium.jl's modular design, GPU-compatibility, and first-class support for AD make it an ideal testbed for cutting edge research on so-called "hybrid" land-atmosphere modeling, which aims to combine data-driven components, such as neural networks, with physics-based simulation.


Terrarium simulations can be plugged into a NumericalEarth `EarthSystemModel`, which provides a unified coupling interface for atmosphere, ocean, land, and sea ice models. In addition, Terrarium directly couples with the SpeedyWeather.jl atmosphere model [@klowerSpeedyWeatherjlReinventingAtmospheric2024], allowing for coupled land-atmosphere simulations using SpeedyWeather's native interface. 

We believe that such a tool will empower geoscientific researchers to move beyond traditional workflows that rely on physical simulations only to generate training and validation data for machine learning models, towards a new generation of scientific models that are capable of learning from data while still respecting physical laws.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# AI usage disclosure

No generative AI was used in the preparation of this manuscript. The software itself is developed with limited AI assistance and strict human oversight. All pull requests are reviewed and signed off by human developers.


# Acknowledgements



# References

