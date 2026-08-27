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

Terrarium.jl is a software package for fast, flexible, and user-friendly simulation of terrestrial heat, water, and carbon dynamics at both global and regional scales. It is based on the "ocean-flavored" fluid dynamics package, Oceananigans.jl [@ramadhanOceananigansjlFastFriendly2020; @wagnerHighlevelHighresolutionOcean2025], which provides a rich set of numerical primitives for GPU-accelerated finite volume computations on both rectangular and spherical grids. It is being developed as part of the NumericalEarth project, which is an international consortium of researchers and software developers dedicated to building modern and accessible tools for multi-scale Earth system modeling, centered around the Oceananigans numerical core. The numerical primitives provided by Oceananigans are used to solve the relevant governing equations for nonlinear transport in porous media (i.e. soil and rock), as well as to facilitate the efficient computation of point-scale processes like the vegetation carbon cycle, surface hydrology, and the exchange of water and energy between the land surface and the atmosphere.

Terrarium.jl is built from the ground up to support multi-architecture parallelization via KernelAbstractions.jl (REF?) and automatic differentiation via Enzyme.jl and Reactant.jl [@mosesInsteadRewritingForeign2020; mosesDJ4EarthDifferentiablePerformancePortable2026]. This makes Terrarium ideally suited to support the implementation of so-called "hybrid" land-atmosphere models that incorporate data-driven components, such as neural networks, into physics-based simulations. In addition, Terrarium is also designed to be highly modular, providing a framework for constructing a multitude of different land model configurations rather than implementing a single, monolithic model design. This allows Terrarium to serve as a common set of numerical tools and physical process implementations which can serve as the basis for a wide range of land and ecosystem models, spanning both highly detailed local-scale simulations with prescribed boundary conditions and intermediate-complexity global-scale simulations as part of a coupled Earth system model.


# Statement of need

Land surface models are key components of Earth System Models (ESMs) which aim to simulate coupled atmosphere, ocean, and land dynamics at global scales. However, most traditional land surface models are implemented in languages such as Fortran and C++, which pose a high barrier of entry for less experienced programmers.

# State of the field


# Software design


Terrarium simulations can be plugged into a NumericalEarth `EarthSystemModel`, which provides a unified coupling interface for atmosphere, ocean, land, and sea ice models. In addition, Terrarium directly couples with the SpeedyWeather.jl atmosphere model [@klowerSpeedyWeatherjlReinventingAtmospheric2024], allowing for coupled land-atmosphere simulations using SpeedyWeather's native interface. 

# Research impact statement


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

