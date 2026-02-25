# Terrarium.jl
<strong> 🌲🌡💧 Fast, differentiable, and GPU-aware land modeling across scales with [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl). </strong>
---

<a href="https://tum-pik-esm.github.io/Terrarium.jl/dev">
<img alt="Development documentation" src="https://img.shields.io/badge/documentation-in%20development-orange?style=flat-square">
</a>
<a href="https://www.repostatus.org/#wip"><img src="https://www.repostatus.org/badges/latest/wip.svg" alt="Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public." /></a>
<a href="https://eupl.eu/1.2/en">
    <img alt="EUPLv1.2 license" src="https://img.shields.io/badge/License-EUPLv1.2-blue.svg?style=flat-square">
</a>
<a href="https://github.com/fredrikekre/Runic.jl">
    <img alt="code style: runic" src="https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black.svg?style=flat-square">
</a>

[Terrarium.jl](https://tum-pik-esm.github.io/Terrarium.jl/dev) is a new and upcoming framework for hybrid physics- and data-driven land modeling across spatial and temporal scales. We envision Terrarium to be part of a new generation of Earth system component models that combine modularity, interactivity, GPU-compability and auto-differentiability (AD) for seamless integration of process-based and data-driven model components in both global and regional scale simulations.

Terrarium is being developed alongside [SpeedyWeather.jl](https://github.com/SpeedyWeather/SpeedyWeather.jl) and [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) as the land component of a new, user-friendly, and fully GPU/AD-compatible Earth System Model in the Julia programming language.

> [!WARNING]
> 🚧🚧 Under Construction! 🚧🚧
>
> Terrarium.jl is still in a prototyping stage and is not production-ready. Expect things to change rapidly and break often. If, however, you share our vision for a new paradigm of Earth system modeling and would like to get involved in the project, please don’t hesitate to reach out by creating an issue on GitHub issues or sending us an email. We are always happy to welcome new collaborators!

## Our goals
We want a land surface (and subsurface) model that is
- **Fast** enough to run global-scale simulations on a laptop at coarse (100 km) resolutions
- **Flexible** enough to scale-up to high resolution simulations in HPC environments
- **Fully GPU-compatible** with the ability to easily switch between running on CPU-based and GPU-based architectures
- **Fully auto-differentiable** with Enzyme.jl to enable systematic parameter estimation and hybrid modeling with neural differential equations
- **Modular** and **extensible** to allow for rapid prototyping of model components with varying levels of complexity
- **Interactive** and **user-friendly** to make land surface modeling fun and accessible for a larger audience of researchers and practitioners, as well as students and educators
- **Open-source** and **community-driven** to foster interdisciplinary collaboration and development

It is important to emphasize, however, what Terrarium is not:
- Terrarium is **not a comprehensive land model**. While we are always open to suggestions and contributions to add new processes, we do not aim to build a state-of-the-art terrestrial ecosystem model rivaling that of, e.g. the [CTSM](https://github.com/ESCOMP/CTSM). Terrarium will always favor simplicity, efficiency, and interoperability over process-complexity.
- Terrarium is **not “just another model”**. We do not intend for users to simply download our model products and cite our papers. We want users to directly interact with our model, ideally running their own simulations and writing their own code.
- Terrarium is **not a monolithic model**. Modularity and extensibility are core to our vision. Terrarium provides a library of models, process implementations, and numerical tools which users can use to build their own simulations. We will provide guidance and a set of well-tested and stable model configurations, but we encourage users to experiment and push the limits of what those models can do.

## Installation

Terrarium is still in an early prototype stage and is not yet registered as a package in the Julia General registry.

However, you can still install the package from the repository via the package manager (type `]` in your REPL):

```
pkg> add https://github.com/TUM-PIK-ESM/Terrarium.jl
```

or clone the repository and start hacking directly!

## Quick start

A natural first step with `Terrarium` is to set up and run your very first `SoilModel`. This represents a standalone model of transient heat, water, and carbon transport over a particular choice of `grid`. We start by chosing a `ColumnGrid` which represents one or more laterally independent vertical columns:

```julia
using Terrarium

# Set up a SoilModel on a ColumnGrid with 10 vertical soil layers that will run on the CPU with 32-bit precision
num_columns = 1
arch = CPU()
grid = ColumnGrid(arch, Float32, ExponentialSpacing(N=10), num_columns)
model = SoilModel(grid)
# Prescribe a constant surface temperature of 1°C
bcs = PrescribedSurfaceTemperature(:T_ub, 1.0)
integrator = initialize(model, ForwardEuler(eltype(grid)), boundary_conditions = bcs)
# Run the simulation forward for 10 model days
@time run!(integrator, period = Day(10))
```

That's it! You already succesfully ran a (very simple) simulation with Terrarium!

Note that setting `num_columns = 1` here corresponds to a point simulation for a single vertical column. However, we can easily scale this up by set `num_columns` to any positive integer (up to the memory limit of your system, of course).

We can also easily adapt this code to run a *global* simulation over a suitable spatial grid. For this, we'll need to have [`RingGrids`](https://github.com/SpeedyWeather/RingGrids.jl) installed (or import it directly from Terrarium with `using Terrarium.RingGrids`). Optionally, if a GPU is available and [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) is installed in the current project or global Julia environment, we can accelerate the global simulation by simply changing `CPU` to `GPU`:

```julia
using RingGrids: FullGaussianGrid
using CUDA # needs to be seaprately installed

rings = FullGaussianGrid(8) # Gaussian grid with 16 lattitudinal rings (512 points, ~9.5˚)
arch = GPU() # run on the GPU!
grid = ColumnRingGrid(arch, Float32, ExponentialSpacing(N=10), rings)
model = SoilModel(grid)
# Prescribe a constant surface temperature of 1°C
bcs = PrescribedSurfaceTemperature(:T_ub, 1.0)
integrator = initialize(model, ForwardEuler(eltype(grid)), boundary_conditions = bcs)
# Run the simulation forward for 10 model days
@time run!(integrator, period = Day(10))
```
and voila! We have just run a GPU-accelerated, global-scale simulation of soil thermal dynamics with minimal additional effort. While more realistic simulations are of course more involved, this simple example demonstrates the core of what we aim to accomplish with Terrarium; a fast, user-friendly, and highly adaptable land model that can easily be configured to run on local, regional, and global scales.

## Why Oceananigans?
It might initially seem strange that a land model would be built on top of a framework for ocean modeling. There are, however, some key advantages in doing so:


1. Firstly, like ocean models, land models are commonly implemented using finite difference and/or finite volume method (FDM/FVM) to approximate spatial gradients in mass and energy conservation laws. Oceananigans provides state-of-the-art tools for FVM simulation in Julia, with a focus on geophysical applications, which aligns well with our goals. Like most land models, Terrarium will initially focus on 1D column modeling; however, using Oceananigans affords us the possibility of very feasibly expanding to 2D and 3D simulations in the future!
2. Secondly, the numerical operators provided by Oceananigans are built to be both auto-differentiable and GPU-compatible out-of-the-box, which means that Terrarium can inherit these capabilities almost “for free”.
3. Finally, and perhaps most importantly, we believe in the vision pioneered by the [Climate Modeling Alliance](https://clima.caltech.edu/) and the [NumericalEarth](https://github.com/NumericalEarth/) projects for the development of a new generation of Earth System Models that are open, accessible, interactive, and capable of learning from data in a multitude of ways that go beyond traditional data assimilation.

## Contributing

An open source project is only as strong as its community of contributors. We're always happy to accept contributions, no matter how big or small!

Terrarium.jl is in a very early stage of development, so this is a golden opportunity for you to get your ideas in on the ground floor. If you have some ideas or code you would like to contribute, please don't hesitate to create an issue and get involved!

## Copyright and license

See the [NOTICE](./NOTICE) file for copyright information.

Terrarium.jl is free and open source licensed under the [European Union Public License v1.2](https://eupl.eu/1.2/en).

What does that mean for you? You are 100% free to
- Copy, modify, and redistribute the code
- Use the software as a package in your own project (regardless of license or copyright status)
- Use the software for both commercial and non-commerical purposes

However, if you make changes or modification to the code, excluding those made purely for the purpose of interoperability, **you are required to re-distribute the modified software under the EUPL v1.2 or a compatible license**. This is vital to ensure the long-term survival of the project and to foster an open, supportive, and diverse community.

## Acknowledgements

Terrarium.jl is a community-oriented project lead by the [FutureLab on Artificial Intelligence](https://www.pik-potsdam.de/en/institute/departments/complexity-science/research/artificial-intelligence) at the Potsdam Institute for Climate Impact Research (PIK) and Earth System Modeling group at the Technical University of Munich (TUM). We acknowledge funding support from
- Horizon Europe grant agreement number 101184070
- Volkswagen Foundation
- Munich Data Science Institute (MDSI) at Technical University of Munich (TUM) via the Linde/MDSI Doctoral Fellowship program

<img width="409" height="91" alt="image" src="https://github.com/user-attachments/assets/d75524fc-2cd4-4951-91c6-d13f0b62229f" />

**Disclaimer**:
*Views and opinions expressed are however those of the author(s)
only and do not necessarily reflect those of the European Union or
European Research Executive Agency (REA). Neither the
European Union nor the granting authority can be held responsible
for them.*
