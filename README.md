# Terra.jl
<strong> 🌲🌡💧 Fast, differentiable, and GPU-aware land modeling across scales with [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl). </strong>
---

<a href="https://www.repostatus.org/#wip"><img src="https://www.repostatus.org/badges/latest/wip.svg" alt="Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public." /></a>
<a href="https://eupl.eu/1.2/en">
    <img alt="EUPLv1.2 license" src="https://img.shields.io/badge/License-EUPLv1.2-blue.svg?style=flat-square">
</a>

Terra.jl is a new and upcoming land model that aims to support hybrid physics- and data-driven land modeling across multiple spatial and temporal scales. We envision Terra to be part of a new generation of land models that combines modularity, interactivity, GPU-compability and differentiability for seamless integration of process-based and data-driven model components in both global and regional scale simulations of the land surface.

Terra is being developed alongside [SpeedyWeather.jl](https://github.com/SpeedyWeather/SpeedyWeather.jl) and [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) as the land component of the proposed DELTA-ESM project which aims to build a new class of Earth System Model in the Julia programming language that enables hybrid geophysical modeling across multiple scales and application domains.

> [!WARNING]
> 🚧🚧 Construction Site! 🚧🚧
>
> This is an early prototype of our model and is not production-ready. Expect things to change rapidly and break often. If you share our vision for a new paradigm of Earth system modeling and would like to get involved in the project, don’t hesitate to reach out by creating an issue on GitHub issues or sending us an email. We are always happy to welcome new collaborators!

## Our goals
We want a land surface (and subsurface) model that is
- Fast enough to run global-scale simulations on a laptop at coarse (100 km) resolutions
- Flexible enough to scale-up to high resolution simulations in HPC environments
- Fully GPU-compatible with the ability to easily switch between running on CPU-based and GPU-based architectures
- Fully auto-differentiable with Enzyme.jl to enable systematic parameter estimation and hybrid modeling with neural differential equations
- Modular and extensible to allow for rapid prototyping of model components with varying levels of complexity
- Interactive and user-friendly to make Earth System Modeling fun and accessible for a larger audience of researchers and practitioners, as well as students and educators
- Open-source and community-driven to foster interdisciplinary collaboration and development

It is important to emphasize, however, what Terra is not:
- Terra is not a comprehensive land model. While we are always open to suggestions and collaborations to implement new processes, we do not aim to build a state-of-the-art terrestrial ecosystem model rivaling that of, e.g. the CTSM. Terra - will always favor simplicity, efficiency, and interoperability over process-complexity.
- Terra is not “just another model”. We do not intend for users to simply download our model products and cite our papers. We want users to directly interact with our model, ideally running their own simulations and writing their own code.
- Terra is not a monolithic model. Modularity and extensibility are core to our vision. Terra provides a library of models, process implementations, and numerical tools which users can use to build their own simulations. We will provide guidance and a set of well-tested and stable model configurations, but we encourage users to experiment and push the limits of what those models can do.

## Why Oceananigans?
It might initially seem strange that a land model would be built on top of a framework for ocean modeling. There are some key advantages in doing so:


1. Firstly, like ocean models, land models are commonly implemented using the finite volume method (FVM) to approximate spatial gradients in mass and energy conservation laws. Oceananigans provides state-of-the-art tools for FVM simulation in Julia, with a focus on geophysical applications, which aligns well with our goals. Although Terra will primarily focus on 1D (vertical) dynamics in the beginning, using Oceananigans affords us the possibility of very feasibly expanding to 2D and 3D simulations of surface and subsurface heat and water transport in the future!
2. Secondly, the numerical operators provided by Oceananigans are built to be both auto-differentiable and GPU-compatible out-of-the-box, which means that Terra can inherit these capabilities almost “for free”.
3. Finally, and perhaps most importantly, we believe in the vision pioneered by the Climate Modeling Alliance and the Oceananigans team for the development of a new generation of Earth System Models that are open, accessible, interactive, and capable of learning from data in ways that go beyond traditional data assimilation.

## Installation

Terra is still in prototype stage and is not yet registered as a package in the Julia General registry.

However, you can still install the package from the repository via the package manager (type `]` in your REPL):

```
pkg> add https://github.com/TUM-PIK-ESM/Terra.jl
```

or clone the repository and start hacking directly!

## Quick start

TODO

## Contributing

An open source project is only as strong as its community of contributors. We're alays happy to accept contributions, no matter how big or small!

Terra.jl is in a very early stage of development, so this is a golden opportunity for you to get your ideas in on the ground floor. If you have some ideas or code you would like to contribute, please don't hesitate to create an issue and get involved!


## Copyright and license

See the [NOTICE](./NOTICE) file for copyright information.

Terra.jl is free and open source licensed under the [European Union Public Lcense v1.2](https://eupl.eu/1.2/en).

What does that mean for you? You are 100% free to
- Copy, modify, and redistribute the code
- Use the software as a package in your own Julia project (regardless of license or copyright status)
- Use the software for both commercial and non-commerical purposes

However, if you make changes or modification to the code (excluding use as a package/library or changes for the purpose of interoperability), **you are required to re-distribute the modified software under the EUPL v1.2 or a compatible license**.
