# Contributing

We gladly welcome any and all contributions to Terrarium.jl. Building a new land model is a huge undertaking that no single developer or scientist can hope to achieve alone. Collaboration always has been and always will be key to building good geoscientific models. Terrarium is no exception to this, and its success will depend on the contributions of the broader community.

There are multiple ways in which you could consider contributing to the project:
- If you have a question or an idea, raise an issue or start a discussion on the GitHub repository.
- If you want to try implementing something, clone the repository and make a pull request. Please take note of our [Software development practices](@ref) below.
- If you're not able to directly contribute yet but would like to support our work, consider sharing the GitHub repository (or this documentation) with others who you think might be interested.

Regardless of how you choose to contribute, we thank you for your participation, and we look forward to working with you!

## Software development practices

### Automated testing and continuous integration

Terrarium.jl adheres to software development standards for automated testing via continuous integration. We write unit tests for every function of our model. In some cases this might appear trivial, but we still want to achieve a near complete coverage of our code in the tests. The majority of the tests should cover the smallest possible units over different input arguments and types (if applicable). Unit tests should typically call the tested functions in a way that is representative for their use in the model, but try to reduce the computational complexity (e.g. by choosing very low dimensional inputs) to keep the overall CI time manageable. Additionally, we have some tests that ensure top-level functionality and stability of the model as well. Every additional proposed feature in a Pull Request has to come with unit tests. Tests verifying differentiability and GPU compatibility will also be required.

### Automatic Differentiation with Enzyme 

For AD, we rely primarily on reverse mode differentiation via [Enzyme.jl](https://enzyme.mit.edu/julia/stable/). In contrast to many other AD systems, Enzyme doesn’t put particularly strong restrictions on coding style. For example, Array mutations are not only allowed, they are even encouraged!

There is one thing however, that is crucial for Enzyme to work: [type stability](https://docs.julialang.org/en/v1/manual/faq/#man-type-stability). The code absolutely has to be type stable, even in parts that are performance non-critical. Every occurrence of type instability may break differentiability with Enzyme. If this happens, we recommend first quickly checking for type instability of core function calls using `@code_warntype` and `@inferred` in unit tests. Then, if the issue is not yet apparent, debug the code using [Cthulhu.jl](https://github.com/JuliaDebug/Cthulhu.jl). Cthulhu is an awesome tool to interactively inspect your code for type instabilities and it is extremely easy to run: just preempt the function that you are trying to differentiate with `@descend` after loading the `Cthulhu` package.

Enzyme does, however, have some disadvantages; it is still not fully mature and bugs do occur. As of the time of writing (August 2025), this is especially the case for Julia 1.11. We currently recommend staying on Julia 1.10.10 (LTS) for the time being. Other cryptic Enzyme error messages have become rarer with time, but they do still occasionally happen. In these cases, we, along with our AD team led by Valentin in the DELTA-ESM project, are happy to offer support to the best of our abilities.

### Code formatting

We use [Runic.jl](https://github.com/fredrikekre/Runic.jl) for automated code formatting. If you submit a PR you probably have seen a comment from our CI that the code is not formatted according to the Runic style. To use runic, you have to run the install Runic install script and then you can choose to either format directly from the command line with `runic --inplace .` (don't forget the `.` at the end!), or configure it in your editor of choice or use a Git hook. All details are described in [Runic's readme](https://github.com/fredrikekre/Runic.jl). We recommend setting up Runic as a Git hook to automatically format your code on commit. For that purpose, we provide the Git hook in the `.githooks` directory of the repository. You can copy the hook from there to your local git hooks directory and make it executable to use like so: 
```
cp .githooks/pre-commit .git/hooks/pre-commit
chmod +x .git/hooks/pre-commit
```

## Documentation

All pull requests that implement new features or modify existing functionality must have associated documentation. At bare minimum, this should consist of docstrings on all of the relevant functions and types. However, in many cases, it can be helpful to add a documentation page or example script that showcases the feature(s) and helps the user understand holistically how it fits into the Terrarium framework.

The docs can be built locally by running

```
julia --project=docs docs/make.jl --local
```

To skip running doctests and example scripts, you can also add `--draft` or `-d` for short.

Preview builds of the documentation associated with pull requests can be reviewed at

https://numericalearth.github.io/Terrarium.jl/previews/PR##

replacing `##` with the pull request ID number.
