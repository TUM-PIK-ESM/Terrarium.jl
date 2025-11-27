using ArgParse
using Documenter
using Literate
using PlutoStaticHTML

using Terrarium

const NOTEBOOK_DIR = joinpath(dirname(@__DIR__), "examples", "notebooks")
const EXAMPLE_DIR = joinpath(@__DIR__, "src", "notebooks")
const EXAMPLE_DIR_RELATIVE = joinpath("notebooks")

# lookup table for all Pluto notebooks to be included 
notebook_lookup = Dict(
    "Model Interface" => "example_model_notebook.md",
    "Differentiating Terrarium" => "differentiate-notebook.md",
)

"""
    build()

Run all Pluto notebooks (".jl" files) in `NOTEBOOK_DIR`.
"""
function build()
    println("Building notebooks in $NOTEBOOK_DIR and moving them to $EXAMPLE_DIR")
    oopts = OutputOptions(; append_build_context=false)
    output_format = documenter_output
    bopts = BuildOptions(NOTEBOOK_DIR; output_format)
    build_notebooks(bopts, oopts)

    # move to docs/src/notebooks because for some reason that's needed
    mkpath(EXAMPLE_DIR)
    for (_, file_name) in notebook_lookup
        mv(joinpath(NOTEBOOK_DIR, file_name), joinpath(EXAMPLE_DIR, file_name))
    end

    return nothing
end

# Build the notebooks; defaults to true.
if get(ENV, "BUILD_DOCS_NOTEBOOKS", "true") == "true"
    build()
end

# Dict for makedocs for notebooks to be included 
notebook_docpages = Pair{String, String}[]
for (title, name) in notebook_lookup
    push!(notebook_docpages, title => joinpath(EXAMPLE_DIR_RELATIVE, name))
end

s = ArgParseSettings()
@add_arg_table! s begin
    "--local", "-l"
       action = :store_true
       help = "Local docs build mode"
    "--draft", "-d"
       action = :store_true
       help = "Whether to build docs in draft mode, i.e. skipping execution of examples and doctests"
end
parsed_args = parse_args(ARGS, s)

IS_LOCAL = parsed_args["local"] || parse(Bool, get(ENV, "LOCALDOCS", "false"))
IS_DRAFT = parsed_args["draft"] || parse(Bool, get(ENV, "DRAFTDOCS", "false"))
if haskey(ENV, "GITHUB_ACTIONS")
    ENV["JULIA_DEBUG"] = "Documenter"
end

makedocs(
    format = Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing)=="true",
        ansicolor=true,
        collapselevel=1,
        canonical = "https://tum-pik-esm.github.io/Terrarium.jl/stable/",
        size_threshold = 600_000,
        # Using MathJax3 since Pluto uses that engine too.
        mathengine=Documenter.MathJax3(),
    ),      # in bytes
    sitename = "Terrarium.jl",
    authors = "Brian Groenke, Maximillian Galbrecht, Maha Badri, and Contributors",
    modules = [Terrarium],
    pages = [
        "Home" => "index.md",
        "Overview" => [
            "Mathematical formulation" => "overview/mathematical_formulation.md",
            "Software architecture" => "overview/software_architecture.md",
        ],
        "Physics" => [
            "Soil physics" => [
                "Energy and water balance" => "physics/soil_energy_water.md",
            ],
            "Vegetation" => "physics/vegetation.md",
        ],
        "Examples" => notebook_docpages,
        "Contributing" => "contributing.md",
        "API Reference" => "api_reference.md",
    ],
    draft=IS_DRAFT,
)

deployconfig = Documenter.auto_detect_deploy_system()

# remove gitignore from build files
# rm(joinpath(@__DIR__, "build", ".gitignore"))

deploydocs(
       repo="github.com/TUM-PIK-ESM/Terrarium.jl.git",
       push_preview = true,
       versions = ["v0" => "v^", "v#.#", "dev" => "dev"],
       deploy_config = deployconfig,
)
