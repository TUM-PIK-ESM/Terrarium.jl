using ArgParse
using Documenter
using Literate

using Terrarium

const SCRIPTS_DIR = joinpath(dirname(@__DIR__), "examples", "scripts")
const EXAMPLES_OUTDIR = joinpath(@__DIR__, "src", "examples")
const EXAMPLES_OUTDIR_RELATIVE = "examples"

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
BUILD_EXAMPLE_DOCS = !IS_DRAFT && parse(Bool, get(ENV, "BUILD_EXAMPLE_DOCS", "true"))
if haskey(ENV, "GITHUB_ACTIONS")
    ENV["JULIA_DEBUG"] = "Documenter"
end

# Literate scripts to be built and included in the docs
# Each entry is (page_title, script_filename)
script_list = [
    "Model Interface" => "model_interface.jl",
    #"Differentiating Terrarium" => "differentiating_terrarium.jl",  # dont include enzyme example for now while it's bugged
]

"""
Convert Literate.jl scripts in `SCRIPTS_DIR` to markdown pages in `EXAMPLES_OUTDIR`.
The differentiation example is never executed during docs builds (Enzyme compile is
too slow); the model interface example is executed unless IS_DRAFT is set.
"""
function build_literate_pages()
    mkpath(EXAMPLES_OUTDIR)
    for (_, filename) in script_list
        ## the differentiation notebook is never auto-executed (Enzyme compile time)
        should_execute = BUILD_EXAMPLE_DOCS && filename != "differentiating_terrarium.jl"
        kwargs = Dict{Symbol, Any}(
            :execute => should_execute,
            :documenter => true,
            :flavor => Literate.DocumenterFlavor(),
        )
        ## For non-executed scripts, use plain julia code fences so Documenter
        ## does not attempt to run them as @example blocks.
        if !should_execute
            kwargs[:codefence] = "```julia" => "```"
        end
        Literate.markdown(
            joinpath(SCRIPTS_DIR, filename),
            EXAMPLES_OUTDIR;
            kwargs...,
        )
    end
    return nothing
end

build_literate_pages()

# Pages vector for makedocs
example_docpages = Pair{String, String}[]
for (title, filename) in script_list
    mdfile = replace(filename, ".jl" => ".md")
    push!(example_docpages, title => joinpath(EXAMPLES_OUTDIR_RELATIVE, mdfile))
end

makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        ansicolor = true,
        collapselevel = 1,
        canonical = "https://tum-pik-esm.github.io/Terrarium.jl/stable/",
        size_threshold = 600_000,
        mathengine = Documenter.MathJax3(),
    ),
    sitename = "Terrarium.jl",
    authors = "Brian Groenke, Maximillian Galbrecht, Maha Badri, and Contributors",
    modules = [Terrarium],
    pages = [
        "Home" => "index.md",
        "Overview" => [
            "Basic concepts" => "introduction/basic_concepts.md",
            "Numerical core" => "introduction/numerical_core.md",
            "Mathematical formulation" => "introduction/mathematical_formulation.md",
        ],
        "Extending Terrarium" => [
            "Core interfaces" => "extending/core_interfaces.md",
        ],
        "Processes" => [
            "Soil" => [
                "Soil hydrology" => "processes/soil/soil_hydrology.md",
                "Soil energy" => "processes/soil/soil_energy.md",
                "Soil stratigraphy" => "processes/soil/soil_stratigraphy.md",
            ],
            "Vegetation" => [
                "Photosynthesis and gas exchange" => "processes/vegetation/photosynthesis.md",
                "Stomatal conductance" => "processes/vegetation/stomatal_conductance.md",
                "Plant available water" => "processes/vegetation/plant_available_water.md",
                "Autotrophic respiration" => "processes/vegetation/autotrophic_respiration.md",
                "Carbon dynamics" => "processes/vegetation/carbon_dynamics.md",
                "Vegetation dynamics" => "processes/vegetation/vegetation_dynamics.md",
                "Vegetation phenology" => "processes/vegetation/vegetation_phenology.md",
                "Root distribution" => "processes/vegetation/root_distribution.md",
            ],
            "Surface hydrology" => [
                "Surface hydrology" => "processes/surface_hydrology/surface_hydrology.md",
                "Canopy interception" => "processes/surface_hydrology/canopy_interception.md",
                "Evapotranspiration" => "processes/surface_hydrology/evapotranspiration.md",
                "Surface runoff" => "processes/surface_hydrology/surface_runoff.md",
            ],
            "Surface energy fluxes" => [
                "Surface energy balance" => "processes/surface_energy/surface_energy_balance.md",
                "Radiative fluxes" => "processes/surface_energy/radiative_fluxes.md",
                "Turbulent fluxes" => "processes/surface_energy/turbulent_fluxes.md",
                "Skin temperature and ground heat" => "processes/surface_energy/skin_temperature.md",
                "Albedo and emissivity" => "processes/surface_energy/albedo.md",
            ],
        ],
        "Models" => [
            "Soil model" => "models/soil_model.md",
        ],
        "Examples" => example_docpages,
        "Contributing" => "contributing.md",
        "API Reference" => "api_reference.md",
    ],
    draft = IS_DRAFT,
)

deployconfig = Documenter.auto_detect_deploy_system()

# remove gitignore from build files
# rm(joinpath(@__DIR__, "build", ".gitignore"))

deploydocs(
    repo = "github.com/NumericalEarth/Terrarium.jl.git",
    push_preview = true,
    versions = ["v0" => "v^", "v#.#", "dev" => "dev"],
    deploy_config = deployconfig,
)
