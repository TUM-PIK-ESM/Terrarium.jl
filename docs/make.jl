using ArgParse
using Documenter
using Literate

using Terra

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
        canonical = "https://tum-pik-esm.github.io/TerraDocumentation/stable/",
        size_threshold = 600_000
    ),      # in bytes
    sitename = "Terra.jl",
    authors = "Brian Groenke, Maximillian Galbrecht, Maha Badri, and Contributors",
    modules = [Terra],
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
        "Contributing" => "contributing.md",
        "API Reference" => "api_reference.md",
    ],
    draft=IS_DRAFT,
)

deployconfig = Documenter.auto_detect_deploy_system()

# remove gitignore from build files
# rm(joinpath(@__DIR__, "build", ".gitignore"))

deploydocs(
       repo="github.com/TUM-PIK-ESM/Terra.jl.git",
       push_preview = true,
       versions = ["v0" => "v^", "v#.#", "dev" => "dev"],
       deploy_config = deployconfig,
)
