abstract type FileFormat end
struct NetCDF <: FileFormat end

"""$(TYPEDSIGNATURES)
File extension (including the leading dot) used by files of the given [`FileFormat`](@ref).
"""
file_extension(::NetCDF) = ".nc"

"""
    $TYPEDEF

Lightweight base type for Terrarium land data assets.
"""
abstract type AbstractLandAsset end

struct ERA5LandInvariants <: AbstractLandAsset end

artifact_name(::ERA5LandInvariants) = "era5-land-invariants"
varnames(::ERA5LandInvariants) = ["cvh", "lsm", "tvl", "cvl", "z", "slt", "dl", "cl", "si10", "tvh"]
format(::ERA5LandInvariants) = NetCDF()
indices(::ERA5LandInvariants) = (:, :, 5)
native_grid(::ERA5LandInvariants) = RingGrids.FullClenshawGrid(900)

struct ERA5LandForcings <: AbstractLandAsset end

artifact_name(::ERA5LandForcings) = "era5-land-forcings-N72"
varnames(::ERA5LandForcings) = ["t2m", "d2m", "tp", "sf", "sp", "ssrd", "strd", "u10", "v10"]
format(::ERA5LandForcings) = NetCDF()
indices(::ERA5LandForcings) = (:, :, :)
native_grid(::ERA5LandForcings) = RingGrids.FullGaussianGrid(72)

"""
    $TYPEDSIGNATURES

Download (if necessary) the given `asset` and return the path to its data file within the installed
artifact directory. The artifact is installed via [`get_artifact`](@ref) and the data file is located
by its [`file_extension`](@ref). Reading the file is left to the caller; use [`load_asset`](@ref) to
read a variable into a data `Field` wrapped on the asset's `native_grid`.
"""
function get_asset(asset::AbstractLandAsset)
    artifact_dir = get_artifact(artifact_name(asset))
    fmt = format(asset)
    path = locate_asset_file(artifact_dir, fmt)
    return path
end

function load_asset(asset::AbstractLandAsset, name::String; NF::Type = Float32, fill_value = NF(NaN))
    path = get_asset(asset)
    fmt = format(asset)
    grid = native_grid(asset)
    return load_asset(path, name, grid, fmt, NF; indices = indices(asset), fill_value)
end

"""
    $TYPEDSIGNATURES

Locate the single data file of the given `format` within an installed artifact directory `dir`,
searching recursively. Errors unless exactly one matching file is found.
"""
function locate_asset_file(dir::String, format::FileFormat)
    ext = file_extension(format)
    matches = String[]
    for (root, _, files) in walkdir(dir)
        for file in files
            endswith(file, ext) && push!(matches, joinpath(root, file))
        end
    end
    length(matches) == 1 || error("expected exactly one '$ext' file in artifact directory $dir, found $(length(matches))")
    return only(matches)
end

function get_artifact(name::String)
    project_root = pkgdir(Terrarium)
    # fallback to @__DIR__ if pkgdir fails
    project_root = isnothing(project_root) ? (@__DIR__) : project_root
    artifact_toml = joinpath(project_root, "Artifacts.toml")
    # Check if the artifact already exists in the TOML
    hash = Pkg.Artifacts.artifact_hash(name, artifact_toml)
    @assert !isnothing(hash) "no artifact with name $name found"
    Pkg.Artifacts.ensure_artifact_installed(name, artifact_toml)
    return Pkg.Artifacts.artifact_path(hash)
end
