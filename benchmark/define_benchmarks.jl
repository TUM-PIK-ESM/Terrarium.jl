# The benchmark suites. Keys are sorted for a deterministic README/documentation section order.
#
# `:bench200` is the headline sweep: its results feed the cross-architecture overview table and the
# comparison figures in the documentation. Keep that in mind before changing it — the stored results
# of other architectures were collected with the resolutions defined here.

benchmarks = Dict{Symbol, BenchmarkSuite}()

## Default resolution used wherever a suite varies something other than resolution.
const DEFAULT_NLAT_HALF = 24        # 4608 columns, ~3.75°
const DEFAULT_NZ = 30

## Horizontal resolution sweep: 128 columns (~22°) up to 165,888 columns (~0.6°).
const RESOLUTION_SWEEP = [4, 8, 16, 24, 32, 48, 72, 96, 144]

## All configurations at one resolution — the "what does each model cost" table.
benchmarks[:bench100] = BenchmarkSuite(
    title = "Model configurations, default resolution",
    nruns = length(CONFIGURATIONS),
    config = collect(CONFIGURATIONS),
)

## Headline sweep: the fully coupled land model across horizontal resolutions.
benchmarks[:bench200] = BenchmarkSuite(
    title = "Land model, horizontal resolution",
    nruns = length(RESOLUTION_SWEEP),
    config = fill(:land, length(RESOLUTION_SWEEP)),
    nlat_half = copy(RESOLUTION_SWEEP),
)

## Same sweep without vegetation, to isolate the cost of the canopy path.
benchmarks[:bench201] = BenchmarkSuite(
    title = "Land model without vegetation, horizontal resolution",
    nruns = length(RESOLUTION_SWEEP),
    config = fill(:land_no_vegetation, length(RESOLUTION_SWEEP)),
    nlat_half = copy(RESOLUTION_SWEEP),
)

## Soil heat conduction only: the Reactant reference configuration (see model_configurations.jl).
benchmarks[:bench202] = BenchmarkSuite(
    title = "Soil heat conduction, horizontal resolution",
    nruns = length(RESOLUTION_SWEEP),
    config = fill(:soil_heat, length(RESOLUTION_SWEEP)),
    nlat_half = copy(RESOLUTION_SWEEP),
)

## Vertical resolution at fixed horizontal resolution.
let nz = [10, 20, 30, 60, 100]
    benchmarks[:bench300] = BenchmarkSuite(
        title = "Number of soil layers",
        nruns = length(nz),
        nz = copy(nz),
    )
end

## Number format.
benchmarks[:bench400] = BenchmarkSuite(
    title = "Number format, Float32 vs Float64",
    nruns = 2,
    NF = [Float32, Float64],
)

## Time steppers. `IMEX` is not benchmarked: it needs an implicit sub-stepper and Terrarium does not
## yet provide a concrete one (`timestepping(…) == Implicit()` has no implementation).
benchmarks[:bench500] = BenchmarkSuite(
    title = "Time stepper",
    nruns = 2,
    model_kwargs = [(timestepper = ForwardEuler,), (timestepper = Heun,)],
    label = ["ForwardEuler", "Heun"],
)