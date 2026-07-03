# CPU-vs-Reactant correctness tests for Terrarium.
#
# Run from the repo root with the dedicated environment (Reactant is NOT a dependency of the
# main test suite):
#
#     julia --project=test/reactant test/reactant/runtests.jl
#
# See PLAN_reactant.md for the overall design.

using Terrarium
using Reactant
using CUDA   # required by Reactant's KernelAbstractions integration, even on CPU
using Test

# Tolerances: XLA reorders floating-point ops, so exact equality is not expected. Start loose
# and tighten empirically per model (pure heat conduction should be far quieter than this).
const DEFAULT_NF = Float32
const NSTEPS = 10
const RTOL = 1.0e-3
const ATOL = 1.0e-6

include("correctness.jl")
include("setup.jl")

@testset "Terrarium CPU vs Reactant" begin
    test_model(:soil_heat_column)
    test_model(:soil_heat_global)
end
