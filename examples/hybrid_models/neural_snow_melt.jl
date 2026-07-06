# # Hybrid modeling: a neural-network snow melt trained with Reactant + Enzyme
#
# This example is a **technical demonstration** of *hybrid* (physics + machine-learning) modeling
# in Terrarium: we take the degree-day snow model from the [snow example](@ref snow_ddm_example)
# and replace its analytic melt law with a small neural network, which we train to reproduce the
# original melt using [Reactant.jl](https://github.com/EnzymeAD/Reactant.jl) (compilation) and
# [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl) (reverse-mode gradients). It is deliberately
# a toy problem — the NN just relearns a function we already know — but it shows the full
# machinery: generating training data from a Terrarium process, training an MLP with Reactant +
# Enzyme, and evaluating that MLP *inside a kernel* as a drop-in replacement process.
#
# ## What we are (and are not) doing
#
# * The degree-day melt is `M(T) = max(0, k·(T − T_melt))` — a rectified linear function of
#   temperature. A tiny MLP can represent it. We train the MLP offline on `(T, M(T))` pairs.
# * The `NeuralSnowMelt` process evaluates the trained MLP **per grid point inside a
#   KernelAbstractions kernel**. Because `Lux.apply` cannot currently be compiled inside a kernel
#   (it pulls in string handling that has no GPU/StableHLO lowering — see the note below), the
#   kernel evaluates the MLP's forward pass *explicitly* from its weight arrays.
# * Training runs on the Reactant device (Reactant + Enzyme). The hybrid model is then *run* on
#   the CPU, because Terrarium's input-driven simulations are not yet traced through Reactant; the
#   same kernel is device-agnostic and compiles under Reactant (that is how it is trained).

using Terrarium
using Reactant, CUDA          # CUDA is required by Reactant's kernel integration, even on CPU
using Lux, Enzyme, Optimisers, Random
using Statistics: mean, std
using KernelAbstractions: @index, @kernel
using Terrarium: launch!, XY, interior, AbstractProcess
using Terrarium: Adapt          # Adapt is re-exported by Terrarium (transitive dep)
import Lux.MLDataDevices as MLDD

# ## The reference degree-day model
#
# This is the analytic process we want to emulate (as in the snow example). We also expose its
# melt law `melt(process, T)` directly, so we can generate training targets from it.

@kwdef struct DegreeDaySnow{NF} <: AbstractProcess{NF}
    "Degree-day factor [m/(K s)]"
    k::NF = 5.0e-3 / (24 * 60 * 60)
    "Melting point of snow on the ground [°C]"
    T_melt::NF = 0.0
end

DegreeDaySnow(::Type{NF}) where {NF} = DegreeDaySnow{NF}()

"Analytic melt rate `M(T) = max(0, k (T − T_melt))` [m/s]."
@inline melt(p::DegreeDaySnow, T) = max(zero(T), p.k * (T - p.T_melt))

# ## The neural snow melt process
#
# `NeuralSnowMelt` mirrors `DegreeDaySnow` — same variables, same tendency `dS/dt = P − M` — but
# computes the melt `M` from an MLP with a single hidden layer (`1 → H → 1`, tanh). It stores the
# MLP's weight arrays and the input/output normalization constants used during training.

struct NeuralSnowMelt{NF, M, V} <: AbstractProcess{NF}
    "Hidden weight matrix (H×1) and bias (H)"
    W1::M
    b1::V
    "Output weight matrix (1×H) and bias (1)"
    W2::M
    b2::V
    "Input normalization for temperature (°C)"
    T_mean::NF
    T_std::NF
    "Output scale mapping the network's O(1) output back to a melt rate [m/s]"
    M_scale::NF
end
Adapt.@adapt_structure NeuralSnowMelt   # so the weight arrays move with the architecture

# Both processes declare exactly the same state variables, so either can be dropped into the model.
const SnowVars{NF} = Union{DegreeDaySnow{NF}, NeuralSnowMelt{NF}}
Terrarium.variables(::SnowVars{NF}) where {NF} = (
    Terrarium.input(:air_temperature, XY(), default = NF(0), units = u"°C", desc = "Near-surface air temperature in °C"),
    Terrarium.input(:snow_fall, XY(), default = NF(0), units = u"m/s", desc = "Snow fall rate in m/s"),
    Terrarium.prognostic(:snow_storage, XY(), units = u"m", desc = "Snow water equivalent in m"),
)

# ## The model
#
# A minimal `AbstractModel` holding one snow-melt process (either kind).

@kwdef struct SnowModel{NF, Grid <: Terrarium.AbstractLandGrid{NF}, Pro, Init, TS} <: Terrarium.AbstractModel{NF, Grid}
    grid::Grid
    snow_melt::Pro = DegreeDaySnow(eltype(grid))
    initializer::Init = DefaultInitializer(eltype(grid))
    timestepper::TS = ForwardEuler(eltype(grid))
end
Terrarium.variables(model::SnowModel) = Terrarium.variables(model.snow_melt)
Terrarium.compute_auxiliary!(state, ::SnowModel) = nothing
Terrarium.compute_tendencies!(state, model::SnowModel) = compute_tendencies!(state, model.grid, model.snow_melt)

# The snow storage is non-negative; clip it after each step (as in the snow example).
function Terrarium.timestep!(state::StateVariables, ::SnowModel, ::Terrarium.AbstractTimeStepper, Δt)
    interior(state.snow_storage) .= max.(interior(state.snow_storage), 0)
    return nothing
end

# ### Tendencies for the analytic process

function Terrarium.compute_tendencies!(state, grid, snow_melt::DegreeDaySnow)
    fields = get_fields(state, snow_melt)
    launch!(grid, XY, ddm_snow_flux_kernel!, state.tendencies, fields, snow_melt)
    return nothing
end
@kernel function ddm_snow_flux_kernel!(tend, grid, fields, snow_melt)
    i, j = @index(Global, NTuple)
    @inbounds tend.snow_storage[i, j, 1] = ddm_snow_flux(i, j, fields, snow_melt)
end
@inline function ddm_snow_flux(i, j, fields, snow_melt::DegreeDaySnow)
    P = @inbounds fields.snow_fall[i, j, 1]
    T = @inbounds fields.air_temperature[i, j, 1]
    return P - melt(snow_melt, T)
end

# ### Tendencies for the neural process
#
# The kernel evaluates the MLP explicitly: normalize `T`, run one hidden `tanh` layer and the
# linear output layer over the `H` hidden units, then un-scale to a melt rate. All operations are
# plain scalar arithmetic and array indexing, so the kernel compiles for CPU, GPU, and Reactant.

function Terrarium.compute_tendencies!(state, grid, snow_melt::NeuralSnowMelt)
    fields = get_fields(state, snow_melt)
    launch!(grid, XY, nn_snow_flux_kernel!, state.tendencies, fields, snow_melt)
    return nothing
end
@kernel function nn_snow_flux_kernel!(tend, grid, fields, snow_melt)
    i, j = @index(Global, NTuple)
    @inbounds tend.snow_storage[i, j, 1] = nn_snow_flux(i, j, fields, snow_melt)
end
@inline function nn_snow_flux(i, j, fields, p::NeuralSnowMelt)
    P = @inbounds fields.snow_fall[i, j, 1]
    T = @inbounds fields.air_temperature[i, j, 1]
    Tn = (T - p.T_mean) / p.T_std
    ## explicit MLP forward pass: y = W2 * tanh.(W1 * Tn .+ b1) .+ b2
    acc = @inbounds p.b2[1]
    for h in 1:length(p.b1)
        z = @inbounds p.W1[h, 1] * Tn + p.b1[h]
        acc += @inbounds p.W2[1, h] * tanh(z)
    end
    M = acc * p.M_scale          # un-scale network output to a melt rate
    return P - M
end

# ## Generate training data from the reference process
#
# We sample temperatures over a realistic range and evaluate the degree-day melt to get targets.
# Inputs and targets are normalized so the network trains on O(1) quantities; the normalization
# constants are stored in the `NeuralSnowMelt` process and undone at inference.

const NF = Float32
reference = DegreeDaySnow(NF)

T_train = reshape(collect(range(NF(-30), NF(30), length = 512)), 1, :)   # (1, N) row of temperatures
M_train = melt.(Ref(reference), T_train)                                  # (1, N) target melt rates

T_mean = mean(T_train)
T_std = std(T_train)
M_scale = maximum(M_train)                     # melt is O(1e-6); scale targets to O(1) for training
T_norm = (T_train .- T_mean) ./ T_std
M_norm = M_train ./ M_scale

# ## Train the MLP with Reactant + Enzyme
#
# We build a small MLP, move it (and the data) to the Reactant device, and train with
# `Lux.Training.single_train_step!` using the Enzyme AD backend. On the Reactant device this
# compiles the train step with XLA; each step computes the loss and its gradient with Enzyme and
# updates the parameters with Adam.

const H = 16
mlp = Lux.Chain(Lux.Dense(1 => H, tanh), Lux.Dense(H => 1))

rdev = MLDD.reactant_device()
rng = Random.Xoshiro(0)
ps, st = Lux.setup(rng, mlp) |> rdev
data = (rdev(T_norm), rdev(M_norm))

function train_mlp(mlp, ps, st, data; epochs = 2000, lr = NF(1.0e-2))
    tstate = Lux.Training.TrainState(mlp, ps, st, Optimisers.Adam(lr))
    local loss
    for epoch in 1:epochs
        _, loss, _, tstate = Lux.Training.single_train_step!(AutoEnzyme(), MSELoss(), data, tstate)
        (epoch == 1 || epoch % 500 == 0) &&
            println("  epoch $epoch   loss = $(Reactant.to_number(loss))")
    end
    return tstate
end

println("Training the neural melt law (Reactant + Enzyme):")
tstate = train_mlp(mlp, ps, st, data)

# ## Assemble the trained neural process
#
# We move the trained weights back to the CPU (we run the hybrid model on the CPU below) and build
# the `NeuralSnowMelt` process. Each Lux `Dense` layer stores a `weight` matrix and `bias` vector.

cdev = MLDD.cpu_device()
ps_trained = cdev(tstate.parameters)
neural = NeuralSnowMelt(
    ps_trained.layer_1.weight, ps_trained.layer_1.bias,   # W1 (H×1), b1 (H)
    ps_trained.layer_2.weight, ps_trained.layer_2.bias,   # W2 (1×H), b2 (1)
    T_mean, T_std, M_scale,
)

# Quick check of the learned melt law against the analytic one (the kernel and this host-side
# forward compute the same thing):
nn_melt(p::NeuralSnowMelt, T) = (Tn = (T - p.T_mean) / p.T_std;
    (p.b2[1] + sum(p.W2[1, h] * tanh(p.W1[h, 1] * Tn + p.b1[h]) for h in 1:length(p.b1))) * p.M_scale)

T_check = range(NF(-30), NF(30), length = 200)
melt_error = maximum(abs, nn_melt.(Ref(neural), T_check) .- melt.(Ref(reference), T_check))
println("Max |M_nn − M| over [-30, 30] °C: ", melt_error, "  (analytic max melt ", M_scale, ")")

# ## Run both models and compare
#
# We build a small global column grid, force it with a latitude-dependent temperature and uniform
# snowfall, start fully snow-covered, and integrate both the analytic and neural snow models. The
# two snow fields should agree to within the network's fit error.

using RingGrids

rings = RingGrids.FullGaussianGrid(8)
grid = ColumnRingGrid(UniformSpacing(N = 1), rings)         # CPU, all-land
_lons, _lats = RingGrids.get_lonlats(rings)                  # per-column lon/lat in radians
air_temperature = NF.(25 .- 30 .* abs.(_lats) ./ (π / 2))    # warm equator (~25°C), cold poles (~-5°C)
snow_fall = fill(NF(1.0e-7), length(_lats))                  # uniform light snowfall [m/s]

inputs = InputSources(
    InputSource(grid, air_temperature, name = :air_temperature, units = u"°C"),
    InputSource(grid, snow_fall, name = :snow_fall, units = u"m/s"),
)
initializers = (snow_storage = NF(0.5),)                      # start with 0.5 m everywhere

function run_snow(process; period = Day(30), Δt = 3600.0)
    model = SnowModel(grid; snow_melt = process, timestepper = ForwardEuler(NF))
    integrator = initialize(model; inputs, initializers)
    run!(integrator; period, Δt)
    return Array(interior(integrator.state.snow_storage))[:, 1, 1]
end

snow_ddm = run_snow(reference)
snow_nn = run_snow(neural)

snow_diff = maximum(abs, snow_ddm .- snow_nn)
println("Ran both snow models for 30 days on $(length(lats)) columns.")
println("Max |S_ddm − S_nn| after 30 days: ", snow_diff, " m  (mean snow ", mean(snow_ddm), " m)")

# The neural snow model reproduces the analytic degree-day model to within the network's fit
# error — the physics-based melt law has been replaced by a trained, kernelized neural network,
# which is the essential building block of hybrid land modeling in Terrarium.
