# Reactant + Enzyme reverse-mode autodiff test.
#
# Compiles a scalar loss and its reverse-mode gradient with Reactant and checks that the
# sensitivity of the loss w.r.t. the initial internal energy is finite and nonzero. Also exercises
# the extension's `checkpointing` argument: the gradient must be identical whether the traced loop
# stores every state (`checkpointing=false`) or only periodic checkpoints (`Periodic(n)`) — the
# checkpointing scheme changes memory use, not the math. Mirrors the Enzyme-through-Reactant
# pattern in Oceananigans' Reactant tests and the CPU example differentiating_terrarium.jl.

using Enzyme
using Statistics: mean

# Scalar loss: mean-square soil temperature after `nsteps`, advanced through `run_timesteps!`
# (the traced stepping loop; the extension's `ReactantState` method) with the given
# `checkpointing` scheme.
function ad_loss(integrator, Δt, nsteps, checkpointing)
    run_timesteps!(integrator, Δt, nsteps, checkpointing)
    return mean(interior(integrator.state.temperature) .^ 2)
end

# Reverse-mode gradient of `ad_loss`; the sensitivity to the initial prognostic internal energy
# accumulates in `dintegrator.state.internal_energy`.
function ad_grad!(integrator, dintegrator, Δt, nsteps, checkpointing)
    _, loss_value = Enzyme.autodiff(
        Enzyme.set_strong_zero(Enzyme.ReverseWithPrimal),
        ad_loss, Enzyme.Active,
        Enzyme.Duplicated(integrator, dintegrator),
        Enzyme.Const(Δt),
        Enzyme.Const(nsteps),
        Enzyme.Const(checkpointing),
    )
    return loss_value
end

# Compile the gradient for a given checkpointing scheme and return (loss, ∂loss/∂U₀ as a CPU array).
function reactant_gradient(config, NF, Δt, nsteps, checkpointing)
    integrator = build_integrator(Val(config), ReactantState(), NF)
    dintegrator = Enzyme.make_zero(integrator)
    compiled_grad! = Reactant.@compile raise = true raise_first = true sync = true ad_grad!(
        integrator, dintegrator, Δt, nsteps, checkpointing)
    loss_value = compiled_grad!(integrator, dintegrator, Δt, nsteps, checkpointing)
    return _scalar(loss_value), Array(interior(dintegrator.state.internal_energy))
end

@testset "Reactant + Enzyme autodiff" begin
    NF = DEFAULT_NF
    Δt = NF(600)
    nsteps = 5

    @testset "gradient is finite and nonzero" begin
        loss_value, dU = reactant_gradient(:soil_heat_column, NF, Δt, nsteps, false)
        println("autodiff: loss=$loss_value  max|∂loss/∂U₀|=$(maximum(abs, dU))")
        @test isfinite(loss_value)
        @test loss_value > 0            # mean-square temperature is positive
        @test all(isfinite, dU)         # gradient finite everywhere
        @test maximum(abs, dU) > 0      # loss actually depends on the initial state
    end

    @testset "checkpointing scheme leaves the gradient unchanged" begin
        _, dU_ref = reactant_gradient(:soil_heat_column, NF, Δt, nsteps, false)
        _, dU_ckpt = reactant_gradient(:soil_heat_column, NF, Δt, nsteps, Reactant.Periodic(2))
        @test all(isfinite, dU_ckpt)
        @test dU_ref ≈ dU_ckpt          # checkpointing changes memory use, not the gradient
    end
end
