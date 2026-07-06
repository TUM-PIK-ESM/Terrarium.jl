# Reactant + Enzyme reverse-mode autodiff smoke test.
#
# Compiles a scalar loss and its reverse-mode gradient with Reactant and checks that the
# sensitivity of the loss w.r.t. the initial internal energy is finite and nonzero. Mirrors the
# Enzyme-through-Reactant pattern in Oceananigans' Reactant tests and the CPU autodiff example
# (examples/autodiff/differentiating_terrarium.jl).

using Enzyme
using Statistics: mean
using Reactant: @trace

# Scalar loss: mean-square soil temperature after `nsteps`. The time loop uses the raw per-step
# API (public `timestep!(integrator, timestepper, Δt)` + `compute_auxiliary!`, NOT the
# auto-compiling `run!`, which would nest `@compile`) inside a `@trace checkpointing=true` loop,
# so Reactant stores periodic checkpoints for the reverse pass instead of every intermediate state.
function ad_loss(integrator, Δt, nsteps)
    timestepper = get_timestepper(integrator.model)
    @trace mincut = true checkpointing = true track_numbers = false for _ in 1:nsteps
        timestep!(integrator, timestepper, Δt)
        compute_auxiliary!(integrator.state, integrator.model)
    end
    return mean(interior(integrator.state.temperature) .^ 2)
end

# Reverse-mode gradient of `ad_loss` w.r.t. the integrator state; the sensitivity to the initial
# prognostic internal energy accumulates in `dintegrator.state.internal_energy`.
function ad_grad!(integrator, dintegrator, Δt, nsteps)
    _, loss_value = Enzyme.autodiff(
        Enzyme.set_strong_zero(Enzyme.ReverseWithPrimal),
        ad_loss, Enzyme.Active,
        Enzyme.Duplicated(integrator, dintegrator),
        Enzyme.Const(Δt),
        Enzyme.Const(nsteps),
    )
    return loss_value
end

@testset "Reactant + Enzyme autodiff" begin
    NF = DEFAULT_NF
    integrator = build_integrator(Val(:soil_heat_column), ReactantState(), NF)
    dintegrator = Enzyme.make_zero(integrator)
    Δt = NF(600)
    nsteps = 5

    compiled_grad! = Reactant.@compile raise = true raise_first = true sync = true ad_grad!(
        integrator, dintegrator, Δt, nsteps)
    loss_value = compiled_grad!(integrator, dintegrator, Δt, nsteps)

    dU = Array(interior(dintegrator.state.internal_energy))
    println("autodiff: loss=$(_scalar(loss_value))  max|∂loss/∂U₀|=$(maximum(abs, dU))")

    @test isfinite(_scalar(loss_value))
    @test _scalar(loss_value) > 0        # mean-square temperature is positive
    @test all(isfinite, dU)              # gradient is finite everywhere
    @test maximum(abs, dU) > 0           # loss actually depends on the initial state
end
