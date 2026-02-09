module TerrariumCheckpointingExt

using DocStringExtensions
using Checkpointing

using Terrarium
using Terrarium: run!

"""
    $SIGNATURES

Run the simulation for `steps` with timestep size `Δt` (in seconds or Dates.Period).

Uses a `checkpointing_scheme` from Checkpointing.jl to checkpoint the integrator state for automatic differentiation with Enzyme.jl. 
"""
function Terrarium.run!(
        integrator::ModelIntegrator,
        checkpointing_scheme::Checkpointing.Scheme,
        steps::Int,
        Δt = default_dt(timestepper(integrator))
    )
    Δt = Terrarium.convert_dt(Δt)

    @ad_checkpoint checkpointing_scheme for _ in 1:steps
        timestep!(integrator, Δt, finalize = false)
    end

    # Update auxiliary variables for final timestep
    compute_auxiliary!(integrator.state, integrator.model)
    return integrator
end

end
