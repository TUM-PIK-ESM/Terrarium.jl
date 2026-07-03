module TerrariumReactantExt

# Reactant support for Terrarium.
#
# Design (see PLAN_reactant.md): a Terrarium model whose grid lives on `ReactantState` is built
# and initialized on the CPU (where eager KernelAbstractions kernel launches work), then the
# initialized state is transferred to the device. Only `timestep!`/`run!` are traced and compiled
# by XLA. Kernel launches inside the compiled step trace fine (this requires `CUDA` to be loaded
# alongside `Reactant`, even on CPU — it provides the KA↔Reactant glue).

using Terrarium
using Reactant
using Oceananigans

using Oceananigans: interior
using Oceananigans.Architectures: ReactantState, CPU, architecture, on_architecture

using Terrarium: Terrarium, AbstractLandGrid, ColumnGrid, ColumnRingGrid, AbstractModel,
    ModelIntegrator, StateVariables, get_field_grid, get_grid, get_timestepper,
    prognostic_names, auxiliary_names, input_names, namespace_names
using Terrarium: getproperties

const RARCH = ReactantState

# Land grids that live on the device.
const ReactantLandGrid{NF} = AbstractLandGrid{NF, <:RARCH}
const ReactantModel{NF} = AbstractModel{NF, <:ReactantLandGrid{NF}}

include("TerrariumReactantExt/grids.jl")
include("TerrariumReactantExt/transfer.jl")
include("TerrariumReactantExt/integrator.jl")

end # module
