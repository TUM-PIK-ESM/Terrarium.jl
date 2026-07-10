module TerrariumReactantExt

# Reactant support for Terrarium.
#
# Design: a Terrarium model whose grid lives on `ReactantState` is built
# and initialized on the CPU (where eager KernelAbstractions kernel launches work), then the
# initialized state is transferred to the device (might be revised in the future).
# Only `timestep!`/`run!` are traced and compiled by Reactant.
# Kernel launches inside the compiled step trace fine (this requires `CUDA` to be loaded
# alongside `Reactant`, even on CPU).

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
