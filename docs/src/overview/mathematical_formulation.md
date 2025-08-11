# Mathematical formulation

Terrarium.jl strictly follows a philosophy of dynamical modeling via discretized differential equations. This means that all physical processes in the model must be written as terms in a set of continuous-time equations of the form:
```math
\frac{\partial u(x,t)}{\partial t} = G(u(x,t)) + F(x,t)
```
where $u(x,t)$ is a field describing the prognostic state of the system, $t$ is time, $F$ is a focing term, and $G$ is a (differentiable) function which computes the tendencies, i.e. the change in the prognostic state at the current time, $t$.

In many areas of science and engineering, this approach to modeling is fairly standard. It can be contrasted, however, with discrete-time dynamical modeling where state changes are computed according to a series of discrete update rules, $\mathbf{D}$, with some fixed time resolution, i.e:
```math
\mathbf{u}_{t+1} = \mathbf{D}(\mathbf{u}_{t})
```

## Practical implications

While the above modeling philosophy has many advantages, it also places some practical restrictions on how we code the model physics. The most important constraints are:
- Prognostic (i.e. time-integrated) state variables must be clearly and coherently distinguished from all other auxiliary variables derived from the prognostic state.
- Prognostic variables of the system should only be updated by the timestepper (or callback) and should not be otherwise modified by the physical processes within a single timestep. Note that this effectively rules out all classical “bucket schemes” for soil hydrology which rely on non-physical, instantaneous routing of water between soil layers.
- For physical coherence, tendencies must be computed only based on the current state of the system. Similarly, all non-prognostic (auxiliary) variables should be derived only from the prognostic state and/or forcings; they should not depend on non-prognostic values from previous timesteps, except in special cases where previous values are used only for computational efficiency (e.g. iterative solvers).
