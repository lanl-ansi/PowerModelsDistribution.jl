# Power Flow Computations

The typical goal of PowerModelsDistribution is to build a JuMP model that is used to solve ditribution power network optimization problems.  The JuMP model abstraction enables PowerModelsDistribution to have state-of-the-art performance on a wide range of problem formulations.  That said, for the specific case of power flow computations, in some specific applications performance gains can be had by avoiding the JuMP model abstraction and solving the problem more directly.  To that end, PowerModelsDistribution includes Julia-native solvers for AC power flow in polar voltage coordinates.  This section provides an overview of the different power flow options that are available in PowerModelsDistribution and under what circumstances they may be beneficial.


## Generic Power Flow

The general purpose power flow solver in PowerModels is,

```@docs
solve_mc_pf
```

This function builds a JuMP model for a wide variety of the power flow formulations supported by PowerModelsDistribution.  For example it supports,
* `ACPUPowerModel` - a non-convex nonlinear AC power flow using complex voltages in polar coordinates
* `IVRUPowerModel` - a non-convex nonlinear AC power flow using current voltage variables in rectangular coordinates
* `AbstractUBFModels` - a non-convex nonlinear branch flow power flow  <!-- TODO type?


The `solve_mc_pf` solution method is both formulation and solver agnostic and can leverage the wide range of solvers that are available in the JuMP ecosystem.  Many of these solvers are commercial-grade, which in turn makes `solve_mc_pf` the most reliable power flow solution method in PowerModelsDistribution.

<!-- !!! note
    Use of `solve_mc_pf` is highly recommended over the other solution methods for increased robustness.  Applications that benefit from the Julia native solution methods are an exception to this general rule. -->


### Warm Starting

In some applications an initial guess of the power flow solution may be available.  In such a case, this information may be able to decrease a solver's time to convergence, especially when solving systems of nonlinear equations.  The `_start` postfix can be used in the network data to initialize the solver's variables and provide a suitable solution guess.  The most common values are as follows,

For each generator,
* `pg_start` - active power injection starting point
* `qg_start` - reactive power injection starting point

For each bus,
* `vm_start` - voltage magnitude starting point for the `ACPUPowerModel` model
* `va_start` - voltage angle starting point for the `ACPUPowerModel` model
* `vr_start` - real voltage starting point for the `IVRUPowerModel` model
* `vi_start` - imaginary voltage starting point for the `IVRUPowerModel` model

The following helper function can be used to use the solution point in the network data as the starting point for `solve_mc_pf`.  <!-- TODO correct?
```@docs
add_start_voltage!
```

!!! warning
    Warm starting a solver is a very delicate task and can easily result in degraded performance.  Using PowerModelsDistribution' default flat-start values is recommended before experimenting with warm starting a solver.


## Native Power Flow

The AC Power Flow problem is ubiquitous in power system analysis. The problem requires solving a system of nonlinear equations, usually via a Newton-Raphson type of algorithm.  In PowerModelsDistribution, the standard Julia library is used for solving this system of nonlinear equations.  The following function is used to solve Power Flow problem with voltages in polar coordinates.
```@docs
compute_pf
```
`compute_pf` will typically provide an identical result to `solve_mc_pf`. However, the existence of solution degeneracy around generator injection assignments and multiple power flow solutions can yield different results.  The primary advantage of `compute_pf` over `solve_mc_pf` is that it does not require building a JuMP model.  If the initial point for the Power Flow solution is near-feasible then `compute_pf` can result in a significant runtime saving by converging quickly and reducing data-wrangling and memory allocation overheads.  This initial guess is provided using the standard `_start` values.  The `add_start_voltage!` function provides a convenient way of setting a suitable starting point.

!!! tip
    If `compute_pf` fails to converge try `solve_mc_pf` instead.


The table below reports the accuracy of the native power flow with respect to OpenDSS native solver tested on three IEEE testcases:
| IEEE testcases | maximum voltage p.u difference with OpenDSS power flow solver |
| ---------------| ------------------------------------------------------------- |
| IEEE13         | 3.765096388188572e-6                                          |
| IEEE34         | 6.805369850332029e-8                                          |
| IEEE123        | 4.021326251365659e-8                                          |


<!-- ## Branch Flow Values

By default none of the Power Flow solvers produce branch flow values. If needed, these can be computed with the network data functions,
```@docs
calc_branch_flow_ac
calc_branch_flow_dc
```
Both of these methods require a complete network data with a valid voltage solution
for computing the branch flows.  For example, one common work flow to recover
branch flow values is,
```julia
result = solve_mc_pf(network, ...)
# check that the solver converged
update_data!(network, result["solution"])
flows = calc_branch_flow_ac(network)
update_data!(network, flows)
``` 
-->


## Network Admittance Matrix

Internally `compute_pf` utilizes an admittance matrix representation of the network data, which may be useful in other contexts.
The foundational type for the admittance matrix representations is `SparseArrays.SparseMatrixCSC`.
<!-- ```@docs
AdmittanceMatrix
``` -->
<!-- In the case of a full admittance matrix and simplified susceptance the type is `AdmittanceMatrix{Complex{Float64}}` and `AdmittanceMatrix{Float64}`, respectively.  -->

<!-- The following functions can be used to compute the admittance matrix and
susceptance matrix from PowerModelsDistribution network data.
```@docs
calc_admittance_matrix
calc_susceptance_matrix
``` 
-->

