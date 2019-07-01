# Quick Start Guide
Once PowerModelsDistribution is installed, Ipopt is installed, and a network data file (e.g. `"case5_c_r_a.m"` or `"case3_unbalanced.dss"` in the package folder under `./test/data`) has been acquired, an unbalanced AC Optimal Power Flow can be executed with,

```julia
using PowerModelsDistribution
using Ipopt

run_ac_tp_opf("case3_unbalanced.dss", IpoptSolver())
```

## Getting Results

The run commands in PowerModelsDistribution return detailed results data in the form of a dictionary. Results dictionaries from either Matpower-style `.m` or OpenDSS' `.dss` files will be identical in format. This dictionary can be saved for further processing as follows,

```julia
result = run_ac_tp_opf("case3_unbalanced.dss", IpoptSolver())
```


## Accessing Different Formulations

The function "run_ac_tp_opf" is a shorthands for a more general formulation-independent OPF execution, "run_tp_opf".
For example, `run_ac_tp_opf` is equivalent to,

```julia
using PowerModels
run_tp_opf("case3_unbalanced.dss", ACPPowerModel, IpoptSolver())
```

Note that PowerModels needs to be loaded to access formulations which are extended by PowerModelsDistribution, here "ACPPowerModel". The PowerModel "ACPPowerModel" indicates an AC formulation in polar coordinates.  This more generic `run_tp_opf()` allows one to solve an OPF problem with any power network formulation implemented in PowerModels or PowerModelsDistribution.  For example, the SDP relaxation of unbalanced Optimal Power Flow (branch flow model) can be run with,

```julia
using SCS
run_tp_opf_bf("case3_unbalanced.dss", SDPUBFPowerModel, SCSSolver())
```
Note that you have to use a SDP-capable solver, e.g. the open-source solver SCS, to solve SDP models.

## Inspecting the Formulation
The following example demonstrates how to break a `run_tp_opf` call into seperate model building and solving steps.  This allows inspection of the JuMP model created by PowerModelsDistribution for the AC-OPF problem,

```julia
data = PowerModelsDistribution.parse_file("case3_unbalanced.dss")
pm = PowerModels.build_model(data, ACPPowerModel, PowerModelsDistribution.post_tp_opf; multiconductor=true)
print(pm.model)
optimize_model!(pm, IpoptSolver())
```
