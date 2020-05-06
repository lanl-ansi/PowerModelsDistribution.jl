# Quick Start Guide

Once PowerModelsDistribution is installed, Ipopt is installed, and a network data file (e.g. `"case3_unbalanced.dss"` in the package folder under `./test/data`) has been acquired, an unbalanced AC Optimal Power Flow can be executed with,

```julia
using PowerModelsDistribution
using Ipopt

run_ac_mc_opf("case3_unbalanced.dss", with_optimizer(Ipopt.Optimizer))
```

## Parsing files

To parse an OpenDSS file into PowerModelsDistribution's default `ENGINEERING` format, use the `parse_file` command

```julia
eng = parse_file("case3_unbalanced.dss")
```

To examine the `MATHEMATICAL` model it is possible to transform the data model using the `transform_data_model` command, but this step is not necessary to run a problem.

```julia
math = transform_data_model(eng)
```

## Getting Results

The run commands in PowerModelsDistribution return detailed results data in the form of a dictionary. This dictionary can be saved for further processing as follows,

```julia
result = run_ac_mc_opf(eng, with_optimizer(Ipopt.Optimizer))
```

Alternatively, you can pass the file path string directly:

```julia
result = run_ac_mc_opf("case3_unbalanced.dss", with_optimizer(Ipopt.Optimizer))
```

## Accessing Different Formulations

The function "run_ac_mc_opf" is a short-hand for a more general formulation-independent OPF execution, "run_mc_opf".
For example, `run_ac_mc_opf` is equivalent to,

```julia
run_mc_opf(eng, ACPPowerModel, with_optimizer(Ipopt.Optimizer))
```

`ACPPowerModel` indicates an AC formulation in polar coordinates.  This more generic `run_mc_opf()` allows one to solve an OPF problem with any power network formulation implemented in PowerModels or PowerModelsDistribution.  For example, the SDPUBFPowerModel relaxation of unbalanced Optimal Power Flow (branch flow model) can be run with,

```julia
using SCS
run_mc_opf(eng, SDPUBFPowerModel, with_optimizer(SCS.Optimizer))
```

Note that you have to use a SDP-capable solver, e.g. the open-source solver SCS, to solve SDP models.

## Inspecting the Formulation

The following example demonstrates how to break a `run_mc_opf` call into seperate model building and solving steps.  This allows inspection of the JuMP model created by PowerModelsDistribution for the AC-OPF problem. Note that the `MATHEMATICAL` model must be passed to `instantiate_model`, so the data model must either be transformed with `transform_data_model` or parsed directly to a `MATHEMATICAL` model using the `data_model` keyword argument:

```julia
math = parse_file("case3_unbalanced.dss"; data_model=MATHEMATICAL)
pm = instantiate_model(math, ACPPowerModel, build_mc_opf; ref_extensions=[ref_add_arcs_trans!])
print(pm.model)
optimize_model!(pm, optimizer=with_optimizer(Ipopt.Optimizer))
```

## Examples

More examples of working with the engineering data model can be found in the `/examples` folder of the PowerModelsDistribution.jl repository.
