# Base

```@docs
ref
var
ids
con
sol
nws
nw_ids
```

## Helper functions

```@docs
set_lower_bound
set_upper_bound
comp_start_value
```

## Ref Creation Functions

```@autodocs
Modules = [PowerModelsDistribution]
Private = false
Order = [:function]
Filter = t -> startswith(string(t), "ref_")
```

## InfrastructureModels Extensions

```@docs
PowerModelsDistribution._IM.solution_preprocessor
PowerModelsDistribution._IM.build_solution_values
```
