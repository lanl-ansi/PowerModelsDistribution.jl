# [Objectives](@id ObjectiveAPI)

```@autodocs
Modules = [PowerModelsDistribution]
Private = false
Order = [:function]
Filter = t -> startswith(string(t), "objective")
```

## Helpers

```@docs
calc_max_cost_index
simplify_cost_terms!
calc_pwl_points
calc_cost_pwl_lines
standardize_cost_terms!
```
