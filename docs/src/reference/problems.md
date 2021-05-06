# [Problems](@id ProblemAPI)

## Solvers

```@autodocs
Modules = [PowerModelsDistribution]
Private = false
Order = [:function]
Filter = t -> startswith(string(t), "solve")
```

## Builders

```@autodocs
Modules = [PowerModelsDistribution]
Private = false
Order = [:function]
Filter = t -> startswith(string(t), "build")
```

## Model Instantiation

```@docs
instantiate_mc_model
```

## Solution Helpers

```@autodocs
Modules = [PowerModelsDistribution]
Private = false
Order = [:function]
Filter = t -> startswith(string(t), "sol_")
```

## DEPRECIATED Solver functions

```@autodocs
Modules = [PowerModelsDistribution]
Private = false
Order = [:function]
Filter = t -> startswith(string(t), "run")
```
