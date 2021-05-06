
# [Data Models](@id DataModelAPI)

## Parsers

```@docs
parse_file
parse_dss
parse_opendss
parse_json
print_file
```

## Constructors

```@autodocs
Modules = [PowerModelsDistribution]
Private = false
Order = [:function]
Pages = ["components.jl"]
```

## Model Transformations

```@docs
transform_data_model
transform_solution
```

## Data Transformations

```@autodocs
Modules = [PowerModelsDistribution]
Private = false
Order = [:function]
Pages = ["transformations.jl"]
```

## Multinetworks

```@autodocs
Modules = [PowerModelsDistribution]
Private = false
Order = [:function]
Pages = ["multinetwork.jl"]
```

## Unit conversions

```@autodocs
Modules = [PowerModelsDistribution]
Private = false
Order = [:function]
Pages = ["units.jl"]
```

## Data Checking and Correction

```@autodocs
Modules = [PowerModelsDistribution]
Private = false
Order = [:function]
Filter = t -> startswith(string(t), "correct") || startswith(string(t), "check")
```

## Statistics and Analysis

```@docs
count_nodes
count_active_connections
count_active_terminals
identify_load_blocks
identify_blocks
identify_islands
calc_connected_components
```

## Helper Functions

```@docs
iseng
ismath
find_conductor_ids!
make_multiconductor!
discover_voltage_zones
calc_voltage_bases
apply_pmd!
get_pmd_data
```
