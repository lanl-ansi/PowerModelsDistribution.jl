# ThreePhasePowerModels.jl Documentation

```@meta
CurrentModule = ThreePhasePowerModels
```

## Overview

ThreePowerModels.jl is a Julia/JuMP extension package to PowerModels.jl for modeling of Multi-Phase (with a focus on three-phase) power grids. 
## Installation

The latest stable release of PowerModels can be installed using the Julia package manager with

```julia
Pkg.add("ThreePhasePowerModels")
```

For the current development version, "checkout" this package with

```julia
Pkg.checkout("ThreePhasePowerModels")
```

At least one solver is required for running ThreePhasePowerModels.  The open-source solver Ipopt is recommended, as it is extremely fast, and can be used to solve a wide variety of the problems and network formulations provided in ThreePhasePowerModels.  The Ipopt solver can be installed via the package manager with

```julia
Pkg.add("Ipopt")
```

Test that the package works by running

```julia
Pkg.test("ThreePhasePowerModels")
```
