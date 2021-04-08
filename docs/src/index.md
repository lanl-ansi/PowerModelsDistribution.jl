# PowerModelsDistribution.jl Documentation

```@meta
CurrentModule = PowerModelsDistribution
```

## Overview

PowerModelsDistribution.jl is a Julia/JuMP extension package to PowerModels.jl for modeling of Multi-Phase (with a focus on three-phase) power grids.

## Installation

The latest stable release of PowerModels can be installed using the Julia package manager with

```julia
Pkg.add("PowerModelsDistribution")
```

For the current development version, add the package with a reference to the master branch of our git repository

```julia
Pkg.add(name="PowerModelsDistribution", rev="master")
```

At least one solver is required for running PowerModelsDistribution.  The open-source solver Ipopt is recommended, as it is extremely fast, and can be used to solve a wide variety of the problems and network formulations provided in PowerModelsDistribution.  The Ipopt solver can be installed via the package manager with

```julia
Pkg.add("Ipopt")
```

Test that the package works by running

```julia
Pkg.test("PowerModelsDistribution")
```
