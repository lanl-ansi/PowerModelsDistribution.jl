# Installation Guide

From Julia, PowerModelsDistribution is installed using the built-in package manager:

```julia
import Pkg
Pkg.add("PowerModelsDistribution")
```

Or, within the Julia REPL:

```julia
]add PowerModelsDistribution
```

## Installing an Optimizer

PowerModelsDistribution depends on optimizers to solve Optimization problems, _e.g._, [`solve_mc_opf`](@ref solve_mc_opf). The table below lists the optimizer packages that have been tested with PowerModelsDistribution, and have been found to work by the team; this list is not exhaustive, there are probably more optimizers that will work.

Install an optimizer using the Julia Package Manager, _e.g._,

```julia
import Pkg
Pkg.add("Ipopt")
```

### Known Working Optimizers

The table below contains a truncated list of optimizers from the JuMP documentation that have been used by the development team and are confirmed to work with our package. There may be other Optimizers that work, and the Optimizers listed below are not guaranteed to work with all problems; they must be selected for the appropriate problems.

| Solver                                                  | Julia Package                                         | Installation | License | Supports                  |
|---------------------------------------------------------|-------------------------------------------------------|--------------|---------|---------------------------|
| [Artelys Knitro](https://www.artelys.com/knitro)        | [KNITRO.jl](https://github.com/jump-dev/KNITRO.jl)    | Manual       | Comm.   | (MI)LP, (MI)SOCP, (MI)NLP |
| [Cbc](https://github.com/coin-or/Cbc)                   | [Cbc.jl](https://github.com/jump-dev/Cbc.jl)          |              | EPL     | (MI)LP                    |
| [CPLEX](https://www.ibm.com/analytics/cplex-optimizer/) | [CPLEX.jl](https://github.com/jump-dev/CPLEX.jl)      | Manual       | Comm.   | (MI)LP, (MI)SOCP          |
| [Gurobi](https://gurobi.com)                            | [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl)    | Manual       | Comm.   | (MI)LP, (MI)SOCP          |
| [Ipopt](https://github.com/coin-or/Ipopt)               | [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl)      |              | EPL     | LP, QP, NLP               |
| [Juniper.jl](https://github.com/lanl-ansi/Juniper.jl)   | [Juniper.jl](https://github.com/lanl-ansi/Juniper.jl) |              | MIT     | (MI)SOCP, (MI)NLP         |
| [SCS](https://github.com/cvxgrp/scs)                    | [SCS.jl](https://github.com/jump-dev/SCS.jl)          |              | MIT     | LP, SOCP, SDP             |

Where:
- LP = Linear programming
- QP = Quadratic programming
- SOCP = Second-order conic programming (including problems with convex quadratic constraints and/or objective)
- NLP = Nonlinear programming
- SDP = Semidefinite programming
- (MI)XXX = Mixed-integer equivalent of problem type `XXX`

For a complete list of JuMP supported optimizers, see [JuMP Documentation](https://jump.dev/JuMP.jl/stable/installation/).

## Unsatisfiable requirements detected

Did you get an error like `Unsatisfiable requirements detected for package D [756980fe]:`?

The Pkg documentation has a [section on how to understand and manage these conflicts](https://julialang.github.io/Pkg.jl/v1/managing-packages/#conflicts).
