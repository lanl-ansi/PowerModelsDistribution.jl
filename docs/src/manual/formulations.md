# Unbalanced Network Formulations

## [`AbstractUnbalancedACPModel`](@ref AbstractUnbalancedACPModel)

Real-valued formulation from:

- Formulation without shunts: Mahdad, B., Bouktir, T., & Srairi, K. (2006). A three-phase power flow modelization: a tool for optimal location and control of FACTS devices in unbalanced power systems. In IEEE Industrial Electronics IECON (pp. 2238–2243).

## [`AbstractUnbalancedDCPModel`](@ref AbstractUnbalancedDCPModel)

Applying all of the standard DC linearization tricks to the [`AbstractUnbalancedACPModel`](@ref AbstractUnbalancedACPModel)

## [`SDPUBFModel`](@ref SDPUBFModel)

The BFM SDP relaxation as described in:

- Gan, L., & Low, S. H. (2014). Convex relaxations and linear approximation for optimal power flow in multiphase radial networks. In PSSC (pp. 1–9). Wroclaw, Poland. [doi:10.1109/PSCC.2014.7038399](https://doi.org/10.1109/PSCC.2014.7038399)

Note that this formulation is complex-valued and additional steps are needed to implement this in JuMP.

## [`SOCNLPUBFModel`](@ref SOCNLPUBFModel)

The starting point is `SDPUBFModel`. The SDP constraint can be relaxed to a set of SOC constraints, starting from either the real or complex form of the matrix on which the PSD-ness constraint is applied.

- Kim, S., Kojima, M., & Yamashita, M. (2003). Second order cone programming relaxation of a positive semidefinite constraint. Optimization Methods and Software, 18(5), 535–541. [doi:10.1080/1055678031000148696](https://doi.org/10.1080/1055678031000148696)
- Andersen, M. S., Hansson, A., & Vandenberghe, L. (2014). Reduced-complexity semidefinite relaxations of optimal power flow problems. IEEE Trans. Power Syst., 29(4), 1855–1863.

## [`SOCConicUBFModel`](@ref SOCConicUBFModel)

See `SOCNLPUBFModel`

## [`LPUBFDiagModel`](@ref LPUBFDiagModel)

This formulation has originally been developed by Sankur et al.

- Sankur, M. D., Dobbe, R., Stewart, E., Callaway, D. S., & Arnold, D. B. (2016). A linearized power flow model for optimization in unbalanced distribution systems. [arXiv:1606.04492v2](https://arxiv.org/abs/1606.04492v2)

This formulation is here cast as only considering the diagonal elements defined in `LPUBFFullModel`, which furthermore leads to the imaginary part of the lifted node voltage variable W being redundant and substituted out.

## [`FBSUBFPowerModel`](@ref FBSUBFPowerModel), [`FOTPUPowerModel`](@ref FOTPUPowerModel), [`FOTRUPowerModel`](@ref FOTRUPowerModel)

The linear FBS and FOT formulations as described in:

- Girigoudar, K., & Roald, L.A. (2021). Linearized  Three-Phase  Optimal  Power  Flow  Models for  Distribution  Grids  with  Voltage  Unbalance. 2021 IEEE Conference on Decision and Control (CDC).

# Unbalanced Network Formulation Type Hierarchy

PowerModelsDistribution has a rich model type hierarchy similiar to PowerModels. At the top of the type hierarchy we can distinguish between conic, active power only, and branch flow models:

```julia
abstract type AbstractUnbalancedConicModel <: AbstractPowerModel end
abstract type AbstractUnbalancedActivePowerModel <: AbstractPowerModel end
abstract type AbstractUBFModel <: AbstractPowerModel end
abstract type AbstractUBFQPModel <: AbstractUBFModel end
abstract type AbstractUBFConicModel <: AbstractUBFModel end
const AbstractUnbalancedConicModels = Union{AbstractUnbalancedConicModel, AbstractUBFConicModel}
```

Several nonlinear (non-convex) models are available at the top level:

```julia
abstract type AbstractUnbalancedACPModel <: AbstractPowerModel end
abstract type AbstractUnbalancedACRModel <: AbstractPowerModel end
abstract type AbstractUnbalancedIVRModel <: AbstractUnbalancedACRModel end
```

The following relaxations are available under these hierarchies:

```julia
abstract type AbstractNLPUBFModel <: AbstractUBFQPModel end
abstract type AbstractConicUBFModel <: AbstractUBFConicModel end
const AbstractUBFModels = Union{AbstractNLPUBFModel, AbstractConicUBFModel}

abstract type SDPUBFModel <: AbstractConicUBFModel end
abstract type SDPUBFKCLMXModel <: SDPUBFModel end
abstract type SOCNLPUBFModel <: AbstractNLPUBFModel end
abstract type SOCConicUBFModel <: AbstractConicUBFModel end
const SOCUBFModels = Union{SOCNLPUBFModel, SOCConicUBFModel}
```

where `UBF` is an unbalanced variant of the __Branch Flow__ models from PowerModels. Models which do not contain `UBF` in their name are __Bus Injection__ Models _e.g._ [`AbstractUnbalancedACPModel`](@ref AbstractUnbalancedACPModel). Finally, some linear unbalanced power flow models are available under the following hierarchy:

```julia
abstract type AbstractUnbalancedDCPModel <: AbstractUnbalancedActivePowerModel end
abstract type AbstractUnbalancedNFAModel <: AbstractUnbalancedDCPModel end
abstract type AbstractLPUBFModel <: AbstractNLPUBFModel end
abstract type LPUBFDiagModel <: AbstractLPUBFModel end
const LinDist3FlowModel = LPUBFDiagModel
abstract type FBSUBFModel <: AbstractLPUBFModel end
```

## Unbalanced Power Models

Each of these Models can be used as the type parameter for an UnbalancedPowerModel:

```julia
mutable struct ACPUPowerModel <: AbstractUnbalancedACPModel @pmd_fields end
mutable struct ACRUPowerModel <: AbstractUnbalancedACRModel @pmd_fields end
mutable struct DCPUPowerModel <: AbstractUnbalancedDCPModel @pmd_fields end
mutable struct NFAUPowerModel <: AbstractUnbalancedNFAModel @pmd_fields end

mutable struct SDPUBFPowerModel <: SDPUBFModel @pmd_fields end
mutable struct SDPUBFKCLMXPowerModel <: SDPUBFKCLMXModel @pmd_fields end

mutable struct SOCNLPUBFPowerModel <: SOCNLPUBFModel @pmd_fields end
mutable struct SOCConicUBFPowerModel <: SOCConicUBFModel @pmd_fields end

mutable struct LPUBFDiagPowerModel <: LPUBFDiagModel @pmd_fields end
const LinDist3FlowPowerModel = LPUBFDiagPowerModel
mutable struct FBSUBFPowerModel <: FBSUBFModel @pmd_fields end
mutable struct FOTPUPowerModel <: AbstractUnbalancedACPModel @pmd_fields end
mutable struct FOTRUPowerModel <: AbstractUnbalancedACRModel @pmd_fields end
```

## Optimization problem classes

- NLP (nonconvex): ACPUPowerModel, ACRUPowerModel, IVRUPowerModel
- SDP: SDPUBFPowerModel, SDPUBFKCLMXPowerModel
- SOC(-representable): SOCNLPUBFPowerModel, SOCConicUBFPowerModel
- Linear: LPUBFDiagPowerModel (LinDist3FlowPowerModel), FBSUBFPowerModel, FOTPUPowerModel, FOTRUPowerModel, DCPUPowerModel, NFAUPowerModel 

## Matrix equations versus scalar equations

JuMP supports vectorized syntax, but not for nonlinear constraints. Therefore, certain formulations must be implemented in a scalar fashion. Other formulations can be written as matrix (in)equalities. The current implementations are categorized as follows:

- Scalar: ACPUPowerModel, ACRUPowerModel, IVRUPowerModel, DCPUPowerModel, NFAPowerModel, FBSUBFPowerModel, FOTPUPowerModel, FOTRUPowerModel
- Matrix: SDPUBFPowerModel, SDPUBFKCLMXPowerModel, SOCNLPUBFPowerModel, SOCConicUBFPowerModel, LPUBFDiagPowerModel
