# Network Formulations

## Type Hierarchy

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

where `UBF` is an unbalanced variant of the __Branch Flow__ models from PowerModels. Models which do not contain `UBF` in their name are __Bus Injection__ Models _e.g._ `AbstractUnbalancedACPModel`. Finally, some linear unbalanced power flow models are available under the following hierarchy:

```julia
abstract type AbstractUnbalancedDCPModel <: AbstractUnbalancedActivePowerModel end
abstract type AbstractUnbalancedNFAModel <: AbstractUnbalancedDCPModel end
abstract type AbstractLPUBFModel <: AbstractNLPUBFModel end
abstract type LPUBFDiagModel <: AbstractLPUBFModel end
const LinDist3FlowModel = LPUBFDiagModel
```

## Power Models

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
```

## Optimization problem classes

- NLP (nonconvex): ACPUPowerModel, ACRUPowerModel, IVRUPowerModel
- SDP: SDPUBFPowerModel, SDPUBFKCLMXPowerModel
- SOC(-representable): SOCNLPUBFPowerModel, SOCConicUBFPowerModel
- Linear: LPUBFDiagPowerModel (LinDist3FlowPowerModel), DCPUPowerModel, NFAUPowerModel

## Matrix equations versus scalar equations

JuMP supports vectorized syntax, but not for nonlinear constraints. Therefore, certain formulations must be implemented in a scalar fashion. Other formulations can be written as matrix (in)equalities. The current implementations are categorized as follows:

- Scalar: ACPUPowerModel, ACRUPowerModel, IVRUPowerModel, DCPUPowerModel, NFAPowerMoel
- Matrix: SDPUBFPowerModel, SDPUBFKCLMXPowerModel, SOCNLPUBFPowerModel, SOCConicUBFPowerModel, LPUBFDiagPowerModel
