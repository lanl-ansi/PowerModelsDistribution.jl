# Network Formulations

## Type Hierarchy

PowerModelsDistribution shares a rich model type hierarchy with PowerModels. The relevant abstract models from PowerModels are documented here for context. At the top of the type hierarchy, starting in PowerModels, we can distinguish between conic, active power only, and branch flow models:

```julia
abstract type PowerModels.AbstractConicModel <: PowerModels.AbstractPowerModel end
abstract type PowerModels.AbstractActivePowerModel <: PowerModels.AbstractPowerModel end
abstract type PowerModels.AbstractBFModel <: PowerModels.AbstractPowerModel end
abstract type PowerModels.AbstractBFQPModel <: PowerModels.AbstractBFModel end
abstract type PowerModels.AbstractBFConicModel <: PowerModels.AbstractBFModel end
const PowerModels.AbstractConicModels = Union{PowerModels.AbstractConicModel, PowerModels.AbstractBFConicModel}
```

Several nonlinear (non-convex) models are available at the top level:

```julia
abstract type PowerModels.AbstractACPModel <: PowerModels.AbstractPowerModel end
abstract type PowerModels.AbstractACRModel <: PowerModels.AbstractPowerModel end
abstract type PowerModels.AbstractIVRModel <: PowerModels.AbstractACRModel end
```

In PowerModelsDistribution, the following relaxations are available under these hierarchies:

```julia
abstract type PowerModelsDistribution.AbstractNLPUBFModel <: PowerModels.AbstractBFQPModel end
abstract type PowerModelsDistribution.AbstractConicUBFModel <: PowerModels.AbstractBFConicModel end
const PowerModelsDistribution.AbstractUBFModels = Union{PowerModelsDistribution.AbstractNLPUBFModel, PowerModelsDistribution.AbstractConicUBFModel}

abstract type PowerModelsDistribution.SDPUBFModel <: PowerModelsDistribution.AbstractConicUBFModel end
abstract type PowerModelsDistribution.SDPUBFKCLMXModel <: PowerModelsDistribution.SDPUBFModel end
abstract type PowerModelsDistribution.SOCNLPUBFModel <: PowerModelsDistribution.AbstractNLPUBFModel end
abstract type PowerModelsDistribution.SOCConicUBFModel <: PowerModelsDistribution.AbstractConicUBFModel end
const PowerModelsDistribution.SOCUBFModels = Union{PowerModelsDistribution.SOCNLPUBFModel, PowerModelsDistribution.SOCConicUBFModel}
```

where `UBF` is an unbalanced variant of the __Branch Flow__ models from PowerModels. Models which do not contain `UBF` in their name are __Bus Injection__ Models _e.g._ `AbstractACPModel`. Finally, some linear unbalanced power flow models are available under the following hierarchy:

```julia
abstract type PowerModels.AbstractDCPModel <: PowerModels.AbstractActivePowerModel end
abstract type PowerModels.AbstractNFAModel <: PowerModels.AbstractDCPModel end
abstract type PowerModelsDistribution.AbstractLPUBFModel <: PowerModelsDistribution.AbstractNLPUBFModel end
abstract type PowerModelsDistribution.LPUBFDiagModel <: PowerModelsDistribution.AbstractLPUBFModel end
const PowerModelsDistribution.LinDist3FlowModel = PowerModelsDistribution.LPUBFDiagModel
```

## Power Models

Each of these Models can be used as the type parameter for a PowerModel:

```julia
mutable struct PowerModels.ACPPowerModel <: PowerModels.AbstractACPModel PowerModels.@pm_fields end
mutable struct PowerModels.ACRPowerModel <: PowerModels.AbstractACRModel PowerModels.@pm_fields end
mutable struct PowerModels.DCPPowerModel <: PowerModels.AbstractDCPModel PowerModels.@pm_fields end
mutable struct PowerModels.NFAPowerModel <: PowerModels.AbstractNFAModel PowerModels.@pm_fields end

mutable struct PowerModelsDistribution.SDPUBFPowerModel <: PowerModelsDistribution.SDPUBFModel PowerModels.@pm_fields end
mutable struct PowerModelsDistribution.SDPUBFKCLMXPowerModel <: PowerModelsDistribution.SDPUBFKCLMXModel PowerModels.@pm_fields end

mutable struct PowerModelsDistribution.SOCNLPUBFPowerModel <: PowerModelsDistribution.SOCNLPUBFModel PowerModels.@pm_fields end
mutable struct PowerModelsDistribution.SOCConicUBFPowerModel <: PowerModelsDistribution.SOCConicUBFModel PowerModels.@pm_fields end

mutable struct PowerModelsDistribution.LPUBFDiagPowerModel <: PowerModelsDistribution.LPUBFDiagModel PowerModels.@pm_fields end
const PowerModelsDistribution.LinDist3FlowPowerModel = PowerModelsDistribution.LPUBFDiagPowerModel
```

## Optimization problem classes

- NLP (nonconvex): ACPPowerModel, ACRPowerModel, IVRPowerModel
- SDP: SDPUBFPowerModel, SDPUBFKCLMXPowerModel
- SOC(-representable): SOCNLPUBFPowerModel, SOCConicUBFPowerModel
- Linear: LPUBFDiagPowerModel (LinDist3FlowPowerModel), DCPPowerModel, NFAPowerModel

## Matrix equations versus scalar equations

JuMP supports vectorized syntax, but not for nonlinear constraints. Therefore, certain formulations must be implemented in a scalar fashion. Other formulations can be written as matrix (in)equalities. The current implementations are categorized as follows:

- Scalar: ACPPowerModel, ACRPowerModel, IVRPowerModel, DCPPowerModel, NFAPowerMoel
- Matrix: SDPUBFPowerModel, SDPUBFKCLMXPowerModel, SOCNLPUBFPowerModel, SOCConicUBFPowerModel, LPUBFDiagPowerModel
