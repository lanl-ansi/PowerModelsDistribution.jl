# Network Formulations

## Type Hierarchy
We begin with the top of the hierarchy, where we can distinguish between conic and non-conic power flow models.
```julia
PowerModels.AbstractConicModels = Union{PowerModels.AbstractConicModel, PowerModels.AbstractBFConicModel}
PowerModels.AbstractConicModel <: PowerModels.AbstractPowerModel
PowerModels.AbstractBFModel <: PowerModels.AbstractPowerModel
PowerModels.AbstractBFQPModel <: PowerModels.AbstractBFModel
PowerModels.AbstractBFConicModel <: PowerModels.AbstractBFModel
```

We begin with the top of the hierarchy, where we can distinguish between AC and DC power flow models.
```julia
PowerModels.AbstractACPModel <: PowerModels.AbstractPowerModel
PowerModels.AbstractDCPModel <: PowerModels.AbstractPowerModel
PowerModelsDistribution.AbstractNLPUBFModel <: PowerModels.AbstractBFQPModel
PowerModelsDistribution.AbstractConicUBFModel <: PowerModels.AbstractBFConicModel
PowerModelsDistribution.AbstractLPUBFModel <: PowerModelsDistribution.AbstractNLPUBFModel
```

From there, different Models are possible:
```julia
#Bus injection models:
PowerModels.AbstractACPModel <: PowerModels.AbstractPowerModel
PowerModels.AbstractDCPModel <: PowerModels.AbstractPowerModel

#Branch flow models:
PowerModelsDistribution.SDPUBFModel <: PowerModelsDistribution.AbstractConicUBFModel
PowerModelsDistribution.SOCNLPUBFModel <: PowerModelsDistribution.AbstractNLPUBFModel
PowerModelsDistribution.SOCConicUBFModel <: PowerModelsDistribution.AbstractConicUBFModel

PowerModelsDistribution.LPLinUBFModel <: PowerModels.AbstractBFModel
PowerModelsDistribution.LPfullUBFModel <: PowerModelsDistribution.AbstractLPUBFModel
PowerModelsDistribution.LPdiagUBFModel <: PowerModelsDistribution.AbstractLPUBFModel
```

## Power Models
Each of these Models can be used as the type parameter for a PowerModel:
```julia
mutable struct PowerModels.ACPPowerModel <: PowerModels.AbstractACPModel PowerModels.@pm_fields end
mutable struct PowerModels.DCPPowerModel <: PowerModels.AbstractDCPModel PowerModels.@pm_fields end

mutable struct PowerModels.SOCWRPowerModel <: PowerModels.SOCWRModel PowerModels.@pm_fields end

mutable struct PowerModelsDistribution.SDPUBFPowerModel <: PowerModelsDistribution.SDPUBFModel PowerModels.@pm_fields end
mutable struct PowerModelsDistribution.SOCNLPUBFPowerModel <: PowerModelsDistribution.SOCNLPUBFModel PowerModels.@pm_fields end
mutable struct PowerModelsDistribution.SOCConicUBFPowerModel <: PowerModelsDistribution.SOCConicUBFModel PowerModels.@pm_fields end

mutable struct PowerModelsDistribution.LPfullUBFPowerModel <: PowerModelsDistribution.LPfullUBFModel PowerModels.@pm_fields end
mutable struct PowerModelsDistribution.LPdiagUBFPowerModel <: PowerModelsDistribution.LPdiagUBFModel PowerModels.@pm_fields end
mutable struct PowerModelsDistribution.LPLinUBFPowerModel <: PowerModelsDistribution.LPLinUBFModel PowerModels.@pm_fields end
```

## Union Types

To support both conic and quadratically-constrained formulation variants for the unbalanced branch flow model, the union type `AbstractUBFModels` is defined. These formulations extend `AbstractBFModel` and are therefore also `AbstractWModels` (as defined in PowerModels proper).

```julia
AbstractUBFModels = Union{AbstractNLPUBFModel, AbstractConicUBFModel}
```

## Optimization problem classes
- NLP (nonconvex): ACPPowerModel
- SDP: SDPUBFPowerModel
- SOC(-representable): SOCWRPowerModel, SOCNLPUBFPowerModel, SOCConicUBFPowerModel
- Linear: LPfullUBFPowerModel, LPdiagUBFPowerModel, LPLinUBFPowerModel, DCPPowerModel


## Matrix equations versus scalar equations
JuMP supports vectorized syntax, but not for nonlinear constraints. Therefore, certain formulations must be implemented in a scalar fashion. Other formulations can be written as matrix (in)equalities. The current implementations are categorized as follows:
- Scalar: ACPPowerModel, DCPPowerModel, LPLinUBFPowerModel, SOCWRPowerModel
- Matrix: SDPUBFPowerModel, SOCNLPUBFPowerModel, SOCConicUBFPowerModel, LPfullUBFPowerModel, LPdiagUBFPowerModel
