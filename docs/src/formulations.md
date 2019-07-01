# Network Formulations

## Type Hierarchy
We begin with the top of the hierarchy, where we can distinguish between conic and non-conic power flow models.
```julia
PowerModels.AbstractConicForms = Union{PowerModels.AbstractConicPowerFormulation, PowerModels.AbstractBFConicForm}
PowerModels.AbstractConicPowerFormulation <: PowerModels.AbstractPowerFormulation
PowerModels.AbstractBFForm <: PowerModels.AbstractPowerFormulation
PowerModels.AbstractBFQPForm <: PowerModels.AbstractBFForm
PowerModels.AbstractBFConicForm <: PowerModels.AbstractBFForm
```

We begin with the top of the hierarchy, where we can distinguish between AC and DC power flow models.
```julia
PowerModels.AbstractACPForm <: PowerModels.AbstractPowerFormulation
PowerModels.AbstractDCPForm <: PowerModels.AbstractPowerFormulation
PowerModels.AbstractWRForm <: PowerModels.AbstractPowerFormulation
PowerModelsDistribution.AbstractNLPUBFForm <: PowerModels.AbstractBFQPForm
PowerModelsDistribution.AbstractConicUBFForm <: PowerModels.AbstractBFConicForm
PowerModelsDistribution.AbstractLPUBFForm <: PowerModelsDistribution.AbstractNLPUBFForm
```

From there, different forms are possible:
```julia
#Bus injection models:
PowerModels.StandardACPForm <: PowerModels.AbstractACPForm
PowerModels.StandardDCPForm <: PowerModels.AbstractDCPForm
PowerModels.SOCWRForm <: PowerModels.AbstractWRForm

#Branch flow models:
PowerModelsDistribution.SDPUBFForm <: PowerModelsDistribution.AbstractConicUBFForm
PowerModelsDistribution.SOCNLPUBFForm <: PowerModelsDistribution.AbstractNLPUBFForm
PowerModelsDistribution.SOCConicUBFForm <: PowerModelsDistribution.AbstractConicUBFForm

PowerModelsDistribution.LPLinUBFForm <: PowerModels.AbstractBFForm
PowerModelsDistribution.LPfullUBFForm <: PowerModelsDistribution.AbstractLPUBFForm
PowerModelsDistribution.LPdiagUBFForm <: PowerModelsDistribution.AbstractLPUBFForm
```

## Power Models
Each of these forms can be used as the type parameter for a PowerModel:
```julia
PowerModels.ACPPowerModel = GenericPowerModel{PowerModels.StandardACPForm}
PowerModels.DCPPowerModel = GenericPowerModel{PowerModels.StandardDCPForm}

PowerModels.SOCWRPowerModel = GenericPowerModel{PowerModels.SOCWRForm}

PowerModelsDistribution.SDPUBFPowerModel = GenericPowerModel{PowerModelsDistribution.SDPUBFForm}
PowerModelsDistribution.SOCNLPUBFPowerModel = GenericPowerModel{PowerModelsDistribution.SOCNLPUBFForm}
PowerModelsDistribution.SOCConicUBFPowerModel = GenericPowerModel{PowerModelsDistribution.SOCConicUBFForm}

PowerModelsDistribution.LPfullUBFPowerModel = GenericPowerModel{PowerModelsDistribution.LPfullUBFForm}
PowerModelsDistribution.LPdiagUBFPowerModel = GenericPowerModel{PowerModelsDistribution.LPdiagUBFForm}
PowerModelsDistribution.LPLinUBFPowerModel = GenericPowerModel{PowerModelsDistribution.LPLinUBFForm}
```

## Union Types

To support both conic and quadratically-constrained formulation variants for the unbalanced branch flow model, the union type `AbstractUBFForm` is defined. These formulations extend `AbstractBFForm` and are therefore also `AbstractWForms` (as defined in PowerModels proper).

```julia
AbstractUBFForm = Union{AbstractNLPUBFForm, AbstractConicUBFForm}
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
