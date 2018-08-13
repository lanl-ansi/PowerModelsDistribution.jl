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
ThreePhasePowerModels.AbstractNLPUBFForm <: PowerModels.AbstractBFQPForm
ThreePhasePowerModels.AbstractConicUBFForm <: PowerModels.AbstractBFConicForm
ThreePhasePowerModels.AbstractLPUBFForm <: ThreePhasePowerModels.AbstractNLPUBFForm
```

From there, different forms are possible:
```julia
#Bus injection models:
PowerModels.StandardACPForm <: PowerModels.AbstractACPForm
PowerModels.StandardDCPForm <: PowerModels.AbstractDCPForm
PowerModels.SOCWRForm <: PowerModels.AbstractWRForm

#Branch flow models:
ThreePhasePowerModels.SDPUBFForm <: ThreePhasePowerModels.AbstractConicUBFForm
ThreePhasePowerModels.SOCNLPUBFForm <: ThreePhasePowerModels.AbstractNLPUBFForm
ThreePhasePowerModels.SOCConicUBFForm <: ThreePhasePowerModels.AbstractConicUBFForm

ThreePhasePowerModels.LPLinUBFForm <: PowerModels.AbstractBFForm
ThreePhasePowerModels.LPfullUBFForm <: ThreePhasePowerModels.AbstractLPUBFForm
ThreePhasePowerModels.LPdiagUBFForm <: ThreePhasePowerModels.AbstractLPUBFForm
```

## Power Models
Each of these forms can be used as the type parameter for a PowerModel:
```julia
PowerModels.ACPPowerModel = GenericPowerModel{PowerModels.StandardACPForm}
PowerModels.DCPPowerModel = GenericPowerModel{PowerModels.StandardDCPForm}

PowerModels.SOCWRPowerModel = GenericPowerModel{PowerModels.SOCWRForm}

ThreePhasePowerModels.SDPUBFPowerModel = GenericPowerModel{ThreePhasePowerModels.SDPUBFForm}
ThreePhasePowerModels.SOCNLPUBFPowerModel = GenericPowerModel{ThreePhasePowerModels.SOCNLPUBFForm}
ThreePhasePowerModels.SOCConicUBFPowerModel = GenericPowerModel{ThreePhasePowerModels.SOCConicUBFForm}

ThreePhasePowerModels.LPfullUBFPowerModel = GenericPowerModel{ThreePhasePowerModels.LPfullUBFForm}
ThreePhasePowerModels.LPdiagUBFPowerModel = GenericPowerModel{ThreePhasePowerModels.LPdiagUBFForm}
ThreePhasePowerModels.LPLinUBFPowerModel = GenericPowerModel{ThreePhasePowerModels.LPLinUBFForm}
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
