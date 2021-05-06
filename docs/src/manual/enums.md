# PowerModelsDistribution Enum Types

Within the PowerModelsDistribution Engineering Model we have included the use of Enums. Here we document the fields for which Enums are expected and the possible Enums available

## Data Model

Any place in PowerModelsDistribution that calls for specifying the `data_model`, either in function calls or the `"data_model"` field inside the data structure itself, will expect a [`DataModel`](@ref DataModel) type.

The `DSS` data model is an output of `parse_dss`, and is an untranslated raw parse of a DSS file. This Enum exists for use by `count_nodes`, where the method to count the number of active nodes is different between all three models

## Component Status

All `"status"` fields in the [`ENGINEERING`](@ref ENGINEERING) model expect a [`Status`](@ref Status) type.

## Connection Configuration

All `"configuration"` fields in the [`ENGINEERING`](@ref ENGINEERING) model expect a [`ConnConfig`](@ref ConnConfig) type:

## Load Model

For `load` objects, the `"model"` field expects a [`LoadModel`](@ref LoadModel) type to specify the type of load model to use, where

- [`POWER`](@ref POWER) indicates constant power,
- [`CURRENT`](@ref CURRENT) indicates constant current,
- [`IMPEDANCE`](@ref IMPEDANCE) indicates constant impedance,
- [`EXPONENTIAL`](@ref EXPONENTIAL) indicates an exponential load model, and
- [`ZIP`](@ref ZIP) indicates a ZIP model

## Shunt Model

For `shunt` objects, the `"model"` field expects a [`ShuntModel`](@ref ShuntModel) type to specify the origin of the shunt object, which is important for transient analysis.

## Switch State

For `switch` objects, the `"state"` field expects a [`SwitchState`](@ref SwitchState) type to specify whether the switch is currently open or closed.

## Dispatchable Component

Some components can be [`Dispatchable`](@ref Dispatchable), _e.g._ if a switch is dispatchable that means it is free to open or close, but if not then it is fixed in place, or if a load is dispatchable it implies that it can be shed in a `run_mc_mld` problem.

## Generator Control Mode

For generator objects, the `"control_mode"` field expects a [`ControlMode`](@ref ControlMode) type to specify whether the generator is operating in an isochronous mode (_i.e._ is frequency forming) or droop mode (_i.e._ is frequency following).
