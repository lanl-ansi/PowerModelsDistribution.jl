# PowerModelsDistribution Enum Types

Within the PowerModelsDistribution Engineering Model we have included the use of Enums. Here we document the fields for which Enums are expected and the possible Enums available

## Data Model

Any place in PowerModelsDistribution that calls for specifying the `data_model`, either in function calls or the `"data_model"` field inside the data structure itself, will expect a `DataModel` type

```julia
@enum DataModel ENGINEERING MATHEMATICAL DSS
```

The `DSS` data model is an output of `parse_dss`, and is an untranslated raw parse of a DSS file. This Enum exists for use by `count_nodes`, where the method to count the number of active nodes is different between all three models

## Component statuses

All `"status"` fields in the `ENGINEERING` model expect a `Status` type:

```julia
@enum Status ENABLED DISABLED
```

## Connection Configurations

All `"configuration"` fields in the `ENGINEERING` model expect a `ConnConfig` type:

```julia
@enum ConnConfig WYE DELTA
```

## Load Models

For `load` objects, the `"model"` field expects a `LoadModel` type to specify the type of load model to use:

```julia
@enum LoadModel POWER CURRENT IMPEDANCE EXPONENTIAL ZIP
```

where `POWER` indicates constant power, `CURRENT` indicates constant current, `IMPEDANCE` indicates constant impedance, `EXPONENTIAL` indicates an exponential load model, and `ZIP` indicates a ZIP model

## Shunt Models

For `shunt` objects, the `"model"` field expects a `ShuntModel` type to specify the origin of the shunt object, which is important for transient analysis:

```julia
@enum ShuntModel GENERIC CAPACITOR REACTOR
```

## Switch States

For `switch` objects, the `"state"` field expects a `SwitchState` type to specify whether the switch is currently open or closed:

```julia
@enum SwitchState OPEN CLOSED
```

## Dispatchable Components

Some components can be dispatchable, _e.g._ if a switch is dispatchable that means it is free to open or close, but if not then it is fixed in place, or if a load is dispatchable it implies that it can be shed in a `run_mc_mld` problem:

```julia
@enum Dispatchable NO YES
```

## Generator Control Modes

For generator objects, the `"control_mode"` field expects a `ControlMode` type to specify whether the generator is operating in an isochronous mode (_i.e._ is frequency forming) or droop mode (_i.e._ is frequency following):

```julia
@enum ControlMode FREQUENCYDROOP ISOCHRONOUS
```
