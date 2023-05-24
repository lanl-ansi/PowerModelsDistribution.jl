# Engineering Data Model

This document describes the [`ENGINEERING`](@ref ENGINEERING) data model type in PowerModelsDistribution, which is transformed at runtime, or at the user's direction into a [`MATHEMATICAL`](@ref MATHEMATICAL) data model for optimization.

In this document,

- `nphases` refers to the number of non-neutral, non-ground active phases connected to a component,
- `nconductors` refers to all active conductors connected to a component, _i.e._ `length(connections)`, and
- `nwindings` refers to the number of windings of a transformer.

The data structure is in the following format

```julia
Dict{String,Any}(
    "data_model" => ENGINEERING,
    "component_type" => Dict{String,Dict{String,Any}}(
        id => Dict{String,Any}(
            "parameter" => value,
            ...
        ),
        ...
    ),
    ...
)
```

Valid component types are those that are documented in the sections below. Each component object is identified by an `id`, which must be a string (`id <: String`), but `id` does not appear inside of the component dictionary, and only appears as keys to the component dictionaries under each component type. Note that this requirement is so that data structures will be JSON serializable.

Each edge or node component (_i.e._ all those that are not data objects or buses), is expected to have `status` fields to specify whether the component is active or disabled, `bus` or `f_bus` and `t_bus`, to specify the buses that are connected to the component, and `connections` or `f_connections` and `t_connections`, to specify the terminals of the buses that are actively connected in an ordered list. **NOTE**: `terminals`, `connections`, `f_connections`, and `t_connections`, must be type `Vector{Int}`.

Parameter values on components are expected to be specified in SI units by default (where applicable) in the engineering data model. Relevant expected units are noted in the sections below. It is possible for the user to select universal scalar factors for power and voltages. For example, if `power_scalar_factor` and `voltage_scalar_factor` are their default values given below, where units is listed as watt or var, real units will be kW and kvar. Where units are listed as volt, real units will be kV (multiplied by `vm_nom`, where that value exists).

The Used column describes the situtations where certain parameters are used. "always" indicates those values are used in all contexts, `opf`, `mld`, or any other problem name abbreviation indicate they are used in particular for those problems. "solution" indicates that those parameters are outputs from the solvers. "multinetwork" indicates these values are only used to build multinetwork problems.

Those parameters that have a default may be omitted by the user from the data model, they will be populated by the specified default values.

Components that support "codes", such as lines, switches, and transformers, behave such that any property on said object that conflicts with a value in the code will override the value given in the code object. This is noted on each object where this is relevant.

## Root-Level Properties

At the root level of the data structure, the following fields can be found.

| Name         | Default       | Type                 | Used   | Description                                                                                                         |
| ------------ | ------------- | -------------------- | ------ | ------------------------------------------------------------------------------------------------------------------- |
| `name`       |               | `String`             |        | Case name                                                                                                           |
| `data_model` | `ENGINEERING` | `DataModel`          | always | `ENGINEERING`, `MATHEMATICAL`, or `DSS`. Type of the data model (this document describes `data_model==ENGINEERING`) |
| `settings`   | `Dict()`      | `Dict{String,<:Any}` | always | Base settings for the data model, see Settings section below for details                                            |

## Settings (`settings`)

At the root-level of the data model a `settings` dictionary object is expected, containing the following fields.

| Name                   | Default | Type                | Units | Used   | Description                                                                  |
| ---------------------- | ------- | ------------------- | ----- | ------ | ---------------------------------------------------------------------------- |
| `voltage_scale_factor` | `1e3`   | `Real`              |       | always | Scalar multiplier for voltage values                                         |
| `power_scale_factor`   | `1e3`   | `Real`              |       | always | Scalar multiplier for power values                                           |
| `vbases_default`       |         | `Dict{String,Real}` |       | always | Instruction to set the vbase at a number of buses for non-dimensionalization |
| `sbase_default`        |         | `Real`              |       | always | Instruction to set the power base for non-dimensionalization                 |
| `base_frequency`       | `60.0`  | `Real`              | Hz    | always | Frequency base, _i.e._ the base frequency of the whole circuit               |

The parameters `voltage_scale_factor` and `power_scale_factor`determine the base
for all voltage and power parameters in this data model. For example,

- `voltage_scale_factor=1E3` and `vm_nom=4.0`: `vm_nom` is `4.0 kV`/`4.0E3 V`,
- `power_scale_factor=1E6` and `pd_nom=2.0`: `pd_nom` is `2.0 MW`/`2.0E6 W`,
- `power_scale_factor=1E6` and `qd_nom=5.0`: `qd_nom` is `5.0 MVAr`/`5.0E6 VAr`,

where the mentioned fields `vm_nom`, `pd_nom` and `qd_nom` are sample voltage and power variables which are defined later.

On the other hand,`vbase_default` and `sbase_default` provide default values for a 'per unit' conversion; these do not affect the interpretation of the parameters in this model, like the scale factors do. Note that `vbase_default` is a `Dict{Any,Real}`, with pairs of bus ids and voltage magnitude levels, since in per unit conversion, the voltage base can change from bus to bus. The power base is the same everywhere, and therefore `sbase_default` has a single value.

## Buses (`bus`)

The data model below allows us to include buses of arbitrary many terminals (_i.e._, more than the usual four). This would be useful for

- underground lines with multiple neutrals which are not joined at every bus;
- distribution lines that carry several conventional lines in parallel (see for example the quad circuits in NEVTestCase).

| Name          | Default     | Type                  | Units  | Used         | Description                                                                          |
| ------------- | ----------- | --------------------- | ------ | ------------ | ------------------------------------------------------------------------------------ |
| `terminals`   | `[1,2,3,4]` | `Vector{Int}`         |        | always       | Terminals for which the bus has active connections                                   |
| `vm_lb`       |             | `Vector{Real}`        | volt   | opf          | Minimum conductor-to-ground voltage magnitude, `size=nphases`                        |
| `vm_ub`       |             | `Vector{Real}`        | volt   | opf          | Maximum conductor-to-ground voltage magnitude, `size=nphases`                        |
| `vm_pair_ub`  |             | `Vector{Tuple}`       |        | opf          | _e.g._ `[(1,2,210)]` means \|U1-U2\|>210                                             |
| `vm_pair_lb`  |             | `Vector{Tuple}`       |        | opf          | _e.g._ `[(1,2,230)]` means \|U1-U2\|<230                                             |
| `grounded`    | `[]`        | `Vector{Int}`         |        | always       | List of terminals which are grounded                                                 |
| `rg`          | `[]`        | `Vector{Real}`        |        | always       | Resistance of each defined grounding, `size=length(grounded)`                        |
| `xg`          | `[]`        | `Vector{Real}`        |        | always       | Reactance of each defined grounding, `size=length(grounded)`                         |
| `vm`          |             | `Vector{Real}`        | volt   | always       | Voltage magnitude at bus. If set, voltage magnitude at bus is fixed                  |
| `va`          |             | `Vector{Real}`        | degree | always       | Voltage angle at bus. If set, voltage angle at bus is fixed                          |
| `status`      | `ENABLED`   | `Status`              |        | always       | `ENABLED` or `DISABLED`. Indicates if component is enabled or disabled, respectively |
| `time_series` |             | `Dict{String,String}` |        | multinetwork | Dictionary containing time series parameters.                                        |

Each terminal `c` of the bus has an associated complex voltage phasor `v[c]`. There are two types of voltage magnitude bounds. The first type bounds the voltage magnitude of each `v[c]` individually,

- `lb <= |v[c]| <= ub`

However, especially in four-wire networks, bounds are more naturally imposed on the difference of two terminal voltages instead, e.g. for terminals `c` and `d`,

- `lb <= |v[c]-v[d]| <= ub`

This is why we introduce the fields `vm_pair_lb` and `vm_pair_ub`, which define bounds for pairs of terminals,

- $\forall$ `(c,d,lb)` $\in$ `vm_pair_lb`: `|v[c]-v[d]| >= lb`
- $\forall$ `(c,d,ub)` $\in$ `vm_pair_ub`: `|v[c]-v[d]| <= ub`

Finally, we give an example of how grounding impedances should be entered. If terminal `4` is grounded through an impedance `Z=1+j2`, we write

- `grounded=[4]`, `rg=[1]`, `xg=[2]`

### Special Case: three-phase bus

For three-phase buses, instead of specifying bounds explicitly for each pair of windings, often we want to specify 'phase-to-phase', 'phase-to-neutral' and 'neutral-to-ground' bounds. This can be done conveniently with a number of additional fields. First, `phases` is a list of the phase terminals, and `neutral` designates a single terminal to be the neutral.

- The bounds `vm_pn_lb` and `vm_pn_ub` specify the same lower and upper bound for the magnitude of the difference of each phase terminal and the neutral.
- The bounds `vm_pp_lb` and `vm_pp_ub` specify the same lower and upper bound for the magnitude of the difference of all phase terminals.
- `vm_ng_ub` specifies an upper bound for the neutral terminal, the lower bound is typically zero.

If all of these are specified, these bounds also imply valid bounds for the individual voltage magnitudes,

- $\forall$ `c` $\in$ `phases`: `vm_pn_lb - vm_ng_ub <= |v[c]| <= vm_pn_ub + vm_ng_ub`
- `0 <= |v[neutral]|<= vm_ng_ub`

Instead of defining the bounds directly, they can be specified through an associated voltage zone.

| Name       | Default | Type          | Units | Used   | Description                                                   |
| ---------- | ------- | ------------- | ----- | ------ | ------------------------------------------------------------- |
| `phases`   |         | `Vector{Int}` |       | always | Identifies the terminal that represents the neutral conductor |
| `neutral`  |         | `Int`         |       | always | Identifies the terminal that represents the neutral conductor |
| `vm_pn_lb` |         | `Real`        |       | opf    | Minimum phase-to-neutral voltage magnitude for all phases     |
| `vm_pn_ub` |         | `Real`        |       | opf    | Maximum phase-to-neutral voltage magnitude for all phases     |
| `vm_pp_lb` |         | `Real`        |       | opf    | Minimum phase-to-phase voltage magnitude for all phases       |
| `vm_pp_ub` |         | `Real`        |       | opf    | Maximum phase-to-phase voltage magnitude for all phases       |
| `vm_ng_ub` |         | `Real`        |       | opf    | Maximum neutral-to-ground voltage magnitude                   |

## Edge Objects

These objects represent edges on the power grid and therefore require `f_bus` and `t_bus` (or `buses` in the case of transformers), and `f_connections` and `t_connections` (or `connections` in the case of transformers).

### Lines (`line`)

This is a pi-model branch. When a `linecode` is given, and any of `rs`, `xs`, `b_fr`, `b_to`, `g_fr` or `g_to` are specified, any of those overwrite the values on the linecode.

| Name            | Default                           | Type           | Units            | Used   | Description                                                                          |
| --------------- | --------------------------------- | -------------- | ---------------- | ------ | ------------------------------------------------------------------------------------ |
| `f_bus`         |                                   | `String`       |                  | always | id of from-side bus connection                                                       |
| `t_bus`         |                                   | `String`       |                  | always | id of to-side bus connection                                                         |
| `f_connections` |                                   | `Vector{Int}`  |                  | always | Indicates for each conductor, to which terminal of the `f_bus` it connects           |
| `t_connections` |                                   | `Vector{Int}`  |                  | always | Indicates for each conductor, to which terminal of the `t_bus` it connects           |
| `linecode`      |                                   | `String`       |                  | always | id of an associated linecode                                                         |
| `rs`            |                                   | `Matrix{Real}` | ohm/meter        | always | Series resistance matrix, `size=(nconductors,nconductors)`                           |
| `xs`            |                                   | `Matrix{Real}` | ohm/meter        | always | Series reactance matrix, `size=(nconductors,nconductors)`                            |
| `g_fr`          | `zeros(nconductors, nconductors)` | `Matrix{Real}` | siemens/meter/Hz | always | From-side conductance, `size=(nconductors,nconductors)`                              |
| `b_fr`          | `zeros(nconductors, nconductors)` | `Matrix{Real}` | siemens/meter/Hz | always | From-side susceptance, `size=(nconductors,nconductors)`                              |
| `g_to`          | `zeros(nconductors, nconductors)` | `Matrix{Real}` | siemens/meter/Hz | always | To-side conductance, `size=(nconductors,nconductors)`                                |
| `b_to`          | `zeros(nconductors, nconductors)` | `Matrix{Real}` | siemens/meter/Hz | always | To-side susceptance, `size=(nconductors,nconductors)`                                |
| `length`        | `1.0`                             | `Real`         | meter            | always | Length of the line                                                                   |
| `cm_ub`         |                                   | `Vector{Real}` | amp              | opf    | Symmetrically applicable current rating, `size=nconductors`                          |
| `sm_ub`         |                                   | `Vector{Real}` | watt             | opf    | Symmetrically applicable power rating, `size=nconductors`                            |
| `vad_lb`        |                                   | `Vector{Real}` | degree           | opf    | Voltage angle difference lower bound                                                 |
| `vad_ub`        |                                   | `Vector{Real}` | degree           | opf    | Voltage angle difference upper bound                                                 |
| `status`        | `ENABLED`                         | `Status`       |                  | always | `ENABLED` or `DISABLED`. Indicates if component is enabled or disabled, respectively |

### Transformers (`transformer`)

These are n-winding (`nwinding`), n-phase (`nphase`), lossy transformers. Note that most properties are now Vectors (or Vectors of Vectors), indexed over the windings.

| Name             | Default                              | Type                   | Units       | Used   | Description                                                                                                                                                    |
| ---------------- | ------------------------------------ | ---------------------- | ----------- | ------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `bus`            |                                      | `Vector{String}`       |             | always | List of bus for each winding, `size=nwindings`                                                                                                                 |
| `connections`    |                                      | `Vector{Vector{Int}}`  |             | always | List of connection for each winding, `size=((nconductors),nwindings)`                                                                                          |
| `configurations` | `fill(WYE, nwindings)`               | `Vector{ConnConfig}`   |             | always | `WYE` or `DELTA`. List of configuration for each winding, `size=nwindings`                                                                                     |
| `xfmrcode`       |                                      | `String`               |             | always | id of                                                                                                                                                          |
| `xsc`            | `zeros(nwindings*(nwindings-1)/2)`   | `Vector{Real}`         | `sm_nom[1]` | always | List of short-circuit reactances between each pair of windings, relative to the VA rating of the first winding; enter as a list of the upper-triangle elements |
| `rw`             | `zeros(nwindings)`                   | `Vector{Real}`         | `sm_nom[1]` | always | Active power lost due to resistance of each winding, relative to the VA rating of each winding winding                                                         |
| `cmag`           | `0.0`                                | `Real`                 | `sm_nom[1]` | always | Total no-load reactive power drawn by the transformer, relative to VA rating of the first winding (magnetizing current)                                        |
| `noloadloss`     | `0.0`                                | `Real`                 | `sm_nom[1]` | always | Total no-load active power drawn by the transformer, relative to VA rating of the first winding                                                                |
| `tm_nom`         | `ones(nwindings)`                    | `Vector{Real}`         |             | always | Nominal tap ratio for the transformer, `size=nwindings` (multiplier)                                                                                           |
| `tm_ub`          |                                      | `Vector{Vector{Real}}` |             | opf    | Maximum tap ratio for each winding and phase, `size=((nphases),nwindings)` (base=`tm_nom`)                                                                     |
| `tm_lb`          |                                      | `Vector{Vector{Real}}` |             | opf    | Minimum tap ratio for for each winding and phase, `size=((nphases),nwindings)` (base=`tm_nom`)                                                                 |
| `tm_set`         | `fill(fill(1.0,nphases),nwindings)`  | `Vector{Vector{Real}}` |             | always | Set tap ratio for each winding and phase, `size=((nphases),nwindings)` (base=`tm_nom`)                                                                         |
| `tm_fix`         | `fill(fill(true,nphases),nwindings)` | `Vector{Vector{Bool}}` |             | oltc   | Indicates for each winding and phase whether the tap ratio is fixed, `size=((nphases),nwindings)`                                                              |
| `polarity`       | `fill(1,nwindings)`                  | `Vector{Int}`          |             | always |                                                                                                                                                                |
| `vm_nom`         |                                      | `Vector{Real}`         | volt        | always |                                                                                                                                                                |
| `sm_nom`         |                                      | `Vector{Real}`         | watt        | always |                                                                                                                                                                |
| `sm_ub`          |                                      | `Real`                 | watt        | opf    | Rating for the total apparent power magnitude at each winding                                                                                                  |
| `status`         | `ENABLED`                            | `Status`               |             | always | `ENABLED` or `DISABLED`. Indicates if component is enabled or disabled, respectively                                                                           |

#### Asymmetric, Lossless, Two-Winding (AL2W) Transformers (`transformer`)

Special case of the Generic transformer, which is still a `transformer` object, but has a simplified method for its definition. These are transformers are asymmetric (A), lossless (L) and two-winding (2W). Asymmetric refers to the fact that the secondary is always has a `WYE` configuration, whilst the primary can be `DELTA`. The table below indicates alternate, more simple ways to specify the special case of an AL2W Transformer. `xsc` and `rw` cannot be specified for an AL2W transformer, because it is lossless. To use this definition format, all of `f_bus`, `t_bus`, `f_connections`, `t_connections`, and `configuration` must be used, and none of `buses`, `connections`, `configurations` may be used. `xfmrcode` is ignored for this component.

| Name            | Default              | Type           | Units | Used   | Description                                                                                                 |
| --------------- | -------------------- | -------------- | ----- | ------ | ----------------------------------------------------------------------------------------------------------- |
| `f_bus`         |                      | `String`       |       | always | Alternative way to specify `buses`, requires both `f_bus` and `t_bus`                                       |
| `t_bus`         |                      | `String`       |       | always | Alternative way to specify `buses`, requires both `f_bus` and `t_bus`                                       |
| `f_connections` |                      | `Vector{Int}`  |       | always | Alternative way to specify `connections`, requires both `f_connections` and `t_connections`, `size=nphases` |
| `t_connections` |                      | `Vector{Int}`  |       | always | Alternative way to specify `connections`, requires both `f_connections` and `t_connections`, `size=nphases` |
| `configuration` | `WYE`                | `ConnConfig`   |       | always | `WYE` or `DELTA`. Alternative way to specify the from-side configuration, to-side is always `WYE`           |
| `tm_nom`        | `1.0`                | `Real`         |       | always | Nominal tap ratio for the transformer (multiplier)                                                          |
| `tm_ub`         |                      | `Vector{Real}` |       | opf    | Maximum tap ratio for each phase (base=`tm_nom`), `size=nphases`                                            |
| `tm_lb`         |                      | `Vector{Real}` |       | opf    | Minimum tap ratio for each phase (base=`tm_nom`), `size=nphases`                                            |
| `tm_set`        | `fill(1.0,nphases)`  | `Vector{Real}` |       | always | Set tap ratio for each phase (base=`tm_nom`), `size=nphases`                                                |
| `tm_fix`        | `fill(true,nphases)` | `Vector{Bool}` |       | oltc   | Indicates for each phase whether the tap ratio is fixed, `size=nphases`                                     |
| `sm_ub`         |                      | `Real`         |       | opf    | Rating for the total apparent power magnitude at each winding                                               |

#### Transformers with voltage regulator control (`controls`)

Special case of the Generic transformer, which is part of the `transformer` object, and emulates a standard utility voltage regulator. The taps of these transformers can be controlled by modelling a line drop compensator.

| Name      | Default | Type                   | Units | Used | Description                                                                                                                                                                |
| --------- | ------- | ---------------------- | ----- | ---- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `vreg`    |         | `Vector{Vector{Real}}` | volt  | oltc | Voltage regulator reference, default value is `120.0` for the controlled winding, `0.0` for winding without regulator control, `size=((nphases),nwindings)`                |
| `band`    |         | `Vector{Vector{Real}}` | volt  | oltc | Voltage bandwidth, default value is `3.0` for the controlled winding, `0.0` for winding without regulator control, `size=((nphases),nwindings)`                            |
| `ptratio` |         | `Vector{Vector{Real}}` |       | oltc | Voltage ratio of the potential transformer, default value is `60.0` for the controlled winding, `0.0` for winding without regulator control, `size=((nphases),nwindings)`  |
| `ctprim`  |         | `Vector{Vector{Real}}` | amp   | oltc | Current transformer rating on primary side, default value is `300.0` for the controlled winding, `0.0` for winding without regulator control, `size=((nphases),nwindings)` |
| `r`       |         | `Vector{Vector{Real}}` | volt  | oltc | Resistance setting on line drop compensator, default value is `0.0` for both controlled winding and winding without regulator control, `size=((nphases),nwindings)`        |
| `x`       |         | `Vector{Vector{Real}}` | volt  | oltc | Reactance setting on line drop compensator, default value is `0.0` for both controlled winding and winding without regulator control, `size=((nphases),nwindings)`         |

### Switches (`switch`)

Switches without `rs`, `xs` or a linecode (conductance/susceptance not considered), defined the switch will be treated as lossless. If lossy parameters are defined, `switch` objects will be decomposed into virtual `branch` & `bus`, and an ideal `switch`.

| Name            | Default                  | Type           | Units | Used   | Description                                                                                      |
| --------------- | ------------------------ | -------------- | ----- | ------ | ------------------------------------------------------------------------------------------------ |
| `f_bus`         |                          | `String`       |       | always | id of from-side bus connection                                                                   |
| `t_bus`         |                          | `String`       |       | always | id of to-side bus connection                                                                     |
| `f_connections` |                          | `Vector{Int}`  |       | always | Indicates for each conductor, to which terminal of the `f_bus` it connects                       |
| `t_connections` |                          | `Vector{Int}`  |       | always | Indicates for each conductor, to which terminal of the `t_bus` it connects                       |
| `cm_ub`         |                          | `Vector{Real}` | amp   | opf    | Symmetrically applicable current rating                                                          |
| `sm_ub`         |                          | `Vector{Real}` | watt  | opf    | Symmetrically applicable power rating                                                            |
| `linecode`      |                          | `String`       |       | always | id of an associated linecode, does not take into account conductance/susceptance                 |
| `rs`            | `zeros(nphases,nphases)` | `Matrix{Real}` | ohm   | always | Series resistance matrix, `size=(nphases,nphases)`                                               |
| `xs`            | `zeros(nphases,nphases)` | `Matrix{Real}` | ohm   | always | Series reactance matrix, `size=(nphases,nphases)`                                                |
| `dispatchable`  | `NO`                     | `Dispatchable` |       |        | `NO` or `YES`, indicates whether switch state can be changed in a switching optimization problem |
| `state`         | `CLOSED`                 | `SwitchState`  |       | always | `CLOSED`: closed or `OPEN`: open, to indicate state of switch                                    |
| `status`        | `ENABLED`                | `Status`       |       | always | `ENABLED` or `DISABLED`. Indicates if component is enabled or disabled, respectively             |

## Node Objects

These are objects that have single bus connections. Every object will have at least `bus`, `connections`, and `status`.

### Shunts (`shunt`)

| Name           | Default   | Type               | Units   | Used         | Description                                                                                                               |
| -------------- | --------- | ------------------ | ------- | ------------ | ------------------------------------------------------------------------------------------------------------------------- |
| `bus`          |           | `String`           |         | always       | id of bus connection                                                                                                      |
| `connections`  |           | `Vector{Int}`      |         | always       | Ordered list of connected conductors, `size=nconductors`                                                                  |
| `gs`           |           | `Matrix{Real}`     | siemens | always       | Conductance, `size=(nconductors,nconductors)`                                                                             |
| `bs`           |           | `Matrix{Real}`     | siemens | always       | Susceptance, `size=(nconductors,nconductors)`                                                                             |
| `model`        | `GENERIC` | `ShuntModel`       |         |              | `GENERIC`, `CAPACITOR`, or `REACTOR`. Indicates the type of shunt which may be necessary for transient stability analysis |
| `dispatchable` | `NO`      | `Dispatchable`     |         | mld          | `NO` or `YES`, indicates whether a shunt can be shed                                                                      |
| `status`       | `ENABLED` | `Status`           |         | always       | `ENABLED` or `DISABLED`. Indicates if component is enabled or disabled, respectively                                      |
| `time_series`  |           | `Dict{String,Any}` |         | multinetwork | Dictionary containing time series parameters.                                                                             |

#### Shunts with capacitor control (`controls`)

Special case of the shunt capacitors, which is part of the `shunt` object, and emulates a typical utility capacitor control (CapControl) by sending switching messages.

| Name           | Default | Type             | Units | Used | Description                                                                                                                                                                                       |
| -------------- | ------- | ---------------- | ----- | ---- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `type`         |         | `Vector{String}` |       | capc | Control type, default is `current` for controlled phase, ` ` for uncontrolled phase,`size=1`for`kvar`type, otherwise`size=(nphases)`                                                              |
| `element`      |         | `String`         |       | capc | `source_id` of element (typically line or transformer) to which CapControl is connected                                                                                                           |
| `terminal`     |         | `Vector{Int}`    |       | capc | Number of the terminal of circuit element to which CapControl is connected, default is `1` for controlled phase, `0` for uncontrolled phase, `size=1` for `kvar` type, otherwise `size=(nphases)` |
| `onsetting`    |         | `Vector{Real}`   |       | capc | Value at which the CapControl switches the capacitor on, default is `300.0` for controlled phase, `0.0` for uncontrolled phase, `size=1` for `kvar` type, otherwise `size=(nphases)`              |
| `offsetting`   |         | `Vector{Real}`   |       | capc | Value at which the CapControl switches the capacitor off, default is `200.0` for controlled phase, `0.0` for uncontrolled phase, `size=1` for `kvar` type, otherwise `size=(nphases)`             |
| `voltoverride` |         | `Vector{Bool}`   |       | capc | Indicate whether voltage over ride is enabled, default is `false` for both controlled and uncontrolled phases, `size=1` for `kvar` type, otherwise `size=(nphases)`                               |
| `ptratio`      |         | `Vector{Real}`   |       | capc | Ratio of potential transformer, default is `60.0` for controlled phase, `0.0` for uncontrolled phase, `size=(nphases)`                                                                            |
| `ctratio`      |         | `Vector{Real}`   |       | capc | Ratio of current transformer, default is `60.0` for controlled phase, `0.0` for uncontrolled phase, `size=(nphases)`                                                                              |
| `vmin`         |         | `Vector{Real}`   | volt  | capc | Minimum voltage below which CapControl switches the capacitor on, default is `115.0` for controlled phase, `0.0` for uncontrolled phase, `size=1` for `kvar` type, otherwise `size=(nphases)`     |
| `vmax`         |         | `Vector{Real}`   | volt  | capc | Maximum voltage above which CapControl switches the capacitor off, default is `126.0` for controlled phase, `0.0` for uncontrolled phase, `size=1` for `kvar` type, otherwise `size=(nphases)`    |

### Loads (`load`)

| Name            | Default   | Type                  | Units | Used           | Description                                                                                        |
| --------------- | --------- | --------------------- | ----- | -------------- | -------------------------------------------------------------------------------------------------- |
| `bus`           |           | `String`              |       | always         | id of bus connection                                                                               |
| `connections`   |           | `Vector{Int}`         |       | always         | Ordered list of connected conductors, `size=nconductors`                                           |
| `configuration` | `WYE`     | `ConnConfig`          |       | always         | `WYE` or `DELTA`. If `WYE`, `connections[end]=neutral`                                             |
| `model`         | `POWER`   | `LoadModel`           |       | always         | `POWER`, `IMPEDANCE`, `CURRENT`, `EXPONENTIAL`, or `ZIP`. Indicates the type of voltage-dependency |
| `pd_nom`        |           | `Vector{Real}`        | watt  | always         | Nominal active load, with respect to `vm_nom`, `size=nphases`                                      |
| `qd_nom`        |           | `Vector{Real}`        | var   | always         | Nominal reactive load, with respect to `vm_nom`, `size=nphases`                                    |
| `vm_nom`        |           | `Real`                | volt  | `model!=POWER` | Nominal voltage (multiplier)                                                                       |
| `dispatchable`  | `NO`      | `Dispatchable`        |       | mld            | `NO` or `YES`, indicates whether a load can be shed                                                |
| `status`        | `ENABLED` | `Status`              |       | always         | `ENABLED` or `DISABLED`. Indicates if component is enabled or disabled, respectively               |
| `time_series`   |           | `Dict{String,String}` |       | multinetwork   | Dictionary containing time series parameters.                                                      |

Multi-phase loads define a number of individual loads connected between two terminals each. How they are connected, is defined both by `configuration` and `connections`. The table below indicates the value of `configuration` and lengths of the other properties for a consistent definition,

| `configuration` | `connections` | `pd_nom \| qd_nom \| pd_exp` |
| --------------- | ------------- | ---------------------------- |
| `DELTA`         | `2`           | `1`                          |
| `DELTA`         | `3`           | `3`                          |
| `WYE`           | `2`           | `1`                          |
| `WYE`           | `3`           | `2`                          |
| `WYE`           | `N`           | `N-1`                        |

Note that for delta loads, only 2 and 3 connections are allowed. Each individual load `i` is connected between two terminals, exposed to a voltage magnitude `v[i]`, which leads to a consumption `pd[i]+j*qd[i]`. The `model` then defines the relationship between these quantities,

| model       | `pd[i]/pd_nom[i]=` | `qd[i]/qd_nom[i]=` |
| ----------- | ------------------ | ------------------ |
| `POWER`     | `1`                | `1`                |
| `CURRENT`   | `(v[i]/vm_nom)`    | `(v[i]/vm_nom)`    |
| `IMPEDANCE` | `(v[i]/vm_nom)^2`  | `(v[i]/vm_nom)^2`  |

Two more model types are supported, which need additional fields and are defined below.

#### `model == EXPONENTIAL`

- `(pd[i]/pd_nom[i]) = (v[i]/vm_nom)^pd_exp[i]`
- `(qd[i]/qd_nom[i]) = (v[i]/vm_nom)^qd_exp[i]`

| Name     | Default | Type   | Units | Used                 | Description |
| -------- | ------- | ------ | ----- | -------------------- | ----------- |
| `pd_exp` |         | `Real` |       | `model==EXPONENTIAL` |             |
| `qd_exp` |         | `Real` |       | `model==EXPONENTIAL` |             |

#### `model == ZIP`

ZIP load models are split into `IMPEDANCE`, `CURRENT`, `POWER` models.

- `(pd[i]/pd_nom) = pd_cz[i]*(v[i]/vm_nom)^2 + pd_ci[i]*(v[i]/vm_nom) + pd_cp[i]`
- `(qd[i]/qd_nom) = qd_cz[i]*(v[i]/vm_nom)^2 + qd_ci[i]*(v[i]/vm_nom) + qd_cp[i]`

| Name     | Default | Type           | Units | Used         | Description                                                                                                                                                                                                             |
| -------- | ------- | ------         | ----- | ------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `zipv`   |         | `Vector{Real}` |       | `model==ZIP` | First 3 are ZIP weighting factors for active power (`pd_cz,pd_ci,pd_cp`), next 3 are ZIP weighting factors for reactive power (`qd_cz,qd_ci,qd_cp`), last 1 is cut-off voltage in p.u. of base kV; load is 0 below this cut-off |

### Generators (`generator`)

| Name            | Default              | Type                  | Units | Used                        | Description                                                                          |
| --------------- | -------------------- | --------------------- | ----- | --------------------------- | ------------------------------------------------------------------------------------ |
| `bus`           |                      | `String`              |       | always                      | id of bus connection                                                                 |
| `connections`   |                      | `Vector{Int}`         |       | always                      | Ordered list of connected conductors, `size=nconductors`                             |
| `configuration` | `WYE`                | `ConnConfig`          |       | always                      | `WYE` or `DELTA`. If `WYE`, `connections[end]=neutral`                               |
| `vg`            |                      | `Vector{Real}`        | volt  | `control_mode==ISOCHRONOUS` | Voltage magnitude setpoint                                                           |
| `pg_lb`         | `zeros(nphases)`     | `Vector{Real}`        | watt  | opf                         | Lower bound on active power generation per phase, `size=nphases`                     |
| `pg_ub`         | `fill(Inf, nphases)` | `Vector{Real}`        | watt  | opf                         | Upper bound on active power generation per phase, `size=nphases`                     |
| `qg_lb`         | `-pg_ub`             | `Vector{Real}`        | var   | opf                         | Lower bound on reactive power generation per phase, `size=nphases`                   |
| `qg_ub`         | `pg_ub`              | `Vector{Real}`        | var   | opf                         | Upper bound on reactive power generation per phase, `size=nphases`                   |
| `pg`            |                      | `Vector{Real}`        | watt  | solution                    | Present active power generation per phase, `size=nphases`                            |
| `qg`            |                      | `Vector{Real}`        | var   | solution                    | Present reactive power generation per phase, `size=nphases`                          |
| `control_mode`  | `FREQUENCYDROOP`     | `ControlMode`         |       |                             | `FREQUENCYDROOP` or `ISOCHRONOUS`                                                    |
| `status`        | `ENABLED`            | `Status`              |       | always                      | `ENABLED` or `DISABLED`. Indicates if component is enabled or disabled, respectively |
| `time_series`   |                      | `Dict{String,String}` |       | multinetwork                | Dictionary containing time series parameters.                                        |

#### `generator` Cost Model

The generator cost model is currently specified by the following fields.

| Name                 | Default           | Type           | Units | Used | Description                                               |
| -------------------- | ----------------- | -------------- | ----- | ---- | --------------------------------------------------------- |
| `cost_pg_model`      | `2`               | `Int`          |       | opf  | Cost model type, `1` = piecewise-linear, `2` = polynomial |
| `cost_pg_parameters` | `[0.0, 1.0, 0.0]` | `Vector{Real}` | $/MVA | opf  | Cost model polynomial                                     |

### Photovoltaic Systems (`solar`)

| Name            | Default   | Type                  | Units | Used         | Description                                                                          |
| --------------- | --------- | --------------------- | ----- | ------------ | ------------------------------------------------------------------------------------ |
| `bus`           |           | `String`              |       | always       | id of bus connection                                                                 |
| `connections`   |           | `Vector{Int}`         |       | always       | Ordered list of connected conductors, `size=nconductors`                             |
| `configuration` | `WYE`     | `ConnConfig`          |       | always       | `WYE` or `DELTA`. If `WYE`, `connections[end]=neutral`                               |
| `pg_lb`         |           | `Vector{Real}`        | watt  | opf          | Lower bound on active power generation per phase, `size=nphases`                     |
| `pg_ub`         |           | `Vector{Real}`        | watt  | opf          | Upper bound on active power generation per phase, `size=nphases`                     |
| `qg_lb`         |           | `Vector{Real}`        | var   | opf          | Lower bound on reactive power generation per phase, `size=nphases`                   |
| `qg_ub`         |           | `Vector{Real}`        | var   | opf          | Upper bound on reactive power generation per phase, `size=nphases`                   |
| `pg`            |           | `Vector{Real}`        | watt  | solution     | Present active power generation per phase, `size=nphases`                            |
| `qg`            |           | `Vector{Real}`        | var   | solution     | Present reactive power generation per phase, `size=nphases`                          |
| `status`        | `ENABLED` | `Status`              |       | always       | `ENABLED` or `DISABLED`. Indicates if component is enabled or disabled, respectively |
| `time_series`   |           | `Dict{String,String}` |       | multinetwork | Dictionary containing time series parameters.                                        |

#### `solar` Cost Model

The cost model for a photovoltaic system currently matches that of generators.

| Name                 | Default           | Type           | Units | Used | Description                                               |
| -------------------- | ----------------- | -------------- | ----- | ---- | --------------------------------------------------------- |
| `cost_pg_model`      | `2`               | `Int`          |       | opf  | Cost model type, `1` = piecewise-linear, `2` = polynomial |
| `cost_pg_parameters` | `[0.0, 1.0, 0.0]` | `Vector{Real}` | $/MVA | opf  | Cost model polynomial                                     |

### Wind Turbine Systems (`wind`)

Wind turbine systems are most closely approximated by induction machines, also known as asynchronous machines. These are not currently supported, but there is plans to support them in the future.

### Storage (`storage`)

A storage object is a flexible component that can represent a variety of energy storage objects, like Li-ion batteries, hydrogen fuel cells, flywheels, etc.

- How to include the inverter model for this? Similar issue as for a PV generator

| Name                   | Default   | Type                  | Units   | Used         | Description                                                                          |
| ---------------------- | --------- | --------------------- | ------- | ------------ | ------------------------------------------------------------------------------------ |
| `bus`                  |           | `String`              |         | always       | id of bus connection                                                                 |
| `connections`          |           | `Vector{Int}`         |         | always       | Ordered list of connected conductors, `size=nconductors`                             |
| `configuration`        | `WYE`     | `ConnConfig`          |         | always       | `WYE` or `DELTA`. If `WYE`, `connections[end]=neutral`                               |
| `energy`               |           | `Real`                | watt-hr | always       | Stored energy                                                                        |
| `energy_ub`            |           | `Real`                |         | opf          | maximum energy rating                                                                |
| `charge_ub`            |           | `Real`                |         | opf          | maximum charge rating                                                                |
| `discharge_ub`         |           | `Real`                |         | opf          | maximum discharge rating                                                             |
| `sm_ub`                |           | `Real`                | watt    | opf          | Power rating,                                                                        |
| `cm_ub`                |           | `Real`                | amp     | opf          | Current rating,                                                                      |
| `charge_efficiency`    |           | `Real`                | percent | always       | charging efficiency (losses)                                                         |
| `discharge_efficiency` |           | `Real`                | percent | always       | disharging efficiency (losses)                                                       |
| `qs_ub`                |           | `Real`                |         | opf          | Maximum reactive power injection,                                                    |
| `qs_lb`                |           | `Real`                |         | opf          | Minimum reactive power injection,                                                    |
| `rs`                   |           | `Real`                | ohm     | always       | converter resistance                                                                 |
| `xs`                   |           | `Real`                | ohm     | always       | converter reactance                                                                  |
| `pex`                  |           | `Real`                |         | always       | Total active power standby exogenous flow (loss)                                     |
| `qex`                  |           | `Real`                |         | always       | Total reactive power standby exogenous flow (loss)                                   |
| `ps`                   |           | `Vector{Real}`        | watt    | solution     | Present active power injection                                                       |
| `qs`                   |           | `Vector{Real}`        | var     | solution     | Present reactive power injection                                                     |
| `status`               | `ENABLED` | `Status`              |         | always       | `ENABLED` or `DISABLED`. Indicates if component is enabled or disabled, respectively |
| `time_series`          |           | `Dict{String,String}` |         | multinetwork | Dictionary containing time series parameters.                                        |

### Voltage Sources (`voltage_source`)

A voltage source is a source of power at a set voltage magnitude and angle connected to a slack bus. If `rs` or `xs` are not specified, the voltage source is assumed to be lossless, otherwise virtual `branch` and `bus` will be created in the mathematical model to represent the internal losses of the voltage source.

| Name            | Default                          | Type                  | Units  | Used         | Description                                                                          |
| --------------- | -------------------------------- | --------------------- | ------ | ------------ | ------------------------------------------------------------------------------------ |
| `bus`           |                                  | `String`              |        | always       | id of bus connection                                                                 |
| `connections`   |                                  | `Vector{Int}`         |        | always       | Ordered list of connected conductors, `size=nconductors`                             |
| `configuration` | `WYE`                            | `ConnConfig`          |        | always       | `WYE` or `DELTA`. If `WYE`, `connections[end]=neutral`                               |
| `vm`            | `ones(nphases)`                  | `Vector{Real}`        | volt   | always       | Voltage magnitude set at slack bus, `size=nphases`                                   |
| `va`            | `zeros(nphases)`                 | `Real`                | degree | always       | Voltage angle offsets at slack bus, applies symmetrically to each phase angle        |
| `rs`            | `zeros(nconductors,nconductors)` | `Matrix{Real}`        | ohm    | always       | Internal series resistance of voltage source, `size=(nconductors,nconductors)`       |
| `xs`            | `zeros(nconductors,nconductors)` | `Matrix{Real}`        | ohm    | always       | Internal series reactance of voltage soure, `size=(nconductors,nconductors)`         |
| `status`        | `ENABLED`                        | `Status`              |        | always       | `ENABLED` or `DISABLED`. Indicates if component is enabled or disabled, respectively |
| `time_series`   |                                  | `Dict{String,String}` |        | multinetwork | Dictionary containing time series parameters.                                        |

## Data Objects (codes, time series, etc.)

These objects are referenced by node and edge objects, but are not part of the network themselves, only containing data.

### Linecodes (`linecode`)

Linecodes are easy ways to specify properties common to multiple lines.

| Name    | Default                          | Type           | Units            | Used   | Description                                             |
| ------- | -------------------------------- | -------------- | ---------------- | ------ | ------------------------------------------------------- |
| `rs`    |                                  | `Matrix{Real}` | ohm/meter        | always | Series resistance, `size=(nconductors,nconductors)`     |
| `xs`    |                                  | `Matrix{Real}` | ohm/meter        | always | Series reactance, `size=(nconductors,nconductors)`      |
| `g_fr`  | `zeros(nconductors,nconductors)` | `Matrix{Real}` | siemens/meter/Hz | always | From-side conductance, `size=(nconductors,nconductors)` |
| `b_fr`  | `zeros(nconductors,nconductors)` | `Matrix{Real}` | siemens/meter/Hz | always | From-side susceptance, `size=(nconductors,nconductors)` |
| `g_to`  | `zeros(nconductors,nconductors)` | `Matrix{Real}` | siemens/meter/Hz | always | To-side conductance, `size=(nconductors,nconductors)`   |
| `b_to`  | `zeros(nconductors,nconductors)` | `Matrix{Real}` | siemens/meter/Hz | always | To-side susceptance, `size=(nconductors,nconductors)`   |
| `cm_ub` | `fill(Inf,nconductors)`          | `Vector{Real}` | ampere           | opf    | maximum current per conductor, symmetrically applicable |
| `sm_ub` | `fill(Inf,nconductors)`          | `Vector{Real}` | watt             | opf    | maximum power per conductor, symmetrically applicable   |

### Transformer Codes (`xfmrcode`)

Transformer codes are easy ways to specify properties common to multiple transformers

| Name             | Default                                | Type                   | Units | Used   | Description                                                                                                                                     |
| ---------------- | -------------------------------------- | ---------------------- | ----- | ------ | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| `configurations` | `fill(WYE, nwindings)`                 | `Vector{ConnConfig}`   |       | always | `WYE` or `DELTA`. List of configuration for each winding, `size=nwindings`                                                                      |
| `xsc`            | `[0.0]`                                | `Vector{Real}`         | ohm   | always | List of short-circuit reactances between each pair of windings; enter as a list of the upper-triangle elements, `size=(nwindings == 2 ? 1 : 3)` |
| `rw`             | `zeros(nwindings)`                     | `Vector{Real}`         | ohm   | always | List of the winding resistance for each winding, `size=nwindings`                                                                               |
| `tm_nom`         | `ones(nwindings)`                      | `Vector{Real}`         |       | always | Nominal tap ratio for the transformer, `size=nwindings` (multiplier)                                                                            |
| `tm_ub`          |                                        | `Vector{Vector{Real}}` |       | opf    | Maximum tap ratio for each winding and phase, `size=((nphases), nwindings)` (base=`tm_nom`)                                                     |
| `tm_lb`          |                                        | `Vector{Vector{Real}}` |       | opf    | Minimum tap ratio for for each winding and phase, `size=((nphases), nwindings)` (base=`tm_nom`)                                                 |
| `tm_set`         | `fill(fill(1.0, nphases), nwindings)`  | `Vector{Vector{Real}}` |       | always | Set tap ratio for each winding and phase, `size=((nphases), nwindings)` (base=`tm_nom`)                                                         |
| `tm_fix`         | `fill(fill(true, nphases), nwindings)` | `Vector{Vector{Bool}}` |       | always | Indicates for each winding and phase whether the tap ratio is fixed, `size=((nphases), nwindings)`                                              |

### Time Series (`time_series`)

Time series objects are used to specify time series for _e.g._ load or generation forecasts.

Some parameters for components specified in this document can support a time series by inserting a referece to a `time_series` object into the `time_series` dictionary inside a component under the relevant parameter name. For example, for a `load`, if `pd_nom` is supposed to be a time series, the user would specify `"time_series" => Dict("pd_nom" => time_series_id)` where `time_series_id` is the `id` of an object in `time_series`, and has type `Any`.

| Name      | Default | Type                                 | Units | Used   | Description                                                                                                   |
| --------- | ------- | ------------------------------------ | ----- | ------ | ------------------------------------------------------------------------------------------------------------- |
| `time`    |         | `Union{Vector{Real},Vector{String}}` | hour  | always | Time points at which values are specified. If time is specified in String, units not required to be in hours. |
| `values`  |         | `Vector{Real}`                       |       | always | Multipers at each time step given in `time`                                                                   |
| `offset`  | `0`     | `Real`                               | hour  | always | Start time offset                                                                                             |
| `replace` | `true`  | `Bool`                               |       | always | Indicates to replace with data, instead of multiply. Will only work on non-Array data                         |

### Fuses (`fuse`)

Fuses can be defined on any terminal of any physical component

| Name                    | Default | Type                    | Units | Used | Description                                          |
| ----------------------- | ------- | ----------------------- | ----- | ---- | ---------------------------------------------------- |
| `component_type`        |         | `String`                |       |      |                                                      |
| `component_id`          |         | `String`                |       |      |                                                      |
| `terminals`             |         | `Vector{Int}`           |       |      |                                                      |
| `fuse_curve`            |         | `Array{Vector{Real},2}` |       |      | specifies the fuse blowing condition                 |
| `minimum_melting_curve` |         | `Array{Vector{Real},2}` |       |      | specifies the minimum melting conditions of the fuse |
