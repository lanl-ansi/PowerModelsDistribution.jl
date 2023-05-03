
# TODO repace with schemas
const _pmd_eng_object_names = String[
    "bus",
    "line",
    "switch",
    "transformer",
    "load",
    "shunt",
    "generator",
    "solar",
    "storage",
    "linecode",
    "xfmrcode",
    "time_series"
]


abstract type Settings <: EngDataObject end

"""
    Settings <: EngDataObject

At the root-level of the data model a `settings` dictionary object is expected, containing the following fields.

| Name                   | Default | Type               | Units | Used   | Description                                                                  |
|------------------------|---------|--------------------|-------|--------|------------------------------------------------------------------------------|
| `voltage_scale_factor` | `1e3`   | `Real`             |       | always | Scalar multiplier for voltage values                                         |
| `power_scale_factor`   | `1e3`   | `Real`             |       | always | Scalar multiplier for power values                                           |
| `vbases_default`       |         | `Dict{<:Any,Real}` |       | always | Instruction to set the vbase at a number of buses for non-dimensionalization |
| `sbase_default`        |         | `Real`             |       | always | Instruction to set the power base for non-dimensionalization                 |
| `base_frequency`       | `60.0`  | `Real`             | Hz    | always | Frequency base, _i.e._ the base frequency of the whole circuit               |

The parameters `voltage_scale_factor` and `power_scale_factor`determine the base
for all voltage and power parameters in this data model. For example,

- `voltage_scale_factor=1E3` and `vm_nom=4.0`: `vm_nom` is `4.0 kV`/`4.0E3 V`,
- `power_scale_factor=1E6` and `pd_nom=2.0`: `pd_nom` is `2.0 MW`/`2.0E6 W`,
- `power_scale_factor=1E6` and `qd_nom=5.0`: `qd_nom` is `5.0 MVAr`/`5.0E6 VAr`,

where the mentioned fields `vm_nom`, `pd_nom` and `qd_nom` are sample voltage and power variables which are defined later.

On the other hand,`vbase_default` and `sbase_default` provide default values for a 'per unit' conversion; these do not affect the interpretation of the parameters in this model, like the scale factors do. Note that `vbase_default` is a `Dict{Any,Real}`, with pairs of bus ids and voltage magnitude levels, since in per unit conversion, the voltage base can change from bus to bus. The power base is the same everywhere, and therefore `sbase_default` has a single value.
"""
Base.@kwdef mutable struct SettingsObj <: Settings
    voltage_scale_factor::Float64 = 1e3
    power_scale_factor::Float64 = 1e3
    sbase_default::Float64 = 1.0
    base_frequency::Float64 = 60.0
    vbases_default::Dict{String,Float64} = Dict{String,Float64}()
    time_elapsed::Float64 = 1.0
    dss::Union{Missing,DssOptions} = missing
end


abstract type Metadata <: EngDataObject end

"""
"""
Base.@kwdef mutable struct MetadataObj <: Metadata
    name::String = ""
    conductors::Vector{Int} = Int[]
    awaiting_ground::Dict{String,Vector{Vector{Int}}} = Dict{String,Vector{Vector{Int}}}()
end


abstract type EngBus <: EngNodeObject end


@doc raw"""
    EngBus <: EngNodeObject

The data model below allows us to include buses of arbitrary many terminals (_i.e._, more than the usual four). This would be useful for

- underground lines with multiple neutrals which are not joined at every bus;
- distribution lines that carry several conventional lines in parallel (see for example the quad circuits in NEVTestCase).

| Name          | Default     | Type                  | Units  | Used         | Description                                                                          |
| ------------- | ----------- | --------------------- | ------ | ------------ | ------------------------------------------------------------------------------------ |
| `terminals`   | `[1,2,3,4]` | `Vector{Int}`         |        | always       | Terminals for which the bus has active connections                                   |
| `vm_lb`       |             | `Vector{Real}`        | volt   | opf          | Minimum conductor-to-ground voltage magnitude, `size=nphases`                        |
| `vm_ub`       |             | `Vector{Real}`        | volt   | opf          | Maximum conductor-to-ground voltage magnitude, `size=nphases`                        |
| `vm_pair_ub`  |             | `Vector{Tuple}`       |        | opf          | _e.g._ `[(1,2,210)]` means \|U1-U2\| > 210                                           |
| `vm_pair_lb`  |             | `Vector{Tuple}`       |        | opf          | _e.g._ `[(1,2,230)]` means \|U1-U2\| < 230                                           |
| `grounded`    | `[]`        | `Vector{Int}`         |        | always       | List of terminals which are grounded                                                 |
| `rg`          | `[]`        | `Vector{Real}`        |        | always       | Resistance of each defined grounding, `size=length(grounded)`                        |
| `xg`          | `[]`        | `Vector{Real}`        |        | always       | Reactance of each defined grounding, `size=length(grounded)`                         |
| `vm`          |             | `Vector{Real}`        | volt   | always       | Voltage magnitude at bus. If set, voltage magnitude at bus is fixed                  |
| `va`          |             | `Vector{Real}`        | degree | always       | Voltage angle at bus. If set, voltage angle at bus is fixed                          |
| `status`      | `ENABLED`   | `Status`              |        | always       | `ENABLED` or `DISABLED`. Indicates if component is enabled or disabled, respectively |
| `time_series` |             | `Dict{String,String}` |        | multinetwork | Dictionary containing time series parameters.                                        |

Each terminal `c` of the bus has an associated complex voltage phasor `v[c]`. There are two types of voltage magnitude bounds. The first type bounds the voltage magnitude of each `v[c]` individually,

- $$lb \leq |v[c]| \leq ub$$

However, especially in four-wire networks, bounds are more naturally imposed on the difference of two terminal voltages instead, e.g. for terminals `c` and `d`,

- $$lb \leq |v[c]-v[d]| \leq ub$$

This is why we introduce the fields `vm_pair_lb` and `vm_pair_ub`, which define bounds for pairs of terminals,

- $$\forall (c,d,lb) \in vm_pair_lb: |v[c]-v[d]| \geq lb$$
- $$\forall (c,d,ub) \in vm_pair_ub: |v[c]-v[d]| \leq ub$$

Finally, we give an example of how grounding impedances should be entered. If terminal `4` is grounded through an impedance `Z=1+j2`, we write

- `grounded=[4]`, `rg=[1]`, `xg=[2]`
"""
Base.@kwdef mutable struct EngBusObj <: EngBus
    name::String
    terminals::Vector{Int} = Int[]
    vm_lb::Union{Missing,Vector{Float64}} = missing
    vm_ub::Union{Missing,Vector{Float64}} = missing
    vm_pair_ub::Union{Missing,Vector{NTuple{3,Float64}}} = missing
    vm_pair_lb::Union{Missing,Vector{NTuple{3,Float64}}} = missing
    grounded::Vector{Int} = Int[]
    rg::Vector{Float64} = Float64[]
    xg::Vector{Float64} = Float64[]
    vm::Union{Missing,Vector{Float64}} = missing
    va::Union{Missing,Vector{Float64}} = missing
    vbase::Union{Missing,Float64} = missing
    status::Status = ENABLED
    source_id::String = "bus.$(name)"
end

abstract type Eng3pBus <: EngNodeObject end

@doc raw"""
    Eng3pBus <: EngNodeObject

For three-phase buses, instead of specifying bounds explicitly for each pair of windings, often we want to specify 'phase-to-phase', 'phase-to-neutral' and 'neutral-to-ground' bounds. This can be done conveniently with a number of additional fields. First, `phases` is a list of the phase terminals, and `neutral` designates a single terminal to be the neutral.

- The bounds `vm_pn_lb` and `vm_pn_ub` specify the same lower and upper bound for the magnitude of the difference of each phase terminal and the neutral.
- The bounds `vm_pp_lb` and `vm_pp_ub` specify the same lower and upper bound for the magnitude of the difference of all phase terminals.
- `vm_ng_ub` specifies an upper bound for the neutral terminal, the lower bound is typically zero.

If all of these are specified, these bounds also imply valid bounds for the individual voltage magnitudes,

- $$\forall c \in phases: vm_pn_lb - vm_ng_ub \leq |v[c]| \leq vm_pn_ub + vm_ng_ub$$
- $$0 \leq |v[neutral]| \leq vm_ng_ub$$

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
"""
Base.@kwdef mutable struct Eng3pBusObj <: Eng3pBus
    name::String
    phases::Vector{Int}
    neutral::Int = 4
    vm_pn_lb::Union{Missing,Float64} = missing
    vm_pn_ub::Union{Missing,Float64} = missing
    vm_pp_lb::Union{Missing,Float64} = missing
    vm_pp_ub::Union{Missing,Float64} = missing
    vm_ng_ub::Union{Missing,Float64} = missing
    source_id::String
end

abstract type EngLine <: EngEdgeObject end
"""
    EngLine <: EngEdgeObject

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
"""
Base.@kwdef mutable struct EngLineObj <: EngLine
    name::String
    f_bus::String
    t_bus::String
    f_connections::Vector{Int}
    t_connections::Vector{Int}
    linecode::String = ""
    rs::Union{Missing,Matrix{Float64}} = missing
    xs::Union{Missing,Matrix{Float64}} = missing
    g_fr::Union{Missing,Matrix{Float64}} = isempty(linecode) ? zeros(Float64, Base.length(f_connections), Base.length(t_connections)) : missing
    b_fr::Union{Missing,Matrix{Float64}} = isempty(linecode) ? zeros(Float64, Base.length(f_connections), Base.length(t_connections)) : missing
    g_to::Union{Missing,Matrix{Float64}} = isempty(linecode) ? zeros(Float64, Base.length(f_connections), Base.length(t_connections)) : missing
    b_to::Union{Missing,Matrix{Float64}} = isempty(linecode) ? zeros(Float64, Base.length(f_connections), Base.length(t_connections)) : missing
    length::Float64 = 1.0
    cm_ub::Union{Missing,Vector{Float64}} = isempty(linecode) ? fill(Inf, Base.length(f_connections)) : missing
    sm_ub::Union{Missing,Vector{Float64}} = isempty(linecode) ? fill(Inf, Base.length(f_connections)) : missing
    vad_lb::Union{Missing,Vector{Float64}} = fill(-5.0, Base.length(f_connections))
    vad_ub::Union{Missing,Vector{Float64}} = fill(5.0, Base.length(f_connections))
    status::Status = ENABLED
    source_id::String = "line.$(name)"
    dss::Union{Missing,DssLine,DssCapacitor,DssReactor} = missing
end

abstract type EngSwitch <: EngEdgeObject end
"""
    EngSwitch <: EngEdgeObject

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
| `length`        | `0.001`                  | `Float64`      | m     | always | Length of switch in meters                                                                       |
| `dispatchable`  | `NO`                     | `Dispatchable` |       |        | `NO` or `YES`, indicates whether switch state can be changed in a switching optimization problem |
| `state`         | `CLOSED`                 | `SwitchState`  |       | always | `CLOSED`: closed or `OPEN`: open, to indicate state of switch                                    |
| `status`        | `ENABLED`                | `Status`       |       | always | `ENABLED` or `DISABLED`. Indicates if component is enabled or disabled, respectively             |
"""
Base.@kwdef mutable struct EngSwitchObj <: EngSwitch
    name::String
    f_bus::String
    t_bus::String
    f_connections::Vector{Int}
    t_connections::Vector{Int}
    linecode::String = ""
    rs::Union{Missing,Matrix{Float64}} = missing
    xs::Union{Missing,Matrix{Float64}} = missing
    g_fr::Union{Missing,Matrix{Float64}} = isempty(linecode) ? zeros(Float64, Base.length(f_connections), Base.length(t_connections)) : missing
    b_fr::Union{Missing,Matrix{Float64}} = isempty(linecode) ? zeros(Float64, Base.length(f_connections), Base.length(t_connections)) : missing
    g_to::Union{Missing,Matrix{Float64}} = isempty(linecode) ? zeros(Float64, Base.length(f_connections), Base.length(t_connections)) : missing
    b_to::Union{Missing,Matrix{Float64}} = isempty(linecode) ? zeros(Float64, Base.length(f_connections), Base.length(t_connections)) : missing
    length::Float64 = 0.001
    cm_ub::Union{Missing,Vector{Float64}} = isempty(linecode) ? fill(Inf, Base.length(f_connections)) : missing
    sm_ub::Union{Missing,Vector{Float64}} = isempty(linecode) ? fill(Inf, Base.length(f_connections)) : missing
    dispatchable::Dispatchable = YES
    state::SwitchState = CLOSED
    status::Status = ENABLED
    source_id::String = "switch.$(name)"
    dss::Union{Missing,DssLine} = missing
end


abstract type EngTransformerControls <: EngControlObject end
"""
    EngTransformerControls <: EngControlObject

Special case of the Generic transformer, which is part of the `transformer` object, and emulates a standard utility voltage regulator. The taps of these transformers can be controlled by modelling a line drop compensator.

| Name      | Default | Type                   | Units | Used | Description                                                                                                                                                                |
| --------- | ------- | ---------------------- | ----- | ---- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `vreg`    |         | `Vector{Vector{Real}}` | volt  | oltc | Voltage regulator reference, default value is `120.0` for the controlled winding, `0.0` for winding without regulator control, `size=((nphases),nwindings)`                |
| `band`    |         | `Vector{Vector{Real}}` | volt  | oltc | Voltage bandwidth, default value is `3.0` for the controlled winding, `0.0` for winding without regulator control, `size=((nphases),nwindings)`                            |
| `ptratio` |         | `Vector{Vector{Real}}` |       | oltc | Voltage ratio of the potential transformer, default value is `60.0` for the controlled winding, `0.0` for winding without regulator control, `size=((nphases),nwindings)`  |
| `ctprim`  |         | `Vector{Vector{Real}}` | amp   | oltc | Current transformer rating on primary side, default value is `300.0` for the controlled winding, `0.0` for winding without regulator control, `size=((nphases),nwindings)` |
| `r`       |         | `Vector{Vector{Real}}` | volt  | oltc | Resistance setting on line drop compensator, default value is `0.0` for both controlled winding and winding without regulator control, `size=((nphases),nwindings)`        |
| `x`       |         | `Vector{Vector{Real}}` | volt  | oltc | Reactance setting on line drop compensator, default value is `0.0` for both controlled winding and winding without regulator control, `size=((nphases),nwindings)`         |
"""
Base.@kwdef mutable struct EngTransformerControlsObj <: EngTransformerControls
    windings::Vector{Int}
    terminals::Vector{Vector{Int}}
    vreg::Vector{Vector{Float64}} = Vector{Float64}[fill(0.0, length(terminals[w])) for w in 1:length(windings)]
    band::Vector{Vector{Float64}} = Vector{Float64}[fill(0.0, length(terminals[w])) for w in 1:length(windings)]
    ptratio::Vector{Vector{Float64}} = Vector{Float64}[fill(0.0, length(terminals[w])) for w in 1:length(windings)]
    ctprim::Vector{Vector{Float64}} = Vector{Float64}[fill(0.0, length(terminals[w])) for w in 1:length(windings)]
    r::Vector{Vector{Float64}} = Vector{Float64}[fill(0.0, length(terminals[w])) for w in 1:length(windings)]
    x::Vector{Vector{Float64}} = Vector{Float64}[fill(0.0, length(terminals[w])) for w in 1:length(windings)]
    dss::Union{Missing,Vector{Vector{DssRegcontrol}}} = missing
end

abstract type EngAl2wTransformer <: EngEdgeObject end
"""
    EngAl2wTransformer <: EngEdgeObject

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
"""
Base.@kwdef mutable struct EngAl2wTransformerObj <: EngAl2wTransformer
    name::String
    f_bus::String
    t_bus::String
    f_connections::Vector{Int}
    t_connections::Vector{Int}
    configuration::ConnConfig = WYE
    tm_nom::Float64 = 1.0
    tm_ub::Union{Missing,Vector{Float64}} = fill( Inf, length(f_connections))
    tm_lb::Union{Missing,Vector{Float64}} = fill(-Inf, length(f_connections))
    tm_set::Vector{Float64} = fill(1.0, length(f_connections))
    tm_fix::Vector{Bool} = fill(true, length(f_connections))
    sm_ub::Union{Missing,Float64} = Inf
    cm_ub::Union{Missing,Float64} = Inf
    status::Status = ENABLED
    source_id::String = "al2wtransformer.$(name)"
end

abstract type EngTransformer <: EngEdgeObject end
"""
    EngTransformer <: EngEdgeObject

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
| `tm_step`        | `fill(fill(1/32,nphases),nwindings)` | `Vector{Vector{Real}}` |             | oltc   |                                                                                                                                                                |
| `tm_fix`         | `fill(fill(true,nphases),nwindings)` | `Vector{Vector{Bool}}` |             | oltc   | Indicates for each winding and phase whether the tap ratio is fixed, `size=((nphases),nwindings)`                                                              |
| `polarity`       | `fill(1,nwindings)`                  | `Vector{Int}`          |             | always |                                                                                                                                                                |
| `vm_nom`         |                                      | `Vector{Real}`         | volt        | always |                                                                                                                                                                |
| `sm_nom`         |                                      | `Vector{Real}`         | watt        | always |                                                                                                                                                                |
| `status`         | `ENABLED`                            | `Status`               |             | always | `ENABLED` or `DISABLED`. Indicates if component is enabled or disabled, respectively                                                                           |

"""
Base.@kwdef mutable struct EngTransformerObj <: EngTransformer
    name::String
    bus::Vector{String}
    connections::Vector{Vector{Int}}
    configurations::Union{Missing,Vector{ConnConfig}} = fill(WYE, length(bus))
    xfmrcode::String = ""
    xsc::Union{Missing,Vector{Union{Missing,Float64}}} = fill(0.0, length(bus) == 2 ? 1 : 3)
    rw::Union{Missing,Vector{Float64}} = fill(0.0, length(bus))
    cmag::Union{Missing,Float64} = 0.0
    noloadloss::Union{Missing,Float64} = 0.0
    tm_nom::Union{Missing,Vector{Float64}} = fill(1.0, length(bus))
    tm_ub::Union{Missing,Vector{Vector{Float64}}} = fill(fill(1.1, length(connections[1])), length(bus))
    tm_lb::Union{Missing,Vector{Vector{Float64}}} = fill(fill(0.9, length(connections[1])), length(bus))
    tm_set::Union{Missing,Vector{Vector{Float64}}} = fill(fill(1.0, length(connections[1])), length(bus))
    tm_step::Union{Missing,Vector{Vector{Float64}}} = fill(fill(1/32, length(connections[1])), length(bus))
    tm_fix::Union{Missing,Vector{Vector{Bool}}} = fill(fill(true, length(connections[1])), length(bus))
    polarity::Vector{Int} = fill(1, length(bus))
    vm_nom::Union{Missing,Vector{Float64}} = missing
    sm_nom::Union{Missing,Vector{Float64}} = missing
    sm_ub::Union{Missing,Float64} = Inf
    cm_ub::Union{Missing,Float64} = Inf
    bank::String = ""
    status::Status = ENABLED
    source_id::String = "transformer.$(name)"
    controls::Union{Missing,EngTransformerControls} = missing
    dss::Union{Missing,DssTransformer} = missing
end

abstract type EngLoad <: EngNodeObject end
@doc raw"""
    EngLoad <: EngNodeObject

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

# `model == EXPONENTIAL`

- `(pd[i]/pd_nom[i]) = (v[i]/vm_nom)^pd_exp[i]`
- `(qd[i]/qd_nom[i]) = (v[i]/vm_nom)^qd_exp[i]`

| Name     | Default | Type   | Units | Used                 | Description |
| -------- | ------- | ------ | ----- | -------------------- | ----------- |
| `pd_exp` |         | `Real` |       | `model==EXPONENTIAL` |             |
| `qd_exp` |         | `Real` |       | `model==EXPONENTIAL` |             |

# `model == ZIP`

- `(pd[i]/pd_nom) = pd_cz[i]*(v[i]/vm_nom)^2 + pd_ci[i]*(v[i]/vm_nom) + pd_cp[i]`
- `(qd[i]/qd_nom) = qd_cz[i]*(v[i]/vm_nom)^2 + qd_ci[i]*(v[i]/vm_nom) + qd_cp[i]`

| Name     | Default | Type   | Units | Used         | Description                  |
| -------- | ------- | ------ | ----- | ------------ | ---------------------------- |
| `vm_nom` |         | `Real` | volt  | `model==ZIP` | Nominal voltage (multiplier) |
| `pd_cz`  |         | `Real` |       | `model==ZIP` |                              |
| `pd_ci`  |         | `Real` |       | `model==ZIP` |                              |
| `pd_cp`  |         | `Real` |       | `model==ZIP` |                              |
| `qd_cz`  |         | `Real` |       | `model==ZIP` |                              |
| `qd_ci`  |         | `Real` |       | `model==ZIP` |                              |
| `qd_cp`  |         | `Real` |       | `model==ZIP` |                              |
"""
Base.@kwdef mutable struct EngLoadObj <: EngLoad
    name::String
    bus::String
    connections::Vector{Int} = Int[]
    configuration::ConnConfig = WYE
    model::LoadModel = POWER
    pd_nom::Vector{Float64} = fill(0.0, length(connections))
    qd_nom::Vector{Float64} = fill(0.0, length(connections))
    vm_nom::Union{Missing,Float64} = missing
    zipv::Union{Missing,Vector{Float64}} = missing
    dispatchable::Dispatchable = NO
    status::Status = ENABLED
    source_id::String = "load.$(name)"
    time_series::Dict{String,String} = Dict{String,String}(
        "pd_nom" => "",
        "qd_nom" => "",
    )
    dss::Union{Missing,DssLoad} = missing
end


abstract type EngShuntControls <: EngControlObject end
"""
    EngShuntControls <: EngControlObject

Special case of the shunt capacitors, which is part of the `shunt` object, and emulates a typical utility capacitor control (CapControl) by sending switching messages.

| Name           | Default | Type             | Units | Used | Description                                                                                                                                                                                       |
| -------------- | ------- | ---------------- | ----- | ---- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `type`         |         | `Vector{String}` |       | capc | Control type, default is `current` for controlled phase, ``for uncontrolled phase,`size=1`for`kvar`type, otherwise`size=(nphases)`                                                                |
| `element`      |         | `String`         |       | capc | `source_id` of element (typically line or transformer) to which CapControl is connected                                                                                                           |
| `terminal`     |         | `Vector{Int}`    |       | capc | Number of the terminal of circuit element to which CapControl is connected, default is `1` for controlled phase, `0` for uncontrolled phase, `size=1` for `kvar` type, otherwise `size=(nphases)` |
| `onsetting`    |         | `Vector{Real}`   |       | capc | Value at which the CapControl switches the capacitor on, default is `300.0` for controlled phase, `0.0` for uncontrolled phase, `size=1` for `kvar` type, otherwise `size=(nphases)`              |
| `offsetting`   |         | `Vector{Real}`   |       | capc | Value at which the CapControl switches the capacitor off, default is `200.0` for controlled phase, `0.0` for uncontrolled phase, `size=1` for `kvar` type, otherwise `size=(nphases)`             |
| `voltoverride` |         | `Vector{Bool}`   |       | capc | Indicate whether voltage over ride is enabled, default is `false` for both controlled and uncontrolled phases, `size=1` for `kvar` type, otherwise `size=(nphases)`                               |
| `ptratio`      |         | `Vector{Real}`   |       | capc | Ratio of potential transformer, default is `60.0` for controlled phase, `0.0` for uncontrolled phase, `size=(nphases)`                                                                            |
| `ctratio`      |         | `Vector{Real}`   |       | capc | Ratio of current transformer, default is `60.0` for controlled phase, `0.0` for uncontrolled phase, `size=(nphases)`                                                                              |
| `vmin`         |         | `Vector{Real}`   | volt  | capc | Minimum voltage below which CapControl switches the capacitor on, default is `115.0` for controlled phase, `0.0` for uncontrolled phase, `size=1` for `kvar` type, otherwise `size=(nphases)`     |
| `vmax`         |         | `Vector{Real}`   | volt  | capc | Maximum voltage above which CapControl switches the capacitor off, default is `126.0` for controlled phase, `0.0` for uncontrolled phase, `size=1` for `kvar` type, otherwise `size=(nphases)`    |
"""
Base.@kwdef mutable struct EngShuntControlsObj <: EngShuntControls
    type::Vector{CapControlType}
    element::String
    terminal::Vector{Int}
    onsetting::Vector{Float64} = fill(0.0, length(terminal))
    offsetting::Vector{Float64} = fill(0.0, length(terminal))
    voltoverride::Vector{Bool} = fill(false, length(terminal))
    ptratio::Vector{Float64} = fill(0.0, length(terminal))
    ctratio::Vector{Float64} = fill(0.0, length(terminal))
    vm_lb::Union{Missing,Vector{Float64}} = fill(-Inf, length(terminal))
    vm_ub::Union{Missing,Vector{Float64}} = fill( Inf, length(terminal))
    dss::Union{Missing,Vector{<:DssObject},DssObject} = missing
end

abstract type EngShunt <: EngNodeObject end
"""
    EngShunt <: EngNodeObject

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

"""
Base.@kwdef mutable struct EngShuntObj <: EngShunt
    name::String
    bus::String
    connections::Vector{Int} = Int[]
    configuration::ConnConfig = WYE
    gs::Matrix{Float64} = fill(0.0, length(connections))
    bs::Matrix{Float64} = fill(0.0, length(connections))
    model::ShuntModel = GENERIC
    dispatchable::Dispatchable = NO
    status::Status = ENABLED
    source_id::String = "shunt.$(name)"
    controls::Union{Missing,EngShuntControls} = missing
    dss::Union{Missing,DssObject} = missing
end

abstract type EngGenerator <: EngNodeObject end
@doc raw"""
    EngGenerator <: EngNodeObject

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

# `generator` Cost Model

The generator cost model is currently specified by the following fields.

| Name                 | Default           | Type           | Units | Used | Description                                               |
| -------------------- | ----------------- | -------------- | ----- | ---- | --------------------------------------------------------- |
| `cost_pg_model`      | `2`               | `Int`          |       | opf  | Cost model type, `1` = piecewise-linear, `2` = polynomial |
| `cost_pg_parameters` | `[0.0, 1.0, 0.0]` | `Vector{Real}` | $/MVA | opf  | Cost model polynomial                                     |
"""
Base.@kwdef mutable struct EngGeneratorObj <: EngGenerator
    name::String
    bus::String
    connections::Vector{Int} = Int[]
    configuration::ConnConfig = WYE
    vg::Union{Missing,Vector{Float64}} = missing
    pg_lb::Union{Missing,Vector{Float64}} = fill(-Inf, length(connections))
    pg_ub::Union{Missing,Vector{Float64}} = fill( Inf, length(connections))
    qg_lb::Union{Missing,Vector{Float64}} = fill(-Inf, length(connections))
    qg_ub::Union{Missing,Vector{Float64}} = fill( Inf, length(connections))
    pg::Vector{Float64} = fill(0.0, length(connections))
    qg::Vector{Float64} = fill(0.0, length(connections))
    control_mode::ControlMode = FREQUENCYDROOP
    cost_pg_model::Int = 2
    cost_pg_parameters::Vector{Float64} = [0.0, 1.0, 0.0]
    status::Status = ENABLED
    source_id::String = "generator.$(name)"
    time_series::Dict{String,String} = Dict{String,String}(
        "pg" => "",
        "qg" => "",
    )
    dss::Union{Missing,DssGenerator} = missing
end


abstract type EngSolar <: EngNodeObject end
@doc raw"""
    EngSolar <: EngNodeObject

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

# `solar` Cost Model

The cost model for a photovoltaic system currently matches that of generators.

| Name                 | Default           | Type           | Units | Used | Description                                               |
| -------------------- | ----------------- | -------------- | ----- | ---- | --------------------------------------------------------- |
| `cost_pg_model`      | `2`               | `Int`          |       | opf  | Cost model type, `1` = piecewise-linear, `2` = polynomial |
| `cost_pg_parameters` | `[0.0, 1.0, 0.0]` | `Vector{Real}` | $/MVA | opf  | Cost model polynomial                                     |
"""
Base.@kwdef mutable struct EngSolarObj <: EngSolar
    name::String
    bus::String
    connections::Vector{Int} = Int[]
    configuration::ConnConfig = WYE
    vg::Union{Missing,Vector{Float64}} = missing
    pg_lb::Union{Missing,Vector{Float64}} = fill(-Inf, length(connections))
    pg_ub::Union{Missing,Vector{Float64}} = fill(Inf, length(connections))
    qg_lb::Union{Missing,Vector{Float64}} = fill(-Inf, length(connections))
    qg_ub::Union{Missing,Vector{Float64}} = fill(Inf, length(connections))
    pg::Vector{Float64} = fill(0.0, length(connections))
    qg::Vector{Float64} = fill(0.0, length(connections))
    cost_pg_model::Int = 2
    cost_pg_parameters::Vector{Float64} = [0.0, 1.0, 0.0]
    status::Status = ENABLED
    source_id::String = "solar.$(name)"
    time_series::Dict{String,String} = Dict{String,String}()
    dss::Union{Missing,DssPvsystem} = missing
end

abstract type EngStorage <: EngNodeObject end
@doc raw"""
    EngStorage <: EngNodeObject

A storage object is a flexible component that can represent a variety of energy storage objects, like Li-ion batteries, hydrogen fuel cells, flywheels, etc.

- TODO: How to include the inverter model for this? Similar issue as for a PV generator

| Name                   | Default   | Type                  | Units   | Used         | Description                                                                          |
| ---------------------- | --------- | --------------------- | ------- | ------------ | ------------------------------------------------------------------------------------ |
| `bus`                  |           | `String`              |         | always       | id of bus connection                                                                 |
| `connections`          |           | `Vector{Int}`         |         | always       | Ordered list of connected conductors, `size=nconductors`                             |
| `configuration`        | `WYE`     | `ConnConfig`          |         | always       | `WYE` or `DELTA`. If `WYE`, `connections[end]=neutral`                               |
| `energy`               |           | `Real`                | watt-hr | always       | Stored energy                                                                        |
| `energy_ub`            |           | `Real`                |         | opf          | maximum energy rating                                                                |
| `charge_ub`            |           | `Real`                |         | opf          | maximum charge rating                                                                |
| `discharge_ub`         |           | `Real`                |         | opf          | maximum discharge rating                                                             |
| `sm_ub`                |           | `Vector{Real}`        | watt    | opf          | Power rating, `size=nphases`                                                         |
| `cm_ub`                |           | `Vector{Real}`        | amp     | opf          | Current rating, `size=nphases`                                                       |
| `charge_efficiency`    |           | `Real`                | percent | always       | charging efficiency (losses)                                                         |
| `discharge_efficiency` |           | `Real`                | percent | always       | disharging efficiency (losses)                                                       |
| `qs_ub`                |           | `Vector{Real}`        |         | opf          | Maximum reactive power injection, `size=nphases`                                     |
| `qs_lb`                |           | `Vector{Real}`        |         | opf          | Minimum reactive power injection, `size=nphases`                                     |
| `rs`                   |           | `Vector{Real}`        | ohm     | always       | converter resistance                                                                 |
| `xs`                   |           | `Vector{Real}`        | ohm     | always       | converter reactance                                                                  |
| `pex`                  |           | `Real`                |         | always       | Total active power standby exogenous flow (loss)                                     |
| `qex`                  |           | `Real`                |         | always       | Total reactive power standby exogenous flow (loss)                                   |
| `ps`                   |           | `Vector{Real}`        | watt    | solution     | Present active power injection                                                       |
| `qs`                   |           | `Vector{Real}`        | var     | solution     | Present reactive power injection                                                     |
| `status`               | `ENABLED` | `Status`              |         | always       | `ENABLED` or `DISABLED`. Indicates if component is enabled or disabled, respectively |
| `time_series`          |           | `Dict{String,String}` |         | multinetwork | Dictionary containing time series parameters.                                        |
"""
Base.@kwdef mutable struct EngStorageObj <: EngStorage
    name::String
    bus::String
    connections::Vector{Int}
    configuration::ConnConfig = WYE
    energy::Float64 = 0.0
    energy_ub::Union{Missing,Float64} = 0.0
    charge_ub::Union{Missing,Float64} = 0.0
    discharge_ub::Float64 = 0.0
    sm_ub::Union{Missing,Float64} = Inf
    cm_ub::Union{Missing,Float64} = Inf
    charge_efficiency::Float64 = 1.0
    discharge_efficiency::Float64 = 1.0
    qs_ub::Union{Missing,Float64} = 0.0
    qs_lb::Union{Missing,Float64} = 0.0
    rs::Float64 = 0.0
    xs::Float64 = 0.0
    pex::Float64 = 0.0
    qex::Float64 = 0.0
    ps::Float64 = 0.0
    qs::Float64 = 0.0
    status::Status = ENABLED
    source_id::String = "storage.$(name)"
    time_series::Dict{String,String} = Dict{String,String}()
    dss::Union{Missing,DssStorage} = missing
end

abstract type EngVoltageSource <: EngNodeObject end
@doc raw"""
    EngVoltageSource <: EngNodeObject

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
"""
Base.@kwdef mutable struct EngVoltageSourceObj <: EngVoltageSource
    name::String
    bus::String
    connections::Vector{Int}
    configuration::ConnConfig = WYE
    vm::Vector{Float64} = fill(1.0, length(connections))
    va::Vector{Float64} = [[0.0, -120.0, 120.0]..., fill(0.0, length(connections))][connections]
    rs::Matrix{Float64} = fill(0.0, length(connections), length(connections))
    xs::Matrix{Float64} = fill(0.0, length(connections), length(connections))
    status::Status = ENABLED
    source_id::String = "voltage_source.$(name)"
    time_series::Dict{String,String} = Dict{String,String}(
        "pg" => "",
        "qg" => "",
        "pg_ub" => "",
        "qg_ub" => "",
    )  # TODO: remove time_series from voltage source
    dss::Union{Missing,DssVsource} = missing
end

abstract type EngLinecode <: EngDataObject end
@doc raw"""
    EngLinecode <: EngDataObject

Linecodes are easy ways to specify properties common to multiple lines.

| Name    | Default                          | Type           | Units            | Used   | Description                                             |
| ------- | -------------------------------- | -------------- | ---------------- | ------ | ------------------------------------------------------- |
| `rs`    |                                  | `Matrix{Real}` | ohm/meter        | always | Series resistance, `size=(nconductors,nconductors)`     |
| `xs`    |                                  | `Matrix{Real}` | ohm/meter        | always | Series reactance, `size=(nconductors,nconductors)`      |
| `g_fr`  | `zeros(nconductors,nconductors)` | `Matrix{Real}` | siemens/meter/Hz | always | From-side conductance, `size=(nconductors,nconductors)` |
| `b_fr`  | `zeros(nconductors,nconductors)` | `Matrix{Real}` | siemens/meter/Hz | always | From-side susceptance, `size=(nconductors,nconductors)` |
| `g_to`  | `zeros(nconductors,nconductors)` | `Matrix{Real}` | siemens/meter/Hz | always | To-side conductance, `size=(nconductors,nconductors)`   |
| `b_to`  | `zeros(nconductors,nconductors)` | `Matrix{Real}` | siemens/meter/Hz | always | To-side susceptance, `size=(nconductors,nconductors)`   |
| `cm_ub` |                                  | `Vector{Real}` | ampere           | opf    | maximum current per conductor, symmetrically applicable |
| `sm_ub` |                                  | `Vector{Real}` | watt             | opf    | maximum power per conductor, symmetrically applicable   |
"""
Base.@kwdef mutable struct EngLinecodeObj <: EngLinecode
    name::String
    rs::Matrix{Float64}
    xs::Matrix{Float64}
    g_fr::Matrix{Float64} = fill(0.0, size(rs))
    b_fr::Matrix{Float64} = fill(0.0, size(rs))
    g_to::Matrix{Float64} = fill(0.0, size(rs))
    b_to::Matrix{Float64} = fill(0.0, size(rs))
    cm_ub::Union{Missing,Vector{Float64}} = fill(Inf, size(rs)[1])
    sm_ub::Union{Missing,Vector{Float64}} = fill(Inf, size(rs)[1])
    source_id::String = "linecode.$(name)"
    dss::Union{Missing,DssLinecode} = missing
end

abstract type EngXfmrcode <: EngDataObject end

@doc raw"""
    EngXfmrcode <: EngDataObject

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
| `tm_step`        | `fill(fill(1/32,nphases),nwindings)`   | `Vector{Vector{Real}}` |       | oltc   |                                                                                                                                                 |
| `tm_fix`         | `fill(fill(true, nphases), nwindings)` | `Vector{Vector{Bool}}` |       | always | Indicates for each winding and phase whether the tap ratio is fixed, `size=((nphases), nwindings)`                                              |
"""
Base.@kwdef mutable struct EngXfmrcodeObj <: EngXfmrcode
    name::String
    configurations::Vector{ConnConfig} = fill(WYE, 2)
    xsc::Vector{Float64} = fill(0.0, length(configurations) == 2 ? 1 : 3)
    rw::Vector{Float64} = fill(0.0, length(configurations))
    noloadloss::Float64 = 0.0
    cmag::Float64 = 0.0
    tm_nom::Union{Missing,Vector{Float64}} = fill(1.0, length(configurations))
    tm_fix::Vector{Vector{Bool}} = fill(fill(true, 3), 2)
    tm_ub::Union{Missing,Vector{Vector{Float64}}} = fill(fill( Inf, length(tm_fix[1])), length(configurations))
    tm_lb::Union{Missing,Vector{Vector{Float64}}} = fill(fill(-Inf, length(tm_fix[1])), length(configurations))
    tm_set::Vector{Vector{Float64}} = fill(fill( Inf, length(tm_fix[1])), length(configurations))
    tm_step::Vector{Vector{Float64}} = fill(fill(1/32, length(tm_fix[1])), length(configurations))
    vm_nom::Vector{Float64} = fill(12.47, length(configurations))
    sm_nom::Vector{Float64} = fill(10.0, length(configurations))
    sm_ub::Union{Missing,Float64} = Inf
    cm_ub::Union{Missing,Float64} = Inf
    source_id::String = "xfmrcode.$(name)"
    dss::Union{Missing,DssXfmrcode} = missing
end

abstract type EngTimeSeries <: EngDataObject end
@doc raw"""
    EngTimeSeries <: EngDataObject

Time series objects are used to specify time series for _e.g._ load or generation forecasts.

    Some parameters for components specified in this document can support a time series by inserting a referece to a `time_series` object into the `time_series` dictionary inside a component under the relevant parameter name. For example, for a `load`, if `pd_nom` is supposed to be a time series, the user would specify `"time_series" => Dict("pd_nom" => time_series_id)` where `time_series_id` is the `id` of an object in `time_series`, and has type `Any`.

    | Name      | Default | Type           | Units | Used   | Description                                                                                                   |
    | --------- | ------- | -------------- | ----- | ------ | ------------------------------------------------------------------------------------------------------------- |
    | `time`    |         | `Vector{Any}`  | hour  | always | Time points at which values are specified. If time is specified in String, units not required to be in hours. |
    | `values`  |         | `Vector{Real}` |       | always | Multipers at each time step given in `time`                                                                   |
    | `offset`  | `0`     | `Real`         | hour  | always | Start time offset                                                                                             |
    | `replace` | `true`  | `Bool`         |       | always | Indicates to replace with data, instead of multiply. Will only work on non-Array data                         |

"""
Base.@kwdef mutable struct EngTimeSeriesObj <: EngTimeSeries
    name::String
    time::Union{Vector{String},Vector{Float64}}
    values::Vector{Float64}
    offset::Float64 = 0.0
    replace::Bool = true
    source_id::String = "time_series.$(name)"
    dss::Union{Missing,DssLoadshape,DssXycurve} = missing
end


"""
"""
Base.@kwdef struct EngineeringDataModel <: EngineeringModel{NetworkModel}
    # metadata
    settings::Settings = SettingsObj()
    metadata::Metadata = MetadataObj()

    # components
    bus::Dict{String,<:Union{EngBus,Eng3pBus}} = Dict{String,Union{EngBus,Eng3pBus}}()

    # edges
    line::Dict{String,<:EngLine} = Dict{String,EngLine}()
    switch::Dict{String,<:EngSwitch} = Dict{String,EngSwitch}()
    transformer::Dict{String,<:Union{EngTransformer,EngAl2wTransformer}} = Dict{String,Union{EngTransformer,EngAl2wTransformer}}()

    # nodes
    load::Dict{String,<:EngLoad} = Dict{String,EngLoad}()
    shunt::Dict{String,<:EngShunt} = Dict{String,EngShunt}()
    generator::Dict{String,<:EngGenerator} = Dict{String,EngGenerator}()
    voltage_source::Dict{String,<:EngVoltageSource} = Dict{String,EngVoltageSource}()
    solar::Dict{String,<:EngSolar} = Dict{String,EngSolar}()
    storage::Dict{String,<:EngStorage} = Dict{String,EngStorage}()

    # data
    linecode::Dict{String,<:EngLinecode} = Dict{String,EngLinecode}()
    xfmrcode::Dict{String,<:EngXfmrcode} = Dict{String,EngXfmrcode}()
    time_series::Dict{String,<:EngTimeSeries} = Dict{String,EngTimeSeries}()

    # additional user data
    extra::Dict{String,Any} = Dict{String,Any}()
end


"""
"""
Base.@kwdef mutable struct EngineeringMultinetworkDataModel <: EngineeringModel{MultinetworkModel}
    # global keys
    metadata::Metadata = MetadataObj()

    # multinetwork
    nw::Dict{String,EngineeringDataModel} = Dict{String,EngineeringDataModel}()
    nw_map::Dict{String,Real} = Dict{String,Real}()
end


""
Base.@kwdef mutable struct UnsupportedEngEdgeObject <: EngEdgeObject
    f_bus::String = ""
    t_bus::String = ""
    f_connections::Vector{Int} = Int[]
    t_connections::Vector{Int} = Int[]
end


""
Base.@kwdef mutable struct UnsupportedEngNodeObject <: EngNodeObject
    bus::String = ""
    connections::Vector{Int} = Int[]
end
