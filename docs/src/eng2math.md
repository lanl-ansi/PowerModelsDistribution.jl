# Engineering to Mathematical Data Model Mapping

In this document we define the mapping from the engineering data model down to the mathematical data model for each physical component.

## `bus`

Buses are parsed into `bus` and potentially `shunt` objects.

The mathematical bus model contains only lossless connections to ground. All other connections to grounds are converted to equivalent shunts at that bus. For example, take a bus defined as

`bus_eng = Dict("grounded"=>[4, 5], "rg"=>[1.0, 0.0], "xg"=>[2.0, 0.0],...)`.

This is equivalent to a shunt `g+im*b = 1/(1.0+im*2.0)` connected to terminal `4`, and a lossless grounding at terminal `5` (since `rg[2]==xg[2]==0.0`). This is mapped to

`bus_math = Dict("grounded"=>[5], ...)`,

`shunt_math = Dict("connections"=>[4], "b"=>[b], "g"=>[g]...)`.

This simplifies the mathematical model, as the modeller does no longer have to consider lossy groundings explicitly.

## `line`

Lines are parsed into `branch` objects with `transformer=false`

## `switch`

Switches are parsed into `switch`. If there are loss parameters provided (_i.e._ `rs` and/or `xs`) then a virtual branch and virtual bus are created to model the impedance

## `transformer`

A transformer can have N windings, each with its own configuration (`delta` or `wye` are supported). This is decomposed to a network of N lossless, two-winding transformers which connect to an internal loss model. The to-winding is always wye-connected, hence we refer to these transformers as 'asymmetric'.

The internal loss model is a function of
- the winding resistance `rw`,
- the short-circuit reactance `xsc`,
- the no-load loss properties `noloadloss` (resistive) and magnetizing current `imag` (reactive).

If all of these are non-zero, this leads to an internal loss model consisting of `N` virtual buses, `(N^2+N)/2` virtual branches, and `1` shunt. These virtual buses and branches are automatically merged and simplified whenever possible; e.g., when all these loss parameters are zero, this simplifies to a single virtual bus, to which all two-winding transformers connect.

For more detail, please refer to [upcoming technical paper]. #TODO add link to paper

## `shunt`

Shunts are parsed directly into `shunt` objects.

## `load`

Loads are parsed into `load` objects. See the discussion under the Load Model documentation on the sidebar, for a detailed discussion of the various load models.

## `generator`

Generators are parsed into `gen` objects.

## `solar`

Solar objects (photovoltaic systems) are parsed into `gen` objects.

## `voltage_source`

Voltage sources are parsed into `gen` objects. If loss parameters are specified (_i.e._ `rs` and/or `xs`) then a virtual bus and branch are created to model the internal impedance.
