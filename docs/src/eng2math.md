# Engineering to Mathematical Data Model Mapping

In this document we define the mapping from the engineering data model down to the mathematical data model for each physical component.

## `bus`

Buses are parsed into `bus` and potentially `shunt` objects

## `line`

Lines are parsed into `branch` objects with `transformer=false`

## `switch`

Switches are parsed into `switch`. If there are loss parameters provided (_i.e._ `rs` and/or `xs`) then a virtual branch and virtual bus are created to model the impedance

## `transformer`

Transformers are parsed into asymmetric lossless 2-winding transformers. When parsing n-winding transformers with n>2 additionally virtual branches and buses are created to connect the new 2-winding transformers. Furthermore, if the loss parameters are non-zero, additional virtual buses and branches to model the transformer impedances

## `shunt`

Shunts are parsed directly into `shunt` objects.

## `load`

Loads are parsed into `load` objects, with a specialized model that can be found in Load Model documentation on the sidebar

## `generator`

Generators are parsed into `gen` objects

## `solar`

Solar objects (photovoltaic systems) are parsed into `gen` objects

## `voltage_source`

Voltage sources are parsed into `gen` objects. If loss parameters are specified (_i.e._ `rs` and/or `xs`) then a virtual bus and branch are created to model the internal impedance
