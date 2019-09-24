# Data Formats

## OpenDSS

PowerModelsDistribution supports parsing OpenDSS format files. In particular, we support a raw parsing of all OpenDSS specified components into a dictionary (serializable) structure, and support converting the data values of a smaller subsection of components into their expected types. Those include

- Line
- Load
- Generator
- Capactior
- Reactor
- Transformer
- Linecode
- Circuit
- VSource
- PVSystem
- Storage

Of those, a subset of configurations are converted into a PowerModelsDistribution internal data model, namely

- Branch (from Lines (incl. Linecodes), Reactors)
- Transformer (arbitrary winding, all connections except zig-zag)
- Generator (from Generators, PVSystems)
- Load (incl. support for Const. Power, Const. Impedance, Const. Current models)
- Shunt (from Capacitors and Reactors)
- Storage

Except for a small subset, in general, commands are not support, e.g. `solve` or `calcvoltagebases` (this is done automatically on parse in PowerModelsDistribution). We support the following commands

- `clear`
- `redirect`
- `compile`
- `set`
- `buscoords`
- `new`

Several notes about the specific design choices w.r.t. OpenDSS are explained below.

### Circuit

The default connection to the transmission system is modeled as an ideal voltage source in OpenDSS; we chosen to model the trunk connection as a loosely bounded generator at a reference bus which is connected to the distribution network via a branch in order to model the inherent impedance of the voltage source.

### Lines

Although Lines and Linecodes are supported, the more generic LineGeometry is not yet supported.

### Transformers

Unfortunately, in the OpenDSS format, multi-phase transformers with different taps for each phase are not explicitly supported, so to work around this limitation multiple single phase transformers should be defined, which are then "banked" together using the `bank` property.

### Capacitors and Reactors

Capacitors and reactors are supported as shunts, although shunts to ground via delta connections are not yet supported. Furthermore, generic reactors are not supported, only those whose second terminal is wye connected to ground (default for unspecified second terminal). Reactors are also supported as a resistanceless line if their second terminal is connected, but only for topological continuity of the network.

## Matlab

We also include a matlab-base format similar in conception to Matpower. This format is in development and details will come later.
