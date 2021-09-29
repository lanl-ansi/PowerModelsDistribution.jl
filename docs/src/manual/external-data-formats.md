# External Data Formats

## OpenDSS

PowerModelsDistribution supports parsing OpenDSS format files. In particular, we support a raw parsing of all OpenDSS specified components into a dictionary (serializable) structure, and support converting the data values of a smaller subsection of components into their expected types. Those include

- Line
- Load
- Generator
- Capactior (shunt capacitors only)
- Reactor
- Transformer
- Linecode
- Xfmrcode
- Loadshape
- XYCurve
- Circuit
- VSource
- PVSystem
- Storage

Of those, a subset of configurations are converted into a PowerModelsDistribution internal data model, namely:

### Edge Elements

- line (from lines and line reactors)
- transformer (arbitrary winding, all connections except zig-zag)
- switch (from lines w/ switch=y)

### Node Elements

- generator
- voltage_source
- solar (from PVSystem)
- load (incl. support for constant [`POWER`](@ref POWER), constant [`IMPEDANCE`](@ref IMPEDANCE), constant [`CURRENT`](@ref CURRENT), and [`EXPONENTIAL`](@ref EXPONENTIAL) models)
- shunt (from shunt capacitors and shunt reactors)
- storage

### Data Elements

- linecode
- xfmrcode
- time_series (from loadshapes)

Several notes about the specific design choices w.r.t. OpenDSS are explained below.

### Circuit

The default connection to the transmission system is modeled as an ideal voltage source named "source" in OpenDSS, which is connected by default to a node named "sourcebus", but this can be changed.

### Lines

Although Lines and Linecodes are supported, the more generic LineGeometry is not yet supported.

### Transformers

Unfortunately, in the OpenDSS format, multi-phase transformers with different taps for each phase are not explicitly supported, so to work around this limitation multiple single phase transformers should be defined, which are then "banked" together using the `bank` property.

### Capacitors and Reactors

Capacitors and reactors are supported as shunts, although shunts to ground via delta connections are not yet supported. Furthermore, generic reactors are not supported, only those whose second terminal is connected to ground (default for unspecified second terminal). Reactors are also supported as a resistanceless line if their second terminal is connected, but only for topological continuity of the network.

## PowerModelsDistribution JSON

You can export a PowerModelsDistribution data structure to a JSON file using the [`print_file`](@ref print_file) command and parse one in using the [`parse_file`](@ref parse_file) command
