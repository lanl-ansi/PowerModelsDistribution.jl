ThreePhasePowerModels.jl Change Log
===================================

### staged
- Allow for arbitrarily named sourcebus
- Add json parser
- Fix bug in OpenDSS parse of Capacitors [zbase factor and wrong sign] (#138)

### v0.3.0
- Update to JuMP v0.19/MathOptInterface

### v0.2.1
- Add transformer support to active power only models
- Fix bug in source_ids of new components created for transformer support

### v0.2.0
- Add transformer to TPPM: wye and delta windings, vector group variations (indirectly) and variable taps per phase with bounds
- Add ACPForm for transformers
- Add transformer parsing to OpenDSS, including mapping of loss model
- Add voltage base calculation, and corresponding change of base
- Add AC rectangular formulation
- Remove support for Julia versions < 1.0
- Switch to using `import` instead of `using` internally

### v0.1.5
- Minor fix to OpenDSS line parsing (length units mismatch with linecode)

### v0.1.4
- Update to PowerModels v0.9

### v0.1.3
- Add opf with multi-phase storage model
- Add support for parsing OpenDSS Storage objects
- Minor fix to branch parsing in matlab format
- Minor fix to OpenDSS parser (parsing ~ lines with preceeding whitespace)
- Fixed parsing OpenDSS files containing redirect/compile/buscoords on case-sensitive filesystems
- Add 'source_id' field to components parsed from OpenDSS, to help determine origin and active phases
- Add conversion of OpenDSS PVSystem objects into generators, using KVA for generator limits
- Add compatibility for Julia v0.7/v1.0

### v0.1.2
- Add support for network flow approximation formulation, NFAPowerModel
- Updates to problem specifications
- Update tests for SCS v0.4
- Minor improvements to OpenDSS parser

### v0.1.1
- Added a variety of matrix-based branch flow formulations
- Updated LPLinUBFForm to more closely match published refrence model
- Update MINLP solvers used in testing
- Added basic docs setup
- Minor improvements to OpenDSS parser #29, #59, #62, #63, #64, #65

### v0.1.0
- Initial release
