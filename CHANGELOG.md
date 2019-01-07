ThreePhasePowerModels.jl Change Log
===================================

### Staged
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
