# PowerModelsDistribution.jl Change Log

## staged

- Add support for Memento v0.13, v1.0
- SDP and SOC relaxations were broken but are fixed again (unit tests added)
- opf_bf_lm removed, functionality available in opf_bf instead
- add BIM SDP (KCL MX and diag) and SOC relaxations

## v0.8.1

- Update to support JuMP v0.21
- Makes bounds optional, turned on by default (#250)
- Updated transformer data model in the mathematical model (#250)
- Add automatic parsing of lon,lat from buscoords file into PMD data structure (#245, #249)
- Updates virtual_sourcebus, which is intended to represent a voltage source, to have a fixed voltage magnitude (#246,#248)
- Add parsing of series data files into array fields in OpenDSS parser
- Add LoadShape parsing to OpenDSS parser (#247)
- The pf and opf problem specifications now contain delta connected, and voltage-dependent load models by default; pf_lm and opf_lm were removed.
- Generators can now also be connected in both delta and wye.
- The delta/voltage-dependent loads in the ACP power models were refactored to be consistent with the new ones.
- The non-linear formulations now support matrix shunts.
- A current-voltage (IVR) formulation was added, and supports all features that ACP supports as well.
- Several power balance constraints now require NLconstraints because the load power/current can contain a NLexpression. This might be optimized further in the future if it leads to performance issues.

## v0.8.0

- Update solution building infrastructure (PMs #77) (breaking). The reported solution is now consistent with the variable space of the formulation.
- Moved multi-conductor support from PowerModels into PowerModelsDistribution. (breaking)
- PMs.var no longer takes conductor as an argument
- Constraints have been (partially) re-written to use vectorized JuMP syntax where possible.
- Bugfixes: generator on-off and storage on-off constraints were incorrect
- Removal of SOCWRPowerModel
- Drop support for matpower format
- Possible regressions in MLD problem types
- Simplified linear UBF formulations. LPLinUBFModel and LPdiagUBFModel are equivalent, and are replaced by LPUBFDiagModel. The more popular name 'LinDist3FlowModel' was added as an alias for LPUBFDiagModel.
- The linearization by Gan & Low is also equivalent to LPUBFDiagModel, though it has redundudant variables and equations. LPfullUBFModel is now implemented as an alias of LPUBFDiagModel, 'LPUBFFullModel'.

## v0.7.0

- Updates function names for PowerModels v0.14 compatibility (breaking) (#194)

## v0.6.1

- Fixed bug with new default SCS settings causing tests to fail (#190)
- Changed unit test '5-bus independent radial different sdpubf opf_bf', testing vm instead of qg
- Added exponential load models, including convex relaxations
- Refactored branch flow relaxations
- Added helper functions for creating matrix variables
- Cleaned up the transformer/oltc methods, including ACP formulation
- added `rate_a` setting to virtual lines (#185, #186)

## v0.6.0

- Adds `count_nodes` function to count the number of nodes in a parsed network (#183)
- Exports `find_bus` and `find_component` functions for better user experience (#183)
- Fixed `solution_bf!` for branch flow solution building (#182)
- Refactored problem definitions to remove any explicit loops over conductors (#181)
- Added data format documentation (#180)
- Moved storage to main MLD and OPF problems (#179)
- Refactor to remove dcline variables and constraints (#179)
- Refactor to genericize `constraint_mc_power_balance` (#179)
- Fix bug in OpenDSS circuit initialization (vsource) (#178)
- Make current rating (c_rating_a|b|c) be the default on branches (breaking)
- Fix bug in transformer `ref` extension where all refs were not built for all `nw` in multinetworks (#171)
- Fix bug in OpenDSS parser where properties were not applied in the order they were received (#170)
- Rename "trans" in data and ref to `transformer` for component naming consistency (breaking) (#169)
- Change internal variable and constraint functions to loop over phases internally (breaking) (#168)
- Fix bug in OpenDSS parser on Lines where the connected phases are listed out of order (#167)
- Add ability to "bank" single phase OpenDSS transformers into a single multiphase transformer (#166)
- Add virtual line to sourcebus to model source impedance (#165)
- Update to JuMP v0.20 / MOI v0.9 (#164)
- Fix bug in OpenDSS parser on Lines / Linecodes related to basefreq (#163)
- Fix bug in OpenDSS parser on Transformers (#162)
- Fix bug in OpenDSS parser on Lines where `switch=y` property is used (#161)
- Update Formulation types to follow PowerModels v0.13 conventions (breaking) (#160)

## v0.5.2

- Fix bug in OpenDSS parser on Capacitors (#158)
- Add support for full matrix line shunts (#153)

## v0.5.1

- Add continuous load shedding problem (mld)

## v0.5.0

- Enforce function naming conventions (starts with `_`: internal function; ends with `!`: transforms data; `correct_`: corrects network data; `check_`: warnings about network data) (breaking)
- Update for PowerModels.jl v0.12 (breaking)
- Enforce constraint/variable naming conventions to include `_tp` (breaking)
- Add automatic export of non-internal functions (all functions not prefixed with `_`)
- Enforce function naming conventions (starts with `_`: internal function; ends with `!`: transforms data; `correct_`: corrects network data; `check_`: warnings about network data)

## v0.4

- First version of PowerModelsDistribution.jl

## v0.3.2

- Final version of ThreePhasePowerModels.jl before name change to PowerModelsDistribution.jl (adds depreciation warnings)

## v0.3.1

- Allow for arbitrarily named sourcebus
- Add json parser
- Add support for additional load models (constant power, constant impedance, constant current; delta or wye connected) (#127)
- Fix bug in OpenDSS parse of Capacitors [zbase factor and wrong sign] (#138)
- Add voltage balance constraints (#129)

## v0.3.0

- Update to JuMP v0.19/MathOptInterface

## v0.2.1

- Add transformer support to active power only models
- Fix bug in source_ids of new components created for transformer support

## v0.2.0

- Add transformer to TPPM: wye and delta windings, vector group variations (indirectly) and variable taps per phase with bounds
- Add ACPForm for transformers
- Add transformer parsing to OpenDSS, including mapping of loss model
- Add voltage base calculation, and corresponding change of base
- Add AC rectangular formulation
- Remove support for Julia versions < 1.0
- Switch to using `import` instead of `using` internally

## v0.1.5

- Minor fix to OpenDSS line parsing (length units mismatch with linecode)

## v0.1.4

- Update to PowerModels v0.9

## v0.1.3

- Add opf with multi-phase storage model
- Add support for parsing OpenDSS Storage objects
- Minor fix to branch parsing in matlab format
- Minor fix to OpenDSS parser (parsing ~ lines with preceeding whitespace)
- Fixed parsing OpenDSS files containing redirect/compile/buscoords on case-sensitive filesystems
- Add 'source_id' field to components parsed from OpenDSS, to help determine origin and active phases
- Add conversion of OpenDSS PVSystem objects into generators, using KVA for generator limits
- Add compatibility for Julia v0.7/v1.0

## v0.1.2

- Add support for network flow approximation formulation, NFAPowerModel
- Updates to problem specifications
- Update tests for SCS v0.4
- Minor improvements to OpenDSS parser

## v0.1.1

- Added a variety of matrix-based branch flow formulations
- Updated LPLinUBFForm to more closely match published refrence model
- Update MINLP solvers used in testing
- Added basic docs setup
- Minor improvements to OpenDSS parser #29, #59, #62, #63, #64, #65

## v0.1.0

- Initial release
