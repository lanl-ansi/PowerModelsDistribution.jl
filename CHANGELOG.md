# PowerModelsDistribution.jl Change Log

## staged

- Fix bug in buscoords parser where comments at the end of the line were not stripped

## v0.10.2

- add `exclude` kwarg to `remove_all_bounds!` transformation, to selectively exclude certain asset types
- fix bug in IVR transformer current variables where bounds needed to be iterated over their connections to be applied
- fix bug in objective function for opf_pbs debugging problem, wrong iteration over bus terminals
- fix typo in `variable_mc_switch_current_imaginary` that overwrote real variables (crsw)
- fix bug in `_map_eng2math_switch!` where vmin/vmax were taken from f_bus instead of t_bus
- fix bug in `_build_eng_multinetwork`, where "dss_options" was missing from const `_pmd_eng_global_keys`
- change enums (SwitchState and Dispatchable) for switches to Reals, was causing problems in loops of OSW problems
- fix bug in `variable_mc_bus_voltage_magnitude_sqr` and `variable_mc_transformer_power_imaginary` where `_start` values were not being iterated over per connection
- depreciate run_ functions in favor of solve_
- add support for `relax_integrality` (InfrastructureModels ~0.5.4)
- fix bug in `variable_mx_real` constructor where it was indexing over terminals instead of enumerates
- added storage variables to automatic unit conversion to si units on math2eng transformation

## v0.10.1

- Fix buspairs ref not getting built
- Fix bug in storage linear constraints and mixed integer variable instantiation
- Add support for ideal (lossless) switches
- Fix bug in handling of dss edit command
- Fix bug in correctly tracking current transformer winding
- Fixes bug in assignment of variables' start values over active conductors/connections

## v0.10.0

- Refactor variables, constraints, objectives to support iterating over arbitrary connections/terminals (breaking)
- Add `ref_add_connections!` that adds lists of connections to ref for each component (breaking)
- Rename constraint and variable functions to better match PowerModels conventions (breaking)
  - variable_mc_gen_power_setpoint -> variable_mc_generator_power
  - variable_mc_gen_power_setpoint_on_off -> variable_mc_generator_power_on_off
  - constraint_mc_gen_setpoint -> constraint_mc_generator_power
  - constraint_mc_slack_power_balance -> constraint_mc_power_balance_slack
  - constraint_mc_shed_power_balance -> constraint_mc_power_balance_shed
  - constraint_mc_load_power_balance -> constraint_mc_power_balance
  - variable_mc_load_setpoint -> variable_mc_load_power
  - constraint_mc_load_setpoint -> constraint_mc_load_power
- Updates objective function for MLD problem (breaking)
- Add `correct_mc_voltage_angle_differences`, `correct_mc_thermal_limits`
- Removed support for JuMP \< v0.21 (breaking)
- Overrides `_objective_min_fuel_cost_polynomial_linquad` from PowerModels to support arbitrary connections on generators
- Updated solution building functions to automatically parse arbitrarily-sized vectors of variables into solutions
- Changed `Array{...,1}` to `Vector{...}` and `Array{...,2}` to `Matrix{...}`
- Removes phase projection by default, but keeps phase projection for delta connected components for now (breaking)
- Adds `apply_phase_projection_delta!` to project phases of delta connected components
- Fixes grounding logic for generator and solar objects
- Fixes bug in parsing of file paths on windows in redirect/compile dss commands
- Adds multiconductor version of the OPF objective function `objective_mc_min_fuel_cost`
- Update publication in README
- Fixes bug in json parser (typo)
- Fixes bug in dss parser that skips some files with same names as previously parsed files

## v0.9.3

- Fix bug in buscoords parser to support more valid buscoords files
- Fix bug in parse_json(::String) which passed additional unused argument to parse_json(::IO)

## v0.9.2

- Add possibility to give vm_start in variable_mc_voltage for ivr and acr formulations
- Remove leftover code that was causing a warning on import, but was not breaking the implementation
- Add multinetwork MLD and OPF problem specifications

## v0.9.1

- Add support for storage components and mld problem in ACR formulation
- Add storage variables and constraints to OPF, PF, and MLD problems
- Fix bug in `con` mapping in power balance constraints where an array was expected
- Fix bug in `constraint_mc_shed_power_balance` where Gt, Bt were not defined
- Fix bug in the status of virtual branches created from switches where the status of the switch was not taken into account, only the state
- Fix bug in transformer parsing where `status` wasn't being included
- Fix constraint for PV buses, which were attempting to be applied in places where `vm` was not specified
- Fix bug in type enforcement of `map` argument in `transform_solution`
- Fix storage parsing including unit conversion and inclusion of time_elapsed from the root level
- Add wye-wye and delta-wye constraints to LPUBFDiagModel formulation
- Fix typo in ivr formulation line 548, was var(pm, nw, :crg_bus, id) now is var(pm, nw, :cig_bus, id)
- Fix missing / incorrect type of some properties on lines in opendss parser (#290)
- Fix connections-determining-code on solar and storage objects to generator object syntax (#291)
- Refactors Kron reduction and padding transformations out of eng2math into their own transformation functions (#287)
- Add functionality of run_mc_mld_bf to run_mc_mld via multiple dispatch
- Fixes inconsistency of connections on MATHEMATICAL components, in particular, virtual objects (#280)
- Add a transformation remove_all_bounds! that removes all fields ending in _ub and _lb (#278)
- Add missing connections for virtual generator at voltage source
- Fix pu conversion bus voltage bounds and add parsing for vm_pair_lb and vm_pair_ub
- Add CONTRIBUTING.md
- Add pull request template
- Add Pull Request Checks (Github Actions), to inform developers of potentially missing aspects of their PR
- Store references for power balance constraints

## v0.9.0

- Add `instantiate_mc_model` to aid in building JuMP model from ENGINEERING data model
- SDP and SOC relaxations were broken but are fixed again (unit tests added)
- Remove `run_mc_opf_iv`, `run_mc_opf_bf`, `run_mc_opf_bf_lm`, `run_mc_pf_bf`, `run_mc_pf_iv`, these can be accessed by using the correct formulation with `run_mc_opf` and `run_mc_pf`
- Add support for Memento 1.1
- Add support for PowerModels v0.17 (breaking)
- Add support for InfrastructureModels v0.5
- Updates JSON parser to handle enum (`"data_model"` values)
- Adds some commonly used InfrastructureModels and PowerModels functions as exports
- Adds model building functions `add_{component}!` to aid in building simple models for testing (experimental)
- Add `run_mc_model` (adds `ref_add_arcs_transformer!` to ref_extensions, and sets `multiconductor=true` by default) (breaking)
- Rename `ref_add_arcs_trans` -> `ref_add_arcs_transformer` (breaking)
- Update `count_nodes`, now counts source nodes as well, excludes \_virtual objects
- Change \_PMs and \_IMs to \_PM, \_IM, respectively
- Add example for PowerModelsDistribution usage (see Jupyter notebooks in `/examples`)
- Update transformer mathematical model
- Introduce new data models: ENGINEERING, MATHEMATICAL (see data model documentation) (breaking)
- Update DSS parser to be more robust, and parse into new format (breaking)
- Updates DSS paser to parse more options/commands, moves these into `"options"` dict (breaking)
- Updates how dss `like` is applied to better match opendss (almost all properties are copied with like) (breaking)
- Add support for new OpenDSS components (loadshape, xfmrcode, xycurve)
- Add support for JuMP v0.22 (exports `optimizer_with_attributtes` by default)
- Add support for PowerModels v0.16 (breaking)
- Add support for Memento v0.13, v1.0

## v0.8.1

- Update to support JuMP v0.21
- Makes bounds optional, turned on by default (#250)
- Updated transformer data model in the mathematical model (#250)
- Add automatic parsing of lon,lat from buscoords file into PowerModelsDistribution data structure (#245, #249)
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

- Update solution building infrastructure (PowerModels #77) (breaking). The reported solution is now consistent with the variable space of the formulation.
- Moved multi-conductor support from PowerModels into PowerModelsDistribution. (breaking)
- PowerModels.var no longer takes conductor as an argument
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
