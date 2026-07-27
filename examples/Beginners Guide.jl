### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 9787aa73-8ffc-4634-bf0f-b70eee0bf377
using CodeTracking, Revise, PlutoUI

# ╔═╡ a1989876-9301-11eb-0783-83b3aa7abfbc
begin
	using PowerModelsDistribution
	import InfrastructureModels
	import JuMP
	import Ipopt
	import JSON
end

# ╔═╡ 2de210f3-75b4-4c5d-89a7-65557b45df2c
using PowerModelsAnalytics

# ╔═╡ c9a5c344-961b-11eb-0458-a746afcb280c
html"""

<style>
main {
	max-width: 1000px;
}
body {
	overflow-x: hidden;
}
</style>
"""


# ╔═╡ b5928052-9616-11eb-14b3-79770b9929df
md"""
# Introduction to PowerModelsDistribution

This Notebook was designed for the following versions:

- `julia = "~1.6"`
- `PowerModelsDistribution = "~0.12"`
- `PowerModelsAnalytics = "~0.4.1"`

This notebook is a begginer's introduction to PowerModelsDistribution, an optimization-focused Julia library for quasi-steady state power distribution modeling, based on JuMP.jl, and part of the larger [InfrastructureModels.jl](https://github.com/lanl-ansi/InfrastructureModels.jl) ecosystem, which notably includes:

- [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) : Transmission (single-phase positive sequence power networks) optimization
- [GasModels.jl](https://github.com/lanl-ansi/GasModels.jl) : Natural Gas pipeline optimization (includes Steady-state and Transient optimization)
- [WaterModels.jl](https://github.com/lanl-ansi/WaterModels.jl) : Water network steady-state optimization

Details about PowerModelsDistribution.jl can be found in our [PSCC Conference Proceedings paper](https://doi.org/10.1016/j.epsr.2020.106664).
"""

# ╔═╡ 62c14531-357a-4669-90cd-2a186df123eb
md"The following packages are used for notebook features only and do not relate to tutorial content"

# ╔═╡ b953c65d-515f-4334-a4d4-b27af1b0e29a
md"""
This notebook will make use of the following packages in various places
"""

# ╔═╡ f4d7c365-ce18-46b8-8ebb-dd2fbc0e73de
IM = InfrastructureModels

# ╔═╡ 749d062d-72a1-4160-afcf-1cdc27d85c84
md"""
## Case Section

This notebook can apply to different data sets, select a case for examples below from the cases included in the PMD unit testing suite:
"""

# ╔═╡ 7e032130-c565-42d1-93ed-87955e1f2334
begin
	pmd_path = joinpath(dirname(pathof(PowerModelsDistribution)), "..")
	@bind case_file Select([
			joinpath(pmd_path, "test/data/opendss/case3_balanced.dss") => "case3_balanced",
			joinpath(pmd_path, "test/data/opendss/case3_unbalanced.dss") => "case3_unbalanced",
			joinpath(pmd_path, "test/data/opendss/case3_balanced_battery.dss") => "case3_balanced_battery",
			joinpath(pmd_path, "test/data/opendss/case5_phase_drop.dss") => "case5_phase_drop",
			joinpath(pmd_path, "test/data/opendss/ut_trans_2w_yy_oltc.dss") => "ut_trans_2w_yy_oltc",
			joinpath(pmd_path, "test/data/opendss/case3_balanced_battery.dss") => "case3_balanced_battery",
		])
end

# ╔═╡ 67029534-961a-11eb-2b06-210f67607d20
begin
	dss = open(case_file, "r") do f
		join(readlines(f),"\n")
	end

	importing_data_md = """
# Importing Data

PMD supports two input formats, __OpenDSS__ and __JSON__. We strongly recommend OpenDSS for new users, as JSON is intended primarily for data models and results portability between colleagues working on the same problem, and OpenDSS is appropriate for specifying new networks.

Below is an example of an OpenDSS specification for feeder $(case_file):

```dss
$(dss)
```

Data is imported via the `parse_file` command, which we will use further down in the tutorial.
"""

	importing_data_md |> Markdown.parse
end

# ╔═╡ 0b9598a8-9618-11eb-1947-4f98dac7129f
md"""
# Data Models

In PMD, there are two data models, an `ENGINEERING` data model, which is meant to be user facing, and to better reflect the engineering realities of the system, and a `MATHEMATICAL` data model, which reflect the mathematical representation of the system.

Data models are identified by a key in the data dictionary, `"data_model"`, whose values are `ENUM`s:

- $(ENGINEERING)
- $(MATHEMATICAL)

## ENGINEERING data model

Full specification of the `ENGINEERING` data model can be found in our [documentation](https://lanl-ansi.github.io/PowerModelsDistribution.jl/stable/eng-data-model/).

The engineering data model supports several broad categories of data:

- metadata
- node objects
- edge objects
- data objects

"""

# ╔═╡ 6d94d2cd-069b-4e48-84da-ba87c554297b
eng = parse_file(case_file)

# ╔═╡ 480f6685-56a0-4fdd-a975-33cdb3459ef5
"""
### Metadata

Metadata is mostly straight-forward, containing meta information about the feeder, on the parse, and what type of data model is currently being represented.

For `case3_balanced_eng` the following metadata fields are available:

- `name`
- `conductor_ids`
- `settings`
- `files`
- `conductors`
- `data_model`
""" |> Markdown.parse

# ╔═╡ b73eedaf-5809-43b9-90e0-69350d6ee2c4
Dict{String,Any}(k => eng[k] for k in ["settings", "conductor_ids", "files", "name", "data_model"])

# ╔═╡ 4243344a-9cdf-46db-96ee-633b0f5e8a71
md"""
`settings` and `data_model` are the most important metadata, required for solving any type of optimization problem.

`data_model` is self-explanatory, but `settings` requires some explanation. First and foremost, `settings` contains the information needed to calculate the voltage bases for all of the buses in the network. When parsing from a DSS file, these are not explicitly stated, and must be extrapolated from the voltage at the source, and the network must be walked-through, and in a case with transformers the new voltage base is adjusted on the other side of the transformer.

Inside `settings`, scalars can also be set, which might be valuable if, _e.g._, the power or voltage values are being scaled into inappropriate per-unit values.

Finally, `settings` also contains `base_frequency`, which by default is 60 Hz, but can be specified in OpenDSS, and is valuable to know if, for example, you are working with data from Europe, which might have a base frequency of 50 Hz.
"""

# ╔═╡ e95a9fba-b6ad-4339-8219-de4ca04536f4
md"""
### Distribution Assets

These consist of the actual physical assests in the distribution feeder, including the buses, which are the connective nodes on the network graph, lines, which are the fundamental edges, and others. In `case3_balanced_eng` there is the following asset types:

- `bus`
- `line`
- `load`
- `voltage_source`
"""

# ╔═╡ 23a3a107-2fc4-4af1-8eb6-e6a0dc29ddbe
Dict{String,Any}(k => eng[k] for k in ["bus", "line", "load", "voltage_source"])

# ╔═╡ ba7c38a6-34f5-40d3-998e-7ea520001454
md"""
Voltage sources are representations of the substation at which the feeder is connected. By default in dss, there is a default voltage source called `"source"`, which has some default values.
"""

# ╔═╡ 47fdba38-56ae-4a68-88a0-26af4e0a907a
md"""
### Data objects

Data objects are things that affect and/or modify other objects, so linecodes, transformer codes (xfmrcode), or time series data.

In `eng`, only the following data object exists:

- `linecode`
"""

# ╔═╡ 4628c72d-32a2-46bb-929e-153f3851f9c0
Dict{String,Any}(k => eng[k] for k in ["linecode"])

# ╔═╡ ca82db50-ff48-411c-bc97-44f5f33c4ddf
md"""
### Enums

In the `ENGINEERING` data model we make heavy use of a Julia data structure called an `Enum`, or an Enumerated Type. This is a type whose values are enumerated, starting with 0. This has the benefit of being much more readable by the user.

If you are familiar with JuMP, you probably are already used to Enums `TerminationStatusCode` and `ResultStatusCode`, which we import explicitly from JuMP.MOI and export, for easy access by the user when using `using PowerModelsDistribution`.

For example, instead of a switch `state` having the possible values 0 or 1, instead we created an enumerated type `SwitchState`, with values `OPEN` (0) and `CLOSED` (1).

We follow the convention that Enum values are all uppercase.

Enums can be cast into their integer values easily:
"""

# ╔═╡ 4f45cac5-3b27-4771-994b-3a9d81ed61c9
Int(OPEN)

# ╔═╡ fb08875c-1d9e-4744-bec2-7b7cad4320f6
md"Integers can be converted back to Enums just as easily..."

# ╔═╡ 5cf93941-3f88-4bd1-98c4-77b924ddb385
SwitchState(0)

# ╔═╡ 8a8222c3-ef01-4d18-97ff-fa43ad89cfcc
md"The following Enum types exist in PMD"

# ╔═╡ 38dfd21e-f156-45be-b9b3-f6eb5008723e
PowerModelsDistributionEnums

# ╔═╡ 0f73a5be-7ee3-4936-95ed-ebb1b913cf4e
md"and the following enum values exist currently in PMD (excluding those imported from JuMP.MOI):"

# ╔═╡ 254a4ed3-263b-4084-b9af-c97eac2d4ab7
[n for n in names(PowerModelsDistribution) if isa(getproperty(PowerModelsDistribution, n), Enum) && !isa(getproperty(PowerModelsDistribution, n),Union{TerminationStatusCode,ResultStatusCode})]

# ╔═╡ 03e27a4b-cf53-4e92-9e88-d07015937b3e
md"""
Some common examples when you will typically see Enums include:

- `status`, on all components,
- `data_model`, at the root level,
- `dispatchable`, on things like switches, loads, and shunts, which indicate an ability to change their "state", like shedding the load, or opening or closing the switch
- `configuration`, which indicates the connection type, `WYE` or `DELTA`
- `model`, e.g., on loads, which can indicate the type of load, like constant `POWER`, `CURRENT`, `IMPEDANCE`, etc.
"""

# ╔═╡ 3e74097f-55ce-4143-938d-543a59581a72
md"""
### Transformations

Transformations are one of the most powerfull aspects of using the engineering model, because items are more simple and self-contained rather than decomposed, editing before transformation into a mathematical model is significantly more straightforward.

The best example of this is Kron reduction, which is still done by default, where it would be too complicated at the mathematical level, requiring significant changes to the transformer models, for example.

Some simple examples that we commonly use involve settings better OPF bounds:
"""

# ╔═╡ 2c11f944-ff4c-478b-9b6c-fba54bfc0afd
begin
	apply_voltage_bounds!(eng; vm_lb=0.9, vm_ub=1.1)
	apply_voltage_angle_difference_bounds!(eng, 1)
end

# ╔═╡ 427bebe9-5a0d-41b0-ace6-6622400e136c
md"""
Some other valuable transformations built into PMD are:

- `make_lossless` (will strip loss models on engineering assets that contain them, e.g., voltage sources or switches)
- `remove_all_bounds` (will remove all bounds, e.g., those parsed in from the raw dss file)
"""

# ╔═╡ 284a471c-5b97-4753-b5be-8896dc096657
"""
```julia
$(@code_string remove_all_bounds!(eng))
```
""" |> Markdown.parse

# ╔═╡ 20c563f0-9304-11eb-16b4-6720b073e911
md"""
## MATHEMATICAL data model

The mathematical data model is a transformation of the engineering components into ones which we can more easily represent in the optimization model.
"""

# ╔═╡ e70c3f4a-9302-11eb-1176-0925031bccc0
math = transform_data_model(eng)

# ╔═╡ 9a65d284-8558-4ddf-9139-e88e3a3625d9
md"""
The mathematical model can also be loaded directly via `parse_file`:

```julia
parse_file(case_file; data_model=MATHEMATICAL)
```

In some cases, the transformations are straight-forward, 1-to-1 type of conversions, where we convert to more optimization-friendly fields and units, like with lines -> branches, or loads -> loads.

But, some other components' transformations are less obvious, like voltage sources.
"""

# ╔═╡ 7edf07c2-208d-47d4-914e-08c09bdfd7f6
@bind source_id_select Select(["$type.$name" for type in pmd_eng_asset_types for name in keys(get(eng, type, Dict()))])

# ╔═╡ 36f01f08-0c19-4a28-a2eb-11a11dec9897
filter(
	x->!isempty(x.second),
	Dict(
		type => Dict(
				name => obj
				for (name,obj) in get(math, type, Dict())
					if source_id_select == obj["source_id"]
				) for type in pmd_math_asset_types
		)
)

# ╔═╡ a3e10f7e-c193-4c02-9fef-e23bde112350
md"""
Note: All objects have the field `source_id`, that is meant to indicate where a `MATHEMATICAL` object originated from with the `ENGINEERING` model.
"""

# ╔═╡ 1cf59273-3f50-41bd-b93c-c5cf5f3cd124
md"""
The reason for this is that some more complex objects can be decomposed into multiple mathematical objects. In the case of this voltage source, there is a non-zero source impedance, and rather than creating an entirely new mathematical object, we can decompose it into a generator with unlimited power bounds, and into an impedance branch, with a connecting bus.
"""

# ╔═╡ 7e36b7b8-198c-4b9e-8b62-d1b481f09eeb
"""
### Componnent ID Format

In the engineering model, all component ids can be arbitrary strings, making it easier to navigate a feeder, but in the mathematical model we use only integers (Strings for dict keys, to maintain JSON compatibility, and Ints for ids within the data properties). For example:
""" |> Markdown.parse

# ╔═╡ b717dab6-30fd-4e4b-b6f8-29ffcff80131
bus_keys = keys(math["bus"])

# ╔═╡ def53588-f511-4dda-8aa9-941d1aa994f7
bus_ids = [bus["bus_i"] for (i, bus) in math["bus"]]

# ╔═╡ f1866650-ca90-4f4b-81ef-1cc3202dc240
md"""

### Additional Metadata

- `map`
- `bus_lookup`
- `basekv`
- `baseMVA`
- `is_projected`
- `per_unit`
- `is_kron_reduced`

Two particular items classified as metadata that are key for understanding how the data model maps between the engineering and mathematical models are `bus_lookup`, which maps bus names into their new integer ids, and `map`, which is an ordered list of actions that were taken to map engineering to mathematical model, so that we can map the solutions back up to the correct components in the engineering model.

"""

# ╔═╡ 27ec55f0-77c2-4cec-b505-976fd86f1004
math["bus_lookup"]

# ╔═╡ 0ef86c0a-6d65-44a3-9804-676d4cc904c3
math["map"]

# ╔═╡ e45adf3e-5b84-4d14-9840-40ba2fbc1573
"""
### Asset Types

For mathematical models we initially adopted the PowerModels data model, which was originally designed based on the Matpower package. We try to largely maintain this parity with PowerModels even though it is no longer a dependency, which is why there there remains some inconsistency in property names compared to the engineering model. The current components in the math model are:

- `bus`
- `load`
- `shunt`
- `gen`
- `branch`
- `transformer`
- `switch`
- `storage`

The one notable divergence is the existance of `transformer`s. In distribution models, phase unbalanced transformers are much more complex than the typical two-winding Pi-branch model commonly utilized in transmission grids.

In the PMD math model, transformers are two-winding lossless transformers that can be either wye-wye or wye-delta connected.
""" |> Markdown.parse

# ╔═╡ 1b17a770-a465-4ba9-a21d-b3d4514ad696
md"""
## `import_all`: Keeping raw DSS properties

The DSS parser included in PMD should parse all raw DSS objects and properties, but by default it will only keep around internal data fields when converting to the `ENGINEERING` or `MATHEMATICAL` data models.

To keep raw DSS properties, use the `import_all` keyword argument when parsing data:
"""

# ╔═╡ 29afef9d-29fe-46fe-9321-a33a488b2ba0
eng_import_all = parse_file(case_file; import_all=true)

# ╔═╡ 2aa3d0ca-115c-476a-be94-f33a2343ea57
md"""
You will notice a couple of extra data structures.

The first you may notice is "dss_options", which will list options that were set at the top-level of the DSS input
"""

# ╔═╡ 81930c47-7f02-4eca-a4bc-df6369691eec
get(eng_import_all, "dss_options", Dict())

# ╔═╡ 394b5c41-1fcc-497d-b598-356e0760da62
md"""
Also, within each asset dictionary will be a `dss` dictionary that will contain the raw dss input about the object from which the asset was derived.
"""

# ╔═╡ 3e2ebbfb-08dd-4099-9367-86c982017e9c
first(eng_import_all["line"]).second["dss"]

# ╔═╡ f090cadb-c87b-4e36-bc28-79736479933b
md"""
This information gets carried through to the `MATHEMATICAL` model as well...
"""

# ╔═╡ 06643acf-6f8e-4bda-a66b-aab42fe4690b
math_import_all = transform_data_model(eng_import_all)

# ╔═╡ 7658ea1b-12e7-475e-898f-d8b31a39665a
first(math_import_all["branch"]).second["dss"]

# ╔═╡ 905bfb0a-2fb5-4b0c-bc62-49bd2888ad30
md"""
# Optimization in PMD

Solving optimization problems in PowerModelsDistribution will feel very familiar to those who use PowerModels.jl for transmission grids (positive sequence representable networks).

Full optimization problems consist of a data model, a mathematical formulation, and a problem specification in which the variables, constraints and objectives are defined.

Additional details can be found in our documentation about the [mathematical problem specifications](https://lanl-ansi.github.io/PowerModelsDistribution.jl/stable/math-model/) and [formulations](https://lanl-ansi.github.io/PowerModelsDistribution.jl/stable/formulation-details/).
"""

# ╔═╡ 3d468a08-c62c-42cb-a096-8ba382bb91e5
md"""
## Formulations

Formulations in PMD are represented by Julia Types, and have a clearly defined hierarchy, starting with our base abstact unbalanced (i.e., multiconductor) PowerModel: `AbstractUnbalancedPowerModel`.

Some useful abbreviations / acronyms in PowerModel and function names:

- `MC`/`_mc_` : multi-conductor, to differentiate from PowerModels.jl name, indicating applies to multiconductor / phase unbalanced / distribution problems.
- `U` : unbalanced
- `BF` : branch flow
- `ACP` : AC polar
- `ACR` : AC rectangular
- `IVR` : IV reectangular
- `LP` : linear program
- `SDP` : semi-definite program
- `SOC` : second-order cone
- `KCL` : Kirchoff's Current Law
- `MX` : matrix

### Non-convex Formulations

- `ACPUPowerModel` : Complex Power-Voltage space polar multiconductor form (bus injection model)
- `ACRUPowerModel` : Complex Power-Voltage space rectangular multiconductor form (bus injection model)
- `IVRUPowerModel` : Complex Current-Voltage space rectangular multiconductor form (bus injection model)

### Linear/Quadratic Formulations

- `LPUBFDiagPowerModel` / `LinDist3Flow` : Diagonal matrix formulation of DistFlow equations (unbalanced branch flow model)
- `NFAUPowerModel <: AbstractUnbalancedActivePowerModel` : Linear, Active-power-only multiconductor form (bus injection model)
- `DCPUPowerModel` : DC polar multiconductor form (bus injection model)

### Semi-definite formulations

- `SDPUBFMCPowerModel` : SDP multiconductor form (unbalanced branch flow model)
- `SDPUBFKCLMXMCPowerModel` : SDP with Matrix KCL constraint multiconductor form (unbalanced branch flow model)

### Second-Order Cone formulations

- `SOCUBFNLPMCPowerModel` : SOC-representable with non-linear ... multiconductor form (unbalanced branch flow model)
- `SOCUBFConicMCPowerModel` : SOC-representable with conic multiconductor form (unbalanced branch flow model)

A more detailed description of the type heirarchies can be found in our [documentation](https://lanl-ansi.github.io/PowerModelsDistribution.jl/stable/formulations/).
"""

# ╔═╡ 114f8b4c-3296-428a-bf72-f5efe6b21b47
filter(x -> endswith(String(x), "PowerModel"), names(PowerModelsDistribution))

# ╔═╡ 0d510b81-6f58-4790-be85-82c9b215322d
"""
## Problem Specifications

Some additional helpful abbreviations for problem specifications:

- `mn` - Multinetwork, i.e., time series problem
- `uc` - unit-commitment, i.e., full load shed only, no partial shed

In PMD, there are two primary problem specifications,

- Optimal Power Flow (OPF) `solve_mc_opf`
  - There is a sub-problem for OPF called On-load Tap Changing (OLTC)
- Maximal Load Delivery (MLD) `solve_mc_mld`
  - There is also a "simple" version, but it will not work in conjunction with a switching optimization (`solve_mc_mld_simple`)
  - There is a caveat that the currently included MLD problem features continuous shedding of individual loads, which is not realistic for real-world distribution operations. While in transmission problems, one may assume some continuous shedding of loads, in real distribution grids, most loads can only be shed by using switches to isolate load blocks.

Recently, with the addition of idealized switches, we have also added a test version of a switching problem (OSW), which has variables for switch states of dispatchable switches, and an additional term in the objective to discourage state changes. This spec is experimental, but ongoing research in this area are expected to yield updates in this problem area.

We also have a debugging problem spec, `solve_mc_opf_pbs`, which will install slacks at every bus, which can be helpful in determining where the issue is in the network.

Power flow (`solve_mc_pf`) is also included, but if should be noted that this uses the same mechanism to solve as all of our optimization problems, and is therefore not efficient or quick, like a Newton-Raphson or Backwards/Forwards method for power flow solving might be.

Below is the primary example for the OPF problem:

```julia
$(@code_string build_mc_opf(instantiate_mc_model(eng, ACPUPowerModel, build_mc_opf)))
```
""" |> Markdown.parse

# ╔═╡ cf3614bc-4c35-4b4d-9ec9-f43307a04c48
filter(x->startswith(String(x), "solve_"), names(PowerModelsDistribution))

# ╔═╡ 6335835e-2b5e-45a0-a5f6-4d97328694bd
md"""
# Examples

Let's start with the AC polar formulation, and solve the OPF problem using Ipopt:
"""

# ╔═╡ bea3a93e-a827-473a-9e19-348ad85dfa25
eng_result = solve_mc_opf(eng, ACPUPowerModel, Ipopt.Optimizer)

# ╔═╡ 86d66ed7-89b1-461e-9487-0e190e6fd02e
md"""
We have designed the `solution` dictionary to be as verbose as possible, including all variables contained in the problem automatically, in the same order and format as the input data.

This means that if an engineering data model is provided, the results will return in the same units and format as that model, unless otherwise instructed...
"""

# ╔═╡ ce39193d-833b-40aa-b532-f9db22f03652
eng_result_pu = solve_mc_opf(eng, ACPUPowerModel, Ipopt.Optimizer; make_si=false)

# ╔═╡ cb632c77-5c8b-4fbb-a867-127a22ce54ad
eng_result_pu["solution"]["bus"]["sourcebus"]

# ╔═╡ b3bd6e16-28dc-44b7-b8ad-1478cdfefecc
eng_result["solution"]["bus"]["sourcebus"]

# ╔═╡ a0e6612f-516e-4be4-b3ae-00c612b981e2
md"""
What if `vm` and `va` are desired, but those variables are not in the model, for example, in `ACRMCPowerModel`?
"""

# ╔═╡ 77899c35-2609-4aef-bcd0-178076237d1d
eng_result_acr2acp = solve_mc_opf(eng, ACRUPowerModel, Ipopt.Optimizer; solution_processors=[sol_data_model!])

# ╔═╡ eb3d6861-0d91-4ae1-9580-7b3182cefef1
eng_result_acr2acp["solution"]["bus"]["sourcebus"]

# ╔═╡ 1b6fd52f-0cdf-4f4f-aa38-31cd8aa9ca7d
"""
```julia
$(@code_string PowerModelsDistribution._sol_data_model_acr!(eng_result))
```
""" |> Markdown.parse

# ╔═╡ e39e89d8-7b3d-4c2d-bbe3-8de202e91c5f
md"""
It is also possible to optimize using the `MATHEMATICAL` data model directly, but it will output results in the same format as the model it is provided:
"""

# ╔═╡ 493d1dfd-719c-4f03-b660-6b0e9e50a222
math_result = solve_mc_opf(math, ACPUPowerModel, Ipopt.Optimizer)

# ╔═╡ b6700c39-8dff-4bcb-b6fa-c67b2ddf5163
md"""
However, if your `MATHEMATICAL` data model contains the `map`, it is possible to manually convert the solution back into the `ENGINEERING` structure...
"""

# ╔═╡ f0cb5faf-0471-4651-9d09-61fe6e5e0fb5
transform_solution(math_result["solution"], math)

# ╔═╡ 89f6282d-9898-4547-a121-028fcdf4f876
"""
# PMD Internals for Specification / Formulation Builders

A problem is formally created using `instantiate_mc_model(data, form, prob)`, and outputs a Julia Struct:

```julia
$(@code_string IM.InitializeInfrastructureModel(NFAUPowerModel, eng, PowerModelsDistribution._pmd_global_keys, pmd_it_sym))
```

The following helper functions are here to help you navigate through the mathematical model.

- `ref`
- `var`
- `con`
- `ids`

""" |> Markdown.parse

# ╔═╡ 7303e40a-56f8-4035-8537-e2eed7c16b8b
pm = instantiate_mc_model(eng, NFAUPowerModel, build_mc_opf)

# ╔═╡ c3dcee8a-edf1-43bb-8733-f2798b29d57a
md"As you can see, the following property names exist on the PowerModel struct:"

# ╔═╡ d7522533-43b9-4b3c-a8e5-1aa4f39700bb
propertynames(pm)

# ╔═╡ f4fb1464-8d33-4265-ab40-8fb72d533466
md"Using the `ref` helper function, for example, we could get all branches in the model..."

# ╔═╡ 213e9ea4-c800-47bc-a020-4eacc334e4fc
math_branches = ref(pm,:branch)

# ╔═╡ e32b7e5e-1d99-482f-b406-460a18cde912
md"Or, we could get only the branch ids by using `ids`..."

# ╔═╡ 6b7e5199-61ed-40b0-a658-1a490298dc6d
branch_ids = ids(pm,:branch)

# ╔═╡ fb48ddde-128a-4a32-ab5f-084561108539
md"""
# Example: Upgrading MLD to use Load Blocks

As mentioned in the section above on Problem Specifications, the MLD problem bundled in PMD represents all loads as individually sheddable, which is not accurate to distribution feeders, where it is unlikely that loads would be sheddable by themselves. Instead, usually loads are only sheddable as a whole block, by opening switches to isolate them.

In this example I am going to get us closer to that more realistic problem by creating indicator variables for loads that apply to the whole load block, instead of variables for each load individually.

To achieve this, first we must be able to calculate the possible load blocks, which we can do with `identify_load_blocks`, which will return all sets of buses that can be isolated with switches...
"""

# ╔═╡ 92fd9826-17a4-44e2-ae8b-8d96d6f21130
identify_load_blocks(math)

# ╔═╡ a7037bcd-e65d-40d4-8070-7fd73811ffad
md"""
Then we should add the loads in each block to a `ref` for easy lookup when we are building our model...
"""

# ╔═╡ 53e1237b-c380-4506-9924-8b7298d0b38a
""
function _ref_add_load_blocks!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
	ref[:load_blocks] = Dict{Int,Set}(i => block for (i,block) in enumerate(identify_load_blocks(data)))

	load_block_map = Dict{Int,Int}()
	for (l,load) in get(data, "load", Dict())
		for (b,block) in ref[:load_blocks]
			if load["load_bus"] in block
				load_block_map[parse(Int,l)] = b
			end
		end
	end
	ref[:load_block_map] = load_block_map

	load_block_switches = Dict{Int,Vector{Int}}(b => Vector{Int}([]) for (b, block) in ref[:load_blocks])
	for (b,block) in ref[:load_blocks]
		for (s,switch) in get(data, "switch", Dict())
			if switch["f_bus"] in block || switch["t_bus"] in block
				if switch["dispatchable"] == 1 && switch["status"] == 1
					push!(load_block_switches[b], parse(Int,s))
				end
			end
		end
	end
	ref[:load_block_switches] = load_block_switches
end

# ╔═╡ 60485fa1-7b4d-4be1-a0c4-50b45a840a66
md"""
Because of a recent upgrade to support multi-infrasture models, we now want to use `apply_pmd!` to help us apply this ref, which will help us apply things to the correct data structure and to each subnetwork, if applicable.
"""

# ╔═╡ 8e2620f3-894f-43ff-9704-637225a6dda3
""
function ref_add_load_blocks!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    apply_pmd!(_ref_add_load_blocks!, ref, data; apply_to_subnetworks=true)
end

# ╔═╡ 516bc915-946d-4216-b5d0-b0bde5db6546
md"""
We will demonstrate how to apply this `add_ref_load_blocks!` later, but in the next steps we will assume we already have these added refs avaiable to us.

Next we need to add the new indicator variables. An indicator variable is a variable z ∈ [0,1] that we can use to shed the loads. This variable gets applied to the real and reactive load power values, so that the load can be dynamically shed in the algorithm.

Because we need to shed whole blocks at a time, there should only be one indicator variable for each block.
"""

# ╔═╡ 87bf0ddd-5678-4aea-9024-26c0a2db7521
"create variables for demand status by load block"
function variable_mc_load_block_indicator(pm::AbstractUnbalancedPowerModel; nw::Int=IM.nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        z_demand = var(pm, nw)[:z_demand_blocks] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :load_blocks)], base_name="$(nw)_z_demand",
            lower_bound = 0,
            upper_bound = 1,
            start = 1.0
        )
    else
        z_demand = var(pm, nw)[:z_demand_blocks] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :load_blocks)], base_name="$(nw)_z_demand",
            binary = true,
            start = 1
        )
    end

    load_block_map = ref(pm, nw, :load_block_map)

    var(pm, nw)[:z_demand] = Dict(l => z_demand[load_block_map[l]] for l in ids(pm, nw, :load))

    # expressions for pd and qd
    pd = var(pm, nw)[:pd] = Dict(i => var(pm, nw)[:z_demand][i].*ref(pm, nw, :load, i)["pd"] for i in ids(pm, nw, :load))
    qd = var(pm, nw)[:qd] = Dict(i => var(pm, nw)[:z_demand][i].*ref(pm, nw, :load, i)["qd"] for i in ids(pm, nw, :load))

    report && IM.sol_component_value(pm, pmd_it_sym, nw, :load, :status, ids(pm, nw, :load), var(pm, nw)[:z_demand])
    report && IM.sol_component_value(pm, pmd_it_sym, nw, :load, :pd, ids(pm, nw, :load), pd)
    report && IM.sol_component_value(pm, pmd_it_sym, nw, :load, :qd, ids(pm, nw, :load), qd)
end

# ╔═╡ f9fc4c70-f091-4ecd-ac70-c1dfcd8a6f93
md"""
The last three lines are how the pd and qd variables get added to the solution.

Finally, we need to update the problem to use this new variable. Lets use the "simple" mld problem as a starting point.

First we need the problem definition (the builder).

Note that in this example I am using a Unbalanced Branch Flow formulation, which obviously has some different constraints in the branch section than would be used with, e.g., the NLP formulations.
"""

# ╔═╡ 33dd7091-5576-4375-b5d9-d9a738656ddb
"Multinetwork load shedding problem for Branch Flow model"
function build_mc_mld_simple_loadblock(pm::AbstractUBFModels)
	variable_mc_bus_voltage(pm)

    variable_mc_branch_power(pm)
	variable_mc_branch_current(pm)
    variable_mc_switch_power(pm)
    variable_mc_transformer_power(pm)

    variable_mc_generator_power(pm)

    variable_mc_load_block_indicator(pm; relax=true)
    variable_mc_shunt_indicator(pm; relax=true)
	variable_mc_storage_power_mi(pm; relax=true)

   	constraint_mc_model_current(pm)

    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    for i in ids(pm, :gen)
        constraint_mc_generator_power(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_mc_power_balance_shed(pm, i)
    end

    for i in ids(pm, :storage)
        constraint_storage_state(pm, i)
        constraint_storage_complementarity_mi(pm, i)
        constraint_mc_storage_losses(pm, i)
        constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_power_losses(pm, i)
        constraint_mc_model_voltage_magnitude_difference(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
        constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    objective_mc_min_load_setpoint_delta_simple(pm)
end

# ╔═╡ 1b18e8ce-f403-4cb6-b833-083e8b9ef5bb
md"""
Next, we need a way to call this problem to solve it:
"""

# ╔═╡ 5369925c-3e0f-43ca-a4d8-c801a90653ff
""
function solve_mc_mld_simple_loadblock(data::Dict{String,<:Any}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mc_mld_simple_loadblock; ref_extensions=[ref_add_load_blocks!], kwargs...)
end

# ╔═╡ 51af76d7-2c62-422d-98d1-b199aba84263
md"""
Note in particular the addition of the `ref_extensions` keyword argument, which takes a vector of function references. This is how we add our custom ref extension `ref_add_load_blocks!`
"""

# ╔═╡ 12ed532a-164c-4a92-a7d1-fc4dd53c3c57
mld_result = solve_mc_mld_simple_loadblock(eng, LPUBFDiagPowerModel, Ipopt.Optimizer)

# ╔═╡ 3207f33b-00e5-4691-a9c9-9461c974b3b0
md"""
But we can see that no loads get shed in the case, because there is no contingency applied to this feeder

What if we instead apply a contingency where power delivery is disabled on one phase, what happens?
"""

# ╔═╡ 8eb99a9b-ee92-4aad-bbfc-c9242dc916cb
begin
	eng_vs_disabled = deepcopy(eng)
	eng_vs_disabled["voltage_source"]["source"]["pg_ub"] = [Inf, Inf, 0]
	eng_vs_disabled["voltage_source"]["source"]["qg_ub"] = [Inf, Inf, 0]
end

# ╔═╡ f70780ef-5ad1-4f69-bb2e-41de737b1bef
mld_vs_disabled_result = solve_mc_mld_simple_loadblock(eng_vs_disabled, LPUBFDiagPowerModel, Ipopt.Optimizer)

# ╔═╡ ab146273-0cf6-4815-b996-dc2e8fc185db
md"""
Because loads are tied together, we must shed all of the load in this block, even though only one phase could not deliver power.

This is a very simplistic example, and therefore the results may not seem interesting in themselves.

Finally, I want to note that in this example, status is not == 0, even though we might expect it to be. This is because we are using the "relaxed" version of the indicator constraints, which will often be not quite zero or one, even when we might expect them to be, especially in the case of the "simple" mld problem.

To guarantee 0 or 1, which is the most realistic for distribution feeders, we should use the unrelaxed indicator variables, but this will require using a different solver that can support mixed-integer variables. In the case of the problem I created above, Juniper, using Cbc and Ipopt might be a good option, but e.g., Gurobi would be better.
"""

# ╔═╡ 3f03084d-f8cc-424e-99a6-78b39ece17ef
md"""
# Multinetwork Data and Problems

Next we will cover multinetwork problems, e.g. time series OPF.

PMD has a lot of tools for multinetwork problems, but really only one constraint that is inherently multinetwork, storage state.

Let's start with constructing a multinetwork data structure, and exploring it.

First, for multinetwork data to be created automatically, we need `time_series` data. If we have chosen one of the data sets above that contains `time_series`, this should have entries...
"""

# ╔═╡ 19e88e22-44e1-4139-bb5b-1247a271e34a
get(eng, "time_series", Dict())

# ╔═╡ 7af68b51-a1b3-4ed8-95bc-43b2523c298c
md"""
`time_series` is one of our "data" objects, in that it does not represent actual assets on the feeder, but represents information about one or more assets. `linecode`s are the most often encountered and most familiar data objects in feeder data.

Each `time_series` object will have "time", "values", and "replace", at a minimum, with optional "source_id", to help find the orginating object, and "offset", for which the intention if for it to add an offset to "time". "offset" is not yet supported, but datetime strings or floats in units of hours in "time" are supported.

"replace" indicates whether the values in "values" will replace the property they are assigned to, if `false`, values will be multiplied by the base value.

To apply a `time_series` object to a property, it must be specified within an asset's specification:
"""

# ╔═╡ 3d769b34-b327-4304-af0a-1b995a4d6e4c
filter(x->!isempty(x.second), Dict(type => Dict(name => obj for (name, obj) in get(eng, type, Dict()) if haskey(obj, "time_series")) for type in pmd_eng_asset_types))

# ╔═╡ 9bb9cbf5-67ef-487b-88a9-830f15378537
md"""
If a input data set with timeseries data has been selected then one should find objects that have their own `"time_series"` dictionary, inside of which the keys are properties to be replaced, and values are references to root-level `time_series` objects.

OpenDSS has multiple ways to specify time series data, most usually through `LoadShapes`, which are specified on `Load` objects via the `daily` or `yearly` properties most often.

By default PMD will parse the `daily` time series data, but you can specify this at parse:

```julia
parse_file(case_file; time_series="yearly")
```

If anything other than `"daily"` is chosen, it might be necessary to adjust `time_elapsed`, which will be discussed below.

If you have `time_series` data specified correctly, building a multinetwork is straightforward:
"""

# ╔═╡ 4bbc5205-4a3e-4e3d-88d6-219195ad8a0f
mn_eng = make_multinetwork(eng)

# ╔═╡ ba8e830f-3c3b-45bc-be25-c440cc736a02
md"""
It is also possible to load from a file directly into a multinetwork data structure:

```julia
parse_file(case_file; multinetwork=true)
```

This transformation changes the root-level of the data model pretty drastically.
"""

# ╔═╡ ea462d05-6325-44b3-84fc-bf34e032bb01
keys(mn_eng)

# ╔═╡ 565385e3-ccb4-4ac7-a0a5-12c2744ceed5
md"Compared to a non-multinetwork structure..."

# ╔═╡ fb05c680-d4ba-46a6-a696-7536aaa0dd57
keys(eng)

# ╔═╡ 5cd2d702-f5fa-44f3-9739-8eb75b0b6472
md"""
However, what is really happening here is that only some information is needed at the root-level, and some information is paired directly with subnetworks (what we can a timestep in the multinetwork structure).

At the top level, we really need `"data_model"`, `"nw"`, and `"multinetwork"`. Even `"mn_lookup"` is only useful for one of our helper functions for manually reorganizing subnetworks `sort_multinetwork!`

The subnetworks live inside `"nw"`:
"""

# ╔═╡ 85d52ceb-2b4f-4c22-880a-290c1a5737a1
mn_eng["nw"]

# ╔═╡ 20f6dfed-5d1a-481c-92f4-440de0f2abca
md"""
Inside `nw`, subnetworks are organized by string integers, corresponding to the "time step". This makes iterating through them consistent...
"""

# ╔═╡ ff6dfa9a-d55e-4bbe-bf41-357a417784a5
sort([parse(Int, n) for n in keys(mn_eng["nw"])])

# ╔═╡ 74a403df-0f94-4a1a-af69-cd8dc5bef430
md"""
Each subnetwork contains all of the information we should need for variables, constraints, etc. at that time step:
"""

# ╔═╡ 28f75060-d0be-47f2-bea9-377d5df89cf6
first(mn_eng["nw"]).second

# ╔═╡ 663e32f5-4b2e-4255-9def-e1c7d383c36d
md"""
The reason for so much duplication of data is that the vision for multinetworks was never just for time series data, but that it could be used more generally in creative ways, like pairing two topologically different networks together with custom problem specifications.

This use case means that things that you might expect to stay the same between time steps, like `conductor_ids`, `linecodes`, `settings`, and the overall topology, could be drastically different, and therefore should be replicated for each subnetwork.

This also explains the origin of the name "multinetwork", in case the standard use case made its name confusing.

As an interesting note, `make_multinetwork` will always return a multinetwork structure, even if there is no `time_series` data, with a single subnetwork with key `"0"`.

Transforming into the `MATHEMATICAL` data model from a multinetwork `ENGINEERING` data model is also the same as for single network data:
"""

# ╔═╡ d3e6486c-4133-404b-91c1-7ac18ec8388f
mn_math = transform_data_model(mn_eng)

# ╔═╡ 8c34cd7f-52c4-440f-a226-afc27bf24799
md"""
Some of the root-level keys will be slightly different, but otherwise you will see what you are already familiar with in the single network `MATHEMATICAL` data model.

It is also possible to transform directly into a multinetwork `MATHEMATICAL` data model from a single network `ENGINEERING` data model:

```julia
mn_math = transform_data_model(eng; multinetwork=true)
```

The one major caveat with these automatic generations of multinetwork data structures is that it **must** be performed before converting to the `MATHEMATICAL` data model. This is because, since we support replacing arbitrary fields with `time_series` data, it is impossible to work out the conversions within the `time_series` objects.

If you have a `MATHEMATICAL` data model and want to convert it to a multinetwork, this is supported, but you must already have a special construction of the `time_series` object that matches the format expected by InfrastructureModels, which has its own, more general, `make_multinetwork` function.

Multinetworks have two key helper functions:

- `sort_multinetwork!`
- `set_time_elapsed!`

The first accepts a Vector of `time` values, which it will use to manually re-sort the subnetworks, and the second accepts either a vector of time deltas or a single time delta, and replaces the `time_elapsed` property within all the subnetworks.

`time_elapsed` a value in hours that indicates how long each time step duration is, which is needed for calculating storage losses.
"""

# ╔═╡ 3562dcc2-b6d7-4bfa-bae7-7a5e992cb422
set_time_elapsed!(mn_eng, 0.5)

# ╔═╡ 805b06be-7ed9-4330-9f1e-a584c4508862
first(mn_eng["nw"]).second["time_elapsed"]

# ╔═╡ e56a2563-45c7-40fd-80ad-b86d7e1c1e44
"""
Solving multinetworks is not anymore difficult than single network cases, but a special problem specification must be used that is multinetwork-aware...

```julia
$(@code_string build_mn_mc_opf(instantiate_mc_model(mn_eng, ACPUPowerModel, build_mn_mc_opf; multinetwork=true)))
```
""" |> Markdown.parse

# ╔═╡ 49808060-ad86-491e-aed2-c69f203d4a0e
md"""
Note that the standard OPF problem loops over each subnetwork (not necessarily in order), with the keyword argument `nw=n`.

Currently the only build-in asset that truly has multinetwork constraints is storage, where you can see the `constraint_storage_state` being called with two nw ids near the bottom of the above specification.

Knowing this, solving a multinetwork OPF problem is straightforward:
"""

# ╔═╡ 59abd602-30de-40af-bce5-7800d11a7108
mn_result = solve_mn_mc_opf(mn_eng, ACRUPowerModel, Ipopt.Optimizer)

# ╔═╡ f40614c0-5aac-491b-bf5b-8aa9d328f640
md"""
# Merging Solution with Data

It is possible to merge your solutions with your data structures, which will make transporting and/or visualizing data easier. This helper function from InfrastructureModels allows you to merge two *nested* dictionaries together:
"""

# ╔═╡ 40e60c21-76fe-44ce-9a33-988030e8d863
eng_copy = deepcopy(eng)

# ╔═╡ 14f71709-35f9-4eb7-8505-a2138b6fb6f5
update_data!(eng_copy, eng_result["solution"])

# ╔═╡ 09f7c681-9654-466d-a1cb-e09eebad62fd
first(eng_copy["bus"]).second

# ╔═╡ 0362394c-9ac9-4978-9b41-bc325464c204
md"""
This will work on both data model types
"""

# ╔═╡ 8361deb2-3baa-4031-927f-e3c034bd7eb9
begin
	math_copy = deepcopy(math)
	update_data!(math_copy, math_result["solution"])
end

# ╔═╡ ba259614-bead-40c4-927f-5f8ad3eb726a
first(math_copy["bus"]).second

# ╔═╡ dee28964-882f-4bad-9368-762f7af99562
md"""
# Exporting and Importing PMD Data Structures

It is possible to export our data structures to JSON, but you may have noticed several items that are not strictly JSON compatible, like Matrix, Enum, Symbol, Inf and NaN. Because most users default to JSON.print to export data structures, we have chosen to create a data model correction helper function that will attempt to fix data structures. For straightforward cases this has shown to work well, but may be fragile in its implementation.
"""

# ╔═╡ 2f10d9b4-c237-4305-ac05-ffe569d70436
begin
	io = PipeBuffer()
	JSON.print(io, eng, 2)
end

# ╔═╡ 26e7b26e-3619-40d0-aef8-0a070af98cbb
raw_from_json = JSON.parse(io)

# ╔═╡ 0c39a147-8587-457e-9a2b-ea38e2bf160a
begin
	parsed_from_json = deepcopy(raw_from_json)
	correct_json_import!(parsed_from_json)
	parsed_from_json
end

# ╔═╡ ec94eeb4-e425-452d-8e46-6ed7f876f463
md"""
This can be easily achieved via `parse_file`:

```julia
parse_file(json_file)
```
"""

# ╔═╡ 292e7528-cf57-4429-995a-571d567b8560
md"""
# Experimental Network Plots with PowerModelsAnalytics

It is possible to quickly create some plots of power networks using `PowerModelsAnalytics.plot_network!`.

Originally we created the plotting functionality in PowerModelsAnalytics primarily for debugging purposes, to look for topological errors, check for errors with load shedding, etc., and had based it on Plots.jl, which is a very popular Julia plotting tool. Unfortunately, plotting graph networks with a lot of nodes is very slow, and we discovered Vega.jl, which is a interface to the Vega visualization grammar.

For the most part, simple plots can be easily achieved with `plot_network!` (best used for Pluto notebooks to produce the plot in the notebook), or `plot_network`, which will return the LightGraphs-based graph representation of the network.
"""

# ╔═╡ 687863f6-1aed-4039-b54e-d9d3a683631b
plot_network!(eng_copy)

# ╔═╡ b448fae8-56d4-40a0-95fa-d502f8a37a40
md"""
While plots won't be publication ready, with some knowledge of Vega, and some tweaking of the plot specifications, it should be possible to produce some nice outputs.

Under the hood, PMA uses Networkx to automatically layout the graph, but `use_coordinates=true` can be used to use any buscoords included in the data set.
"""

# ╔═╡ 53f3cd12-05a9-4495-937a-781ce3907174
try
	plot_network!(eng_copy; use_coordinates=true)
catch
	md"**no buscoords exist for this case_file**"
end

# ╔═╡ 232a4bde-3f85-46c4-9e57-ffc64ecafbb6
md"""
# Development

PowerModelsDistribution is subject to active, ongoing development, and is used internally by various high-profile projects, making its improvement and maintanence high priority.

If you find bugs while using PMD, we encourage you to submit bug reports on our [GitHub Issues](https://github.com/lanl-ansi/PowerModelsDistribution.jl/issues).

If you have questions about using PMD, [JuliaLang Discourse](https://discourse.julialang.org/) is a great place, which several of our developers regularly watch, particularly in the [Optimization Category](https://discourse.julialang.org/c/domain/opt/13).

We always welcome [Pull Requests](https://github.com/lanl-ansi/PowerModelsDistribution.jl/pulls) for new features and bug fixes as well.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CodeTracking = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
InfrastructureModels = "2030c09a-7f63-5d83-885d-db604e0e9cc0"
Ipopt = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PowerModelsAnalytics = "7425d11c-6b70-11e9-2bc6-9bcd3191d4b9"
PowerModelsDistribution = "d7431456-977f-11e9-2de3-97ff7677985e"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"

[compat]
CodeTracking = "~1.2.2"
InfrastructureModels = "~0.7.6"
Ipopt = "~1.2.0"
JSON = "~0.21.3"
JuMP = "~1.9.0"
PlutoUI = "~0.7.50"
PowerModelsAnalytics = "~0.4.1"
PowerModelsDistribution = "~0.14.5"
Revise = "~3.5.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "d91d298db48f96d95fb3b99572823ead1906af35"

[[deps.ASL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6252039f98492252f9e47c312c8ffda0e3b9e78d"
uuid = "ae81ac8f-d209-56e5-92de-9978fef736f9"
version = "0.1.3+0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "SnoopPrecompile", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "c700cce799b51c9045473de751e9319bdd1c6e94"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.9"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "485193efd2176b88e6622a39a246f8c5b600e74e"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.6"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "d57c99cc7e637165c81b30eb268eabe156a45c49"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.2.2"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "550e3f8a429fd2f8a9da59f1589c5e268ddc97b3"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.46.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "e32a90da027ca45d84678b826fffd3110bb3fc90"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.8.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "89a9db8d28102b094992472d333674bd1a83ce2a"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "a4ad7ef19d2cdc2eff57abbbe68032b1cd0bd8f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.13.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Glob]]
git-tree-sha1 = "4df9f7e06108728ebf00a0a11edee4b29a482bb2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "37e4657cd56b11abe3d10cd4a1ec5fbdb4180263"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.7.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.InfrastructureModels]]
deps = ["JuMP", "Memento"]
git-tree-sha1 = "88da90ad5d8ca541350c156bea2715f3a23836ce"
uuid = "2030c09a-7f63-5d83-885d-db604e0e9cc0"
version = "0.7.6"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.Ipopt]]
deps = ["Ipopt_jll", "LinearAlgebra", "MathOptInterface", "OpenBLAS32_jll", "SnoopPrecompile"]
git-tree-sha1 = "7690de6bc4eb8d8e3119dc707b5717326c4c0536"
uuid = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
version = "1.2.0"

[[deps.Ipopt_jll]]
deps = ["ASL_jll", "Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "MUMPS_seq_jll", "OpenBLAS32_jll", "Pkg", "libblastrampoline_jll"]
git-tree-sha1 = "563b23f40f1c83f328daa308ce0cdf32b3a72dc4"
uuid = "9cc047cb-c261-5740-88fc-0cf96f7bdcc7"
version = "300.1400.403+1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JSONSchema]]
deps = ["HTTP", "JSON", "URIs"]
git-tree-sha1 = "8d928db71efdc942f10e751564e6bbea1e600dfe"
uuid = "7d188eb4-7ad8-530c-ae41-71a32a6d4692"
version = "1.0.1"

[[deps.JuMP]]
deps = ["LinearAlgebra", "MathOptInterface", "MutableArithmetics", "OrderedCollections", "Printf", "SnoopPrecompile", "SparseArrays"]
git-tree-sha1 = "611b9f12f02c587d860c813743e6cec6264e94d8"
uuid = "4076af6c-e467-56ae-b986-b466b2749572"
version = "1.9.0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "d9ae7a9081d9b1a3b2a5c1d3dac5e2fdaafbd538"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.22"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LightGraphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "432428df5f360964040ed60418dd5601ecd240b6"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.3.5"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "60168780555f3e663c536500aa790b6368adc02a"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.3.0"

[[deps.METIS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1fd0a97409e418b78c53fac671cf4622efdf0f21"
uuid = "d00139f3-1899-568f-a2f0-47f597d42d70"
version = "5.1.2+0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MUMPS_seq_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "METIS_jll", "OpenBLAS32_jll", "Pkg", "libblastrampoline_jll"]
git-tree-sha1 = "f429d6bbe9ad015a2477077c9e89b978b8c26558"
uuid = "d7ed1dd3-d0ae-5e8e-bfb4-87a502085b8d"
version = "500.500.101+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "DataStructures", "ForwardDiff", "JSON", "LinearAlgebra", "MutableArithmetics", "NaNMath", "OrderedCollections", "Printf", "SnoopPrecompile", "SparseArrays", "SpecialFunctions", "Test", "Unicode"]
git-tree-sha1 = "f219b62e601c2f2e8adb7b6c48db8a9caf381c82"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.13.1"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Memento]]
deps = ["Dates", "Distributed", "Requires", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "bb2e8f4d9f400f6e90d57b34860f6abdc51398e5"
uuid = "f28f55f0-a522-5efc-85c2-fe41dfb9b2d9"
version = "1.4.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3295d296288ab1a0a2528feb424b854418acff57"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.2.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NodeJS]]
deps = ["Pkg"]
git-tree-sha1 = "905224bbdd4b555c69bb964514cfa387616f0d3a"
uuid = "2bd173c7-0d6d-553b-b6af-13a54713934c"
version = "1.3.0"

[[deps.OpenBLAS32_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c6c2ed4b7acd2137b878eb96c68e63b76199d0f"
uuid = "656ef2d0-ae68-5445-9ca0-591084a874a2"
version = "0.3.17+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "6503b77492fd7fcb9379bf73cd31035670e3c509"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.PolyhedralRelaxations]]
deps = ["DataStructures", "ForwardDiff", "JuMP", "Logging", "LoggingExtras"]
git-tree-sha1 = "05f2adc696ae9a99be3de99dd8970d00a4dccefe"
uuid = "2e741578-48fa-11ea-2d62-b52c946f73a0"
version = "0.3.5"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PowerModelsAnalytics]]
deps = ["ColorVectorSpace", "Colors", "Compat", "LightGraphs", "LinearAlgebra", "PyCall", "Random", "Setfield", "SpecialFunctions", "Statistics", "Vega"]
git-tree-sha1 = "e222c67fc3d63250e9af159857df8beab55c00bb"
uuid = "7425d11c-6b70-11e9-2bc6-9bcd3191d4b9"
version = "0.4.1"

[[deps.PowerModelsDistribution]]
deps = ["CSV", "Dates", "FilePaths", "Glob", "InfrastructureModels", "JSON", "JuMP", "LinearAlgebra", "Logging", "LoggingExtras", "PolyhedralRelaxations", "Statistics"]
git-tree-sha1 = "497384070b047898b763a7e9c883c5bc196582ad"
uuid = "d7431456-977f-11e9-2de3-97ff7677985e"
version = "0.14.5"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "4ba3651d33ef76e24fef6a598b63ffd1c5e1cd17"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.92.5"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "90cb983381a9dc7d3dff5fb2d1ee52cd59877412"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "77d3c4726515dca71f6d80fbb5e251088defe305"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.18"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "fca29e68c5062722b5b4435594c3d1ba557072a3"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.7.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "6da46b16e6bca4abe1b6c6fa40b94beb0c87f4ac"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.8"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "6aa098ef1012364f2ede6b17bf358c7f1fbe90d4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.17"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.TableTraitsUtils]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Missings", "TableTraits"]
git-tree-sha1 = "78fecfe140d7abb480b53a44f3f85b6aa373c293"
uuid = "382cd787-c1b6-5bf2-a167-d5b971a19bda"
version = "1.0.2"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIParser]]
deps = ["Unicode"]
git-tree-sha1 = "53a9f49546b8d2dd2e688d216421d050c9a31d0d"
uuid = "30578b45-9adc-5946-b283-645ec420af67"
version = "0.4.1"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Vega]]
deps = ["DataStructures", "DataValues", "Dates", "FileIO", "FilePaths", "IteratorInterfaceExtensions", "JSON", "JSONSchema", "MacroTools", "NodeJS", "Pkg", "REPL", "Random", "Setfield", "TableTraits", "TableTraitsUtils", "URIParser"]
git-tree-sha1 = "c6bd0c396ce433dce24c4a64d5a5ab6dc8e40382"
uuid = "239c3e63-733f-47ad-beb7-a12fde22c578"
version = "2.3.1"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─c9a5c344-961b-11eb-0458-a746afcb280c
# ╟─b5928052-9616-11eb-14b3-79770b9929df
# ╟─62c14531-357a-4669-90cd-2a186df123eb
# ╠═9787aa73-8ffc-4634-bf0f-b70eee0bf377
# ╟─b953c65d-515f-4334-a4d4-b27af1b0e29a
# ╠═a1989876-9301-11eb-0783-83b3aa7abfbc
# ╠═f4d7c365-ce18-46b8-8ebb-dd2fbc0e73de
# ╟─749d062d-72a1-4160-afcf-1cdc27d85c84
# ╟─7e032130-c565-42d1-93ed-87955e1f2334
# ╟─67029534-961a-11eb-2b06-210f67607d20
# ╟─0b9598a8-9618-11eb-1947-4f98dac7129f
# ╠═6d94d2cd-069b-4e48-84da-ba87c554297b
# ╟─480f6685-56a0-4fdd-a975-33cdb3459ef5
# ╠═b73eedaf-5809-43b9-90e0-69350d6ee2c4
# ╟─4243344a-9cdf-46db-96ee-633b0f5e8a71
# ╟─e95a9fba-b6ad-4339-8219-de4ca04536f4
# ╠═23a3a107-2fc4-4af1-8eb6-e6a0dc29ddbe
# ╟─ba7c38a6-34f5-40d3-998e-7ea520001454
# ╟─47fdba38-56ae-4a68-88a0-26af4e0a907a
# ╟─4628c72d-32a2-46bb-929e-153f3851f9c0
# ╟─ca82db50-ff48-411c-bc97-44f5f33c4ddf
# ╠═4f45cac5-3b27-4771-994b-3a9d81ed61c9
# ╟─fb08875c-1d9e-4744-bec2-7b7cad4320f6
# ╠═5cf93941-3f88-4bd1-98c4-77b924ddb385
# ╟─8a8222c3-ef01-4d18-97ff-fa43ad89cfcc
# ╠═38dfd21e-f156-45be-b9b3-f6eb5008723e
# ╟─0f73a5be-7ee3-4936-95ed-ebb1b913cf4e
# ╠═254a4ed3-263b-4084-b9af-c97eac2d4ab7
# ╟─03e27a4b-cf53-4e92-9e88-d07015937b3e
# ╟─3e74097f-55ce-4143-938d-543a59581a72
# ╠═2c11f944-ff4c-478b-9b6c-fba54bfc0afd
# ╟─427bebe9-5a0d-41b0-ace6-6622400e136c
# ╟─284a471c-5b97-4753-b5be-8896dc096657
# ╟─20c563f0-9304-11eb-16b4-6720b073e911
# ╠═e70c3f4a-9302-11eb-1176-0925031bccc0
# ╟─9a65d284-8558-4ddf-9139-e88e3a3625d9
# ╟─7edf07c2-208d-47d4-914e-08c09bdfd7f6
# ╟─36f01f08-0c19-4a28-a2eb-11a11dec9897
# ╟─a3e10f7e-c193-4c02-9fef-e23bde112350
# ╟─1cf59273-3f50-41bd-b93c-c5cf5f3cd124
# ╟─7e36b7b8-198c-4b9e-8b62-d1b481f09eeb
# ╠═b717dab6-30fd-4e4b-b6f8-29ffcff80131
# ╠═def53588-f511-4dda-8aa9-941d1aa994f7
# ╟─f1866650-ca90-4f4b-81ef-1cc3202dc240
# ╠═27ec55f0-77c2-4cec-b505-976fd86f1004
# ╠═0ef86c0a-6d65-44a3-9804-676d4cc904c3
# ╟─e45adf3e-5b84-4d14-9840-40ba2fbc1573
# ╟─1b17a770-a465-4ba9-a21d-b3d4514ad696
# ╠═29afef9d-29fe-46fe-9321-a33a488b2ba0
# ╟─2aa3d0ca-115c-476a-be94-f33a2343ea57
# ╠═81930c47-7f02-4eca-a4bc-df6369691eec
# ╟─394b5c41-1fcc-497d-b598-356e0760da62
# ╠═3e2ebbfb-08dd-4099-9367-86c982017e9c
# ╟─f090cadb-c87b-4e36-bc28-79736479933b
# ╠═06643acf-6f8e-4bda-a66b-aab42fe4690b
# ╠═7658ea1b-12e7-475e-898f-d8b31a39665a
# ╟─905bfb0a-2fb5-4b0c-bc62-49bd2888ad30
# ╟─3d468a08-c62c-42cb-a096-8ba382bb91e5
# ╠═114f8b4c-3296-428a-bf72-f5efe6b21b47
# ╟─0d510b81-6f58-4790-be85-82c9b215322d
# ╠═cf3614bc-4c35-4b4d-9ec9-f43307a04c48
# ╟─6335835e-2b5e-45a0-a5f6-4d97328694bd
# ╠═bea3a93e-a827-473a-9e19-348ad85dfa25
# ╟─86d66ed7-89b1-461e-9487-0e190e6fd02e
# ╠═ce39193d-833b-40aa-b532-f9db22f03652
# ╠═cb632c77-5c8b-4fbb-a867-127a22ce54ad
# ╠═b3bd6e16-28dc-44b7-b8ad-1478cdfefecc
# ╟─a0e6612f-516e-4be4-b3ae-00c612b981e2
# ╠═77899c35-2609-4aef-bcd0-178076237d1d
# ╠═eb3d6861-0d91-4ae1-9580-7b3182cefef1
# ╟─1b6fd52f-0cdf-4f4f-aa38-31cd8aa9ca7d
# ╟─e39e89d8-7b3d-4c2d-bbe3-8de202e91c5f
# ╠═493d1dfd-719c-4f03-b660-6b0e9e50a222
# ╟─b6700c39-8dff-4bcb-b6fa-c67b2ddf5163
# ╠═f0cb5faf-0471-4651-9d09-61fe6e5e0fb5
# ╟─89f6282d-9898-4547-a121-028fcdf4f876
# ╠═7303e40a-56f8-4035-8537-e2eed7c16b8b
# ╟─c3dcee8a-edf1-43bb-8733-f2798b29d57a
# ╠═d7522533-43b9-4b3c-a8e5-1aa4f39700bb
# ╟─f4fb1464-8d33-4265-ab40-8fb72d533466
# ╠═213e9ea4-c800-47bc-a020-4eacc334e4fc
# ╟─e32b7e5e-1d99-482f-b406-460a18cde912
# ╠═6b7e5199-61ed-40b0-a658-1a490298dc6d
# ╟─fb48ddde-128a-4a32-ab5f-084561108539
# ╠═92fd9826-17a4-44e2-ae8b-8d96d6f21130
# ╟─a7037bcd-e65d-40d4-8070-7fd73811ffad
# ╠═53e1237b-c380-4506-9924-8b7298d0b38a
# ╟─60485fa1-7b4d-4be1-a0c4-50b45a840a66
# ╠═8e2620f3-894f-43ff-9704-637225a6dda3
# ╟─516bc915-946d-4216-b5d0-b0bde5db6546
# ╠═87bf0ddd-5678-4aea-9024-26c0a2db7521
# ╟─f9fc4c70-f091-4ecd-ac70-c1dfcd8a6f93
# ╠═33dd7091-5576-4375-b5d9-d9a738656ddb
# ╟─1b18e8ce-f403-4cb6-b833-083e8b9ef5bb
# ╠═5369925c-3e0f-43ca-a4d8-c801a90653ff
# ╟─51af76d7-2c62-422d-98d1-b199aba84263
# ╠═12ed532a-164c-4a92-a7d1-fc4dd53c3c57
# ╠═3207f33b-00e5-4691-a9c9-9461c974b3b0
# ╠═8eb99a9b-ee92-4aad-bbfc-c9242dc916cb
# ╠═f70780ef-5ad1-4f69-bb2e-41de737b1bef
# ╟─ab146273-0cf6-4815-b996-dc2e8fc185db
# ╟─3f03084d-f8cc-424e-99a6-78b39ece17ef
# ╠═19e88e22-44e1-4139-bb5b-1247a271e34a
# ╟─7af68b51-a1b3-4ed8-95bc-43b2523c298c
# ╠═3d769b34-b327-4304-af0a-1b995a4d6e4c
# ╟─9bb9cbf5-67ef-487b-88a9-830f15378537
# ╠═4bbc5205-4a3e-4e3d-88d6-219195ad8a0f
# ╟─ba8e830f-3c3b-45bc-be25-c440cc736a02
# ╠═ea462d05-6325-44b3-84fc-bf34e032bb01
# ╟─565385e3-ccb4-4ac7-a0a5-12c2744ceed5
# ╠═fb05c680-d4ba-46a6-a696-7536aaa0dd57
# ╟─5cd2d702-f5fa-44f3-9739-8eb75b0b6472
# ╠═85d52ceb-2b4f-4c22-880a-290c1a5737a1
# ╟─20f6dfed-5d1a-481c-92f4-440de0f2abca
# ╠═ff6dfa9a-d55e-4bbe-bf41-357a417784a5
# ╟─74a403df-0f94-4a1a-af69-cd8dc5bef430
# ╠═28f75060-d0be-47f2-bea9-377d5df89cf6
# ╟─663e32f5-4b2e-4255-9def-e1c7d383c36d
# ╠═d3e6486c-4133-404b-91c1-7ac18ec8388f
# ╟─8c34cd7f-52c4-440f-a226-afc27bf24799
# ╠═3562dcc2-b6d7-4bfa-bae7-7a5e992cb422
# ╠═805b06be-7ed9-4330-9f1e-a584c4508862
# ╟─e56a2563-45c7-40fd-80ad-b86d7e1c1e44
# ╟─49808060-ad86-491e-aed2-c69f203d4a0e
# ╠═59abd602-30de-40af-bce5-7800d11a7108
# ╟─f40614c0-5aac-491b-bf5b-8aa9d328f640
# ╠═40e60c21-76fe-44ce-9a33-988030e8d863
# ╠═14f71709-35f9-4eb7-8505-a2138b6fb6f5
# ╠═09f7c681-9654-466d-a1cb-e09eebad62fd
# ╟─0362394c-9ac9-4978-9b41-bc325464c204
# ╠═8361deb2-3baa-4031-927f-e3c034bd7eb9
# ╠═ba259614-bead-40c4-927f-5f8ad3eb726a
# ╟─dee28964-882f-4bad-9368-762f7af99562
# ╠═2f10d9b4-c237-4305-ac05-ffe569d70436
# ╠═26e7b26e-3619-40d0-aef8-0a070af98cbb
# ╠═0c39a147-8587-457e-9a2b-ea38e2bf160a
# ╟─ec94eeb4-e425-452d-8e46-6ed7f876f463
# ╟─292e7528-cf57-4429-995a-571d567b8560
# ╠═2de210f3-75b4-4c5d-89a7-65557b45df2c
# ╠═687863f6-1aed-4039-b54e-d9d3a683631b
# ╟─b448fae8-56d4-40a0-95fa-d502f8a37a40
# ╠═53f3cd12-05a9-4495-937a-781ce3907174
# ╟─232a4bde-3f85-46c4-9e57-ffc64ecafbb6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
