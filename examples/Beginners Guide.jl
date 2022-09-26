### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
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
begin
	eng_copy = deepcopy(eng)
	update_data!(eng_copy, eng_result["solution"])
end

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
CodeTracking = "~1.0.6"
InfrastructureModels = "~0.6.1"
Ipopt = "~0.7.0"
JSON = "~0.21.2"
JuMP = "~0.21.10"
PlutoUI = "~0.7.16"
PowerModelsAnalytics = "~0.4.1"
PowerModelsDistribution = "~0.12.0"
Revise = "~3.1.20"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ASL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "370cafc70604b2522f2c7cf9915ebcd17b4cd38b"
uuid = "ae81ac8f-d209-56e5-92de-9978fef736f9"
version = "0.1.2+0"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "61adeb0823084487000600ef8b1c00cc2474cd47"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.2.0"

[[BinaryProvider]]
deps = ["Libdl", "Logging", "SHA"]
git-tree-sha1 = "ecdec412a9abc8db54c0efc5548c64dfce072058"
uuid = "b99e7846-7c00-51b0-8f62-c81ae34c0232"
version = "0.5.10"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "c0a735698d1a0a388c5c7ae9c7fb3da72fd5424e"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.9.9"

[[Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "3533f5a691e60601fe60c90d8bc47a27aa2907ec"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.0"

[[CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "9aa8a5ebb6b5bf469a7e0e2b5202cf6f8c291104"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.0.6"

[[CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "45efb332df2e86f2cb2e992239b6267d97c9e0b6"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.7"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "31d0151f5716b655421d9d75b7fa74cc4e744df2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.39.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Conda]]
deps = ["JSON", "VersionParsing"]
git-tree-sha1 = "299304989a5e6473d985212c28928899c74e9421"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.5.2"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "7220bc21c33e990c14f4a9a319b1d242ebc5b269"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.3.1"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "3c041d2ac0a52a12a27af2782b34900d9c3ee68c"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.11.1"

[[FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[FilePathsBase]]
deps = ["Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "7fb0eaac190a7a68a56d2407a6beff1142daf844"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.12"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "63777916efbcb0ab6173d09a658fb7f2783de485"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.21"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "14eece7a3308b4d8be910e265c724a6ba51a9798"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.16"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "5efcf53d798efede8fee5b2c8b09284be359bf24"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.2"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[InfrastructureModels]]
deps = ["JuMP", "MathOptInterface", "Memento"]
git-tree-sha1 = "14ad999045d7ea263cb909b1942814850227d457"
uuid = "2030c09a-7f63-5d83-885d-db604e0e9cc0"
version = "0.6.1"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "19cb49649f8c41de7fea32d089d37de917b553da"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.0.1"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "f0c6489b12d28fb4c2103073ec7452f3423bd308"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.1"

[[Ipopt]]
deps = ["BinaryProvider", "Ipopt_jll", "Libdl", "LinearAlgebra", "MathOptInterface", "MathProgBase"]
git-tree-sha1 = "380786b4929b8d18d76e909c6b2eca355b7c3bd6"
uuid = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
version = "0.7.0"

[[Ipopt_jll]]
deps = ["ASL_jll", "Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "MUMPS_seq_jll", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "82124f27743f2802c23fcb05febc517d0b15d86e"
uuid = "9cc047cb-c261-5740-88fc-0cf96f7bdcc7"
version = "3.13.4+2"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JSONSchema]]
deps = ["HTTP", "JSON", "URIs"]
git-tree-sha1 = "2f49f7f86762a0fbbeef84912265a1ae61c4ef80"
uuid = "7d188eb4-7ad8-530c-ae41-71a32a6d4692"
version = "0.3.4"

[[JuMP]]
deps = ["Calculus", "DataStructures", "ForwardDiff", "JSON", "LinearAlgebra", "MathOptInterface", "MutableArithmetics", "NaNMath", "Printf", "Random", "SparseArrays", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "4358b7cbf2db36596bdbbe3becc6b9d87e4eb8f5"
uuid = "4076af6c-e467-56ae-b986-b466b2749572"
version = "0.21.10"

[[JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "e273807f38074f033d94207a201e6e827d8417db"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.8.21"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LightGraphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "432428df5f360964040ed60418dd5601ecd240b6"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.3.5"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "6193c3815f13ba1b78a51ce391db8be016ae9214"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.4"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "dfeda1c1130990428720de0024d4516b1902ce98"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.7"

[[LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "491a883c4fef1103077a7f648961adbf9c8dd933"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.1.2"

[[METIS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "2dc1a9fc87e57e32b1fc186db78811157b30c118"
uuid = "d00139f3-1899-568f-a2f0-47f597d42d70"
version = "5.1.0+5"

[[MUMPS_seq_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "METIS_jll", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "1a11a84b2af5feb5a62a820574804056cdc59c39"
uuid = "d7ed1dd3-d0ae-5e8e-bfb4-87a502085b8d"
version = "5.2.1+4"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "5a5bc6bf062f0f95e62d0fe0a2d99699fed82dd9"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.8"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "JSON", "JSONSchema", "LinearAlgebra", "MutableArithmetics", "OrderedCollections", "SparseArrays", "Test", "Unicode"]
git-tree-sha1 = "575644e3c05b258250bb599e57cf73bbf1062901"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "0.9.22"

[[MathProgBase]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9abbe463a1e9fc507f12a69e7f29346c2cdc472c"
uuid = "fdba3010-5040-5b88-9595-932c9decdf73"
version = "0.7.8"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Memento]]
deps = ["Dates", "Distributed", "JSON", "Serialization", "Sockets", "Syslogs", "Test", "TimeZones", "UUIDs"]
git-tree-sha1 = "19650888f97362a2ae6c84f0f5f6cda84c30ac38"
uuid = "f28f55f0-a522-5efc-85c2-fe41dfb9b2d9"
version = "1.2.0"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "29714d0a7a8083bba8427a4fbfb00a540c681ce7"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.3"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "8d9496b2339095901106961f44718920732616bb"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.22"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[NodeJS]]
deps = ["Pkg"]
git-tree-sha1 = "905224bbdd4b555c69bb964514cfa387616f0d3a"
uuid = "2bd173c7-0d6d-553b-b6af-13a54713934c"
version = "1.3.0"

[[OpenBLAS32_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ba4a8f683303c9082e84afba96f25af3c7fb2436"
uuid = "656ef2d0-ae68-5445-9ca0-591084a874a2"
version = "0.3.12+1"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "f19e978f81eca5fd7620650d7dbea58f825802ee"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlutoUI]]
deps = ["Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "4c8a7d080daca18545c56f1cac28710c362478f3"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.16"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a193d6ad9c45ada72c14b731a318bedd3c2f00cf"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.3.0"

[[PowerModelsAnalytics]]
deps = ["ColorVectorSpace", "Colors", "Compat", "LightGraphs", "LinearAlgebra", "PyCall", "Random", "Setfield", "SpecialFunctions", "Statistics", "Vega"]
git-tree-sha1 = "e222c67fc3d63250e9af159857df8beab55c00bb"
uuid = "7425d11c-6b70-11e9-2bc6-9bcd3191d4b9"
version = "0.4.1"

[[PowerModelsDistribution]]
deps = ["CSV", "Dates", "InfrastructureModels", "JSON", "JuMP", "LinearAlgebra", "Logging", "LoggingExtras", "MathOptInterface", "Statistics"]
git-tree-sha1 = "7103d919119554bae66cf3677f8fad3373446a78"
uuid = "d7431456-977f-11e9-2de3-97ff7677985e"
version = "0.12.0"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "4ba3651d33ef76e24fef6a598b63ffd1c5e1cd17"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.92.5"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "41deb3df28ecf75307b6e492a738821b031f8425"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.1.20"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "f45b34656397a1f6e729901dc9ef679610bd12b5"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.8"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "fca29e68c5062722b5b4435594c3d1ba557072a3"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.7.1"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f0bccf98e16759818ffc5d97ac3ebf87eb950150"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.1"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[Syslogs]]
deps = ["Printf", "Sockets"]
git-tree-sha1 = "46badfcc7c6e74535cc7d833a91f4ac4f805f86d"
uuid = "cea106d9-e007-5e6c-ad93-58fe2094e9c4"
version = "0.3.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[TableTraitsUtils]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Missings", "TableTraits"]
git-tree-sha1 = "78fecfe140d7abb480b53a44f3f85b6aa373c293"
uuid = "382cd787-c1b6-5bf2-a167-d5b971a19bda"
version = "1.0.2"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "b4c6460412b1db0b4f1679ab2d5ef72568a14a57"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.6.1"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[URIParser]]
deps = ["Unicode"]
git-tree-sha1 = "53a9f49546b8d2dd2e688d216421d050c9a31d0d"
uuid = "30578b45-9adc-5946-b283-645ec420af67"
version = "0.4.1"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Vega]]
deps = ["DataStructures", "DataValues", "Dates", "FileIO", "FilePaths", "IteratorInterfaceExtensions", "JSON", "JSONSchema", "MacroTools", "NodeJS", "Pkg", "REPL", "Random", "Setfield", "TableTraits", "TableTraitsUtils", "URIParser"]
git-tree-sha1 = "43f83d3119a868874d18da6bca0f4b5b6aae53f7"
uuid = "239c3e63-733f-47ad-beb7-a12fde22c578"
version = "2.3.0"

[[VersionParsing]]
git-tree-sha1 = "e575cf85535c7c3292b4d89d89cc29e8c3098e47"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.1"

[[WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "c69f9da3ff2f4f02e811c3323c22e5dfcb584cfa"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.1"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─c9a5c344-961b-11eb-0458-a746afcb280c
# ╟─b5928052-9616-11eb-14b3-79770b9929df
# ╟─62c14531-357a-4669-90cd-2a186df123eb
# ╠═9787aa73-8ffc-4634-bf0f-b70eee0bf377
# ╟─b953c65d-515f-4334-a4d4-b27af1b0e29a
# ╠═a1989876-9301-11eb-0783-83b3aa7abfbc
# ╟─f4d7c365-ce18-46b8-8ebb-dd2fbc0e73de
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
