### A Pluto.jl notebook ###
# v0.14.2

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

# ╔═╡ 9e6f21ee-5caf-4bcc-aca5-4d9b55a4acc9
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.Registry.update()
	Pkg.add([
			Pkg.PackageSpec(;name="Revise"),
			Pkg.PackageSpec(;name="CodeTracking"),
			Pkg.PackageSpec(;name="PlutoUI"),
			Pkg.PackageSpec(;name="PowerModelsDistribution", version="0.11"),
			Pkg.PackageSpec(;name="PowerModelsAnalytics", version="0.4.1"),
			Pkg.PackageSpec(;name="InfrastructureModels", version="0.6"),
			Pkg.PackageSpec(;name="JuMP", version="0.21.7"),
			Pkg.PackageSpec(;name="Ipopt", version="0.6.5"),
			Pkg.PackageSpec(;name="JSON", version="0.21"),
			])
end

# ╔═╡ 9787aa73-8ffc-4634-bf0f-b70eee0bf377
using CodeTracking, Revise, PlutoUI

# ╔═╡ a1989876-9301-11eb-0783-83b3aa7abfbc
begin
	using PowerModelsDistribution
	using PowerModelsAnalytics
	import InfrastructureModels as IM
	import JuMP
	import Ipopt
	import JSON
end

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
- `PowerModelsDistribution = "~0.11"`
- `PowerModelsAnalytics = "~0.4.1"`

This notebook is a begginer's introduction to PowerModelsDistribution, an optimization-focused Julia library for quasi-steady state power distribution modeling, based on JuMP.jl, and part of the larger [InfrastructureModels.jl](https://github.com/lanl-ansi/InfrastructureModels.jl) ecosystem, which notably includes:

- [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) : Transmission (single-phase positive sequence power networks) optimization
- [GasModels.jl](https://github.com/lanl-ansi/GasModels.jl) : Natural Gas pipeline optimization (includes Steady-state and Transient optimization)
- [WaterModels.jl](https://github.com/lanl-ansi/WaterModels.jl) : Water network steady-state optimization

Details about PowerModelsDistribution.jl can be found in our [PSCC Conference Proceedings paper](https://doi.org/10.1016/j.epsr.2020.106664).

## Julia Environment Setup

The following code block will setup a Julia environment for you with the correct versions of packages for this Pluto notebook
"""

# ╔═╡ 62c14531-357a-4669-90cd-2a186df123eb
md"The following packages are used for notebook features only and do not relate to tutorial content"

# ╔═╡ b953c65d-515f-4334-a4d4-b27af1b0e29a
md"""
This notebook will make use of the following packages in various places
"""

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

If you are familiar with JuMP, you probably are already used to Enums `TerminationStatusCode` and `ResultStatusCode`, which we import explicitly from MathOptInterface and export, for easy access by the user when using `using PowerModelsDistribution`.

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
md"and the following enum values exist currently in PMD (excluding those imported from MathOptInterface):"

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

# ╔═╡ Cell order:
# ╟─c9a5c344-961b-11eb-0458-a746afcb280c
# ╟─b5928052-9616-11eb-14b3-79770b9929df
# ╠═9e6f21ee-5caf-4bcc-aca5-4d9b55a4acc9
# ╟─62c14531-357a-4669-90cd-2a186df123eb
# ╠═9787aa73-8ffc-4634-bf0f-b70eee0bf377
# ╟─b953c65d-515f-4334-a4d4-b27af1b0e29a
# ╠═a1989876-9301-11eb-0783-83b3aa7abfbc
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
# ╠═687863f6-1aed-4039-b54e-d9d3a683631b
# ╟─b448fae8-56d4-40a0-95fa-d502f8a37a40
# ╠═53f3cd12-05a9-4495-937a-781ce3907174
# ╟─232a4bde-3f85-46c4-9e57-ffc64ecafbb6
