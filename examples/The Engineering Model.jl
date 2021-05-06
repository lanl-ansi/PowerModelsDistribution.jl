### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ 2661286b-2d18-42e3-b309-7974fc2db425
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

# ╔═╡ f30cd0d0-b0da-4f63-a245-568a763a93d8
begin
	using PowerModelsDistribution
	using Ipopt
end

# ╔═╡ c55b2c42-9d27-11eb-24ca-e90a5472ffbb
md"""
# Introduction to the PowerModelsDistribution Data Models

In this notebook we introduce the engineering data model added to PowerModelsDistribution in version v0.9.0. We will give several examples of how to use this new data model directly, including new transformations that have become easier with its introduction, how to convert it to the the lower-level mathematical model that was previously the only user interface we offered, and how to get various types of results using this new model.

## Julia Environment Setup

The following code block will setup a Julia environment for you with the correct versions of packages for this Pluto notebook.
"""

# ╔═╡ de92bc20-b125-4f3d-930a-b1da63d5cef5
md"""
## Imports

All commands in this document with no package namespace specified are directly exported by PowerModelsDistribution or already available in Julia base. Any commands that are only avaiable via an external package will be specified by including by using `import`, which will require specifying the originating package before the command, _e.g._ `Ipopt.Optimizer` as you will see below.
"""

# ╔═╡ ac94e556-b544-4ca6-87bd-ff2d7a7414e7
pmd_path = joinpath(dirname(pathof(PowerModelsDistribution)), "..")

# ╔═╡ 1e791262-261d-4756-bec4-edebe4732700
md"In these examples we will use the following optimization solvers, specified using `optimizer_with_attributes` from JuMP v0.21"

# ╔═╡ 89de80df-1dd0-4f94-a27d-74693a978059
ipopt_solver = optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>0)

# ╔═╡ e22a7d3e-3f21-41b0-abf2-2c723e87c57e
md"""
## Parsing Data

Here we give the first example of how to parse data into the `ENGINEERING` data model structure, which is the default data structure type that the user will see without passing additional arguments, as we demonstrate later.

We start with a 3 bus unbalanced load case provided as a dss file in the `test` folder of the PowerModelsDistribution.jl repository
"""

# ╔═╡ 5fe2d186-40c6-46c2-823e-9401cb3b6d6c
eng = parse_file("$pmd_path/test/data/opendss/case3_unbalanced.dss")

# ╔═╡ c884fd5b-4ecb-402b-b228-01b9a61db8bf
md"""
Different information and warning messages will be given depending on the input file. In the case above, these messages all related to various parse notifications that arise during a parse of a dss file, and can be safely ignored

The resulting data structure is a Julia dictionary. The first notable field is `"data_model"` which specifies which data model this data structure corresponds to, in this case `ENGINEERING`. This value is expected to be an `Enum` of type `DataModel`

The next notable field is `"settings"`, which contains some important default/starting values for the distribution network
"""

# ╔═╡ a5b736f9-2776-4760-b073-d02027baef13
eng["settings"]

# ╔═╡ 29a8a560-7d5e-4929-877a-2dab17309968
md"""
- `"sbase_default"` is the starting value for the power base,
- `"vbases_default"` is the starting voltage base for the case, and multiple voltage bases can be specified, which would be useful in cases where there are multiple isolated islands with their own generation,
- `"voltage_scale_factor"` is a scaling factor for all voltage values, which in the case of OpenDSS is in kV by default
- `"power_scale_factor"` is a scaling factor for all power values
- `"base_frequency"` is the base frequency of the network in Hz, which is useful to know for mixed frequency networks

Next we look at the `"bus"` components
"""

# ╔═╡ 9ac6b0ed-58a8-4079-aa9f-d70024f4d4b4
eng["bus"]

# ╔═╡ 89846878-ebfc-49a5-91ef-762f692ba5ea
md"""
We can see there are three buses in this system, identified by ids `"primary"`, `"sourcebus"`, and `"loadbus"`. 

__NOTE__: In Julia, order of Dictionary keys is not fixed, nor does it retain the order in which it was parsed like _e.g._ `Vectors`. 

Identifying components by non-integer names is a new feature of the `ENGINEERING` model, and makes network debugging more straightforward. 

__NOTE__: all names are converted to lowercase on parse from the originating dss file.

Each bus component has the following properties in the `ENGINEERING` model
"""

# ╔═╡ 8ba081b4-8444-48e0-afce-20649f7fdd01
eng["bus"]["sourcebus"]

# ╔═╡ 2a8c1536-ef5f-4e37-b30e-1bb73982d7b0
md"""
- `"terminals"` indicates which terminals on the bus have active connections
- `"grounded"` indicates which terminals are grounded
- `"rg"` and `"xg"` indicate the grounding resistance and reactance of the ground
- `"status"` indicates whether a bus is `ENABLED` or `DISABLED`, and is specified for every component in the engineering model

Next, we look at the `"line"` components, which is a generic name for both overhead lines and underground cables, which we do not differentiate between in the nomenclature
"""

# ╔═╡ 83ce434c-0f92-46a8-8648-322084044600
eng["line"]

# ╔═╡ 2f9c3f64-3e5d-423f-939e-8d96bd19e7fe
eng["line"]["quad"]

# ╔═╡ f06b9592-0e16-4fd1-bf6e-624ec0cdb8fb
md"""
Again, we see components identified by their OpenDSS names. A `"line"` is an edge object, which will always have the following properties:

- `"f_bus"`
- `"t_bus"`
- `"f_connections"` - list of terminals to which the line is connected on the from-side
- `"t_connections"` - list of terminals to which the line is connected on the to-side

Here we are also introduced to two important concepts, the `"source_id"`, which is an easy way to identify from where an object originates in the dss file, and a data type element, pointed to by `"linecode"` in this case.

A data type element is an element that does not represent a real engineering object, but only contains data that one of those real objects can refer to, in this case a linecode, which contains information like line resistance/reactance and conductance/susceptance.
"""

# ╔═╡ b15593b9-f775-47af-aa9b-4980a4028faa
eng["linecode"]["4/0quad"]

# ╔═╡ 2e4fe000-449c-4fc0-8707-4c51ac50ab03
md"""
We can see that the length of the Vectors for `"pd_nom"` and `"qd_nom"` are only one, although the number of terminals listed in `"connections"` is two. This is because the connection is WYE, and therefore the final connection is a grounded neutral

Here we are also introduced to two new Enums, `WYE`, which gives the connection configuration, and `NO` under dispatchable, which indicates that if this case were used in an MLD problem, _i.e._ with `run_mc_mld` that this load would not be sheddable.

Finally, we show the generation source for this case, which in opendss is a voltage source named `"source"`
"""

# ╔═╡ 36ca17bc-b560-48c9-988f-255780479d05
eng["voltage_source"]["source"]

# ╔═╡ 7e5811cf-86a2-43fb-8c94-44091f33c031
md"""
- `"vm"` - specifies the fixed voltage magnitudes per phase at the bus
- `"va"` - specifies the fixed reference angles per phases at the bus
- `"rs"` and `"xs"` specifies internal impedances of the voltage source

### Importing raw dss properties

In case there are additional properties that you want to use from dss, it is possible to import those directly into the `ENGINEERING` (and `MATHEMATICAL`) data structure with the `import_all` keyword argument
"""

# ╔═╡ 026088f7-b1df-4647-ad81-52c6c8ea7944
eng_all = parse_file("../test/data/opendss/case3_unbalanced.dss"; import_all=true)

# ╔═╡ 5f72ecec-bb3f-4231-93c1-719e38e17be4
md"""
You will note the presence of `"dss"` dictionaries under components, and `"dss_options"` at the root level
"""

# ╔═╡ f070dbbc-1a55-48f2-aa21-2119eb573b5b
eng_all["line"]

# ╔═╡ f3458a5c-70c6-4d0c-8253-4feb0c87ee76
md"""
### Time Series Parsing Example

In the `ENGINEERING` model, we have included the `time_series` data type, which holds all time series data and can be referred to similar to `"linecode"` as demonstrated above.

Below we can see an example of a parse that includes some time_series components
"""

# ╔═╡ 542e592f-8f75-41f9-8f98-5373b431fcc6
eng_ts = parse_file("../test/data/opendss/case3_balanced.dss"; time_series="daily")

# ╔═╡ d0e02ee5-f57c-4286-940d-b01097c840be
md"""
You can see that under the actual component, in this case a `"load"`, that there is a `"time_series"` dictionary that contains `ENGINEERING` model variable names and references to the identifiers of a root-level `time_series` object, 
"""

# ╔═╡ e852712c-d09a-41e1-9785-308c219c7ad8
eng_ts["time_series"]["ls1"]

# ╔═╡ 0eba403e-c2ba-4af3-8dcc-fe14a3c22bd7
md"""
This feature is useful for building multinetwork data structures, which will be described below in the section on the `MATHEMATICAL` model
"""

# ╔═╡ 8d5701d1-9e91-43b5-975b-fc8e0e306e99
md"""
## Running Optimal Power Flow

In this section we introduce how to run an optimal power flow (opf) in PowerModelsDistribution on an engineering data model

In order to run an OPF problem you will need

1. a data model
2. a formulation
3. a solver

In these examples we will use the `eng` model we worked with above, the `ACPUPowerModel`, which is a AC power flow formulation in polar coordinates, and the `ipopt_solver` we already defined above
"""

# ╔═╡ 745f33bf-98d7-4b27-95fd-ec6bf5e3d5cf
result = solve_mc_opf(eng, ACPUPowerModel, ipopt_solver)

# ╔═╡ fb81d3c1-a47b-46d9-b463-1766022be114
md"""
The result of `solve_mc_opf` will be very familiar to those who are already familiar with PowerModels and PowerModelsDistribution. The notable difference will be in the `"solution"` dictionary
"""

# ╔═╡ 27d5a237-4c2e-4b14-a082-f4847153c8d5
result["solution"]

# ╔═╡ 564200e0-2080-44da-8dbf-37c5e3a38412
md"""
Here you can see that the solution comes back out by default into the same data model as is provided by the user to the run_ command, as well as being in SI units, as opposed to per unit, which is used during the solve. For example,
"""

# ╔═╡ 8b066039-ed65-4617-b404-86952d3fc78e
result["solution"]["bus"]["loadbus"]

# ╔═╡ a09de8ab-c1e1-4c9f-9d86-4268ce503cae
md"""
If for some reason you want to return the result in per-unit rather than SI, you can specify this in the `solve_` command by
"""

# ╔═╡ 763a1504-f927-4dfc-a042-ac66324708a9
result_pu = solve_mc_opf(eng, ACPUPowerModel, ipopt_solver; make_si=false)

# ╔═╡ e7b588e9-fb5f-4383-b5c6-bce9ed8f5a60
result_pu["solution"]["bus"]["loadbus"]

# ╔═╡ f882f070-74af-4f3b-9ad1-bd7ad453778f
md"""
### Branch Flow formulations

Previously, to use a branch flow formulation, such as `SOCNLPUBFPowerModel`, it was required to use a different `solve_` command, but now, by using multiple dispatch we have simplified this for the user
"""

# ╔═╡ f2fa779e-f62d-4053-b055-e9f0aee6c3f7
result_bf = solve_mc_opf(eng, SOCNLPUBFPowerModel, ipopt_solver)

# ╔═╡ 544108d1-d1a7-4dc6-9a3c-6a10025cb8a3
md"""
### Running Time Series Models

By default, `time_series` object will be ignored when running a model. To use the time series information you will need to have a multinetwork problem specification

In the example below we use a test case, which is not exported by default, and therefore requires the specification of the PowerModelsDistribution namespace
"""

# ╔═╡ ae8daf6c-4055-497e-905c-d6d5a6373ad3
result_mn = PowerModelsDistribution._solve_mn_mc_opb(eng_ts, NFAUPowerModel, ipopt_solver)

# ╔═╡ 3d86989a-a4db-4766-aec7-abc9a25ffde8
md"""
## Engineering Model Transformations

One of the power things about the engineering model is that data transformations are much more simple. Here we illustrate two examples that are currently included in PowerModelsDistribution, but writing your own data transformation functions will be trivial, as we will show.

__Note__: In v0.9, `apply_kron_reduction!` and `apply_phase_projection!` are applied by default, but can be disabled with the keyword arguments `kron_reduced=false` and `project_phases=false`, respectively in `parse_file` or `transform_data_model`.

First, there are several objects that have loss models by default when parsing from dss files, such as voltage sources, transformers, and switches. To remove these loss models, therefore making these components lossless, we can use the included `make_lossess!` function. Here we use a basic 2-winding wye-wye connected transformer case from `test` to illustrate this
"""

# ╔═╡ 7a86e26d-fdf6-41b6-b1dd-dc0c4a8e9d50
eng_ut = parse_file("../test/data/opendss/ut_trans_2w_yy.dss")

# ╔═╡ 82e9330b-2e94-42a1-b56d-430674b11241
eng_ut["transformer"]["tx1"]

# ╔═╡ 2c358d7a-5d84-4f62-a147-577bb775ac7f
md"""
We can see that "noloadloss", "rw", and "imag" are non-zero, but if we apply the make_lossless! function we can see these parameters are set to zero, effectively eliminating the losses
"""

# ╔═╡ 5e5dd24c-95b5-4951-8547-13e20bb0feec
make_lossless!(eng_ut)

# ╔═╡ 45c320dd-e0b2-46ec-af09-d83e5652026d
eng_ut["transformer"]["tx1"]

# ╔═╡ 723f4e26-8065-4e22-8c5c-de32d2e42b47
md"Alternatively, we can apply this function at parse"

# ╔═╡ 58f80651-4a44-4781-a6c0-a056ec76b07e
parse_file("../test/data/opendss/ut_trans_2w_yy.dss"; transformations=[make_lossless!])

# ╔═╡ b1e55263-9c76-43b2-b91c-115b2f080182
md"""
Another transformation function included in PowerModelsDistribution is the `apply_voltage_bounds!` function, which will apply some voltage bounds in SI units, given some percent value, _e.g._ if we want the lower bound on voltage to be `0.9` and upper bound `1.1` after per-unit conversion
"""

# ╔═╡ 2cdd682a-b97f-49dc-8829-47e7146497b9
apply_voltage_bounds!(eng_ut; vm_lb=0.9, vm_ub=1.1)

# ╔═╡ 74836b2c-db67-4435-a552-ba5c3e93e43f
eng_ut["bus"]["2"]

# ╔═╡ 2bf802e2-f769-48db-adef-a5a066c565c2
md"Alternatively, this can be specified at parse by"

# ╔═╡ 4020d67a-4a9e-4030-b7cd-e5d60e78495c
parse_file("../test/data/opendss/ut_trans_2w_yy.dss"; transformations=[make_lossless!, (apply_voltage_bounds!, "vm_lb"=>0.9, "vm_ub"=>1.1)])

# ╔═╡ dde2da59-63fa-4549-8829-491e5f3c675e
md"""
### Transformations on Multinetworks

Transformations on Multinetworks should happen __before__ the network is converted into a `MATHEMATICAL` data model, so that they can generally follow the same pattern as shown above and can be seen in the `make_lossless!` and `apply_voltage_bounds!` functions already in PowerModelsDistribution
"""

# ╔═╡ 2c375878-f16b-4833-9aaa-b826fbe056a2
md"""
## Mathematical Model

In this section we introduce the mathematical model, which was the previous user-facing model in PowerModelsDistribution, explain how conversions between the model happen in practice, and give an example of how to do this conversion manually

In practice, unless the user is interested, the conversion between the `ENGINEERING` and `MATHEMATICAL` models should be seemless and invisible to the user. By providing an `ENGINEERING` model to a `run_` command the `run_mc_model` command will know to convert the model to `MATHEMATICAL`, which will be used to the generate the JuMP model that will actually be optimized. Similarly, the solution generated by this optimization will be automatically converted back to the format of the `ENGINEERING` model.

Let's first take a look at how to convert to the `MATHEMATICAL` model
"""

# ╔═╡ a7348fc8-b25e-4ecf-9fe6-314036b7df9c
math = transform_data_model(eng)

# ╔═╡ 54280696-beec-47ee-8114-ea9311aca9eb
md"""
There are a couple of things to notice right away. First, the data model transform automatically converts the model to per-unit. Second, there are a lot of empty component sets, whereas in the `ENGINEERING` model, only component types that had components in them were listed. In the `MATHEMATICAL` model certain component dictionaries are always expected to exist, and the `eng2math` conversion functions will automatically populate these.

Next, there are a few unusal fields, such as `"settings"`, which previously didn't exist in the `MATHEMATICAL` model. This is used for the per-unit conversion specifically in PowerModelsDistribution. Also, is the `"map"` field, which is a `Vector` of Dictionaries that enable the conversion back to `ENGINEERING` from `MATHEMATICAL`. Without this it would be impossible to convert back, and in fact only the solution can be converted, because some properties are combined destructively during the conversion to the `MATHEMATICAL` model, and therefore cannot be reverse engineered. However, since the conversion to `MATHEMATICAL` is not in-place, you will always have a copy of `eng` alongside `math`.

Here is an example of one of the `"map"` entries
"""

# ╔═╡ 1dc06007-b413-4cbb-b43b-2c7955ecb6d2
math["map"][end]

# ╔═╡ 96a65bdc-7ee4-4e9f-94f5-8a83f7512fe9
md"""
Alternatively, the `MATHEMATICAL` model can be returned directly from the `parse_file` command with the `data_model` keyword argument 
"""

# ╔═╡ 545c9270-257a-46e8-970d-6a9fcf09dc17
parse_file("../test/data/opendss/case3_unbalanced.dss"; data_model=MATHEMATICAL)

# ╔═╡ fda420a1-e563-4a41-84d4-2a570bf163ae
md"""
### Multinetworks

In this subsection we cover parsing into a multinetwork data structure, which is a structure that only exists in the `MATHEMATICAL` model

For those unfamiliar, the InfrastructureModels family of packages has a feature called multinetworks, which is useful for, among other things, running optimization problems on time series type problems. 

Multinetwork data structures are formatted like so

```julia
mn = Dict{String,Any}(
    "multinetwork" => true,
    "nw" => Dict{String,Any}(
        "1" => Dict{String,Any}(
            "bus" => Dict{String,Any}(),
            ...
        ),
        ...
    ),
    ...
)
```

To automatically create a multinetwork structure from an engineering model that contains `time_series` elements, we can use the `multinetwork` keyword argument in `transform_data_model`
"""

# ╔═╡ 785a1c39-aa56-4244-b7b9-31d7a31b73c3
math_mn = transform_data_model(eng_ts; multinetwork=true)

# ╔═╡ ad48b3ad-ede1-4109-9647-007c1d8d196c
parse_file("../test/data/opendss/case3_balanced.dss"; multinetwork=true, data_model=MATHEMATICAL)

# ╔═╡ 666223d5-802a-4164-882c-8ab1a1319d6d
md"""
### Running `MATHEMATICAL` models

There is very little difference from the user point-of-view in running `MATHEMATICAL` models other than the results will not be automatically converted back to the the format of the `ENGINEERING` model
"""

# ╔═╡ e3fe0a5e-9431-453d-b7e5-821de624a4ef
result_math = solve_mc_opf(math, ACPUPowerModel, ipopt_solver)

# ╔═╡ 5ab6cd16-3c24-4f60-bdf2-7f48b3e25085
md"""
It is also possible to manually convert the solution back to the `ENGINEERING` format, provided you have the __map__
"""

# ╔═╡ 05360163-4a6f-41da-8fae-070e327c2605
sol_eng = transform_solution(result_math["solution"], math)

# ╔═╡ b8e01ce2-d7f7-4c73-8230-9fdeecab2c4e
md"""
#### Running `MATHEMATICAL` Multinetworks

As with the `ENGINEERING` example of running a multinetwork problem, you will need a multinetwork problem specification, and as with the previous single `MATHEMATICAL` network example above, we only obtain the `MATHEMATICAL` solution, and can transform the solution in the same manner as before
"""

# ╔═╡ 96757dde-52be-4930-9bdd-9303bb715e55
result_math_mn = PowerModelsDistribution._solve_mn_mc_opb(math_mn, NFAUPowerModel, ipopt_solver)

# ╔═╡ 867d8d1d-042a-42ec-843e-56b8ba2db96a
sol_eng_mn = transform_solution(result_math_mn["solution"], math_mn)

# ╔═╡ db4d9eb9-e120-4489-a3c2-48c43f2a3c3b
md"""
## Building the JuMP Model

In some cases the user will want to directly build the JuMP model, which would traditionally be done with `instantiate_model` from PowerModels. In order to facilitate using the `ENGINEERING` model we have introduced `instantiate_mc_model` to aid in the generation of the JuMP model. `instantiate_mc_model` will automatically convert the data model to MATHEMATICAL if necessary (notifying the user of the conversion), and pass the MATHEMATICAL model off to PowerModels' `instantiate_model` with `ref_add_arcs_transformer!` in `ref_extensions`, which is a required ref extension for PowerModelsDistribution.
"""

# ╔═╡ 4603bc56-c402-48af-9715-0722e42d90f7
pm_eng = instantiate_mc_model(eng, NFAUPowerModel, build_mc_opf)

# ╔═╡ 12ba02ff-91dd-4355-ad2b-f3ca3b61e318
"```$(pm_eng.model)```" |> Markdown.parse

# ╔═╡ 07fc5a1b-4a9e-4a10-bf7a-d2682769ac97
md"""
## Conclusion

This concludes the introduction to the `ENGINEERING` data model and conversion to the `MATHEMATICAL` model. We hope that you will find this new data model abstraction easy to use and simple to understand
"""

# ╔═╡ Cell order:
# ╟─c55b2c42-9d27-11eb-24ca-e90a5472ffbb
# ╠═2661286b-2d18-42e3-b309-7974fc2db425
# ╟─de92bc20-b125-4f3d-930a-b1da63d5cef5
# ╠═f30cd0d0-b0da-4f63-a245-568a763a93d8
# ╠═ac94e556-b544-4ca6-87bd-ff2d7a7414e7
# ╟─1e791262-261d-4756-bec4-edebe4732700
# ╠═89de80df-1dd0-4f94-a27d-74693a978059
# ╟─e22a7d3e-3f21-41b0-abf2-2c723e87c57e
# ╠═5fe2d186-40c6-46c2-823e-9401cb3b6d6c
# ╟─c884fd5b-4ecb-402b-b228-01b9a61db8bf
# ╠═a5b736f9-2776-4760-b073-d02027baef13
# ╟─29a8a560-7d5e-4929-877a-2dab17309968
# ╠═9ac6b0ed-58a8-4079-aa9f-d70024f4d4b4
# ╟─89846878-ebfc-49a5-91ef-762f692ba5ea
# ╠═8ba081b4-8444-48e0-afce-20649f7fdd01
# ╟─2a8c1536-ef5f-4e37-b30e-1bb73982d7b0
# ╠═83ce434c-0f92-46a8-8648-322084044600
# ╠═2f9c3f64-3e5d-423f-939e-8d96bd19e7fe
# ╟─f06b9592-0e16-4fd1-bf6e-624ec0cdb8fb
# ╠═b15593b9-f775-47af-aa9b-4980a4028faa
# ╟─2e4fe000-449c-4fc0-8707-4c51ac50ab03
# ╠═36ca17bc-b560-48c9-988f-255780479d05
# ╟─7e5811cf-86a2-43fb-8c94-44091f33c031
# ╠═026088f7-b1df-4647-ad81-52c6c8ea7944
# ╟─5f72ecec-bb3f-4231-93c1-719e38e17be4
# ╠═f070dbbc-1a55-48f2-aa21-2119eb573b5b
# ╟─f3458a5c-70c6-4d0c-8253-4feb0c87ee76
# ╠═542e592f-8f75-41f9-8f98-5373b431fcc6
# ╟─d0e02ee5-f57c-4286-940d-b01097c840be
# ╠═e852712c-d09a-41e1-9785-308c219c7ad8
# ╟─0eba403e-c2ba-4af3-8dcc-fe14a3c22bd7
# ╟─8d5701d1-9e91-43b5-975b-fc8e0e306e99
# ╠═745f33bf-98d7-4b27-95fd-ec6bf5e3d5cf
# ╟─fb81d3c1-a47b-46d9-b463-1766022be114
# ╠═27d5a237-4c2e-4b14-a082-f4847153c8d5
# ╟─564200e0-2080-44da-8dbf-37c5e3a38412
# ╠═8b066039-ed65-4617-b404-86952d3fc78e
# ╟─a09de8ab-c1e1-4c9f-9d86-4268ce503cae
# ╠═763a1504-f927-4dfc-a042-ac66324708a9
# ╠═e7b588e9-fb5f-4383-b5c6-bce9ed8f5a60
# ╟─f882f070-74af-4f3b-9ad1-bd7ad453778f
# ╠═f2fa779e-f62d-4053-b055-e9f0aee6c3f7
# ╟─544108d1-d1a7-4dc6-9a3c-6a10025cb8a3
# ╠═ae8daf6c-4055-497e-905c-d6d5a6373ad3
# ╟─3d86989a-a4db-4766-aec7-abc9a25ffde8
# ╠═7a86e26d-fdf6-41b6-b1dd-dc0c4a8e9d50
# ╠═82e9330b-2e94-42a1-b56d-430674b11241
# ╟─2c358d7a-5d84-4f62-a147-577bb775ac7f
# ╠═5e5dd24c-95b5-4951-8547-13e20bb0feec
# ╠═45c320dd-e0b2-46ec-af09-d83e5652026d
# ╟─723f4e26-8065-4e22-8c5c-de32d2e42b47
# ╠═58f80651-4a44-4781-a6c0-a056ec76b07e
# ╟─b1e55263-9c76-43b2-b91c-115b2f080182
# ╠═2cdd682a-b97f-49dc-8829-47e7146497b9
# ╠═74836b2c-db67-4435-a552-ba5c3e93e43f
# ╟─2bf802e2-f769-48db-adef-a5a066c565c2
# ╠═4020d67a-4a9e-4030-b7cd-e5d60e78495c
# ╟─dde2da59-63fa-4549-8829-491e5f3c675e
# ╟─2c375878-f16b-4833-9aaa-b826fbe056a2
# ╠═a7348fc8-b25e-4ecf-9fe6-314036b7df9c
# ╟─54280696-beec-47ee-8114-ea9311aca9eb
# ╠═1dc06007-b413-4cbb-b43b-2c7955ecb6d2
# ╟─96a65bdc-7ee4-4e9f-94f5-8a83f7512fe9
# ╠═545c9270-257a-46e8-970d-6a9fcf09dc17
# ╟─fda420a1-e563-4a41-84d4-2a570bf163ae
# ╠═785a1c39-aa56-4244-b7b9-31d7a31b73c3
# ╠═ad48b3ad-ede1-4109-9647-007c1d8d196c
# ╟─666223d5-802a-4164-882c-8ab1a1319d6d
# ╠═e3fe0a5e-9431-453d-b7e5-821de624a4ef
# ╟─5ab6cd16-3c24-4f60-bdf2-7f48b3e25085
# ╠═05360163-4a6f-41da-8fae-070e327c2605
# ╟─b8e01ce2-d7f7-4c73-8230-9fdeecab2c4e
# ╠═96757dde-52be-4930-9bdd-9303bb715e55
# ╠═867d8d1d-042a-42ec-843e-56b8ba2db96a
# ╟─db4d9eb9-e120-4489-a3c2-48c43f2a3c3b
# ╠═4603bc56-c402-48af-9715-0722e42d90f7
# ╟─12ba02ff-91dd-4355-ad2b-f3ca3b61e318
# ╟─07fc5a1b-4a9e-4a10-bf7a-d2682769ac97
