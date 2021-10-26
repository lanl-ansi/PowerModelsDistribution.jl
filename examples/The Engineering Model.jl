### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# ╔═╡ f30cd0d0-b0da-4f63-a245-568a763a93d8
begin
	using PowerModelsDistribution
	using Ipopt
end

# ╔═╡ c55b2c42-9d27-11eb-24ca-e90a5472ffbb
md"""
# Introduction to the PowerModelsDistribution Data Models

In this notebook we introduce the engineering data model added to PowerModelsDistribution in version v0.9.0. We will give several examples of how to use this new data model directly, including new transformations that have become easier with its introduction, how to convert it to the the lower-level mathematical model that was previously the only user interface we offered, and how to get various types of results using this new model.
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Ipopt = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
PowerModelsDistribution = "d7431456-977f-11e9-2de3-97ff7677985e"

[compat]
Ipopt = "~0.7.0"
PowerModelsDistribution = "~0.12.0"
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

[[FilePathsBase]]
deps = ["Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "7fb0eaac190a7a68a56d2407a6beff1142daf844"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.12"

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

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a193d6ad9c45ada72c14b731a318bedd3c2f00cf"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.3.0"

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

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "f45b34656397a1f6e729901dc9ef679610bd12b5"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.8"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

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

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

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

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

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
# ╟─c55b2c42-9d27-11eb-24ca-e90a5472ffbb
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
