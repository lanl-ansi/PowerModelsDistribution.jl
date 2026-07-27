### A Pluto.jl notebook ###
# v0.19.22

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

__Note__: In v0.9, `apply_kron_reduction!` and `apply_phase_projection!` are applied by default, but can be disabled with the keyword arguments `kron_reduce=false` and `phase_project=false`, respectively in `parse_file` or `transform_data_model`.

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
Ipopt = "~1.2.0"
PowerModelsDistribution = "~0.14.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ASL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6252039f98492252f9e47c312c8ffda0e3b9e78d"
uuid = "ae81ac8f-d209-56e5-92de-9978fef736f9"
version = "0.1.3+0"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "SnoopPrecompile", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "c700cce799b51c9045473de751e9319bdd1c6e94"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.9"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "485193efd2176b88e6622a39a246f8c5b600e74e"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.6"

[[CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "a4ad7ef19d2cdc2eff57abbbe68032b1cd0bd8f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.13.0"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[Glob]]
git-tree-sha1 = "4df9f7e06108728ebf00a0a11edee4b29a482bb2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.0"

[[InfrastructureModels]]
deps = ["JuMP", "Memento"]
git-tree-sha1 = "88da90ad5d8ca541350c156bea2715f3a23836ce"
uuid = "2030c09a-7f63-5d83-885d-db604e0e9cc0"
version = "0.7.6"

[[InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[Ipopt]]
deps = ["Ipopt_jll", "LinearAlgebra", "MathOptInterface", "OpenBLAS32_jll", "SnoopPrecompile"]
git-tree-sha1 = "7690de6bc4eb8d8e3119dc707b5717326c4c0536"
uuid = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
version = "1.2.0"

[[Ipopt_jll]]
deps = ["ASL_jll", "Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "MUMPS_seq_jll", "OpenBLAS32_jll", "Pkg", "libblastrampoline_jll"]
git-tree-sha1 = "563b23f40f1c83f328daa308ce0cdf32b3a72dc4"
uuid = "9cc047cb-c261-5740-88fc-0cf96f7bdcc7"
version = "300.1400.403+1"

[[IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[JuMP]]
deps = ["LinearAlgebra", "MathOptInterface", "MutableArithmetics", "OrderedCollections", "Printf", "SnoopPrecompile", "SparseArrays"]
git-tree-sha1 = "611b9f12f02c587d860c813743e6cec6264e94d8"
uuid = "4076af6c-e467-56ae-b986-b466b2749572"
version = "1.9.0"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[METIS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1fd0a97409e418b78c53fac671cf4622efdf0f21"
uuid = "d00139f3-1899-568f-a2f0-47f597d42d70"
version = "5.1.2+0"

[[MUMPS_seq_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "METIS_jll", "OpenBLAS32_jll", "Pkg", "libblastrampoline_jll"]
git-tree-sha1 = "f429d6bbe9ad015a2477077c9e89b978b8c26558"
uuid = "d7ed1dd3-d0ae-5e8e-bfb4-87a502085b8d"
version = "500.500.101+0"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "DataStructures", "ForwardDiff", "JSON", "LinearAlgebra", "MutableArithmetics", "NaNMath", "OrderedCollections", "Printf", "SnoopPrecompile", "SparseArrays", "SpecialFunctions", "Test", "Unicode"]
git-tree-sha1 = "f219b62e601c2f2e8adb7b6c48db8a9caf381c82"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.13.1"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[Memento]]
deps = ["Dates", "Distributed", "Requires", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "bb2e8f4d9f400f6e90d57b34860f6abdc51398e5"
uuid = "f28f55f0-a522-5efc-85c2-fe41dfb9b2d9"
version = "1.4.1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3295d296288ab1a0a2528feb424b854418acff57"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.2.3"

[[NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[OpenBLAS32_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c6c2ed4b7acd2137b878eb96c68e63b76199d0f"
uuid = "656ef2d0-ae68-5445-9ca0-591084a874a2"
version = "0.3.17+0"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

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
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[PolyhedralRelaxations]]
deps = ["DataStructures", "ForwardDiff", "JuMP", "Logging", "LoggingExtras"]
git-tree-sha1 = "05f2adc696ae9a99be3de99dd8970d00a4dccefe"
uuid = "2e741578-48fa-11ea-2d62-b52c946f73a0"
version = "0.3.5"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[PowerModelsDistribution]]
deps = ["CSV", "Dates", "FilePaths", "Glob", "InfrastructureModels", "JSON", "JuMP", "LinearAlgebra", "Logging", "LoggingExtras", "PolyhedralRelaxations", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "fd2a5efc06acb1b449a985c48d4b3d8004a3b371"
uuid = "d7431456-977f-11e9-2de3-97ff7677985e"
version = "0.14.7"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

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
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "77d3c4726515dca71f6d80fbb5e251088defe305"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.18"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "6aa098ef1012364f2ede6b17bf358c7f1fbe90d4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.17"

[[StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
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
