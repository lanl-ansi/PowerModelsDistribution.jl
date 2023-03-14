### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ a3d1557a-000b-4bc8-a390-571c08e55586
begin
	using PowerModelsDistribution
	import Ipopt
	import DataFrames: DataFrame
end

# ╔═╡ e12e9ffa-1d5e-11ec-0ef4-435a3b4ee5d5
md"""
# Introduction to Explicit Neutral Models

This notebook illustrates how you can use PMD to optimize network models with an explicit neutral (EN) conductor representation. It will go through the full workflow, consisting of
- importing OpenDSS network data (and applying transformations as needed);
- adding OPF-specific data to the model;
- optimizing the model;
- inspecting the results.
"""

# ╔═╡ 947fccff-8647-48c5-9bd1-c483ee589a33
md"""
This notebook will make use of the following packages in various places
"""

# ╔═╡ 1842fe02-c874-42c7-926e-d843868eaa5f
pmd_path = joinpath(dirname(pathof(PowerModelsDistribution)), "..")

# ╔═╡ 5c82c6c2-b1b5-4914-97ec-fe4b9214d5b0
md"""
## Building a case

### Importing network data
"""

# ╔═╡ 40e870aa-ca5a-4d4d-a904-c79c419aab87
md"""
OpenDSS cases with explicit neutral conductors tend to use components in a way which is not supported directly by PMD. These cases model the grounding as a 'reactor', connected between different terminals of the same bus (i.e. from the 4th terminal to ground).

In PMD, a reactor is mapped by default to a line, which will then be connected on both ends to the same bus, which is not allowed. This problem is solved by applying the data model transformation `transform_loops!`, which will map the reactor to a shunt instead, or merge terminals if the reactor is in fact a short-circuit.
"""

# ╔═╡ 52eeb3df-1fed-404a-a2ca-e94b19e0df4d
begin
	case_path = pmd_path*"/test/data/en_validation_case_data/test_grounding.dss"
	data_eng = parse_file(case_path, transformations=[transform_loops!])
end

# ╔═╡ 08672065-d643-4769-adf3-7d640bd70f19
md"""
Note that this test case is *very* unbalanced. To obtain a more realistic scenario, we reduce the loading by a factor 3.
"""

# ╔═╡ 2b583476-e5d7-48fb-a9d2-778890a5656a
begin
	for (_,load) in data_eng["load"]
		load["pd_nom"] *= 1/3
		load["qd_nom"] *= 1/3
	end
end

# ╔═╡ 7a2e018b-8bbf-4583-8191-b5a3e70ca25d
md"""
### Adding OPF specific data

First of all, it is a good practice to remove any bounds that may have been imported from the OpenDSS network data. These can be default bounds, which might not make sense for the specific case and can lead to infeasibility. For example, this occurs for IEEE13, where the default OpenDSS line ratings are too tight for the base case power flows.
"""

# ╔═╡ 57778559-f76e-438f-863e-78fb255f9eb1
remove_all_bounds!(data_eng)

# ╔═╡ a7a2c1e9-8477-4f39-9fb5-6320a8542fce
md"""
#### Voltage bounds

Voltage bounds are naturally expressed *between* two terminals.

For example, if a wye-connected load is connected between phase `a` and the neutral `n`, we want to apply the bounds `lb <= |Ua-Un| <= ub`, in order to ensure the voltage across the load does not get too low or high. If the neutral is not modeled and assumed to be grounded everywhere, i.e. `Un=0`, then this can be done equivalently by bounding `Ua` directtly.

PMD comes with several data model transformations which allow you to quickly apply voltage bounds. The bounds are specified in per unit, relative to the local voltage base (phase-to-ground). If this was not the case, it would be difficult to specify bounds for multiple voltage zones at once.

First consider `add_bus_absolute_vbounds!`, which allows you to specify *absolute* bounds, i.e. for each terminal individually. This will set the bus properties `vm_lb` and `vm_ub` appropriately.
"""

# ╔═╡ 1ebfb921-817b-4615-8cab-a6e8af7cb0e3
add_bus_absolute_vbounds!(
	data_eng,
	phase_lb_pu = 0.8,
	phase_ub_pu = 1.2,
	neutral_ub_pu = 0.1,
)

# ╔═╡ 5b8def6d-f786-439d-ac89-50054b2591dc
md"""
The `ENGINEERING` data model includes special properties to easily apply symmetrical bounds for three-phase buses (refer to the [documentation](https://lanl-ansi.github.io/PowerModelsDistribution.jl/stable/manual/eng-data-model.html#Special-Case:-three-phase-bus) for more details). These properties can be populated easily with `add_bus_pn_pp_ng_vbounds!`, which specifies
- phase-to-neutral (pn)
- phase-to-phase (pp)
- neutral-to-ground (ng)
bounds. Note that since the voltage base is in phase-to-ground, the `pp` bounds have to multiplied by `sqrt(3)` for a three-phase network.
"""

# ╔═╡ 4d512fe7-c114-4a9a-973e-77eb680a8745
add_bus_pn_pp_ng_vbounds!(
	data_eng, [1:3...], 4,
	pn_lb_pu = 0.9,
	pn_ub_pu = 1.1,
	pp_lb_pu = 0.9*sqrt(3),
	pp_ub_pu = 1.1*sqrt(3),
	ng_ub_pu = 0.1
)

# ╔═╡ 2779e6c1-fe5e-45b9-a0a7-5b4bbfef6953
md"""
Often, the true goal of voltage bounds is to protect the connected 'units', i.e. loads, generators etc. Therefore, it can suffice to apply voltage bounds only to those buses and terminals which actually have a unit connected to them. This is what `add_unit_vbounds!` is for. The `delta_multiplier` is the correction factor which is applied to the specified bounds if a unit has `configuration=DELTA`.
"""

# ╔═╡ ecaa0ad3-80ca-4121-a1d0-b00117d6cebc
add_unit_vbounds!(
	data_eng,
	lb_pu = 0.91,
	ub_pu = 1.09,
	delta_multiplier = sqrt(3),
	unit_comp_types = ["load"]
)

# ╔═╡ 259b0a46-f47a-46ea-b31b-a81ee19ac040
md"""
Whether it is desirable to apply as many bounds as possible, or only ones which are really needed, will depend on the chosen formulation later on. In any case, these data model transformations make it easy for the user to apply any desired bounds.

You might have noticed that applying all these transformations will lead to redundant constraints. This is resolved in the data model transformation, which determines a set of *absolute* and *pairwise* voltage constraints which imply all other ones.
"""

# ╔═╡ 73547628-d5dc-4ac6-8a14-51e67e87a048
md"""
For example, let's inspect bus `b2`. This bus has a single-phase load between terminals 1 and 4, so only for that pair the lower bound should be 0.91 instead of 0.9. For an explanation of how the `vm_pair_lb` property is structured, refer to the [documentation](https://lanl-ansi.github.io/PowerModelsDistribution.jl/stable/manual/eng-data-model.html#Buses-(bus)).
"""

# ╔═╡ 05faff5b-f590-4e9b-9959-3770469bbf96
begin
	data_math_tmp = transform_data_model(data_eng, kron_reduce=false, phase_project=false)
	b2_index = data_math_tmp["bus_lookup"]["b2"]
	b2_math = data_math_tmp["bus"]["$b2_index"]
	b2_math["vm_pair_lb"]
end

# ╔═╡ 579a151f-680d-431a-8db1-b35e783f16ce
md"""
### Adding a generator

So far, we imported a power flow case from OpenDSS. This means that the problem is fully determined, i.e. there is no degrees of freedom left over which to optimize. Therefore, we will add an additional generator to the problem.
"""

# ╔═╡ 71cc07b7-b232-42fa-8f4f-2c95feda04b8
begin
	data_eng["generator"] = Dict{String,Any}()
	data_eng["generator"]["g1"] = Dict{String,Any}(
		"status" => ENABLED,
		"bus" => "b2",
		"configuration" => WYE,
		"connections" => [2,4],
		"pg_lb" => [0.0],
		"pg_ub" => [20.0],
		"qg_lb" => [0.0],
		"qg_ub" => [0.0],
		"cost_pg_parameters" => fill(0.0, 3)
	)
end

# ╔═╡ fafe3980-f88b-4ad5-8297-2a97011f6eeb
md"""
## Solving a simple OPF problem
"""

# ╔═╡ 54bb64b1-4544-416b-94fb-f7500d60dd9e
md"""
First of all, we apply the data model transformation to go from `ENGINEERING` to `MATHEMATICAL`. We need to pass the flags `kron_reduce=false` and `phase_project=false`; these activate data model modifications which are required for some formulations, but not for the EN ones.
"""

# ╔═╡ 3d31b6c8-63a3-466b-8ab0-03164b9ec43c
data_math = transform_data_model(data_eng, kron_reduce=false, phase_project=false)

# ╔═╡ c7bbb931-c01f-41cd-8d71-48b5c6d51b7b
md"""
Before solving the problem, it is important to add initialization values for the voltage variables. Failing to do so will almost always result in solver issues.

For a single-phase equivalent network, this is very easy to do; simply initialize each voltage variable to `1.0+im*0.0`, which usually corresponds to the voltage profile without any network load and ignoring all shunts and linecharging (refered to as 'no-load voltage').

The equivalent in a general, multi-conductor network is less trivial; for example, transformers can shift the no-load voltage angle of individual terminals in complicated ways. Therefore, we developed the method `add_start_vrvi!`, which will infer the no-load voltage for each terminal in the network and add initialization properties to the `MATHEMATICAL` data model.
"""

# ╔═╡ 74cd872c-3c68-43f5-97a0-03b1d4ff2058
add_start_vrvi!(data_math)

# ╔═╡ d3de88eb-0ef3-4253-a167-869d6d6340e8
md"""
Now we are ready to solve the OPF problem.
"""

# ╔═╡ 8b7e8f12-ef0c-4dad-8e16-a338ea09cee6
res = solve_mc_opf(data_math, IVRENPowerModel, Ipopt.Optimizer)

# ╔═╡ c7c0c3bc-ecca-46f1-be97-f07a66efd1c5
md"""
## Inspecting results
"""

# ╔═╡ e998a8a2-c3ec-4652-95af-f6c56dab5e93
md"""
The result dictionary contains the solutions.
"""

# ╔═╡ a029214f-62e3-4b86-ad7d-ddd55d6556a1
sol_math = res["solution"]

# ╔═╡ 54038baf-974c-4d37-920e-845e7346bd3c
md"""
Next, we can transform the solution back to the `ENGINEERING` data model.
"""

# ╔═╡ 6068bcbd-8217-4baa-97e4-5b25c2f5383a
sol_eng = transform_solution(sol_math, data_math)

# ╔═╡ f849592e-d744-467b-94ff-fa3bca1c93f7
md"""
For example, let's inspect the dispatched active generator power.
"""

# ╔═╡ 28f5437b-3c61-4315-a8dd-98f3c893cfca
sol_eng["generator"]["g1"]["pg"]

# ╔═╡ db43ce4b-a1fa-4e8c-aece-8656fc29d1e0
md"""
The active generator power bounds are not active, `0 <= 5.15 <= 10.0`. This means some other bound is active, because since `g1` has zero cost, it is normally more optimal to dispatch more active power, which the source bus generator will pay the default price for.

So, let's inspect the voltage at bus `b2`. We will obtain the voltage base through 'data_math', which will allow us to inspect the voltage in pu. This will make it easier to compare the values against the bounds we specified before.
"""

# ╔═╡ 79a57586-4b18-4722-9cbf-a47bf514ae62
begin
	vbase_b2 = data_math["bus"][string(data_math["bus_lookup"]["b1"])]["vbase"]
	v_b2_pu = (sol_eng["bus"]["b2"]["vr"]+im*sol_eng["bus"]["b2"]["vi"])./vbase_b2
	vm_b2_pu = abs.(v_b2_pu)
end

# ╔═╡ b0764438-e6f7-43f8-9a5a-daaaa0d5d1b0
md"""
As it turns out, the voltage magnitude constraint om the neutral terminal is binding, `vm_b2_pu[4]=0.1`.
"""

# ╔═╡ c8a9adc2-0f0d-473d-a339-bcc9b54088ad
md"""
## Available formulations
"""

# ╔═╡ 353386cd-be39-4d27-ad3a-83e110d266ab
md"""
We create a `results` dictionary which will collect the result dictionary generated by each formulation. We will end this section with a comparison of them.

There are several formulations available which support explicit neutrals. The main distinction is whether the flow variables are current or power variables.
"""

# ╔═╡ b29de9b9-3346-4676-91fa-4f9147d4eb4d
results = Dict{String,Any}()

# ╔═╡ d41a36ed-b8a5-4e3f-941b-8431e0b6a5cd
md"""
### Current flow variables

The preferred class of exact formulations are `IVR`, i.e. with current flow variables (I) and rectangular voltage variables (VR).

`IVRENPowerModel` is a non-linear formulation.
"""

# ╔═╡ b99d3796-740c-4f8b-ba99-1dae4f774770
results["IVRENPowerModel"] = solve_mc_opf(data_math, IVRENPowerModel, Ipopt.Optimizer)

# ╔═╡ 5264da2d-5ad5-4b1a-801b-506dd90f615a
md"""
`IVRQuadraticENPowerModel` is an equivalent quadratic formulation.
"""

# ╔═╡ 90b9f912-913f-4211-b658-5266a429e9cc
results["IVRQuadraticENPowerModel"] = solve_mc_opf(data_math, IVRQuadraticENPowerModel, Ipopt.Optimizer)

# ╔═╡ c607215e-bbd4-4b31-af74-a744a48f0467
md"""
`Reduced` models only create explicit series current variables, and create the total current variables as linear expressions of those. Since branches tend to be the dominate component in number, this can lead to a big reduction in the number of variables.

`IVRReducedENPowerModel` is the branch-reduced version of `IVRENPowerModel`.
"""

# ╔═╡ af801bae-3c39-488d-a1a1-35b75acd4e92
results["IVRReducedENPowerModel"] = solve_mc_opf(data_math, IVRReducedENPowerModel, Ipopt.Optimizer)

# ╔═╡ 4586c0b6-d2d7-4696-a95d-dd19fe04f33c
md"""
`IVRReducedQuadraticENPowerModel` is the branch-reduced version of `IVRQuadraticENPowerModel`.
"""

# ╔═╡ ad755f03-eb64-4db0-b861-305e8d52ec48
results["IVRReducedQuadraticENPowerModel"] = solve_mc_opf(data_math, IVRReducedQuadraticENPowerModel, Ipopt.Optimizer)

# ╔═╡ 63893e92-e646-4fa0-a001-a12966f25206
md"""
### Power flow variables

We also included a single formulation with power flow variables. Because it also has rectangular voltage variables and legacy reasons, it is referred to as `ACRENPowerModel`.
"""

# ╔═╡ 1705a2ad-c69f-4728-9d94-14a5c259f95c
results["ACRENPowerModel"] = solve_mc_opf(data_math, ACRENPowerModel, Ipopt.Optimizer)

# ╔═╡ 7fc51410-113b-45d7-b802-e90019ff522d
md"""
We mentioned before that the `IVR` formulations are preferred for `EN` models.
This is because in `EN` models, the voltage magnitude cannot be bounded below for some terminals (the ones belonging to the neutral conductor).
In a formulation with power flow variables, this means that KCL cannot be enforced for those terminals.

In short, ACR allows non-physical groundings of the neutral conductor, and is therefore a relaxation of the original problem. To demonstrate this, let's create a summary of the obtained objective values and generator set points.
"""

# ╔═╡ b926e396-e2c2-4282-a83b-befa1eea0a49
begin
	forms = sort([keys(results)...])
	sol_engs = Dict(f=>transform_solution(results[f]["solution"], data_math) for f in forms)
	DataFrame(
		"formulation" => forms,
		"objective value" => [results[f]["objective"] for f in forms],
		"g1 pg" => [sol_engs[f]["generator"]["g1"]["pg"][1] for f in forms],
	)
end

# ╔═╡ 3f803a29-caf4-4d33-a25d-5f972377c778
md"""
This table illustrates that `ACRENPowerModel` is a relaxation  of the other ones. However, note that it is not guaranteed that the objective value will be lower, because all problems are only solved to local optimality.  And in fact, the `ACR` solution is very sensitive to changes in the initialization. The optional virtual groundings seem to introduce many potential local optima.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Ipopt = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
PowerModelsDistribution = "d7431456-977f-11e9-2de3-97ff7677985e"

[compat]
DataFrames = "~1.5.0"
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

[[Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "aa51303df86f8626a962fccb878430cdb0a97eee"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.5.0"

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

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

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

[[InvertedIndices]]
git-tree-sha1 = "82aec7a3dd64f4d9584659dc0b62ef7db2ef3e19"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.2.0"

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

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

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

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

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

[[PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "96f6db03ab535bdb901300f88335257b0018689d"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.2"

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

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

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

[[StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

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
# ╟─e12e9ffa-1d5e-11ec-0ef4-435a3b4ee5d5
# ╟─947fccff-8647-48c5-9bd1-c483ee589a33
# ╠═a3d1557a-000b-4bc8-a390-571c08e55586
# ╠═1842fe02-c874-42c7-926e-d843868eaa5f
# ╟─5c82c6c2-b1b5-4914-97ec-fe4b9214d5b0
# ╟─40e870aa-ca5a-4d4d-a904-c79c419aab87
# ╠═52eeb3df-1fed-404a-a2ca-e94b19e0df4d
# ╟─08672065-d643-4769-adf3-7d640bd70f19
# ╠═2b583476-e5d7-48fb-a9d2-778890a5656a
# ╟─7a2e018b-8bbf-4583-8191-b5a3e70ca25d
# ╠═57778559-f76e-438f-863e-78fb255f9eb1
# ╟─a7a2c1e9-8477-4f39-9fb5-6320a8542fce
# ╠═1ebfb921-817b-4615-8cab-a6e8af7cb0e3
# ╟─5b8def6d-f786-439d-ac89-50054b2591dc
# ╠═4d512fe7-c114-4a9a-973e-77eb680a8745
# ╟─2779e6c1-fe5e-45b9-a0a7-5b4bbfef6953
# ╠═ecaa0ad3-80ca-4121-a1d0-b00117d6cebc
# ╟─259b0a46-f47a-46ea-b31b-a81ee19ac040
# ╟─73547628-d5dc-4ac6-8a14-51e67e87a048
# ╠═05faff5b-f590-4e9b-9959-3770469bbf96
# ╟─579a151f-680d-431a-8db1-b35e783f16ce
# ╠═71cc07b7-b232-42fa-8f4f-2c95feda04b8
# ╟─fafe3980-f88b-4ad5-8297-2a97011f6eeb
# ╟─54bb64b1-4544-416b-94fb-f7500d60dd9e
# ╠═3d31b6c8-63a3-466b-8ab0-03164b9ec43c
# ╟─c7bbb931-c01f-41cd-8d71-48b5c6d51b7b
# ╠═74cd872c-3c68-43f5-97a0-03b1d4ff2058
# ╟─d3de88eb-0ef3-4253-a167-869d6d6340e8
# ╠═8b7e8f12-ef0c-4dad-8e16-a338ea09cee6
# ╟─c7c0c3bc-ecca-46f1-be97-f07a66efd1c5
# ╟─e998a8a2-c3ec-4652-95af-f6c56dab5e93
# ╠═a029214f-62e3-4b86-ad7d-ddd55d6556a1
# ╟─54038baf-974c-4d37-920e-845e7346bd3c
# ╠═6068bcbd-8217-4baa-97e4-5b25c2f5383a
# ╟─f849592e-d744-467b-94ff-fa3bca1c93f7
# ╠═28f5437b-3c61-4315-a8dd-98f3c893cfca
# ╟─db43ce4b-a1fa-4e8c-aece-8656fc29d1e0
# ╠═79a57586-4b18-4722-9cbf-a47bf514ae62
# ╟─b0764438-e6f7-43f8-9a5a-daaaa0d5d1b0
# ╟─c8a9adc2-0f0d-473d-a339-bcc9b54088ad
# ╟─353386cd-be39-4d27-ad3a-83e110d266ab
# ╠═b29de9b9-3346-4676-91fa-4f9147d4eb4d
# ╟─d41a36ed-b8a5-4e3f-941b-8431e0b6a5cd
# ╠═b99d3796-740c-4f8b-ba99-1dae4f774770
# ╟─5264da2d-5ad5-4b1a-801b-506dd90f615a
# ╠═90b9f912-913f-4211-b658-5266a429e9cc
# ╟─c607215e-bbd4-4b31-af74-a744a48f0467
# ╠═af801bae-3c39-488d-a1a1-35b75acd4e92
# ╟─4586c0b6-d2d7-4696-a95d-dd19fe04f33c
# ╠═ad755f03-eb64-4db0-b861-305e8d52ec48
# ╟─63893e92-e646-4fa0-a001-a12966f25206
# ╠═1705a2ad-c69f-4728-9d94-14a5c259f95c
# ╟─7fc51410-113b-45d7-b802-e90019ff522d
# ╠═b926e396-e2c2-4282-a83b-befa1eea0a49
# ╟─3f803a29-caf4-4d33-a25d-5f972377c778
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
