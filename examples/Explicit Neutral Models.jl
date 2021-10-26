### A Pluto.jl notebook ###
# v0.16.4

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
	data_math_tmp = transform_data_model(data_eng, kron_reduced=false, phase_projected=false)
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
First of all, we apply the data model transformation to go from `ENGINEERING` to `MATHEMATICAL`. We need to pass the flags `kron_reduced=false` and `phase_projected=false`; these activate data model modifications which are required for some formulations, but not for the EN ones.
"""

# ╔═╡ 3d31b6c8-63a3-466b-8ab0-03164b9ec43c
data_math = transform_data_model(data_eng, kron_reduced=false, phase_projected=false)

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
DataFrames = "~1.2.2"
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

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d785f42445b63fc86caa08bb9a9351008be9b765"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.2.2"

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

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

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

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

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

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "d940010be611ee9d67064fe559edbb305f8cc0eb"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.2.3"

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

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

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

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

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
