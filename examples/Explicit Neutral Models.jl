### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 22cd4c84-7d9d-43dd-9b5c-519016b4ff41
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.Registry.update()
	Pkg.add([
			Pkg.PackageSpec(;name="PowerModelsDistribution", version="0.12"),
			Pkg.PackageSpec(;name="Ipopt", version="0.6.5"),
			Pkg.PackageSpec(;name="DataFrames", version="1.1.1"),
			])
end


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

## Julia environment setup
The following code block will setup a Julia environment for you with the correct versions of packages for this Pluto notebook
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

# ╔═╡ Cell order:
# ╟─e12e9ffa-1d5e-11ec-0ef4-435a3b4ee5d5
# ╠═22cd4c84-7d9d-43dd-9b5c-519016b4ff41
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
