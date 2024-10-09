### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ a3d1557a-000b-4bc8-a390-571c08e55586
begin
	import JuMP
	import PlutoUI
	using PowerModelsDistribution
	import InfrastructureModels
	import Ipopt
	import Plots
end

# ╔═╡ e12e9ffa-1d5e-11ec-0ef4-435a3b4ee5d5
md"""
# Extending PowerModelsDistribution.jl

PowerModelsDistribution.jl (PMD) can be extended in many ways, most of which would require knowledge of the code base. This tutorial is not about those extensions. Instead, we will use generators as an interface to connect device models to the network model of PMD. This can be done without modiyfing PMD itself, and is actually a very common use case.

For this tutorial, we will develop an EV charging model. The aim is to optimize the EV charging schedules, minimizing the user dissatisfaction.

First, we will develop a model that does not rely on PMD, but simply limits each EV to charging at a rate of at most 5 kW. We refer to this model as 'copperplate', since it does not represent the limitations of the network. This can also serve as a simple tutorial on how you can use `JuMP` to create and solve optimization models.

Secondly, we leverage PMD to model the distribution network and its technical constraints. We lift the conservative 5 kW charging constraint, and allow any charging rate as long as the voltage quality is maintained.
"""

# ╔═╡ 947fccff-8647-48c5-9bd1-c483ee589a33
md"""
This notebook will make use of the following packages in various places. Note that we do not apply `using` to any packages; we avoid this so you do not get confused where certain methods are coming from. Also, we introduce the shorthand `PMD` for `PowerModelsDistribution`.
"""

# ╔═╡ f1abaf46-59bc-4ad3-8054-0509e37c959c
const PMD = PowerModelsDistribution

# ╔═╡ 6e7eae9b-de28-4069-bdad-a8d4033c533c
md"""
The variable `NOTEBOOK_DIR` contains the path to the folder in which you placed this Pluto script. The script expects the network data file at `$ROOT_DIR/resources/lvtestcase_notrans.dss`.
"""

# ╔═╡ 1842fe02-c874-42c7-926e-d843868eaa5f
NOTEBOOK_DIR = join(split(split(@__FILE__, "#==#")[1], "/")[1:end-1], "/")


# ╔═╡ 550efc60-e19e-46d3-bce3-090b17ccbb30
md"""
Note that in what follows, we surpress the output of code blocks by adding a ';' at the end. The only purpose of this is to not clutter this notebook. Doing this in a regular script is not a good practice.

"""

# ╔═╡ 7569cb6f-1ad2-465b-b0cf-9edccb0e617a
md"""
## Problem description

Before building any optimization model, let's create a dummy problem.

We have a set of EVs for which we want to optimize the charging schedule.
"""

# ╔═╡ 6f96134f-fa1e-47e9-8047-9814e5b911ec
EVs = collect(1:10); # set of EVs

# ╔═╡ f2c56f02-f6c3-4576-ac6d-81dc6e7b20b3
md"""
Define a set of `timestamps` for the problem, defined in hours. Next, we convert this to a set of indexed timesteps, `K`, with an associated duration `D_k`. Note the notation `_k`, which indicates that this variable is indexed over K.
"""

# ╔═╡ da345583-51b6-4537-a540-e1100ff90d06
begin
	timestamps = collect(0:1:10) # timestamps
	D_k = timestamps[2:end].-timestamps[1:end-1] # duration of each timestep
	K = 1:length(D_k) # set of timesteps
end;

# ╔═╡ 38924354-4977-46e8-80d8-256d11a10242
md"""
Define a set of EV parameters, indexed over e $\in$ EVs.
"""

# ╔═╡ 78ba1c9e-bb56-4158-a67b-e77770f73645
begin
	Emax_e = fill(82.0, length(EVs)) # maximum SoC in KWh
	E0_e = Emax_e.*(1/length(EVs)*[0:length(EVs)-1...]) # initial SoC in KWh
	Pmax = 5.0 # maximum charge in KW
	Pmax_e = fill(Pmax, length(EVs)) # all EVs have the same maximum charge
end;

# ╔═╡ bb8be0dc-d8b3-4026-a9ea-1b7d1b04b5df
md"""
The objective is to minimize the dissatisfaction (or maximize satisfaction) of the EV owners. Per timestep `k` and EV `e`, the dissatisfaction goes linearly from 1 to 0 as the state-of-charge (SoC) goes from epty to full. See below for an illustration.
"""

# ╔═╡ bfaf5d7d-1042-45fd-bd22-0d2a78bf4d78
begin
	Plots.plot(xticks=([0,1],["empty","full"]), xlabel="SoC")
	Plots.plot!(ylabel="dissatisfaction")
	Plots.plot!([0,1], [1,0], label="", linewidth=2)
end

# ╔═╡ 7847a75f-28f2-472f-9f0a-0d40e7a40cf9
md"""
## Copperplate optimization model

First, create a blank JuMP optimization model.

We prefix all variables with `cp_`, because Pluto requires that variables are only defined once. For the network-constrained model, we will use `nc_`.
"""

# ╔═╡ 8a7de573-2e02-4332-b416-0e08dcb8c6d1
cp_model = JuMP.Model();

# ╔═╡ 72c124c1-fc0b-40f1-a7c9-7e475c3f74d7
md"""
Create variables for the charge rate and SoC.
"""

# ╔═╡ e6ac4be0-e524-4341-a1ec-28e16a3b2bfb
begin
	# charge for EV e and timestep k, in kW
	JuMP.@variable(cp_model, 0 <= cp_P_ek[e in EVs, k in K])
	# SoC for EV e at end of timestep k, in kWh
	JuMP.@variable(cp_model, 0 <= cp_E_ek[e in EVs, k in K] <= Emax_e[e]);
end;

# ╔═╡ cd654189-90fe-49e2-a0b9-c9f4e3c3f97c
md"""
Link the SoC to the charge rate.
"""

# ╔═╡ 6ff8cb3b-d153-4639-8eda-3a298860a905
begin
	# for the first timestemp, use initial SoC E0_e
	JuMP.@constraint(cp_model, [e in EVs],
		cp_E_ek[e,1] == E0_e[e] + D_k[1]*cp_P_ek[e,1]
	)
	# for the other timestemps, use the SoC at the preceding timestep
	JuMP.@constraint(cp_model, [e in EVs, k in K[2:end]],
		cp_E_ek[e,k] == cp_E_ek[e,k-1] + D_k[k]*cp_P_ek[e,k]
	)
end;

# ╔═╡ 8e0fff25-e6fd-439b-96ec-e5a13c9cff24
md"""
In the copperplate model, we do not model the network. Instead, we limit the charging of each EV to at most 5 kW, in the hope this will ensure the voltage quality.
"""

# ╔═╡ 5c3b70e4-99d0-4dd0-a7ef-d62bfa904bb7
JuMP.@constraint(cp_model, [e in EVs, k in K], cp_P_ek[e,k]<=Pmax_e[e]);

# ╔═╡ eb8479f7-99c6-49c9-b767-43a8d31cc502
md"""
Add the objective to the model.
"""

# ╔═╡ c40088d6-78b7-44dc-ac4f-ab5357490818
begin
	# define the dissatisfaction for each EV e at each timestep k
	cp_dissatisfaction_ek = [
		(Emax_e[e]-cp_E_ek[e,k])/Emax_e[e]
		for e in EVs, k in K]
	# the objective is to minimize the total dissatisfaction across EVs and timesteps,
	# taking into account the duration of each step
	JuMP.@objective(cp_model, Min,
		sum(
			D_k[k]*cp_dissatisfaction_ek[e,k]
			for e in EVs, k in K)
		)
end;

# ╔═╡ 09ae6158-3bd9-4d95-8e1e-029f308c7a91
md"""
Set a solver and solve the problem.
"""

# ╔═╡ fb7f7904-a219-46a5-ae0b-f728dc4a3fd9
begin
	# set the optimizer
	JuMP.set_optimizer(cp_model, Ipopt.Optimizer)
	# solve the problem
	JuMP.optimize!(cp_model)
	# inspect the termination status
	JuMP.termination_status(cp_model)
end

# ╔═╡ 5988c7bd-6419-4d46-9a0c-6c83f5bc4c25
md"""
The figure below shows how the SoC of the EVs will evolve with the obtained charging schedule. As could be expected, each EV will simply charge at 5 kW until the maximum energy capapcity is reached.
"""

# ╔═╡ 49539850-e2b1-4976-95dc-1a39a160747a
begin
	cp_E_ek_vals = JuMP.value.(cp_E_ek.data)
	Plots.plot(legend=:none, title="", xlabel="time [h]", ylabel="SoC[kWh]", ylim=[0, maximum(Emax_e)])
	for e in EVs
		Plots.plot!(timestamps, [E0_e[e], cp_E_ek_vals[e,:]...], markershape=:circle, markersize=3)
	end
	Plots.plot!()
end

# ╔═╡ ce8ac568-c3b5-4413-936f-b42ae2b9234d
md"""
## Network-constrained optimization model

Now we will use PMD to include a full network model in the optimization model. Let's first discuss conceptually how this works.

PMD is unaware of our EV model, so how can we get it to connect our EVs? This is achieved by using generators as a generic interface for external devices. For each EV, create a generator with the same connection settings. By linking the generator power variables to our EV model with a constraint, a combined model is obtained. This approach is illustrated in the figure below.
"""

# ╔═╡ 2b2d5e60-f4ce-4813-9608-5090bf6c5213
PlutoUI.LocalResource("resources/extension_summary.svg", "width"=>"700px")

# ╔═╡ 423c6e71-9978-4845-89e4-803f32098af4
md"""
**Build PMD data model**

Let's start by importing a `dss` network file as a PMD data model. We remove any bounds that might be imported from the `dss` file, because we will set our own anyways (otherwise default bounds might get imported, which can cause unexpected feasibility issues).
"""

# ╔═╡ 0be0abd5-f0ad-4d15-b6fb-d59d74e6058c
data_eng = PMD.parse_file(
	"$NOTEBOOK_DIR/resources/lvtestcase_notrans.dss",
	transformations=[remove_all_bounds!]
);

# ╔═╡ 04ef243c-2586-406e-9528-140d732c28f7
md"""
This network has about 900 lines. However, it is possible to obtain a reduced network model which is very similar and sometimes equivalent (when linecharging is negligible, which is often the case).
"""

# ╔═╡ 51d2dae6-a59c-4014-bdd6-886864847683
PMD.reduce_lines!(data_eng);

# ╔═╡ 2bd488f2-1a73-4c8e-a2bd-42a6fbd76d94
md"""
PMD uses two data models: a high-level `ENGINEERING` one, and a low-level `MATHEMATICAL` one. `ENGINEERING` is in SI-units, whilst `MATHEMATICAL` is in pu.

Since the optimization model is generated from the `MATHEMATICAL` model, we want to specify explicitly what the power base should be. Set it to 1 kW, so the unit is the same as in our EV charging model.
"""

# ╔═╡ 1d878788-13bd-4090-9ed0-b9eb77a8575d
data_eng["settings"]["sbase_default"] = 1.0*1E3/data_eng["settings"]["power_scale_factor"];

# ╔═╡ 548ffacb-af94-46c1-bcea-3695f98b4516
md"""
We require that in per unit, the phase voltage magnitude $|U_p|$ and neutral voltage magnitude $|U_n|$ should obey

$0.9 \leq |U_p| \leq 1.1, \hspace{3em} |U_n| \leq 0.1$.

We can easily add these bounds to the data model with `PMD.add_bus_absolute_vbounds!`. Note that PMD can also constrain phase-to-neutral voltages instead of only absolute ones, but we ommit that here for simplicity.
"""

# ╔═╡ 0f97d7aa-cdbe-454c-83f0-978964c83b2a
PMD.add_bus_absolute_vbounds!(
	data_eng,
	phase_lb_pu=0.9,
	phase_ub_pu=1.1,
	neutral_ub_pu=0.1
);

# ╔═╡ 161c2052-7fbd-4993-a21f-3ed14b750619
md"""
So far, the network model contains only single-period data. Since we actually need a multi-period model, we add a time series to it which has the same amount of steps as our EV charging model.

Whilst we are at it, we apply this time series to modulate the consumption of the loads in the network.
"""

# ╔═╡ 8d24de5c-749c-4b0c-acaa-25d010a33844
begin
	# add a new time series to the data model
	data_eng["time_series"] = Dict{String, Any}()
	data_eng["time_series"]["normalized_load_profile"] = Dict{String, Any}(
		"replace" => false,
		"time" => K,
		"values" => 0.2*cos.((pi/2/maximum(K)).*K)
	)
	# attach a reference to each load, so that the consumption will be scaled
	# by the profile we created
	for (_,load) in data_eng["load"]
		load["time_series"] = Dict(
			"pd_nom"=>"normalized_load_profile",
			"qd_nom"=>"normalized_load_profile"
		)
	end
end

# ╔═╡ eb8c0d08-5ac0-446f-9447-944399eacdb0
md"""
We need to add a generator for each EV, and specify the connection settings. In the test case we imported, LVTestCase, each load represents a household with a single-phase connection. We now associate each EV with a household, and give it the same bus and phase connection.
"""

# ╔═╡ 8cf4316e-1cfa-4cf0-abff-1438edcb6d36
begin
	# load to which each EV belongs
	load_e = "load".*string.(EVs)
	# bus to which each EV is connected (same as associated load)
	bus_e = [data_eng["load"][id]["bus"] for id in load_e]
	# phase terminal for each EV (same as associated load)
	phase_e = [data_eng["load"][id]["connections"][1] for id in load_e]
end;

# ╔═╡ f4344910-d95e-4f33-9b16-3f52c9ced4ac
md"""
Now we are ready to add the generators to the data model.
"""

# ╔═╡ 97017a8b-f86b-49b4-ac83-66303df1f63c
begin
	data_eng["generator"] = Dict{String, Any}()
	for e in EVs
		data_eng["generator"]["EV_gen_$e"] = Dict{String, Any}(
			"status" => ENABLED,
			"bus" => bus_e[e],
			"connections" => [phase_e[e], 4],
			"configuration" => WYE,
		)
	end
end;

# ╔═╡ 3a1abed7-d111-49e1-bcfb-ee9af48d4a6d
md"""
Transform the `ENGINEERING` data model to a `MATHEMATICAL ones`, and dont forget the `multinetwork=true` flag.
"""

# ╔═╡ 52c40603-42fe-45d3-840c-f530fc3951f2
data_math_mn = transform_data_model(data_eng, multinetwork=true);

# ╔═╡ b2e5c25e-e497-4080-ba74-2cfa8e6c05d4
md"""
Before solving the problem, it is important to add initialization values for the voltage variables. Failing to do so will almost always result in solver issues.
"""

# ╔═╡ 49e89e0f-480d-4dc1-850a-4ad9c4ae0ebd
add_start_vrvi!(data_math_mn);

# ╔═╡ 8ac591c5-e3d4-42bd-994c-7bc885536824
md"""
**Build PMD optimization model**

Generate the PMD optimization model based on the data model.
"""

# ╔═╡ 11ba8af0-50d3-4f5d-8c2b-4a1af9e5f5d5
pm = instantiate_mc_model(data_math_mn, IVRUPowerModel, build_mn_mc_opf);

# ╔═╡ 5ecc3cc4-0f8e-4ec5-8e2f-47a5871f1304
md"""
**Add EV charging model**

Start by extracting the JuMP model itself.
"""

# ╔═╡ bbba2730-31a7-4c27-ab88-c05f358b99b6
nc_model = pm.model;

# ╔═╡ ff40f2b8-1f77-4e2b-a0f8-27415a8e3978
md"""
Add the EV charging model to it. The code below is identical to the model in the previous section, except for the prefix `nc_` and the omission of the charge rate limit (`nc_P_ek[e,k]<=Pmax_e[e]`).
"""

# ╔═╡ 12d8fc45-a242-4c53-ac31-161af6331f0a
begin
	# charge for EV e and timestep k, in kW
	JuMP.@variable(nc_model, 0 <= nc_P_ek[e in EVs, k in K])
	# SoC for EV e at end of timestep k, in kWh
	JuMP.@variable(nc_model, 0 <= nc_E_ek[e in EVs, k in K] <= Emax_e[e]);

	# relate SoC to charge
	# for the first timestemp, use initial SoC E0_e
	JuMP.@constraint(nc_model, [e in EVs],
		nc_E_ek[e,1] == E0_e[e] + D_k[1]*nc_P_ek[e,1]
	)
	# for the other timestemps, use the SoC at the preceding timestep
	JuMP.@constraint(nc_model, [e in EVs, k in K[2:end]],
		nc_E_ek[e,k] == nc_E_ek[e,k-1] + D_k[k]*nc_P_ek[e,k]
	)

	# define the dissatisfaction for each EV e at each timestep k
	nc_dissatisfaction_ek = [
		(Emax_e[e]-nc_E_ek[e,k])/Emax_e[e]
		for e in EVs, k in K]
	# the objective is to minimize the total dissatisfaction across EVs and timesteps,
	# taking into account the duration of each step
	JuMP.@objective(nc_model, Min,
		sum(
			D_k[k]*nc_dissatisfaction_ek[e,k]
			for e in EVs, k in K)
		)
end;

# ╔═╡ 496df139-97ec-4d2a-9c13-1ecd61ffa64f
md"""
**Establish link between PMD and EV charging**
"""

# ╔═╡ 565a33a7-3789-4a22-b6e1-0333cc5b7324
gen_name2ind = Dict(gen["name"] => gen["index"] for (_, gen) in data_math_mn["nw"]["1"]["gen"]);


# ╔═╡ ea6237f8-d9f8-4cc5-a2f9-86f045254a6c
ev_gen_ind_e = [gen_name2ind["EV_gen_$e"] for e in EVs];

# ╔═╡ 3d42555f-99df-41d7-8726-9f8e990aedc0
begin
	nc_Pg_ek = [var(pm, k, :pg, ev_gen_ind_e[e]).data[1] for e in EVs, k in K]
	nc_Qg_ek = [var(pm, k, :qg, ev_gen_ind_e[e]).data[1] for e in EVs, k in K]
end;

# ╔═╡ f7b500a6-7dbf-44cb-9433-3b534f13cd6b
begin
	# link charge to generator models
	JuMP.@constraint(nc_model, [e in EVs, k in K],
		nc_Pg_ek[e,k] == -nc_P_ek[e,k]
		)
	JuMP.@constraint(nc_model, [e in EVs, k in K],
		nc_Qg_ek[e,k] == 0.0
		)
end;

# ╔═╡ 284972a8-0e9f-495d-b778-812d04febba1
md"""
**Solve**

As before, we could solve the JuMP model directly, i.e. `JuMP.optimize!(nc_model)`. However, it is better to use the PMD wrapper, because this will also generate a solution dictionary for the network variables.
"""

# ╔═╡ d07996ac-67de-41bd-a8b5-35c8e513b145
begin
	res = optimize_model!(pm, optimizer=Ipopt.Optimizer)
	res["termination_status"]
end

# ╔═╡ ef198d40-77c8-4aa4-8af7-f3a535c7c74b
md"""
**Inspect EV charging variables**

Finally, let's explore the new solution through a series of figures.

Below you can see how the SoC evolves for the EVs with the new charging schedule. By the end, all EVs are now fully charged.
"""

# ╔═╡ f39c50ef-2251-4249-b41b-0a2c87fe18e1
begin
	nc_E_ek_vals = JuMP.value.(nc_E_ek.data)
	Plots.plot(legend=:none, title="", xlabel="time [h]", ylabel="SoC[kWh]", ylim=[0, maximum(Emax_e)])
	for e in EVs
		Plots.plot!(timestamps, [E0_e[e], nc_E_ek_vals[e,:]...], markershape=:circle, markersize=3)
	end
	Plots.plot!()
end

# ╔═╡ 4cf28814-0a1c-4df6-9ab3-ecb4b89cd66f
md"""
Below you see the optimal charging rate setpoints for all EVs. These regularly exceed the conservative 5 kW limit (indicated with a gray line).
"""

# ╔═╡ 0d75f5db-16c0-4345-86ba-b72491ee82d1
begin
	nc_P_ek_vals = JuMP.value.(nc_P_ek.data)
	Plots.plot(legend=:none, title="", xlabel="time [h]", ylabel="charge [kW]")
	Plots.plot!([timestamps[1], timestamps[end-1]], [Pmax, Pmax], color=:gray, linewidth=2)
	for e in EVs
		Plots.scatter!(timestamps[1:end-1], nc_P_ek_vals[e,:], markershape=:circle, markersize=3)
	end
	Plots.plot!()
end

# ╔═╡ 37ab2c49-7016-417c-9715-a13a78ed50b8
md"""
**Inspect network variables**

First, we convert the `MATHEMATICAL` solution back to the `ENGINEERING` representation. This will allow us to inspect the results with the familiar component names and in SI-units. The `MATHEMATICAL` solution uses indices instead of string identifiers, and is in pu.
"""

# ╔═╡ 8795bf04-cdcb-445c-8916-f5ac6497cb79
begin
	sol_math = res["solution"]
	sol_eng  = transform_solution(sol_math, data_math_mn)
end;

# ╔═╡ 19d23af9-27e6-4a54-bdfd-593835de0861
md"""
There are a total of 117x3 phase voltages. To keep the plot light, we extract below the voltage for only those terminals which have a load connected to them. Furthermore, we convert th valus again to pu, because that is how we specified the bounds as well.
"""

# ╔═╡ 4fdd43cd-dc38-4e42-b003-dfb9f2c87739
begin
	vm_pu_lk = fill(NaN, length(data_eng["load"]), length(K))
	for k in K, l in 1:length(data_eng["load"])
		bus_id = data_eng["load"]["load$l"]["bus"]
		bus_ind = data_math_mn["bus_lookup"][bus_id]
		sol_bus = sol_eng["nw"]["$k"]["bus"][bus_id]
		data_bus = data_eng["bus"][bus_id]
		vbase = data_math_mn["nw"]["$k"]["bus"]["$bus_ind"]["vbase"]
		phase = data_eng["load"]["load$l"]["connections"][1]
		ind = findfirst(data_bus["terminals"].==phase)
		vm_pu_lk[l,k] = abs(sol_bus["vr"][ind]+im*sol_bus["vi"][ind])/vbase
	end
end

# ╔═╡ da7d6303-490a-4073-9b10-440cb9e23097
md"""
Below you can find scatter plot of the load phase voltages at each timestep. The voltage bounds are indicated by the red lines.
"""

# ╔═╡ 4b41e9dc-5c29-4979-a43a-714b1eeb804d
begin
	Plots.plot(xlabel="time [h]", ylabel="load phase voltage [pu]", legend=:none)
	Plots.plot!([timestamps[K[1]], timestamps[K[end]]], [0.9, 0.9], color=:red, linewidth=3)
	Plots.plot!([timestamps[K[1]], timestamps[K[end]]], [1.1, 1.1], color=:red, linewidth=3)
	for k in K
		Plots.scatter!(fill(timestamps[k], length(data_eng["load"])), vm_pu_lk[:,k], markershape=:circle, markersize=3, label="")
	end
	Plots.
	Plots.plot!()
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
InfrastructureModels = "2030c09a-7f63-5d83-885d-db604e0e9cc0"
Ipopt = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PowerModelsDistribution = "d7431456-977f-11e9-2de3-97ff7677985e"

[compat]
InfrastructureModels = "~0.7.6"
Ipopt = "~1.2.0"
JuMP = "~1.9.0"
Plots = "~1.38.7"
PlutoUI = "~0.7.50"
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

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

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

[[BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

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

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

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

[[ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

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

[[Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

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

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

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

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

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

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

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

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "660b2ea2ec2b010bb02823c6d0ff6afd9bdc5c16"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.71.7"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d5e1fd17ac7f3aa4c5287a61ee28d4f8b8e98873"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.71.7+0"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[Glob]]
git-tree-sha1 = "4df9f7e06108728ebf00a0a11edee4b29a482bb2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.0"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "37e4657cd56b11abe3d10cd4a1ec5fbdb4180263"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.7.4"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[InfrastructureModels]]
deps = ["JuMP", "Memento"]
git-tree-sha1 = "88da90ad5d8ca541350c156bea2715f3a23836ce"
uuid = "2030c09a-7f63-5d83-885d-db604e0e9cc0"
version = "0.7.6"

[[IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

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

[[JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

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

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[JuMP]]
deps = ["LinearAlgebra", "MathOptInterface", "MutableArithmetics", "OrderedCollections", "Printf", "SnoopPrecompile", "SparseArrays"]
git-tree-sha1 = "611b9f12f02c587d860c813743e6cec6264e94d8"
uuid = "4076af6c-e467-56ae-b986-b466b2749572"
version = "1.9.0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

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

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

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

[[MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

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

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

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

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

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

[[OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "6503b77492fd7fcb9379bf73cd31035670e3c509"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.3"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "c95373e73290cf50a8a22c3375e4625ded5c5280"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.4"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "cfcd24ebf8b066b4f8e42bade600c8558212ed83"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.7"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

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

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "e974477be88cb5e3040009f3767611bc6357846f"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.11"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "77d3c4726515dca71f6d80fbb5e251088defe305"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.18"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

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

[[StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

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

[[TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c6edfe154ad7b313c01aceca188c05c835c67360"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.4+0"

[[fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─e12e9ffa-1d5e-11ec-0ef4-435a3b4ee5d5
# ╟─947fccff-8647-48c5-9bd1-c483ee589a33
# ╠═a3d1557a-000b-4bc8-a390-571c08e55586
# ╠═f1abaf46-59bc-4ad3-8054-0509e37c959c
# ╟─6e7eae9b-de28-4069-bdad-a8d4033c533c
# ╟─1842fe02-c874-42c7-926e-d843868eaa5f
# ╟─550efc60-e19e-46d3-bce3-090b17ccbb30
# ╟─7569cb6f-1ad2-465b-b0cf-9edccb0e617a
# ╠═6f96134f-fa1e-47e9-8047-9814e5b911ec
# ╟─f2c56f02-f6c3-4576-ac6d-81dc6e7b20b3
# ╠═da345583-51b6-4537-a540-e1100ff90d06
# ╟─38924354-4977-46e8-80d8-256d11a10242
# ╠═78ba1c9e-bb56-4158-a67b-e77770f73645
# ╟─bb8be0dc-d8b3-4026-a9ea-1b7d1b04b5df
# ╟─bfaf5d7d-1042-45fd-bd22-0d2a78bf4d78
# ╟─7847a75f-28f2-472f-9f0a-0d40e7a40cf9
# ╠═8a7de573-2e02-4332-b416-0e08dcb8c6d1
# ╟─72c124c1-fc0b-40f1-a7c9-7e475c3f74d7
# ╠═e6ac4be0-e524-4341-a1ec-28e16a3b2bfb
# ╟─cd654189-90fe-49e2-a0b9-c9f4e3c3f97c
# ╠═6ff8cb3b-d153-4639-8eda-3a298860a905
# ╟─8e0fff25-e6fd-439b-96ec-e5a13c9cff24
# ╠═5c3b70e4-99d0-4dd0-a7ef-d62bfa904bb7
# ╟─eb8479f7-99c6-49c9-b767-43a8d31cc502
# ╠═c40088d6-78b7-44dc-ac4f-ab5357490818
# ╟─09ae6158-3bd9-4d95-8e1e-029f308c7a91
# ╠═fb7f7904-a219-46a5-ae0b-f728dc4a3fd9
# ╟─5988c7bd-6419-4d46-9a0c-6c83f5bc4c25
# ╟─49539850-e2b1-4976-95dc-1a39a160747a
# ╟─ce8ac568-c3b5-4413-936f-b42ae2b9234d
# ╟─2b2d5e60-f4ce-4813-9608-5090bf6c5213
# ╟─423c6e71-9978-4845-89e4-803f32098af4
# ╠═0be0abd5-f0ad-4d15-b6fb-d59d74e6058c
# ╟─04ef243c-2586-406e-9528-140d732c28f7
# ╠═51d2dae6-a59c-4014-bdd6-886864847683
# ╟─2bd488f2-1a73-4c8e-a2bd-42a6fbd76d94
# ╠═1d878788-13bd-4090-9ed0-b9eb77a8575d
# ╟─548ffacb-af94-46c1-bcea-3695f98b4516
# ╠═0f97d7aa-cdbe-454c-83f0-978964c83b2a
# ╟─161c2052-7fbd-4993-a21f-3ed14b750619
# ╠═8d24de5c-749c-4b0c-acaa-25d010a33844
# ╟─eb8c0d08-5ac0-446f-9447-944399eacdb0
# ╠═8cf4316e-1cfa-4cf0-abff-1438edcb6d36
# ╟─f4344910-d95e-4f33-9b16-3f52c9ced4ac
# ╠═97017a8b-f86b-49b4-ac83-66303df1f63c
# ╟─3a1abed7-d111-49e1-bcfb-ee9af48d4a6d
# ╠═52c40603-42fe-45d3-840c-f530fc3951f2
# ╟─b2e5c25e-e497-4080-ba74-2cfa8e6c05d4
# ╠═49e89e0f-480d-4dc1-850a-4ad9c4ae0ebd
# ╟─8ac591c5-e3d4-42bd-994c-7bc885536824
# ╠═11ba8af0-50d3-4f5d-8c2b-4a1af9e5f5d5
# ╟─5ecc3cc4-0f8e-4ec5-8e2f-47a5871f1304
# ╠═bbba2730-31a7-4c27-ab88-c05f358b99b6
# ╟─ff40f2b8-1f77-4e2b-a0f8-27415a8e3978
# ╠═12d8fc45-a242-4c53-ac31-161af6331f0a
# ╟─496df139-97ec-4d2a-9c13-1ecd61ffa64f
# ╠═565a33a7-3789-4a22-b6e1-0333cc5b7324
# ╠═ea6237f8-d9f8-4cc5-a2f9-86f045254a6c
# ╠═3d42555f-99df-41d7-8726-9f8e990aedc0
# ╠═f7b500a6-7dbf-44cb-9433-3b534f13cd6b
# ╟─284972a8-0e9f-495d-b778-812d04febba1
# ╠═d07996ac-67de-41bd-a8b5-35c8e513b145
# ╟─ef198d40-77c8-4aa4-8af7-f3a535c7c74b
# ╟─f39c50ef-2251-4249-b41b-0a2c87fe18e1
# ╟─4cf28814-0a1c-4df6-9ab3-ecb4b89cd66f
# ╟─0d75f5db-16c0-4345-86ba-b72491ee82d1
# ╟─37ab2c49-7016-417c-9715-a13a78ed50b8
# ╠═8795bf04-cdcb-445c-8916-f5ac6497cb79
# ╟─19d23af9-27e6-4a54-bdfd-593835de0861
# ╠═4fdd43cd-dc38-4e42-b003-dfb9f2c87739
# ╟─da7d6303-490a-4073-9b10-440cb9e23097
# ╟─4b41e9dc-5c29-4979-a43a-714b1eeb804d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
