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
			Pkg.PackageSpec(;name="PowerModelsDistribution", version="0.12", rev="main"),
			Pkg.PackageSpec(;name="Ipopt", version="0.6.5"),
			])
end


# ╔═╡ a3d1557a-000b-4bc8-a390-571c08e55586
begin
	import JuMP
	import PlutoUI
	using PowerModelsDistribution; const PMD = PowerModelsDistribution
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

Before any of this, we will first setup the Julia environment for this tutorial.
"""

# ╔═╡ 52144adc-3c7f-408b-a53b-3ef14c6ea135
md"""

## Julia environment setup
The following code block will setup a Julia environment for you with the correct versions of packages for this Pluto notebook. It is possible that more recent versions will work as well, but the notebook has only been tested with this combination of package versions.
"""

# ╔═╡ 947fccff-8647-48c5-9bd1-c483ee589a33
md"""
This notebook will make use of the following packages in various places. Note that we do not apply `using` to any packages; we avoid this so you do not get confused where certain methods are coming from. Also, we introduce the shorthand `PMD` for `PowerModelsDistribution`.
"""

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
	"$NOTEBOOK_DIR/data/lvtestcase_notrans.dss", 
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
	JuMP.@NLconstraint(nc_model, [e in EVs, k in K], 
		nc_Pg_ek[e,k] == -nc_P_ek[e,k]
		)
	JuMP.@NLconstraint(nc_model, [e in EVs, k in K], 
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

# ╔═╡ Cell order:
# ╟─e12e9ffa-1d5e-11ec-0ef4-435a3b4ee5d5
# ╟─52144adc-3c7f-408b-a53b-3ef14c6ea135
# ╠═22cd4c84-7d9d-43dd-9b5c-519016b4ff41
# ╟─947fccff-8647-48c5-9bd1-c483ee589a33
# ╠═a3d1557a-000b-4bc8-a390-571c08e55586
# ╠═6e7eae9b-de28-4069-bdad-a8d4033c533c
# ╠═1842fe02-c874-42c7-926e-d843868eaa5f
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
# ╠═3a1abed7-d111-49e1-bcfb-ee9af48d4a6d
# ╠═52c40603-42fe-45d3-840c-f530fc3951f2
# ╟─b2e5c25e-e497-4080-ba74-2cfa8e6c05d4
# ╠═49e89e0f-480d-4dc1-850a-4ad9c4ae0ebd
# ╟─8ac591c5-e3d4-42bd-994c-7bc885536824
# ╠═11ba8af0-50d3-4f5d-8c2b-4a1af9e5f5d5
# ╟─5ecc3cc4-0f8e-4ec5-8e2f-47a5871f1304
# ╠═bbba2730-31a7-4c27-ab88-c05f358b99b6
# ╠═ff40f2b8-1f77-4e2b-a0f8-27415a8e3978
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
