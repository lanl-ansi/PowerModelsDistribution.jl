### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 3da55c61-2439-4282-91ff-af8c794872e4
begin
    using PlutoUI
    import JSON
    import PowerModelsDistribution as PMD
end

# ╔═╡ 1faefc14-e9f3-4728-98b8-55b75b79b350
md"""
### Setting up the notebook environment

In order to install a specific package branch or a commit, the notebook environment needs to be activated before Pluto.run(), as follows:

 - `julia> import Pluto`

 - `julia> Pluto.activate_notebook_environment("~/notebook.jl")`

 - `]`

 - `(notebook.jl) > add specific_package_branch_commit_url`

 - `Pluto.run()`

"""

# ╔═╡ 107ba92b-3042-4121-81ea-f5b5697bb950
md"""
This notebook uses the following packages.
"""

# ╔═╡ e12e9ffa-1d5e-11ec-0ef4-435a3b4ee5d5
md"""
# Native Power Flow Solver

PowerModelsDistribution.jl (PMD) uses a fixed point iteration current injection method to solve power flow, similar to OpenDSS built-in solver.

This notebook illustrates accuracy of this embedded native power flow solver across diverse networks with and without an explicit neutral (EN) conductor representation. It will go through the full workflow, consisting of
- importing OpenDSS network data (and applying transformations as needed);
- updating network data;
- optimizing the model and inspecting the results;

"""


# ╔═╡ 95a94c06-3652-4fde-b496-087fb2d25948
md"""
## Importing Network Data and Inspecting Power Flow Results

We illustrate the native power flow accuracy by choosing from the list of networks below, each with different network elements, number of conductors, and configurations:
"""

# ╔═╡ 9cfb3edd-af18-4cf7-8310-0bfeb3d81826
md"""
## Updating network data

The network data need to be updated based on representation of the the explicit neitral conductor.

### Source bus

The vector sizes in source bus data should be compatible with whether or not the neutral conductor is explicitly represented or not.

### Other network elements

Vector and matrix sizes of other elements should also be compatible with whether or not the neutral conductor is explicitly represented or not.
"""


# ╔═╡ ce0d2924-94c6-484e-9d0a-c805e8a5d5df
function sourcebus_voltage_vector_correction!(data_math::Dict{String,<:Any}; explicit_neutral::Bool=true)
    if haskey(data_math, "multinetwork")
        for (n,nw) in data_math["nw"]
            for (i, bus) in data_math["nw"]["bus"]
                if bus["bus_type"] == 3
                    if explicit_neutral
                        bus["vm"] = bus["vm"][1:length(bus["terminals"])]
                        bus["va"] = bus["va"][1:length(bus["terminals"])]
                        bus["vmin"] = bus["vmin"][1:length(bus["terminals"])]
                        bus["vmax"] = bus["vmax"][1:length(bus["terminals"])]
                        bus["grounded"] = bus["grounded"][1:length(bus["terminals"])]
                    else
                        if PMD.neutral_idx ∈ bus["terminals"]
                            bus["terminals"] = bus["terminals"][1:end-1]
                        end
                        bus["vm"] = bus["vm"][1:length(bus["terminals"])]
                        bus["va"] = bus["va"][1:length(bus["terminals"])]
                        bus["vmin"] = bus["vmin"][1:length(bus["terminals"])]
                        bus["vmax"] = bus["vmax"][1:length(bus["terminals"])]
                        bus["grounded"] = bus["grounded"][1:length(bus["terminals"])]
                    end
                end
            end
        end
    else
        for (i, bus) in data_math["bus"]
            if bus["bus_type"] == 3
                if explicit_neutral
                    bus["vm"] = bus["vm"][1:length(bus["terminals"])]
                    bus["va"] = bus["va"][1:length(bus["terminals"])]
                    bus["vmin"] = bus["vmin"][1:length(bus["terminals"])]
                    bus["vmax"] = bus["vmax"][1:length(bus["terminals"])]
                    bus["grounded"] = bus["grounded"][1:length(bus["terminals"])]
                else
                    if PMD.neutral_idx ∈ bus["terminals"]
                        bus["terminals"] = bus["terminals"][1:end-1]
                    end
                    bus["vm"] = bus["vm"][1:length(bus["terminals"])]
                    bus["va"] = bus["va"][1:length(bus["terminals"])]
                    bus["vmin"] = bus["vmin"][1:length(bus["terminals"])]
                    bus["vmax"] = bus["vmax"][1:length(bus["terminals"])]
                    bus["grounded"] = bus["grounded"][1:length(bus["terminals"])]
                end
            end
        end
    end
    return nothing
end


# ╔═╡ 8fa26a09-0352-45e8-bc31-c5e442a56cc2
function update_math_model_3wire!(math)
    math["conductor_ids"] = math["conductor_ids"][1:3]

    for (i,bus) in math["bus"]
        explicit_neutral = false
        if haskey(bus, "terminals") && PMD.neutral_idx ∈ bus["terminals"]
            explicit_neutral = true
        end

        if explicit_neutral
            idx = findall(x->x==PMD.neutral_idx, bus["terminals"])
            if haskey(bus, "terminals")
                deleteat!(bus["terminals"], bus["terminals"].==bus["terminals"][idx])
                # bus["terminals"] = bus["terminals"][1:end-1]
            end
            if haskey(bus, "grounded")
                deleteat!(bus["grounded"], bus["grounded"].==bus["grounded"][idx])
                # bus["grounded"] = bus["grounded"][1:end-1]
            end
            if haskey(bus, "vmax")
                deleteat!(bus["vmax"], bus["vmax"].==bus["vmax"][idx])
                # bus["vmax"] = bus["vmax"][1:end-1]
            end
            if haskey(bus, "vmin")
                deleteat!(bus["vmin"], bus["vmin"].==bus["vmin"][idx])
                # bus["vmin"] = bus["vmin"][1:end-1]
            end
            bus["vmin"] = 0.9 * ones(length(bus["terminals"]))
            bus["vmax"] = 1.1 * ones(length(bus["terminals"]))
        end
    end

    for (l,branch) in math["branch"]
        explicit_neutral = false
        if haskey(branch, "t_connections") && PMD.neutral_idx ∈ branch["t_connections"]
            explicit_neutral = true
            deleteat!(branch["t_connections"], branch["t_connections"] .== PMD.neutral_idx)
            # branch["t_connections"] = branch["t_connections"][1:end-1]
        end
        if haskey(branch, "f_connections") && PMD.neutral_idx ∈ branch["f_connections"]
            explicit_neutral = true
            deleteat!(branch["f_connections"], branch["f_connections"] .== PMD.neutral_idx)
            # branch["f_connections"] = branch["f_connections"][1:end-1]
        end
        if haskey(branch, "br_r") && explicit_neutral
            branch["br_r"] = branch["br_r"][1:end-1,1:end-1]
        end
        if haskey(branch, "br_x") && explicit_neutral
            branch["br_x"] = branch["br_x"][1:end-1,1:end-1]
        end
        if haskey(branch, "g_to") && explicit_neutral
            branch["g_to"] = branch["g_to"][1:end-1,1:end-1]
        end
        if haskey(branch, "g_fr") && explicit_neutral
            branch["g_fr"] = branch["g_fr"][1:end-1,1:end-1]
        end
        if haskey(branch, "b_to") && explicit_neutral
            branch["b_to"] = branch["b_to"][1:end-1,1:end-1]
        end
        if haskey(branch, "b_fr") && explicit_neutral
            branch["b_fr"] = branch["b_fr"][1:end-1,1:end-1]
        end
        if haskey(branch, "c_rating_a") && explicit_neutral
            branch["c_rating_a"] = branch["c_rating_a"][1:end-1]
        end
    end

    for (t,transformer) in math["transformer"]
        if haskey(transformer, "t_connections") && PMD.neutral_idx ∈ transformer["t_connections"]
            if transformer["t_connections"][end] !== PMD.neutral_idx
                transformer["polarity"] = -1
                deleteat!(transformer["t_connections"], transformer["t_connections"] .== PMD.neutral_idx)
            else
                transformer["t_connections"] = transformer["t_connections"][1:end-1]
            end
        end
        if haskey(transformer, "f_connections") && PMD.neutral_idx ∈ transformer["f_connections"]
            if transformer["f_connections"][end] !== PMD.neutral_idx
                transformer["polarity"] = -1
                deleteat!(transformer["f_connections"], transformer["f_connections"] .== PMD.neutral_idx)
            else
                transformer["f_connections"] = transformer["f_connections"][1:end-1]
            end
        end
    end

    for (g,gen) in math["gen"]
        if PMD.neutral_idx in gen["connections"]
            gen["connections"] = gen["connections"][1:end-1]
            gen["vg"] = gen["vg"][1:end-1]
            gen["pg"] = gen["pg"][1:end-1]
            gen["qg"] = gen["qg"][1:end-1]
            gen["pmax"] = gen["pmax"][1:end-1]
            gen["pmin"] = gen["pmin"][1:end-1]
            gen["qmax"] = gen["qmax"][1:end-1]
            gen["qmin"] = gen["qmin"][1:end-1]
            gen["cost"] = 1000 .* gen["cost"]
        end
    end

    for (l,load) in math["load"]
        if load["configuration"] == PMD.WYE && PMD.neutral_idx ∈ load["connections"]
            load["connections"] = load["connections"][1:end-1]
        end
    end

    return nothing
end


# ╔═╡ 179a598f-2efc-4e7f-9287-1dbc9a68cf4f
md"""
## Optimizing the Model and Inspecting the Results

The updated network model is then solved and the results are compared with OpenDSS solutions to illustrate the native power solver accuracy.
"""

# ╔═╡ 0aedaa8d-ab6e-4255-8535-ad6ae3ab4ecd
function compare_sol_dss_pmd(sol_dss::Dict{String,<:Any}, sol_pmd::Dict{String,<:Any}, data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}; compare_math::Bool=false, verbose::Bool=true, floating_buses::Vector=[], skip_buses::Vector=[], v_err_print_tol::Real=1E-6)
    max_v_err_pu = 0.0

    # voltage base for ENGINEERING buses in [V]
    vbase = Dict(id=>data_math["bus"]["$ind"]["vbase"]*data_math["settings"]["voltage_scale_factor"] for (id,ind) in data_math["bus_lookup"])

    buses_intersected = intersect(keys(sol_dss["bus"]), keys(sol_pmd["bus"]))
    for id in setdiff(buses_intersected, skip_buses)
        pmd_bus = sol_pmd["bus"][id]
        dss_bus = sol_dss["bus"][id]

        terminals = data_eng["bus"][id]["terminals"]
        if compare_math
            ts = filter(x->haskey(dss_bus["vm"], "$x"), terminals)
            v_dss = [dss_bus["vm"]["$t"]*exp(im*dss_bus["va"]["$t"]) for t in ts]
            # convert to V instead of usual kV
            v_pmd = [pmd_bus["vm"][idx]*exp(im*deg2rad(pmd_bus["va"][idx]))*data_eng["settings"]["voltage_scale_factor"] for (idx,t) in enumerate(ts)]
        else
            ts = filter(x->haskey(dss_bus["vm"], x), terminals)
            v_dss = [dss_bus["vm"]["$t"]*exp(im*dss_bus["va"]["$t"]) for t in ts]
            # convert to V instead of usual kV
            v_pmd = [(pmd_bus["vr"][idx]+im*pmd_bus["vi"][idx])*data_eng["settings"]["voltage_scale_factor"] for (idx,t) in enumerate(ts)]
        end

        # convert to pu
        v_dss_pu = v_dss/vbase[id]
        v_pmd_pu = v_pmd/vbase[id]

        # convert to diffs if floating
        N = length(v_dss)
        if id in floating_buses && N>1
            v_dss_pu = [v_dss_pu[i]-v_dss_pu[j] for i in 1:N for j in i:N if i!=j]/vbase[id]
            v_pmd_pu = [v_pmd_pu[i]-v_pmd_pu[j] for i in 1:N for j in i:N if i!=j]/vbase[id]
            labels = ["$(ts[i])-$(ts[j])" for i in 1:N for j in i:N if i!=j]
        else
            labels = string.(ts)
        end

        for i in eachindex(v_pmd_pu)
            v_err_pu = abs.(v_dss_pu[i]-v_pmd_pu[i]); max_v_err_pu = max(max_v_err_pu, v_err_pu)

            if v_err_pu>v_err_print_tol && verbose
                println("terminal $id.$(labels[i])")
                println("\t |U| dss: $(abs(v_dss_pu[i]))")
                println("\t     pmd: $(abs(v_pmd_pu[i]))")
                println("\t  ∠U dss:  $(angle(v_dss_pu[i]))")
                println("\t     pmd:  $(angle(v_pmd_pu[i]))")
            end
        end
    end

    return max_v_err_pu
end


# ╔═╡ c35141d1-df28-4260-aef6-a5a7e50149c5
begin
	function solve_compute_mc_pf(dss_file, solution_file; explicit_neutral=true, max_iter=100)
	    if explicit_neutral
	        data_eng = PMD.parse_file(dss_file, transformations=[PMD.transform_loops!])

	        data_math = PMD.transform_data_model(data_eng;kron_reduce=false)
	        res = PMD.compute_mc_pf(data_math; explicit_neutral=true, max_iter=max_iter)
	    else
	        data_eng = PMD.parse_file(dss_file, transformations=[PMD.transform_loops!]);
	        data_eng["is_kron_reduced"] = true
	        data_eng["settings"]["sbase_default"] = 1

	        data_math = PMD.transform_data_model(data_eng;kron_reduce=false, phase_project=false);
	        sourcebus_voltage_vector_correction!(data_math, explicit_neutral=false);
	        update_math_model_3wire!(data_math);
	        res = PMD.compute_mc_pf(data_math; explicit_neutral=false, max_iter=max_iter)
	    end

	    # obtain solution from dss
	    sol_dss = open(solution_file, "r") do f
	        JSON.parse(f)
	    end
	    sol_pmd = PMD.transform_solution(res["solution"], data_math, make_si=true);

	    v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)

	    return data_eng, data_math, res, v_maxerr_pu
	end
end

# ╔═╡ 987398ad-4a31-47f7-af20-4f2efe7eb36f
begin
	@bind case Select([
			"test_gen_1ph_delta" => "single phase EN network with delta generator",
			"test_gen_1ph_wye" => "single phase EN network with wye generator",
			"test_gen_3ph_delta" => "three phase EN network with delta generator",
			"test_gen_3ph_wye" => "three phase EN network with wye generator",
			"test_load_1ph_delta_cp" => "single phase EN network with delta constant power load",
			"test_load_1ph_wye_cp" => "single phase EN network with wye constant power load",
			"test_load_3ph_delta_cz" => "three phase EN network with delta constant impedance load",
			"test_load_3ph_delta_ci" => "three phase EN network with delta constant current load",
            "test_load_3ph_delta_cp" => "three phase EN network with delta constant power load",
			"test_load_3ph_wye_cz" => "three phase EN network with wye constant impedance load",
            "test_load_3ph_wye_ci" => "three phase EN network with wye constant current load",
            "test_load_3ph_wye_cp" => "three phase EN network with wye constant power load",
            "test_switch" => "four wire network with switch",
			"test_switch_3w" => "three wire network with switch",
			"test_switch_1w" => "single wire network with switch",
			"ut_trans_3w_yyy_1" => "three wire network with three winding (yyy) transformer",
			"ut_trans_3w_dyy_1" => "three wire network with three winding (dyy) transformer",
			"test_trans_dy_3w" => "three wire network with two winding (dy) transformer",
			"test_trans_yy_3w" => "three wire network with two winding (yy) transformer",
            "test_trans_dy" => "four wire network with two winding (dy) transformer",
			"test_trans_yy" => "four wire network with two winding (yy) transformer",
		])
end


# ╔═╡ b322c8f7-b761-44bf-ae18-afe8842c7b78
begin
	pmd_path = joinpath(dirname(pathof(PMD)), "..");
	case_path = joinpath(pmd_path, "test/data/en_validation_case_data/$case.dss");
    solution_path = joinpath(pmd_path, "test/data/en_validation_case_solutions");
	solution1 = joinpath(solution_path, "$case.json");

    if case ∈ ["test_switch_3w", "test_switch_1w", "ut_trans_3w_yyy_1", "ut_trans_3w_dyy_1", "test_trans_dy_3w", "test_trans_yy_3w"]
        explicit_neutral = false
    else
        explicit_neutral = true
    end

    data_eng, data_math, res, v_maxerr_pu = solve_compute_mc_pf(case_path, solution1; explicit_neutral=explicit_neutral);

	"maximum voltage error p.u is $v_maxerr_pu, total time is $(res["time_total"])"

end

# ╔═╡ b5d99973-ea63-41dc-90c8-a744104437fe
md"""
### Larger networks

The native power flow solver is tested on larger networks. and the results are reported in table below.

To interactively test the larger networks, the reader is encouraged to download the network OpenDSS files from the links provided and place them in the same folder as the notebook file, name the folder/file as network/network.dss, and update the case paths accordingly. The native power flow results are validated against the OpenDSS solution using OpenDSSDirect.jl.

- IEEE networks: (https://github.com/sanderclaeys/DistributionTestCases.jl)
  - IEEE 13  (13 buses, 13 branches and switches, 2 transformers)
  - IEEE 34  (34 buses, 51 branches and switches, 4 transformers)
  - IEEE 123 (123 buses, 126 branches and switches, 3 transformers)

- Egrid networks:
  - Egrid GreensBoro Industrial network (8460+ buses, 7750+ branches and switches, 1620+ transformers) (https://egriddata.org/dataset/greensboro-synthetic-network)
  - Egrid SantaFe uhs0_1247/uhs0_1247--udt4776 network (3280+ buses, 2860+ branches and switches, 485 transformers) (https://egriddata.org/dataset/santa-fe-synthetic-network)

The table below shows the maximum per unit voltage error for the tested networks:

| `Network`                                     | `maximum voltage error pu`|
| --------------------------------------------  | ------------------------  |
| "IEEE 13 three wire network"                  |  3.765097382122667e-6     |
| "IEEE 34 three wire network"                  |  6.801818003195945e-8     |
| "IEEE 123 three wire network"                 |  4.044977697179336e-8     |
| "Egrid GreensBoro Industrial network"         |  0.0018373325064614753    |
| "Egrid SantaFe/urban-suburban/uhs0_1247/uhs0_1247--udt4776 network" |  0.00011818425933292976   |
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PowerModelsDistribution = "d7431456-977f-11e9-2de3-97ff7677985e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "SnoopPrecompile", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "c700cce799b51c9045473de751e9319bdd1c6e94"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.9"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "485193efd2176b88e6622a39a246f8c5b600e74e"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.6"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "a4ad7ef19d2cdc2eff57abbbe68032b1cd0bd8f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.13.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Glob]]
git-tree-sha1 = "97285bbd5230dd766e9ef6749b80fc617126d496"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InfrastructureModels]]
deps = ["JuMP", "Memento"]
git-tree-sha1 = "dc1e2eba1a98aa457b629fe1723d9078ecb74340"
uuid = "2030c09a-7f63-5d83-885d-db604e0e9cc0"
version = "0.7.7"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuMP]]
deps = ["LinearAlgebra", "MathOptInterface", "MutableArithmetics", "OrderedCollections", "Printf", "SnoopPrecompile", "SparseArrays"]
git-tree-sha1 = "4ec0e68fecbbe1b78db2ddf1ac573963ed5adebc"
uuid = "4076af6c-e467-56ae-b986-b466b2749572"
version = "1.10.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "DataStructures", "ForwardDiff", "JSON", "LinearAlgebra", "MutableArithmetics", "NaNMath", "OrderedCollections", "PrecompileTools", "Printf", "SparseArrays", "SpecialFunctions", "Test", "Unicode"]
git-tree-sha1 = "8e054675d393ce5866dcdd6a071075e25e21a39c"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.15.1"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Memento]]
deps = ["Dates", "Distributed", "Requires", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "bb2e8f4d9f400f6e90d57b34860f6abdc51398e5"
uuid = "f28f55f0-a522-5efc-85c2-fe41dfb9b2d9"
version = "1.4.1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3295d296288ab1a0a2528feb424b854418acff57"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.2.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.PolyhedralRelaxations]]
deps = ["DataStructures", "ForwardDiff", "JuMP", "Logging", "LoggingExtras"]
git-tree-sha1 = "05f2adc696ae9a99be3de99dd8970d00a4dccefe"
uuid = "2e741578-48fa-11ea-2d62-b52c946f73a0"
version = "0.3.5"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PowerModelsDistribution]]
deps = ["CSV", "Dates", "FilePaths", "Glob", "InfrastructureModels", "JSON", "JuMP", "LinearAlgebra", "Logging", "LoggingExtras", "PolyhedralRelaxations", "SparseArrays", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "4cdf5575752da0238e5f8b1bcf382e31a79e2392"
repo-rev = "four-wire-native-pf"
repo-url = "../../../../../../../../Users/hei06j/Documents/repositories/remote/PowerModelsDistribution.jl"
uuid = "d7431456-977f-11e9-2de3-97ff7677985e"
version = "0.14.8"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "bc2bda41d798c2e66e7c44a11007bb329b15941b"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.0.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "77d3c4726515dca71f6d80fbb5e251088defe305"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.18"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "63e84b7fdf5021026d0f17f76af7c57772313d99"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.21"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "0b829474fed270a4b0ab07117dce9b9a2fa7581a"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.12"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─1faefc14-e9f3-4728-98b8-55b75b79b350
# ╟─107ba92b-3042-4121-81ea-f5b5697bb950
# ╠═3da55c61-2439-4282-91ff-af8c794872e4
# ╟─e12e9ffa-1d5e-11ec-0ef4-435a3b4ee5d5
# ╟─95a94c06-3652-4fde-b496-087fb2d25948
# ╟─9cfb3edd-af18-4cf7-8310-0bfeb3d81826
# ╠═ce0d2924-94c6-484e-9d0a-c805e8a5d5df
# ╠═8fa26a09-0352-45e8-bc31-c5e442a56cc2
# ╟─179a598f-2efc-4e7f-9287-1dbc9a68cf4f
# ╠═0aedaa8d-ab6e-4255-8535-ad6ae3ab4ecd
# ╠═c35141d1-df28-4260-aef6-a5a7e50149c5
# ╟─987398ad-4a31-47f7-af20-4f2efe7eb36f
# ╟─b322c8f7-b761-44bf-ae18-afe8842c7b78
# ╟─b5d99973-ea63-41dc-90c8-a744104437fe
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
