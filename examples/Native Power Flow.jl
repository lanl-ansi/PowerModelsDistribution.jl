### A Pluto.jl notebook ###
# v0.19.16

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
	import Pkg
	Pkg.activate("./native_pf_testcases")
	# Pkg.activate(mktempdir())
	Pkg.instantiate()
	Pkg.add(url="https://github.com/hei06j/PowerModelsDistribution.jl.git#four-wire-native-pf")
	Pkg.add("PlutoUI")
	Pkg.add("JSON")
	
	Pkg.status()

	using PowerModelsDistribution, PlutoUI, JSON
	const PMD = PowerModelsDistribution	
end

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

### Voltage source

Voltage source impedance matrix size should be compatible with whether or not the neutral conductor is explicitly represented or not. 

### Source bus 

The vector sizes in source bus data should also be compatible with whether or not the neutral conductor is explicitly represented or not. 

### Other network elements

Vector and matrix sizes of other elements should also be compatible with whether or not the neutral conductor is explicitly represented or not. 
"""

# ╔═╡ b343dd2b-5624-4836-b12f-357b29025772
function vsource_correction!(data_eng; explicit_neutral=true)
    if explicit_neutral
        if haskey(data_eng, "multinetwork")
            for (n,nw) in data_eng["nw"]
                nw["voltage_source"]["source"]["rs"][4,4] = nw["voltage_source"]["source"]["rs"][1,1]
                nw["voltage_source"]["source"]["rs"][1:3,4] .= nw["voltage_source"]["source"]["rs"][1,2]
                nw["voltage_source"]["source"]["rs"][4,1:3] .= nw["voltage_source"]["source"]["rs"][1,2]
                nw["voltage_source"]["source"]["xs"][4,4] = nw["voltage_source"]["source"]["xs"][1,1]
                nw["voltage_source"]["source"]["xs"][1:3,4] .= nw["voltage_source"]["source"]["xs"][1,2]
                nw["voltage_source"]["source"]["xs"][4,1:3] .= nw["voltage_source"]["source"]["xs"][1,2]
            end
        else
            data_eng["voltage_source"]["source"]["rs"][4,4] = data_eng["voltage_source"]["source"]["rs"][1,1]
            data_eng["voltage_source"]["source"]["rs"][1:3,4] .= data_eng["voltage_source"]["source"]["rs"][1,2]
            data_eng["voltage_source"]["source"]["rs"][4,1:3] .= data_eng["voltage_source"]["source"]["rs"][1,2]
            data_eng["voltage_source"]["source"]["xs"][4,4] = data_eng["voltage_source"]["source"]["xs"][1,1]
            data_eng["voltage_source"]["source"]["xs"][1:3,4] .= data_eng["voltage_source"]["source"]["xs"][1,2]
            data_eng["voltage_source"]["source"]["xs"][4,1:3] .= data_eng["voltage_source"]["source"]["xs"][1,2]        
        end
    else
        if haskey(data_eng, "multinetwork")
            for (n,nw) in data_eng["nw"]
                nw["voltage_source"]["source"]["rs"] = nw["voltage_source"]["source"]["rs"][1:3,1:3]
                nw["voltage_source"]["source"]["xs"] = nw["voltage_source"]["source"]["xs"][1:3,1:3]
                nw["voltage_source"]["source"]["connections"] = nw["voltage_source"]["source"]["connections"][1:3]
            end
        else
            data_eng["voltage_source"]["source"]["rs"] = data_eng["voltage_source"]["source"]["rs"][1:3,1:3]
            data_eng["voltage_source"]["source"]["xs"] = data_eng["voltage_source"]["source"]["xs"][1:3,1:3]
            data_eng["voltage_source"]["source"]["connections"] = data_eng["voltage_source"]["source"]["connections"][1:3]
        end
    end
    return nothing
end

# ╔═╡ ce0d2924-94c6-484e-9d0a-c805e8a5d5df
function sourcebus_voltage_vector_correction!(data_math::Dict{String, Any}; explicit_neutral=true)
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
                        if neutral_idx ∈ bus["terminals"]
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
                    if neutral_idx ∈ bus["terminals"]
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
        if haskey(bus, "terminals") && neutral_idx ∈ bus["terminals"]
            explicit_neutral = true
        end

        if explicit_neutral
            idx = findall(x->x==neutral_idx, bus["terminals"])
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
        if haskey(branch, "t_connections") && neutral_idx ∈ branch["t_connections"]
            explicit_neutral = true
            deleteat!(branch["t_connections"], branch["t_connections"] .== neutral_idx)
            # branch["t_connections"] = branch["t_connections"][1:end-1]
        end
        if haskey(branch, "f_connections") && neutral_idx ∈ branch["f_connections"]
            explicit_neutral = true
            deleteat!(branch["f_connections"], branch["f_connections"] .== neutral_idx)
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
        if haskey(transformer, "t_connections") && neutral_idx ∈ transformer["t_connections"]
            if transformer["t_connections"][end] !== neutral_idx
                transformer["polarity"] = -1
                deleteat!(transformer["t_connections"], transformer["t_connections"] .== neutral_idx)
            else
                transformer["t_connections"] = transformer["t_connections"][1:end-1]
            end
        end
        if haskey(transformer, "f_connections") && neutral_idx ∈ transformer["f_connections"]
            if transformer["f_connections"][end] !== neutral_idx
                transformer["polarity"] = -1
                deleteat!(transformer["f_connections"], transformer["f_connections"] .== neutral_idx)
            else
                transformer["f_connections"] = transformer["f_connections"][1:end-1]
            end
        end
    end

    for (g,gen) in math["gen"]
        if neutral_idx in gen["connections"]
            idx = findall(x->x==neutral_idx, gen["connections"])
            deleteat!(gen["connections"], gen["connections"].==gen["connections"][idx])
            deleteat!(gen["vg"], gen["vg"].==gen["vg"][idx])
            deleteat!(gen["pg"], gen["pg"].==gen["pg"][idx])
            deleteat!(gen["qg"], gen["qg"].==gen["qg"][idx])
            deleteat!(gen["pmax"], gen["pmax"].==gen["pmax"][idx])
            deleteat!(gen["pmin"], gen["pmin"].==gen["pmin"][idx])
            deleteat!(gen["qmax"], gen["qmax"].==gen["qmax"][idx])
            deleteat!(gen["qmin"], gen["qmin"].==gen["qmin"][idx])
            # gen["connections"] = gen["connections"][1:end-1]
            # gen["vg"] = gen["vg"][1:end-1]
            # gen["pg"] = gen["pg"][1:end-1]
            # gen["qg"] = gen["qg"][1:end-1]
            # gen["pmax"] = gen["pmax"][1:end-1]
            # gen["pmin"] = gen["pmin"][1:end-1]
            # gen["qmax"] = gen["qmax"][1:end-1]
            # gen["qmin"] = gen["qmin"][1:end-1]
            gen["cost"] = 1000 .* gen["cost"]
        end
    end

    for (l,load) in math["load"]
        if load["configuration"] == WYE && neutral_idx ∈ load["connections"]
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

# ╔═╡ c35141d1-df28-4260-aef6-a5a7e50149c5
begin
	function solve_compute_mc_pf(dss_file, solution_file; explicit_neutral=true, max_iter=100)
	    if explicit_neutral
	        data_eng = PMD.parse_file(dss_file, transformations=[transform_loops!])
	        vsource_correction!(data_eng, explicit_neutral=true)
	
	        data_math = PMD.transform_data_model(data_eng;kron_reduce=false)
	        res = PMD.compute_mc_pf(data_math; explicit_neutral=true, max_iter=max_iter)
	    else
	        data_eng = PMD.parse_file(dss_file, transformations=[transform_loops!]);
	        data_eng["is_kron_reduced"] = true
	        data_eng["settings"]["sbase_default"] = 1
	        vsource_correction!(data_eng, explicit_neutral=false)
	
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
    pmd_path = joinpath(dirname(pathof(PowerModelsDistribution)), "..")
	@bind case PlutoUI.Select([
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
            "test_switch" => "three wire network with switch",
			"test_switch_3w" => "three wire network with switch",
			"test_switch_1w" => "single wire network with switch",
			"ut_trans_3w_yyy_1" => "three wire network with three winding (yyy) transformer",
			"ut_trans_3w_dyy_1" => "three wire network with three winding (dyy) transformer",
			"test_trans_dy_3w" => "three wire network with two winding (dy) transformer",
			"test_trans_yy_3w" => "three wire network with two winding (yy) transformer",
            "test_trans_dy" => "four wire network with two winding (dy) transformer",
			"test_trans_yy" => "four wire network with two winding (dy) transformer",
		])
    
    case_path = joinpath(pmd_path, "test/data/en_validation_case_data/$case.dss")
    solution_path = joinpath(pmd_path, "test/data/en_validation_case_solutions/$case.json")

    if case ∈ ["test_switch_3w", "test_switch_1w", "ut_trans_3w_yyy_1", "ut_trans_3w_dyy_1", "test_trans_dy_3w", "test_trans_yy_3w"]
        explicit_neutral = false
    else
        explicit_neutral = true
    end

    data_eng, data_math, res, v_maxerr_pu = solve_compute_mc_pf(case_path, solution_path; explicit_neutral=explicit_neutral);
    
    @show v_maxerr_pu
end



# ╔═╡ b5d99973-ea63-41dc-90c8-a744104437fe
md"""
### Larger networks

The native power flow solver is also tested on larger networks 

- IEEE 13, IEEE 34, IEEE 123 (https://github.com/sanderclaeys/DistributionTestCases.jl)

- Egrid GreensBoro and SantaFe networks (https://egriddata.org/dataset/greensboro-synthetic-network, https://egriddata.org/dataset/santa-fe-synthetic-network)

The table below shows the maximum per unit voltage error for the tested networks:

| `Network`                                     | `maximum voltage error pu`|
| --------------------------------------------  | ------------------------  |
| "IEEE 13 three wire network"                  |  3.765097382122667e-6     |
| "IEEE 34 three wire network"                  |  6.801818003195945e-8     |
| "IEEE 123 three wire network"                 |  4.044977697179336e-8     |
| "Egrid GreensBoro Industrial"                 |  0.0018373325064614753    |
| "Egrid SantaFe Urban/Suberban uhs0 1"         |  0.00011703535416486209   |
"""

# ╔═╡ Cell order:
# ╠═3da55c61-2439-4282-91ff-af8c794872e4
# ╟─e12e9ffa-1d5e-11ec-0ef4-435a3b4ee5d5
# ╟─95a94c06-3652-4fde-b496-087fb2d25948
# ╟─9cfb3edd-af18-4cf7-8310-0bfeb3d81826
# ╟─b343dd2b-5624-4836-b12f-357b29025772
# ╟─ce0d2924-94c6-484e-9d0a-c805e8a5d5df
# ╟─8fa26a09-0352-45e8-bc31-c5e442a56cc2
# ╟─179a598f-2efc-4e7f-9287-1dbc9a68cf4f
# ╟─c35141d1-df28-4260-aef6-a5a7e50149c5
# ╟─987398ad-4a31-47f7-af20-4f2efe7eb36f
# ╟─b5d99973-ea63-41dc-90c8-a744104437fe
