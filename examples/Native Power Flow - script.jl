using Pkg
using PowerModelsDistribution
const PMD = PowerModelsDistribution

using JSON

data_dir = "examples/native_pf_testcases"
solution_dir = "examples/native_pf_testcases/solutions"

###
function compare_sol_dss_pmd(sol_dss::Dict{String,Any}, sol_pmd::Dict{String,Any}, data_eng::Dict{String,Any}, data_math::Dict{String,Any}; compare_math=false, verbose=true, floating_buses=[], skip_buses=[], v_err_print_tol=1E-6)
    max_v_err_pu = 0.0

    # voltage base for ENGINEERING buses in [V]
    vbase = Dict(id=>data_math["bus"]["$ind"]["vbase"]*data_math["settings"]["voltage_scale_factor"] for (id,ind) in data_math["bus_lookup"])

    buses_intersected = intersect(keys(sol_dss["bus"]), keys(sol_pmd["bus"]))
    for id in setdiff(buses_intersected, skip_buses)
        pmd_bus = sol_pmd["bus"][id]
        dss_bus = sol_dss["bus"][id]

        terminals = data_eng["bus"][id]["terminals"]
        if compare_math
            ts = [t for t in string.(terminals) if haskey(dss_bus["vm"], t)]
            v_dss = [dss_bus["vm"][t]*exp(im*dss_bus["va"][t]) for t in ts]
            # convert to V instead of usual kV
            v_pmd = [pmd_bus["vm"][t]*exp(im*deg2rad(pmd_bus["va"][t]))*data_eng["settings"]["voltage_scale_factor"] for t in ts]
        else
            ts = [t for t in terminals if haskey(dss_bus["vm"], t)]
            v_dss = [dss_bus["vm"][t]*exp(im*dss_bus["va"][t]) for t in ts]
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

        for i in 1:length(v_pmd_pu)
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



function solve_compute_mc_pf(dss_file, solution_file; explicit_neutral=true, max_iter=100)
    if explicit_neutral
        data_eng = parse_file(dss_file, transformations=[transform_loops!])
        vsource_correction!(data_eng, explicit_neutral=true)

        data_math = transform_data_model(data_eng;kron_reduce=false)
        res = compute_mc_pf(data_math; explicit_neutral=true, max_iter=max_iter)
    else
        data_eng = parse_file(dss_file, transformations=[transform_loops!]);
        data_eng["is_kron_reduced"] = true
        data_eng["settings"]["sbase_default"] = 1
        vsource_correction!(data_eng, explicit_neutral=false)

        data_math = transform_data_model(data_eng;kron_reduce=false, phase_project=false);
        sourcebus_voltage_vector_correction!(data_math, explicit_neutral=false);
        update_math_model_3wire!(data_math);
        res = compute_mc_pf(data_math; explicit_neutral=false, max_iter=max_iter)
    end

    # obtain solution from dss
    sol_dss = open(solution_file, "r") do f
        JSON.parse(f)
    end
    sol_pmd = transform_solution(res["solution"], data_math, make_si=true);

    if occursin("ieee123_pmd", dss_file)
        delete!(sol_pmd["bus"], "87")
    end

    v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)

    return data_eng, data_math, res, v_maxerr_pu
end


## ############## Generators ##############
### 1-phase PV generator - Delta
case = "test_gen_1ph_delta"
data_eng, data_math, res, v_maxerr_pu_g1d = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);

### 1-phase PV generator - Wye
case = "test_gen_1ph_wye"
data_eng, data_math, res, v_maxerr_pu_g1y = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);

### 3-phase PV generator - Delta
case = "test_gen_3ph_delta"
data_eng, data_math, res, v_maxerr_pu_g3d = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);

### 3-phase PV generator - Wye
case = "test_gen_3ph_wye"
data_eng, data_math, res, v_maxerr_pu_g3y = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);

names_gen = ["test_gen_1ph_delta", "test_gen_1ph_wye", "test_gen_3ph_delta", "test_gen_3ph_wye"]
v_maxerr_pu_g = [v_maxerr_pu_g1d, v_maxerr_pu_g1y, v_maxerr_pu_g3d, v_maxerr_pu_g3y]
@show v_maxerr_pu_g


## ############## Loads ##############
### 1-phase load - Delta - Constant P
case = "test_load_1ph_delta_cp"
data_eng, data_math, res, v_maxerr_pu_l1dcp = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);

### 1-phase load - Wye - Constant P
case = "test_load_1ph_wye_cp"
data_eng, data_math, res, v_maxerr_pu_l1ycp = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);

### 3-phase load - Delta - Constant Z
case = "test_load_3ph_delta_cz"
data_eng, data_math, res, v_maxerr_pu_l3dcz = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);

### 3-phase load - Delta - Constant I
case = "test_load_3ph_delta_ci"
data_eng, data_math, res, v_maxerr_pu_l3dci = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);

### 3-phase load - Delta - Constant P
case = "test_load_3ph_delta_cp"
data_eng, data_math, res, v_maxerr_pu_l3dcp = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);

### 3-phase load - Wye - Constant Z
case = "test_load_3ph_wye_cz"
data_eng, data_math, res, v_maxerr_pu_l3ycz = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);

### 3-phase load - Wye - Constant I
case = "test_load_3ph_wye_ci"
data_eng, data_math, res, v_maxerr_pu_l3yci = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);

### 3-phase load - Wye - Constant P
case = "test_load_3ph_wye_cp"
data_eng, data_math, res, v_maxerr_pu_l3ycp = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);


names_load = ["test_load_1ph_delta_cp", "test_load_1ph_wye_cp", "test_load_3ph_delta_cz", "test_load_3ph_delta_ci", "test_load_3ph_delta_cp", "test_load_3ph_wye_cz", "test_load_3ph_wye_ci", "test_load_3ph_wye_cp"]
v_maxerr_pu_l = [v_maxerr_pu_l1dcp, v_maxerr_pu_l1ycp, v_maxerr_pu_l3dcz, v_maxerr_pu_l3dci, v_maxerr_pu_l3dcp, v_maxerr_pu_l3ycz, v_maxerr_pu_l3yci, v_maxerr_pu_l3ycp]
@show v_maxerr_pu_l


## ############## SWITCH ##############
### 4wire - switch
case = "test_switch"
data_eng, data_math, res, v_maxerr_pu_s4w = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);

### 3wire - switch
case = "test_switch_3w"
data_eng, data_math, res, v_maxerr_pu_s3w = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=false);

### 1wire - switch
case = "test_switch_1w"
data_eng, data_math, res, v_maxerr_pu_s1w = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=false);

names_switch = ["test_switch", "test_switch_3w", "test_switch_1w"]
v_maxerr_pu_s = [v_maxerr_pu_s4w, v_maxerr_pu_s3w, v_maxerr_pu_s1w]
@show v_maxerr_pu_s

## ############## TRANSFORMER ##############
### 3wire - 3 winding 1 phase -> wye-wye-wye
case = "ut_trans_3w_yyy_1"
data_eng, data_math, res, v_maxerr_pu_t3wyyy = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=false);

### 3wire - 3 winding 3 phase -> delta-wye-wye
case = "ut_trans_3w_dyy_1"
data_eng, data_math, res, v_maxerr_pu_t3wdyy = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=false);

### 3wire - 2 winding -> delta-wye
case = "test_trans_dy_3w"
data_eng, data_math, res, v_maxerr_pu_t3wdy = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=false);

### 3wire - 2 winding -> wye-wye
case = "test_trans_yy_3w"
data_eng, data_math, res, v_maxerr_pu_t3wyy = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=false);

### 4wire - 2 winding -> delta-wye
case = "test_trans_dy"
data_eng, data_math, res, v_maxerr_pu_t4wdy = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);

### 4wire - 2 winding -> wye-wye
case = "test_trans_yy"
data_eng, data_math, res, v_maxerr_pu_t4wyy = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=true);

names_transformer = ["ut_trans_3w_yyy_1", "ut_trans_3w_dyy_1", "test_trans_dy_3w", "test_trans_yy_3w", "test_trans_dy", "test_trans_yy"]
v_maxerr_pu_t = [v_maxerr_pu_t3wyyy, v_maxerr_pu_t3wdyy, v_maxerr_pu_t3wdy, v_maxerr_pu_t3wyy, v_maxerr_pu_t4wdy, v_maxerr_pu_t4wyy]
@show v_maxerr_pu_t

## ############## 3 wire IEEE test cases  ##############  get these files from  https://github.com/sanderclaeys/DistributionTestCases.jl
case = "ieee13_pmd"
data_eng, data_math, res, v_maxerr_pu_IEEE13 = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=false);

case = "ieee34_pmd"
data_eng, data_math, res, v_maxerr_pu_IEEE34 = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=false);

case = "ieee123_pmd"
data_eng, data_math, res, v_maxerr_pu_IEEE123 = solve_compute_mc_pf("$data_dir/$case.dss", "$solution_dir/$case.json"; explicit_neutral=false);

names_IEEE = ["ieee13_pmd", "ieee34_pmd", "ieee123_pmd"]
v_maxerr_pu_IEEE = [v_maxerr_pu_IEEE13, v_maxerr_pu_IEEE34, v_maxerr_pu_IEEE123]
@show v_maxerr_pu_IEEE

## ############## Egrid - Greensboro ##############
egrid_data_dir = "/Users/hei06j/Documents/repositories/remote/PowerModelsDistribution.jl/examples/native_pf_testcases/Industrial"
egrid_solution_file = "../solutions/egrid_GreensBoro.json"
cd(egrid_data_dir)
data_eng, data_math, res, v_maxerr_pu_GreensBoro_industrial = solve_compute_mc_pf("master.dss", egrid_solution_file; explicit_neutral=false, max_iter=1000);
@show v_maxerr_pu_GreensBoro_industrial
cd("../../../")


## ############## Egrid - SantaFe ##############
egrid_data_dir = "/Users/hei06j/Documents/repositories/remote/PowerModelsDistribution.jl/examples/native_pf_testcases/EgridData - SantaFe/urban-suburban/uhs0_1247/uhs0_1247--udt4776"
egrid_solution_file = "../../../../solutions/egrid_SantaFe.json"
cd(egrid_data_dir)
data_eng, data_math, res, v_maxerr_pu_SantaFE_uh0_1 = solve_compute_mc_pf("Master.dss", egrid_solution_file; explicit_neutral=false, max_iter=1000);
@show v_maxerr_pu_SantaFE_uh0_1
cd("../../../../../../")


names = vcat(names_gen, names_load, names_switch, names_transformer, names_IEEE, "GreensBoro_industrial", "SantaFE_uh0_1")
v_maxerr_pu = vcat(v_maxerr_pu_g, v_maxerr_pu_l, v_maxerr_pu_s, v_maxerr_pu_t, v_maxerr_pu_IEEE, v_maxerr_pu_GreensBoro_industrial, v_maxerr_pu_SantaFE_uh0_1)