using Pkg
using PowerModelsDistribution

const PMD = PowerModelsDistribution

using JSON

##
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


function vsource_correction_to_4w!(data_eng)
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
    return nothing
end

function vsource_correction_to_3w!(data_eng)
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
            if haskey(bus, "terminals")
                bus["terminals"] = bus["terminals"][1:end-1]
            end
            if haskey(bus, "grounded")
                bus["grounded"] = bus["grounded"][1:end-1]
            end
            if haskey(bus, "vmax")
                bus["vmax"] = bus["vmax"][1:end-1]
            end
            if haskey(bus, "vmin")
                bus["vmin"] = bus["vmin"][1:end-1]
            end
            bus["vmin"] = 0.9 * ones(length(bus["terminals"]))
            bus["vmax"] = 1.1 * ones(length(bus["terminals"]))
        end
    end

    for (l,branch) in math["branch"]
        explicit_neutral = false
        if haskey(branch, "t_connections") && neutral_idx ∈ branch["t_connections"]
            explicit_neutral = true
            branch["t_connections"] = branch["t_connections"][1:end-1]
        end
        if haskey(branch, "f_connections") && neutral_idx ∈ branch["f_connections"]
            explicit_neutral = true
            branch["f_connections"] = branch["f_connections"][1:end-1]
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
            transformer["t_connections"] = transformer["t_connections"][1:end-1]
        end
        if haskey(transformer, "f_connections") && neutral_idx ∈ transformer["f_connections"]
            transformer["f_connections"] = transformer["f_connections"][1:end-1]
        end
    end

    for (g,gen) in math["gen"]
        if neutral_idx in gen["connections"]
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
        if load["configuration"] == WYE && neutral_idx ∈ load["connections"]
            load["connections"] = load["connections"][1:end-1]
        end
    end

    return nothing
end



## 3 wire   ---  get these files from  https://github.com/sanderclaeys/DistributionTestCases.jl
data_dir = "examples/native_pf_testcases"
solution_dir = "examples/native_pf_testcases/solutions"

case = "ieee13_pmd"
# case = "ieee34_pmd"
# case = "ieee123_pmd"

data_eng = parse_file("$data_dir/$case.dss", transformations=[transform_loops!])
data_eng["is_kron_reduced"] = true
data_eng["settings"]["sbase_default"] = 1
vsource_correction_to_3w!(data_eng)
data_math = transform_data_model(data_eng;kron_reduce=false, phase_project=false)
sourcebus_voltage_vector_correction!(data_math, explicit_neutral=false)
update_math_model_3wire!(data_math)
res = compute_pf(data_math; explicit_neutral=false)
# obtain solution from dss
sol_dss = open("$solution_dir/$case.json", "r") do f
    JSON.parse(f)
end
sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)



## 3wire and 4wire
case = "test_trans_dy"
data_eng = parse_file("$data_dir/$case.dss", transformations=[transform_loops!])
vsource_correction_to_4w!(data_eng)
data_math = transform_data_model(data_eng;kron_reduce=false)
res = compute_pf(data_math; explicit_neutral=true)
# obtain solution from dss
sol_dss = open("$solution_dir/$case.json", "r") do f
    JSON.parse(f)
end
sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)

##
case = "test_trans_yy"
data_eng = parse_file("$data_dir/$case.dss", transformations=[transform_loops!])
vsource_correction_to_4w!(data_eng)
data_math = transform_data_model(data_eng;kron_reduce=false)
res = compute_pf(data_math; explicit_neutral=true)
# obtain solution from dss
sol_dss = open("$solution_dir/$case.json", "r") do f
    JSON.parse(f)
end
sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)

##
case = "test_trans_dy_3w"
data_eng = parse_file("$data_dir/$case.dss", transformations=[transform_loops!])
data_eng["is_kron_reduced"] = true
data_eng["settings"]["sbase_default"] = 1
vsource_correction_to_3w!(data_eng)
data_math = transform_data_model(data_eng;kron_reduce=false, phase_project=false)
sourcebus_voltage_vector_correction!(data_math, explicit_neutral=false)
update_math_model_3wire!(data_math)
res = compute_pf(data_math; explicit_neutral=false)
# obtain solution from dss
sol_dss = open("$solution_dir/$case.json", "r") do f
    JSON.parse(f)
end
sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)

##
case = "test_trans_yy_3w"
data_eng = parse_file("$data_dir/$case.dss", transformations=[transform_loops!])
data_eng["is_kron_reduced"] = true
data_eng["settings"]["sbase_default"] = 1
vsource_correction_to_3w!(data_eng)
data_math = transform_data_model(data_eng;kron_reduce=false, phase_project=false)
sourcebus_voltage_vector_correction!(data_math, explicit_neutral=false)
update_math_model_3wire!(data_math)
res = compute_pf(data_math; explicit_neutral=false)
# obtain solution from dss
sol_dss = open("$solution_dir/$case.json", "r") do f
    JSON.parse(f)
end
sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)


##
# run the unit tests with KLU, and with original factorize
# accessing OpenDSSDirect for results consistency
# validating the results are consistent with OpenDSS
# double check what is actually happeind with the aux current and voltage variables in ideal transformers
# fixing the unit tests
# double checking KLU settings (default or not?) with respect to OpenDSS settings


## solution builder - > make sure load currents are in dict (phase->cr, ci) format

