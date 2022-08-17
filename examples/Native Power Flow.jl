using Pkg
using PowerModelsDistribution

Pkg.activate("./examples")

using PowerModelsDistribution
const PMD = PowerModelsDistribution

##
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


function multinetwork_data_math_correction!(data_math::Dict{String, Any})
    @assert data_math["multinetwork"]
    @assert data_math["data_model"]==MATHEMATICAL
    for (nw, dm) in data_math["nw"]
        dm["data_model"] = MATHEMATICAL
        dm["map"] = data_math["map"]
        dm["bus_lookup"] = data_math["bus_lookup"][nw]
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
    data_math["conductor_ids"] = data_math["conductor_ids"][1:3]

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


function find_source_bus_id(math)
    for (b, bus) in math["bus"]
        if bus["name"] == "sourcebus"
            return parse(Int, b)
        end
    end
    return error()
end



## 3 wire

### get these files from  https://github.com/sanderclaeys/DistributionTestCases.jl
# case_path = "./examples/IEEE testcases Sander/ieee13_pmd.dss"
# case_path = "./examples/IEEE testcases Sander/ieee34_pmd.dss"
case_path = "./examples/IEEE testcases Sander/ieee123_pmd.dss"

data_eng = parse_file(case_path, transformations=[transform_loops!, remove_all_bounds!])
data_eng["is_kron_reduced"] = true
data_eng["settings"]["sbase_default"] = 1
vsource_correction_to_3w!(data_eng)

data_math = transform_data_model(data_eng;kron_reduce=false, phase_project=false)
sourcebus_voltage_vector_correction!(data_math, explicit_neutral=false)
update_math_model_3wire!(data_math)

res = PowerModelsDistribution.compute_pf(data_math; explicit_neutral=false)


## 4 wire

case_path = "test/data/en_validation_case_data/test_trans_dy.dss"

data_eng = parse_file(case_path, transformations=[transform_loops!, remove_all_bounds!])
vsource_correction_to_4w!(data_eng)

data_math = transform_data_model(data_eng;kron_reduce=false, phase_project=false)

res = PowerModelsDistribution.compute_pf(data_math; explicit_neutral=true)


## adapt the script by comparing a 4wire and 3wire testcase -> add_start_voltage!(dm, coordinates=:rectangular, epsilon=0) ???


case_path = "./test/data/opendss/ut_trans_2w_dy_lag.dss"
data_eng = parse_file(case_path, transformations=[transform_loops!])
vsource_correction_to_4w!(data_eng)
data_math = transform_data_model(data_eng;kron_reduce=false)
res = PowerModelsDistribution.compute_pf(data_math; explicit_neutral=true)




##
case_path_4w = "./examples/ENWL testcases/four-wire with explicit_neutral"
cd(case_path_4w)
data_eng = parse_file("Master.dss", transformations=[transform_loops!])
vsource_correction_to_4w!(data_eng)
data_math = transform_data_model(data_eng;kron_reduce=false)
res = PowerModelsDistribution.compute_pf(data_math; explicit_neutral=true)

cd("../three-wire")
data_eng = parse_file("Master.dss", transformations=[transform_loops!])
data_eng["is_kron_reduced"] = true
vsource_correction_to_3w!(data_eng)
data_math = transform_data_model(data_eng;kron_reduce=false)
sourcebus_voltage_vector_correction!(data_math, explicit_neutral=false)
update_math_model_3wire!(data_math)
res = PowerModelsDistribution.compute_pf(data_math; explicit_neutral=false)


##
# accessing OpenDSSDirect for results consistency
# validating the results are consistent with OpenDSS
# double check what is actually happeind with the aux current and voltage variables in ideal transformers
# fixing the unit tests
# double checking KLU settings (default or not?) with respect to OpenDSS settings