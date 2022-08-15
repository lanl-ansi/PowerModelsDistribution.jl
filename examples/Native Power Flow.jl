using Pkg
Pkg.activate("./examples")

using PowerModelsDistribution
const PMD = PowerModelsDistribution

# case_path = "./test/data/opendss/ut_trans_2w_yy.dss"

### get these files from  https://github.com/sanderclaeys/DistributionTestCases.jl
# case_path = "./examples/Sander's IEEE testcases/ieee13_pmd.dss"
# case_path = "./examples/Sander's IEEE testcases/ieee34_pmd.dss"
case_path = "./examples/Sander's IEEE testcases/ieee123_pmd.dss"

##
function vsource_correction!(data_eng)
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


function sourcebus_voltage_vector_correction!(data_math::Dict{String, Any})
    if haskey(data_math, "multinetwork")
        for (n,nw) in data_math["nw"]
            for (i, bus) in data_math["nw"]["bus"]
                if bus["bus_type"] == 3 && length(bus["terminals"]) != length(bus["vm"])
                    bus["vm"] = bus["vm"][1:length(bus["terminals"])]
                    bus["va"] = bus["va"][1:length(bus["terminals"])]
                end
            end
        end
    else
        for (i, bus) in data_math["bus"]
            if bus["bus_type"] == 3 && length(bus["terminals"]) != length(bus["vm"])
                bus["vm"] = bus["vm"][1:length(bus["terminals"])]
                bus["va"] = bus["va"][1:length(bus["terminals"])]
            end
        end
    end
    return nothing
end

##

data_eng = parse_file(case_path, transformations=[transform_loops!])
vsource_correction!(data_eng)

data_math = transform_data_model(data_eng;kron_reduce=false)
sourcebus_voltage_vector_correction!(data_math)
# data_math["bus"]["58"]["vm"] = data_math["bus"]["58"]["vm"][1:4]
# data_math["bus"]["58"]["va"] = data_math["bus"]["58"]["va"][1:4]

res = PowerModelsDistribution.compute_pf(data_math; explicit_neutral=true)



## adapt the script by comparing a 4wire and 3wire testcase -> add_start_voltage!(dm, coordinates=:rectangular, epsilon=0) ???
