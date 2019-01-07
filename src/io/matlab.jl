current_version = v"1"

""
function parse_matlab(io::IOStream)
    data_string = read(io, String)
    ml_data = parse_matlab_string(data_string)

    pm_data = matlab_to_tppm(ml_data)

    return pm_data
end


""
function parse_matlab(file_string::String)
    pm_data = open(file_string) do io
        parse_matlab(io)
    end
    return pm_data
end


### Data and functions specific to Matlab format ###
tppm_data_names = [
    "tppmc.baseMVA", "tppmc.baseKV", "tppmc.bus", "tppmc.load", "tppmc.shunt",
    "tppmc.gen", "tppmc.branch", "tppmc.bus_name", "tppmc.gencost"
]

tppm_bus_columns = [
    ("bus_i", Int),
    ("bus_type", Int),
    ("vmin_1", Float64), ("vmax_1", Float64),
    ("vmin_2", Float64), ("vmax_2", Float64),
    ("vmin_3", Float64), ("vmax_3", Float64),
    ("vm_1", Float64), ("va_1", Float64),
    ("vm_2", Float64), ("va_2", Float64),
    ("vm_3", Float64), ("va_3", Float64)
]

tppm_bus_name_columns = [
    ("bus_name", AbstractString)
]

tppm_load_columns = [
    ("load_bus", Int),
    ("pd_1", Float64), ("qd_1", Float64),
    ("pd_2", Float64), ("qd_2", Float64),
    ("pd_3", Float64), ("qd_3", Float64),
    ("status", Int)
]

tppm_shunt_columns = [
    ("shunt_bus", Int),
    ("gs_1", Float64), ("bs_1", Float64),
    ("gs_2", Float64), ("bs_2", Float64),
    ("gs_3", Float64), ("bs_3", Float64),
    ("status", Int)
]

tppm_gen_columns = [
    ("gen_bus", Int),
    ("pmin_1", Float64), ("pmax_1", Float64),
    ("qmin_1", Float64), ("qmax_1", Float64),
    ("pmin_2", Float64), ("pmax_2", Float64),
    ("qmin_2", Float64), ("qmax_2", Float64),
    ("pmin_3", Float64), ("pmax_3", Float64),
    ("qmin_3", Float64), ("qmax_3", Float64),
    ("pg_1", Float64), ("qg_1", Float64),
    ("pg_2", Float64), ("qg_2", Float64),
    ("pg_3", Float64), ("qg_3", Float64),
    ("gen_status", Int)
]

tppm_branch_columns = [
    ("f_bus", Int),
    ("t_bus", Int),
    ("r_11", Float64), ("x_11", Float64),
    ("r_12", Float64), ("x_12", Float64),
    ("r_13", Float64), ("x_13", Float64),
    ("r_22", Float64), ("x_22", Float64),
    ("r_23", Float64), ("x_23", Float64),
    ("r_33", Float64), ("x_33", Float64),
    ("b_1", Float64),
    ("b_2", Float64),
    ("b_3", Float64),
    ("rate_a", Float64),
    ("rate_b", Float64),
    ("rate_c", Float64),
    ("angmin", Float64), ("angmax", Float64),
    ("br_status", Int)
]


""
function parse_matlab_string(data_string::String)
    matlab_data, func_name, colnames = InfrastructureModels.parse_matlab_string(data_string, extended=true)

    #println(matlab_data)

    case = Dict{String,Any}()

    if func_name != nothing
        case["name"] = func_name
    else
        warn(LOGGER, string("no case name found in matlab file.  The file seems to be missing \"function mpc = ...\""))
        case["name"] = "no_name_found"
    end

    case["source_type"] = "matlab"
    if haskey(matlab_data, "tppmc.version")
        case["source_version"] = VersionNumber(matlab_data["tppmc.version"])
    else
        warn(LOGGER, "No version number found, file may not be compatible with parser.")
        case["source_version"] = v"0"
    end

    if haskey(matlab_data, "tppmc.baseMVA")
        case["baseMVA"] = matlab_data["tppmc.baseMVA"]
    else
        warn(LOGGER, string("no baseMVA found in matlab file.  The file seems to be missing \"tppmc.baseMVA = ...\""))
        case["baseMVA"] = 1.0
    end

    if haskey(matlab_data, "tppmc.baseKV")
        case["baseKV"] = matlab_data["tppmc.baseKV"]
    else
        warn(LOGGER, string("no baseKV found in matlab file.  The file seems to be missing \"tppmc.baseKV = ...\""))
        case["baseKV"] = 1.0
    end

    if haskey(matlab_data, "tppmc.bus")
        buses = []
        for bus_row in matlab_data["tppmc.bus"]
            bus_data = InfrastructureModels.row_to_typed_dict(bus_row, tppm_bus_columns)
            bus_data["index"] = InfrastructureModels.check_type(Int, bus_row[1])
            push!(buses, bus_data)
        end
        case["bus"] = buses
    else
        error(string("no bus table found in matlab file.  The file seems to be missing \"tppmc.bus = [...];\""))
    end

    if haskey(matlab_data, "tppmc.load")
        loads = []
        for (i,load_row) in enumerate(matlab_data["tppmc.load"])
            load_data = InfrastructureModels.row_to_typed_dict(load_row, tppm_load_columns)
            load_data["index"] = i
            push!(loads, load_data)
        end
        case["load"] = loads
    else
        error(string("no load table found in matlab file.  The file seems to be missing \"tppmc.load = [...];\""))
    end

    if haskey(matlab_data, "tppmc.shunt")
        shunts = []
        for (i,shunt_row) in enumerate(matlab_data["tppmc.shunt"])
            shunt_data = InfrastructureModels.row_to_typed_dict(shunt_row, tppm_shunt_columns)
            shunt_data["index"] = i
            push!(shunts, shunt_data)
        end
        case["shunt"] = shunts
    end

    if haskey(matlab_data, "tppmc.gen")
        gens = []
        for (i, gen_row) in enumerate(matlab_data["tppmc.gen"])
            gen_data = InfrastructureModels.row_to_typed_dict(gen_row, tppm_gen_columns)
            gen_data["index"] = i
            push!(gens, gen_data)
        end
        case["gen"] = gens
    else
        error(string("no gen table found in matlab file.  The file seems to be missing \"tppmc.gen = [...];\""))
    end

    if haskey(matlab_data, "tppmc.branch")
        branches = []
        for (i, branch_row) in enumerate(matlab_data["tppmc.branch"])
            branch_data = InfrastructureModels.row_to_typed_dict(branch_row, tppm_branch_columns)
            branch_data["index"] = i
            push!(branches, branch_data)
        end
        case["branch"] = branches
    else
        error(string("no branch table found in matlab file.  The file seems to be missing \"tppmc.branch = [...];\""))
    end


    if haskey(matlab_data, "tppmc.bus_name")
        bus_names = []
        for (i, bus_name_row) in enumerate(matlab_data["tppmc.bus_name"])
            bus_name_data = InfrastructureModels.row_to_typed_dict(bus_name_row, tppm_bus_name_columns)
            bus_name_data["index"] = i
            push!(bus_names, bus_name_data)
        end
        case["bus_name"] = bus_names

        if length(case["bus_name"]) != length(case["bus"])
            error("incorrect Matpower file, the number of bus names ($(length(case["bus_name"]))) is inconsistent with the number of buses ($(length(case["bus"]))).\n")
        end
    end

    if haskey(matlab_data, "tppmc.gencost")
        gencost = []
        for (i, gencost_row) in enumerate(matlab_data["tppmc.gencost"])
            gencost_data = PMs.mp_cost_data(gencost_row)
            gencost_data["index"] = i
            push!(gencost, gencost_data)
        end
        case["gencost"] = gencost

        if length(case["gencost"]) != length(case["gen"])
            error("incorrect matlab file, the number of generator cost functions ($(length(case["gencost"]))) is inconsistent with the number of generators ($(length(case["gen"]))).\n")
        end
    end

    for k in keys(matlab_data)
        if !in(k, tppm_data_names) && startswith(k, "tppmc.")
            case_name = k[7:length(k)]
            value = matlab_data[k]
            if isa(value, Array)
                column_names = []
                if haskey(colnames, k)
                    column_names = colnames[k]
                end
                tbl = []
                for (i, row) in enumerate(matlab_data[k])
                    row_data = PMs.row_to_dict(row, column_names)
                    row_data["index"] = i
                    push!(tbl, row_data)
                end
                case[case_name] = tbl
                info(LOGGER, "extending matlab format with data: $(case_name) $(length(tbl))x$(length(tbl[1])-1)")
            else
                case[case_name] = value
                info(LOGGER, "extending matlab format with constant data: $(case_name)")
            end
        end
    end

    #println("Case:")
    #println(case)

    return case
end


"Translates legacy versions into current version format"
function translate_version!(ml_data::Dict{String,Any})
    # Future Version translation here
    if ml_data["source_version"] == current_version
        return ml_data
    else
        warn(LOGGER, "matlab source data has unrecognized version $(ml_data["source_version"]), cannot translate to version $current_version, parse may be invalid")
        return ml_data
    end
end


"""
Converts a Matlab dict into a ThreePhasePowerModels dict
"""
function matlab_to_tppm(ml_data::Dict{String,Any})
    ml_data = deepcopy(ml_data)

    translate_version!(ml_data)

    ml_data["multinetwork"] = false
    ml_data["per_unit"] = false
    ml_data["conductors"] = 3
    ml_data["dcline"] = []
    ml_data["storage"] = []

    # required default values
    if !haskey(ml_data, "shunt")
        ml_data["shunt"] = []
    end

    ml2pm_bus(ml_data)
    ml2pm_load(ml_data)
    ml2pm_shunt(ml_data)
    ml2pm_gen(ml_data)
    ml2pm_branch(ml_data)

    PMs.merge_bus_name_data(ml_data)
    PMs.merge_generator_cost_data(ml_data)

    PMs.merge_generic_data(ml_data)

    InfrastructureModels.arrays_to_dicts!(ml_data)

    for optional in ["dcline", "load", "shunt", "storage"]
        if length(ml_data[optional]) == 0
            ml_data[optional] = Dict{String,Any}()
        end
    end

    return ml_data
end


"convert raw bus data into arrays"
function ml2pm_bus(data::Dict{String,Any})
    for bus in data["bus"]
        make_mpv!(bus, "vmin", ["vmin_1", "vmin_2", "vmin_3"])
        make_mpv!(bus, "vmax", ["vmax_1", "vmax_2", "vmax_3"])
        make_mpv!(bus,   "vm", ["vm_1", "vm_2", "vm_3"])
        make_mpv!(bus,   "va", ["va_1", "va_2", "va_3"])
    end
end


"convert raw load data into arrays"
function ml2pm_load(data::Dict{String,Any})
    for load in data["load"]
        make_mpv!(load, "pd", ["pd_1", "pd_2", "pd_3"])
        make_mpv!(load, "qd", ["qd_1", "qd_2", "qd_3"])
    end
end


"convert raw shunt data into arrays"
function ml2pm_shunt(data::Dict{String,Any})
    for load in data["shunt"]
        make_mpv!(load, "gs", ["gs_1", "gs_2", "gs_3"])
        make_mpv!(load, "bs", ["bs_1", "bs_2", "bs_3"])
    end
end


"convert raw generator data into arrays"
function ml2pm_gen(data::Dict{String,Any})
    for gen in data["gen"]
        make_mpv!(gen, "pmin", ["pmin_1", "pmin_2", "pmin_3"])
        make_mpv!(gen, "pmax", ["pmax_1", "pmax_2", "pmax_3"])
        make_mpv!(gen, "qmin", ["qmin_1", "qmin_2", "qmin_3"])
        make_mpv!(gen, "qmax", ["qmax_1", "qmax_2", "qmax_3"])
        make_mpv!(gen, "pg", ["pg_1", "pg_2", "pg_3"])
        make_mpv!(gen, "qg", ["qg_1", "qg_2", "qg_3"])
    end
end


"convert raw branch data into arrays"
function ml2pm_branch(data::Dict{String,Any})
    for branch in data["branch"]
        branch["rate_a"] = PMs.MultiConductorVector(branch["rate_a"], 3)
        branch["rate_b"] = PMs.MultiConductorVector(branch["rate_b"], 3)
        branch["rate_c"] = PMs.MultiConductorVector(branch["rate_c"], 3)

        branch["angmin"] = PMs.MultiConductorVector(branch["angmin"], 3)
        branch["angmax"] = PMs.MultiConductorVector(branch["angmax"], 3)

        set_default(branch, "g_fr_1", 0.0)
        set_default(branch, "g_fr_2", 0.0)
        set_default(branch, "g_fr_3", 0.0)

        set_default(branch, "g_to_1", 0.0)
        set_default(branch, "g_to_2", 0.0)
        set_default(branch, "g_to_3", 0.0)

        make_mpv!(branch, "g_fr", ["g_fr_1", "g_fr_2", "g_fr_3"])
        make_mpv!(branch, "g_to", ["g_to_1", "g_to_2", "g_to_3"])

        branch["b_fr"] = PMs.MultiConductorVector([branch["b_1"], branch["b_2"], branch["b_3"]]) / 2.0
        branch["b_to"] = PMs.MultiConductorVector([branch["b_1"], branch["b_2"], branch["b_3"]]) / 2.0

        delete!(branch, "b_1")
        delete!(branch, "b_2")
        delete!(branch, "b_3")

        branch["tap"] = PMs.MultiConductorVector(1.0, 3)
        branch["shift"] = PMs.MultiConductorVector(0.0, 3)
        branch["transformer"] = false

        branch["br_r"] = PMs.MultiConductorMatrix([
            branch["r_11"]     branch["r_12"]/2.0 branch["r_13"]/2.0;
            branch["r_12"]/2.0 branch["r_22"]     branch["r_23"]/2.0;
            branch["r_13"]/2.0 branch["r_23"]/2.0 branch["r_33"];
        ])

        branch["br_x"] = PMs.MultiConductorMatrix([
            branch["x_11"]     branch["x_12"]/2.0 branch["x_13"]/2.0;
            branch["x_12"]/2.0 branch["x_22"]     branch["x_23"]/2.0;
            branch["x_13"]/2.0 branch["x_23"]/2.0 branch["x_33"];
        ])

        for k in ["r_11", "r_12", "r_13", "r_22", "r_23", "r_33",
            "x_11", "x_12", "x_13", "x_22", "x_23", "x_33"]
            delete!(branch, k)
        end
    end
end


"collects several from_keys in an array and sets it to the to_key, removes from_keys"
function make_mpv!(data::Dict{String,Any}, to_key::String, from_keys::Array{String,1})
    @assert !(haskey(data, to_key))
    data[to_key] = PMs.MultiConductorVector([data[k] for k in from_keys])
    for k in from_keys
        delete!(data, k)
    end
end


"checks if the given dict has a value, if not, sets a default value"
function set_default(data::Dict{String,Any}, key::String, default_value)
    if !(haskey(data, key))
        data[key] = default_value
    end
end


