""
function _map_math2eng!(data_math)
    @assert get(data_math, "data_model", "mathematical") == "mathematical" "Cannot map data to engineering model: provided data is not a mathematical model"
    @assert haskey(data_math, "map") "Cannot map data to engineering model: no mapping from mathematical to engineering data model is provided"

    data_eng = Dict{String,Any}()

    unmap_keys = sort(keys(data_math["map"]), rev=true)
    reverse_bus_lookup = Dict{Int,Any}((bus_math, bus_eng) for (bus_eng, bus_math) in data_math["bus_lookup"])
    for key in unmap_keys
        getfield(PowerModelsDistribution, map[:unmap_function])(data_eng, data_math, data_math["map"][key], reverse_bus_lookup)
    end

    return data_eng
end


# MAP SOLUTION UP
""
function solution_math2eng(solution_math::Dict, data_math::Dict; make_si::Bool=true, make_deg::Bool=true)
    solution_eng = Dict{String, Any}()

    if make_deg || make_si
        solution_math = solution_make_si(solution_math, data_math; mult_vbase=make_si, mult_sbase=make_si, convert_rad2deg=make_deg)
    end

    map_keys = sort(collect(keys(data_math["map"])); rev=true)
    for map_key in map_keys
        map = data_math["map"][map_key]
        reverse_bus_lookup = Dict{Int,Any}((bus_math, bus_eng) for (bus_eng, bus_math) in data_math["bus_lookup"])
        getfield(PowerModelsDistribution, map[:unmap_function])(solution_eng, solution_math, map, reverse_bus_lookup; map_solution=true)
    end

    return solution_eng
end


function _map_math2eng_voltage_source!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any}; map_solution::Bool=false)
    if !haskey(data_eng, "voltage_source")
        data_eng["voltage_source"] = Dict{Any,Any}()
    end

    eng_obj = _init_unmap_eng_obj!(data_eng, "voltage_source", map; map_solution=map_solution)

    if map_solution
        for to_id in map[:to]
            math_obj = _get_math_obj(data_math, to_id)
            if startswith(to_id, "gen")
                for property in ["pg", "qg", "pg_bus", "qg_bus"]
                    if haskey(math_obj, property)
                        eng_obj[property] = math_obj[property]
                    end
                end
            end
        end
    else
        for to_id in map[:to]
            math_obj = _get_math_obj(data_math, to_id)
            if startswith(to_id, "bus")
                eng_obj["vm"] = math_obj["vm"]
                eng_obj["va"] = math_obj["va"]

            elseif startswith(to_id, "gen")
                eng_obj["source_id"] = strip(math_obj["source_id"], "_virtual_gen.")

            elseif startswith(to_id, "branch")
                eng_obj["bus"] = bus_lookup[math_obj["t_bus"]]
                eng_obj["rs"] = math_obj["br_r"]
                eng_obj["xs"] = math_obj["br_x"]
            else
                Memento.warn(_LOGGER, "transforming from $to_id to $(map[:from]) is not supported")
            end
        end
    end

    data_eng["voltage_source"][map[:from]] = eng_obj
end


function _map_math2eng_bus!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any}; map_solution::Bool=false)
    eng_obj = _init_unmap_eng_obj!(data_eng, "bus", map; map_solution=map_solution)
    math_obj = _get_math_obj(data_math, map[:to])

    if map_solution
        eng_obj = math_obj
    else
        eng_obj["status"] = math_obj["bus_type"] == 4 ? 0 : 1

        for key in ["vm", "va", "vmin", "vmax"]
            if haskey(math_obj, key)
                eng_obj[key] = math_obj[key]
            end
        end
    end
    data_eng["bus"][map[:from]] = eng_obj
end


function _map_math2eng_load!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any}; map_solution::Bool=false)
    eng_obj = _init_unmap_eng_obj!(data_eng, "load", map; map_solution=map_solution)
    math_obj = _get_math_obj(data_math, map[:to])

    if map_solution
        eng_obj = math_obj
    else
        eng_obj["vnom"] = math_obj["vnom_kv"]
        eng_obj["bus"] = bus_lookup[math_obj["load_bus"]]

        eng_obj["pd"] = math_obj["pd"]
        eng_obj["qd"] = math_obj["qd"]

        eng_obj["status"] = math_obj["status"]

        eng_obj["configuration"] = math_obj["configuration"]
    end

    data_eng["load"][map[:from]] = eng_obj
end


function _map_math2eng_shunt_capacitor!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any}; map_solution::Bool=false)
    eng_obj = _init_unmap_eng_obj!(data_eng, "shunt_capacitor", map; map_solution=map_solution)
    math_obj = _get_math_obj(data_math, map[:to])

    if map_solution
        eng_obj = math_obj
    else
        eng_obj["bus"] = bus_lookup["shunt_bus"]
        eng_obj["status"] = math_obj["status"]

        eng_obj["configuration"] = eng_obj["configuration"]

        B = math_obj["bs"]

        if math_obj["configuration"] == "wye"
            b_cap_pu = diag(B)
        else
            # TODO
        end

        # TODO how to pull kv, kvar back out, this is a destructive transformation
    end

    data_eng["shunt_capacitor"][map[:from]] = eng_obj
end


function _map_math2eng_shunt!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any}; map_solution::Bool=false)
end


function _map_math2eng_generator!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any}; map_solution::Bool=false)
    eng_obj = _init_unmap_eng_obj!(data_eng, "generator", map)
    math_obj = _get_math_obj(data_math, map[:to])

    eng_obj["bus"] = bus_lookup[math_obj["gen_bus"]]
    eng_obj["configuration"] = math_obj["configuration"]
    eng_obj["status"] = math_obj["gen_status"]

    eng_obj["kw"] = math_obj["pg"]
    eng_obj["kvar"] = math_obj["qg"]
    eng_obj["kv"] = math_obj["vg"]

    eng_obj["kvar_min"] = math_obj["qmin"]
    eng_obj["kvar_max"] = math_obj["qmax"]

    eng_obj["kw_min"] = math_obj["pmin"]
    eng_obj["kw_max"] = math_obj["pmax"]

    gen_bus = data_math["bus"]["$(math_obj["gen_bus"])"]
    if gen_bus["bus_type"] == 2
        eng_obj["control_mode"] = 3
    else
        eng_obj["control_mode"] = 1
    end

    data_eng["generator"][map[:from]] = eng_obj
end


function _map_math2eng_solar!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any}; map_solution::Bool=false)
end


function _map_math2eng_storage!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any}; map_solution::Bool=false)
end


function _map_math2eng_line!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any}; map_solution::Bool=false)
end


function _map_math2eng_line_reactor!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any}; map_solution::Bool=false)
end


function _map_math2eng_switch!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any}; map_solution::Bool=false)
end


function _map_math2eng_transformer!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any}; map_solution::Bool=false)
    eng_obj = _init_unmap_eng_obj!(data_eng, "transformer", map; map_solution=map_solution)

    if map_solution
        trans_2wa_ids = [index for (comp_type, index) in split.(map[:to], ".", limit=2) if comp_type=="transformer"]

        prop_map = Dict("pf"=>"p", "qf"=>"q")
        for (prop_from, prop_to) in prop_map
            eng_obj[prop_to] = [get(data_math["transformer"][id], prop_from, NaN) for id in trans_2wa_ids]
        end
    else
    end

    data_eng["transformer"][map[:from]] = eng_obj
end


function _map_math2eng_root!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any}; map_solution::Bool=false)
    if map_solution
    else
    end
end