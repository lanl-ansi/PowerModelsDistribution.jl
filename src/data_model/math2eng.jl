""
function solution_math2eng(solution_math::Dict{String,<:Any}, data_math::Dict{String,<:Any}; make_si::Bool=true, make_deg::Bool=true)::Dict{String,Any}
    solution_eng = Dict{String, Any}()

    if make_deg || make_si
        solution_math = solution_make_si(solution_math, data_math; mult_vbase=make_si, mult_sbase=make_si, convert_rad2deg=make_deg)
    end

    map_keys = sort(collect(keys(data_math["map"])); rev=true)
    for map_key in map_keys
        map = data_math["map"][map_key]
        reverse_bus_lookup = Dict{Int,Any}((bus_math, bus_eng) for (bus_eng, bus_math) in data_math["bus_lookup"])
        getfield(PowerModelsDistribution, map[:unmap_function])(solution_eng, solution_math, map, reverse_bus_lookup)
    end

    for (k,v) in solution_eng
        if isempty(v)
            delete!(solution_eng, k)
        end
    end

    return solution_eng
end


""
function _map_math2eng_voltage_source!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any})
    if !haskey(data_eng, "voltage_source")
        data_eng["voltage_source"] = Dict{Any,Any}()
    end

    eng_obj = _init_unmap_eng_obj!(data_eng, "voltage_source", map)

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

    if !isempty(eng_obj)
        data_eng["voltage_source"][map[:from]] = eng_obj
    end
end


""
function _map_math2eng_bus!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "bus", map)
    math_obj = _get_math_obj(data_math, map[:to])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["bus"][map[:from]] = eng_obj
    end
end


""
function _map_math2eng_load!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "load", map)
    math_obj = _get_math_obj(data_math, map[:to])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["load"][map[:from]] = eng_obj
    end
end


function _map_math2eng_shunt_capacitor!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "shunt_capacitor", map)
    math_obj = _get_math_obj(data_math, map[:to])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["shunt_capacitor"][map[:from]] = eng_obj
    end
end


function _map_math2eng_shunt!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any})
end


function _map_math2eng_generator!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "generator", map)
    math_obj = _get_math_obj(data_math, map[:to])

    @warn "gen" math_obj eng_obj

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["generator"][map[:from]] = eng_obj
    end
end


function _map_math2eng_solar!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any})
end


function _map_math2eng_storage!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any})
end


function _map_math2eng_line!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "line", map)
    math_obj = _get_math_obj(data_math, map[:to])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["line"][map[:from]] = eng_obj
    end
end


function _map_math2eng_line_reactor!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any})
end


function _map_math2eng_switch!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any})
end


function _map_math2eng_transformer!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "transformer", map)

    trans_2wa_ids = [index for (comp_type, index) in split.(map[:to], ".", limit=2) if comp_type=="transformer"]

    prop_map = Dict("pf"=>"p", "qf"=>"q", "crt_fr"=>"crt", "cit_fr"=>"cit")
    for (prop_from, prop_to) in prop_map
        if any(haskey(data_math["transformer"][id], prop_from) for id in trans_2wa_ids)
            eng_obj[prop_to] = [get(data_math["transformer"][id], prop_from, NaN) for id in trans_2wa_ids]
        end
    end

    if !isempty(eng_obj)
        data_eng["transformer"][map[:from]] = eng_obj
    end
end


function _map_math2eng_root!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{Symbol,<:Any}, bus_lookup::Dict{Int,<:Any})
    data_eng["per_unit"] = data_math["per_unit"]

    data_eng["settings"] = Dict{String,Any}("sbase" => data_math["baseMVA"])
end
