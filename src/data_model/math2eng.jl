""
function transform_solution(solution_math::Dict{String,<:Any}, data_math::Dict{String,<:Any}; map::Union{Vector{<:Dict{String,<:Any}},Missing}=missing, make_si::Bool=true, convert_rad2deg::Bool=true)::Dict{String,Any}
    @assert get(data_math, "data_model", MATHEMATICAL) == MATHEMATICAL "provided solution cannot be converted to an engineering model"
    if ismultinetwork(data_math)
        solution_eng = Dict{String,Any}(
            "nw" => Dict{String,Any}(
                k => Dict{String,Any}() for k in keys(data_math["nw"])
            )
        )
    else
        solution_eng = Dict{String,Any}()
    end

    solution_math = solution_make_si(solution_math, data_math; mult_vbase=make_si, mult_sbase=make_si, convert_rad2deg=convert_rad2deg)

    map = ismissing(map) ? get(data_math, "map", Vector{Dict{String,Any}}()) : map
    @assert !isempty(map) "Map is empty, cannot map solution up to engineering model"

    for map_item in reverse(map)
        if ismultinetwork(data_math) && map_item["unmap_function"] != "_map_math2eng_root!"
            for (n, nw) in solution_math["nw"]
                getfield(PowerModelsDistribution, Symbol(map_item["unmap_function"]))(solution_eng["nw"][n], nw, map_item)
            end
        else
            getfield(PowerModelsDistribution, Symbol(map_item["unmap_function"]))(solution_eng, solution_math, map_item)
        end
    end

    if ismultinetwork(data_math)
        for (n,nw) in solution_eng["nw"]
            for (k,v) in nw
                if isempty(v)
                    delete!(nw, k)
                end
            end
        end
    else
        for (k,v) in solution_eng
            if isempty(v)
                delete!(solution_eng, k)
            end
        end
    end

    return solution_eng
end


""
function _map_math2eng_voltage_source!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    if !haskey(data_eng, "voltage_source")
        data_eng["voltage_source"] = Dict{String,Any}()
    end

    eng_obj = _init_unmap_eng_obj!(data_eng, "voltage_source", map)

    map["to"] = isa(map["to"], Vector) ? map["to"] : [map["to"]]

    for to_id in map["to"]
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
        data_eng["voltage_source"][map["from"]] = eng_obj
    end
end


""
function _map_math2eng_bus!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "bus", map)
    math_obj = _get_math_obj(data_math, map["to"])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["bus"][map["from"]] = eng_obj
    end
end


""
function _map_math2eng_load!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "load", map)
    math_obj = _get_math_obj(data_math, map["to"])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["load"][map["from"]] = eng_obj
    end
end


function _map_math2eng_shunt!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "shunt", map)
    math_obj = _get_math_obj(data_math, map["to"])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["shunt"][map["from"]] = eng_obj
end
end


function _map_math2eng_generator!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "generator", map)
    math_obj = _get_math_obj(data_math, map["to"])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["generator"][map["from"]] = eng_obj
    end
end


function _map_math2eng_solar!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "solar", map)
    math_obj = _get_math_obj(data_math, map["to"])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["solar"][map["from"]] = eng_obj
    end
end


function _map_math2eng_storage!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "storage", map)
    math_obj = _get_math_obj(data_math, map["to"])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["storage"][map["from"]] = eng_obj
    end
end


function _map_math2eng_line!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "line", map)
    math_obj = _get_math_obj(data_math, map["to"])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["line"][map["from"]] = eng_obj
    end
end


function _map_math2eng_switch!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "switch", map)

    prop_map = Dict{String,String}(
        "pf" => "psw_fr",
        "qf" => "qsw_fr",
        "cr_fr" => "crsw_fr",
        "ci_fr" => "cisw_fr"
    )

    if isa(map["to"], String)
        math_obj = _get_math_obj(data_math, map["to"])
        merge!(eng_obj, math_obj)
    else
        for to_id in map["to"]
            if startswith(to_id, "switch")
                math_obj = _get_math_obj(data_math, to_id)
                merge!(eng_obj, math_obj)
            elseif startswith(to_id, "branch")
                math_obj = _get_math_obj(data_math, to_id)
                for k in ["pf", "qf", "cr_fr", "ci_fr"]
                    if haskey(eng_obj, prop_map[k]) && haskey(math_obj, k)
                        eng_obj[prop_map[k]] = math_obj[k]
                    end
                end
            end
        end
    end

    if haskey(eng_obj, "state")
        eng_obj["state"] = SwitchState(Int(round(eng_obj["state"])))
    end

    if !isempty(eng_obj)
        data_eng["switch"][map["from"]] = eng_obj
    end
end


function _map_math2eng_transformer!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "transformer", map)

    trans_2wa_ids = [index for (comp_type, index) in split.(map["to"], ".", limit=2) if comp_type=="transformer"]

    prop_map = Dict("pf"=>"p", "qf"=>"q", "crt_fr"=>"crt", "cit_fr"=>"cit")
    for (prop_from, prop_to) in prop_map
        if haskey(data_math, "transformer")
            if any(haskey(data_math["transformer"][id], prop_from) for id in trans_2wa_ids)
                eng_obj[prop_to] = [get(data_math["transformer"][id], prop_from, NaN) for id in trans_2wa_ids]
            end
        end
    end

    if !isempty(eng_obj)
        data_eng["transformer"][map["from"]] = eng_obj
    end
end


function _map_math2eng_root!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    data_eng["per_unit"] = data_math["per_unit"]

    if !ismultinetwork(data_math)
        data_eng["settings"] = Dict{String,Any}("sbase" => data_math["baseMVA"])
    else
        for (n,nw) in data_eng["nw"]
            nw["settings"] = Dict{String,Any}("sbase" => data_math["nw"][n]["baseMVA"])
        end
    end
end
