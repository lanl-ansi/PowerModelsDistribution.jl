"""
    transform_solution(
        solution_math::Dict{String,<:Any},
        data_math::Dict{String,<:Any};
        map::Union{Vector{<:Dict{String,<:Any}},Missing}=missing,
        make_si::Bool=true,
        convert_rad2deg::Bool=true,
        map_math2eng_extensions::Dict{String,<:Function}=Dict{String,Function}(),
        make_si_extensions::Vector{<:Function}=Function[],
        dimensionalize_math_extensions::Dict{String,<:Dict{String,<:Vector{<:String}}}=Dict{String,Dict{String,Vector{String}}}(),
        )::Dict{String,Any}

Transforms solutions from MATHEMATICAL data structures, back to an ENGINEERING data structure, given
a `map::Vector{Dict{String,Any}}`, typically which was produced automatically by
[`transform_data_model`](@ref transform_data_model).

# Notes

If `make_si==false`, the solution will remain in per-unit, rather than being converted back to SI
units (default). Angles will still be converted to degrees unless `convert_rad2deg` is utilized.

If `convert_rad2deg==false`, angles will remain in radians, instead of getting converted to
degrees (default).

## Custom SI unit conversions

See [`solution_make_si`](@ref solution_make_si)

## Custom math2eng transformations

To enable automatically mapping back custom components solutions' to the ENGINEERING structure,
`eng2math_extensions` added in [`transform_data_model`](@ref transform_data_model) should include
a push of an item to the `map` dictionary in the `data_math` structure. These items should have
the structure:

    Dict{String,Any}(
        "from" => String,
        "to" => Union{String,Vector{String}},
        "unmap_function" => PowerModelsDistribution.function!,
        "apply_to_subnetworks" => Bool
    )

Important things to note are that

1. The function must be included in `map_math2eng_extensions`, which has the form:

    ```julia
    Dict{String,Function}(
        "_map_math2eng_func!" => _map_math2eng_func!,
    )
    ```

1. `"apply_to_subnetworks"` is optional, and is true by default.
1. `"from"` needs to be a single object
1. `"to"` can be multiple objects or a single object
"""
function transform_solution(
    solution_math::Dict{String,<:Any},
    data_math::Dict{String,<:Any};
    map::Union{Vector{<:Dict{String,<:Any}},Missing}=missing,
    make_si::Bool=true,
    convert_rad2deg::Bool=true,
    map_math2eng_extensions::Dict{String,<:Function}=Dict{String,Function}(),
    make_si_extensions::Vector{<:Function}=Function[],
    dimensionalize_math_extensions::Dict{String,<:Dict{String,<:Vector{<:String}}}=Dict{String,Dict{String,Vector{String}}}(),
    )::Dict{String,Any}

    @assert ismath(data_math) "cannot be converted to an engineering model"

    # convert solution to si?
    solution_math = solution_make_si(
        solution_math,
        data_math;
        mult_vbase=make_si,
        mult_sbase=make_si,
        convert_rad2deg=convert_rad2deg,
        make_si_extensions=make_si_extensions,
        dimensionalize_math_extensions=dimensionalize_math_extensions,
    )

    if ismultinetwork(data_math)
        nws_math_sol = solution_math["nw"]
        nws_math_data = data_math["nw"]
    else
        nws_math_sol = Dict("0" => solution_math)
        nws_math_data = Dict("0" => data_math)
    end

    nws_eng_sol = Dict(k => Dict{String,Any}() for k in keys(nws_math_sol))
    solution_eng = Dict{String,Any}("nw" => nws_eng_sol)

    map = ismissing(map) ? get(data_math, "map", Vector{Dict{String,Any}}()) : map
    @assert !isempty(map) "Map is empty, cannot map solution up to engineering model"

    # apply unmap functions
    for map_item in reverse(map)
        if map_item["unmap_function"] != "_map_math2eng_root!" && get(map_item, "apply_to_subnetworks", true)
            for (n, nw_math_sol) in nws_math_sol
                if haskey(map_math2eng_extensions, map_item["unmap_function"])
                    map_math2eng_extensions[map_item["unmap_function"]](nws_eng_sol[n], nw_math_sol, map_item)
                else
                    getfield(PowerModelsDistribution, Symbol(map_item["unmap_function"]))(nws_eng_sol[n], nw_math_sol, map_item)
                end
            end
        else
            if haskey(map_math2eng_extensions, map_item["unmap_function"])
                map_math2eng_extensions[map_item["unmap_function"]](solution_eng, solution_math, map_item)
            else
                getfield(PowerModelsDistribution, Symbol(map_item["unmap_function"]))(solution_eng, solution_math, map_item)
            end
        end
    end

    # cleanup empty solution dicts
    for (n,nw_eng_sol) in nws_eng_sol
        for (k,v) in nw_eng_sol
            if isempty(v)
                delete!(nw_eng_sol, k)
            end
        end
    end

    # if !multinetwork, correct solution dict
    if !ismultinetwork(data_math)
        merge!(solution_eng, nws_eng_sol["0"])
        delete!(solution_eng, "nw")
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


""
function _map_math2eng_shunt!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "shunt", map)
    math_obj = _get_math_obj(data_math, map["to"])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["shunt"][map["from"]] = eng_obj
    end
end


""
function _map_math2eng_generator!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "generator", map)
    math_obj = _get_math_obj(data_math, map["to"])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["generator"][map["from"]] = eng_obj
    end
end


""
function _map_math2eng_solar!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "solar", map)
    math_obj = _get_math_obj(data_math, map["to"])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["solar"][map["from"]] = eng_obj
    end
end


""
function _map_math2eng_storage!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "storage", map)
    math_obj = _get_math_obj(data_math, map["to"])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["storage"][map["from"]] = eng_obj
    end
end


""
function _map_math2eng_line!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "line", map)
    math_obj = _get_math_obj(data_math, map["to"])

    merge!(eng_obj, math_obj)

    if !isempty(eng_obj)
        data_eng["line"][map["from"]] = eng_obj
    end
end


""
function _map_math2eng_switch!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    eng_obj = _init_unmap_eng_obj!(data_eng, "switch", map)

    prop_map = Dict{String,String}(
        "pt" => "psw_to",
        "qt" => "qsw_to",
        "cr_to" => "crsw_to",
        "ci_to" => "cisw_to",
        "csr_fr" => "csrsw_fr",
        "csi_fr" => "csisw_fr",
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
                for k in keys(prop_map)
                    if haskey(math_obj, k)
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


""
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


""
function _map_math2eng_root!(data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}, map::Dict{String,<:Any})
    data_eng["per_unit"] = data_math["per_unit"]

    if !ismultinetwork(data_math)
        data_eng["settings"] = data_math["settings"]
    else
        for (n,nw) in data_eng["nw"]
            nw["settings"] = Dict{String,Any}("sbase" => data_math["nw"][n]["settings"]["sbase"])
        end
    end
end
