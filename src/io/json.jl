# Parse PowerModels data from JSON exports of PowerModels data structures.
# Customize JSON printing of matrices

import JSON.Serializations: CommonSerialization, StandardSerialization
import JSON.Writer: StructuralContext, show_json
struct PMDSerialization <: CommonSerialization end


"converts julia Version into serializable structure"
function _jsonver2juliaver!(data::Dict{String,<:Any})
    if haskey(data, "source_version") && isa(data["source_version"], Dict)
        data["source_version"] = "$(data["source_version"]["major"]).$(data["source_version"]["minor"]).$(data["source_version"]["patch"])"
    end
end


"recursive helper function for parsing serialized matrices"
function _parse_mats_recursive!(data::Dict{String,<:Any})::Dict{String,Any}
    parse =  haskey(data, "type") && data["type"]=="Matrix" && haskey(data, "eltype") && haskey(data, "value")
    if parse
        return _parse_matrix_value(data["value"], data["eltype"])
    else
        for (k,v) in data
            if isa(v, Dict)
                data[k] = _parse_mats!(v)
            elseif isa(v, Vector) && Dict <: eltype(v)
                data[k] = [isa(x, Dict) ? _parse_mats!(x) : x for x in v]
            end
        end

        return data
    end
end


"parser function for serialized matrices"
function _parse_mats!(data::Dict{String,<:Any})
    stack = Array{Tuple{Any, Any, Any}}([(data, k, v) for (k, v) in data])
    while !isempty(stack)
        (store, k, v) = pop!(stack)

        parse =  isa(v, Dict) && haskey(v, "type") && v["type"]=="Matrix" && haskey(v, "eltype") && haskey(v, "value")
        if parse
            store[k] = _parse_matrix_value(v["value"], v["eltype"])
        elseif isa(v, Dict)
            append!(stack, [(v, xk, xv) for (xk, xv) in v])
        elseif isa(v, Vector)
            store[k] = [x for x in v]
            if Dict <: eltype(v)
                append!(stack, [(v, i, xv) for (i, xv) in enumerate(v) if isa(xv, Dict)])
            end
        end
    end
end


"parses enums from json"
function _parse_enums!(data::Dict{String,<:Any})
    data["data_model"] = DataModel(get(data, "data_model", 1))

    for (root_type, root_value) in data
        if isa(root_value, Dict)
            for (component_id, component) in root_value
                if isa(component, Dict)
                    if haskey(component, "status")
                        component["status"] = Status(component["status"])
                    end

                    if haskey(component, "dispatchable")
                        component["dispatchable"] = Dispatchable(component["dispatchable"])
                    end

                    if haskey(component, "configuration")
                        if isa(component["configuration"], Vector)
                            component["configuration"] = Vector{ConnConfig}([ConnConfig(el) for el in component["configuration"]])
                        else
                            component["configuration"] = ConnConfig(component["configuration"])
                        end
                    end

                    if root_type == "switch" && haskey(component, "state")
                        component["state"] = SwitchState(component["state"])
                    end

                    if root_type == "generator" && haskey(component, "control_mode")
                        component["generator"] = ControlMode(component["control_model"])
                    end

                    if root_type == "load" && haskey(component, "model")
                        component["model"] = LoadModel(component["model"])
                    end

                    if root_type == "shunt" && haskey(component, "model")
                        component["model"] = ShuntModel(component["model"])
                    end
                end
            end
        end
    end
end


"Parses a JSON file into a PMD data structure"
function parse_json(file::String; validate::Bool=false)
    data = open(file) do io
        parse_json(io; filetype=split(lowercase(file), '.')[end], validate=validate)
    end
    return data
end


"Parses a JSON file into a PMD data structure"
function parse_json(io::IO; validate::Bool=false)::Dict{String,Any}
    data = JSON.parse(io)

    _jsonver2juliaver!(data)

    _parse_mats!(data)

    _parse_enums!(data)

    if validate
        correct_network_data!(data)
    end

    return data
end


"prints a PMD data structure into a JSON file"
function print_file(path::String, data::Dict{String,<:Any}; indent::Int=2)
    open(path, "w") do io
        print_file(io, data; indent=indent)
    end
end


"prints a PMD data structure into a JSON file"
function print_file(io::IO, data::Dict{String,<:Any}; indent::Int=2)
    if indent == 0
        JSON.print(io, JSON.parse(sprint(show_json, PMDSerialization(), data)))
    else
        JSON.print(io, JSON.parse(sprint(show_json, PMDSerialization(), data)), indent)
    end
end


"turns a matrix into a serializable structure"
function show_json(io::StructuralContext, ::PMDSerialization, f::Matrix{<:Any})
    N, M = size(f)
    value = string("[", join([join([f[i,j] for j in 1:M], " ") for i in 1:N], "; "), "]")
    eltyp = isempty(f) ? eltype(f) : typeof(f[1,1])
    out = Dict(:type=>:Matrix, :eltype=>eltyp, :value=>value)
    return show_json(io, StandardSerialization(), out)
end


"custom handling for enums output to json"
function show_json(io::StructuralContext, ::CommonSerialization, f::PowerModelsDistributionEnums)
    return show_json(io, StandardSerialization(), Int(f))
end


"custom handling for enums output to json"
JSON.lower(p::PowerModelsDistributionEnums) = Int(p)


"parses in a serialized matrix"
function _parse_matrix_value(value::String, eltyp::String)
    if value=="[]"
        eltyp =
        return Array{eval(Symbol(eltyp)), 2}(undef, 0, 0)
    else
        tmp = [split(strip(row_str), " ") for row_str in split(value[2:end-1], ";")]
        tmp = [tmp[i][j] for i in 1:length(tmp), j in 1:length(tmp[1])]
        M = parse.(eval(Symbol(eltyp)), string.(tmp))
        return M
    end
end
