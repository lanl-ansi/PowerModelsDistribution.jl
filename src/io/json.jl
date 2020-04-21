# Parse PowerModels data from JSON exports of PowerModels data structures.
# Customize JSON printing of matrices

import JSON.Serializations: CommonSerialization, StandardSerialization
import JSON.Writer: StructuralContext, show_json
struct PMDSerialization <: CommonSerialization end


function _jsonver2juliaver!(data::Dict{String,<:Any})
    if haskey(data, "source_version") && isa(data["source_version"], Dict)
        data["source_version"] = "$(data["source_version"]["major"]).$(data["source_version"]["minor"]).$(data["source_version"]["patch"])"
    end
end


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


""
function parse_json(file::String; kwargs...)
    data = open(file) do io
        parse_json(io; filetype=split(lowercase(file), '.')[end], kwargs...)
    end
    return data
end


"Parses json from iostream or string"
function parse_json(io::IO; kwargs...)::Dict{String,Any}
    data = JSON.parse(io)

    _jsonver2juliaver!(data)

    _parse_mats!(data)

    if haskey(data, "files")
        data["files"] = Set(data["files"])
    end

    if get(kwargs, :validate, true)
        PowerModels.correct_network_data!(data)
    end

    return data
end


function print_file(path::String, data; kwargs...)
    open(path, "w") do io
        print_file(io, data; kwargs...)
    end
end


function print_file(io::IO, data; indent=false)
    if indent
        JSON.print(io, JSON.parse(sprint(show_json, PMDSerialization(), data)), 4)
    else
        JSON.print(io, JSON.parse(sprint(show_json, PMDSerialization(), data)))
    end
end


function show_json(io::StructuralContext, ::PMDSerialization, f::Matrix{<:Any})
    N, M = size(f)
    value = string("[", join([join([f[i,j] for j in 1:M], " ") for i in 1:N], "; "), "]")
    eltyp = isempty(f) ? eltype(f) : typeof(f[1,1])
        out = Dict(:type=>:Matrix, :eltype=>eltyp, :value=>value)
    return show_json(io, StandardSerialization(), out)
end


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
