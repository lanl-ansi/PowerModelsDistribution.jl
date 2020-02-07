# Parse PowerModels data from JSON exports of PowerModels data structures.
# Customize JSON printing of matrices

import JSON.Serializations: CommonSerialization, StandardSerialization
import JSON.Writer: StructuralContext, show_json
struct PMDSerialization <: CommonSerialization end

function _jsonver2juliaver!(pm_data)
    if haskey(pm_data, "source_version") && isa(pm_data["source_version"], Dict)
        pm_data["source_version"] = "$(pm_data["source_version"]["major"]).$(pm_data["source_version"]["minor"]).$(pm_data["source_version"]["patch"])"
    end
end


function _parse_mats_recursive!(dict)
    parse =  haskey(dict, "type") && dict["type"]=="Matrix" && haskey(dict, "eltype") && haskey(dict, "value")
    if parse
        return _parse_matrix_value(dict["value"], dict["eltype"])
    else
        for (k,v) in dict
            if isa(v, Dict)
                dict[k] = _parse_mats!(v)
            elseif isa(v, Vector) && Dict <: eltype(v)
                dict[k] = [isa(x, Dict) ? _parse_mats!(x) : x for x in v]
            end
        end

        return dict
    end
end


function _parse_mats!(root_dict)
    stack = [(root_dict, k, v) for (k, v) in root_dict]
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
    pmd_data = open(file) do io
        parse_json(io; filetype=split(lowercase(file), '.')[end], kwargs...)
    end
    return pmd_data
end


"Parses json from iostream or string"
function parse_json(io::IO; kwargs...)::Dict{String,Any}
    pm_data = JSON.parse(io)

    _jsonver2juliaver!(pm_data)

    _parse_mats!(pm_data)

    if get(kwargs, :validate, true)
        PowerModels.correct_network_data!(pm_data)
    end

    return pm_data
end


function print_file(path::String, pmd_data)
    open(path, "w") do io
        print_file(io, pmd_data)
    end
end


function print_file(io::IO, pmd_data)
    JSON.print(io, JSON.parse(sprint(show_json, PMDSerialization(), pmd_data)))
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
