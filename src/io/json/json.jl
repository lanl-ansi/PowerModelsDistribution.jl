"""
    parse_json(file::String)

parses json files that were dumped via JSON.print (or PMD.print_file)
"""
function parse_json(file::String)
    open(file, "r") do io
        parse_json(io)
    end
end


"""
    parse_json(io::IO)

parses json files that were dumped via JSON.print (or PMD.print_file)
"""
function parse_json(io::IO)
    data = JSON.parse(io)
    correct_json_import!(data)

    return data
end


"""
    correct_json_import!(data::Dict{String,<:Any})

helper function to correct data imported from json
"""
function correct_json_import!(data::Dict{String,<:Any})
    _fix_dtypes!(data)
end


"recursive function to fix data types from data imported from json"
function _fix_dtypes!(data::Dict)
    for (k, v) in data
        if isa(v, Dict)
            _fix_dtypes!(v)
        else
            _fix_enums!(data, k, data[k])
            _fix_arrays!(data, k, data[k])
            _fix_nulls!(data, k, data[k])
        end
    end
end


"helper function to fix matrices (from vector of vectors) and vector dtypes"
function _fix_arrays!(obj, prop, val)
    if isa(val, Vector)
        dtypes = unique(Vector{Any}([typeof(v) for v in val]))
        if isempty(val)
        elseif all([isa(v, Vector) for v in val])
            obj[prop] = hcat(obj[prop]...)
        else
            if length(dtypes) == 1
                dtype = dtypes[1]
                obj[prop] = Vector{dtype}(val)
            end
        end
    end
end


"helper function to convert stringified enums"
function _fix_enums!(obj, prop, val)
    if isa(val, String) && uppercase(val) == val && Symbol(val) in names(PowerModelsDistribution)
        obj[prop] = getfield(PowerModelsDistribution, Symbol(val))
    end
end


"helper function to fix null values from json (usually Inf or NaN)"
function _fix_nulls!(obj, prop, val)
    if endswith(prop, "_ub") || endswith(prop, "max")
        fill_val = Inf
    elseif endswith(prop, "_lb") || endswith(prop, "min")
        fill_val = -Inf
    else
        fill_val = NaN
    end

    if isa(val, Matrix) && any(val .=== nothing)
        @debug "a 'null' was encountered in the json import, making an assumption that null values in $prop = $fill_val"
        valdtype = valtype(val)
        if isa(valdtype, Union)
            dtype = [getproperty(valdtype, n) for n in propertynames(valdtype) if getproperty(valdtype, n) != Nothing][end]
        else
            dtype = valdtype
        end
        val[val .=== nothing] .= fill_val
        obj[prop] = Matrix{valtype(val) == Nothing ? typeof(fill_val) : valtype(val)}(val)
    elseif isa(val, Vector) && any(v === nothing for v in val)
        @debug "a 'null' was encountered in the json import, making an assumption that null values in $prop = $fill_val"
        obj[prop] = Vector{valtype(val) == Nothing ? typeof(fill_val) : valtype(val)}([v === nothing ? fill_val : v for v in val])
    elseif val === nothing
        @debug "a 'null' was encountered in the json import, making an assumption that null values in $prop = $fill_val"
        obj[prop] = fill_val
    end
end


"""
    print_file(path::String, data::Dict{String,<:Any}; indent::Int=2)

prints a PowerModelsDistribution data structure into a JSON file
"""
function print_file(path::String, data::Dict{String,<:Any}; indent::Int=2)
    open(path, "w") do io
        print_file(io, data; indent=indent)
    end
end


"""
    print_file(io::IO, data::Dict{String,<:Any}; indent::Int=2)

prints a PowerModelsDistribution data structure into a JSON file
"""
function print_file(io::IO, data::Dict{String,<:Any}; indent::Int=2)
    JSON.print(io, data, indent)
end
