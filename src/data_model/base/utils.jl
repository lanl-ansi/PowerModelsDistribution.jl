"Helper function for Base.parse for dss matrices"
function _fill_matrix!(matrix::Matrix{T}, rows::Vector{Vector{String}})::Nothing where T
    nphases = size(matrix)[1]

    if length(rows) == 1
        for i in 1:nphases
            matrix[i, i] = parse(T, rows[1][1])
        end
    elseif all([length(row) for row in rows] .== [i for i in 1:nphases])
        for (i, row) in enumerate(rows)
            for (j, col) in enumerate(row)
                matrix[i, j] = matrix[j, i] = parse(T, col)
            end
        end
    elseif all([length(row) for row in rows] .== nphases)
        for (i, row) in enumerate(rows)
            for (j, col) in enumerate(row)
                matrix[i, j] = parse(T, col)
            end
        end
    end

    nothing
end


"Helper function for Base.parse for dss matrices"
function _parse_dss_matrix(data::String)::Vector{Vector{String}}
    rows = Vector{String}[]
    for line in split(strip(strip(data, _array_delimiters)), '|')
        cols = String[]
        for item in split(strip(line), _dss_matrix_regex)
            push!(cols, string(item))
        end
        push!(rows, cols)
    end

    return rows
end


"Helper function for Base.parse for dss vectors"
function _parse_dss_vector(data::String)::Vector{String}
    if _isa_rpn(data)
        matches = collect((m.match for m = eachmatch(_dss_rpn_array_sep_regex, data, overlap=false)))
        if length(matches) == 2
            return data
        else
            return string.(strip.(getfield.(collect(eachmatch(_dss_cmd_new_regex, data)), :match)))
        end
    else

        for delim in _array_delimiters
            data = replace(data, delim => "")
        end

        return filter(x->!isempty(x), strip.(split(strip(data), _dss_array_regex)))
    end
end


"Helper function for Base.parse for dss vectors"
function _fill_vector!(vector::Vector{T}, elements::Vector{String})::Nothing where T
    for (i, element) in enumerate(elements)
        if _isa_rpn(element)
            vector[i] = _parse_rpn(T, element)
        else
            vector[i] = parse(T, element)
        end
    end

    nothing
end


"Helper function to convert InfrastructureModel,InfrastructureObject into Dict{String,Any}"
function _convert_model_to_dict(data::Union{InfrastructureModel,InfrastructureObject})::Dict{String,Any}
    out = Dict{String,Any}(
    )

    for property in propertynames(data)
        item = getproperty(data, property)

        if isa(item, Dict)
            out["$property"] = Dict{String,Any}()
            for (id, obj) in item
                if isa(obj, InfrastructureObject)
                    out["$property"]["$id"] = _convert_model_to_dict(obj)
                else
                    out["$property"]["$id"] = obj
                end
            end
        elseif isa(item, InfrastructureObject)
            out["$property"] = _convert_model_to_dict(item)
        else
            out["$property"] = item
        end
    end

    return filter(x->!ismissing(x.second), out)
end


"Helper function to pull the specified properties from dss property pairs"
function _get_raw_fields(property_pairs::Vector{Pair{String,String}})::Vector{Symbol}
    collect(Symbol(x.first) for x in property_pairs)
end
