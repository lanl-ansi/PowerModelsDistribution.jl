
"""
"""
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


"""
"""
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


"""
"""
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


"""
"""
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
