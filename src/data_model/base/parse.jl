
"""
"""
function Base.parse(::Type{T}, data::String) where {subtype, T <: Vector{subtype}}
    elements = _parse_dss_vector(data)
    nphases = length(elements)
    vector = zeros(subtype, nphases)
    _fill_vector!(vector, elements)

    return vector
end


"""
"""
function Base.parse(::Type{T}, data::String) where {subtype, T <: Matrix{subtype}}
    rows = _parse_dss_matrix(data)
    nphases = maximum([length(row) for row in rows])
    matrix = zeros(subtype, nphases, nphases)
    _fill_matrix!(matrix, rows)

    return matrix
end


Base.parse(::Type{T}, data::String) where T <: String = data


"""
"""
function Base.parse(::Type{T}, conn::String)::T where T <: ConnConfig
    if conn ∈ ["wye", "y", "ln"]
        return WYE
    elseif conn ∈ ["delta", "ll"]
        return DELTA
    end

    @warn "Unsupported connection $conn, defaulting to WYE"
    return WYE
end


"""
"""
function Base.parse(::Type{T}, model::String)::T where T <: LoadModel
    model = parse(Int, model)
    if model ∈ [3, 4, 7, 8]
        @warn "dss load model '$model' not supported; treating as constant POWER model"
        model = 1
    elseif model == 6
        @warn "dss load model '$model' identical to model '1' in current feature set; treating as constant POWER model"
        model = 1
    end

    return _dss2pmd_load_model[model]
end


"""
"""
function Base.parse(::Type{T}, model::String)::T where T <: CapControlType

    if isempty(model)
        return CAP_DISABLED
    elseif model == "kvar"
        return CAP_REACTIVE_POWER
    elseif model == "current"
        return CAP_CURRENT
    elseif model == "voltage"
        return CAP_VOLTAGE
    elseif model == "time"
        return CAP_TIME
    else
        @warn "cap control type '$(model)' unrecognized, returning CAP_DISABLED"

        return CAP_DISABLED
    end
end


"""
"""
function Base.parse(::Type{T}, status::String)::T where T <: Status
    if status ∈ ["y", "yes", "true"]
        return ENABLED
    elseif status ∈ ["n", "no", "false"]
        return DISABLED
    end

    @warn "enabled code '$status' not recognized, defaulting to ENABLED"
    return ENABLED
end
