"returns number of phases implied by a two-bus (edge) object"
function _get_implied_nphases(bus1::AbstractString, bus2::AbstractString; default::Int=3)
    f_conds = _get_conductors_ordered(bus1; default=collect(1:default), check_length=false)
    t_conds = _get_conductors_ordered(bus2; default=collect(1:default), check_length=false)

    if !isempty(f_conds) || !isempty(t_conds)
        return maximum([length(f_conds), length(t_conds)])
    else
        return default
    end
end


"returns number of phases implied by a transformer object"
function _get_implied_nphases(buses::Vector{<:AbstractString}; default::Int=3)
    nphases = Int[]
    for bus in buses
        conds = _get_conductors_ordered(bus; default=collect(1:default), check_length=false)
        if !isempty(conds)
            push!(nphases, length(conds))
        end
    end

    if !isempty(nphases)
        return maximum(nphases)
    else
        return default
    end
end


"returns number of phases implied by a single-bus (node) object"
function _get_implied_nphases(bus1::AbstractString; default::Int=3)::Int
    conds = _get_conductors_ordered(bus1; default=collect(1:default), check_length=false)

    if !isempty(conds)
        return length(conds)
    else
        return default
    end
end


""
function _infer_partial_property_name(pn::AbstractString, dss_obj::T)::String where T <: DssObject
    try
        pn = _dss_short_prop_names_map[typeof(dss_obj)][pn]
    catch KeyError
    end

    return pn
end
