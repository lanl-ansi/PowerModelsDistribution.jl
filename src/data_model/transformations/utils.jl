"helper function to determine if `property1` appears after `property2` in the property pair list"
function _is_after(prop_order::Vector{Pair{String,String}}, property1::String, property2::String)::Bool
    property1_idx = 0
    property2_idx = 0

    for (i, (prop,_)) in enumerate(prop_order)
        if prop == property1
            property1_idx = i
        elseif prop == property2
            property2_idx = i
        end
    end

    return property1_idx > property2_idx
end


"helper function to determine if `property1` appears after `property2` in the property pair list for objects that have windings"
function _is_after(prop_order::Vector{Pair{String,String}}, property1::String, property2::String, winding::Int)::Bool
    property1_idx = 0
    property2_idx = 0
    wdg_idx = 0

    for (i, (prop,value)) in enumerate(prop_order)
        if prop == "wdg" && parse(Int, value) == winding
            wdg_idx = i
        end

        if prop == property1
            property1_idx = i
        elseif prop == property2
            property2_idx = i
        end
    end

    return property1_idx > property2_idx && property1_idx > wdg_idx
end


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


"Helper function to infer a full property name from a partial one"
function _infer_partial_property_name(pn::AbstractString, dss_obj::T)::String where T <: DssObject
    try
        pn = _dss_short_prop_names_map[typeof(dss_obj)][pn]
    catch KeyError
    end

    return pn
end
