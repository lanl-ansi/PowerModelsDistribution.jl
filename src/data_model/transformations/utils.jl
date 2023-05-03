"""
"""
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


"""
"""
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



"""
"""
function _register_awaiting_ground!(eng::EngineeringDataModel, bus::EngBus, connections::Vector{Int})
    if !haskey(eng.metadata["awaiting_ground"], bus.name)
        eng.metadata.awaiting_ground[bus.name] = Vector{Int}[]
    end

    push!(eng.metadata.awaiting_ground[bus.name], connections)
end


"""
"""
function _parse_dss_pairs(property_pairs::Vector{Pair{String,String}}, dss_obj_type::String)::Dict{Symbol,Any}
    dss_struct = getfield(PowerModelsDistribution, Symbol("Dss$(titlecase(dss_obj_type))"))
    dtypes = Dict{String,Type}(string(k) => v for (k,v) in zip(fieldnames(dss_struct), fieldtypes(dss_struct)))

    Dict{Symbol,Any}(Symbol(k) => parse(dtypes[k], v) for (k,v) in property_pairs)
end


""
function _correct_ptphase!(controller::EngTransformerControls, transformer::EngTransformer)
    controller["terminals"][1][1] = transformer["connections"][1][1]
end


"""
"""
function apply_linecode!(line::EngLine, linecode::EngLinecode)::EngLine
    for pn in [:rs, :xs, :g_fr, :g_to, :b_fr, :b_to, :cm_ub, :sm_ub]
        if ismissing(line[pn]) && !ismissing(linecode[pn])
            line[pn] = linecode[pn]
        end
    end

    return line
end


"""
"""
function apply_linecodes!(eng::EngineeringModel{T})::EngineeringModel{T} where T <: NetworkModel
    for line in eng.line
        if !ismissing(line.linecode)
            apply_linecode!(line, eng.linecode[line.linecode])
        end
    end

    return eng
end
