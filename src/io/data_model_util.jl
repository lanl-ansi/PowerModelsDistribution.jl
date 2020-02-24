
#import PowerModelsDistribution
#get = PowerModelsDistribution.get

function scale(dict, key, scale)
    if haskey(dict, key)
        dict[key] *= scale
    end
end


function add_virtual!(data_model, comp_type, comp)
    if !haskey(data_model, comp_type)
        data_model[comp_type] = Dict{Any, Any}()
    end
    comp_dict = data_model[comp_type]
    virtual_ids = [parse(Int, x[1]) for x in [match(r"_virtual_([1-9]{1}[0-9]*)", id) for id in keys(comp_dict) if isa(id, AbstractString)] if x!=nothing]
    if isempty(virtual_ids)
        id = "_virtual_1"
    else
        id = "_virtual_$(maximum(virtual_ids)+1)"
    end
    comp["id"] = id
    comp_dict[id] = comp
    return comp
end

add_virtual_get_id!(data_model, comp_type, comp) = add_virtual!(data_model, comp_type, comp)["id"]

function delete_component!(data_model, comp_type, comp::Dict)
    delete!(data_model[comp_type], comp["id"])
    if isempty(data_model[comp_type])
        delete!(data_model, comp_type)
    end
end

function delete_component!(data_model, comp_type, id::Any)
    delete!(data_model[comp_type], id)
    if isempty(data_model[comp_type])
        delete!(data_model, comp_type)
    end
end

function add_mappings!(data_model::Dict{String, Any}, mapping_type::String, mappings::Vector)
    if !haskey(data_model, "mappings")
        data_model["mappings"] = []
    end

    append!(data_model["mappings"], [(mapping_type, mapping) for mapping in mappings])
end


function _get_next_index(last_index, presets)
    new_index = last_index+1
    while new_index in presets
        new_index += 1
    end
    return new_index
end


function data_model_index!(data_model; components=["bus", "line", "shunt", "generator", "load", "transformer_2wa"], index_presets=Dict())
    comp_id2ind = Dict()

    # bus should be the first component, because we want to
    for comp_type in components
        comp_dict = Dict{String, Any}()

        if !haskey(index_presets, comp_type)
            for (i,(id,comp)) in enumerate(data_model[comp_type])
                @assert(!haskey(comp, "index"), "$comp_type $id: component already has an index.")
                comp["index"] = i
                comp["id"] = id
                comp_dict["$i"] = comp
            end
        else
            last_index = 0

            for (id, comp) in data_model[comp_type]
                @assert(!haskey(comp, "index"), "$comp_type $id: component already has an index.")
                if haskey(index_presets[comp_type], id)
                    println("$comp_type: $id found, ->$(index_presets[comp_type][id])")
                    comp["index"] = index_presets[comp_type][id]
                else
                    comp["index"] = _get_next_index(last_index, values(index_presets[comp_type]))
                    last_index = comp["index"]
                end

                comp["id"] = id
                comp_dict["$(comp["index"])"] = comp
            end
        end

        data_model[comp_type] = comp_dict
        comp_id2ind[comp_type] = Dict(comp["id"]=>comp["index"] for comp in values(comp_dict))
    end

    # update bus references
    for comp_type in components
        for (_, comp) in data_model[comp_type]
            for bus_key in ["f_bus", "t_bus", "bus"]
                if haskey(comp, bus_key)
                    comp[bus_key] = comp_id2ind["bus"][comp[bus_key]]
                end
            end
        end
    end

    return data_model
end


function solution_identify!(solution, data_model; id_prop="id")
    for comp_type in keys(solution)
        if isa(solution[comp_type], Dict)
            comp_dict = Dict{Any, Any}()
            for (ind, comp) in solution[comp_type]
                id = data_model[comp_type][ind][id_prop]
                comp_dict[id] = comp
            end
            solution[comp_type] = comp_dict
        end
    end

    return solution
end

function add_solution!(solution, comp_type, id, data)
    if !haskey(solution, comp_type)
        solution[comp_type] = Dict()
    end

    if !haskey(solution[comp_type], id)
        solution[comp_type][id] = Dict{String, Any}()
    end

    for (key, prop) in data
        solution[comp_type][id][key] = prop
    end
end


function delete_solution!(solution, comp_type, id, props)
    if haskey(solution, comp_type)
        if haskey(solution[comp_type], id)
            for prop in props
                delete!(solution[comp_type][id], prop)
            end
        end
    end
end


function delete_solution!(solution, comp_type, id)
    if haskey(solution, comp_type)
        if haskey(solution[comp_type], id)
            delete!(solution[comp_type], id)
        end
    end
end

##
function _get_new_ground(terminals)
    if isa(terminals, Vector{Int})
        return maximum(terminals)+1
    else
        nrs = [parse(Int, x[1]) for x in [match(r"n([1-9]{1}[0-9]*)", string(t)) for t in terminals] if x!=nothing]
        new = isempty(nrs) ? 1 : maximum(nrs)+1
        if isa(terminals, Vector{Symbol})
            return Symbol("g$new")
        else
            return "g$new"
        end
    end
end


function _get_ground!(bus)
    # find perfect groundings (true ground)
    grounded_perfect = []
    for i in 1:length(bus["grounded"])
        if bus["rg"][i]==0 && bus["xg"][i]==0
            push!(grounded_perfect, bus["grounded"][i])
        end
    end

    if !isempty(grounded_perfect)
        return grounded_perfect[1]
    else
        g = _get_new_ground(bus["terminals"])
        push!(bus["terminals"], g)
        push!(bus["rg"], 0.0)
        push!(bus["xg"], 0.0)
        return g
    end
end
