
#import PowerModelsDistribution
#get = PowerModelsDistribution.get

function scale(dict, key, scale)
    if haskey(dict, key)
        dict[key] *= scale
    end
end


function add_virtual_get_id!(data_model, comp_type, comp)
    if !haskey(data_model, comp_type)
        data_model[comp_type] = Dict{Any, Any}()
    end
    comp_dict = data_model[comp_type]
    virtual_ids = [parse(Int, x[1]) for x in [match(r"_virtual_([1-9]{1}[0-9]*)", id) for id in keys(comp_dict) if isa(id, AbstractString)] if !isnothing(x)]
    if isempty(virtual_ids)
        id = "_virtual_1"
    else
        id = "_virtual_$(maximum(virtual_ids)+1)"
    end
    comp["id"] = id
    comp_dict[id] = comp
    return id
end

function delete_component!(data_model, comp_type, comp)
    delete!(data_model[comp_type], comp["id"])
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


function data_model_index!(data_model; components=["line", "shunt", "generator", "load", "transformer_2wa"])
    bus_id2ind = Dict()

    for (i, id) in enumerate(keys(data_model["bus"]))
        data_model["bus"][id]["index"] = i
        bus_id2ind[id] = i
    end
    data_model["bus"] = Dict{String, Any}(string(bus_id2ind[id])=>bus for (id, bus) in data_model["bus"])

    for comp_type in components
        comp_dict = Dict{String, Any}()
        for (i,(id,comp)) in enumerate(data_model[comp_type])
            @assert(!haskey(comp, "index"), "$comp_type $id: component already has an index.")
            comp["index"] = i
            comp["id"] = id
            comp_dict["$i"] = comp
            for bus_key in ["f_bus", "t_bus", "bus"]
                if haskey(comp, bus_key)
                    comp[bus_key] = bus_id2ind[comp[bus_key]]
                end
            end
        end
        data_model[comp_type] = comp_dict
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
