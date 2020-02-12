
#import PowerModelsDistribution
#get = PowerModelsDistribution.get

function scale(dict, key, scale)
    if haskey(dict, key)
        dict[key] *= scale
    end
end


function add_component!(data_model, comp_type, comp)
    if !haskey(data_model, comp_type)
        data_model[comp_type] = Dict{String, Any}()
    end
    comp_dict = data_model[comp_type]
    virtual_ids = [parse(Int, x[1]) for x in [match(r"_virtual_([1-9]{1}[0-9]*)", id) for id in keys(comp_dict)] if !isnothing(x)]
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

function add_mapping!(data_model::Dict{String, Any}, mapping_type::String, mapping::Dict{String, <:Any})
    if !haskey(data_model, "mapping")
        data_model["mapping"] = Dict{String, Any}()
    end

    if !haskey(data_model["mapping"], mapping_type)
        data_model["mapping"][mapping_type] = Array{Dict{String, Any}, 1}()
    end

    push!(data_model["mapping"][mapping_type], mapping)
end


function data_model_index!(data_model; components=["line", "bus", "shunt", "transformer_2w_ideal", "load", "generator"])
    for comp_type in components
        comp_dict = Dict{String, Any}()
        for (i,(id,comp)) in enumerate(data_model[comp_type])
            @assert(!haskey(comp, "index"), "$comp_type $id: component already has an index.")
            comp["index"] = i
            comp["id"] = id
            comp_dict["$i"] = comp
        end
        data_model[comp_type] = comp_dict
    end
    return data_model
end


function solution_ind2id(solution, data_model; id_prop="id")
    for comp_type in keys(solution)
        if isa(solution[comp_type], Dict)
            comp_dict = Dict{String, Any}()
            for (ind, comp) in solution[comp_type]
                id = data_model[comp_type][ind][id_prop]
                comp_dict[id] = comp
            end
            solution[comp_type] = comp_dict
        end
    end

    return solution
end
