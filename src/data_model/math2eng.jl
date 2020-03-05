""
function _map_math2eng!(data_math)
    @assert get(data_math, "data_model", "mathematical") == "mathematical" "Cannot map data to engineering model: provided data is not a mathematical model"
    @assert haskey(data_math, "map") "Cannot map data to engineering model: no mapping from mathematical to engineering data model is provided"

    data_eng = Dict{<:Any,<:Any}()

    map_keys = sort(keys(data_math["map"]); reverse=true)
    for map in map_keys
        # TODO
    end

end


# MAP SOLUTION UP
""
function solution_unmap!(solution::Dict, data_model::Dict)
    for (name, data) in reverse(data_model["mappings"])
        if name=="decompose_transformer_nw"
            for bus_id in values(data["vbuses"])
                delete!(solution["bus"], bus_id)
            end

            for line_id in values(data["vlines"])
                delete!(solution["branch"], line_id)
            end

            pt = [solution["transformer"][tr_id]["pf"] for tr_id in data["trans_2wa"]]
            qt = [solution["transformer"][tr_id]["qf"] for tr_id in data["trans_2wa"]]
            for tr_id in data["trans_2wa"]
                delete!(solution["transformer"], tr_id)
            end

            add_solution!(solution, "transformer_nw", data["trans"]["id"], Dict("pt"=>pt, "qt"=>qt))
        elseif name=="capacitor_to_shunt"
            # shunt has no solutions defined
            delete_solution!(solution, "shunt", data["shunt_id"])
            add_solution!(solution, "capacitor", data["capacitor"]["id"], Dict())
        elseif name=="load_to_shunt"
            # shunt has no solutions, but a load should have!
            delete!(solution, "shunt", data["shunt_id"])
            add_solution!(solution, "load", data["load"]["id"], Dict())
        elseif name=="decompose_voltage_source"
            gen = solution["gen"][data["gen_id"]]
            delete_solution!(solution, "gen", data["gen_id"])
            add_solution!(solution, "voltage_source", data["voltage_source"]["id"], Dict("pg"=>gen["pg"], "qg"=>gen["qg"]))

        end
    end

    # remove component dicts if empty
    for (comp_type, comp_dict) in solution
        if isa(comp_dict, Dict) && isempty(comp_dict)
            delete!(solution, comp_type)
        end
    end
end

