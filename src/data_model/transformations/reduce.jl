"""
    reduce_line_series!(eng::Dict{String,<:Any}; remove_original_lines::Bool=false)::Dict{String,<:Any}

This is a function to merge series of lines which only connect to buses with no other connections
(i.e., string of buses with no loads, generators, transformers, etc.). This function will preserve
the total length of the merged lines.

If `remove_original_lines`, the original lines and eliminated buses will be deleted from the data structure,
otherwise the lines and buses will be `DISABLED`.
"""
function reduce_line_series!(eng::Dict{String,<:Any}; remove_original_lines::Bool=false)
    G = Graphs.Graph(sum([length(get(eng, k, Dict())) for k in ["bus", PowerModelsDistribution._eng_node_elements...]]))

    id_map = Dict((k,i) => n for (n,(k,i)) in enumerate([(k,i) for k in ["bus", PowerModelsDistribution._eng_node_elements...] for i in keys(get(eng, k, Dict()))]))

    for k in PowerModelsDistribution._eng_node_elements
        for (i,obj) in get(eng, k, Dict())
            Graphs.add_edge!(G, id_map[(k,i)], id_map[("bus",obj["bus"])])
        end
    end

    lines_by_id = Set()
    line_lookup_by_bus = Dict(i => [] for i in keys(eng["bus"]))
    for k in PowerModelsDistribution._eng_edge_elements
        for (i,obj) in get(eng, k, Dict())
            if k == "transformer" && haskey(obj, "bus")
                for bus1 in obj["bus"]
                    for bus2 in obj["bus"]
                        if bus1 != bus2
                            Graphs.add_edge!(G, id_map[("bus", bus1)], id_map[("bus", bus2)])
                        end
                    end
                end
            else
                Graphs.add_edge!(G, id_map[("bus", obj["f_bus"])], id_map[("bus", obj["t_bus"])])
                if k == "line"
                    push!(lines_by_id, Set([id_map[("bus", obj["f_bus"])], id_map[("bus", obj["t_bus"])]]))
                    push!(line_lookup_by_bus[obj["f_bus"]], i)
                    push!(line_lookup_by_bus[obj["t_bus"]], i)
                end
            end
        end
    end

    id_lookup = Dict(n => (k,i) for ((k,i), n) in id_map)

    b2r = Set()
    for node in Graphs.vertices(G)
        n = Graphs.all_neighbors(G, node)
        if length(n) == 2 && all([id_lookup[i][1] == "bus" for i in n]) && all([Set([node, i]) in lines_by_id for i in n])
            push!(b2r, node)
        end
    end
    buses2remove = Set([id_lookup[n][2] for n in b2r])
    lines2merge = [Set(line_lookup_by_bus[n]) for n in buses2remove if length(line_lookup_by_bus[n])==2]

    while true
        to_combine = []
        for (i, l1) in enumerate(lines2merge)
            for (j, l2) in enumerate(lines2merge)
                if i<j && !isempty(intersect(l1, l2))
                    push!(to_combine, (i,j))
                end
            end
        end

        if isempty(to_combine)
            break
        else
            touched = []
            merged = []
            for (i,j) in to_combine
                if (j ∉ touched) && (i ∉ touched)
                    lines2merge[i] = union(lines2merge[i], lines2merge[j])
                    push!(touched, j)
                    push!(touched, i)
                    push!(merged, j)
                end
            end

            lines2merge = [lines2merge[i] for i in 1:length(lines2merge) if (i ∉ merged)]
        end
    end

    for line_group in lines2merge
        if length(Set([Set(eng["line"][line_id][connection_type]) for line_id in line_group for connection_type in ["f_connections", "t_connections"]])) == 1
            buses = [eng["line"][line_id][bus_end] for line_id in line_group for bus_end in ["f_bus", "t_bus"]]
            buses = Dict(b => count(isequal(b), buses) for b in buses)
            ft_buses = [k for (k,v) in buses if v == 1]
            new_line = Dict{String,Any}(
                "f_bus" => ft_buses[1],
                "t_bus" => ft_buses[2],
                "f_connections" => eng["line"][first(line_group)]["f_connections"],
                "t_connections" => eng["line"][first(line_group)]["t_connections"],
                "source_id" => "line._virtual.$(join(line_group, "+"))",
                "status" => all(eng["line"][line_id]["status"] == ENABLED for line_id in line_group) ? ENABLED : DISABLED
            )
        end

        if length(Set([eng["line"][line_id]["linecode"] for line_id in line_group])) == 1
            new_line["length"] = sum(eng["line"][line_id]["length"] for line_id in line_group)
            new_line["linecode"] = eng["line"][first(line_group)]["linecode"]
        else
            for line_id in line_group
                _apply_linecode!(eng["line"][line_id], eng)
            end

            new_line["length"] = sum(eng["line"][line_id]["length"] for line_id in line_group)

            for prop_id in ["rs", "xs"]
                new_line[prop_id] = sum(eng["line"][line_id][prop_id]*eng["line"][line_id]["length"] for line_id in line_group) / new_line["length"]
            end

            for prop_id in ["b_fr", "b_to", "g_fr", "g_to"]
                new_line[prop_id] = LinearAlgebra.pinv(sum(LinearAlgebra.pinv(eng["line"][line_id][prop_id]) for line_id in line_group)) / new_line["length"]
            end
        end

        eng["line"]["_virtual.$(join(line_group, "+"))"] = new_line

        for line_id in line_group
            if remove_original_lines
                for bus in keys(buses)
                    if bus ∉ ft_buses
                        delete!(eng["bus"], bus)
                    end
                end
                delete!(eng["line"], line_id)
            else
                for bus in keys(buses)
                    if bus ∉ ft_buses
                        eng["bus"][bus]["status"] = DISABLED
                    end
                end

                eng["line"][line_id]["status"] = DISABLED
            end
        end
    end

    return eng["line"]
end
