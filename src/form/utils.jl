"Merges flow variables that enter the same terminals, i.e. multiple neutrals of an underground cable connected to same neutral terminal"
function _merge_bus_flows(pm::ExplicitNeutralModels, flows::Vector, connections::Vector)::JuMP.Containers.DenseAxisArray
    flows_merged = []
    conns_unique = unique(connections)
    for t in conns_unique
        idxs = findall(connections.==t)
        flows_t = flows[idxs]
        if length(flows_t)==1
            flows_merged_t = flows_t[1]
        else
            flows_merged_t = sum(flows_t)
        end
        push!(flows_merged, flows_merged_t)
    end
    JuMP.Containers.DenseAxisArray(flows_merged, conns_unique)
end
