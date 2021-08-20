"Merges flow variables that enter the same terminals, i.e. multiple neutrals of an underground cable connected to same neutral terminal"
function _merge_bus_flows(pm::ExplicitNeutralModels, flows::Vector, connections::Vector)::JuMP.Containers.DenseAxisArray
    flows_merged = []
    conns_unique = unique(connections)
    for t in conns_unique
        idxs = findall(connections.==t)
        flows_t = flows[idxs]
        if length(flows_t)==1
            flows_merged_t = flows_t[1]
        elseif any(isa(a, JuMP.NonlinearExpression) for a in flows_t)
            flows_merged_t = JuMP.@NLexpression(pm.model, sum(flows_t[i] for i in 1:length(flows_t)))
        else
            flows_merged_t = sum(flows_t)
        end
        push!(flows_merged, flows_merged_t)
    end
    JuMP.Containers.DenseAxisArray(flows_merged, conns_unique)
end
