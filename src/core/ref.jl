"""
    apply_pmd!(func!::Function, ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any}; apply_to_subnetworks::Bool=true)

PowerModelsDistribution wrapper for the InfrastructureModels `apply!` function
"""
function apply_pmd!(func!::Function, ref::Dict{Symbol, <:Any}, data::Dict{String, <:Any}; apply_to_subnetworks::Bool = true)
    _IM.apply!(func!, ref, data, pmd_it_sym; apply_to_subnetworks = apply_to_subnetworks)
end


"""
    ref_calc_storage_injection_bounds(storage, buses)

Computes storage bounds
"""
function ref_calc_storage_injection_bounds(storage, buses)
    injection_lb = Dict()
    injection_ub = Dict()

    for (i, strg) in storage
        connections = strg["connections"]
        ncnds = length(connections)
        injection_lb[i] = fill(-Inf, ncnds)
        injection_ub[i] = fill( Inf, ncnds)

        if haskey(strg, "thermal_rating")
            injection_lb[i] = max.(injection_lb[i], -strg["thermal_rating"])
            injection_ub[i] = min.(injection_ub[i],  strg["thermal_rating"])
        end

        if haskey(strg, "current_rating")
            for (j, t) in connections
                vmax = buses[strg["storage_bus"]]["vmax"][findfirst(isequal(t), buses[strg["storage_bus"]]["terminals"])]

                injection_lb[i][j] = max(injection_lb[i][j], -strg["current_rating"][j]*vmax)
                injection_ub[i][j] = min(injection_ub[i][j],  strg["current_rating"][j]*vmax)
            end
        end
    end

    return injection_lb, injection_ub
end


"""
    calc_buspair_parameters(buses, branches)

Computes the buspair parameters for the network
"""
function calc_buspair_parameters(buses, branches)
    bus_lookup = Dict(bus["index"] => bus for (i,bus) in buses if bus["bus_type"] != 4)
    branch_lookup = Dict(branch["index"] => branch for (i,branch) in branches if branch["br_status"] == 1 && haskey(bus_lookup, branch["f_bus"]) && haskey(bus_lookup, branch["t_bus"]))

    buspair_indexes = Set((branch["f_bus"], branch["t_bus"], fc, tc) for (i,branch) in branch_lookup for (fc, tc) in zip(branch["f_connections"], branch["t_connections"]))

    bp_branch = Dict((bp, typemax(Int)) for bp in buspair_indexes)

    bp_angmin = Dict((bp, -Inf) for bp in buspair_indexes)
    bp_angmax = Dict((bp,  Inf) for bp in buspair_indexes)

    for (l,branch) in branch_lookup
        i = branch["f_bus"]
        j = branch["t_bus"]

        f_connections = branch["f_connections"]
        t_connections = branch["t_connections"]

        for (idx, (fc, tc)) in enumerate(zip(f_connections, t_connections))
            bp_angmin[(i,j,fc,tc)] = max(bp_angmin[(i,j,fc,tc)], branch["angmin"][idx])
            bp_angmax[(i,j,fc,tc)] = min(bp_angmax[(i,j,fc,tc)], branch["angmax"][idx])

            bp_branch[(i,j,fc,tc)] = min(bp_branch[i,j,fc,tc], l)
        end
    end

    buspairs = Dict(
        (i,j,fc,tc) => Dict(
            "branch" => bp_branch[(i,j,fc,tc)],
            "angmin" => bp_angmin[(i,j,fc,tc)],
            "angmax" => bp_angmax[(i,j,fc,tc)],
            "vm_fr_min" => bus_lookup[i]["vmin"][findfirst(isequal(fc), bus_lookup[i]["terminals"])],
            "vm_fr_max" => bus_lookup[i]["vmax"][findfirst(isequal(fc), bus_lookup[i]["terminals"])],
            "vm_to_min" => bus_lookup[j]["vmax"][findfirst(isequal(tc), bus_lookup[j]["terminals"])],
            "vm_to_max" => bus_lookup[j]["vmax"][findfirst(isequal(tc), bus_lookup[j]["terminals"])]
        ) for (i, j, fc, tc) in buspair_indexes
    )

    # add optional parameters
    for bp in buspair_indexes
        branch = branch_lookup[bp_branch[bp]]
        if haskey(branch, "rate_a")
            buspairs[bp]["rate_a"] = branch["rate_a"][findfirst(isequal(bp[end]), branch["t_connections"])]
        end
        if haskey(branch, "c_rating_a")
            buspairs[bp]["c_rating_a"] = branch["c_rating_a"][findfirst(isequal(bp[end]), branch["t_connections"])]
        end
    end
    return buspairs
end


"""
    find_conductor_ids!(data::Dict{String,Any})

Finds all conductor ids and puts a list of them under "conductor_ids" at the root level
"""
function find_conductor_ids!(data::Dict{String,<:Any})
    conductor_ids = []

    for (_,bus) in get(data, "bus", Dict())
        for t in get(bus, "terminals", [])
            if !(t in conductor_ids)
                push!(conductor_ids, t)
            end
        end
    end

    data["conductor_ids"] = [c for c in sort(conductor_ids)]
end


"""
    ref_add_core!(ref::Dict{Symbol,Any})

Returns a dict that stores commonly used pre-computed data from of the data dictionary,
primarily for converting data-types, filtering out deactivated components, and storing
system-wide values that need to be computed globally.
Some of the common keys include:
* `:off_angmin` and `:off_angmax` (see `calc_theta_delta_bounds(data)`),
* `:bus` -- the set `{(i, bus) in ref[:bus] : bus["bus_type"] != 4}`,
* `:gen` -- the set `{(i, gen) in ref[:gen] : gen["gen_status"] == 1 && gen["gen_bus"] in keys(ref[:bus])}`,
* `:branch` -- the set of branches that are active in the network (based on the component status values),
* `:arcs_branch_from` -- the set `[(i,b["f_bus"],b["t_bus"]) for (i,b) in ref[:branch]]`,
* `:arcs_branch_to` -- the set `[(i,b["t_bus"],b["f_bus"]) for (i,b) in ref[:branch]]`,
* `:arcs_branch` -- the set of arcs from both `arcs_from` and `arcs_to`,
* `:arcs_switch_from` -- the set `[(i,b["f_bus"],b["t_bus"]) for (i,b) in ref[:switch]]`,
* `:arcs_switch_to` -- the set `[(i,b["t_bus"],b["f_bus"]) for (i,b) in ref[:switch]]`,
* `:arcs_switch` -- the set of arcs from both `arcs_switch_from` and `arcs_switch_to`,
* `:arcs_transformer_from` -- the set `[(i,b["f_bus"],b["t_bus"]) for (i,b) in ref[:transformer]]`,
* `:arcs_transformer_to` -- the set `[(i,b["t_bus"],b["f_bus"]) for (i,b) in ref[:transformer]]`,
* `:arcs_transformer` -- the set of arcs from both `arcs_transformer_from` and `arcs_transformer_to`,
* `:bus_arcs_branch` -- the mapping `Dict(i => [(l,i,j) for (l,i,j) in ref[:arcs_branch]])`,
* `:bus_arcs_transformer` -- the mapping `Dict(i => [(l,i,j) for (l,i,j) in ref[:arcs_transformer]])`,
* `:bus_arcs_switch` -- the mapping `Dict(i => [(l,i,j) for (l,i,j) in ref[:arcs_switch]])`,
* `:buspairs` -- (see `buspair_parameters(ref[:arcs_branch_from], ref[:branch], ref[:bus])`),
* `:bus_gens` -- the mapping `Dict(i => [gen["gen_bus"] for (i,gen) in ref[:gen]])`.
* `:bus_loads` -- the mapping `Dict(i => [load["load_bus"] for (i,load) in ref[:load]])`.
* `:bus_shunts` -- the mapping `Dict(i => [shunt["shunt_bus"] for (i,shunt) in ref[:shunt]])`.
"""
function ref_add_core!(ref::Dict{Symbol,Any})
    for (nw, nw_ref) in ref[:it][pmd_it_sym][:nw]
        ### filter out inactive components ###

        nw_ref[:bus] = Dict(x for x in nw_ref[:bus] if (x.second["bus_type"] != pmd_math_component_status_inactive["bus"]))

        for (type, status) in pmd_math_component_status
            if type in _math_node_elements
                nw_ref[Symbol(type)] = Dict(x for x in get(nw_ref, Symbol(type), Dict()) if (x.second[status] != pmd_math_component_status_inactive[type] && x.second["$(type)_bus"] in keys(nw_ref[:bus])))
            elseif type in _math_edge_elements
                nw_ref[Symbol(type)] = Dict(x for x in get(nw_ref, Symbol(type), Dict()) if (x.second[status] != pmd_math_component_status_inactive[type] && x.second["f_bus"] in keys(nw_ref[:bus]) && x.second["t_bus"] in keys(nw_ref[:bus])))
            end
        end

        # dispatchable types ref
        for type in _math_dispatchable_elements
            nw_ref[Symbol("$(type)_dispatchable")] = filter(x->get(x.second, "dispatchable", 0) == 1, nw_ref[Symbol(type)])
        end

        ### setup arcs from edges ###
        nw_ref[:arcs_branch_from] = [(i,branch["f_bus"],branch["t_bus"]) for (i,branch) in nw_ref[:branch]]
        nw_ref[:arcs_branch_to]   = [(i,branch["t_bus"],branch["f_bus"]) for (i,branch) in nw_ref[:branch]]
        nw_ref[:arcs_branch] = [nw_ref[:arcs_branch_from]; nw_ref[:arcs_branch_to]]

        nw_ref[:arcs_switch_from] = [(i,switch["f_bus"],switch["t_bus"]) for (i,switch) in nw_ref[:switch]]
        nw_ref[:arcs_switch_to]   = [(i,switch["t_bus"],switch["f_bus"]) for (i,switch) in nw_ref[:switch]]
        nw_ref[:arcs_switch] = [nw_ref[:arcs_switch_from]; nw_ref[:arcs_switch_to]]

        nw_ref[:arcs_transformer_from] = [(i, transformer["f_bus"], transformer["t_bus"]) for (i,transformer) in nw_ref[:transformer]]
        nw_ref[:arcs_transformer_to] = [(i, transformer["t_bus"], transformer["f_bus"]) for (i,transformer) in nw_ref[:transformer]]
        nw_ref[:arcs_transformer] = [nw_ref[:arcs_transformer_from]; nw_ref[:arcs_transformer_to]]

        ### bus connected component lookups ###
        for type in _math_node_elements
            bus_objs = Dict((i, Int[]) for (i,bus) in nw_ref[:bus])
            for (i, obj) in nw_ref[Symbol(type)]
                push!(bus_objs[obj["$(type)_bus"]], i)
            end
            nw_ref[Symbol("bus_$(type)s")] = bus_objs
        end

        for type in _math_edge_elements
            bus_arcs = Dict((i, Tuple{Int,Int,Int}[]) for (i,bus) in nw_ref[:bus])
            for (l,i,j) in nw_ref[Symbol("arcs_$type")]
                push!(bus_arcs[i], (l,i,j))
            end
            nw_ref[Symbol("bus_arcs_$(type)")] = bus_arcs
        end

        ### connections
        for (type, status) in pmd_math_component_status
            if type in _math_node_elements
                conns = Dict{Int,Vector{Tuple{Int,Vector{Int}}}}([(i, []) for (i, bus) in nw_ref[:bus]])
                for (i, obj) in nw_ref[Symbol(type)]
                    if obj[status] != pmd_math_component_status_inactive[type]
                        push!(conns[obj["$(type)_bus"]], (i, obj["connections"]))
                    end
                end
                nw_ref[Symbol("bus_conns_$(type)")] = conns
            elseif type in _math_edge_elements
                conns = Dict{Int,Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}}([(i, []) for (i, bus) in nw_ref[:bus]])
                for (i, obj) in nw_ref[Symbol(type)]
                    if obj[status] != pmd_math_component_status_inactive[type]
                        push!(conns[obj["f_bus"]], ((obj["index"], obj["f_bus"], obj["t_bus"]), obj["f_connections"]))
                        if obj["f_bus"] != obj["t_bus"]
                            push!(conns[obj["t_bus"]], ((obj["index"], obj["t_bus"], obj["f_bus"]), obj["t_connections"]))
                        end
                    end
                end
                nw_ref[Symbol("bus_arcs_conns_$(type)")] = conns
            end
        end

        ### reference bus lookup (a set to support multiple connected components) ###
        ref_buses = Dict{Int,Any}()
        for (k,v) in nw_ref[:bus]
            if v["bus_type"] == 3
                ref_buses[k] = v
            end
        end

        nw_ref[:ref_buses] = ref_buses

        if length(ref_buses) > 1
            @debug "multiple reference buses found, $(keys(ref_buses)), this can cause infeasibility if they are in the same connected component"
        end

        ### aggregate info for pairs of connected buses ###
        if !haskey(nw_ref, :buspairs)
            nw_ref[:buspairs] = calc_buspair_parameters(nw_ref[:bus], nw_ref[:branch])
        end
    end
end
