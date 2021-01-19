import LinearAlgebra: diagm


""
function _calc_mc_voltage_product_bounds(pm::_PM.AbstractPowerModel, buspairs; nw::Int=pm.cnw)
    wr_min = Dict([(bp, -Inf) for bp in buspairs])
    wr_max = Dict([(bp,  Inf) for bp in buspairs])
    wi_min = Dict([(bp, -Inf) for bp in buspairs])
    wi_max = Dict([(bp,  Inf) for bp in buspairs])

    for (i, j, c, d) in buspairs
        if i == j
            bus = ref(pm, nw, :bus)[i]
            vm_fr_max = bus["vmax"][c]
            vm_to_max = bus["vmax"][d]
            vm_fr_min = bus["vmin"][c]
            vm_to_min = bus["vmin"][d]
        else
            buspair = ref(pm, nw, :buspairs)[(i, j)]
            vm_fr_max = buspair["vm_fr_max"][c]
            vm_to_max = buspair["vm_to_max"][d]
            vm_fr_min = buspair["vm_fr_min"][c]
            vm_to_min = buspair["vm_to_min"][d]
        end

        wr_max[(i, j, c, d)] =  vm_fr_max * vm_to_max
        wr_min[(i, j, c, d)] = -vm_fr_max * vm_to_max
        wi_max[(i, j, c, d)] =  vm_fr_max * vm_to_max
        wi_min[(i, j, c, d)] = -vm_fr_max * vm_to_max
    end

    return wr_min, wr_max, wi_min, wi_max
end


""
function _find_ref_buses(pm::_PM.AbstractPowerModel, nw)
    buses = ref(pm, nw, :bus)
    return [b for (b,bus) in buses if bus["bus_type"]==3]
    # return [bus for (b,bus) in buses ]
end


"Adds arcs for PowerModelsDistribution transformers; for dclines and branches this is done in PowerModels"
function ref_add_arcs_transformer!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if _IM.ismultinetwork(data)
        nws_data = data["nw"]
    else
        nws_data = Dict("0" => data)
    end

    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        nw_ref = ref[:nw][nw_id]

        if !haskey(nw_ref, :transformer)
            # this might happen when parsing data from matlab format
            # the OpenDSS parser always inserts a trans dict
            nw_ref[:transformer] = Dict{Int, Any}()
        end

        nw_ref[:arcs_from_trans] = [(i, trans["f_bus"], trans["t_bus"]) for (i,trans) in nw_ref[:transformer]]
        nw_ref[:arcs_to_trans] = [(i, trans["t_bus"], trans["f_bus"]) for (i,trans) in nw_ref[:transformer]]
        nw_ref[:arcs_trans] = [nw_ref[:arcs_from_trans]..., nw_ref[:arcs_to_trans]...]
        nw_ref[:bus_arcs_trans] = Dict{Int64, Array{Any, 1}}()

        for (i,bus) in nw_ref[:bus]
            nw_ref[:bus_arcs_trans][i] = [e for e in nw_ref[:arcs_trans] if e[2]==i]
        end
    end
end


"Adds arcs for PowerModelsDistribution transformers; for dclines and branches this is done in PowerModels"
function ref_add_arcs_switch!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if _IM.ismultinetwork(data)
        nws_data = data["nw"]
    else
        nws_data = Dict("0" => data)
    end

    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        nw_ref = ref[:nw][nw_id]

        if !haskey(nw_ref, :switch)
            # this might happen when parsing data from matlab format
            # the OpenDSS parser always inserts a trans dict
            nw_ref[:switch] = Dict{Int, Any}()
        end

        nw_ref[:arcs_from_sw] = [(i, switch["f_bus"], switch["t_bus"]) for (i,switch) in nw_ref[:switch]]
        nw_ref[:arcs_to_sw] = [(i, switch["t_bus"], switch["f_bus"]) for (i,switch) in nw_ref[:switch]]
        nw_ref[:arcs_sw] = [nw_ref[:arcs_from_sw]..., nw_ref[:arcs_to_sw]...]
        nw_ref[:bus_arcs_sw] = Dict{Int64, Array{Any, 1}}()

        for (i,bus) in nw_ref[:bus]
            nw_ref[:bus_arcs_sw][i] = [e for e in nw_ref[:arcs_sw] if e[2]==i]
        end
    end
end



""
function _calc_mc_transformer_Tvi(pm::_PM.AbstractPowerModel, i::Int; nw=pm.cnw)
    trans = ref(pm, nw, :transformer,  i)
    # transformation matrices
    # Tv and Ti will be compositions of these
    Tbr = [0 0 1; 1 0 0; 0 1 0]                             # barrel roll
    Tdelt  = [1 -1 0; 0 1 -1; -1 0 1]                       # delta transform
    # grounding disregarded
    for config in [trans["config_to"], trans["config_fr"]]
        if haskey(config, "grounded") && config["grounded"]==false
            Memento.warning(_LOGGER, "The wye winding is considered to be grounded instead of ungrounded.")
        end
    end
    # make sure the secondary is y+123
    if trans["config_to"]["type"]!=WYE
        Memento.error(_LOGGER, "Secondary should always be of wye type.")
    end
    if trans["config_to"]["cnd"]!=[1,2,3]
        Memento.error(_LOGGER, "Secondary should always be connected in 123.")
    end
    if trans["config_to"]["polarity"]!='+'
        Memento.error(_LOGGER, "Secondary should always be of positive polarity.")
    end
    # connection transformers
    perm = trans["config_fr"]["cnd"]
    if !(perm in [[1,2,3], [3,1,2], [2,3,1]])
        Memento.error(_LOGGER, "Only the permutations \"123\", \"312\" and \"231\" are supported, but got \"$perm\".")
    end
    polarity = trans["config_fr"]["polarity"]
    if !(polarity in ['+', '-'])
        Memento.error(_LOGGER, "The polarity should be either \'+\' or \'-\', but got \'$polarity\'.")
    end
    dyz = trans["config_fr"]["type"]
    if !(dyz in [DELTA, WYE])
        Memento.error(_LOGGER, "The winding type should be either delta or wye, but got \'$dyz\'.")
    end
    # for now, grounded by default
    #grounded = length(trans["configuration"])>5 && trans["configuration"][6]=='n'
    # Tw will contain transformations related to permutation and polarity
    perm_to_trans = Dict(
        [1,2,3]=>diagm(0=>ones(Float64, 3)),
        [3,1,2]=>Tbr,
        [2,3,1]=>Tbr*Tbr
    )
    Tw = perm_to_trans[perm]
    Tw = (polarity=='+') ? Tw : -Tw
    #Tw = diagm(0=>ones(Float64, 3))
    vmult = 1.0 # compensate for change in LN
    if dyz==WYE
        Tv_fr = Tw
        Tv_im = diagm(0=>ones(Float64, 3))
        Ti_fr = Tw
        Ti_im = diagm(0=>ones(Float64, 3))
        # if !grounded
        #     # if not grounded, phase currents should sum to zero
        #     Ti_fr = [Ti_fr; ones(1,3)]
        #     Ti_im = [Ti_im; zeros(1,3)]
        # end
    elseif dyz==DELTA
        Tv_fr = Tdelt*Tw
        Tv_im = diagm(0=>ones(Float64, 3))
        Ti_fr = Tw
        Ti_im = Tdelt'
        vmult = sqrt(3)
    elseif dyz=="zig-zag"
        #TODO zig-zag here
    end
    # make equations dimensionless
    # if vbase across a transformer scales according to the ratio of vnom_kv,
    # this will simplify to 1.0
    bkv_fr = ref(pm, nw, :bus, trans["f_bus"], "base_kv")
    bkv_to = ref(pm, nw, :bus, trans["t_bus"], "base_kv")
    Cv_to = trans["config_fr"]["vm_nom"]/trans["config_to"]["vm_nom"]*bkv_to/bkv_fr
    # compensate for change of LN voltage of a delta winding
    Cv_to *= vmult
    return (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to)
end


""
function ref_add_connections!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if ismultinetwork(data)
        nws_data = data["nw"]
    else
        nws_data = Dict("0" => data)
    end

    for (n, nw) in nws_data
        n = parse(Int, n)
        nw_ref = ref[:nw][n]

        for (type, status) in [("gen", "gen_status"), ("load", "status"), ("shunt", "status"), ("storage", "status")]
            nw_ref[Symbol("bus_conns_$(type)")] = Dict{Int,Vector{Tuple{Int,Vector{Int}}}}([(bus["index"], []) for (_,bus) in nw["bus"]])
            for (_,obj) in get(nw, type, Dict())
                if obj[status] != 0
                    push!(nw_ref[Symbol("bus_conns_$(type)")][obj["$(type)_bus"]], (obj["index"], obj["connections"]))
                end
            end
        end

        for (type, status) in [("transformer", "status"), ("branch", "br_status"), ("switch", "status")]
            nw_ref[Symbol("bus_arcs_conns_$(type)")] = Dict{Int,Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}}([(bus["index"], []) for (_,bus) in nw["bus"]])
            for (_,obj) in get(nw, type, Dict())
                if obj[status] != 0
                    push!(nw_ref[Symbol("bus_arcs_conns_$(type)")][obj["f_bus"]], ((obj["index"], obj["f_bus"], obj["t_bus"]), obj["f_connections"]))
                    push!(nw_ref[Symbol("bus_arcs_conns_$(type)")][obj["t_bus"]], ((obj["index"], obj["t_bus"], obj["f_bus"]), obj["t_connections"]))
                end
            end
        end
    end
end


"computes storage bounds"
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


"overwrites PowerModels buspairs ref by conductor"
function _PM.calc_buspair_parameters(buses, branches, conductor_ids, is_multiconductor::Bool)
    bus_lookup = Dict(bus["index"] => bus for (i,bus) in buses if bus["bus_type"] != 4)
    branch_lookup = Dict(branch["index"] => branch for (i,branch) in branches if branch["br_status"] == 1 && haskey(bus_lookup, branch["f_bus"]) && haskey(bus_lookup, branch["t_bus"]))

    if is_multiconductor
        buspair_indexes = Set((branch["f_bus"], branch["t_bus"], fc, tc) for (i,branch) in branch_lookup for (fc, tc) in zip(branch["f_connections"], branch["t_connections"]))
    else
        buspair_indexes = Set((branch["f_bus"], branch["t_bus"]) for (i,branch) in branch_lookup)
    end

    bp_branch = Dict((bp, typemax(Int64)) for bp in buspair_indexes)

    bp_angmin = Dict((bp, -Inf) for bp in buspair_indexes)
    bp_angmax = Dict((bp,  Inf) for bp in buspair_indexes)

    for (l,branch) in branch_lookup
        i = branch["f_bus"]
        j = branch["t_bus"]

        if is_multiconductor
            f_connections = branch["f_connections"]
            t_connections = branch["t_connections"]

            for (idx, (fc, tc)) in enumerate(zip(f_connections, t_connections))
                bp_angmin[(i,j,fc,tc)] = max(bp_angmin[(i,j,fc,tc)], branch["angmin"][idx])
                bp_angmax[(i,j,fc,tc)] = min(bp_angmax[(i,j,fc,tc)], branch["angmax"][idx])

                bp_branch[(i,j,fc,tc)] = min(bp_branch[i,j,fc,tc], l)
            end
        else
            bp_angmin[(i,j)] = max(bp_angmin[(i,j)], branch["angmin"])
            bp_angmax[(i,j)] = min(bp_angmax[(i,j)], branch["angmax"])

            bp_branch[(i,j)] = min(bp_branch[(i,j)], l)
        end
    end

    if is_multiconductor
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
    else
        buspairs = Dict((i,j) => Dict(
            "branch"=>bp_branch[(i,j)],
            "angmin"=>bp_angmin[(i,j)],
            "angmax"=>bp_angmax[(i,j)],
            "tap"=>get(branch_lookup[bp_branch[(i,j)]], "tap", 1.0),
            "vm_fr_min"=>bus_lookup[i]["vmin"],
            "vm_fr_max"=>bus_lookup[i]["vmax"],
            "vm_to_min"=>bus_lookup[j]["vmin"],
            "vm_to_max"=>bus_lookup[j]["vmax"]
            ) for (i,j) in buspair_indexes
        )
    end

    # add optional parameters
    for bp in buspair_indexes
        branch = branch_lookup[bp_branch[bp]]
        if haskey(branch, "rate_a")
            if is_multiconductor
                buspairs[bp]["rate_a"] = branch["rate_a"][findfirst(isequal(bp[end]), branch["t_connections"])]
            else
                buspairs[bp]["rate_a"] = branch["rate_a"]
            end
        end
        if haskey(branch, "c_rating_a")
            if is_multiconductor
                buspairs[bp]["c_rating_a"] = branch["c_rating_a"][findfirst(isequal(bp[end]), branch["t_connections"])]
            else
                buspairs[bp]["c_rating_a"] = branch["c_rating_a"]
            end
        end
    end

end


""
function find_conductor_ids!(data::Dict{String,<:Any})
    conductor_ids = []

    for (_,bus) in get(data, "bus", Dict())
        for t in get(bus, "terminals", [])
            if !(t in conductor_ids)
                push!(conductor_ids, t)
            end
        end
    end

    data["conductors"] = length(conductor_ids)
    data["conductor_ids"] = [c for c in conductor_ids]
end
