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

"Adds connectivity for storage subsystems to converters"
function ref_add_converter_storage!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if _IM.ismultinetwork(data)
        nws_data = data["nw"]
    else
        nws_data = Dict("0" => data)
    end

    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        nw_ref = ref[:nw][nw_id]

        if !haskey(nw_ref, :converter)
            # this might happen when parsing data from matlab format
            # the OpenDSS parser always inserts a trans dict
            nw_ref[:converter] = Dict{Int, Any}()
        end

        converter_storage = Dict((i, Int[]) for (i,converter) in nw_ref[:converter])
        for (i,storage) in nw_ref[:storage]
            push!(converter_storage[storage["converter"]], i)
        end
        nw_ref[:converter_storage] = converter_storage
    end
end

"Adds connectivity for converters to buses"
function ref_add_bus_converter!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if _IM.ismultinetwork(data)
        nws_data = data["nw"]
    else
        nws_data = Dict("0" => data)
    end

    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        nw_ref = ref[:nw][nw_id]

        if !haskey(nw_ref, :converter)
            # this might happen when parsing data from matlab format
            # the OpenDSS parser always inserts a trans dict
            nw_ref[:converter] = Dict{Int, Any}()
        end

        bus_converter = Dict((i, Int[]) for (i,bus) in nw_ref[:bus])
        for (i,converter) in nw_ref[:converter]
            push!(bus_converter[converter["converter_bus"]], i)
        end
        nw_ref[:bus_converter] = bus_converter
    end
end
