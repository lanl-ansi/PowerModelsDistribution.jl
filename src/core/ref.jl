import LinearAlgebra: diagm


""
function _calc_tp_voltage_product_bounds(pm::_PMs.GenericPowerModel, buspairs; nw::Int=pm.cnw)
    wr_min = Dict([(bp, -Inf) for bp in buspairs])
    wr_max = Dict([(bp,  Inf) for bp in buspairs])
    wi_min = Dict([(bp, -Inf) for bp in buspairs])
    wi_max = Dict([(bp,  Inf) for bp in buspairs])

    for (i, j, c, d) in buspairs
        if i == j
            bus = _PMs.ref(pm, nw, :bus)[i]
            vm_fr_max = bus["vmax"][c]
            vm_to_max = bus["vmax"][d]
            vm_fr_min = bus["vmin"][c]
            vm_to_min = bus["vmin"][d]
        else
            buspair = _PMs.ref(pm, nw, :buspairs)[(i, j)]
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
function _find_ref_buses(pm::_PMs.GenericPowerModel, nw)
    buses = _PMs.ref(pm, nw, :bus)
    return [b for (b,bus) in buses if bus["bus_type"]==3]
    # return [bus for (b,bus) in buses ]
end


"Adds arcs for PMD transformers; for dclines and branches this is done in PMs"
function ref_add_arcs_trans!(pm::_PMs.GenericPowerModel)
    if !haskey(_PMs.ref(pm, pm.cnw), :trans)
        # this might happen when parsing data from matlab format
        # the OpenDSS parser always inserts a trans dict
        _PMs.ref(pm, pm.cnw)[:trans] = Dict{Int, Any}()
    end
    # dirty fix add arcs_from/to_trans and bus_arcs_trans
    pm.ref[:nw][0][:arcs_from_trans] = [(i, trans["f_bus"], trans["t_bus"]) for (i,trans) in _PMs.ref(pm, :trans)]
    pm.ref[:nw][0][:arcs_to_trans] = [(i, trans["t_bus"], trans["f_bus"]) for (i,trans) in _PMs.ref(pm, :trans)]
    pm.ref[:nw][0][:arcs_trans] = [pm.ref[:nw][0][:arcs_from_trans]..., pm.ref[:nw][0][:arcs_to_trans]...]
    pm.ref[:nw][0][:bus_arcs_trans] = Dict{Int64, Array{Any, 1}}()
    for i in _PMs.ids(pm, :bus)
        pm.ref[:nw][0][:bus_arcs_trans][i] = [e for e in pm.ref[:nw][0][:arcs_trans] if e[2]==i]
    end
end


""
function _calc_tp_trans_Tvi(pm::_PMs.GenericPowerModel, i::Int; nw=pm.cnw)
    trans = _PMs.ref(pm, nw, :trans,  i)
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
    if trans["config_to"]["type"]!="wye"
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
    if !(dyz in ["delta", "wye"])
        Memento.error(_LOGGER, "The winding type should be either delta or wye, but got \'$dyz\'.")
    end
    # for now, grounded by default
    #grounded = length(trans["conn"])>5 && trans["conn"][6]=='n'
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
    if dyz=="wye"
        Tv_fr = Tw
        Tv_im = diagm(0=>ones(Float64, 3))
        Ti_fr = Tw
        Ti_im = diagm(0=>ones(Float64, 3))
        # if !grounded
        #     # if not grounded, phase currents should sum to zero
        #     Ti_fr = [Ti_fr; ones(1,3)]
        #     Ti_im = [Ti_im; zeros(1,3)]
        # end
    elseif dyz=="delta"
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
    bkv_fr = _PMs.ref(pm, nw, :bus, trans["f_bus"], "base_kv")
    bkv_to = _PMs.ref(pm, nw, :bus, trans["t_bus"], "base_kv")
    Cv_to = trans["config_fr"]["vm_nom"]/trans["config_to"]["vm_nom"]*bkv_to/bkv_fr
    # compensate for change of LN voltage of a delta winding
    Cv_to *= vmult
    return (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to)
end


function _bus_vm_ll_bounds(bus::Dict; eps=0.1)
    vmax = bus["vmax"].values
    vmin = bus["vmin"].values
    if haskey(bus, "vm_ll_max")
        vdmax = bus["vm_ll_max"].values
    else
        vdmax = 2*vmax
    end
    if haskey(bus, "vm_ll_min")
        vdmin = bus["vm_ll_min"].values
    else
        vdmin = ones(3)*eps
    end
    return (vdmin, vdmax)
end


function _load_pq_bounds(load::Dict, bus::Dict)
    a, α, b, β = _load_expmodel_params(load, bus)
    vmin, vmax = _load_vbounds(load, bus)
    # get bounds
    pmin = min.(a.*vmin.^α, a.*vmax.^α)
    pmax = max.(a.*vmin.^α, a.*vmax.^α)
    qmin = min.(b.*vmin.^β, b.*vmax.^β)
    qmax = max.(b.*vmin.^β, b.*vmax.^β)
    return (pmin, pmax, qmin, qmax)
end

function _load_curr_max(load::Dict, bus::Dict)
    pmin, pmax, qmin, qmax = _load_pq_bounds(load, bus)
    pabsmax = max.(abs.(pmin), abs.(pmax))
    qabsmax = max.(abs.(qmin), abs.(qmax))
    smax = sqrt.(pabsmax.^2 + qabsmax.^2)

    vmin, vmax = _load_vbounds(load, bus)

    return smax./vmin
end


function _load_expmodel_params(load::Dict, bus::Dict)
    pd = load["pd"].values
    qd = load["qd"].values
    ncnds = length(pd)
    if load["model"]=="constant_power"
        return (pd, zeros(ncnds), qd, zeros(ncnds))
    else
        vmin, vmax = _load_vbounds(load, bus)
        # get exponents
        if load["model"]=="constant_current"
            α = ones(ncnds)
            β  =ones(ncnds)
        elseif load["model"]=="constant_impedance"
            α = ones(ncnds)*2
            β  =ones(ncnds)*2
        elseif load["model"]=="exponential"
            α = load["alpha"]
            β = load["beta"]
        end
        # calculate proportionality constants
        v0 = load["vnom_kv"]/(bus["base_kv"]/sqrt(3))
        a = pd./v0.^α
        b = qd./v0.^β
        # get bounds
        return (a, α, b, β)
    end
end

function _load_vbounds(load::Dict, bus::Dict)
    if load["conn"]=="wye"
        vmin = bus["vmin"].values
        vmax = bus["vmax"].values
    elseif load["conn"]=="delta"
        vmin, vmax = _bus_vm_ll_bounds(bus)
    end
    return vmin, vmax
end

"""
Returns a Bool, indicating whether the convex hull of the voltage-dependent
relationship needs a cone inclusion constraint.
"""
function _load_needs_cone(load::Dict)
    if load["model"]=="constant_current"
        return true
    elseif load["model"]=="exponential"
        return true
    else
        return false
    end
end
