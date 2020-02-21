
function _find_zones(data_model)
    unused_line_ids = Set(keys(data_model["line"]))
    bus_lines = Dict([(id,Set()) for id in keys(data_model["bus"])])
    for (line_id,line) in data_model["line"]
        f_bus = line["f_bus"]
        t_bus = line["t_bus"]
        push!(bus_lines[f_bus], (line_id,t_bus))
        push!(bus_lines[t_bus], (line_id,f_bus))
    end
    zones = []
    buses = Set(keys(data_model["bus"]))
    while !isempty(buses)
        stack = [pop!(buses)]
        zone = Set()
        while !isempty(stack)
            bus = pop!(stack)
            delete!(buses, bus)
            push!(zone, bus)
            for (line_id,bus_to) in bus_lines[bus]
                if line_id in unused_line_ids && bus_to in buses
                    delete!(unused_line_ids, line_id)
                    push!(stack, bus_to)
                end
            end
        end
        append!(zones, [zone])
    end
    zones = Dict(enumerate(zones))

    return zones
end


function _calc_vbase(data_model, vbase_sources::Dict{String,<:Real})
    # find zones of buses connected by lines
    zones = _find_zones(data_model)
    bus_to_zone = Dict([(bus,zone) for (zone, buses) in zones for bus in buses])

    # assign specified vbase to corresponding zones
    zone_vbase = Dict{Int, Union{Missing,Real}}([(zone,missing) for zone in keys(zones)])
    for (bus,vbase) in vbase_sources
        if !ismissing(zone_vbase[bus_to_zone[bus]])
            Memento.warn("You supplied multiple voltage bases for the same zone; ignoring all but the last one.")
        end
        zone_vbase[bus_to_zone[bus]] = vbase
    end

    # transformers form the edges between these zones
    zone_edges = Dict([(zone,[]) for zone in keys(zones)])
    edges = Set()
    for (i,(_,trans)) in enumerate(data_model["transformer_2wa"])
        push!(edges,i)
        f_zone = bus_to_zone[trans["f_bus"]]
        t_zone = bus_to_zone[trans["t_bus"]]
        tm_nom = trans["configuration"]=="delta" ? trans["tm_nom"]/sqrt(3) : trans["tm_nom"]
        push!(zone_edges[f_zone], (i, t_zone, 1/tm_nom))
        push!(zone_edges[t_zone], (i, f_zone, tm_nom))
    end

    # initialize the stack with all specified zones
    stack = [zone for (zone,vbase) in zone_vbase if !ismissing(vbase)]

    while !isempty(stack)

        zone = pop!(stack)

        for (edge_id, zone_to, scale) in zone_edges[zone]
            delete!(edges, edge_id)

            if ismissing(zone_vbase[zone_to])
                zone_vbase[zone_to] = zone_vbase[zone]*scale
                push!(stack, zone_to)
            end
        end
    end

    bus_vbase = Dict([(bus,zone_vbase[zone]) for (bus,zone) in bus_to_zone])
    line_vbase = Dict([(id, bus_vbase[line["f_bus"]]) for (id,line) in data_model["line"]])
    return (bus_vbase, line_vbase)
end


function data_model_make_pu!(data_model; sbase=missing, vbases=missing)
    v_var_scalar = data_model["settings"]["v_var_scalar"]

    if haskey(data_model["settings"], "sbase")
        sbase_old = data_model["settings"]["sbase"]
    else
        sbase_old = 1
    end

    if ismissing(sbase)
        if haskey(data_model["settings"], "set_sbase_val")
            sbase = data_model["settings"]["set_sbase_val"]
        else
            sbase = 1
        end
    end

    # automatically find a good vbase
    if ismissing(vbases)
        if haskey(data_model["settings"], "set_vbase_val") && haskey(data_model["settings"], "set_vbase_bus")
            vbases = Dict(data_model["settings"]["set_vbase_bus"]=>data_model["settings"]["set_vbase_val"])
        else
            buses_type_3 = [(id, sum(bus["vm"])/length(bus["vm"])) for (id,bus) in data_model["bus"] if haskey(bus, "bus_type") && bus["bus_type"]==3]
            if !isempty(buses_type_3)
                vbases = Dict([buses_type_3[1]])
            else
                Memento.error("Please specify vbases manually; cannot make an educated guess for this data model.")
            end
        end
    end

    bus_vbase, line_vbase = _calc_vbase(data_model, vbases)

    for (id, bus) in data_model["bus"]
        _rebase_pu_bus!(bus, bus_vbase[id], sbase, sbase_old, v_var_scalar)
    end

    for (id, line) in data_model["line"]
        vbase = line_vbase[id]
        _rebase_pu_line!(line, line_vbase[id], sbase, sbase_old, v_var_scalar)
    end

    for (id, shunt) in data_model["shunt"]
        _rebase_pu_shunt!(shunt, bus_vbase[shunt["bus"]], sbase, sbase_old, v_var_scalar)
    end

    for (id, load) in data_model["load"]
        _rebase_pu_load!(load, bus_vbase[load["bus"]], sbase, sbase_old, v_var_scalar)
    end

    for (id, gen) in data_model["generator"]
        _rebase_pu_generator!(gen, bus_vbase[gen["bus"]], sbase, sbase_old, v_var_scalar)
    end

    for (id, trans) in data_model["transformer_2wa"]
        # voltage base across transformer does not have to be consistent with the ratio!
        f_vbase = bus_vbase[trans["f_bus"]]
        t_vbase = bus_vbase[trans["t_bus"]]
        _rebase_pu_transformer_2w_ideal!(trans, f_vbase, t_vbase, sbase_old, sbase, v_var_scalar)
    end

    data_model["settings"]["sbase"] = sbase

    return data_model
end


function _rebase_pu_bus!(bus, vbase, sbase, sbase_old, v_var_scalar)

    # if not in p.u., these are normalized with respect to vnom
    prop_vnom = ["vm", "vmax", "vmin", "vm_set", "vm_ln_min", "vm_ln_max", "vm_lg_min", "vm_lg_max", "vm_ng_min", "vm_ng_max", "vm_ll_min", "vm_ll_max"]

    if !haskey(bus, "vbase")

        # if haskey(bus, "vnom")
        #     vnom = bus["vnom"]
        #     _scale_props!(bus, ["vnom"], 1/vbase)
        # end
        _scale_props!(bus, prop_vnom, 1/vbase)

        z_old = 1.0
    else
        vbase_old = bus["vbase"]
        _scale_props!(bus, [prop_vnom..., "vnom"], vbase_old/vbase)

        z_old = vbase_old^2*sbase_old*v_var_scalar
    end

    # rebase grounding resistance
    z_new = vbase^2/sbase*v_var_scalar
    z_scale = z_old/z_new
    _scale_props!(bus, ["rg", "xg"], z_scale)

    # save new vbase
    bus["vbase"] = vbase
end


function _rebase_pu_line!(line, vbase, sbase, sbase_old, v_var_scalar)

    if !haskey(line, "vbase")
        z_old = 1
    else
        z_old = vbase_old^2/sbase_old*v_var_scalar
    end

    z_new = vbase^2/sbase*v_var_scalar
    z_scale = z_old/z_new
    y_scale = 1/z_scale

    _scale_props!(line, ["rs", "xs"], z_scale)
    _scale_props!(line, ["b_fr", "g_fr", "b_to", "g_to"], y_scale)

    # save new vbase
    line["vbase"] = vbase
end


function _rebase_pu_shunt!(shunt, vbase, sbase, sbase_old, v_var_scalar)

    if !haskey(shunt, "vbase")
        z_old = 1
    else
        z_old = vbase_old^2/sbase_old*v_var_scalar
    end

    # rebase grounding resistance
    z_new = vbase^2/sbase*v_var_scalar
    z_scale = z_old/z_new
    y_scale = 1/z_scale
    scale(shunt, "gs", y_scale)
    scale(shunt, "bs", y_scale)

    # save new vbase
    shunt["vbase"] = vbase
end


function _rebase_pu_load!(load, vbase, sbase, sbase_old, v_var_scalar)

    if !haskey(load, "vbase")
        vbase_old = 1
        sbase_old = 1
    else
        vbase_old = load["vbase"]
    end

    vbase_old = get(load, "vbase", 1.0)
    vbase_scale = vbase_old/vbase
    scale(load, "vnom", vbase_scale)

    sbase_scale = sbase_old/sbase
    scale(load, "pd", sbase_scale)
    scale(load, "qd", sbase_scale)

    # save new vbase
    load["vbase"] = vbase
end


function _rebase_pu_generator!(gen, vbase, sbase, sbase_old, v_var_scalar)
    vbase_old = get(gen, "vbase", 1.0/v_var_scalar)
    vbase_scale = vbase_old/vbase
    sbase_scale = sbase_old/sbase

    for key in ["pd_set", "qd_set", "pd_min", "qd_min", "pd_max", "qd_max"]
        scale(gen, key, sbase_scale)
    end

    scale(gen, "cost", 1/sbase_scale)

    # save new vbase
    gen["vbase"] = vbase
end


function _rebase_pu_transformer_2w_ideal!(trans, f_vbase_new, t_vbase_new, sbase_old, sbase_new, v_var_scalar)
    f_vbase_old = get(trans, "f_vbase", 1.0)
    t_vbase_old = get(trans, "t_vbase", 1.0)
    f_vbase_scale = f_vbase_old/f_vbase_new
    t_vbase_scale = t_vbase_old/t_vbase_new

    scale(trans, "tm_nom", f_vbase_scale/t_vbase_scale)

    # save new vbase
    trans["f_vbase"] = f_vbase_new
    trans["t_vbase"] = t_vbase_new
end


function _scale_props!(comp::Dict{String, Any}, prop_names::Array{String, 1}, scale::Real)
    for name in prop_names
        if haskey(comp, name)
            comp[name] *= scale
        end
    end
end

#data_model_user = make_test_data_model()
#data_model_base = map_down_data_model(data_model_user)
#bus_vbase, line_vbase = get_vbase(data_model_base, Dict("1"=>230.0))

#make_pu!(data_model_base)

function add_big_M!(data_model; kwargs...)
    big_M = Dict{String, Any}()

    big_M["v_phase_pu_min"] = add_kwarg!(big_M, kwargs, :v_phase_pu_min, 0.5)
    big_M["v_phase_pu_max"] = add_kwarg!(big_M, kwargs, :v_phase_pu_max, 2.0)
    big_M["v_neutral_pu_min"] = add_kwarg!(big_M, kwargs, :v_neutral_pu_min, 0.0)
    big_M["v_neutral_pu_max"] = add_kwarg!(big_M, kwargs, :v_neutral_pu_max, 0.5)

    data_model["big_M"] = big_M
end

function solution_unmake_pu!(solution, data_model)
    sbase = data_model["sbase"]
    for (comp_type, comp_dict) in [(x,y) for (x,y) in solution if isa(y, Dict)]
        for (id, comp) in comp_dict
            for (prop, val) in comp
                if any([occursin(x, prop) for x in ["p", "q"]])
                    comp[prop] = val*sbase
                elseif any([occursin(x, prop) for x in ["vm", "vr", "vi"]])
                    comp[prop] = val*data_model[comp_type][id]["vbase"]
                end
            end
        end
    end

    return solution
end
