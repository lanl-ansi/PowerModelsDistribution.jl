const _dimensionalize_math = Dict(
    "bus"  => Dict("rad2deg"=>["va"], "vbase"=>["vm", "vr", "vi"]),
    "gen"  => Dict("sbase"=>["pg", "qg", "pg_bus", "qg_bus"]),
    "load" => Dict("sbase"=>["pd", "qd", "pd_bus", "qd_bus"]),
    "line" => Dict("sbase"=>["pf", "qf", "pt", "qt"]),
)


"converts data model between per-unit and SI units"
function make_per_unit!(data::Dict{String,<:Any})
    data_model_type = get(data, "data_model", "mathematical")

    if  data_model_type == "mathematical"
        if !get(data, "per_unit", false)
            bus_indexed_id = string(data["bus_lookup"][data["settings"]["base_bus"]])
            vbases = Dict(bus_indexed_id=>data["settings"]["vbase"])
            sbase = data["settings"]["sbase"]

            _make_math_per_unit!(data, vbases=vbases, sbase=sbase, v_var_scalar=data["settings"]["v_var_scalar"])
        else
            # make math model si units
        end
    elseif data_model_type == "engineering"
        if !get(data, "per_unit", false)
            # make eng model per unit
        else
            # make eng model si units
        end
    else
        Memento.warn(_LOGGER, "Data model '$data_model_type' is not recognized, no per-unit transformation performed")
    end
end


"finds voltage zones"
function _find_zones(data_model::Dict{String,<:Any})
    unused_line_ids = Set(keys(data_model["branch"]))
    bus_lines = Dict([(id,Set()) for id in keys(data_model["bus"])])
    for (line_id,line) in data_model["branch"]
        f_bus = string(line["f_bus"])
        t_bus = string(line["t_bus"])
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


"calculates voltage bases for each voltage zone"
function _calc_vbase(data_model::Dict{String,<:Any}, vbase_sources::Dict{String,<:Real})
    # find zones of buses connected by lines
    zones = _find_zones(data_model)
    bus_to_zone = Dict([(bus,zone) for (zone, buses) in zones for bus in buses])

    # assign specified vbase to corresponding zones
    zone_vbase = Dict{Int, Union{Missing,Real}}([(zone,missing) for zone in keys(zones)])
    for (bus,vbase) in vbase_sources
        if !ismissing(zone_vbase[bus_to_zone[bus]])
            Memento.warn(_LOGGER, "You supplied multiple voltage bases for the same zone; ignoring all but the last one.")
        end
        zone_vbase[bus_to_zone[bus]] = vbase
    end

    # transformers form the edges between these zones
    zone_edges = Dict([(zone,[]) for zone in keys(zones)])
    edges = Set()
    for (i,(_,transformer)) in enumerate(get(data_model, "transformer", Dict{Any,Dict{String,Any}}()))
        push!(edges,i)
        f_zone = bus_to_zone[string(transformer["f_bus"])]
        t_zone = bus_to_zone[string(transformer["t_bus"])]
        tm_nom = transformer["configuration"]=="delta" ? transformer["tm_nom"]/sqrt(3) : transformer["tm_nom"]
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
    line_vbase = Dict([(id, bus_vbase[string(line["f_bus"])]) for (id,line) in data_model["branch"]])
    return (bus_vbase, line_vbase)
end


"converts to per unit from SI"
function _make_math_per_unit!(data_model::Dict{String,<:Any}; settings::Union{Missing,Dict{String,<:Any}}=missing, sbase::Union{Real,Missing}=1.0, vbases::Union{Dict{String,<:Real},Missing}=missing, v_var_scalar::Union{Missing,Real}=missing)
    if ismissing(sbase)
        if !ismissing(settings) && haskey(settings, "set_sbase")
            sbase = settings["set_sbase"]
        else
            sbase = 1.0
        end
    end

    if haskey(data_model, "sbase")
        sbase_old = data_model["sbase"]
    else
        sbase_old = 1.0
    end

    # automatically find a good vbase
    if ismissing(vbases)
        if !ismissing(settings) && haskey(settings, "set_vbases")
            vbases = settings["vbases"]
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

    for (id, line) in data_model["branch"]
        vbase = line_vbase[id]
        _rebase_pu_branch!(line, line_vbase[id], sbase, sbase_old, v_var_scalar)
    end

    for (id, shunt) in data_model["shunt"]
        _rebase_pu_shunt!(shunt, bus_vbase[string(shunt["shunt_bus"])], sbase, sbase_old, v_var_scalar)
    end

    for (id, load) in data_model["load"]
        _rebase_pu_load!(load, bus_vbase[string(load["load_bus"])], sbase, sbase_old, v_var_scalar)
    end

    for (id, gen) in data_model["gen"]
        _rebase_pu_generator!(gen, bus_vbase[string(gen["gen_bus"])], sbase, sbase_old, v_var_scalar, data_model)
    end

    for (id, storage) in data_model["storage"]
        # TODO
    end

    for (id, switch) in data_model["switch"]
        # TODO are there any properties that need to be converted to pu?
    end

    if haskey(data_model, "transformer") # transformers are not required by PowerModels
        for (id, trans) in data_model["transformer"]
            # voltage base across transformer does not have to be consistent with the ratio!
            f_vbase = bus_vbase[string(trans["f_bus"])]
            t_vbase = bus_vbase[string(trans["t_bus"])]
            _rebase_pu_transformer_2w_ideal!(trans, f_vbase, t_vbase, sbase_old, sbase, v_var_scalar)
        end
    end

    data_model["sbase"] = sbase
    data_model["per_unit"] = true

    return data_model
end


"per-unit conversion for buses"
function _rebase_pu_bus!(bus::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real, v_var_scalar::Real)
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

    # TODO fix
    if haskey(bus ,"va")
        bus["va"] = deg2rad.(bus["va"])
    end

    # save new vbase
    bus["vbase"] = vbase
end


"per-unit conversion for branches"
function _rebase_pu_branch!(branch::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real, v_var_scalar::Real)
    if !haskey(branch, "vbase")
        z_old = 1
    else
        z_old = vbase_old^2/sbase_old*v_var_scalar
    end

    z_new = vbase^2/sbase*v_var_scalar
    z_scale = z_old/z_new
    y_scale = 1/z_scale

    _scale_props!(branch, ["br_r", "br_x"], z_scale)
    _scale_props!(branch, ["b_fr", "g_fr", "b_to", "g_to"], y_scale)

    branch["angmin"] = deg2rad.(branch["angmin"])
    branch["angmax"] = deg2rad.(branch["angmax"])

    # save new vbase
    branch["vbase"] = vbase
end


"per-unit conversion for shunts"
function _rebase_pu_shunt!(shunt::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real, v_var_scalar::Real)
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


"per-unit conversion for loads"
function _rebase_pu_load!(load::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real, v_var_scalar::Real)
    if !haskey(load, "vbase")
        vbase_old = 1
        sbase_old = 1
    else
        vbase_old = load["vbase"]
    end

    vbase_old = get(load, "vbase", 1.0)
    vbase_scale = vbase_old/vbase
    scale(load, "vnom_kv", vbase_scale)

    sbase_scale = sbase_old/sbase
    scale(load, "pd", sbase_scale)
    scale(load, "qd", sbase_scale)

    # save new vbase
    load["vbase"] = vbase
end


"per-unit conversion for generators"
function _rebase_pu_generator!(gen::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real, v_var_scalar::Real, data_model::Dict{String,<:Any})
    vbase_old = get(gen, "vbase", 1.0/v_var_scalar)
    vbase_scale = vbase_old/vbase
    sbase_scale = sbase_old/sbase

    for key in ["pg", "qg", "pmin", "qmin", "pmax", "qmax"]
        scale(gen, key, sbase_scale)
    end

    # if not in per unit yet, the cost has is in $/MWh
    if !haskey(data_model, "sbase")
        sbase_old_cost = 1E6/v_var_scalar
        sbase_scale_cost = sbase_old_cost/sbase
    else
        sbase_scale_cost = sbase_scale
    end

    _PM._rescale_cost_model!(gen, 1/sbase_scale_cost)

    # save new vbase
    gen["vbase"] = vbase
end


"per-unit conversion for ideal 2-winding transformers"
function _rebase_pu_transformer_2w_ideal!(transformer::Dict{String,<:Any}, f_vbase_new::Real, t_vbase_new::Real, sbase_old::Real, sbase_new::Real, v_var_scalar::Real)
    f_vbase_old = get(transformer, "f_vbase", 1.0)
    t_vbase_old = get(transformer, "t_vbase", 1.0)
    f_vbase_scale = f_vbase_old/f_vbase_new
    t_vbase_scale = t_vbase_old/t_vbase_new

    scale(transformer, "tm_nom", f_vbase_scale/t_vbase_scale)

    # save new vbase
    transformer["f_vbase"] = f_vbase_new
    transformer["t_vbase"] = t_vbase_new
end


"helper function to apply a scale factor to given properties"
function _scale_props!(comp::Dict{String,<:Any}, prop_names::Vector{String}, scale::Real)
    for name in prop_names
        if haskey(comp, name)
            comp[name] *= scale
        end
    end
end


""
function add_big_M!(data_model::Dict{String,<:Any}; kwargs...)
    big_M = Dict{String, Any}()

    big_M["v_phase_pu_min"] = add_kwarg!(big_M, kwargs, :v_phase_pu_min, 0.5)
    big_M["v_phase_pu_max"] = add_kwarg!(big_M, kwargs, :v_phase_pu_max, 2.0)
    big_M["v_neutral_pu_min"] = add_kwarg!(big_M, kwargs, :v_neutral_pu_min, 0.0)
    big_M["v_neutral_pu_max"] = add_kwarg!(big_M, kwargs, :v_neutral_pu_max, 0.5)

    data_model["big_M"] = big_M
end


""
function solution_make_si(solution::Dict{String,<:Any}, math_model::Dict{String,<:Any}; mult_sbase::Bool=true, mult_vbase::Bool=true, convert_rad2deg::Bool=true)
    solution_si = deepcopy(solution)

    sbase = math_model["sbase"]

    for (comp_type, comp_dict) in [(x,y) for (x,y) in solution_si if isa(y, Dict)]
        dimensionalize_math_comp = get(_dimensionalize_math, comp_type, Dict())
        vbase_props   = mult_vbase      ? get(dimensionalize_math_comp, "vbase", [])   : []
        sbase_props   = mult_sbase      ? get(dimensionalize_math_comp, "sbase", [])   : []
        rad2deg_props = convert_rad2deg ? get(dimensionalize_math_comp, "rad2deg", []) : []

        for (id, comp) in comp_dict
            if !isempty(vbase_props)
                vbase = math_model[comp_type][id]["vbase"]
            end

            for (prop, val) in comp
                if prop in vbase_props
                    comp[prop] = val*vbase
                elseif prop in sbase_props
                    comp[prop] = val*sbase
                elseif prop in rad2deg_props
                    comp[prop] = _wrap_to_180(rad2deg.(val))
                end
            end
        end
    end

    solution_si["per_unit"] = false

    return solution_si
end
