"lists of scaling factors and what they apply to"
const _dimensionalize_math = Dict{String,Dict{String,Vector{String}}}(
    "bus"  => Dict{String,Vector{String}}(
        "rad2deg"=>Vector{String}(["va", "va_pp", "va_pn"]),
        "vbase"=>Vector{String}(["vm", "vr", "vi", "vm_pp", "vm_pn"])
    ),
    "gen"  => Dict{String,Vector{String}}(
        "sbase"=>Vector{String}(["pg", "qg", "pg_bus", "qg_bus"]),
        "ibase"=>Vector{String}(["crg", "cig", "crg_bus", "cig_bus"])
    ),
    "load" => Dict{String,Vector{String}}(
        "sbase"=>Vector{String}(["pd", "qd", "pd_bus", "qd_bus"]),
        "ibase"=>Vector{String}(["crd", "cid", "crd_bus", "cid_bus"])
    ),
    "branch" => Dict{String,Vector{String}}(
        "sbase"=>Vector{String}(["pf", "qf", "pt", "qt"]),
        "ibase"=>Vector{String}(["cr_fr", "ci_fr", "cr_to", "cr_to"])
    ),
    "transformer" => Dict{String,Vector{String}}(
        "ibase_fr"=>Vector{String}(["crt_fr", "cit_fr"]),
        "ibase_to"=>Vector{String}(["crt_to", "cit_to"])
    ),
)


"converts data model between per-unit and SI units"
function make_per_unit!(data::Dict{String,<:Any}; vbases::Union{Dict{<:Any,<:Real},Missing}=missing, sbase::Union{Real,Missing}=missing, data_model_type::DataModel=get(data, "data_model", MATHEMATICAL))
    if data_model_type == MATHEMATICAL
        if !get(data, "per_unit", false)
            if ismissing(vbases)
                vbases = Dict{String,Real}("$(data["bus_lookup"][id])"=>vbase for (id, vbase) in data["settings"]["vbases_default"])
            end

            if ismissing(sbase)
                sbase = data["settings"]["sbase_default"]
            end

            if ismultinetwork(data)
                for (n, nw) in data["nw"]
                    _make_math_per_unit!(nw, data; sbase=sbase, vbases=vbases)
                end
            else
                _make_math_per_unit!(data, data; sbase=sbase, vbases=vbases)
            end
        else
            # TODO make math model si units
        end
    else
        Memento.warn(_LOGGER, "Data model '$data_model_type' is not recognized, no per-unit transformation performed")
    end
end


"finds voltage zones"
function _find_zones(data_model::Dict{String,<:Any})::Dict{Int,Set{String}}
    unused_components = Set("$comp_type.$id" for comp_type in ["branch", "switch"] for id in keys(data_model[comp_type]))
    bus_connectors = Dict([(id,Set()) for id in keys(data_model["bus"])])
    for comp_type in ["branch", "switch"]
        for (id,obj) in data_model[comp_type]
            f_bus = string(obj["f_bus"])
            t_bus = string(obj["t_bus"])
            push!(bus_connectors[f_bus], ("$comp_type.$id",t_bus))
            push!(bus_connectors[t_bus], ("$comp_type.$id",f_bus))
        end
    end

    zones = []
    buses = Set(keys(data_model["bus"]))
    while !isempty(buses)
        stack = [pop!(buses)]
        zone = Set{String}()
        while !isempty(stack)
            bus = pop!(stack)
            delete!(buses, bus)
            push!(zone, bus)
            for (id,bus_to) in bus_connectors[bus]
                if id in unused_components && bus_to in buses
                    delete!(unused_components, id)
                    push!(stack, bus_to)
                end
            end
        end
        append!(zones, [zone])
    end
    zones = Dict{Int,Set{String}}(enumerate(zones))

    return zones
end


"calculates voltage bases for each voltage zone"
function _calc_vbase(data_model::Dict{String,<:Any}, vbase_sources::Dict{String,<:Real})::Tuple{Dict,Dict}
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
    zone_edges = Dict{Int,Vector{Tuple{Int,Real}}}([(zone,[]) for zone in keys(zones)])
    for (_,transformer) in get(data_model, "transformer", Dict{Any,Dict{String,Any}}())
        f_zone = bus_to_zone["$(transformer["f_bus"])"]
        t_zone = bus_to_zone["$(transformer["t_bus"])"]
        tm_nom = transformer["configuration"]==DELTA ? transformer["tm_nom"]/sqrt(3) : transformer["tm_nom"]
        push!(zone_edges[f_zone], (t_zone, 1/tm_nom))
        push!(zone_edges[t_zone], (f_zone,   tm_nom))
    end

    # initialize the stack with all specified zones
    stack = [zone for (zone,vbase) in zone_vbase if !ismissing(vbase)]
    while !isempty(stack)
        f_zone = pop!(stack)
        for (t_zone, scale) in zone_edges[f_zone]
            if ismissing(zone_vbase[t_zone])
                zone_vbase[t_zone] = zone_vbase[f_zone]*scale
                push!(stack, t_zone)
            end
        end
    end

    bus_vbase = Dict([(bus,zone_vbase[zone]) for (bus,zone) in bus_to_zone])
    line_vbase = Dict([(id, bus_vbase[string(line["f_bus"])]) for (id,line) in data_model["branch"]])
    return (bus_vbase, line_vbase)
end


"converts to per unit from SI"
function _make_math_per_unit!(data_model::Dict{String,<:Any}, data_math::Dict{String,<:Any}; sbase::Union{Real,Missing}=missing, vbases::Union{Dict{String,<:Real},Missing}=missing)
    if ismissing(sbase)
        if haskey(data_math["settings"], "sbase_default")
            sbase = data_math["settings"]["sbase_default"]
        else
            sbase = 1.0
        end
    end

    if haskey(data_math["settings"], "sbase")
        sbase_old = data_math["settings"]["sbase"]
    else
        sbase_old = 1.0
    end

    # automatically find a good vbase if not provided
    if ismissing(vbases)
        if haskey(data_math["settings"], "vbases_default")
            vbases = Dict{String,Real}("$(data_math["bus_lookup"][id])" => vbase for (id, vbase) in data_math["settings"]["vbases_default"])
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
    voltage_scale_factor = data_math["settings"]["voltage_scale_factor"]

    for (id, bus) in data_model["bus"]
        _rebase_pu_bus!(bus, bus_vbase[id], sbase, sbase_old, voltage_scale_factor)
    end

    for (id, line) in data_model["branch"]
        vbase = line_vbase[id]
        _rebase_pu_branch!(line, line_vbase[id], sbase, sbase_old, voltage_scale_factor)
    end

    for (id, shunt) in data_model["shunt"]
        _rebase_pu_shunt!(shunt, bus_vbase[string(shunt["shunt_bus"])], sbase, sbase_old, voltage_scale_factor)
    end

    for (id, load) in data_model["load"]
        _rebase_pu_load!(load, bus_vbase[string(load["load_bus"])], sbase, sbase_old, voltage_scale_factor)
    end

    for (id, gen) in data_model["gen"]
        _rebase_pu_generator!(gen, bus_vbase[string(gen["gen_bus"])], sbase, sbase_old, data_math)
    end

    for (id, storage) in data_model["storage"]
        # TODO
    end

    for (id, switch) in data_model["switch"]
        # TODO
    end

    if haskey(data_model, "transformer") # transformers are not required by PowerModels
        for (id, trans) in data_model["transformer"]
            # voltage base across transformer does not have to be consistent with the ratio!
            f_vbase = bus_vbase[string(trans["f_bus"])]
            t_vbase = bus_vbase[string(trans["t_bus"])]
            _rebase_pu_transformer_2w_ideal!(trans, f_vbase, t_vbase, sbase_old, sbase, voltage_scale_factor)
        end
    end

    data_math["settings"]["sbase"] = sbase
    data_math["per_unit"] = true
end


"per-unit conversion for buses"
function _rebase_pu_bus!(bus::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real, voltage_scale_factor::Real)
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

        z_old = vbase_old^2*sbase_old*voltage_scale_factor
    end

    # rebase grounding resistance
    z_new = vbase^2/sbase*voltage_scale_factor
    z_scale = z_old/z_new
    _scale_props!(bus, ["rg", "xg"], z_scale)

    if haskey(bus ,"va")
        bus["va"] = deg2rad.(bus["va"])
    end

    # save new vbase
    bus["vbase"] = vbase
end


"per-unit conversion for branches"
function _rebase_pu_branch!(branch::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real, voltage_scale_factor::Real)
    if !haskey(branch, "vbase")
        z_old = 1
    else
        z_old = vbase_old^2/sbase_old*voltage_scale_factor
    end

    z_new = vbase^2/sbase*voltage_scale_factor
    z_scale = z_old/z_new
    y_scale = 1/z_scale
    sbase_scale = sbase_old/sbase

    _scale_props!(branch, ["br_r", "br_x"], z_scale)
    _scale_props!(branch, ["b_fr", "g_fr", "b_to", "g_to"], y_scale)
    _scale_props!(branch, ["c_rating_a", "c_rating_b", "c_rating_c", "rate_a", "rate_b", "rate_c"], sbase_scale)

    branch["angmin"] = deg2rad.(branch["angmin"])
    branch["angmax"] = deg2rad.(branch["angmax"])

    # save new vbase
    branch["vbase"] = vbase
end


"per-unit conversion for shunts"
function _rebase_pu_shunt!(shunt::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real, voltage_scale_factor::Real)
    if !haskey(shunt, "vbase")
        z_old = 1
    else
        z_old = vbase_old^2/sbase_old*voltage_scale_factor
    end

    # rebase grounding resistance
    z_new = vbase^2/sbase*voltage_scale_factor

    z_scale = z_old/z_new
    y_scale = 1/z_scale
    _scale(shunt, "gs", y_scale)
    _scale(shunt, "bs", y_scale)

    # save new vbase
    shunt["vbase"] = vbase
end


"per-unit conversion for loads"
function _rebase_pu_load!(load::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real, voltage_scale_factor::Real)
    if !haskey(load, "vbase")
        vbase_old = 1
        sbase_old = 1
    else
        vbase_old = load["vbase"]
    end

    vbase_old = get(load, "vbase", 1.0)
    vbase_scale = vbase_old/vbase
    _scale(load, "vnom_kv", vbase_scale)

    sbase_scale = sbase_old/sbase
    _scale(load, "pd", sbase_scale)
    _scale(load, "qd", sbase_scale)

    # save new vbase
    load["vbase"] = vbase
end


"per-unit conversion for generators"
function _rebase_pu_generator!(gen::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real, data_model::Dict{String,<:Any})
    vbase_old = get(gen, "vbase", 1.0/data_model["settings"]["voltage_scale_factor"])
    vbase_scale = vbase_old/vbase
    sbase_scale = sbase_old/sbase

    for key in ["pg", "qg", "pmin", "qmin", "pmax", "qmax"]
        _scale(gen, key, sbase_scale)
    end

    # if not in per unit yet, the cost has is in $/MWh
    if !haskey(data_model["settings"], "sbase")
        sbase_old_cost = 1E6/data_model["settings"]["power_scale_factor"]
        sbase_scale_cost = sbase_old_cost/sbase
    else
        sbase_scale_cost = sbase_scale
    end

    _PM._rescale_cost_model!(gen, 1/sbase_scale_cost)

    # save new vbase
    gen["vbase"] = vbase
end


"per-unit conversion for ideal 2-winding transformers"
function _rebase_pu_transformer_2w_ideal!(transformer::Dict{String,<:Any}, f_vbase_new::Real, t_vbase_new::Real, sbase_old::Real, sbase_new::Real, voltage_scale_factor::Real)
    f_vbase_old = get(transformer, "f_vbase", 1.0)
    t_vbase_old = get(transformer, "t_vbase", 1.0)
    f_vbase_scale = f_vbase_old/f_vbase_new
    t_vbase_scale = t_vbase_old/t_vbase_new

    _scale(transformer, "tm_nom", f_vbase_scale/t_vbase_scale)

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


_apply_func_vals(x, f) = isa(x, Dict) ? Dict(k=>f(v) for (k,v) in x) : f.(x)


""
function solution_make_si(solution, math_model; mult_sbase=true, mult_vbase=true, mult_ibase=true, convert_rad2deg=true)
    solution_si = deepcopy(solution)

    sbase = math_model["settings"]["sbase"]

    for (comp_type, comp_dict) in [(x,y) for (x,y) in solution_si if isa(y, Dict)]
        dimensionalize_math_comp = get(_dimensionalize_math, comp_type, Dict())
        vbase_props   = mult_vbase      ? get(dimensionalize_math_comp, "vbase", [])   : []
        sbase_props   = mult_sbase      ? get(dimensionalize_math_comp, "sbase", [])   : []
        ibase_props   = mult_ibase      ? get(dimensionalize_math_comp, "ibase", [])   : []
        rad2deg_props = convert_rad2deg ? get(dimensionalize_math_comp, "rad2deg", []) : []



        for (id, comp) in comp_dict
            if !isempty(vbase_props) || !isempty(ibase_props)
                vbase = math_model[comp_type][id]["vbase"]
                ibase = sbase/vbase
            end

            for (prop, val) in comp
                if prop in vbase_props
                    comp[prop] = _apply_func_vals(comp[prop], x->x*vbase)
                elseif prop in sbase_props
                    comp[prop] = _apply_func_vals(comp[prop], x->x*sbase)
                elseif prop in ibase_props
                    comp[prop] = _apply_func_vals(comp[prop], x->x*ibase)
                elseif prop in rad2deg_props
                    comp[prop] = _apply_func_vals(comp[prop], x->_wrap_to_180(rad2deg(x)))
                end
            end

            if comp_type=="transformer"
                # transformers have different vbase/ibase on each side
                f_bus = math_model["transformer"][id]["f_bus"]
                f_vbase = math_model["bus"]["$f_bus"]["vbase"]
                t_bus = math_model["transformer"][id]["t_bus"]
                t_vbase = math_model["bus"]["$t_bus"]["vbase"]
                f_ibase = sbase/f_vbase
                t_ibase = sbase/t_vbase

                for (prop, val) in comp
                    if prop in dimensionalize_math_comp["ibase_fr"]
                        comp[prop] = _apply_func_vals(comp[prop], x->x*f_ibase)
                    elseif prop in dimensionalize_math_comp["ibase_to"]
                        comp[prop] = _apply_func_vals(comp[prop], x->x*t_ibase)
                    end
                end
            end
        end
    end

    solution_si["per_unit"] = false

    return solution_si
end
