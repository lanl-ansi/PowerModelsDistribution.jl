"lists of scaling factors and what they apply to"
const dimensionalize_math = Dict{String,Dict{String,Vector{String}}}(
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
        "ibase"=>Vector{String}(["cr_fr", "ci_fr", "cr_to", "ci_to", "csr_fr", "csi_fr"])
    ),
    "transformer" => Dict{String,Vector{String}}(
        "sbase"=>Vector{String}(["pf", "qf", "pt", "qt"]),
        "ibase_fr"=>Vector{String}(["cr_fr", "ci_fr"]),
        "ibase_to"=>Vector{String}(["cr_to", "ci_to"])
    ),
    "switch" => Dict{String,Vector{String}}(
        "sbase" => Vector{String}(["pf", "pt", "qf", "qt"]),
        "ibase" => Vector{String}(["cr_fr", "ci_fr", "cr_to", "ci_to"])
    ),
    "storage" => Dict{String,Vector{String}}(
        "sbase"=>Vector{String}(["ps", "qs", "energy", "se", "sd", "sc"]),
    ),
    "shunt" => Dict{String,Vector{String}}(
        "ibase" => Vector{String}(["crsh", "cish"])
    ),
)


"""
    make_per_unit!(
        data::Dict{String,Any};
        vbases::Union{Missing,Dict{String,Real}}=missing,
        sbase::Union{Missing,Real}=missing,
        make_pu_extensions::Vector{<:Function}=Function[],
    )

Converts units of properties to per-unit from SI units

# `make_pu_extensions`

To add additional per-unit transformations, a user can supply custom functions to
`make_pu_extensions::Vector{<:Function}`, which will only be used if `make_pu==true`.

For example, if custom properties are added to the MATHEMATICAL model via
`eng2math_passthrough` or `eng2math_extensions`, those properties will not be converted
to per-unit by default, and custom rules will need to be added with functions with the
signature:

    rebase_pu_func!(
        nw::Dict{String,Any},
        data_math::Dict{String,Any},
        bus_vbase::Dict{String,Real},
        line_vbase::Dict{String,Real},
        sbase::Real,
        sbase_old::Real,
        voltage_scale_factor::Real
    )

where,

- `nw` is equivalent to the a single subnetwork in a multinetwork data structure (which may be the same as `data_math`, in the case of a single network),
- `data_math` is the complete data structure with the global keys,
- `bus_vbase` is a dictionary of the voltage bases of each bus indexed by their MATHEMATICAL model indices,
- `line_vbase` is a dictionary of the voltage bases of each branch indexed by their MATHEMATICAL model indices,
- `sbase` is the new power base,
- `sbase_old` is the power base the data structure started with, and
- `voltage_scale_factor` is the scaling factor for voltage.

"""
function make_per_unit!(
    data::Dict{String,<:Any};
    vbases::Union{Dict{<:Any,<:Real},Missing}=missing,
    sbase::Union{Real,Missing}=missing,
    make_pu_extensions::Vector{<:Function}=Function[],
    )

    data_model_type = get(data, "data_model", MATHEMATICAL)

    if ismath(data)
        if !ismultinetwork(data)
            nw_data = Dict("0" => data)
        else
            nw_data = data["nw"]
        end

        for (n, nw) in nw_data
            vbases = ismissing(vbases) ? nw["settings"]["vbases_default"] : vbases
            sbase  = ismissing(sbase) ? nw["settings"]["sbase_default"] : sbase

            nw["data_model"] = data["data_model"]
            !get(nw, "per_unit", false) && _make_math_per_unit!(nw, data; sbase=sbase, vbases=vbases, make_pu_extensions=make_pu_extensions)
            if ismultinetwork(data)
                delete!(nw, "data_model")
            end
        end
    else
        @warn "Data model '$data_model_type' is not recognized, no per-unit transformation performed"
    end
end


"""
    discover_voltage_zones(data_model::Dict{String,<:Any})::Dict{Int,Set{Any}}

finds voltage zones by walking through the network and analyzing the transformers, attempting to decern the type of `data_model`
"""
function discover_voltage_zones(data_model::Dict{String,<:Any})::Dict{Int,Set{Any}}
    @assert iseng(data_model) || ismath(data_model) "unsupported data model"

    return ismath(data_model) ? discover_math_voltage_zones(data_model) : discover_eng_voltage_zones(data_model)
end


"""
    discover_math_voltage_zones(data_model::Dict{String,<:Any})::Dict{Int,Set{Any}}

finds voltage zones by walking through the network and analyzing the transformers for a MATHEMATICAL `data_model`
"""
discover_math_voltage_zones(data_model::Dict{String,Any})::Dict{Int,Set{Any}} = _discover_voltage_zones(data_model, _math_edge_elements)


"""
    discover_voltage_zones(data_model::Dict{String,<:Any})::Dict{Int,Set{Any}}

finds voltage zones by walking through the network and analyzing the transformers for a ENGINEERING `data_model`
"""
discover_eng_voltage_zones(data_model::Dict{String,Any})::Dict{Int,Set{Any}} = _discover_voltage_zones(data_model, _eng_edge_elements)


"""
    discover_voltage_zones(data_model::Dict{String,<:Any}, edge_elements::Vector{String})::Dict{Int,Set{Any}}

finds voltage zones by walking through the network and analyzing the transformers
"""
function _discover_voltage_zones(data_model::Dict{String,<:Any}, edge_elements::Vector{String})::Dict{Int,Set{Any}}
    unused_components = Set("$comp_type.$id" for comp_type in edge_elements[edge_elements .!= "transformer"] for id in keys(get(data_model, comp_type, Dict())))
    bus_connectors = Dict([(id,Set()) for id in keys(get(data_model, "bus", Dict()))])
    for comp_type in edge_elements[edge_elements .!= "transformer"]
        for (id,obj) in get(data_model, comp_type, Dict())
            f_bus = string(obj["f_bus"])
            t_bus = string(obj["t_bus"])
            push!(bus_connectors[f_bus], ("$comp_type.$id",t_bus))
            push!(bus_connectors[t_bus], ("$comp_type.$id",f_bus))
        end
    end

    zones = []
    buses = Set(keys(get(data_model, "bus", Dict())))
    while !isempty(buses)
        stack = [pop!(buses)]
        zone = Set{Any}()
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
    zones = Dict{Int,Set{Any}}(enumerate(zones))

    return zones
end

"""
    calc_math_voltage_bases(data_model::Dict{String,<:Any}, vbase_sources::Dict{String,<:Real})::Tuple{Dict,Dict}

Calculates voltage bases for each voltage zone for buses and branches for a MATHEMATICAL `data_model`
"""
calc_math_voltage_bases(data_model::Dict{String,<:Any}, vbase_sources::Dict{String,<:Real})::Tuple{Dict,Dict} = _calc_voltage_bases(data_model, vbase_sources, _math_edge_elements)


"""
    calc_eng_voltage_bases(data_model::Dict{String,<:Any}, vbase_sources::Dict{String,<:Real})::Tuple{Dict,Dict}

Calculates voltage bases for each voltage zone for buses and branches for a ENGINEERING `data_model`
"""
calc_eng_voltage_bases(data_model::Dict{String,<:Any}, vbase_sources::Dict{String,<:Real})::Tuple{Dict,Dict} = _calc_voltage_bases(data_model, vbase_sources, _eng_edge_elements)


"""
    calc_voltage_bases(data_model::Dict{String,<:Any}, vbase_sources::Dict{String,<:Real})::Tuple{Dict,Dict}

Calculates voltage bases for each voltage zone for buses and branches, attempting to automatically decern the `data_model` type
"""
function calc_voltage_bases(data_model::Dict{String,<:Any}, vbase_sources::Dict{String,<:Real})::Tuple{Dict,Dict}
    return ismath(data_model) ? calc_math_voltage_bases(data_model, vbase_sources) : calc_eng_voltage_bases(data_model, vbase_sources)
end


"""
    _calc_voltage_bases(data_model::Dict{String,<:Any}, vbase_sources::Dict{String,<:Real}, edge_elements::Vector{String})::Tuple{Dict,Dict}

Calculates voltage bases for each voltage zone for buses and branches given a list of `edge_elements`
"""
function _calc_voltage_bases(data_model::Dict{String,<:Any}, vbase_sources::Dict{String,<:Real}, edge_elements::Vector{String})::Tuple{Dict,Dict}
    # find zones of buses connected by lines
    zones = _discover_voltage_zones(data_model, edge_elements)
    bus_to_zone = Dict([(bus,zone) for (zone, buses) in zones for bus in buses])

    # assign specified vbase to corresponding zones
    zone_vbase = Dict{Int, Union{Missing,Real}}([(zone,missing) for zone in keys(zones)])
    for (bus,vbase) in vbase_sources
        if !ismissing(zone_vbase[bus_to_zone[bus]])
            @warn "You supplied multiple voltage bases for the same zone; ignoring all but the last one."
        end
        zone_vbase[bus_to_zone[bus]] = vbase
    end

    # transformers form the edges between these zones
    zone_edges = Dict{Int,Vector{Tuple{Int,Real}}}([(zone,[]) for zone in keys(zones)])
    for (_,transformer) in get(data_model, "transformer", Dict{Any,Dict{String,Any}}())
        if ismath(data_model)
            f_zone = bus_to_zone["$(transformer["f_bus"])"]
            t_zone = bus_to_zone["$(transformer["t_bus"])"]
            tm_nom = transformer["configuration"]==DELTA ? transformer["tm_nom"]/sqrt(3) : transformer["tm_nom"]
            push!(zone_edges[f_zone], (t_zone, 1/tm_nom))
            push!(zone_edges[t_zone], (f_zone,   tm_nom))
        else
            if !isempty(get(transformer, "xfmrcode", ""))
                _apply_xfmrcode!(transformer, data_model)
            end

            if haskey(transformer, "f_bus")
                f_zone = bus_to_zone["$(transformer["f_bus"])"]
                t_zone = bus_to_zone["$(transformer["f_bus"])"]
                tm_nom = transformer["configuration"] == DELTA ? transformer["tm_nom"]/sqrt(3) : transformer["tm_nom"]
                push!(zone_edges[f_zone], (t_zone, 1/tm_nom))
                push!(zone_edges[t_zone], (f_zone,   tm_nom))
            else
                nrw = length(transformer["bus"])
                f_zone = bus_to_zone[transformer["bus"][1]]
                f_vnom = transformer["configuration"][1] == DELTA ? transformer["vm_nom"][1]/sqrt(3) : transformer["vm_nom"][1]
                for w in 2:nrw
                    t_zone = bus_to_zone[transformer["bus"][w]]
                    t_vnom = transformer["configuration"][1] == DELTA ? transformer["vm_nom"][w]/sqrt(3) : transformer["vm_nom"][w]
                    tm_nom = f_vnom / t_vnom
                    push!(zone_edges[f_zone], (t_zone, 1/tm_nom))
                    push!(zone_edges[t_zone], (f_zone,   tm_nom))
                end
            end
        end
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
    edge_vbase = Dict([("$edge_type.$id", bus_vbase["$(obj["f_bus"])"]) for edge_type in edge_elements[edge_elements .!= "transformer"] if haskey(data_model, edge_type) for (id,obj) in data_model[edge_type]])

    return (bus_vbase, edge_vbase)
end


"converts the MATHEMATICAL model to per unit from SI"
function _make_math_per_unit!(
    nw::Dict{String,<:Any},
    data_math::Dict{String,<:Any};
    sbase::Union{Real,Missing}=missing,
    vbases::Union{Dict{String,<:Real},Missing}=missing,
    make_pu_extensions::Vector{<:Function}=Function[],
    )
    if ismissing(sbase)
        if haskey(nw["settings"], "sbase_default")
            sbase = nw["settings"]["sbase_default"]
        else
            sbase = 1.0
        end
    end

    if haskey(nw["settings"], "sbase")
        sbase_old = nw["settings"]["sbase"]
    else
        sbase_old = 1.0
    end

    # automatically find a good vbase if not provided
    if ismissing(vbases)
        if haskey(nw["settings"], "vbases_default")
            vbases = Dict{String,Real}("$(data_math["bus_lookup"][id])" => vbase for (id, vbase) in nw["settings"]["vbases_default"])
        else
            buses_type_3 = [(id, sum(bus["vm"])/length(bus["vm"])) for (id,bus) in nw["bus"] if haskey(bus, "bus_type") && bus["bus_type"]==3]
            if !isempty(buses_type_3)
                vbases = Dict([buses_type_3[1]])
            else
                error("Please specify vbases manually; cannot make an educated guess for this data model.")
            end
        end
    end

    bus_vbase, line_vbase = calc_voltage_bases(nw, vbases)
    voltage_scale_factor = nw["settings"]["voltage_scale_factor"]

    for (id, bus) in nw["bus"]
        if ismissing(bus_vbase[id])
            error("calc_voltage_bases was unabled to discover the voltage base for bus $id (source_id: $(get(bus, "source_id", ""))), there may be islands in the network model")
        end
        _rebase_pu_bus!(bus, bus_vbase[id], sbase, sbase_old, voltage_scale_factor)
    end

    for (id, line) in nw["branch"]
        vbase = line_vbase["branch.$id"]
        _rebase_pu_branch!(line, vbase, sbase, sbase_old, voltage_scale_factor)
    end

    for (id, shunt) in nw["shunt"]
        _rebase_pu_shunt!(shunt, bus_vbase[string(shunt["shunt_bus"])], sbase, sbase_old, voltage_scale_factor)
    end

    for (id, load) in nw["load"]
        _rebase_pu_load!(load, bus_vbase[string(load["load_bus"])], sbase, sbase_old, voltage_scale_factor)
    end

    for (id, gen) in nw["gen"]
        _rebase_pu_generator!(gen, bus_vbase[string(gen["gen_bus"])], sbase, sbase_old, nw)
    end

    for (id, storage) in nw["storage"]
        _rebase_pu_storage!(storage, bus_vbase[string(storage["storage_bus"])], sbase, sbase_old)
    end

    for (id, switch) in nw["switch"]
        vbase = line_vbase["switch.$id"]
        _rebase_pu_switch!(switch, vbase, sbase, sbase_old, voltage_scale_factor)
    end

    if haskey(nw, "transformer")
        for (id, trans) in nw["transformer"]
            # voltage base across transformer does not have to be consistent with the ratio!
            f_vbase = bus_vbase[string(trans["f_bus"])]
            t_vbase = bus_vbase[string(trans["t_bus"])]
            _rebase_pu_transformer_2w_ideal!(trans, f_vbase, t_vbase, sbase_old, sbase, voltage_scale_factor)
        end
    end

    for rebase_pu_func! in make_pu_extensions
        rebase_pu_func!(nw, data_math, bus_vbase, line_vbase, sbase, sbase_old, voltage_scale_factor)
    end

    nw["settings"]["sbase"] = sbase
    nw["per_unit"] = true
end


"per-unit conversion for buses"
function _rebase_pu_bus!(bus::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real, voltage_scale_factor::Real)
    # if not in p.u., these are normalized with respect to vnom
    prop_vnom = ["vm", "vm_start", "vmax", "vmin", "vm_set", "vm_ln_min", "vm_pn_lb", "vm_pn_ub", "vm_pp_lb", "vm_pp_ub", "vm_ng_ub", "vr_start", "vi_start"]

    if !haskey(bus, "vbase")

        # if haskey(bus, "vm_nom")
        #     vnom = bus["vm_nom"]
        #     _scale_props!(bus, ["vm_nom"], 1/vbase)
        # end
        _scale_props!(bus, prop_vnom, 1/vbase)
        # tupples need special treatment
        for field in ["vm_pair_ub", "vm_pair_lb"]
            # the empty check is needed because otherwise type specialization is lost,
            # the type would then become Vector{Any,Any,Any}
            if !isempty(get(bus, field, []))
                bus[field] = [(c,d,b*1/vbase) for (c,d,b) in bus[field]]
            end
        end

        z_old = 1.0
    else
        vbase_old = bus["vbase"]
        _scale_props!(bus, [prop_vnom..., "vm_nom"], vbase_old/vbase)
        # tupples need special treatment
        for field in ["vm_pair_ub", "vm_pair_lb"]
            if haskey(bus, field)
                bus[field] = [(c,d,b*vbase_old/vbase) for (c,d,b) in bus[field]]
            end
        end

        z_old = vbase_old^2*sbase_old*voltage_scale_factor
    end

    # rebase grounding resistance
    z_new = vbase^2/sbase*voltage_scale_factor
    z_scale = z_old/z_new
    _scale_props!(bus, ["rg", "xg"], z_scale)

    for prop in ["va", "va_start"]
        if haskey(bus, prop)
            bus[prop] = deg2rad.(bus[prop])
        end
    end

    # save new vbase
    bus["vbase"] = vbase
end


"per-unit conversion for branches"
function _rebase_pu_branch!(branch::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real, voltage_scale_factor::Real)
    if !haskey(branch, "vbase")
        z_old = 1
        vbase_old = 1.0
    else
        vbase_old = branch["vbase"]
        z_old = vbase_old^2/sbase_old*voltage_scale_factor
    end

    z_new = vbase^2/sbase*voltage_scale_factor
    z_scale = z_old/z_new
    y_scale = 1/z_scale
    sbase_scale = sbase_old/sbase
    ibase_scale = sbase_scale/(vbase_old/vbase)

    _scale_props!(branch, ["br_r", "br_x"], z_scale)
    _scale_props!(branch, ["b_fr", "g_fr", "b_to", "g_to"], y_scale)
    _scale_props!(branch, ["c_rating_a", "c_rating_b", "c_rating_c"], ibase_scale)
    _scale_props!(branch, ["rate_a", "rate_b", "rate_c"], sbase_scale)

    branch["angmin"] = deg2rad.(branch["angmin"])
    branch["angmax"] = deg2rad.(branch["angmax"])

    # save new vbase
    branch["vbase"] = vbase
end


"per-unit conversion for switches"
function _rebase_pu_switch!(switch::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real, voltage_scale_factor::Real)
    vbase_old = !haskey(switch, "vbase") ? 1.0 : switch["vbase"]

    sbase_scale = sbase_old / sbase
    ibase_scale = sbase_scale/(vbase_old/vbase)

    _scale_props!(switch, ["current_rating", "c_rating_b", "c_rating_c"], ibase_scale)
    _scale_props!(switch, ["thermal_rating", "rate_b", "rate_c"], sbase_scale)

    switch["vbase"] = vbase
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

    # convert capcontrol items to per unit
    if haskey(shunt,"controls")
        if shunt["controls"]["type"] == CAP_REACTIVE_POWER
            power_scale = length(shunt["connections"]) == 1 ? 1 : 3
            shunt["controls"]["onsetting"] = shunt["controls"]["onsetting"] /(sbase*power_scale)
            shunt["controls"]["offsetting"] = shunt["controls"]["offsetting"]/(sbase*power_scale)
            if shunt["controls"]["voltoverride"]
                shunt["controls"]["vmin"] = shunt["controls"]["vmin"]*shunt["controls"]["ptratio"]/(vbase*voltage_scale_factor)
                shunt["controls"]["vmax"] = shunt["controls"]["vmax"]*shunt["controls"]["ptratio"]/(vbase*voltage_scale_factor)
            end
        elseif shunt["controls"]["type"] == CAP_TIME
            # do nothing
        else
            for (idx,val) in enumerate(shunt["controls"]["type"])
                if shunt["controls"]["voltoverride"][idx]
                    shunt["controls"]["vmin"][idx] = shunt["controls"]["vmin"][idx]*shunt["controls"]["ptratio"][idx]/(vbase*voltage_scale_factor)
                    shunt["controls"]["vmax"][idx] = shunt["controls"]["vmax"][idx]*shunt["controls"]["ptratio"][idx]/(vbase*voltage_scale_factor)
                end
                if val == CAP_VOLTAGE
                    shunt["controls"]["onsetting"][idx]  = shunt["controls"]["onsetting"][idx] *shunt["controls"]["ptratio"][idx]/(vbase*voltage_scale_factor)
                    shunt["controls"]["offsetting"][idx] = shunt["controls"]["offsetting"][idx]*shunt["controls"]["ptratio"][idx]/(vbase*voltage_scale_factor)
                elseif val == CAP_CURRENT
                    shunt["controls"]["onsetting"][idx]  = shunt["controls"]["onsetting"][idx] *shunt["controls"]["ctratio"][idx]/(sbase*1e3)*(vbase*voltage_scale_factor)
                    shunt["controls"]["offsetting"][idx] = shunt["controls"]["offsetting"][idx]*shunt["controls"]["ctratio"][idx]/(sbase*1e3)*(vbase*voltage_scale_factor)
                end
            end
        end
    end
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
function _rebase_pu_generator!(gen::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real, nw::Dict{String,<:Any})
    vbase_old = get(gen, "vbase", 1.0/nw["settings"]["voltage_scale_factor"])
    vbase_scale = vbase_old/vbase
    sbase_scale = sbase_old/sbase

    for key in ["pg", "qg", "pmin", "qmin", "pmax", "qmax"]
        _scale(gen, key, sbase_scale)
    end

    for key in ["vg"]
        _scale(gen, key, vbase_scale)
    end

    # if not in per unit yet, the cost has is in $/MWh
    if !haskey(nw["settings"], "sbase")
        sbase_old_cost = 1E6/nw["settings"]["power_scale_factor"]
        sbase_scale_cost = sbase_old_cost/sbase
    else
        sbase_scale_cost = sbase_scale
    end

    _rescale_cost_model!(gen, 1/sbase_scale_cost)

    # save new vbase
    gen["vbase"] = vbase
end


"per-unit conversion for storage"
function _rebase_pu_storage!(gen::Dict{String,<:Any}, vbase::Real, sbase::Real, sbase_old::Real)
    sbase_scale = sbase_old/sbase

    for key in ["energy", "energy_rating", "charge_rating", "discharge_rating", "thermal_rating", "current_rating", "qmin", "qmax", "p_loss", "q_loss", "ps", "qs"]
        _scale(gen, key, sbase_scale)
    end
end


"per-unit conversion for ideal 2-winding transformers"
function _rebase_pu_transformer_2w_ideal!(transformer::Dict{String,<:Any}, f_vbase_new::Real, t_vbase_new::Real, sbase_old::Real, sbase_new::Real, voltage_scale_factor::Real)
    f_vbase_old = get(transformer, "f_vbase", 1.0)
    t_vbase_old = get(transformer, "t_vbase", 1.0)
    f_vbase_scale = f_vbase_old/f_vbase_new
    t_vbase_scale = t_vbase_old/t_vbase_new
    sbase_scale = sbase_old/sbase_new

    _scale(transformer, "tm_nom", f_vbase_scale/t_vbase_scale)
    _scale(transformer, "sm_ub", sbase_scale)

    # save new vbase
    transformer["f_vbase"] = f_vbase_new
    transformer["t_vbase"] = t_vbase_new

    # convert regcontrol items to per unit
    if haskey(transformer,"controls")
        # convert reference voltage and band from volts to per unit
        transformer["controls"]["vreg"] = transformer["controls"]["vreg"].*transformer["controls"]["ptratio"]/(f_vbase_new*voltage_scale_factor)
        transformer["controls"]["band"] = transformer["controls"]["band"].*transformer["controls"]["ptratio"]/(f_vbase_new*voltage_scale_factor)
        # convert regulator impedance from volts to per unit
        baseZ = (f_vbase_new*voltage_scale_factor)^2/(sbase_new*1e3)
        transformer["controls"]["r"] = transformer["controls"]["r"].*transformer["controls"]["ptratio"]./transformer["controls"]["ctprim"]/baseZ
        transformer["controls"]["x"] = transformer["controls"]["x"].*transformer["controls"]["ptratio"]./transformer["controls"]["ctprim"]/baseZ
    end
end


"helper function to apply a scale factor to given properties"
function _scale_props!(comp::Dict{String,<:Any}, prop_names::Vector{String}, scale::Real)
    for name in prop_names
        if haskey(comp, name)
            comp[name] *= scale
        end
    end
end


"helper to apply function values"
_apply_func_vals(x, f) = isa(x, Dict) ? Dict(k=>f(v) for (k,v) in x) : f.(x)


"""
    solution_make_si(
        solution::Dict{String,<:Any},
        math_model::Dict{String,<:Any};
        mult_sbase::Bool=true,
        mult_vbase::Bool=true,
        mult_ibase::Bool=true,
        convert_rad2deg::Bool=true,
        make_si_extensions::Vector{<:Function}=Function[],
        dimensionalize_math_extensions::Dict{String,<:Dict{String,<:Vector{<:String}}}=Dict{String,Dict{String,Vector{String}}}()
    )::Dict{String,Any}

Transforms solution dictionaries `solution` from per-unit back to SI units, requiring the original MATHEMATICAL model `math_model`
to perform the transformation.

If `mult_sbase` is false, sbase variables will not be multiplied, thus remaining in per-unit

If `mult_vbase` is false, vbase variables will not be multiplied, thus remaining in per-unit

If `mult_ibase` is false, ibase variables will not be multiplied, thus remaining in per-unit

If `convert_rad2deg` is false, angle variables will not be multiplied, thus remaining in radians

# Custom SI unit conversions

To convert custom properties not part of formulations already included within PowerModelsDistribution,
users will need to either specify multiplicative factors via `dimensionalize_math_extensions`, or pass
user functions via `make_si_extensions`.

The latter case requires functions with the signature

    make_si_func!(nw_solution, nw_data, solution, data)

where `nw_solution` and `nw_data` are equivalent to a single subnetwork of a multinetwork structure of
the solution and the data in the MATHEMATICAL format, respectively, and `solution` and `data` are the
full data structures, which may be equivalent to `nw_solution` and `nw_data`, if the data is not
multinetwork. Changes should be applied to `nw_solution` in the user functions.

For `dimensionalize_math_extensions`, it is possible to easily extended the SI conversions if they are
straightforward conversions using `vbase`, `sbase`, `ibase`, or `rad2deg`. For example, if a custom
variable `cfr` is added to branches, and is scaled by `ibase`, the following dictionary would be passed:

    Dict{String,Dict{String,Vector{String}}}(
        "branch" => Dict{String,Vector{String}}(
            "ibase" => String["cfr"]
        )
    )

which would ensure that this variable gets converted back to SI units upon transformation.
"""
function solution_make_si(
    solution::Dict{String,<:Any},
    math_model::Dict{String,<:Any};
    mult_sbase::Bool=true,
    mult_vbase::Bool=true,
    mult_ibase::Bool=true,
    convert_rad2deg::Bool=true,
    make_si_extensions::Vector{<:Function}=Function[],
    dimensionalize_math_extensions::Dict{String,<:Dict{String,<:Vector{<:String}}}=Dict{String,Dict{String,Vector{String}}}()
    )::Dict{String,Any}

    solution_si = deepcopy(solution)

    if ismultinetwork(math_model)
        nw_data = math_model["nw"]
        nw_sol = get(solution_si, "nw", Dict{String,Any}())  # in case solution is not found
    else
        nw_data = Dict("0" => math_model)
        nw_sol = Dict("0" => solution_si)
    end

    for (n,nw) in nw_sol
        if !isempty(nw) && nw["per_unit"]
            sbase = nw["settings"]["sbase"]
            for (comp_type, comp_dict) in [(x,y) for (x,y) in nw if isa(y, Dict) && x != "settings"]
                dimensionalize_math_comp = get(dimensionalize_math, comp_type, Dict())
                ext_comp = get(dimensionalize_math_extensions, comp_type, Dict())

                vbase_props   = mult_vbase      ? [get(dimensionalize_math_comp, "vbase", []); get(ext_comp, "vbase", [])]   : []
                sbase_props   = mult_sbase      ? [get(dimensionalize_math_comp, "sbase", []); get(ext_comp, "sbase", [])]   : []
                ibase_props   = mult_ibase      ? [get(dimensionalize_math_comp, "ibase", []); get(ext_comp, "ibase", [])]   : []
                rad2deg_props = convert_rad2deg ? [get(dimensionalize_math_comp, "rad2deg", []); get(ext_comp, "rad2deg", [])] : []

                for (id, comp) in comp_dict
                    if !isempty(vbase_props) || !isempty(ibase_props)
                        vbase = nw_data[n][comp_type][id]["vbase"]
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
                        f_bus = nw_data[n]["transformer"][id]["f_bus"]
                        f_vbase = nw_data[n]["bus"]["$f_bus"]["vbase"]
                        t_bus = nw_data[n]["transformer"][id]["t_bus"]
                        t_vbase = nw_data[n]["bus"]["$t_bus"]["vbase"]
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
            nw["per_unit"] = false
        end

        for make_si_func! in make_si_extensions
            make_si_func!(nw, nw_data[n], solution, math_model)
        end
    end

    return solution_si
end
