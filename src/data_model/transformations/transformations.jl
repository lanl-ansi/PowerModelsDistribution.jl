# This file contains useful transformation functions for the engineering data model

"collect of components and their properties that define loss models when transforming to the MATHEMATICAL model"
const _loss_model_objects = Dict{String,Vector{String}}(
    "linecode" => String["rs", "xs", "g_fr", "g_to", "b_fr", "b_to"],
    "line" => String["rs", "xs", "g_fr", "g_to", "b_fr", "b_to"],
    "switch" => String["rs", "xs", "g_fr", "g_to", "b_fr", "b_to"],
    "xfmrcode" => String["rw", "xsc", "cmag", "noloadloss"],
    "transformer" => String["rw", "xsc", "cmag", "noloadloss"],
    "voltage_source" => String["rs", "xs"],
    "generator" => String["rs", "xs"],
    "solar" => String["rs","xs"],
    "storage" => String["rs", "xs", "pex", "qex"],
)


"""
    make_lossless!(data_eng::Dict{String,<:Any})

Remove parameters from objects with loss models to make them lossless. This includes linecodes,
lines, switches, xfmrcodes, transformers, voltage sources, generators, solar, and storage, which
all have (or will have in the future), loss model parameters that can be omitted.
"""
function make_lossless!(data::Dict{String,<:Any}; exclude::Vector{String}=String[])
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_make_lossless!, data; apply_to_subnetworks=true, exclude=exclude)
end


"""
    _make_lossless!(data_eng::Dict{String,<:Any})

Remove parameters from objects with loss models to make them lossless. This includes linecodes,
lines, switches, xfmrcodes, transformers, voltage sources, generators, solar, and storage, which
all have (or will have in the future), loss model parameters that can be omitted.
"""
function _make_lossless!(data_eng::Dict{String,<:Any}; exclude::Vector{String}=String[])
    for (object_type, parameters) in _loss_model_objects
        if haskey(data_eng, object_type) && !(object_type in exclude)
            for (id, eng_obj) in data_eng[object_type]
                if object_type == "switch"
                    delete!(eng_obj, "linecode")
                end
                for parameter in parameters
                    if haskey(eng_obj, parameter)
                        eng_obj[parameter] = 0.0 .* eng_obj[parameter]
                    end
                end
            end
        end
    end
end


"""
    apply_phase_projection!(data::Dict{String,<:Any})

pad matrices and vectors to max number of conductors
"""
function apply_phase_projection!(data::Dict{String,<:Any})
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_apply_phase_projection!, data; apply_to_subnetworks=true)
end


"""
    _apply_phase_projection!(data_eng::Dict{String,<:Any})

pad matrices and vectors to max number of conductors
"""
function _apply_phase_projection!(data_eng::Dict{String,<:Any})
    @assert get(data_eng, "is_kron_reduced", false)

    all_conductors = _get_complete_conductor_set(data_eng)

    if haskey(data_eng, "bus")
        for (_,eng_obj) in data_eng["bus"]
            if !haskey(eng_obj, "vm_lb")
                eng_obj["vm_lb"] = fill(0.0, length(eng_obj["terminals"]))
            end

            if !haskey(eng_obj, "vm_ub")
                eng_obj["vm_ub"] = fill(Inf, length(eng_obj["terminals"]))
            end

            _pad_properties!(eng_obj, ["vm", "va", "vm_lb", "vm_ub", "vm_start", "va_start"], eng_obj["terminals"], all_conductors)

            _pad_connections!(eng_obj, "terminals", all_conductors)
        end
    end

    if haskey(data_eng, "line")
        for (_,eng_obj) in data_eng["line"]
            _apply_linecode!(eng_obj, data_eng)

            _pad_properties!(eng_obj, ["rs", "xs", "g_fr", "g_to", "b_fr", "b_to"], eng_obj["f_connections"], all_conductors)
            _pad_properties!(eng_obj, ["vad_lb"], eng_obj["f_connections"], all_conductors; pad_value=-10.0)
            _pad_properties!(eng_obj, ["vad_ub"], eng_obj["f_connections"], all_conductors; pad_value=10.0)
            for key in ["cm_ub", "cm_ub_b", "cm_ub_c", "sm_ub", "sm_ub_b", "sm_ub_c"]
                if haskey(eng_obj, key)
                    _pad_properties!(eng_obj, [key], eng_obj["f_connections"], all_conductors)
                end
            end

            _pad_connections!(eng_obj, "f_connections", all_conductors)
            _pad_connections!(eng_obj, "t_connections", all_conductors)
        end
    end

    if haskey(data_eng, "transformer")
        for (_,eng_obj) in data_eng["transformer"]
            _apply_xfmrcode!(eng_obj, data_eng)

            if haskey(eng_obj, "f_connections")
                _pad_properties!(eng_obj, ["tm_lb", "tm_ub", "tm_set"], eng_obj["f_connections"], all_conductors; pad_value=1.0)
                _pad_properties!(eng_obj, ["tm_fix"], eng_obj["f_connections"], all_conductors; pad_value=true)

                _pad_connections!(eng_obj, "f_connections", all_conductors)
                _pad_connections!(eng_obj, "t_connections", all_conductors)
            else
                for (w, connections) in enumerate(eng_obj["connections"])
                    _pad_properties!(eng_obj, ["tm_lb", "tm_ub", "tm_set"], w, connections, all_conductors; pad_value=1.0)
                    _pad_properties!(eng_obj, ["tm_fix"], w, connections, all_conductors; pad_value=true)

                    _pad_connections!(eng_obj, "connections", w, all_conductors)
                end
            end
        end
    end

    if haskey(data_eng, "switch")
        for (_,eng_obj) in data_eng["switch"]
            _apply_linecode!(eng_obj, data_eng)

            _pad_properties!(eng_obj, ["rs", "xs"], eng_obj["f_connections"], all_conductors)
            _pad_properties!(eng_obj, ["vad_lb"], eng_obj["f_connections"], all_conductors; pad_value=-10.0)
            _pad_properties!(eng_obj, ["vad_ub"], eng_obj["f_connections"], all_conductors; pad_value=10.0)
            for key in ["cm_ub", "cm_ub_b", "cm_ub_c", "sm_ub", "sm_ub_b", "sm_ub_c"]
                if haskey(eng_obj, key)
                    _pad_properties!(eng_obj, [key], eng_obj["f_connections"], all_conductors)
                end
            end

            _pad_connections!(eng_obj, "f_connections", all_conductors)
            _pad_connections!(eng_obj, "t_connections", all_conductors)
        end
    end

    if haskey(data_eng, "shunt")
        for (_,eng_obj) in data_eng["shunt"]
            _pad_properties!(eng_obj, ["gs", "bs"], eng_obj["connections"], all_conductors)

            _pad_connections!(eng_obj, "connections", all_conductors)
        end
    end

    if haskey(data_eng, "load")
        for (_,eng_obj) in data_eng["load"]
            if eng_obj["configuration"] == WYE
                _pad_properties!(eng_obj, ["pd_nom", "qd_nom"], eng_obj["connections"], all_conductors)
            else
                _pad_properties_delta!(eng_obj, ["pd_nom", "qd_nom"], eng_obj["connections"], all_conductors)
            end

            _pad_connections!(eng_obj, "connections", all_conductors)
        end
    end

    if haskey(data_eng, "generator")
        for (_,eng_obj) in data_eng["generator"]
            if eng_obj["configuration"]==WYE
                _pad_properties!(eng_obj, ["pg", "qg", "vg", "pg_lb", "pg_ub", "qg_lb", "qg_ub"], eng_obj["connections"], all_conductors)
            else
                _pad_properties_delta!(eng_obj, ["pg", "qg", "vg", "pg_lb", "pg_ub", "qg_lb", "qg_ub"], eng_obj["connections"], all_conductors)
            end

            _pad_connections!(eng_obj, "connections", all_conductors)
        end
    end

    if haskey(data_eng, "solar")
        for (_,eng_obj) in data_eng["solar"]
            if eng_obj["configuration"]==WYE
                _pad_properties!(eng_obj, ["pg", "qg", "vg", "pg_lb", "pg_ub", "qg_lb", "qg_ub"], eng_obj["connections"], all_conductors)
            else
                _pad_properties_delta!(eng_obj, ["pg", "qg", "vg", "pg_lb", "pg_ub", "qg_lb", "qg_ub"], eng_obj["connections"], all_conductors)
            end

            _pad_connections!(eng_obj, "connections", all_conductors)
        end
    end

    if haskey(data_eng, "storage")
        for (_,eng_obj) in data_eng["storage"]
            _pad_properties!(eng_obj, ["cm_ub", "qs_lb", "qs_ub", "rs", "xs", "ps", "qs"], eng_obj["connections"], all_conductors)

            _pad_connections!(eng_obj, "connections", all_conductors)
        end
    end

    if haskey(data_eng, "voltage_source")
        for (_,eng_obj) in data_eng["voltage_source"]
            _pad_properties!(eng_obj, ["pg", "qg", "pg_lb", "pg_ub", "qg_lb", "qg_ub", "vm", "va"], eng_obj["connections"], all_conductors)

            _pad_connections!(eng_obj, "connections", all_conductors)
        end
    end

    data_eng["is_projected"] = true
end


"""
    apply_phase_projection_delta!(data::Dict{String,<:Any})

phase projection for components where unprojected states are not yet supported (delta configurations).

See [`apply_phase_projection!`](@ref apply_phase_projection!)
"""
function apply_phase_projection_delta!(data::Dict{String,<:Any})
    apply_pmd!(_apply_phase_projection_delta!, data; apply_to_subnetworks=true)
end


"""
    _apply_phase_projection_delta!(data_eng::Dict{String,<:Any})

phase projection for components where unprojected states are not yet supported (delta configurations).

See [`apply_phase_projection!`](@ref apply_phase_projection!)
"""
function _apply_phase_projection_delta!(data_eng::Dict{String,<:Any})
    @assert get(data_eng, "is_kron_reduced", false)

    bus_terminals = Dict{String,Vector{Int}}()

    all_conductors = _get_complete_conductor_set(data_eng)

    if haskey(data_eng, "transformer")
        for (_,eng_obj) in data_eng["transformer"]
            _apply_xfmrcode!(eng_obj, data_eng)

            if haskey(eng_obj, "f_connections")
                if eng_obj["configuration"] == DELTA
                    _pad_properties!(eng_obj, ["tm_lb", "tm_ub", "tm_set"], eng_obj["f_connections"], all_conductors; pad_value=1.0)
                    _pad_properties!(eng_obj, ["tm_fix"], eng_obj["f_connections"], all_conductors; pad_value=true)

                    _pad_connections!(eng_obj, "f_connections", all_conductors)
                    _pad_connections!(eng_obj, "t_connections", all_conductors)

                    bus_terminals[eng_obj["f_bus"]] = haskey(bus_terminals, eng_obj["f_bus"]) ? _pad_connections!(bus_terminals, eng_obj["f_bus"], eng_obj["f_connections"]) : eng_obj["f_connections"]
                    bus_terminals[eng_obj["t_bus"]] = haskey(bus_terminals, eng_obj["t_bus"]) ? _pad_connections!(bus_terminals, eng_obj["t_bus"], eng_obj["t_connections"]) : eng_obj["t_connections"]
                end
            else
                for (w, connections) in enumerate(eng_obj["connections"])
                    if eng_obj["configuration"][w] == DELTA && length(eng_obj["tm_set"][w]) == 1
                        _pad_properties_delta!(eng_obj, ["tm_lb", "tm_ub", "tm_set"], connections, w, all_conductors)
                        _pad_properties_delta!(eng_obj, ["tm_fix"], connections, w, all_conductors)
                        for i in 1:length(eng_obj["tm_set"][w])
                            eng_obj["tm_set"][w][i] = eng_obj["tm_set"][w][i] == 0 ? 1 : eng_obj["tm_set"][w][i]
                            eng_obj["tm_fix"][w][i] = eng_obj["tm_fix"][w][i] == 0 ? true : eng_obj["tm_fix"][w][i]
                        end
                        _pad_connections!(eng_obj, "connections", w, all_conductors)
                        bus_terminals[eng_obj["bus"][w]] = haskey(bus_terminals, eng_obj["bus"][w]) ? _pad_connections!(bus_terminals, eng_obj["bus"][w], eng_obj["connections"][w]) : eng_obj["connections"][w]
                    end
                end
            end
        end
    end

    if haskey(data_eng, "load")
        for (_,eng_obj) in data_eng["load"]
            if eng_obj["configuration"] == DELTA
                if eng_obj["connections"]==[1, 2] && eng_obj["vm_nom"]==0.24 # check if load is connected between split-phase terminals of triplex node (nominal line-line voltage=240V), TODO: better generalization
                    bus_terminals[eng_obj["bus"]] = eng_obj["connections"] = [1]
                else
                    _pad_properties_delta!(eng_obj, ["pd_nom", "qd_nom"], eng_obj["connections"], all_conductors)
                    _pad_connections!(eng_obj, "connections", all_conductors)
                    bus_terminals[eng_obj["bus"]] = haskey(bus_terminals, eng_obj["bus"]) ? _pad_connections!(bus_terminals, eng_obj["bus"], eng_obj["connections"]) : eng_obj["connections"]
                end
            end
        end
    end

    if haskey(data_eng, "generator")
        for (_,eng_obj) in data_eng["generator"]
            if eng_obj["configuration"]==DELTA
                if eng_obj["connections"]==[1, 2] && eng_obj["vg"][1]==0.24 # check if generator is connected between split-phase terminals of triplex node (nominal line-line voltage=240V), TODO: better generalization
                    bus_terminals[eng_obj["bus"]] = eng_obj["connections"] = [1]
                else
                    _pad_properties_delta!(eng_obj, ["pg", "qg", "vg", "pg_lb", "pg_ub", "qg_lb", "qg_ub"], eng_obj["connections"], all_conductors)
                    _pad_connections!(eng_obj, "connections", all_conductors)
                    bus_terminals[eng_obj["bus"]] = haskey(bus_terminals, eng_obj["bus"]) ? _pad_connections!(bus_terminals, eng_obj["bus"], eng_obj["connections"]) : eng_obj["connections"]
                end
            end
        end
    end

    if haskey(data_eng, "solar")
        for (_,eng_obj) in data_eng["solar"]
            if eng_obj["configuration"]==DELTA
                if eng_obj["connections"]==[1, 2] && eng_obj["vg"][1]==0.24 # check if solar is connected between split-phase terminals of triplex node (nominal line-line voltage=240V), TODO: better generalization
                    bus_terminals[eng_obj["bus"]] = eng_obj["connections"] = [1]
                else
                    _pad_properties_delta!(eng_obj, ["pg", "qg", "vg", "pg_lb", "pg_ub", "qg_lb", "qg_ub"], eng_obj["connections"], all_conductors)
                    _pad_connections!(eng_obj, "connections", all_conductors)
                    bus_terminals[eng_obj["bus"]] = haskey(bus_terminals, eng_obj["bus"]) ? _pad_connections!(bus_terminals, eng_obj["bus"], eng_obj["connections"]) : eng_obj["connections"]
                end
            end
        end
    end

    _update_bus_terminal_projections!(data_eng, bus_terminals)

    data_eng["is_projected"] = true
end


"helper function to update the terminals on projected buses"
function _update_bus_terminal_projections!(data_eng::Dict{String,<:Any}, bus_terminals::Dict{String,<:Vector{Int}})
    for (id,terminals) in bus_terminals
        eng_obj = data_eng["bus"][id]

        if !haskey(eng_obj, "vm_lb")
            eng_obj["vm_lb"] = fill(0.0, length(eng_obj["terminals"]))
        end

        if !haskey(eng_obj, "vm_ub")
            eng_obj["vm_ub"] = fill(Inf, length(eng_obj["terminals"]))
        end

        old_terms = deepcopy(eng_obj["terminals"])
        new_terms = _pad_connections!(eng_obj, "terminals", terminals)
        _pad_properties!(eng_obj, ["vm", "va", "vm_lb", "vm_ub", "vm_start", "va_start"], old_terms, new_terms)
    end
end


"Indicates whether the passed component has a oneport structure (e.g. loads and generators)."
_is_oneport_component(comp::Dict{String,Any})::Bool = haskey(comp, "bus") && isa(comp["bus"], AbstractString)


"Indicates whether the passed component has a twoport structure (e.g. lines and switches)."
_is_twoport_component(comp::Dict{String,Any})::Bool = haskey(comp, "f_bus") && haskey(comp, "t_bus")


"Indicates whether the passed component has a multiport structure (e.g. transformers)."
_is_multiport_component(comp::Dict{String,Any})::Bool = haskey(comp, "bus") && isa(comp["bus"], Vector{<:AbstractString})


"Obtain impedance parameters, directly or from linecode."
function _get_line_impedance_parameters(data_eng::Dict{String,Any}, line::Dict{String,Any})
    @assert data_eng["data_model"]==ENGINEERING

    if haskey(line, "rs")
        z_s = line["rs"].+im*line["xs"]
        y_fr = line["g_fr"].+im*line["b_fr"]
        y_to = line["g_to"].+im*line["b_to"]
    else
        lc = data_eng["linecode"][line["linecode"]]
        z_s = (lc["rs"].+im*lc["xs"])*line["length"]
        y_fr = (lc["g_fr"].+im*lc["b_fr"])*line["length"]
        y_to = (lc["g_to"].+im*lc["b_to"])*line["length"]
    end

    return z_s, y_fr, y_to
end


"Create an equivalent shunt for a line which connects to a single bus."
function _loop_line_to_shunt(data_eng::Dict{String,Any}, line_id::AbstractString)
    @assert data_eng["data_model"]==ENGINEERING

    # only possible when the line is a 'loop' with respect to its bus
    line = data_eng["line"][line_id]
    @assert line["f_bus"]==line["t_bus"]

    # obtain impedance parameters, directly or from linecode
    z_s, y_fr, y_to = _get_line_impedance_parameters(data_eng, line)
    y_s = inv(z_s)

    # build shunt addmittance
    Yb = [y_s+y_fr -y_s; -y_s y_s+y_to]
    conns = [line["f_connections"]..., line["t_connections"]...]

    # simplify to a unique set of of connections and equivalent addmittance
    conns_unique = unique(conns)
    M = [conns[i]==conns_unique[j] ? 1 : 0 for i in 1:length(conns), j in 1:length(conns_unique)]
    Yb_unique = M'*Yb*M

    # build shunt dict
    shunt =  Dict{String, Any}(
        "status" => line["status"],
        "dispatchable" => NO,
        "bus" => line["f_bus"],
        "connections" => conns_unique,
        "gs" => real.(Yb_unique),
        "bs" => imag.(Yb_unique),
    )
    if haskey(line, "source_id")
        shunt["source_id"] = line["source_id"]
    end

    return shunt
end


"Merge a terminal into another for a specified bus, i.e. as if they are short-ciruited."
function _merge_terminals!(data_eng::Dict{String,Any}, bus_id::String, t_fr, t_to)
    @assert data_eng["data_model"]==ENGINEERING

    bus = data_eng["bus"][bus_id]
    old_terminals = bus["terminals"]
    # exclude t_fr from the bus terminals
    merged_terminals = bus["terminals"] = [x for x in old_terminals if x!=t_fr]
    # find position of merged terminals in old terminal vector
    old_idxs = Dict(t=>findall(t.==old_terminals) for t in merged_terminals)
    # assign idxs of t_fr to idxs of t_to
    append!(old_idxs[t_to], findall(t_fr.==old_terminals))


    # resolve 'vmin' and 'vmax' properties
    if haskey(bus, "vmin")
        bus["vmin"] = [maximum(bus["vmin"][old_idxs[t]]) for t in merged_terminals]
    end
    if haskey(bus, "vmax")
        bus["vmax"] = [minimum(bus["vmax"][old_idxs[t]]) for t in merged_terminals]
    end

    # resolve properties which should be the same for the merged terminals
    for prop in ["vm", "va", "vr_start", "vi_start", "vm_start", "va_start"]
        if haskey(bus, prop)
            vals = bus[prop][old_idx[t_to]]
            @assert vals.==vals[1] "Cannot merge bus property $prop because the merged terminals differ."
            bus[prop] = [bus[prop][findfirst(old_terminals.==t)] for t in merged_terminals]
        end
    end

    # merge grounding properties (i.e. equivalent impedance of parallel impedances to ground if multiple groundings)
    if haskey(bus, "grounded")
        grounded = bus["grounded"]
        if t_fr in grounded
            idxs_fr = findall(grounded.==t_fr)
            idxs_to = findall(grounded.==t_to)
            idxs = [idxs_fr..., idxs_to...]
            zgs = [bus["rg"][i]+im*bus["xg"][i] for i in idxs]
            zg = any(iszero.(zgs)) ? im*0.0 : 1/sum(1/x for x in zgs)

            idxs_other = setdiff(1:length(grounded), idxs)
            bus["grounded"] = [grounded[idxs_other]..., t_to]
            bus["rg"] = [bus["rg"][idxs_other]..., real.(zg)]
            bus["xg"] = [bus["xg"][idxs_other]..., imag.(zg)]
        end
    end

    #TODO take care of other bounds in the bus data model

    # update connection properties to reflect the merge
    for type in setdiff(intersect(pmd_eng_asset_types, keys(data_eng)), ["bus"])
        for (_,comp) in data_eng[type]
            # one-port components
            if _is_oneport_component(comp) && comp["bus"]==bus_id
                comp["connections"] = [t==t_fr ? t_to : t for t in comp["connections"]]
            # two-port components
            elseif _is_twoport_component(comp)
                if comp["f_bus"]==bus_id
                    comp["f_connections"] = [t==t_fr ? t_to : t for t in comp["f_connections"]]
                end
                if comp["t_bus"]==bus_id
                    comp["t_connections"] = [t==t_fr ? t_to : t for t in comp["t_connections"]]
                end
            # multi-port components
            elseif _is_multiport_component(comp)
                for (w,bus_w) in enumerate(comp["bus"])
                    if bus_w==bus_id
                        comp["connections"][w] = [t==t_fr ? t_to : t for t in comp["connections"][w]]
                    end
                end
            end
        end
    end
end


"""
    transform_loops!(
        data::Dict{String,Any};
        zero_series_impedance_threshold::Real=1E-8,
        shunt_id_prefix::AbstractString="line_loop"
    )::Dict{String,Any}

Transform line loops (connected to a single bus), which are not allowed in the mathematical model.
Lossy line loops are converted to equivalent shunts, and lossless ones (i.e. short-circuits) are represented by merging the short-circuited terminals.
The argument 'zero_series_impedance_threshold' controls the threshold below which the series impedance is considered to be a short-ciruit.
This is useful because OpenDSS modelers have to insert tiny impedances to represent short-circuit reactors.
The addmittance to ground should be zero to trigger the short-circuit handling.
"""
function transform_loops!(
    data::Dict{String,Any};
    zero_series_impedance_threshold::Real=1E-8,
    shunt_id_prefix::AbstractString="line_loop"
    )::Dict{String,Any}

    @assert iseng(data) "wrong data model type"

    apply_pmd!(_transform_loops!, data; apply_to_subnetworks=true, zero_series_impedance_threshold=zero_series_impedance_threshold, shunt_id_prefix=shunt_id_prefix)
end


"""
    _transform_loops!(
        data_eng::Dict{String,Any};
        zero_series_impedance_threshold::Real=1E-8,
        shunt_id_prefix::AbstractString="line_loop"
    )::Dict{String,Any}

Transform line loops (connected to a single bus), which are not allowed in the mathematical model.
Lossy line loops are converted to equivalent shunts, and lossless ones (i.e. short-circuits) are represented by merging the short-circuited terminals.
The argument 'zero_series_impedance_threshold' controls the threshold below which the series impedance is considered to be a short-ciruit.
This is useful because OpenDSS modelers have to insert tiny impedances to represent short-circuit reactors.
The addmittance to ground should be zero to trigger the short-circuit handling.
"""
function _transform_loops!(
    data_eng::Dict{String,Any};
    zero_series_impedance_threshold::Real=1E-8,
    shunt_id_prefix::AbstractString="line_loop"
    )::Dict{String,Any}

    for (id,line) in get(data_eng, "line", Dict())
        if line["f_bus"]==line["t_bus"]
            # obtain impedance parameters, directly or from linecode
            z_s, y_fr, y_to = _get_line_impedance_parameters(data_eng, line)

            # remove short-circuit line and merge terminals
            if all(iszero.(y_fr)) && all(iszero.(y_to)) && all(abs.(z_s).<=zero_series_impedance_threshold)
                for (t_fr, t_to) in zip(line["f_connections"], line["t_connections"])
                    _merge_terminals!(data_eng, line["f_bus"], sort([t_fr, t_to], rev=true)...)
                end
            # convert line to a shunt
            else
                shunt = _loop_line_to_shunt(data_eng, id)
                if !haskey(data_eng, "shunt")
                    data_eng["shunt"] = Dict{String, Any}()
                end
                data_eng["shunt"]["$shunt_id_prefix.$id"] = shunt
            end
            # delete the line now that it has been handled
            delete!(data_eng["line"], id)
        end
    end

    # update the 'conductor_ids' property
    find_conductor_ids!(data_eng)

    return data_eng
end


"""
    remove_unconnected_terminals!(data::Dict{String,Any})::Dict{String,Any}

Remove all terminals which are unconnected (not considering a grounding as a connection).
"""
function remove_unconnected_terminals!(data::Dict{String,Any})::Dict{String,Any}
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_remove_unconnected_terminals!, data; apply_to_subnetworks=true)
end


"""
    _remove_unconnected_terminals!(data_eng::Dict{String,Any})::Dict{String,Any}

Remove all terminals which are unconnected (not considering a grounding as a connection).
"""
function _remove_unconnected_terminals!(data_eng::Dict{String,Any})::Dict{String,Any}
    # find all connected bts (a 'bt' is a bus-terminal pair)
    connected_bts = []
    for type in setdiff(intersect(pmd_eng_asset_types, keys(data_eng)), ["bus"])
        for (id,comp) in data_eng[type]
            if _is_oneport_component(comp)
                connected_bts = append!(connected_bts, [(comp["bus"],t) for t in comp["connections"]])
            elseif _is_twoport_component(comp)
                connected_bts = append!(connected_bts, [(comp["f_bus"],t) for t in comp["f_connections"]])
                connected_bts = append!(connected_bts, [(comp["t_bus"],t) for t in comp["t_connections"]])
            elseif _is_multiport_component(comp)
                for (b_w,conns_w) in zip(comp["bus"], comp["connections"])
                    connected_bts = append!(connected_bts, [(b_w,t) for t in conns_w])
                end
            end
        end
    end
    connected_bts = unique(connected_bts)

    # remove all unconnected bts
    bts = [(b,t) for (b,bus) in data_eng["bus"] for t in bus["terminals"]]
    for (b,t) in setdiff(bts, connected_bts)
        bus = data_eng["bus"][b]
        keep_ts_order = (!).(bus["terminals"].==t)
        _apply_filter!(bus, ["terminals", "vm", "va", "vmin", "vmax"], keep_ts_order)
        if t in bus["grounded"]
            keep_gr_order = (!).(bus["grounded"].==t)
            _apply_filter!(bus, ["grounded", "rg", "xg"], keep_gr_order)
        end
        #TODO other bounds
    end

    return data_eng
end


"""
    reduce_lines!(data_eng::Dict{String,Any})::Dict{String,Any}

Reduces line models, by removing trailing lines and merging lines in series with the same linecode.
"""
function reduce_lines!(data_eng::Dict{String,Any})::Dict{String,Any}
    delete_trailing_lines!(data_eng)
    join_lines!(data_eng)

    return data_eng
end


"Returns a reduced data model, by removing trailing lines and merging lines in series with the same linecode."
reduce_lines(data_eng::Dict{String,Any}) = reduce_lines!(deepcopy(data_eng))


"Reverse the direction of a line."
function _line_reverse!(line::Dict{String,Any})
    prop_pairs = [("f_bus", "t_bus"), ("f_connections", "t_connections"), ("g_fr", "g_to"), ("b_fr","b_to")]

    for (x,y) in prop_pairs
        if haskey(line, x)
            tmp = line[x]
            line[x] = line[y]
            line[y] = tmp
        end
    end
end


"""
    get_defined_buses(data_eng::Dict{String,Any}; comp_types=pmd_eng_asset_types)::Vector{String}

Returns a unique list of all buses specified in the data model.
The argument 'comp_types' specifies which component types are searched to build the list.
"""
function get_defined_buses(data_eng::Dict{String,Any}; comp_types=pmd_eng_asset_types)::Vector{String}
    @assert data_eng["data_model"]==ENGINEERING

    buses_exclude = Vector{String}()
    for comp_type in intersect(comp_types, keys(data_eng))
        for (id,comp) in data_eng[comp_type]
            if haskey(comp, "f_bus")
                buses = [comp["f_bus"], comp["t_bus"]]
            elseif haskey(comp, "bus")
                # works for a vector of buses and a single bus as a string
                buses = isa(comp["bus"], Vector) ? comp["bus"] : [comp["bus"]]
            else
                buses = []
            end
            append!(buses_exclude, buses)
        end
    end

    return unique(buses_exclude)
end


"""
    delete_trailing_lines!(data_eng::Dict{String,Any})::Dict{String,Any}

Deletes trailing lines,
i.e. lines connected to a bus with no other connected components and which is not grounded.
"""
function delete_trailing_lines!(data_eng::Dict{String,Any})::Dict{String,Any}
    @assert data_eng["data_model"]==ENGINEERING

    # exclude buses that appear in components other than lines
    comp_types_exclude = setdiff(pmd_eng_asset_types, ["line"])
    buses_exclude = get_defined_buses(data_eng, comp_types=comp_types_exclude)

    # build auxiliary variables
    line_has_shunt = Dict()
    bus_lines = Dict(k=>[] for k in keys(data_eng["bus"]))
    for (id,line) in data_eng["line"]
        _, y_fr, y_to = _get_line_impedance_parameters(data_eng, line)
        line_has_shunt[id] = !iszero(y_fr) || !iszero(y_to)
        push!(bus_lines[line["f_bus"]], id)
        push!(bus_lines[line["t_bus"]], id)
    end

    # eligible buses connect to a single line, and that line should have zero addmittance to ground
    eligible_buses = [bus_id for (bus_id,line_ids) in bus_lines if length(line_ids)==1 && !line_has_shunt[line_ids[1]]]
    # exclude buses that connect to a component other than lines
    eligible_buses = setdiff(eligible_buses, buses_exclude)
    # exclude buses that are grounded
    grounded_buses = [bus_id for (bus_id,bus) in data_eng["bus"] if !isempty(bus["grounded"])]
    eligible_buses = setdiff(eligible_buses, grounded_buses)

    # remove trailing lines and buses
    while !isempty(eligible_buses)
        for bus_id in eligible_buses
            # this trailing bus has one associated line
            line_id = bus_lines[bus_id][1]
            line = data_eng["line"][line_id]

            delete!(data_eng["line"], line_id)
            delete!(data_eng["bus"],  bus_id)

            other_end_bus = line["f_bus"]==bus_id ? line["t_bus"] : line["f_bus"]
            bus_lines[other_end_bus] = setdiff(bus_lines[other_end_bus], [line_id])
            delete!(bus_lines,  bus_id)
        end

        eligible_buses = [bus_id for (bus_id, line_ids) in bus_lines if length(line_ids)==1 && !(bus_id in buses_exclude) && !line_has_shunt[line_ids[1]]]
    end

    return data_eng
end


"Join lines which are connected in series, and of which the intermediate bus is ungrounded and does not connect to any other components."
function join_lines!(data_eng::Dict{String,Any})
    @assert data_eng["data_model"]==ENGINEERING

    # a bus is eligible for reduction if it only appears in exactly two lines
    buses_all = collect(keys(data_eng["bus"]))

    # exclude buses that appear in components other than lines
    comp_types_exclude = setdiff(pmd_eng_asset_types, ["line"])
    buses_exclude = get_defined_buses(data_eng, comp_types=comp_types_exclude)
    # exclude buses that are grounded
    grounded_buses = [bus_id for (bus_id,bus) in data_eng["bus"] if !isempty(bus["grounded"])]
    buses_exclude = union(buses_exclude, grounded_buses)

    # per bus, list all inbound or outbound lines
    bus_lines = Dict(bus=>[] for bus in buses_all)
    for (id,line) in data_eng["line"]
        push!(bus_lines[line["f_bus"]], id)
        push!(bus_lines[line["t_bus"]], id)
    end

    # exclude all buses that do not have exactly two lines connected to it
    buses_exclude = union(buses_exclude, [bus for (bus,lines) in bus_lines if length(lines)!=2])

    # now loop over remaining buses
    candidates = setdiff(buses_all, buses_exclude)
    for bus in candidates
        line1_id, line2_id = bus_lines[bus]
        line1 = data_eng["line"][line1_id]
        line2 = data_eng["line"][line2_id]

        # reverse lines if needed to get the order
        # (x)--fr-line1-to--(bus)--to-line2-fr--(x)
        if line1["f_bus"]==bus
            _line_reverse!(line1)
        end
        if line2["f_bus"]==bus
            _line_reverse!(line2)
        end

        reducable = true
        reducable = reducable && line1["linecode"]==line2["linecode"]
        reducable = reducable && all(line1["t_connections"].==line2["t_connections"])
        if reducable

            line1["length"] += line2["length"]
            line1["t_bus"] = line2["f_bus"]
            line1["t_connections"] = line2["f_connections"]

            delete!(data_eng["line"], line2_id)
            delete!(data_eng["bus"], bus)
            for x in candidates
                if line2_id in bus_lines[x]
                    bus_lines[x] = [setdiff(bus_lines[x], [line2_id])..., line1_id]
                end
            end
        end
    end

    return data_eng
end