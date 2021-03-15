# This file contains useful transformation functions for the engineering data model

const _loss_model_objects = Dict{String,Vector{String}}(
    "switch" => Vector{String}(["linecode", "rs", "xs"]),
    "voltage_source" => Vector{String}(["rs", "xs"]),
    "transformer" => Vector{String}(["rw", "xsc", "cmag", "noloadloss"])
)


"remove parameters from objects with loss models to make them lossless"
function make_lossless!(data_eng::Dict{String,<:Any})
    @assert data_eng["data_model"] == ENGINEERING "incorrect data model type"

    for (object_type, parameters) in _loss_model_objects
        if haskey(data_eng, object_type)
            for (id, eng_obj) in data_eng[object_type]
                for parameter in parameters
                    if haskey(eng_obj, parameter)
                        if parameter == "linecode"
                            delete!(eng_obj, parameter)
                        else
                            eng_obj[parameter] = 0 .* eng_obj[parameter]
                        end
                    end
                end
            end
        end
    end
end


"add voltage bounds"
function apply_voltage_bounds!(data_eng::Dict{String,<:Any}; vm_lb::Union{Real,Missing}=0.9, vm_ub::Union{Real,Missing}=1.1)
    @assert data_eng["data_model"] == ENGINEERING "incorrect data model type"

    (bus_vbases, edge_vbases) = calc_voltage_bases(data_eng, data_eng["settings"]["vbases_default"])
    if haskey(data_eng, "bus")
        for (id,bus) in data_eng["bus"]
            vbase = bus_vbases[id]
            if !ismissing(vm_lb) && !haskey(bus, "vm_lb")
                bus["vm_lb"] = vbase .* fill(vm_lb, length(bus["terminals"]))
            end

            if !ismissing(vm_ub) && !haskey(bus, "vm_ub")
                bus["vm_ub"] = vbase .* fill(vm_ub, length(bus["terminals"]))
            end
        end
    end
end


"removes all fields ending in '_ub' or '_lb'"
function remove_all_bounds!(data_eng; exclude::Vector{<:String}=Vector{String}([]))
    for (k,v) in data_eng
        if isa(v, Dict) && k!="settings" && !(k in exclude)
            for (id, comp) in v
                for field in keys(comp)
                    if endswith(field, "_lb") || endswith(field, "_ub")
                        delete!(comp, field)
                    end
                end
            end
        end
    end
end


"kron reduction"
function apply_kron_reduction!(data_eng::Dict{String,<:Any}; kr_phases::Union{Vector{Int},Vector{String}}=[1,2,3], kr_neutral::Union{Int,String}=4)
    if haskey(data_eng, "bus")
        for (id,eng_obj) in data_eng["bus"]
            filter = eng_obj["terminals"] .!= kr_neutral
            terminals_kr = eng_obj["terminals"][filter]

            @assert all(t in kr_phases for t in terminals_kr) "bus $id has terminals $(eng_obj["terminals"]), outside of $kr_phases, cannot be kron reduced"

            _apply_filter!(eng_obj, ["vm", "va", "vm_lb", "vm_ub"], filter)
            eng_obj["terminals"] = terminals_kr

            gr_filter = eng_obj["grounded"] .!= kr_neutral
            _apply_filter!(eng_obj, ["grounded", "rg", "xg"], gr_filter)
        end
    end

    if haskey(data_eng, "line")
        for (_,eng_obj) in data_eng["line"]
            @assert all(eng_obj["f_connections"].==eng_obj["t_connections"]) "Kron reduction is only supported if f_connections == t_connections"

            _apply_linecode!(eng_obj, data_eng)

            filter = _kron_reduce_branch!(eng_obj, ["rs", "xs"], ["g_fr", "b_fr", "g_to", "b_to"], eng_obj["f_connections"], kr_neutral)
            _apply_filter!(eng_obj, ["vad_lb", "vad_ub", "cm_ub", "cm_ub_b", "cm_ub_c", "sm_ub", "sm_ub_b", "sm_ub_c", "f_connections", "t_connections"], filter)
        end
    end

    if haskey(data_eng, "transformer")
        for (id, eng_obj) in data_eng["transformer"]
            _apply_xfmrcode!(eng_obj, data_eng)

            if haskey(eng_obj, "f_connections")
                @assert all(eng_obj["f_connections"].==eng_obj["t_connections"]) "Kron reduction is only supported if f_connections == t_connections"

                if eng_obj["configuration"] == WYE
                    filter = eng_obj["f_connections"] .!= kr_neutral
                    _apply_filter!(eng_obj, ["f_connections", "t_connections"], filter)
                end
            else
                for (w, connections) in enumerate(eng_obj["connections"])
                    if eng_obj["configuration"][w] == WYE
                        filter = connections .!= kr_neutral
                        _apply_filter!(eng_obj, ["connections"], w, filter)
                    end
                end
            end
        end
    end

    if haskey(data_eng, "switch")
        for (_,eng_obj) in data_eng["switch"]
            @assert all(eng_obj["f_connections"].==eng_obj["t_connections"]) "Kron reduction is only supported if f_connections == t_connections"

            _apply_linecode!(eng_obj, data_eng)

            filter = _kron_reduce_branch!(eng_obj, Vector{String}([k for k in ["rs", "xs"] if haskey(eng_obj, k)]), Vector{String}([]), eng_obj["f_connections"], kr_neutral)
            _apply_filter!(eng_obj, ["vad_lb", "vad_ub", "cm_ub", "cm_ub_b", "cm_ub_c", "sm_ub", "sm_ub_b", "sm_ub_c", "f_connections", "t_connections"], filter)
        end
    end

    if haskey(data_eng, "shunt")
        for (_,eng_obj) in data_eng["shunt"]
            filter = _kron_reduce_branch!(eng_obj, Vector{String}([]), ["gs", "bs"], eng_obj["connections"], kr_neutral)
            _apply_filter!(eng_obj, ["connections"], filter)
        end
    end

    if haskey(data_eng, "load")
        for (_,eng_obj) in data_eng["load"]
            if eng_obj["configuration"]==WYE
                @assert eng_obj["connections"][end] == kr_neutral "for wye-connected loads, to kron reduce the connections list should end with a neutral"

                filter = eng_obj["connections"] .!= kr_neutral
                _apply_filter!(eng_obj, ["connections"], filter)
            end
        end
    end

    if haskey(data_eng, "generator")
        for (_,eng_obj) in data_eng["generator"]
            if eng_obj["configuration"]==WYE
                @assert eng_obj["connections"][end] == kr_neutral "for wye-connected generators, to kron reduce the connections list should end with a neutral"

                filter = eng_obj["connections"] .!= kr_neutral
                _apply_filter!(eng_obj, ["connections"], filter)
            end
        end
    end

    if haskey(data_eng, "solar")
        for (_,eng_obj) in data_eng["solar"]
            if eng_obj["configuration"]==WYE
                @assert eng_obj["connections"][end] == kr_neutral "for wye-connected solar, to kron reduce the connections list should end with a neutral"

                filter = eng_obj["connections"] .!= kr_neutral
                _apply_filter!(eng_obj, ["connections"], filter)
            end
        end
    end

    if haskey(data_eng, "storage")
        for (_,eng_obj) in data_eng["storage"]
            filter = eng_obj["connections"] .!= kr_neutral
            _apply_filter!(eng_obj, ["connections"], filter)
        end
    end

    if haskey(data_eng, "voltage_source")
        for (_,eng_obj) in data_eng["voltage_source"]
            filter = eng_obj["connections"] .!= kr_neutral
            _apply_filter!(eng_obj, ["vm", "va", "vm_lb", "vm_ub", "rs", "xs", "connections"], filter)
        end
    end

    delete!(data_eng, "linecode")  # kron reduction moves linecode properties directly to lines
    delete!(data_eng, "xfmrcode")  # kron reduction moves xfmrcode properties directly to transformers

    data_eng["is_kron_reduced"] = true
end


"pad matrices to max number of conductors"
function apply_phase_projection!(data_eng::Dict{String,<:Any})
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

            _pad_properties!(eng_obj, ["vm", "va", "vm_lb", "vm_ub"], eng_obj["terminals"], all_conductors)

            _pad_connections!(eng_obj, "terminals", all_conductors)
        end
    end

    if haskey(data_eng, "line")
        for (_,eng_obj) in data_eng["line"]
            _apply_linecode!(eng_obj, data_eng)

            _pad_properties!(eng_obj, ["rs", "xs", "g_fr", "g_to", "b_fr", "b_to"], eng_obj["f_connections"], all_conductors)
            _pad_properties!(eng_obj, ["vad_lb"], eng_obj["f_connections"], all_conductors; pad_value=-60.0)
            _pad_properties!(eng_obj, ["vad_ub"], eng_obj["f_connections"], all_conductors; pad_value=60.0)
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
            _pad_properties!(eng_obj, ["vad_lb"], eng_obj["f_connections"], all_conductors; pad_value=-60.0)
            _pad_properties!(eng_obj, ["vad_ub"], eng_obj["f_connections"], all_conductors; pad_value=60.0)
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


"phase projection for components where unprojected states are not yet supported"
function apply_phase_projection_delta!(data_eng::Dict{String,<:Any})
    @assert get(data_eng, "is_kron_reduced", false)

    bus_terminals = Dict{Any,Vector{Int}}()

    all_conductors = _get_complete_conductor_set(data_eng)

    if haskey(data_eng, "transformer")
        for (_,eng_obj) in data_eng["transformer"]
            _apply_xfmrcode!(eng_obj, data_eng)

            if haskey(eng_obj, "f_connections")
                _pad_properties!(eng_obj, ["tm_lb", "tm_ub", "tm_set"], eng_obj["f_connections"], all_conductors; pad_value=1.0)
                _pad_properties!(eng_obj, ["tm_fix"], eng_obj["f_connections"], all_conductors; pad_value=true)

                _pad_connections!(eng_obj, "f_connections", all_conductors)
                _pad_connections!(eng_obj, "t_connections", all_conductors)

                bus_terminals[eng_obj["f_bus"]] = haskey(bus_terminals, eng_obj["f_bus"]) ? unique([bus_terminals[eng_obj["f_bus"]]..., all_conductors...]) : all_conductors
                bus_terminals[eng_obj["t_bus"]] = haskey(bus_terminals, eng_obj["t_bus"]) ? unique([bus_terminals[eng_obj["t_bus"]]..., all_conductors...]) : all_conductors
            else
                for (w, connections) in enumerate(eng_obj["connections"])
                    if eng_obj["configuration"][w] == DELTA && length(eng_obj["tm_set"][w]) == 1
                        _pad_properties_delta!(eng_obj, ["tm_lb", "tm_ub", "tm_set"], connections, w, all_conductors)
                        _pad_properties_delta!(eng_obj, ["tm_fix"], connections, w, all_conductors)
                        for i in 1:length(eng_obj["tm_set"][w])
                            eng_obj["tm_set"][w][i] = eng_obj["tm_set"][w][i] == 0 ? 1 : eng_obj["tm_set"][w][i]
                            eng_obj["tm_fix"][w][i] = eng_obj["tm_fix"][w][i] == 0 ? true : eng_obj["tm_fix"][w][i]
                        end
                    else
                        _pad_properties!(eng_obj, ["tm_lb", "tm_ub", "tm_set"], w, connections, all_conductors; pad_value=1.0)
                        _pad_properties!(eng_obj, ["tm_fix"], w, connections, all_conductors; pad_value=true)
                    end
                    _pad_connections!(eng_obj, "connections", w, all_conductors)
                    bus_terminals[eng_obj["bus"][w]] = haskey(bus_terminals, eng_obj["bus"][w]) ? unique([bus_terminals[eng_obj["bus"][w]]..., all_conductors...]) : all_conductors
                end
            end
        end
    end

    if haskey(data_eng, "load")
        for (_,eng_obj) in data_eng["load"]
            if eng_obj["configuration"] == DELTA
                _pad_properties_delta!(eng_obj, ["pd_nom", "qd_nom"], eng_obj["connections"], all_conductors)
                _pad_connections!(eng_obj, "connections", all_conductors)
                bus_terminals[eng_obj["bus"]] = haskey(bus_terminals, eng_obj["bus"]) ? unique([bus_terminals[eng_obj["bus"]]..., all_conductors...]) : all_conductors
            end
        end
    end

    if haskey(data_eng, "generator")
        for (_,eng_obj) in data_eng["generator"]
            if eng_obj["configuration"]==DELTA
                _pad_properties_delta!(eng_obj, ["pg", "qg", "vg", "pg_lb", "pg_ub", "qg_lb", "qg_ub"], eng_obj["connections"], all_conductors)
                _pad_connections!(eng_obj, "connections", all_conductors)
                bus_terminals[eng_obj["bus"]] = haskey(bus_terminals, eng_obj["bus"]) ? unique([bus_terminals[eng_obj["bus"]]..., all_conductors...]) : all_conductors
            end
        end
    end

    if haskey(data_eng, "solar")
        for (_,eng_obj) in data_eng["solar"]
            if eng_obj["configuration"]==DELTA
                _pad_properties_delta!(eng_obj, ["pg", "qg", "vg", "pg_lb", "pg_ub", "qg_lb", "qg_ub"], eng_obj["connections"], all_conductors)
                _pad_connections!(eng_obj, "connections", all_conductors)
                bus_terminals[eng_obj["bus"]] = haskey(bus_terminals, eng_obj["bus"]) ? unique([bus_terminals[eng_obj["bus"]]..., all_conductors...]) : all_conductors
            end
        end
    end

    _update_bus_terminal_projections!(data_eng, bus_terminals)

    data_eng["is_projected"] = true
end


function _update_bus_terminal_projections!(data_eng::Dict{String,<:Any}, bus_terminals::Dict{Any,<:Vector{Int}})
    for (id,terminals) in bus_terminals
        eng_obj = data_eng["bus"][id]

        if !haskey(eng_obj, "vm_lb")
            eng_obj["vm_lb"] = fill(0.0, length(eng_obj["terminals"]))
        end

        if !haskey(eng_obj, "vm_ub")
            eng_obj["vm_ub"] = fill(Inf, length(eng_obj["terminals"]))
        end

        _pad_properties!(eng_obj, ["vm", "va", "vm_lb", "vm_ub"], eng_obj["terminals"], terminals)

        _pad_connections!(eng_obj, "terminals", terminals)
    end
end


"voltage angle bounds"
function apply_voltage_angle_bounds!(eng::Dict{String,<:Any}, vad::Real)
    if haskey(eng, "line")
        for (_,line) in eng["line"]
            line["vad_lb"] = fill(-vad, length(line["f_connections"]))
            line["vad_ub"] = fill( vad, length(line["f_connections"]))
        end
    end
end
