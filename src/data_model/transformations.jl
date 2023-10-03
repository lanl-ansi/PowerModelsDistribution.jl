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
    "storage_ne" => String["rs", "xs", "pex", "qex"],
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
    adjust_small_line_impedances!(data::Dict{String,<:Any}; min_impedance_val::Float64=1e-2, replace_impedance_val::Float64=0.0)

Replaces impedances (rs, xs) on lines, linecodes, and switches lower than `min_impedance_val` with `replace_impedance_val`.
"""
function adjust_small_line_impedances!(data::Dict{String,<:Any}; min_impedance_val::Float64=1e-2, replace_impedance_val::Float64=0.0)
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_adjust_small_line_impedances!, data; apply_to_subnetworks=true, min_impedance_val=min_impedance_val, replace_impedance_val=replace_impedance_val)
end


"""
    _adjust_small_line_impedances!(data_eng::Dict{String,<:Any}; min_impedance_val::Float64=1e-2, replace_impedance_val::Float64=0.0)

Replaces impedances (rs, xs) on lines, linecodes, and switches lower than `min_impedance_val` with `replace_impedance_val`.
"""
function _adjust_small_line_impedances!(data_eng::Dict{String,<:Any}; min_impedance_val::Float64=1e-2, replace_impedance_val::Float64=0.0)
    for type in ["line", "linecode", "switch"]
        if haskey(data_eng, type)
            for (id,obj) in data_eng[type]
                for k in ["rs", "xs"]
                    if haskey(obj, k)
                        obj[k][obj[k].<min_impedance_val] .= replace_impedance_val
                    end
                end
            end
        end
    end
end


"""
    adjust_small_line_admittances!(data::Dict{String,<:Any}; min_admittance_val::Float64=1e-2, replace_admittance_val::Float64=0.0)

Replaces admittances (g_fr, g_to, b_fr, b_to) on lines, linecodes, and switches lower than `min_admittance_val` with `replace_admittance_val`.
"""
function adjust_small_line_admittances!(data::Dict{String,<:Any}; min_admittance_val::Float64=1e-2, replace_admittance_val::Float64=0.0)
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_adjust_small_line_admittances!, data; apply_to_subnetworks=true, min_admittance_val=min_admittance_val, replace_admittance_val=replace_admittance_val)
end


"""
    _adjust_small_line_admittances!(data_eng::Dict{String,<:Any}; min_admittance_val::Float64=1e-2, replace_admittance_val::Float64=0.0)

Replaces admittances (g_fr, g_to, b_fr, b_to) on lines, linecodes, and switches lower than `min_admittance_val` with `replace_admittance_val`.
"""
function _adjust_small_line_admittances!(data_eng::Dict{String,<:Any}; min_admittance_val::Float64=1e-2, replace_admittance_val::Float64=0.0)
    for type in ["line", "linecode", "switch"]
        if haskey(data_eng, type)
            for (id,obj) in data_eng[type]
                for k in ["g_fr", "b_fr", "g_to", "b_to"]
                    if haskey(obj, k)
                        obj[k][obj[k].<min_admittance_val] .= replace_admittance_val
                    end
                end
            end
        end
    end
end


"""
    adjust_small_line_lengths!(data::Dict{String,<:Any}; min_length_val::Float64=25.0, replace_length_val::Float64=0.0)

Replaces length on lines, switches lower than `min_length_val` with `replace_length_val`.
"""
function adjust_small_line_lengths!(data::Dict{String,<:Any}; min_length_val::Float64=25.0, replace_length_val::Float64=0.0)
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_adjust_small_line_lengths!, data; apply_to_subnetworks=true, min_length_val=min_length_val, replace_length_val=replace_length_val)
end


"""
    _adjust_small_line_lengths!(data_eng::Dict{String,<:Any}; min_length_val::Float64=25.0, replace_length_val::Float64=0.0)

Replaces length on lines, switches lower than `min_length_val` with `replace_length_val`.
"""
function _adjust_small_line_lengths!(data_eng::Dict{String,<:Any}; min_length_val::Float64=25.0, replace_length_val::Float64=0.0)
    for type in ["line", "switch"]
        if haskey(data_eng, type)
            for (id,obj) in data_eng[type]
                if haskey(obj, "length") &&  obj["length"] < min_length_val
                    obj["length"] = replace_length_val
                end
            end
        end
    end
end


"""
    apply_voltage_bounds!(data::Dict{String,<:Any}; vm_lb::Union{Real,Missing}=0.9, vm_ub::Union{Real,Missing}=1.1)

add voltage bounds to all buses based on per-unit upper (`vm_ub`) and lower (`vm_lb`), scaled by the bus's voltage based
"""
function apply_voltage_bounds!(data::Dict{String,<:Any}; vm_lb::Union{Real,Missing}=0.9, vm_ub::Union{Real,Missing}=1.1, exclude::Vector{String}=!isempty(get(data, "voltage_source", Dict())) ? String[x.second["bus"] for x in data["voltage_source"]] : String[])
    # @assert iseng(data) "wrong data model type"
    if iseng(data)
        apply_pmd!(_apply_voltage_bounds!, data; apply_to_subnetworks=true, vm_lb=vm_lb, vm_ub=vm_ub)
    else
        apply_pmd!(_apply_voltage_bounds_math!, data; apply_to_subnetworks=true, vm_lb=vm_lb, vm_ub=vm_ub)
    end
end


"""
    _apply_voltage_bounds!(data_eng::Dict{String,<:Any}; vm_lb::Union{Real,Missing}=0.9, vm_ub::Union{Real,Missing}=1.1)

add voltage bounds to all buses based on per-unit upper (`vm_ub`) and lower (`vm_lb`), scaled by the bus's voltage based
"""
function _apply_voltage_bounds!(data_eng::Dict{String,<:Any}; vm_lb::Union{Real,Missing}=0.9, vm_ub::Union{Real,Missing}=1.1, exclude::Vector{String}=!isempty(get(data_eng, "voltage_source", Dict())) ? String[x.second["bus"] for x in data_eng["voltage_source"]] : String[])
    (bus_vbases, edge_vbases) = calc_voltage_bases(data_eng, data_eng["settings"]["vbases_default"])
    for (id, bus) in filter(x->!(x.first in exclude), get(data_eng, "bus", Dict{String,Any}()))
        vbase = bus_vbases[id]
        if !ismissing(vm_lb) && !ismissing(vbase)
            data_eng["bus"][id]["vm_lb"] = vbase .* fill(vm_lb, length(bus["terminals"]))
            data_eng["bus"][id]["vm_lb"][any.(bus["grounded"] .== t for t in bus["terminals"])] .= 0.0
        end

        if !ismissing(vm_ub) && !ismissing(vbase)
            data_eng["bus"][id]["vm_ub"] = vbase .* fill(vm_ub, length(bus["terminals"]))
            data_eng["bus"][id]["vm_ub"][any.(bus["grounded"] .== t for t in bus["terminals"])] .= Inf
        end
    end
end


"""
    _apply_voltage_bounds!(data_math::Dict{String,<:Any}; vm_lb::Union{Real,Missing}=0.9, vm_ub::Union{Real,Missing}=1.1)

add voltage bounds to all buses based on per-unit upper (`vm_ub`) and lower (`vm_lb`)
"""
function _apply_voltage_bounds_math!(data_math::Dict{String,<:Any}; vm_lb::Union{Real,Missing}=0.9, vm_ub::Union{Real,Missing}=1.1, exclude::Vector{String}=!isempty(get(data_eng, "voltage_source", Dict())) ? String[x.second["bus"] for x in data_eng["voltage_source"]] : String[])
    # (bus_vbases, edge_vbases) = calc_voltage_bases(data_math, data_eng["settings"]["vbases_default"])
    data = data_math["nw"]
    for network in data
        for (id, bus) in filter(x->!(x.first in exclude), get(network, "bus", Dict{String,Any}()))
            # vbase = bus_vbases[id]
            if !ismissing(vm_lb)
                network["bus"][id]["vm_lb"] = fill(vm_lb, length(bus["terminals"]))
                network["bus"][id]["vm_lb"][any.(bus["grounded"] .== t for t in bus["terminals"])] .= 0.0
            end

            if !ismissing(vm_ub)
                network["bus"][id]["vm_ub"] = fill(vm_ub, length(bus["terminals"]))
                network["bus"][id]["vm_ub"][any.(bus["grounded"] .== t for t in bus["terminals"])] .= Inf
            end
        end
    end
end

"""
    remove_all_bounds!(data; exclude::Vector{<:String}=String["energy_ub"], exclude_asset_type::Vector{String}=String[])

Removes all fields ending in '_ub' or '_lb' that aren't required by the math model. Properties
can be excluded from this removal with `exclude::Vector{String}`

Whole asset types (e.g., "line") can be excluded using the keyword argument `exclude_asset_type::Vector{String}`

By default, `"energy_ub"` is excluded from this removal, since it is a required properly on storage.
"""
function remove_all_bounds!(data; exclude::Vector{<:String}=String["energy_ub"], exclude_asset_type::Vector{String}=String[])
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_remove_all_bounds!, data; apply_to_subnetworks=true, exclude=exclude, exclude_asset_type=exclude_asset_type)
end


"""
    _remove_all_bounds!(data_eng; exclude::Vector{<:String}=String["energy_ub"], exclude_asset_type::Vector{String}=String[])

Removes all fields ending in '_ub' or '_lb' that aren't required by the math model. Properties
can be excluded from this removal with `exclude::Vector{String}`

Whole asset types (e.g., "line") can be excluded using the keyword argument `exclude_asset_type::Vector{String}`

By default, `"energy_ub"` is excluded from this removal, since it is a required properly on storage.
"""
function _remove_all_bounds!(data_eng; exclude::Vector{<:String}=String["energy_ub"], exclude_asset_type::Vector{String}=String[])
    for (k,v) in data_eng
        if isa(v, Dict) && k!="settings" && !(k in exclude_asset_type)
            for (_, comp) in v
                for field in keys(comp)
                    if !(field in exclude) && (endswith(field, "_lb") || endswith(field, "_ub"))
                        delete!(comp, field)
                    end
                end
            end
        end
    end
end


"""
    adjust_line_limits!(data::Dict{String,<:Any}, mult::Float64)

Multiplies limits (`sm_ub` and/or `cm_ub`) on line objects (`line`, `linecode`, `switch`) by a multiplier `mult`
"""
function adjust_line_limits!(data::Dict{String,<:Any}, mult::Float64)
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_adjust_line_limits!, data, mult; apply_to_subnetworks=true)
end


"""
    _adjust_line_limits!(data_eng::Dict{String,<:Any}, mult::Float64)

Multiplies limits (`sm_ub` and/or `cm_ub`) on line objects (`line`, `linecode`, `switch`) by a multiplier `mult`
"""
function _adjust_line_limits!(data_eng::Dict{String,<:Any}, mult::Float64)
    for type in ["linecode", "line", "switch"]
        if haskey(data_eng, type)
            for (_,obj) in data_eng[type]
                if haskey(obj, "cm_ub")
                    obj["cm_ub"] .*= mult
                end
                if haskey(obj, "sm_ub")
                    obj["sm_ub"] .*= mult
                end
            end
        end
    end
end


"""
    remove_line_limits!(data::Dict{String,<:Any})

Removes fields `cm_ub` and `sm_ub` from lines, switches, and linecodes
"""
function remove_line_limits!(data::Dict{String,<:Any})
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_remove_line_limits!, data; apply_to_subnetworks=true)
end


"""
    _remove_line_limits!(data_eng::Dict{String,<:Any})

Removes fields `cm_ub` and `sm_ub` from lines, switches, and linecodes
"""
function _remove_line_limits!(data_eng::Dict{String,<:Any})
    for type in ["linecode", "line", "switch"]
        if haskey(data_eng, type)
            for (_,obj) in data_eng[type]
                delete!(obj, "cm_ub")
                delete!(obj, "sm_ub")
            end
        end
    end
end


"""
    adjust_transformer_limits!(data::Dict{String,<:Any}, mult::Float64)

Multiplies limits (`sm_ub` and/or `cm_ub`) on transformer objects by a multiplier `mult`
"""
function adjust_transformer_limits!(data::Dict{String,<:Any}, mult::Float64)
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_adjust_transformer_limits!, data, mult; apply_to_subnetworks=true)
end


"""
    _adjust_transformer_limits!(data_eng::Dict{String,<:Any}, mult::Float64)

Multiplies limits (`sm_ub` and/or `cm_ub`) on transformer objects by a multiplier `mult`
"""
function _adjust_transformer_limits!(data_eng::Dict{String,<:Any}, mult::Float64)
    if haskey(data_eng, "transformer")
        for (_,obj) in data_eng["transformer"]
            if haskey(obj, "cm_ub")
                obj["cm_ub"] *= mult
            end
            if haskey(obj, "sm_ub")
                obj["sm_ub"] *= mult
            end
        end
    end
end


"""
    remove_transformer_limits!(data_eng::Dict{String,<:Any})

Removes field `sm_ub` from transformers, xfmrcodes
"""
function remove_transformer_limits!(data::Dict{String,<:Any})
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_remove_transformer_limits!, data; apply_to_subnetworks=true)
end


"""
    _remove_transformer_limits!(data_eng::Dict{String,<:Any})

Removes field `sm_ub` from transformers, xfmrcodes
"""
function _remove_transformer_limits!(data_eng::Dict{String,<:Any})
    for type in ["xfmrcode", "transformer"]
        if haskey(data_eng, type)
            for (_,obj) in data_eng[type]
                delete!(obj, "sm_ub")
            end
        end
    end
end


"""
    apply_kron_reduction!(data::Dict{String,<:Any}; kr_phases::Union{Vector{Int},Vector{String}}=[1,2,3], kr_neutral::Union{Int,String}=4)

Applies a Kron Reduction to the network, reducing out the `kr_neutral`, leaving only the `kr_phases`
"""
function apply_kron_reduction!(data::Dict{String,<:Any}; kr_phases::Union{Vector{Int},Vector{String}}=[1,2,3], kr_neutral::Union{Int,String}=4)
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_apply_kron_reduction!, data; apply_to_subnetworks=true, kr_phases=kr_phases, kr_neutral=kr_neutral)
end


"""
    _apply_kron_reduction!(data_eng::Dict{String,<:Any}; kr_phases::Union{Vector{Int},Vector{String}}=[1,2,3], kr_neutral::Union{Int,String}=4)

Applies a Kron Reduction to the network, reducing out the `kr_neutral`, leaving only the `kr_phases`
"""
function _apply_kron_reduction!(data_eng::Dict{String,<:Any}; kr_phases::Union{Vector{Int},Vector{String}}=[1,2,3], kr_neutral::Union{Int,String}=4)
    if !get(data_eng, "is_kron_reduced", false)
        if haskey(data_eng, "bus")
            for (id,eng_obj) in data_eng["bus"]
                filter = eng_obj["terminals"] .!= kr_neutral
                terminals_kr = eng_obj["terminals"][filter]

                @assert all(t in kr_phases for t in terminals_kr) "bus $id has terminals $(eng_obj["terminals"]), outside of $kr_phases, cannot be kron reduced"

                _apply_filter!(eng_obj, ["vm", "va", "vm_lb", "vm_ub", "vm_start", "va_start"], filter)
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
                _apply_filter!(eng_obj, ["vad_lb", "vad_ub", "cm_ub", "sm_ub", "f_connections", "t_connections"], filter)
            end
        end

        if haskey(data_eng, "transformer")
            for (id, eng_obj) in data_eng["transformer"]
                _apply_xfmrcode!(eng_obj, data_eng)

                if haskey(eng_obj, "f_connections")
                    @assert all(eng_obj["f_connections"].==eng_obj["t_connections"]) "Kron reduction is only supported if f_connections == t_connections"

                    filter = eng_obj["f_connections"] .!= kr_neutral
                    _apply_filter!(eng_obj, ["f_connections", "t_connections"], filter)
                else
                    for (w, connections) in enumerate(eng_obj["connections"])
                        filter = connections .!= kr_neutral
                        _apply_filter!(eng_obj, ["connections"], w, filter)
                    end
                end
            end
        end

        if haskey(data_eng, "switch")
            for (_,eng_obj) in data_eng["switch"]
                @assert all(eng_obj["f_connections"].==eng_obj["t_connections"]) "Kron reduction is only supported if f_connections == t_connections"

                _apply_linecode!(eng_obj, data_eng)

                filter = _kron_reduce_branch!(eng_obj, Vector{String}([k for k in ["rs", "xs"] if haskey(eng_obj, k)]), String[], eng_obj["f_connections"], kr_neutral)
                _apply_filter!(eng_obj, ["vad_lb", "vad_ub", "cm_ub", "sm_ub", "f_connections", "t_connections"], filter)
            end
        end

        if haskey(data_eng, "shunt")
            for (_,eng_obj) in data_eng["shunt"]
                filter = _kron_reduce_branch!(eng_obj, String[], ["gs", "bs"], eng_obj["connections"], kr_neutral)
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

        if haskey(data_eng, "solar_ne")
            for (_,eng_obj) in data_eng["solar_ne"]
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

        if haskey(data_eng, "storage_ne")
            for (_,eng_obj) in data_eng["storage_ne"]
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

        find_conductor_ids!(data_eng)
        data_eng["is_kron_reduced"] = true
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

    if haskey(data_eng, "solar_ne")
        for (_,eng_obj) in data_eng["solar_ne"]
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
                _pad_properties_delta!(eng_obj, ["pd_nom", "qd_nom"], eng_obj["connections"], all_conductors)
                _pad_connections!(eng_obj, "connections", all_conductors)
                bus_terminals[eng_obj["bus"]] = haskey(bus_terminals, eng_obj["bus"]) ? _pad_connections!(bus_terminals, eng_obj["bus"], eng_obj["connections"]) : eng_obj["connections"]
            end
        end
    end

    if haskey(data_eng, "generator")
        for (_,eng_obj) in data_eng["generator"]
            if eng_obj["configuration"]==DELTA
                _pad_properties_delta!(eng_obj, ["pg", "qg", "vg", "pg_lb", "pg_ub", "qg_lb", "qg_ub"], eng_obj["connections"], all_conductors)
                _pad_connections!(eng_obj, "connections", all_conductors)
                bus_terminals[eng_obj["bus"]] = haskey(bus_terminals, eng_obj["bus"]) ? _pad_connections!(bus_terminals, eng_obj["bus"], eng_obj["connections"]) : eng_obj["connections"]
            end
        end
    end

    if haskey(data_eng, "solar")
        for (_,eng_obj) in data_eng["solar"]
            if eng_obj["configuration"]==DELTA
                _pad_properties_delta!(eng_obj, ["pg", "qg", "vg", "pg_lb", "pg_ub", "qg_lb", "qg_ub"], eng_obj["connections"], all_conductors)
                _pad_connections!(eng_obj, "connections", all_conductors)
                bus_terminals[eng_obj["bus"]] = haskey(bus_terminals, eng_obj["bus"]) ? _pad_connections!(bus_terminals, eng_obj["bus"], eng_obj["connections"]) : eng_obj["connections"]
            end
        end
    end

    if haskey(data_eng, "solar_ne")
        for (_,eng_obj) in data_eng["solar_ne"]
            if eng_obj["configuration"]==DELTA
                _pad_properties_delta!(eng_obj, ["pg", "qg", "vg", "pg_lb", "pg_ub", "qg_lb", "qg_ub"], eng_obj["connections"], all_conductors)
                _pad_connections!(eng_obj, "connections", all_conductors)
                bus_terminals[eng_obj["bus"]] = haskey(bus_terminals, eng_obj["bus"]) ? _pad_connections!(bus_terminals, eng_obj["bus"], eng_obj["connections"]) : eng_obj["connections"]
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


"""
    apply_voltage_angle_difference_bounds!(data::Dict{String,<:Any}, vad::Real=5.0)

Applies voltage angle difference bound given by `vad::Real` in degrees (_i.e._, the allowed drift of angle from one end
of a line to another) to all lines. By default, `vad=5.0`.
"""
function apply_voltage_angle_difference_bounds!(data::Dict{String,<:Any}, vad::Real=5.0)
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_apply_voltage_angle_difference_bounds!, data, vad; apply_to_subnetworks=true)
end


"""
    _apply_voltage_angle_difference_bounds!(eng::Dict{String,<:Any}, vad::Real=5.0)

Applies voltage angle difference bound given by `vad::Real` in degrees (_i.e._, the allowed drift of angle from one end
of a line to another) to all lines. By default, `vad=5.0`.
"""
function _apply_voltage_angle_difference_bounds!(eng::Dict{String,<:Any}, vad::Real=5.0)
    if haskey(eng, "line")
        for (_,line) in eng["line"]
            line["vad_lb"] = fill(-vad, length(line["f_connections"]))
            line["vad_ub"] = fill( vad, length(line["f_connections"]))
        end
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


"Return the Kron-reduction of the specified neutral conductors of a series impedance matrix."
function _kron_reduce_series_impedance(Z::Matrix, neutral_conductors::Vector{Int})
    # neutral conductor idxs
    N = neutral_conductors
    # phase conductor idxs (complement)
    P = setdiff(1:size(Z)[1], N)

    Zkr = Z[P,P] - Z[P,N]*inv(Z[N,N])*Z[N,P]

    return Zkr
end


"Return the Kron-reduction of the specified neutral conductors of a shunt addmittance matrix."
function _kron_reduce_shunt_addmittance(Y::Matrix, neutral_conductors::Vector{Int})
    # phase conductor idxs (complement)
    P = setdiff(1:size(Y)[1], neutral_conductors)

    Ykr = Y[P,P]

    return Ykr
end


"Kron-reduce specified neutral conductors of a linecode."
function _kron_reduce_linecode!(l, neutral_conductors::Vector{Int})
    z_s = _kron_reduce_series_impedance(l["rs"].+im*l["xs"], neutral_conductors)
    y_fr = _kron_reduce_shunt_addmittance(l["g_fr"].+im*l["b_fr"], neutral_conductors)
    y_to = _kron_reduce_shunt_addmittance(l["g_to"].+im*l["b_to"], neutral_conductors)
    l["rs"] = real.(z_s)
    l["xs"] = imag.(z_s)
    l["g_fr"] = real.(y_fr)
    l["b_fr"] = imag.(y_fr)
    l["g_to"] = real.(y_to)
    l["b_to"] = imag.(y_to)
end


"""
    kron_reduce_implicit_neutrals!(data::Dict{String,Any})::Dict{String,Any}

Kron-reduce all (implied) neutral conductors of lines, switches and shunts, and remove any terminals which become unconnected.
A line or switch conductor is considered as a neutral conductor if it is connected between two neutral terminals.
A terminal is a neutral terminals if it is galvanically connected (i.e. through a line or switch)
to a grounded terminal, or the neutral conductor of a wye-connected component.
"""
function kron_reduce_implicit_neutrals!(data::Dict{String,Any})::Dict{String,Any}
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_kron_reduce_implicit_neutrals!, data; apply_to_subnetworks=true)
end


"""
    _kron_reduce_implicit_neutrals!(data_eng::Dict{String,Any})::Dict{String,Any}

Kron-reduce all (implied) neutral conductors of lines, switches and shunts, and remove any terminals which become unconnected.
A line or switch conductor is considered as a neutral conductor if it is connected between two neutral terminals.
A terminal is a neutral terminals if it is galvanically connected (i.e. through a line or switch)
to a grounded terminal, or the neutral conductor of a wye-connected component.
"""
function _kron_reduce_implicit_neutrals!(data_eng::Dict{String,Any})::Dict{String,Any}
    # store the original linecode ids to detect naming clashes
    orig_lc_ids = [keys(data_eng["linecode"])...]

    # obtain all (implicit) neutral terminals
    nbts = _infer_neutral_terminals(data_eng)

    # Kron-reduce each line if eligible
    for (id,line) in get(data_eng, "line", Dict())
        doubly_grounded = [(line["f_bus"],t_fr) in nbts && (line["t_bus"],t_to) in nbts for (t_fr, t_to) in zip(line["f_connections"], line["t_connections"])]
        if any(doubly_grounded)
            keep = (!).(doubly_grounded)
            neutral_conductors = findall(doubly_grounded)
            _apply_filter!(line, ["f_connections", "t_connections", "cm_ub", "cm_ub_b", "cm_ub_c"], keep)
            if haskey(line, "linecode")
                suffix = "_kr_"*join(findall(doubly_grounded), ".")
                lc_orig_id = line["linecode"]
                lc_kr_id = lc_orig_id*suffix
                @assert !(lc_kr_id in orig_lc_ids) "Kron-reduced linecode naming clashes with original linecode names."
                line["linecode"] = lc_kr_id
                if !haskey(data_eng["linecode"], lc_kr_id)
                    data_eng["linecode"][lc_kr_id] = lc = deepcopy(data_eng["linecode"][lc_orig_id])
                    _kron_reduce_linecode!(lc, neutral_conductors)
                end
            end
            if haskey(line, "rs")
                _kron_reduce_linecode!(line, neutral_conductors)
            end
        end

        return data_eng
    end

    # Kron-reduce each shunt if eligible
    for (id,shunt) in get(data_eng, "shunt", Dict())
        cond_is_neutral = [(shunt["bus"],t) in nbts for t in shunt["connections"]]
        cond_keep = (!).(cond_is_neutral)
        _apply_filter!(shunt, ["connections"], cond_keep)
        Ys = _kron_reduce_shunt_addmittance(shunt["gs"].+im*shunt["bs"], findall(cond_is_neutral))
        shunt["gs"] = real.(Ys)
        shunt["bs"] = imag.(Ys)
    end

    # Kron-reduce each switch if eligible
    for (id,switch) in get(data_eng, "switch", Dict())
        doubly_grounded = [(switch["f_bus"],t_fr) in nbts && (switch["t_bus"],t_to) in nbts for (t_fr, t_to) in zip(switch["f_connections"], switch["t_connections"])]
        if any(doubly_grounded)
            keep = (!).(doubly_grounded)
            neutral_conductors = findall(doubly_grounded)
            _apply_filter!(switch, ["f_connections", "t_connections", "cm_ub", "cm_ub_b", "cm_ub_c"], keep)
            if haskey(switch, "linecode")
                suffix = "_kr_"*join(findall(doubly_grounded), ".")
                lc_orig_id = switch["linecode"]
                lc_kr_id = lc_orig_id*suffix
                @assert !(lc_kr_id in orig_lc_ids) "Kron-reduced linecode naming clashes with original linecode names."
                line["linecode"] = lc_kr_id
                if !haskey(data_eng["linecode"], lc_kr_id)
                    data_eng["linecode"][lc_kr_id] = _kron_reduce_linecode(data_eng["linecode"][lc_orig_id], neutral_conductors)
                end
            end
            if haskey(switch, "rs")
                Zs = _kron_reduce_series_impedance(switch["rs"].+im*switch["xs"], findall(doubly_grounded))
                switch["rs"] = real.(Zs)
                switch["xs"] = imag.(Zs)
            end
        end
    end

    # remove unconnected terminals (likely to be caused by earlier Kron-reductions)
    remove_unconnected_terminals!(data_eng)

    # ground remaining neutral terminals
    remaining_bts = [(b, t) for (b,bus) in data_eng["bus"] for t in bus["terminals"]]
    for (b,t) in intersect(nbts, remaining_bts)
        bus = data_eng["bus"][b]
        perfectly_grounded = bus["grounded"][iszero.(bus["rg"].+im*bus["xg"])]
        if !(t in perfectly_grounded)
            push!(bus["grounded"], t)
            push!(bus["rg"], 0.0)
            push!(bus["xg"], 0.0)
        end
    end

    # update 'conductor_ids' property
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
Return a list of all implicit neutrals as a list of bus-terminal pairs.
This is done by starting from a list of all terminals which are either
    a.connected to the neutral of wye-connected components;
    b. or are grounded.
This initial list is then expanded to all terminals which are
galvanically connected to terminals in the initial list.
"""
function _infer_neutral_terminals(data_eng::Dict{String,Any})
    # find all galvanic edges through lines and switches
    bts = [(b,t) for (b,bus) in data_eng["bus"] for t in bus["terminals"]]
    bt_neighbours = Dict(bt=>[] for bt in bts)
    for type in intersect(["line", "switch"], keys(data_eng))
        for (_,comp) in data_eng[type]
            f_bus = comp["f_bus"]
            t_bus = comp["t_bus"]
            for (t_fr, t_to) in zip(comp["f_connections"], comp["t_connections"])
                push!(bt_neighbours[(f_bus, t_fr)], (t_bus, t_to))
                push!(bt_neighbours[(t_bus, t_to)], (f_bus, t_fr))
            end
        end
    end

    # find starpoints of wye-connected devices
    starpoints = []
    for type in setdiff(intersect(pmd_eng_asset_types, keys(data_eng)), ["bus", "line", "switch"])
        for (_,comp) in data_eng[type]
            if haskey(comp, "configuration")
                conf = comp["configuration"]
                if isa(conf, Vector)
                    for w in 1:length(conf)
                        if conf[w]==WYE
                            push!(starpoints, (comp["bus"][w], comp["connections"][w][end]))
                        end
                    end
                else
                    if conf==WYE
                        push!(starpoints, (comp["bus"], comp["connections"][end]))
                    end
                end
            end
        end
    end

    # find grounded terminals
    grounded_terminals = [(b,t) for (b,bus) in data_eng["bus"] for t in bus["grounded"]]

    # assign initial neutral terminals
    nbts = [starpoints..., grounded_terminals...]

    # propagate neutral terminals through galvanic connections
    stack = deepcopy(nbts)
    while !isempty(stack)
        bt = pop!(stack)
        new_nbts = setdiff(bt_neighbours[bt], nbts)
        append!(stack, new_nbts)
        append!(nbts, new_nbts)
    end

    return nbts
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
                buses = isa(comp["bus"], String) ? [comp["bus"]] : comp["bus"]
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


"""
    add_bus_absolute_vbounds!(
        data_eng::Dict{String,Any};
        phase_lb_pu::Real=0.9,
        phase_ub_pu::Real=1.1,
        neutral_ub_pu::Real=0.3
    )::Dict{String,Any}

Adds absolute (i.e. indivdially, not between a pair of terminals) voltage bounds through the 'vm_lb' and 'vm_ub' property.
Bounds are specified in per unit, and automatically converted to SI units by calculating the voltage base.
If you change data_eng["settings"]["vbases_default"], the data model transformation will however produce inconsistent bounds in per unit.
Neutral terminals are automatically detected, and set to [0,phase_ub_pu*vbase].
"""
function add_bus_absolute_vbounds!(
    data_eng::Dict{String,Any};
    phase_lb_pu::Real=0.9,
    phase_ub_pu::Real=1.1,
    neutral_ub_pu::Real=0.3
    )::Dict{String,Any}

    nbts = _infer_neutral_terminals(data_eng)
    bus_vbase, _ = calc_voltage_bases(data_eng, data_eng["settings"]["vbases_default"])
    for (id, bus) in data_eng["bus"]
        vbase = bus_vbase[id]
        bus["vm_lb"] = [(id,t) in nbts ? 0.0 : phase_lb_pu*vbase for t in bus["terminals"]]
        bus["vm_ub"] = [(id,t) in nbts ? neutral_ub_pu*vbase : phase_ub_pu*vbase for t in bus["terminals"]]
    end

    return data_eng
end


"Generates pairwise bounds for oneport components."
function _generate_vm_pairs(connections::Vector, model::ConnConfig, lb::Real, ub::Real; delta_multiplier::Real=sqrt(3))
    vm_pair_ub = Tuple{Any,Any,Real}[]
    vm_pair_lb = Tuple{Any,Any,Real}[]

    if model==WYE
        n = connections[end]
        for p in connections[1:end-1]
            push!(vm_pair_lb, (p,n,lb))
            push!(vm_pair_ub, (p,n,ub))
        end
    elseif model==DELTA
        nphases = length(connections)
        if nphases==2
            push!(vm_pair_lb, (connections[1],connections[2],lb*delta_multiplier))
            push!(vm_pair_ub, (connections[1],connections[2],ub*delta_multiplier))
        elseif nphases==3
            for (p,q) in zip(connections, [connections[end], connections[1:end-1]...])
                push!(vm_pair_lb, (p,q,lb*delta_multiplier))
                push!(vm_pair_ub, (p,q,ub*delta_multiplier))
            end
        else
            error("Only 2-phase and 3-phase delta-connections components are supported.")
        end
    else
        error("The configuration $model is not supported.")
    end

    return vm_pair_lb, vm_pair_ub
end


"""
    add_unit_vbounds!(
        data_eng::Dict{String,Any};
        lb_pu::Real=0.9,
        ub_pu::Real=1.1,
        delta_multiplier::Real=sqrt(3),
        unit_comp_types::Vector{<:AbstractString}=["load", "generator", "storage", "pv"],
    )::Dict{String,Any}

Adds voltage bounds to the bus terminals to which units are connected.
'Units' in this context are all oneport component types specified by the argument 'unit_comp_types'.
Bounds are specified in per unit, and automatically converted to SI units by calculating the voltage base.
If you change data_eng["settings"]["vbases_default"], the data model transformation will however produce inconsistent bounds in per unit.
The delta multiplier controls the scaling of bounds of delta-connected units.
"""
function add_unit_vbounds!(
    data_eng::Dict{String,Any};
    lb_pu::Real=0.9,
    ub_pu::Real=1.1,
    delta_multiplier::Real=sqrt(3),
    unit_comp_types::Vector{<:AbstractString}=["load", "generator", "storage", "pv"],
    )::Dict{String,Any}

    bus_vbase, line_vbase = calc_voltage_bases(data_eng, data_eng["settings"]["vbases_default"])
    for comp_type in intersect(unit_comp_types, keys(data_eng))
        for (id, unit) in data_eng[comp_type]
            @assert _is_oneport_component(unit) "The $comp_type $id is not a oneport component."
            bus_id = unit["bus"]
            bus = data_eng["bus"][bus_id]
            unit_vm_pair_lb, unit_vm_pair_ub = _generate_vm_pairs(unit["connections"], unit["configuration"], lb_pu*bus_vbase[bus_id], ub_pu*bus_vbase[bus_id], delta_multiplier=delta_multiplier)
            bus_vm_pair_lb = haskey(bus, "vm_pair_lb") ? bus["vm_pair_lb"] : []
            bus_vm_pair_ub = haskey(bus, "vm_pair_ub") ? bus["vm_pair_ub"] : []
            bus["vm_pair_lb"] = unique(vcat(bus_vm_pair_lb, unit_vm_pair_lb))
            bus["vm_pair_ub"] = unique(vcat(bus_vm_pair_ub, unit_vm_pair_ub))
        end
    end

    return data_eng
end


"""
    add_bus_pn_pp_ng_vbounds!(data_eng::Dict{String,Any}, phase_terminals::Vector, neutral_terminal;
        pn_lb_pu::Union{Real,Missing}=missing,
        pn_ub_pu::Union{Real,Missing}=missing,
        pp_lb_pu::Union{Real,Missing}=missing,
        pp_ub_pu::Union{Real,Missing}=missing,
        ng_ub_pu::Union{Real,Missing}=missing,
    )::Dict{String,Any}

Adds symmetric phase-to-neutral and phase-to-phase voltage bounds when possible
for each bus through the three-phase bus syntax.
"""
function add_bus_pn_pp_ng_vbounds!(data_eng::Dict{String,Any}, phase_terminals::Vector, neutral_terminal;
    pn_lb_pu::Union{Real,Missing}=missing,
    pn_ub_pu::Union{Real,Missing}=missing,
    pp_lb_pu::Union{Real,Missing}=missing,
    pp_ub_pu::Union{Real,Missing}=missing,
    ng_ub_pu::Union{Real,Missing}=missing,
    )::Dict{String,Any}

    bus_vbase, _ = calc_voltage_bases(data_eng, data_eng["settings"]["vbases_default"])

    for (id,bus) in data_eng["bus"]
        vbase = bus_vbase[id]
        terminals = bus["terminals"]

        has_neutral = neutral_terminal in terminals
        if has_neutral
            bus["neutral"] = neutral_terminal
        end
        bus_phases = intersect(phase_terminals, terminals)
        has_phases = !isempty(bus_phases)
        if has_phases
            bus["phases"] = bus_phases
        end

        if has_neutral && !ismissing(ng_ub_pu)
            bus["vm_ng_ub"] = ng_ub_pu*vbase
        end
        if has_phases && has_neutral && !ismissing(pn_lb_pu)
            bus["vm_pn_lb"] = pn_lb_pu*vbase
        end
        if has_phases && has_neutral && !ismissing(pn_ub_pu)
            bus["vm_pn_ub"] = pn_ub_pu*vbase
        end
        if has_phases && length(bus_phases)>1 && !ismissing(pp_lb_pu)
            bus["vm_pp_lb"] = pp_lb_pu*vbase
        end
        if has_phases && length(bus_phases)>1 && !ismissing(pp_ub_pu)
            bus["vm_pp_ub"] = pp_ub_pu*vbase
        end
    end

    return data_eng
end


"""
    calc_start_voltage(
        data_math::Dict{String,Any};
        max_iter=Inf,
        epsilon::Number=1E-3
    )::Dict{Tuple{Int,Any},Union{Complex,Missing}}

Calculate no-load starting values for all bus-terminals pairs.
"""
function calc_start_voltage(
    data_math::Dict{String,Any};
    max_iter=Inf,
    epsilon::Number=1E-3
    )::Dict{Tuple{Int,Any},Union{Complex,Missing}}

    if haskey(data_math, "multinetwork")
        @assert !data_math["multinetwork"] "This method should be called on individual networks."
    end
    if haskey(data_math, "data_model")
        @assert data_math["data_model"]==MATHEMATICAL
    end

    node_links = Dict(vcat([[((bus["index"],t), []) for t in bus["terminals"]] for (_, bus) in data_math["bus"]]...))
    for (id, comp) in [data_math["branch"]..., data_math["switch"]...]
        f_bus   = comp["f_bus"]
        f_conns = comp["f_connections"]
        t_bus   = comp["t_bus"]
        t_conns = comp["t_connections"]

        arcs = []
        append!(arcs, [((f_bus,f_conns[k]), (t_bus,t_conns[k]), 1.0) for k in 1:length(f_conns)])
        append!(arcs, [((t_bus,t_conns[k]), (f_bus,f_conns[k]), 1.0) for k in 1:length(f_conns)])
        for (fr, to, scale) in arcs
            push!(node_links[fr], (to, scale))
        end
    end

    # initialize v_start for all bts to missing
    v_start = Dict{Tuple{Int,Any},Union{Complex,Missing}}((bus["index"],t)=>missing for (_, bus) in data_math["bus"] for t in bus["terminals"])

    for (_, bus) in data_math["bus"]
        # set fixed nodes
        if haskey(bus, "vm") && haskey(bus, "va")
            for (i,t) in enumerate(bus["terminals"])
                v_start[(bus["index"],t)] = bus["vm"][i]*exp(im*bus["va"][i])
            end

        end
        # set grounded nodes to zero
        for t in bus["terminals"][bus["grounded"]]
            v_start[(bus["index"],t)] = 0.0+im*0.0
        end
    end

    progress = round(sum(ismissing.(values(v_start)))/length(v_start)*100, digits=2)
    @debug "it. 0:\t$(100-progress)% specified at start"

    # propogate within zones and then between them across transformers as long as progress is made
    # a 'zone' in this context is a collection of buses which are galvanically connected
    # zones are connected through transformers
    count = 0
    stack = collect(keys(filter(kv->!ismissing(kv[2]), v_start)))
    while !isempty(stack) && count<max_iter
        count += 1

        # propogate all nodes in stack through switches/branches
        while !isempty(stack)
            node = pop!(stack)
            for (to, scale) in node_links[node]
                if ismissing(v_start[to])
                    v_start[to] = v_start[node]*scale
                    push!(stack, to)
                end
            end
        end

        # cross zones through transformers
        for (_, tr) in data_math["transformer"]
            f_bus = tr["f_bus"]; f_conns = tr["f_connections"]
            t_bus = tr["t_bus"]; t_conns = tr["t_connections"]
            tm_scale = calculate_tm_scale(tr, data_math["bus"]["$f_bus"], data_math["bus"]["$t_bus"])
            scale = (tm_scale*tr["polarity"]).*tr["tm_set"]
            if tr["configuration"]==WYE || (tr["configuration"]==DELTA && length(tr["tm_set"])==1)
                v_fr = Array{Union{Complex,Missing}}([v_start[(f_bus, t)] for t in f_conns])
                v_to = Array{Union{Complex,Missing}}([v_start[(t_bus, t)] for t in t_conns])
                # forward propagation
                if all((!).(ismissing.(v_fr))) && any(ismissing.(v_to))
                    N = length(v_fr)-1
                    Mpn = [LinearAlgebra.diagm(0=>ones(N)) fill(-1.0, N, 1)]
                    v_fr_pn = Mpn*v_fr
                    if all(ismissing.(v_to))
                        anchor_ind = length(v_to)
                        anchor_val = 0.0*im
                    else
                        (anchor_ind, anchor_val) = [(i, v) for (i, v) in enumerate(v_to) if !ismissing(v)][1]
                    end
                    v_to_prop = inv([Mpn; [i==anchor_ind ? 1 : 0 for i in 1:N+1]'])*[v_fr_pn./scale..., 0]

                    for (i,t) in enumerate(t_conns)
                        if ismissing(v_start[(t_bus,t)])
                            v_start[(t_bus,t)] = v_to_prop[i]
                            push!(stack, (t_bus,t))
                        end
                    end
                end
                # backward propagation
                if all((!).(ismissing.(v_to))) && any(ismissing.(v_fr))
                    N = length(v_fr)-1
                    Mpn = [LinearAlgebra.diagm(0=>ones(N)) fill(-1.0, N, 1)]
                    v_to_pn = Mpn*v_to
                    if all(ismissing.(v_fr))
                        anchor_ind = length(v_fr)
                        anchor_val = 0.0*im
                    else
                        (anchor_ind, anchor_val) = [(i, v) for (i, v) in enumerate(v_fr) if !ismissing(v)][1]
                    end
                    v_fr_prop = inv([Mpn; [i==anchor_ind ? 1 : 0 for i in 1:N+1]'])*[scale.*v_to_pn..., 0]

                    for (i,t) in enumerate(f_conns)
                        if ismissing(v_start[(f_bus,t)])
                            v_start[(f_bus,t)] = v_fr_prop[i]
                            push!(stack, (f_bus,t))
                        end
                    end
                end
            elseif tr["configuration"]==DELTA
                v_fr = Array{Union{Complex,Missing}}([v_start[(f_bus, t)] for t in f_conns])
                v_to = Array{Union{Complex,Missing}}([v_start[(t_bus, t)] for t in t_conns])
                # forward propagation
                if all((!).(ismissing.(v_fr))) && any(ismissing.(v_to))
                    if ismissing(v_to[end])
                        v_to[end] = 0.0*im
                    end
                    v_to_n = v_to[end]
                    M = _get_delta_transformation_matrix(length(v_fr))
                    v_to = [(M*v_fr)./scale.+v_to_n..., v_to_n]
                    for (i,t) in enumerate(t_conns)
                        if ismissing(v_start[(t_bus,t)])
                            v_start[(t_bus,t)] = v_to[i]
                            push!(stack, (t_bus,t))
                        end
                    end
                end
                # backward propagation
                if all((!).(ismissing.(v_to))) && any(ismissing.(v_fr))
                    v_to_p = v_to[1:end-1]; v_to_n = v_to[end]
                    M = _get_delta_transformation_matrix(length(v_fr))
                    Mp = [M[1:end-1,:]; fill(1.0, 1, length(v_fr))]
                    v_to_pn_scaled = (v_to_p.-v_to_n).*scale
                    v_fr = inv(Mp)*[v_to_pn_scaled[1:end-1]..., 0.0]
                    for (i,t) in enumerate(f_conns)
                        if ismissing(v_start[(f_bus,t)])
                            v_start[(f_bus,t)] = v_fr[i]
                            push!(stack, (f_bus,t))
                        end
                    end
                end
            end
        end
        progress = round(sum(ismissing.(values(v_start)))/length(v_start)*100, digits=2)
        @debug "it. $count:\t$progress% left to initialize at end"
    end

    # increment non-grounded zero values with epsilon
    for (k,v) in v_start
        (b,t) = k
        data_bus = data_math["bus"][string(b)]
        grounded_t = Dict(data_bus["terminals"].=>data_bus["grounded"])
        if !ismissing(v) && iszero(v) && !grounded_t[t]
            v_start[k] = convert(Complex, epsilon)
        end
    end

    return v_start
end


"""
    add_start_voltage!(
        data_math::Dict{String,Any};
        coordinates=:rectangular,
        uniform_v_start=missing,
        vr_default=0.0,
        vi_default=0.0,
        vm_default=0.0,
        va_default=0.0,
        epsilon::Number=1E-3,
    )::Dict{String,Any}

Adds start values for the voltage to the buses.
For a multinetwork data model, you can calculate the start voltages for a representative network through 'calc_start_voltage',
and pass the result as 'uniform_v_start' to use the same values for all networks and avoid recalculating it for each network.
The argument 'epsilon' controls the offset added to ungrounded terminals which would otherwise be set to zero.
"""
function add_start_voltage!(
    data_math::Dict{String,Any};
    coordinates=:rectangular,
    uniform_v_start=missing,
    vr_default=0.0,
    vi_default=0.0,
    vm_default=0.0,
    va_default=0.0,
    epsilon::Number=1E-3,
    )::Dict{String,Any}

    @assert data_math["data_model"]==MATHEMATICAL
    @assert coordinates in [:polar, :rectangular] "Legal values for the 'coordinates' argument are [:polar,:rectangular], not :$coordinates."

    is_mn = haskey(data_math, "multinetwork")

    for (nw,dm) in (is_mn ? data_math["nw"] : [("", data_math)])
        if ismissing(uniform_v_start)
            v_start = calc_start_voltage(dm, epsilon=epsilon)
        else
            v_start = uniform_v_start
        end
        for (_, bus) in dm["bus"]
            index = bus["index"]
            if coordinates==:rectangular
                bus["vr_start"] = [ismissing(v_start[(index, t)]) ? vr_default : real(v_start[(index, t)]) for t in bus["terminals"]]
                bus["vi_start"] = [ismissing(v_start[(index, t)]) ? vi_default : imag(v_start[(index, t)]) for t in bus["terminals"]]
            elseif coordinates==:polar
                bus["vm_start"] = [ismissing(v_start[(index, t)]) ? vm_default : abs(v_start[(index, t)]) for t in bus["terminals"]]
                bus["va_start"] = [ismissing(v_start[(index, t)]) ? va_default : angle(v_start[(index, t)]) for t in bus["terminals"]]
            end
        end
    end

    return data_math
end


"""
    add_start_vrvi!(data_math::Dict{String,Any}; kwargs...)

    Short-hand for [`add_start_voltage`](@ref add_start_voltage) with rectangular coordinates (coordinates=:rectangular).
"""
add_start_vrvi!(data_math::Dict{String,Any}; kwargs...) = add_start_voltage!(data_math, coordinates=:rectangular, kwargs...)


"""
    add_start_vmva!(data_math::Dict{String,Any}; kwargs...)

Short-hand for [`add_start_voltage`](@ref add_start_voltage) with polar coordinates (coordinates=:polar).
"""
add_start_vmva!(data_math::Dict{String,Any}; kwargs...) = add_start_voltage!(data_math, coordinates=:polar, kwargs...)
