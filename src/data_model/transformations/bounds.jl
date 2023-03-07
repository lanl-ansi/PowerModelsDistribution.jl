"""
    adjust_small_line_impedances!(data::Dict{String,<:Any}; min_impedance_val::Real=1e-2, replace_impedance_val::Real=0.0)

Replaces impedances (rs, xs) on lines, linecodes, and switches lower than `min_impedance_val` with `replace_impedance_val`.
"""
function adjust_small_line_impedances!(data::Dict{String,<:Any}; min_impedance_val::Real=1e-2, replace_impedance_val::Real=0.0)
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_adjust_small_line_impedances!, data; apply_to_subnetworks=true, min_impedance_val=min_impedance_val, replace_impedance_val=replace_impedance_val)
end


"""
    _adjust_small_line_impedances!(data_eng::Dict{String,<:Any}; min_impedance_val::Real=1e-2, replace_impedance_val::Real=0.0)

Replaces impedances (rs, xs) on lines, linecodes, and switches lower than `min_impedance_val` with `replace_impedance_val`.
"""
function _adjust_small_line_impedances!(data_eng::Dict{String,<:Any}; min_impedance_val::Real=1e-2, replace_impedance_val::Real=0.0)
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
    adjust_small_line_admittances!(data::Dict{String,<:Any}; min_admittance_val::Real=1e-2, replace_admittance_val::Real=0.0)

Replaces admittances (g_fr, g_to, b_fr, b_to) on lines, linecodes, and switches lower than `min_admittance_val` with `replace_admittance_val`.
"""
function adjust_small_line_admittances!(data::Dict{String,<:Any}; min_admittance_val::Real=1e-2, replace_admittance_val::Real=0.0)
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_adjust_small_line_admittances!, data; apply_to_subnetworks=true, min_admittance_val=min_admittance_val, replace_admittance_val=replace_admittance_val)
end


"""
    _adjust_small_line_admittances!(data_eng::Dict{String,<:Any}; min_admittance_val::Real=1e-2, replace_admittance_val::Real=0.0)

Replaces admittances (g_fr, g_to, b_fr, b_to) on lines, linecodes, and switches lower than `min_admittance_val` with `replace_admittance_val`.
"""
function _adjust_small_line_admittances!(data_eng::Dict{String,<:Any}; min_admittance_val::Real=1e-2, replace_admittance_val::Real=0.0)
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
    adjust_small_line_lengths!(data::Dict{String,<:Any}; min_length_val::Real=25.0, replace_length_val::Real=0.0)

Replaces length on lines, switches lower than `min_length_val` with `replace_length_val`.
"""
function adjust_small_line_lengths!(data::Dict{String,<:Any}; min_length_val::Real=25.0, replace_length_val::Real=0.0)
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_adjust_small_line_lengths!, data; apply_to_subnetworks=true, min_length_val=min_length_val, replace_length_val=replace_length_val)
end


"""
    _adjust_small_line_lengths!(data_eng::Dict{String,<:Any}; min_length_val::Real=25.0, replace_length_val::Real=0.0)

Replaces length on lines, switches lower than `min_length_val` with `replace_length_val`.
"""
function _adjust_small_line_lengths!(data_eng::Dict{String,<:Any}; min_length_val::Real=25.0, replace_length_val::Real=0.0)
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
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_apply_voltage_bounds!, data; apply_to_subnetworks=true, vm_lb=vm_lb, vm_ub=vm_ub)
end



"""
    _apply_voltage_bounds!(data_eng::Dict{String,<:Any}; vm_lb::Union{Real,Missing}=0.9, vm_ub::Union{Real,Missing}=1.1)

add voltage bounds to all buses based on per-unit upper (`vm_ub`) and lower (`vm_lb`), scaled by the bus's voltage based
"""
function _apply_voltage_bounds!(data_eng::Dict{String,<:Any}; vm_lb::Union{Real,Missing}=0.9, vm_ub::Union{Real,Missing}=1.1, exclude::Vector{String}=!isempty(get(data_eng, "voltage_source", Dict())) ? String[x.second["bus"] for x in data_eng["voltage_source"]] : String[])
    (bus_vbases, _) = calc_eng_voltage_bases(data_eng, data_eng["settings"]["vbases_default"])
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
    adjust_line_limits!(data::Dict{String,<:Any}, mult::Real)

Multiplies limits (`sm_ub` and/or `cm_ub`) on line objects (`line`, `linecode`, `switch`) by a multiplier `mult`
"""
function adjust_line_limits!(data::Dict{String,<:Any}, mult::Real)
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_adjust_line_limits!, data, mult; apply_to_subnetworks=true)
end


"""
    _adjust_line_limits!(data_eng::Dict{String,<:Any}, mult::Real)

Multiplies limits (`sm_ub` and/or `cm_ub`) on line objects (`line`, `linecode`, `switch`) by a multiplier `mult`
"""
function _adjust_line_limits!(data_eng::Dict{String,<:Any}, mult::Real)
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
    adjust_transformer_limits!(data::Dict{String,<:Any}, mult::Real)

Multiplies limits (`sm_ub` and/or `cm_ub`) on transformer objects by a multiplier `mult`
"""
function adjust_transformer_limits!(data::Dict{String,<:Any}, mult::Real)
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_adjust_transformer_limits!, data, mult; apply_to_subnetworks=true)
end


"""
    _adjust_transformer_limits!(data_eng::Dict{String,<:Any}, mult::Real)

Multiplies limits (`sm_ub` and/or `cm_ub`) on transformer objects by a multiplier `mult`
"""
function _adjust_transformer_limits!(data_eng::Dict{String,<:Any}, mult::Real)
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
