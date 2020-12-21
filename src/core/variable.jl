
function _get_conductor_indicator(comp::Dict{String,<:Any})
    if haskey(comp, "terminals")
        return "terminals"
    elseif haskey(comp, "connections")
        return "connections"
    elseif haskey(comp, "f_connections")
        return "f_connections"
    else
        return ""
    end
end


function comp_start_value(comp::Dict{String,<:Any}, key::String, conductor::Int, default)
    cond_ind = _get_conductor_indicator(comp)
    if haskey(comp, key) && !isempty(cond_ind)
        return comp[key][findfirst(isequal(conductor), comp[cond_ind])]
    else
        return default
    end
end


function comp_start_value(comp::Dict{String,<:Any}, key::String, default)
    return _PM.comp_start_value(comp, key, default)
end


""
function variable_mc_bus_voltage_angle(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    va = var(pm, nw)[:va] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_va_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "va_start", t, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    report && _IM.sol_component_value(pm, nw, :bus, :va, ids(pm, nw, :bus), va)
end


""
function variable_mc_bus_voltage_magnitude_only(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    vm = var(pm, nw)[:vm] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_vm_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "vm_start", t, 1.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in ref(pm, nw, :bus)
            for (idx, t) in enumerate(terminals[i])
                if haskey(bus, "vmin")
                    set_lower_bound(vm[i][t], bus["vmin"][idx])
                end
                if haskey(bus, "vmax")
                    set_upper_bound(vm[i][t], bus["vmax"][idx])
                end
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :bus, :vm, ids(pm, nw, :bus), vm)
end

""
function variable_mc_bus_voltage_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    vr = var(pm, nw)[:vr] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_vr_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "vr_start", t, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in ref(pm, nw, :bus)
            if haskey(bus, "vmax")
                for (idx,t) in enumerate(terminals[i])
                    set_lower_bound(vr[i][t], -bus["vmax"][idx])
                    set_upper_bound(vr[i][t],  bus["vmax"][idx])
                end
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :bus, :vr, ids(pm, nw, :bus), vr)
end

""
function variable_mc_bus_voltage_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    vi = var(pm, nw)[:vi] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_vi_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "vi_start", t, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in ref(pm, nw, :bus)
            if haskey(bus, "vmax")
                for (idx,t) in enumerate(terminals[i])
                    set_lower_bound(vi[i][t], -bus["vmax"][idx])
                    set_upper_bound(vi[i][t],  bus["vmax"][idx])
                end
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :bus, :vi, ids(pm, nw, :bus), vi)
end


"branch flow variables, delegated back to PowerModels"
function variable_mc_branch_power(pm::_PM.AbstractPowerModel; kwargs...)
    variable_mc_branch_power_real(pm; kwargs...)
    variable_mc_branch_power_imaginary(pm; kwargs...)
end


"variable: `p[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_power_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)
    p = var(pm, nw)[:p] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_p_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "p_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs)
            smax = _calc_branch_power_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
            for (idx, c) in enumerate(connections[(l,i,j)])
                set_upper_bound(p[(l,i,j)][c],  smax[idx])
                set_lower_bound(p[(l,i,j)][c], -smax[idx])
            end
        end
    end

    for (l,branch) in ref(pm, nw, :branch)
        if haskey(branch, "pf_start")
            f_idx = (l, branch["f_bus"], branch["t_bus"])
            for (idx,c) in enumerate(connections[f_idx])
                JuMP.set_start_value(p[f_idx][c], branch["pf_start"][idx])
            end
        end
        if haskey(branch, "pt_start")
            t_idx = (l, branch["t_bus"], branch["f_bus"])
            for (idx,c) in enumerate(connections[t_idx])
                JuMP.set_start_value(p[t_idx][c], branch["pt_start"][idx])
            end
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :branch, :pf, :pt, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), p)
end

"variable: `q[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)
    q = var(pm, nw)[:q] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_q_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "q_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs)
            smax = _calc_branch_power_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
            for (idx, c) in enumerate(connections[(l,i,j)])
                set_upper_bound(q[(l,i,j)][c],  smax[idx])
                set_lower_bound(q[(l,i,j)][c], -smax[idx])
            end
        end
    end

    for (l,branch) in ref(pm, nw, :branch)
        if haskey(branch, "qf_start")
            f_idx = (l, branch["f_bus"], branch["t_bus"])
            for (idx,c) in enumerate(connections[f_idx])
                JuMP.set_start_value(q[f_idx][c], branch["qf_start"][idx])
            end
        end
        if haskey(branch, "qt_start")
            t_idx = (l, branch["t_bus"], branch["f_bus"])
            for (idx,c) in enumerate(connections[t_idx])
                JuMP.set_start_value(q[t_idx][c], branch["qt_start"][idx])
            end
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :branch, :qf, :qt, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), q)
end


"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_current_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)
    cr = var(pm, nw)[:cr] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_cr_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "cr_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs)
            cmax = _calc_branch_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
            for (idx,c) in enumerate(connections[(l,i,j)])
                set_upper_bound(cr[(l,i,j)][c],  cmax[idx])
                set_lower_bound(cr[(l,i,j)][c], -cmax[idx])
            end
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :branch, :cr_fr, :cr_to, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), cr)
end


"variable: `ci[l,i,j] ` for `(l,i,j)` in `arcs`"
function variable_mc_branch_current_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)
    ci = var(pm, nw)[:ci] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_ci_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "ci_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs)
            cmax = _calc_branch_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
            for (idx,c) in enumerate(connections[(l,i,j)])
                set_upper_bound(ci[(l,i,j)][c],  cmax[idx])
                set_lower_bound(ci[(l,i,j)][c], -cmax[idx])
            end
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :branch, :ci_fr, :ci_to, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), ci)
end


"variable: `csr[l]` for `l` in `branch`"
function variable_mc_branch_current_series_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)
    csr = var(pm, nw)[:csr] = Dict(l => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_csr_$(l)",
            start = comp_start_value(ref(pm, nw, :branch, l), "csr_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_from)
            cmax = _calc_branch_series_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
            for (idx,c) in enumerate(connections[(l,i,j)])
                set_upper_bound(csr[l][c],  cmax[idx])
                set_lower_bound(csr[l][c], -cmax[idx])
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :branch, :csr_fr, ids(pm, nw, :branch), csr)
end


"variable: `csi[l]` for `l` in `branch`"
function variable_mc_branch_current_series_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)
    csi = var(pm, nw)[:csi] = Dict(l => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_csi_$(l)",
            start = comp_start_value(ref(pm, nw, :branch, l), "csi_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_from)
            cmax = _calc_branch_series_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
            for (idx,c) in enumerate(connections[(l,i,j)])
                set_upper_bound(csi[l][c],  cmax[idx])
                set_lower_bound(csi[l][c], -cmax[idx])
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :branch, :csi_fr, ids(pm, nw, :branch), csi)
end


"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_transformer_current_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_transformer) for ((l,i,j), connections) in entry)
    cr = var(pm, nw)[:crt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_crt_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :transformer, l), "cr_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_trans)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_from_trans)
            trans = ref(pm, nw, :transformer, l)
            f_bus = ref(pm, nw, :bus, i)
            t_bus = ref(pm, nw, :bus, j)
            cmax_fr, cmax_to = _calc_transformer_current_max_frto(trans, f_bus, t_bus)
            set_lower_bound(cr[(l,i,j)], -cmax_fr)
            set_upper_bound(cr[(l,i,j)],  cmax_fr)
            set_lower_bound(cr[(l,j,i)], -cmax_to)
            set_upper_bound(cr[(l,j,i)],  cmax_to)
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :transformer, :cr_fr, :cr_to, ref(pm, nw, :arcs_from_trans), ref(pm, nw, :arcs_to_trans), cr)
end


"variable: `ci[l,i,j] ` for `(l,i,j)` in `arcs`"
function variable_mc_transformer_current_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_transformer) for ((l,i,j), connections) in entry)
    ci = var(pm, nw)[:cit] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_cit_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :transformer, l), "ci_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_trans)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_from_trans)
            trans = ref(pm, nw, :transformer, l)
            f_bus = ref(pm, nw, :bus, i)
            t_bus = ref(pm, nw, :bus, j)
            cmax_fr, cmax_to = _calc_transformer_current_max_frto(trans, f_bus, t_bus)
            set_lower_bound(ci[(l,i,j)], -cmax_fr)
            set_upper_bound(ci[(l,i,j)],  cmax_fr)
            set_lower_bound(ci[(l,j,i)], -cmax_to)
            set_upper_bound(ci[(l,j,i)],  cmax_to)
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :transformer, :ci_fr, :ci_to, ref(pm, nw, :arcs_from_trans), ref(pm, nw, :arcs_to_trans), ci)
end


"variable: `w[i] >= 0` for `i` in `buses"
function variable_mc_bus_voltage_magnitude_sqr(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    w = var(pm, nw)[:w] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_w_$(i)",
            lower_bound = 0.0,
            start = comp_start_value(ref(pm, nw, :bus, i), "w_start", 1.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for i in ids(pm, nw, :bus)
            bus = ref(pm, nw, :bus, i)
            for (idx, t) in enumerate(terminals[i])
                set_upper_bound(w[i][t], max(bus["vmin"][idx]^2, bus["vmax"][idx]^2))
                if bus["vmin"][idx] > 0
                    set_lower_bound(w[i][t], bus["vmin"][idx]^2)
                end
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :bus, :w, ids(pm, nw, :bus), w)
end


"variables for modeling storage units, includes grid injection and internal variables"
function variable_mc_storage_power(pm::_PM.AbstractPowerModel; kwargs...)
    variable_mc_storage_power_real(pm; kwargs...)
    variable_mc_storage_power_imaginary(pm; kwargs...)
    variable_mc_storage_power_control_imaginary(pm; kwargs...)
    variable_mc_storage_current(pm; kwargs...)
    _PM.variable_storage_energy(pm; kwargs...)
    _PM.variable_storage_charge(pm; kwargs...)
    _PM.variable_storage_discharge(pm; kwargs...)
end


""
function variable_mc_storage_power_mi(pm::_PM.AbstractPowerModel; relax::Bool=false, kwargs...)
    variable_mc_storage_power_on_off(pm; kwargs...)
    variable_mc_storage_power_control_imaginary(pm; kwargs...)
    variable_mc_storage_current(pm; kwargs...)
    variable_mc_storage_indicator(pm; relax=relax, kwargs...)
    _PM.variable_storage_energy(pm; kwargs...)
    _PM.variable_storage_charge(pm; kwargs...)
    _PM.variable_storage_discharge(pm; kwargs...)
    _PM.variable_storage_complementary_indicator(pm; relax=relax, kwargs...)
end


""
function variable_mc_storage_power_mi_on_off(pm::_PM.AbstractPowerModel; relax::Bool=false, kwargs...)
    variable_mc_storage_power_real_on_off(pm; kwargs...)
    variable_mc_storage_power_imaginary_on_off(pm; kwargs...)
    variable_mc_storage_current(pm; kwargs...)
    variable_mc_storage_power_control_imaginary_on_off(pm; kwargs...)
    _PM.variable_storage_energy(pm; kwargs...)
    _PM.variable_storage_charge(pm; kwargs...)
    _PM.variable_storage_discharge(pm; kwargs...)
    _PM.variable_storage_complementary_indicator(pm; relax=relax, kwargs...)
end


"do nothing by default but some formulations require this"
function variable_mc_storage_current(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
end


"""
a reactive power slack variable that enables the storage device to inject or
consume reactive power at its connecting bus, subject to the injection limits
of the device.
"""
function variable_mc_storage_power_control_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    qsc = var(pm, nw)[:qsc] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :storage)], base_name="$(nw)_qsc_$(i)",
        start = _PM.comp_start_value(ref(pm, nw, :storage, i), "qsc_start")
    )

    if bounded
        inj_lb, inj_ub = _PM.ref_calc_storage_injection_bounds(ref(pm, nw, :storage), ref(pm, nw, :bus))
        for (i,storage) in ref(pm, nw, :storage)
            if !isinf(sum(inj_lb[i])) || haskey(storage, "qmin")
                set_lower_bound(qsc[i], max(sum(inj_lb[i]), sum(get(storage, "qmin", -Inf))))
            end
            if !isinf(sum(inj_ub[i])) || haskey(storage, "qmax")
                set_upper_bound(qsc[i], min(sum(inj_ub[i]), sum(get(storage, "qmax",  Inf))))
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :storage, :qsc, ids(pm, nw, :storage), qsc)
end


"""
a reactive power slack variable that enables the storage device to inject or
consume reactive power at its connecting bus, subject to the injection limits
of the device.
"""
function variable_mc_storage_power_control_imaginary_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    qsc = var(pm, nw)[:qsc] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :storage)], base_name="$(nw)_qsc_$(i)",
        start = _PM.comp_start_value(ref(pm, nw, :storage, i), "qsc_start")
    )

    if bounded
        inj_lb, inj_ub = _PM.ref_calc_storage_injection_bounds(ref(pm, nw, :storage), ref(pm, nw, :bus))
        for (i,storage) in ref(pm, nw, :storage)
            if !isinf(sum(inj_lb[i])) || haskey(storage, "qmin")
                lb = max(sum(inj_lb[i]), sum(get(storage, "qmin", -Inf)))
                set_lower_bound(qsc[i], min(lb, 0.0))
            end
            if !isinf(sum(inj_ub[i])) || haskey(storage, "qmax")
                ub = min(sum(inj_ub[i]), sum(get(storage, "qmax", Inf)))
                set_upper_bound(qsc[i], max(ub, 0.0))
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :storage, :qsc, ids(pm, nw, :storage), qsc)
end



""
function variable_mc_storage_power_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i,strg) in ref(pm, nw, :storage))
    ps = var(pm, nw)[:ps] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_ps_$(i)",
            start = comp_start_value(ref(pm, nw, :storage, i), "ps_start", c, 0.0)
        ) for i in ids(pm, nw, :storage)
    )

    if bounded
        flow_lb, flow_ub = ref_calc_storage_injection_bounds(ref(pm, nw, :storage), ref(pm, nw, :bus))
        for i in ids(pm, nw, :storage)
            for (idx, c) in enumerate(connections[i])
                if !isinf(flow_lb[i][idx])
                    set_lower_bound(ps[i][c], flow_lb[i][idx])
                end
                if !isinf(flow_ub[i][idx])
                    set_upper_bound(ps[i][c], flow_ub[i][idx])
                end
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :storage, :ps, ids(pm, nw, :storage), ps)
end


""
function variable_mc_storage_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i,strg) in ref(pm, nw, :storage))
    qs = var(pm, nw)[:qs] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_qs_$(i)",
            start = comp_start_value(ref(pm, nw, :storage, i), "qs_start", c, 0.0)
        ) for i in ids(pm, nw, :storage)
    )

    if bounded
        flow_lb, flow_ub = ref_calc_storage_injection_bounds(ref(pm, nw, :storage), ref(pm, nw, :bus))
        for i in ids(pm, nw, :storage)
            for (idx, c) in enumerate(connections[i])
                if !isinf(flow_lb[i][idx])
                    set_lower_bound(qs[i][c], flow_lb[i][idx])
                end
                if !isinf(flow_ub[i][idx])
                    set_upper_bound(qs[i][c], flow_ub[i][idx])
                end
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :storage, :qs, ids(pm, nw, :storage), qs)
end


"generates variables for both `active` and `reactive` slack at each bus"
function variable_mc_slack_bus_power(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, kwargs...)
    variable_mc_slack_bus_power_real(pm; nw=nw, kwargs...)
    variable_mc_slack_bus_power_imaginary(pm; nw=nw, kwargs...)
end


""
function variable_mc_slack_bus_power_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, report::Bool=true)
    terminals = Dict(i => ref(pm, nw, :bus, i)["terminals"] for i in ids(pm, nw, :bus))
    p_slack = var(pm, nw)[:p_slack] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_p_slack_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "p_slack_start", t, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    report && _IM.sol_component_value(pm, nw, :bus, :p_slack, ids(pm, nw, :bus), p_slack)
end


""
function variable_mc_slack_bus_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, report::Bool=true)
    terminals = Dict(i => ref(pm, nw, :bus, i)["terminals"] for i in ids(pm, nw, :bus))
    q_slack = var(pm, nw)[:q_slack] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_q_slack_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "q_slack_start", t, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    report && _IM.sol_component_value(pm, nw, :bus, :q_slack, ids(pm, nw, :bus), q_slack)
end


"Creates variables for both `active` and `reactive` power flow at each transformer."
function variable_mc_transformer_power(pm::_PM.AbstractPowerModel; kwargs...)
    variable_mc_transformer_power_real(pm; kwargs...)
    variable_mc_transformer_power_imaginary(pm; kwargs...)
end


"Create variables for the active power flowing into all transformer windings."
function variable_mc_transformer_power_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_transformer) for ((l,i,j), connections) in entry)
    pt = var(pm, nw)[:pt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_pt_$((l,i,j))",
        ) for (l,i,j) in ref(pm, nw, :arcs_trans)
    )

    if bounded
        for arc in ref(pm, nw, :arcs_from_trans)
            (l,i,j) = arc
            rate_a_fr, rate_a_to = _calc_transformer_power_ub_frto(ref(pm, nw, :transformer, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
            for (idx, (fc, tc)) in enumerate(zip(connections[(l,i,j)], connections[(l,j,i)]))
                set_lower_bound(pt[(l,i,j)][fc], -rate_a_fr[idx])
                set_upper_bound(pt[(l,i,j)][fc],  rate_a_fr[idx])
                set_lower_bound(pt[(l,j,i)][tc], -rate_a_to[idx])
                set_upper_bound(pt[(l,j,i)][tc],  rate_a_to[idx])
            end
        end
    end

    for (l,transformer) in ref(pm, nw, :transformer)
        if haskey(transformer, "pf_start")
            f_idx = (l, transformer["f_bus"], transformer["t_bus"])
            for (idx, c) in enumerate(connections[f_idx])
                JuMP.set_start_value(pt[f_idx][c], transformer["pf_start"][idx])
            end
        end
        if haskey(transformer, "pt_start")
            t_idx = (l, transformer["t_bus"], transformer["f_bus"])
            for (idx, c) in enumerate(connections[t_idx])
                JuMP.set_start_value(pt[t_idx][c], transformer["pt_start"][idx])
            end
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :transformer, :pf, :pt, ref(pm, nw, :arcs_from_trans), ref(pm, nw, :arcs_to_trans), pt)
end


"Create variables for the reactive power flowing into all transformer windings."
function variable_mc_transformer_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_transformer) for ((l,i,j), connections) in entry)
    qt = var(pm, nw)[:qt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_qt_$((l,i,j))",
            start = 0.0
        ) for (l,i,j) in ref(pm, nw, :arcs_trans)
    )

    if bounded
        for arc in ref(pm, nw, :arcs_from_trans)
            (l,i,j) = arc
            rate_a_fr, rate_a_to = _calc_transformer_power_ub_frto(ref(pm, nw, :transformer, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))

            for (idx, (fc,tc)) in enumerate(zip(connections[(l,i,j)], connections[(l,j,i)]))
                set_lower_bound(qt[(l,i,j)][fc], -rate_a_fr[idx])
                set_upper_bound(qt[(l,i,j)][fc],  rate_a_fr[idx])
                set_lower_bound(qt[(l,j,i)][tc], -rate_a_to[idx])
                set_upper_bound(qt[(l,j,i)][tc],  rate_a_to[idx])
            end
        end
    end

    for (l,transformer) in ref(pm, nw, :transformer)
        if haskey(transformer, "qf_start")
            f_idx = (l, transformer["f_bus"], transformer["t_bus"])
            JuMP.set_start_value(qt[f_idx], transformer["qf_start"])
        end
        if haskey(transformer, "qt_start")
            t_idx = (l, transformer["t_bus"], transformer["f_bus"])
            JuMP.set_start_value(qt[t_idx], transformer["qt_start"])
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :transformer, :qf, :qt, ref(pm, nw, :arcs_from_trans), ref(pm, nw, :arcs_to_trans), qt)
end


"Create tap variables."
function variable_mc_oltc_transformer_tap(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    # when extending to 4-wire, this should iterate only over the phase conductors
    p_oltc_ids = [id for (id,trans) in ref(pm, nw, :transformer) if !all(trans["tm_fix"])]
    tap = var(pm, nw)[:tap] = Dict(i => JuMP.@variable(pm.model,
        [p in 1:length(ref(pm,nw,:transformer,i,"f_connections"))],
        base_name="$(nw)_tm_$(i)",
        start=ref(pm, nw, :transformer, i, "tm_set")[p]
    ) for i in p_oltc_ids)
    if bounded
        for tr_id in p_oltc_ids, p in 1:length(ref(pm,nw,:transformer,tr_id,"f_connections"))
            set_lower_bound(var(pm, nw)[:tap][tr_id][p], ref(pm, nw, :transformer, tr_id, "tm_lb")[p])
            set_upper_bound(var(pm, nw)[:tap][tr_id][p], ref(pm, nw, :transformer, tr_id, "tm_ub")[p])
        end
    end

    report && _IM.sol_component_value(pm, nw, :transformer, :tap, ids(pm, nw, :transformer), tap)
end


"""
Create a dictionary with values of type Any for the load.
Depending on the load model, this can be a parameter or a NLexpression.
These will be inserted into KCL.
"""
function variable_mc_load_power(pm::_PM.AbstractPowerModel; nw=pm.cnw, bounded::Bool=true, report::Bool=true)
    var(pm, nw)[:pd] = Dict{Int, Any}()
    var(pm, nw)[:qd] = Dict{Int, Any}()
    var(pm, nw)[:pd_bus] = Dict{Int, Any}()
    var(pm, nw)[:qd_bus] = Dict{Int, Any}()
end


"Create variables for demand status"
function variable_mc_load_indicator(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, relax::Bool=false, report::Bool=true)
    if relax
        z_demand = var(pm, nw)[:z_demand] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :load)], base_name="$(nw)_z_demand",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(ref(pm, nw, :load, i), "z_demand_on_start", 1.0)
        )
    else
        z_demand = var(pm, nw)[:z_demand] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :load)], base_name="$(nw)_z_demand",
            binary = true,
            start = comp_start_value(ref(pm, nw, :load, i), "z_demand_on_start", 1.0)
        )
    end

    # expressions for pd and qd
    pd = var(pm, nw)[:pd] = Dict(i => var(pm, nw)[:z_demand][i].*ref(pm, nw, :load, i)["pd"] for i in ids(pm, nw, :load))
    qd = var(pm, nw)[:qd] = Dict(i => var(pm, nw)[:z_demand][i].*ref(pm, nw, :load, i)["qd"] for i in ids(pm, nw, :load))

    report && _IM.sol_component_value(pm, nw, :load, :status, ids(pm, nw, :load), z_demand)
    report && _IM.sol_component_value(pm, nw, :load, :pd, ids(pm, nw, :load), pd)
    report && _IM.sol_component_value(pm, nw, :load, :qd, ids(pm, nw, :load), qd)
end


"Create variables for shunt status"
function variable_mc_shunt_indicator(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, relax=false, report::Bool=true)
    if relax
        z_shunt = var(pm, nw)[:z_shunt] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :shunt)], base_name="$(nw)_z_shunt",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(ref(pm, nw, :shunt, i), "z_shunt_on_start", 1.0)
        )
    else
        z_shunt = var(pm, nw)[:z_shunt] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :shunt)], base_name="$(nw)_z_shunt",
            binary=true,
            start = comp_start_value(ref(pm, nw, :shunt, i), "z_shunt_on_start", 1.0)
        )
    end

    report && _IM.sol_component_value(pm, nw, :shunt, :status, ids(pm, nw, :shunt), z_shunt)
end


"Create variables for bus status"
function variable_mc_bus_voltage_indicator(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_voltage = var(pm, nw)[:z_voltage] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :bus)], base_name="$(nw)_z_voltage",
            binary = true,
            start = comp_start_value(ref(pm, nw, :bus, i), "z_voltage_start", 1.0)
        )
    else
        z_voltage =var(pm, nw)[:z_voltage] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :bus)], base_name="$(nw)_z_voltage",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(ref(pm, nw, :bus, i), "z_voltage_start", 1.0)
        )
    end

    report && _IM.sol_component_value(pm, nw, :bus, :status, ids(pm, nw, :bus), z_voltage)
end


"Create variables for generator status"
function variable_mc_gen_indicator(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_gen = var(pm, nw)[:z_gen] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :gen)], base_name="$(nw)_z_gen",
            binary = true,
            start = comp_start_value(ref(pm, nw, :gen, i), "z_gen_start", 1.0)
        )
    else
        z_gen = var(pm, nw)[:z_gen] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :gen)], base_name="$(nw)_z_gen",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(ref(pm, nw, :gen, i), "z_gen_start", 1.0)
        )
    end

    report && _IM.sol_component_value(pm, nw, :gen, :gen_status, ids(pm, nw, :gen), z_gen)
end


"Create variables for storage status"
function variable_mc_storage_indicator(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_storage = var(pm, nw)[:z_storage] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :storage)], base_name="$(nw)-z_storage",
            binary = true,
            start = comp_start_value(ref(pm, nw, :storage, i), "z_storage_start", 1.0)
        )
    else
        z_storage = var(pm, nw)[:z_storage] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :storage)], base_name="$(nw)_z_storage",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(ref(pm, nw, :storage, i), "z_storage_start", 1.0)
        )
    end

    report && _IM.sol_component_value(pm, nw, :storage, :status, ids(pm, nw, :storage), z_storage)
end


"Create variables for `active` and `reactive` storage injection"
function variable_mc_storage_power_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, kwargs...)
    variable_mc_storage_power_real_on_off(pm; nw=nw, kwargs...)
    variable_mc_storage_power_imaginary_on_off(pm; nw=nw, kwargs...)
end


"Create variables for `active` storage injection"
function variable_mc_storage_power_real_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i,strg) in ref(pm, nw, :storage))
    ps = var(pm, nw)[:ps] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_ps_$(i)",
        start = comp_start_value(ref(pm, nw, :storage, i), "ps_start", c, 0.0)
    ) for i in ids(pm, nw, :storage))

    if bounded
        for (i, strg) in ref(pm, nw, :storage)
            for (idx, c) in enumerate(connections[i])
                inj_lb, inj_ub = _PM.ref_calc_storage_injection_bounds(ref(pm, nw, :storage), ref(pm, nw, :bus), idx)

                set_lower_bound(ps[i][c], inj_lb[i])
                set_upper_bound(ps[i][c], inj_ub[i])
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :storage, :ps, ids(pm, nw, :storage), ps)
end


"Create variables for `reactive` storage injection"
function variable_mc_storage_power_imaginary_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i,strg) in ref(pm, nw, :storage))
    qs = var(pm, nw)[:qs] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_qs_$(i)",
        start = comp_start_value(ref(pm, nw, :storage, i), "qs_start", c, 0.0)
    ) for i in ids(pm, nw, :storage))

    if bounded
        for (i, strg) in ref(pm, nw, :storage)
            if haskey(strg, "qmin")
                for (idx, c) in enumerate(connections[i])
                    set_lower_bound(qs[i][c], strg["qmin"][idx])
                end
            end

            if haskey(strg, "qmax")
                for (idx, c) in enumerate(connections[i])
                    set_upper_bound(qs[i][c], strg["qmax"][idx])
                end
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :storage, :qs, ids(pm, nw, :storage), qs)
end


"voltage variable magnitude squared (relaxed form)"
function variable_mc_bus_voltage_magnitude_sqr_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    w = var(pm, nw)[:w] = Dict(i => JuMP.@variable(pm.model,
        [t in terminals[i]], base_name="$(nw)_w_$(i)",
        start = comp_start_value(ref(pm, nw, :bus, i), "w_start", t, 1.001)
    ) for i in ids(pm, nw, :bus))

    if bounded
        for (i, bus) in ref(pm, nw, :bus)
            for (idx, t) in enumerate(terminals[i])
                set_lower_bound(w[i][t], 0.0)

                if haskey(bus, "vmax")
                    set_upper_bound(w[i][t], bus["vmax"][idx]^2)
                end
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :bus, :w, ids(pm, nw, :bus), w)
end


"on/off voltage magnitude variable"
function variable_mc_bus_voltage_magnitude_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    vm = var(pm, nw)[:vm] = Dict(i => JuMP.@variable(pm.model,
        [t in terminals[i]], base_name="$(nw)_vm_$(i)",
        start = comp_start_value(ref(pm, nw, :bus, i), "vm_start", t, 1.0)
    ) for i in ids(pm, nw, :bus))

    if bounded
        for (i, bus) in ref(pm, nw, :bus)
            for (idx,t) in enumerate(terminals[i])
                set_lower_bound(vm[i][t], 0.0)

                if haskey(bus, "vmax")
                    set_upper_bound(vm[i][t], bus["vmax"][idx])
                end
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :bus, :vm, ids(pm, nw, :bus), vm)

end


"create variables for generators, delegate to PowerModels"
function variable_mc_generator_power(pm::_PM.AbstractPowerModel; kwargs...)
    variable_mc_generator_power_real(pm; kwargs...)
    variable_mc_generator_power_imaginary(pm; kwargs...)
end


""
function variable_mc_generator_power_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => gen["connections"] for (i,gen) in ref(pm, nw, :gen))
    pg = var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_pg_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "pg_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )

    if bounded
        for (i,gen) in ref(pm, nw, :gen)
            if haskey(gen, "pmin")
                for (idx,c) in enumerate(connections[i])
                    set_lower_bound(pg[i][c], gen["pmin"][idx])
                end
            end
            if haskey(gen, "pmax")
                for (idx,c) in enumerate(connections[i])
                    set_upper_bound(pg[i][c], gen["pmax"][idx])
                end
            end
        end
    end

    var(pm, nw)[:pg_bus] = Dict{Int, Any}()

    report && _IM.sol_component_value(pm, nw, :gen, :pg, ids(pm, nw, :gen), pg)
end


""
function variable_mc_generator_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => gen["connections"] for (i,gen) in ref(pm, nw, :gen))
    qg = var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_qg_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "qg_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )

    if bounded
        for (i,gen) in ref(pm, nw, :gen)
            if haskey(gen, "qmin")
                for (idx,c) in enumerate(connections[i])
                    set_lower_bound(qg[i][c], gen["qmin"][idx])
                end
            end
            if haskey(gen, "qmax")
                for (idx,c) in enumerate(connections[i])
                    set_upper_bound(qg[i][c], gen["qmax"][idx])
                end
            end
        end
    end

    var(pm, nw)[:qg_bus] = Dict{Int, Any}()

    report && _IM.sol_component_value(pm, nw, :gen, :qg, ids(pm, nw, :gen), qg)
end


"variable: `crg[j]` for `j` in `gen`"
function variable_mc_generator_current_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => gen["connections"] for (i,gen) in ref(pm, nw, :gen))
    crg = var(pm, nw)[:crg] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_crg_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "crg_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )
    if bounded
        for (i, g) in ref(pm, nw, :gen)
            cmax = _calc_gen_current_max(g, ref(pm, nw, :bus, g["gen_bus"]))
            for (idx,c) in enumerate(connections[i])
                set_lower_bound(crg[i][c], -cmax[idx])
                set_upper_bound(crg[i][c],  cmax[idx])
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :gen, :crg, ids(pm, nw, :gen), crg)
end

"variable: `cig[j]` for `j` in `gen`"
function variable_mc_generator_current_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => gen["connections"] for (i,gen) in ref(pm, nw, :gen))
    cig = var(pm, nw)[:cig] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_cig_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "cig_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )
    if bounded
        for (i, g) in ref(pm, nw, :gen)
            cmax = _calc_gen_current_max(g, ref(pm, nw, :bus, g["gen_bus"]))
            for (idx,c) in enumerate(connections[i])
                set_lower_bound(cig[i][c], -cmax[idx])
                set_upper_bound(cig[i][c],  cmax[idx])
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :gen, :cig, ids(pm, nw, :gen), cig)
end


function variable_mc_generator_power_on_off(pm::_PM.AbstractPowerModel; kwargs...)
    variable_mc_generator_power_real_on_off(pm; kwargs...)
    variable_mc_generator_power_imaginary_on_off(pm; kwargs...)
end


function variable_mc_generator_power_real_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => gen["connections"] for (i,gen) in ref(pm, nw, :gen))
    pg = var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_pg_$(i)",
        start = comp_start_value(ref(pm, nw, :gen, i), "pg_start", c, 0.0)
    ) for i in ids(pm, nw, :gen))

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            if haskey(gen, "pmin")
                for (idx, c) in enumerate(connections[i])
                    set_lower_bound(pg[i][c], gen["pmin"][idx])
                end
            end

            if haskey(gen, "pmax")
                for (idx, c) in enumerate(connections[i])
                    set_upper_bound(pg[i][c], gen["pmax"][idx])
                end
            end
        end
    end

    var(pm, nw)[:pg_bus] = Dict{Int, Any}()

    report && _IM.sol_component_value(pm, nw, :gen, :pg, ids(pm, nw, :gen), pg)
end


function variable_mc_generator_power_imaginary_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => gen["connections"] for (i,gen) in ref(pm, nw, :gen))
    qg = var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_qg_$(i)",
        start = comp_start_value(ref(pm, nw, :gen, i), "qg_start", c, 0.0)
    ) for i in ids(pm, nw, :gen))

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            if haskey(gen, "qmin")
                for (idx, c) in enumerate(connections[i])
                    set_lower_bound(qg[i][c], gen["qmin"][idx])
                end
            end

            if haskey(gen, "qmax")
                for (idx, c) in enumerate(connections[i])
                    set_upper_bound(qg[i][c], gen["qmax"][idx])
                end
            end
        end
    end

    var(pm, nw)[:qg_bus] = Dict{Int, Any}()

    report && _IM.sol_component_value(pm, nw, :gen, :qg, ids(pm, nw, :gen), qg)
end
