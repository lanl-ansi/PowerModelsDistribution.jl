
function comp_start_value(comp::Dict{String,<:Any}, key::String, conductor::Int, default)
    if haskey(comp, key)
        return comp[key][conductor]
    else
        return default
    end
end


function comp_start_value(comp::Dict{String,<:Any}, key::String, default)
    return _PM.comp_start_value(comp, key, default)
end


""
function variable_mc_bus_voltage_angle(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    va = var(pm, nw)[:va] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_va_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "va_start", 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    report && _IM.sol_component_value(pm, nw, :bus, :va, ids(pm, nw, :bus), va)
end

""
function variable_mc_bus_voltage_magnitude_only(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    vm = var(pm, nw)[:vm] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_vm_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "vm_start", c, 1.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in ref(pm, nw, :bus)
            if haskey(bus, "vmin")
                set_lower_bound.(vm[i], bus["vmin"])
            end
            if haskey(bus, "vmax")
                set_upper_bound.(vm[i], bus["vmax"])
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :bus, :vm, ids(pm, nw, :bus), vm)
end

""
function variable_mc_bus_voltage_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    vr = var(pm, nw)[:vr] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_vr_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "vr_start", c, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in ref(pm, nw, :bus)
            if haskey(bus, "vmax")
                set_lower_bound.(vr[i], -bus["vmax"])
                set_upper_bound.(vr[i],  bus["vmax"])
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :bus, :vr, ids(pm, nw, :bus), vr)
end

""
function variable_mc_bus_voltage_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    vi = var(pm, nw)[:vi] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_vi_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "vi_start", c, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in ref(pm, nw, :bus)
            if haskey(bus, "vmax")
                set_lower_bound.(vi[i], -bus["vmax"])
                set_upper_bound.(vi[i],  bus["vmax"])
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
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    p = var(pm, nw)[:p] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_p_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "p_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs)
            smax = _calc_branch_power_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
            set_upper_bound.(p[(l,i,j)],  smax)
            set_lower_bound.(p[(l,i,j)], -smax)
        end
    end

    for (l,branch) in ref(pm, nw, :branch)
        if haskey(branch, "pf_start")
            f_idx = (l, branch["f_bus"], branch["t_bus"])
            set_start_value(p[f_idx], branch["pf_start"])
        end
        if haskey(branch, "pt_start")
            t_idx = (l, branch["t_bus"], branch["f_bus"])
            set_start_value(p[t_idx], branch["pt_start"])
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :branch, :pf, :pt, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), p)
end

"variable: `q[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    q = var(pm, nw)[:q] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_q_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "q_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs)
            smax = _calc_branch_power_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
            set_upper_bound.(q[(l,i,j)],  smax)
            set_lower_bound.(q[(l,i,j)], -smax)
        end
    end

    for (l,branch) in ref(pm, nw, :branch)
        if haskey(branch, "qf_start")
            f_idx = (l, branch["f_bus"], branch["t_bus"])
            set_start_value(q[f_idx], branch["qf_start"])
        end
        if haskey(branch, "qt_start")
            t_idx = (l, branch["t_bus"], branch["f_bus"])
            set_start_value(q[t_idx], branch["qt_start"])
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :branch, :qf, :qt, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), q)
end


"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_current_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    branch = ref(pm, nw, :branch)
    bus = ref(pm, nw, :bus)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    cr = var(pm, nw)[:cr] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_cr_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "cr_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs)
            cmax = _calc_branch_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
            set_upper_bound.(cr[(l,i,j)],  cmax)
            set_lower_bound.(cr[(l,i,j)], -cmax)
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :branch, :cr_fr, :cr_to, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), cr)
end


"variable: `ci[l,i,j] ` for `(l,i,j)` in `arcs`"
function variable_mc_branch_current_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    branch = ref(pm, nw, :branch)
    bus = ref(pm, nw, :bus)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    ci = var(pm, nw)[:ci] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_ci_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "ci_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs)
            cmax = _calc_branch_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
            set_upper_bound.(ci[(l,i,j)],  cmax)
            set_lower_bound.(ci[(l,i,j)], -cmax)
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :branch, :ci_fr, :ci_to, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), ci)
end


"variable: `csr[l]` for `l` in `branch`"
function variable_mc_branch_current_series_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    branch = ref(pm, nw, :branch)
    bus = ref(pm, nw, :bus)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    csr = var(pm, nw)[:csr] = Dict(l => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_csr_$(l)",
            start = comp_start_value(ref(pm, nw, :branch, l), "csr_start", c, 0.0)
        ) for l in ids(pm, nw, :branch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_from)
            cmax = _calc_branch_series_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
            set_upper_bound.(csr[l],  cmax)
            set_lower_bound.(csr[l], -cmax)
        end
    end

    report && _IM.sol_component_value(pm, nw, :branch, :csr_fr, ids(pm, nw, :branch), csr)
end


"variable: `csi[l]` for `l` in `branch`"
function variable_mc_branch_current_series_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    branch = ref(pm, nw, :branch)
    bus = ref(pm, nw, :bus)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    csi = var(pm, nw)[:csi] = Dict(l => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_csi_$(l)",
            start = comp_start_value(ref(pm, nw, :branch, l), "csi_start", c, 0.0)
        ) for l in ids(pm, nw, :branch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_from)
            cmax = _calc_branch_series_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
            set_upper_bound.(csi[l],  cmax)
            set_lower_bound.(csi[l], -cmax)
        end
    end

    report && _IM.sol_component_value(pm, nw, :branch, :csi_fr, ids(pm, nw, :branch), csi)
end


"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_transformer_current_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    #trans = ref(pm, nw, :transformer)
    #bus = ref(pm, nw, :bus)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    cr = var(pm, nw)[:crt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_crt_$((l,i,j))",
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
    #trans = ref(pm, nw, :transformer)
    #bus = ref(pm, nw, :bus)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    ci = var(pm, nw)[:cit] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_cit_$((l,i,j))",
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
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    w = var(pm, nw)[:w] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_w_$(i)",
            lower_bound = 0.0,
            start = comp_start_value(ref(pm, nw, :bus, i), "w_start", 1.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for i in ids(pm, nw, :bus)
            bus = ref(pm, nw, :bus, i)
            vmax=bus["vmax"]
            vmin=bus["vmin"]
            for c in 1:ncnds
                set_upper_bound.((w[i])[c], max(vmin[c]^2, vmax[c]^2))
                if(vmin[c]>0)
                    set_lower_bound.((w[i])[c], vmin[c]^2)
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
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

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
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

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
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    ps = var(pm, nw)[:ps] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_ps_$(i)",
            start = comp_start_value(ref(pm, nw, :storage, i), "ps_start", c, 0.0)
        ) for i in ids(pm, nw, :storage)
    )

    if bounded
        for c in cnds
            flow_lb, flow_ub = _PM.ref_calc_storage_injection_bounds(ref(pm, nw, :storage), ref(pm, nw, :bus), c)

            for i in ids(pm, nw, :storage)
                if !isinf(flow_lb[i])
                    set_lower_bound(ps[i][c], flow_lb[i])
                end
                if !isinf(flow_ub[i])
                    set_upper_bound(ps[i][c], flow_ub[i])
                end
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :storage, :ps, ids(pm, nw, :storage), ps)
end


""
function variable_mc_storage_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    qs = var(pm, nw)[:qs] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_qs_$(i)",
            start = comp_start_value(ref(pm, nw, :storage, i), "qs_start", c, 0.0)
        ) for i in ids(pm, nw, :storage)
    )

    if bounded
        for c in cnds
            flow_lb, flow_ub = _PM.ref_calc_storage_injection_bounds(ref(pm, nw, :storage), ref(pm, nw, :bus), c)

            for i in ids(pm, nw, :storage)
                if !isinf(flow_lb[i])
                    set_lower_bound(qs[i][c], flow_lb[i])
                end
                if !isinf(flow_ub[i])
                    set_upper_bound(qs[i][c], flow_ub[i])
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
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    p_slack = var(pm, nw)[:p_slack] = Dict(i => JuMP.@variable(pm.model,
            [cnd in 1:ncnds], base_name="$(nw)_p_slack_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "p_slack_start", cnd, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    report && _IM.sol_component_value(pm, nw, :bus, :p_slack, ids(pm, nw, :bus), p_slack)
end


""
function variable_mc_slack_bus_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    q_slack = var(pm, nw)[:q_slack] = Dict(i => JuMP.@variable(pm.model,
            [cnd in 1:ncnds], base_name="$(nw)_q_slack_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "q_slack_start", cnd, 0.0)
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
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    pt = var(pm, nw)[:pt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_pt_$((l,i,j))",
        ) for (l,i,j) in ref(pm, nw, :arcs_trans)
    )

    if bounded
        for arc in ref(pm, nw, :arcs_from_trans)
            (t,i,j) = arc
            rate_a_fr, rate_a_to = _calc_transformer_power_ub_frto(ref(pm, nw, :transformer, t), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
            set_lower_bound.(pt[(t,i,j)], -rate_a_fr)
            set_upper_bound.(pt[(t,i,j)],  rate_a_fr)
            set_lower_bound.(pt[(t,j,i)], -rate_a_to)
            set_upper_bound.(pt[(t,j,i)],  rate_a_to)
        end
    end

    for (l,transformer) in ref(pm, nw, :transformer)
        if haskey(transformer, "pf_start")
            f_idx = (l, transformer["f_bus"], transformer["t_bus"])
            set_start_value(pt[f_idx], branch["pf_start"])
        end
        if haskey(transformer, "pt_start")
            t_idx = (l, transformer["t_bus"], transformer["f_bus"])
            set_start_value(pt[t_idx], transformer["pt_start"])
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :transformer, :pf, :pt, ref(pm, nw, :arcs_from_trans), ref(pm, nw, :arcs_to_trans), pt)
end


"Create variables for the reactive power flowing into all transformer windings."
function variable_mc_transformer_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    qt = var(pm, nw)[:qt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_qt_$((l,i,j))",
            start = 0.0
        ) for (l,i,j) in ref(pm, nw, :arcs_trans)
    )

    if bounded
        for arc in ref(pm, nw, :arcs_from_trans)
            (t,i,j) = arc
            rate_a_fr, rate_a_to = _calc_transformer_power_ub_frto(ref(pm, nw, :transformer, t), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))

            set_lower_bound.(qt[(t,i,j)], -rate_a_fr)
            set_upper_bound.(qt[(t,i,j)],  rate_a_fr)
            set_lower_bound.(qt[(t,j,i)], -rate_a_to)
            set_upper_bound.(qt[(t,j,i)],  rate_a_to)
        end
    end

    for (l,transformer) in ref(pm, nw, :transformer)
        if haskey(transformer, "qf_start")
            f_idx = (l, transformer["f_bus"], transformer["t_bus"])
            set_start_value(qt[f_idx], branch["qf_start"])
        end
        if haskey(transformer, "qt_start")
            t_idx = (l, transformer["t_bus"], transformer["f_bus"])
            set_start_value(qt[t_idx], transformer["qt_start"])
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :transformer, :qf, :qt, ref(pm, nw, :arcs_from_trans), ref(pm, nw, :arcs_to_trans), qt)
end


"Create tap variables."
function variable_mc_oltc_transformer_tap(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    # when extending to 4-wire, this should iterate only over the phase conductors
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    nph = 3
    p_oltc_ids = [id for (id,trans) in ref(pm, nw, :transformer) if !all(trans["tm_fix"])]
    tap = var(pm, nw)[:tap] = Dict(i => JuMP.@variable(pm.model,
        [p in 1:nph],
        base_name="$(nw)_tm_$(i)",
        start=ref(pm, nw, :transformer, i, "tm_set")[p]
    ) for i in p_oltc_ids)
    if bounded
        for tr_id in p_oltc_ids, p in 1:nph
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
function variable_mc_load_setpoint(pm::_PM.AbstractPowerModel; nw=pm.cnw, bounded::Bool=true, report::Bool=true)
    var(pm, nw)[:pd] = Dict{Int, Any}()
    var(pm, nw)[:qd] = Dict{Int, Any}()
    var(pm, nw)[:pd_bus] = Dict{Int, Any}()
    var(pm, nw)[:qd_bus] = Dict{Int, Any}()
end


"Create variables for demand status"
function variable_mc_load_indicator(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, relax::Bool=false, report::Bool=true)
    # this is not indexedon cnd; why used in start value?
    cnd = 1
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
    pd = var(pm, nw)[:pd] = Dict(i => var(pm, nw)[:z_demand][i].*ref(pm, nw, :load, i)["pd"]
     for i in ids(pm, nw, :load))
    qd = var(pm, nw)[:qd] = Dict(i => var(pm, nw)[:z_demand][i].*ref(pm, nw, :load, i)["qd"]
     for i in ids(pm, nw, :load))

    report && _IM.sol_component_value(pm, nw, :load, :status, ids(pm, nw, :load), z_demand)
    report && _IM.sol_component_value(pm, nw, :load, :pd, ids(pm, nw, :load), pd)
    report && _IM.sol_component_value(pm, nw, :load, :qd, ids(pm, nw, :load), qd)
end


"Create variables for shunt status"
function variable_mc_shunt_indicator(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, relax=false, report::Bool=true)
    # this is not indexedon cnd; why used in start value?
    cnd = 1
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
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    ps = var(pm, nw)[:ps] = Dict(i => JuMP.@variable(pm.model,
        [cnd in 1:ncnds], base_name="$(nw)_ps_$(i)",
        start = comp_start_value(ref(pm, nw, :storage, i), "ps_start", cnd, 0.0)
    ) for i in ids(pm, nw, :storage))

    if bounded
        for cnd in 1:ncnds
            inj_lb, inj_ub = _PM.ref_calc_storage_injection_bounds(ref(pm, nw, :storage), ref(pm, nw, :bus), cnd)

            for (i, strg) in ref(pm, nw, :storage)
                set_lower_bound.(ps[i][cnd], inj_lb[i])
                set_upper_bound.(ps[i][cnd], inj_ub[i])
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :storage, :ps, ids(pm, nw, :storage), ps)
end


"Create variables for `reactive` storage injection"
function variable_mc_storage_power_imaginary_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    qs = var(pm, nw)[:qs] = Dict(i => JuMP.@variable(pm.model,
        [cnd in 1:ncnds], base_name="$(nw)_qs_$(i)",
        start = comp_start_value(ref(pm, nw, :storage, i), "qs_start", cnd, 0.0)
    ) for i in ids(pm, nw, :storage))

    if bounded
        for (i, strg) in ref(pm, nw, :storage)
            if haskey(strg, "qmin")
                set_lower_bound.(qs[i], strg["qmin"])
            end

            if haskey(strg, "qmax")
                set_upper_bound.(qs[i], strg["qmax"])
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :storage, :qs, ids(pm, nw, :storage), qs)
end


"voltage variable magnitude squared (relaxed form)"
function variable_mc_bus_voltage_magnitude_sqr_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    w = var(pm, nw)[:w] = Dict(i => JuMP.@variable(pm.model,
        [c in 1:ncnds], base_name="$(nw)_w_$(i)",
        start = comp_start_value(ref(pm, nw, :bus, i), "w_start", c, 1.001)
    ) for i in ids(pm, nw, :bus))

    if bounded
        for (i, bus) in ref(pm, nw, :bus)
            set_lower_bound.(w[i], 0.0)

            if haskey(bus, "vmax")
                set_upper_bound.(w[i], bus["vmax"].^2)
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :bus, :w, ids(pm, nw, :bus), w)
end


"on/off voltage magnitude variable"
function variable_mc_bus_voltage_magnitude_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    vm = var(pm, nw)[:vm] = Dict(i => JuMP.@variable(pm.model,
        [c in 1:ncnds], base_name="$(nw)_vm_$(i)",
        start = comp_start_value(ref(pm, nw, :bus, i), "vm_start", c, 1.0)
    ) for i in ids(pm, nw, :bus))

    if bounded
        for (i, bus) in ref(pm, nw, :bus)
            set_lower_bound.(vm[i], 0.0)

            if haskey(bus, "vmax")
                set_upper_bound.(vm[i], bus["vmax"])
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :bus, :vm, ids(pm, nw, :bus), vm)

end


"create variables for generators, delegate to PowerModels"
function variable_mc_gen_power_setpoint(pm::_PM.AbstractPowerModel; kwargs...)
    variable_mc_gen_power_setpoint_real(pm; kwargs...)
    variable_mc_gen_power_setpoint_imaginary(pm; kwargs...)
end


function variable_mc_gen_power_setpoint_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    pg = var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_pg_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "pg_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )

    if bounded
        for (i,gen) in ref(pm, nw, :gen)
            if haskey(gen, "pmin")
                set_lower_bound.(pg[i], gen["pmin"])
            end
            if haskey(gen, "pmax")
                set_upper_bound.(pg[i], gen["pmax"])
            end
        end
    end

    var(pm, nw)[:pg_bus] = Dict{Int, Any}()

    report && _IM.sol_component_value(pm, nw, :gen, :pg, ids(pm, nw, :gen), pg)
end


function variable_mc_gen_power_setpoint_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    qg = var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_qg_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "qg_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )

    if bounded
        for (i,gen) in ref(pm, nw, :gen)
            if haskey(gen, "qmin")
                set_lower_bound.(qg[i], gen["qmin"])
            end
            if haskey(gen, "qmax")
                set_upper_bound.(qg[i], gen["qmax"])
            end
        end
    end

    var(pm, nw)[:qg_bus] = Dict{Int, Any}()

    report && _IM.sol_component_value(pm, nw, :gen, :qg, ids(pm, nw, :gen), qg)
end


"variable: `crg[j]` for `j` in `gen`"
function variable_mc_gen_current_setpoint_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    gen = ref(pm, nw, :gen)
    bus = ref(pm, nw, :bus)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    crg = var(pm, nw)[:crg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_crg_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "crg_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )
    if bounded
        for (i, g) in gen
            cmax = _calc_gen_current_max(g, ref(pm, nw, :bus, g["gen_bus"]))
            set_lower_bound.(crg[i], -cmax)
            set_upper_bound.(crg[i],  cmax)
        end
    end

    report && _IM.sol_component_value(pm, nw, :gen, :crg, ids(pm, nw, :gen), crg)
end

"variable: `cig[j]` for `j` in `gen`"
function variable_mc_gen_current_setpoint_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    gen = ref(pm, nw, :gen)
    bus = ref(pm, nw, :bus)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    cig = var(pm, nw)[:cig] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_cig_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "cig_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )
    if bounded
        for (i, g) in gen
            cmax = _calc_gen_current_max(g, ref(pm, nw, :bus, g["gen_bus"]))
            set_lower_bound.(cig[i], -cmax)
            set_upper_bound.(cig[i],  cmax)
        end
    end

    report && _IM.sol_component_value(pm, nw, :gen, :cig, ids(pm, nw, :gen), cig)
end


function variable_mc_gen_power_setpoint_on_off(pm::_PM.AbstractPowerModel; kwargs...)
    variable_mc_gen_power_setpoint_real_on_off(pm; kwargs...)
    variable_mc_gen_power_setpoint_imaginary_on_off(pm; kwargs...)
end


function variable_mc_gen_power_setpoint_real_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(conductor_ids(pm, nw))

    pg = var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
        [cnd in 1:ncnds], base_name="$(nw)_pg_$(i)",
        start = comp_start_value(ref(pm, nw, :gen, i), "pg_start", cnd, 0.0)
    ) for i in ids(pm, nw, :gen))

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            if haskey(gen, "pmin")
                set_lower_bound.(pg[i], gen["pmin"])
            end

            if haskey(gen, "pmax")
                set_upper_bound.(pg[i], gen["pmax"])
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :gen, :pg, ids(pm, nw, :gen), pg)
end


function variable_mc_gen_power_setpoint_imaginary_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    qg = var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
        [cnd in 1:ncnds], base_name="$(nw)_qg_$(i)",
        start = comp_start_value(ref(pm, nw, :gen, i), "qg_start", cnd, 0.0)
    ) for i in ids(pm, nw, :gen))

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            if haskey(gen, "qmin")
                set_lower_bound.(qg[i], gen["qmin"])
            end

            if haskey(gen, "qmax")
                set_upper_bound.(qg[i], gen["qmax"])
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :gen, :qg, ids(pm, nw, :gen), qg)
end
