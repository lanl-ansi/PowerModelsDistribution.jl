# voltage variables

""
function variable_mc_bus_voltage_angle(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    va_start_defaults = Dict(i => [0.0, -120.0, 120.0, fill(0.0, length(terms))...][terms] for (i, terms) in terminals)
    va = var(pm, nw)[:va] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_va_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), ["va_start", "va"], t, va_start_defaults[i][findfirst(isequal(t), terminals[i])]),
        ) for i in ids(pm, nw, :bus)
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :va, ids(pm, nw, :bus), va)
end


""
function variable_mc_bus_voltage_magnitude_only(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    vm = var(pm, nw)[:vm] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_vm_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), ["vm_start", "vm", "vmin"], t, 1.0)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :vm, ids(pm, nw, :bus), vm)
end


""
function variable_mc_bus_voltage_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i, bus) in ref(pm, nw, :bus))

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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :vr, ids(pm, nw, :bus), vr)
end


""
function variable_mc_bus_voltage_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :vi, ids(pm, nw, :bus), vi)
end


"variable: `w[i] >= 0` for `i` in `buses"
function variable_mc_bus_voltage_magnitude_sqr(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    w = var(pm, nw)[:w] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_w_$(i)",
            lower_bound = 0.0,
            start = comp_start_value(ref(pm, nw, :bus, i), "w_start", t, comp_start_value(ref(pm, nw, :bus, i), ["vm_start", "vm", "vmin"], t, 1.0)^2)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :w, ids(pm, nw, :bus), w)
end

## voltage on/off variables

"on/off voltage magnitude variable"
function variable_mc_bus_voltage_magnitude_on_off(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    vm = var(pm, nw)[:vm] = Dict(i => JuMP.@variable(pm.model,
        [t in terminals[i]], base_name="$(nw)_vm_$(i)",
        start = comp_start_value(ref(pm, nw, :bus, i), ["vm_start", "vm", "vmin"], t, 1.0)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :vm, ids(pm, nw, :bus), vm)

end


"voltage variable magnitude squared (relaxed form)"
function variable_mc_bus_voltage_magnitude_sqr_on_off(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    w = var(pm, nw)[:w] = Dict(i => JuMP.@variable(pm.model,
        [t in terminals[i]], base_name="$(nw)_w_$(i)",
        start = comp_start_value(ref(pm, nw, :bus, i), "w_start", t, comp_start_value(ref(pm, nw, :bus, i), ["vm_start", "vm", "vmin"], t, 1.0)^2)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :w, ids(pm, nw, :bus), w)
end


# branch variables

## branch power variables

"branch flow variables, delegated back to PowerModels"
function variable_mc_branch_power(pm::AbstractUnbalancedPowerModel; kwargs...)
    variable_mc_branch_power_real(pm; kwargs...)
    variable_mc_branch_power_imaginary(pm; kwargs...)
end


"variable: `p[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_power_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)
    p = var(pm, nw)[:p] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_p_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "p_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_branch)
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

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch, :pf, :pt, ref(pm, nw, :arcs_branch_from), ref(pm, nw, :arcs_branch_to), p)
end


"variable: `q[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_power_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)
    q = var(pm, nw)[:q] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_q_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "q_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_branch)
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

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch, :qf, :qt, ref(pm, nw, :arcs_branch_from), ref(pm, nw, :arcs_branch_to), q)
end

## branch current variables

"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_current_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)
    cr = var(pm, nw)[:cr] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_cr_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "cr_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_branch)
            cmax = _calc_branch_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
            for (idx,c) in enumerate(connections[(l,i,j)])
                set_upper_bound(cr[(l,i,j)][c],  cmax[idx])
                set_lower_bound(cr[(l,i,j)][c], -cmax[idx])
            end
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch, :cr_fr, :cr_to, ref(pm, nw, :arcs_branch_from), ref(pm, nw, :arcs_branch_to), cr)
end


"variable: `ci[l,i,j] ` for `(l,i,j)` in `arcs`"
function variable_mc_branch_current_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)
    ci = var(pm, nw)[:ci] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_ci_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "ci_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_branch)
            cmax = _calc_branch_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
            for (idx,c) in enumerate(connections[(l,i,j)])
                set_upper_bound(ci[(l,i,j)][c],  cmax[idx])
                set_lower_bound(ci[(l,i,j)][c], -cmax[idx])
            end
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch, :ci_fr, :ci_to, ref(pm, nw, :arcs_branch_from), ref(pm, nw, :arcs_branch_to), ci)
end


"variable: `csr[l]` for `l` in `branch`"
function variable_mc_branch_current_series_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)
    csr = var(pm, nw)[:csr] = Dict(l => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_csr_$(l)",
            start = comp_start_value(ref(pm, nw, :branch, l), "csr_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_branch_from)
            cmax = _calc_branch_series_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
            for (idx,c) in enumerate(connections[(l,i,j)])
                set_upper_bound(csr[l][c],  cmax[idx])
                set_lower_bound(csr[l][c], -cmax[idx])
            end
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch, :csr_fr, ids(pm, nw, :branch), csr)
end


"variable: `csi[l]` for `l` in `branch`"
function variable_mc_branch_current_series_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)
    csi = var(pm, nw)[:csi] = Dict(l => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_csi_$(l)",
            start = comp_start_value(ref(pm, nw, :branch, l), "csi_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_branch_from)
            cmax = _calc_branch_series_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
            for (idx,c) in enumerate(connections[(l,i,j)])
                set_upper_bound(csi[l][c],  cmax[idx])
                set_lower_bound(csi[l][c], -cmax[idx])
            end
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch, :csi_fr, ids(pm, nw, :branch), csi)
end


# transformer variables

## transformer power variables

"Creates variables for both `active` and `reactive` power flow at each transformer."
function variable_mc_transformer_power(pm::AbstractUnbalancedPowerModel; kwargs...)
    variable_mc_transformer_power_real(pm; kwargs...)
    variable_mc_transformer_power_imaginary(pm; kwargs...)
end


"Create variables for the active power flowing into all transformer windings."
function variable_mc_transformer_power_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_transformer) for ((l,i,j), connections) in entry)
    pt = var(pm, nw)[:pt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]],
            base_name="$(nw)_pt_$((l,i,j))",
            start = 0.0,
        ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
    )

    if bounded
        for arc in ref(pm, nw, :arcs_transformer_from)
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

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :pf, :pt, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), pt)
end


"Create variables for the reactive power flowing into all transformer windings."
function variable_mc_transformer_power_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_transformer) for ((l,i,j), connections) in entry)
    qt = var(pm, nw)[:qt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_qt_$((l,i,j))",
            start = 0.0
        ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
    )

    if bounded
        for arc in ref(pm, nw, :arcs_transformer_from)
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
            for (idx, fc) in enumerate(connections[f_idx])
                JuMP.set_start_value(qt[f_idx][fc], transformer["qf_start"][idx])
            end
        end
        if haskey(transformer, "qt_start")
            t_idx = (l, transformer["t_bus"], transformer["f_bus"])
            for (idx, tc) in enumerate(connections[t_idx])
                JuMP.set_start_value(qt[t_idx][tc], transformer["qt_start"][idx])
            end
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :qf, :qt, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), qt)
end

## transformer current variables

"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_transformer_current_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_transformer) for ((l,i,j), connections) in entry)
    cr = var(pm, nw)[:crt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_crt_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :transformer, l), "cr_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_transformer_from)
            trans = ref(pm, nw, :transformer, l)
            f_bus = ref(pm, nw, :bus, i)
            t_bus = ref(pm, nw, :bus, j)
            cmax_fr, cmax_to = _calc_transformer_current_max_frto(trans, f_bus, t_bus)
            for (idx, (fc,tc)) in enumerate(zip(trans["f_connections"], trans["t_connections"]))
                set_lower_bound(cr[(l,i,j)][fc], -cmax_fr[idx])
                set_upper_bound(cr[(l,i,j)][fc],  cmax_fr[idx])
                set_lower_bound(cr[(l,j,i)][tc], -cmax_to[idx])
                set_upper_bound(cr[(l,j,i)][tc],  cmax_to[idx])
            end
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :cr_fr, :cr_to, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), cr)
end


"variable: `ci[l,i,j] ` for `(l,i,j)` in `arcs`"
function variable_mc_transformer_current_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_transformer) for ((l,i,j), connections) in entry)
    ci = var(pm, nw)[:cit] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_cit_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :transformer, l), "ci_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_transformer_from)
            trans = ref(pm, nw, :transformer, l)
            f_bus = ref(pm, nw, :bus, i)
            t_bus = ref(pm, nw, :bus, j)
            cmax_fr, cmax_to = _calc_transformer_current_max_frto(trans, f_bus, t_bus)
            for (idx, (fc,tc)) in enumerate(zip(trans["f_connections"], trans["t_connections"]))
                set_lower_bound(ci[(l,i,j)][fc], -cmax_fr[idx])
                set_upper_bound(ci[(l,i,j)][fc],  cmax_fr[idx])
                set_lower_bound(ci[(l,j,i)][tc], -cmax_to[idx])
                set_upper_bound(ci[(l,j,i)][tc],  cmax_to[idx])
            end
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :ci_fr, :ci_to, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), ci)
end

## transformer tap variables

"Create tap variables."
function variable_mc_oltc_transformer_tap(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # when extending to 4-wire, this should iterate only over the phase conductors
    p_oltc_ids = [id for (id,trans) in ref(pm, nw, :transformer) if !all(trans["tm_fix"])]
    tap = var(pm, nw)[:tap] = Dict(i => JuMP.@variable(pm.model,
        [p in 1:length(ref(pm,nw,:transformer,i,"f_connections"))],
        base_name="$(nw)_tm_$(i)",
        start=comp_start_value(ref(pm, nw, :transformer, i), "tm_set", ref(pm, nw, :transformer, i, "f_connections")[p], 1.0),
    ) for i in p_oltc_ids)

    if bounded
        for tr_id in p_oltc_ids, p in 1:length(ref(pm,nw,:transformer,tr_id,"f_connections"))
            set_lower_bound(var(pm, nw)[:tap][tr_id][p], ref(pm, nw, :transformer, tr_id, "tm_lb")[p])
            set_upper_bound(var(pm, nw)[:tap][tr_id][p], ref(pm, nw, :transformer, tr_id, "tm_ub")[p])
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :transformer, :tap, p_oltc_ids, tap)
end


# switch variables

## switch power variables

""
function variable_mc_switch_power(pm::AbstractUnbalancedPowerModel; kwargs...)
    variable_mc_switch_power_real(pm; kwargs...)
    variable_mc_switch_power_imaginary(pm; kwargs...)
end


""
function variable_mc_switch_power_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_switch) for ((l,i,j), connections) in entry)
    psw = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_psw_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :switch, l), "psw_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_switch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_switch)
            smax = _calc_branch_power_max(ref(pm, nw, :switch, l), ref(pm, nw, :bus, i))
            for (idx, c) in enumerate(connections[(l,i,j)])
                set_upper_bound(psw[(l,i,j)][c],  smax[idx])
                set_lower_bound(psw[(l,i,j)][c], -smax[idx])
            end
        end
    end

    # this explicit type erasure is necessary
    psw_expr = Dict{Any,Any}( (l,i,j) => psw[(l,i,j)] for (l,i,j) in ref(pm, nw, :arcs_switch_from) )
    psw_expr = merge(psw_expr, Dict( (l,j,i) => -1.0.*psw[(l,i,j)] for (l,i,j) in ref(pm, nw, :arcs_switch_from)))

    # This is needed to get around error: "unexpected affine expression in nlconstraint"
    psw_auxes = Dict{Any,Any}(
        (l,i,j) => JuMP.@variable(
            pm.model, [c in connections[(l,i,j)]],
            base_name="$(nw)_psw_aux_$((l,i,j))"
        ) for (l,i,j) in ref(pm, nw, :arcs_switch)
    )
    for ((l,i,j), psw_aux) in psw_auxes
        for (idx, c) in enumerate(connections[(l,i,j)])
            JuMP.@constraint(pm.model, psw_expr[(l,i,j)][c] == psw_aux[c])
        end
    end

    var(pm, nw)[:psw] = psw_auxes

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :switch, :pf, :pt, ref(pm, nw, :arcs_switch_from), ref(pm, nw, :arcs_switch_to), psw_expr)
end


""
function variable_mc_switch_power_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_switch) for ((l,i,j), connections) in entry)
    qsw = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_qsw_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :switch, l), "qsw_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_switch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_switch)
            smax = _calc_branch_power_max(ref(pm, nw, :switch, l), ref(pm, nw, :bus, i))
            for (idx, c) in enumerate(connections[(l,i,j)])
                set_upper_bound(qsw[(l,i,j)][c],  smax[idx])
                set_lower_bound(qsw[(l,i,j)][c], -smax[idx])
            end
        end
    end

    # this explicit type erasure is necessary
    qsw_expr = Dict{Any,Any}( (l,i,j) => qsw[(l,i,j)] for (l,i,j) in ref(pm, nw, :arcs_switch_from) )
    qsw_expr = merge(qsw_expr, Dict( (l,j,i) => -1.0*qsw[(l,i,j)] for (l,i,j) in ref(pm, nw, :arcs_switch_from)))

    # This is needed to get around error: "unexpected affine expression in nlconstraint"
    qsw_auxes = Dict{Any,Any}(
        (l,i,j) => JuMP.@variable(
            pm.model, [c in connections[(l,i,j)]],
            base_name="$(nw)_qsw_aux_$((l,i,j))"
        ) for (l,i,j) in ref(pm, nw, :arcs_switch)
    )
    for ((l,i,j), qsw_aux) in qsw_auxes
        for (idx, c) in enumerate(connections[(l,i,j)])
            JuMP.@constraint(pm.model, qsw_expr[(l,i,j)][c] == qsw_aux[c])
        end
    end

    var(pm, nw)[:qsw] = qsw_auxes

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :switch, :qf, :qt, ref(pm, nw, :arcs_switch_from), ref(pm, nw, :arcs_switch_to), qsw_expr)
end

## switch current variables

""
function variable_mc_switch_current(pm::AbstractUnbalancedPowerModel; kwargs...)
    variable_mc_switch_current_real(pm; kwargs...)
    variable_mc_switch_current_imaginary(pm; kwargs...)
end


""
function variable_mc_switch_current_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_switch) for ((l,i,j), connections) in entry)
    crsw = var(pm, nw)[:crsw] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_crsw_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :switch, l), "crsw_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_switch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_switch)
            cmax = _calc_branch_current_max(ref(pm, nw, :switch, l), ref(pm, nw, :bus, i))
            for (idx,c) in enumerate(connections[(l,i,j)])
                set_upper_bound(crsw[(l,i,j)][c],  cmax[idx])
                set_lower_bound(crsw[(l,i,j)][c], -cmax[idx])
            end
        end
    end

    # this explicit type erasure is necessary
    crsw_expr = Dict{Any,Any}( (l,i,j) => crsw[(l,i,j)] for (l,i,j) in ref(pm, nw, :arcs_switch_from) )
    crsw_expr = merge(crsw_expr, Dict( (l,j,i) => -crsw[(l,i,j)] for (l,i,j) in ref(pm, nw, :arcs_switch_from)))

    # This is needed to get around error: "unexpected affine expression in nlconstraint"
    crsw_auxes = Dict{Any,Any}(
        (l,i,j) => JuMP.@variable(
            pm.model, [c in connections[(l,i,j)]],
            base_name="$(nw)_crsw_aux_$((l,i,j))"
        ) for (l,i,j) in ref(pm, nw, :arcs_switch)
    )
    for ((l,i,j), crsw_aux) in crsw_auxes
        for (idx, c) in enumerate(connections[(l,i,j)])
            JuMP.@constraint(pm.model, crsw_expr[(l,i,j)][c] == crsw_aux[c])
        end
    end

    var(pm, nw)[:crsw] = crsw_auxes

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :switch, :cr_fr, :cr_to, ref(pm, nw, :arcs_switch_from), ref(pm, nw, :arcs_switch_to), crsw_expr)
end


""
function variable_mc_switch_current_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_switch) for ((l,i,j), connections) in entry)
    cisw = var(pm, nw)[:cisw] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_cisw_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :switch, l), "cisw_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_switch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_switch)
            cmax = _calc_branch_current_max(ref(pm, nw, :switch, l), ref(pm, nw, :bus, i))
            for (idx,c) in enumerate(connections[(l,i,j)])
                set_upper_bound(cisw[(l,i,j)][c],  cmax[idx])
                set_lower_bound(cisw[(l,i,j)][c], -cmax[idx])
            end
        end
    end

    # this explicit type erasure is necessary
    cisw_expr = Dict{Any,Any}( (l,i,j) => cisw[(l,i,j)] for (l,i,j) in ref(pm, nw, :arcs_switch_from) )
    cisw_expr = merge(cisw_expr, Dict( (l,j,i) => -cisw[(l,i,j)] for (l,i,j) in ref(pm, nw, :arcs_switch_from)))

    # This is needed to get around error: "unexpected affine expression in nlconstraint"
    cisw_auxes = Dict{Any,Any}(
        (l,i,j) => JuMP.@variable(
            pm.model, [c in connections[(l,i,j)]],
            base_name="$(nw)_cisw_aux_$((l,i,j))"
        ) for (l,i,j) in ref(pm, nw, :arcs_switch)
    )
    for ((l,i,j), cisw_aux) in cisw_auxes
        for (idx, c) in enumerate(connections[(l,i,j)])
            JuMP.@constraint(pm.model, cisw_expr[(l,i,j)][c] == cisw_aux[c])
        end
    end

    var(pm, nw)[:cisw] = cisw_auxes

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :switch, :ci_fr, :ci_to, ref(pm, nw, :arcs_switch_from), ref(pm, nw, :arcs_switch_to), cisw_expr)
end

# switch state variables

"switch state (open/close) variables"
function variable_mc_switch_state(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, report::Bool=true, relax::Bool=false)
    if relax
        state = var(pm, nw)[:switch_state] = JuMP.@variable(
            pm.model,
            [l in ids(pm, nw, :switch_dispatchable)],
            base_name="$(nw)_switch_state_$(l)",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(ref(pm, nw, :switch, l), ["state_start", "state"], 0)
        )
    else
        state = var(pm, nw)[:switch_state] = JuMP.@variable(
            pm.model,
            [l in ids(pm, nw, :switch_dispatchable)],
            base_name="$(nw)_switch_state_$(l)",
            binary = true,
            start = comp_start_value(ref(pm, nw, :switch, l), ["state_start", "state"], 0)
        )
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :switch, :state, ids(pm, nw, :switch_dispatchable), state)
end


"""
    variable_mc_capacitor_switch_state(pm::AbstractUnbalancedPowerModel, relax::Bool; nw::Int=nw_id_default, report::Bool=true)

Capacitor (with capcontrol) switch state (open/close) variables
"""
function variable_mc_capacitor_switch_state(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=true, report::Bool=true)
    cap_switch_ids = [id for (id,cap) in ref(pm, nw, :shunt) if haskey(cap,"controls")]
    if relax
        cap_state = var(pm, nw)[:capacitor_state] = Dict(i => JuMP.@variable(pm.model,
        [p in ref(pm, nw, :shunt, i, "connections")],
        base_name="$(nw)_cap_sw_$(i)",
        lower_bound = 0,
        upper_bound = 1,
        start = 0.5
        ) for i in cap_switch_ids)
    else
        cap_state = var(pm, nw)[:capacitor_state] = Dict(i => JuMP.@variable(pm.model,
        [p in ref(pm, nw, :shunt, i, "connections")],
        base_name="$(nw)_cap_sw_$(i)",
        binary = true,
        start = 1
        ) for i in cap_switch_ids)
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :shunt, :cap_state, cap_switch_ids, cap_state)
end


"""
    variable_mc_capacitor_reactive_power(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default)

Capacitor (with capcontrol) relaxed power variables for AbstractLPUBFModel (using McCormick envelopes)
"""
function variable_mc_capacitor_reactive_power(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default)
    cap_switch_ids = [id for (id,cap) in ref(pm, nw, :shunt) if haskey(cap,"controls")]
    cap_reactive_power = var(pm, nw)[:capacitor_reactive_power] = Dict(i => JuMP.@variable(pm.model,
    [p in ref(pm, nw, :shunt, i, "connections")],
    base_name="$(nw)_cap_cur_$(i)"
    ) for i in cap_switch_ids)
end


# generator variables

## generator power variables

"create variables for generators, delegate to PowerModels"
function variable_mc_generator_power(pm::AbstractUnbalancedPowerModel; kwargs...)
    variable_mc_generator_power_real(pm; kwargs...)
    variable_mc_generator_power_imaginary(pm; kwargs...)
end


""
function variable_mc_generator_power_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => gen["connections"] for (i,gen) in ref(pm, nw, :gen))
    pg = var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_pg_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), ["pg_start", "pg", "pmin"], c, 0.0)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :pg, ids(pm, nw, :gen), pg)
end


""
function variable_mc_generator_power_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => gen["connections"] for (i,gen) in ref(pm, nw, :gen))
    qg = var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_qg_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), ["qg_start", "qg", "qmin"], c, 0.0)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :qg, ids(pm, nw, :gen), qg)
end

### generator power on/off variables

""
function variable_mc_generator_power_on_off(pm::AbstractUnbalancedPowerModel; kwargs...)
    variable_mc_generator_power_real_on_off(pm; kwargs...)
    variable_mc_generator_power_imaginary_on_off(pm; kwargs...)
end


""
function variable_mc_generator_power_real_on_off(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => gen["connections"] for (i,gen) in ref(pm, nw, :gen))
    pg = var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_pg_$(i)",
        start = comp_start_value(ref(pm, nw, :gen, i), ["pg_start", "pg", "pmin"], c, 0.0)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :pg, ids(pm, nw, :gen), pg)
end


""
function variable_mc_generator_power_imaginary_on_off(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => gen["connections"] for (i,gen) in ref(pm, nw, :gen))
    qg = var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_qg_$(i)",
        start = comp_start_value(ref(pm, nw, :gen, i), ["qg_start", "qg", "qmin"], c, 0.0)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :qg, ids(pm, nw, :gen), qg)
end

## generator current variables

"variable: `crg[j]` for `j` in `gen`"
function variable_mc_generator_current_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :crg, ids(pm, nw, :gen), crg)
end


"variable: `cig[j]` for `j` in `gen`"
function variable_mc_generator_current_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :cig, ids(pm, nw, :gen), cig)
end


# storage variables

## storage state variables
""
function variable_storage_energy(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    se = var(pm, nw)[:se] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :storage)], base_name="$(nw)_se",
        start = comp_start_value(ref(pm, nw, :storage, i), "se_start", 1)
    )

    if bounded
        for (i, storage) in ref(pm, nw, :storage)
            set_lower_bound(se[i], 0)
            set_upper_bound(se[i], storage["energy_rating"])
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :se, ids(pm, nw, :storage), se)
end


""
function variable_storage_charge(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    sc = var(pm, nw)[:sc] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :storage)], base_name="$(nw)_sc",
        start = comp_start_value(ref(pm, nw, :storage, i), "sc_start", 1)
    )

    if bounded
        for (i, storage) in ref(pm, nw, :storage)
            set_lower_bound(sc[i], 0)
            set_upper_bound(sc[i], storage["charge_rating"])
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :sc, ids(pm, nw, :storage), sc)
end


""
function variable_storage_discharge(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    sd = var(pm, nw)[:sd] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :storage)], base_name="$(nw)_sd",
        start = comp_start_value(ref(pm, nw, :storage, i), "sd_start", 1)
    )

    if bounded
        for (i, storage) in ref(pm, nw, :storage)
            set_lower_bound(sd[i], 0)
            set_upper_bound(sd[i], storage["discharge_rating"])
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :sd, ids(pm, nw, :storage), sd)
end


""
function variable_storage_complementary_indicator(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if !relax
        sc_on = var(pm, nw)[:sc_on] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :storage)], base_name="$(nw)_sc_on",
            binary = true,
            start = comp_start_value(ref(pm, nw, :storage, i), "sc_on_start", 0)
        )
        sd_on = var(pm, nw)[:sd_on] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :storage)], base_name="$(nw)_sd_on",
            binary = true,
            start = comp_start_value(ref(pm, nw, :storage, i), "sd_on_start", 0)
        )
    else
        sc_on = var(pm, nw)[:sc_on] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :storage)], base_name="$(nw)_sc_on",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(ref(pm, nw, :storage, i), "sc_on_start", 0)
        )
        sd_on = var(pm, nw)[:sd_on] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :storage)], base_name="$(nw)_sd_on",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(ref(pm, nw, :storage, i), "sd_on_start", 0)
        )
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :sc_on, ids(pm, nw, :storage), sc_on)
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :sd_on, ids(pm, nw, :storage), sd_on)
end

## storage power variables

"variables for modeling storage units, includes grid injection and internal variables"
function variable_mc_storage_power(pm::AbstractUnbalancedPowerModel; kwargs...)
    variable_mc_storage_power_real(pm; kwargs...)
    variable_mc_storage_power_imaginary(pm; kwargs...)
    variable_mc_storage_power_control_imaginary(pm; kwargs...)
    variable_mc_storage_current(pm; kwargs...)
    variable_storage_energy(pm; kwargs...)
    variable_storage_charge(pm; kwargs...)
    variable_storage_discharge(pm; kwargs...)
end


""
function variable_mc_storage_power_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i,strg) in ref(pm, nw, :storage))
    ps = var(pm, nw)[:ps] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_ps_$(i)",
            start = comp_start_value(ref(pm, nw, :storage, i), ["ps_start", "ps"], c, 0.0)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :ps, ids(pm, nw, :storage), ps)
end


""
function variable_mc_storage_power_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i,strg) in ref(pm, nw, :storage))
    qs = var(pm, nw)[:qs] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_qs_$(i)",
            start = comp_start_value(ref(pm, nw, :storage, i), ["qs_start", "qs"], c, 0.0)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :qs, ids(pm, nw, :storage), qs)
end


"""
a reactive power slack variable that enables the storage device to inject or
consume reactive power at its connecting bus, subject to the injection limits
of the device.
"""
function variable_mc_storage_power_control_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    qsc = var(pm, nw)[:qsc] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :storage)], base_name="$(nw)_qsc_$(i)",
        start = comp_start_value(ref(pm, nw, :storage, i), "qsc_start")
    )

    if bounded
        inj_lb, inj_ub = ref_calc_storage_injection_bounds(ref(pm, nw, :storage), ref(pm, nw, :bus))
        for (i,storage) in ref(pm, nw, :storage)
            if !isinf(sum(inj_lb[i])) || haskey(storage, "qmin")
                set_lower_bound(qsc[i], max(sum(inj_lb[i]), sum(get(storage, "qmin", -Inf))))
            end
            if !isinf(sum(inj_ub[i])) || haskey(storage, "qmax")
                set_upper_bound(qsc[i], min(sum(inj_ub[i]), sum(get(storage, "qmax",  Inf))))
            end
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :qsc, ids(pm, nw, :storage), qsc)
end


""
function variable_mc_storage_power_mi(pm::AbstractUnbalancedPowerModel; relax::Bool=false, kwargs...)
    variable_mc_storage_power_real(pm; kwargs...)
    variable_mc_storage_power_imaginary(pm; kwargs...)
    variable_mc_storage_power_control_imaginary(pm; kwargs...)
    variable_mc_storage_current(pm; kwargs...)
    variable_mc_storage_indicator(pm; relax=relax, kwargs...)
    variable_storage_energy(pm; kwargs...)
    variable_storage_charge(pm; kwargs...)
    variable_storage_discharge(pm; kwargs...)
    variable_storage_complementary_indicator(pm; relax=relax, kwargs...)
end

### storage power on/off variables


"Create variables for `active` and `reactive` storage injection"
function variable_mc_storage_power_on_off(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, kwargs...)
    variable_mc_storage_power_real_on_off(pm; nw=nw, kwargs...)
    variable_mc_storage_power_imaginary_on_off(pm; nw=nw, kwargs...)
end


"Create variables for `active` storage injection"
function variable_mc_storage_power_real_on_off(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i,strg) in ref(pm, nw, :storage))
    ps = var(pm, nw)[:ps] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_ps_$(i)",
        start = comp_start_value(ref(pm, nw, :storage, i), ["ps_start", "ps"], c, 0.0)
    ) for i in ids(pm, nw, :storage))

    if bounded
        inj_lb, inj_ub = ref_calc_storage_injection_bounds(ref(pm, nw, :storage), ref(pm, nw, :bus))
        for (i, strg) in ref(pm, nw, :storage)
            for (idx, c) in enumerate(connections[i])
                set_lower_bound(ps[i][c], inj_lb[i][idx])
                set_upper_bound(ps[i][c], inj_ub[i][idx])
            end
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :ps, ids(pm, nw, :storage), ps)
end


"Create variables for `reactive` storage injection"
function variable_mc_storage_power_imaginary_on_off(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i,strg) in ref(pm, nw, :storage))
    qs = var(pm, nw)[:qs] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_qs_$(i)",
        start = comp_start_value(ref(pm, nw, :storage, i), ["qs_start", "qs"], c, 0.0)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :qs, ids(pm, nw, :storage), qs)
end


"""
a reactive power slack variable that enables the storage device to inject or
consume reactive power at its connecting bus, subject to the injection limits
of the device.
"""
function variable_mc_storage_power_control_imaginary_on_off(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    qsc = var(pm, nw)[:qsc] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :storage)], base_name="$(nw)_qsc_$(i)",
        start = comp_start_value(ref(pm, nw, :storage, i), "qsc_start")
    )

    if bounded
        inj_lb, inj_ub = ref_calc_storage_injection_bounds(ref(pm, nw, :storage), ref(pm, nw, :bus))
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :qsc, ids(pm, nw, :storage), qsc)
end


""
function variable_mc_storage_power_mi_on_off(pm::AbstractUnbalancedPowerModel; relax::Bool=false, kwargs...)
    variable_mc_storage_power_real_on_off(pm; kwargs...)
    variable_mc_storage_power_imaginary_on_off(pm; kwargs...)
    variable_mc_storage_current(pm; kwargs...)
    variable_mc_storage_power_control_imaginary_on_off(pm; kwargs...)
    variable_storage_energy(pm; kwargs...)
    variable_storage_charge(pm; kwargs...)
    variable_storage_discharge(pm; kwargs...)
    variable_storage_complementary_indicator(pm; relax=relax, kwargs...)
end


## storage current variables

"do nothing by default but some formulations require this"
function variable_mc_storage_current(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
end


# load variables

"""
Create a dictionary with values of type Any for the load.
Depending on the load model, this can be a parameter or a NLexpression.
These will be inserted into KCL.
"""
function variable_mc_load_power(pm::AbstractUnbalancedPowerModel; nw=nw_id_default, bounded::Bool=true, report::Bool=true)
    var(pm, nw)[:pd] = Dict{Int, Any}()
    var(pm, nw)[:qd] = Dict{Int, Any}()
    var(pm, nw)[:pd_bus] = Dict{Int, Any}()
    var(pm, nw)[:qd_bus] = Dict{Int, Any}()
end


# indicator variables (mld)

"Create variables for demand status"
function variable_mc_load_indicator(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :status, ids(pm, nw, :load), z_demand)
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :pd, ids(pm, nw, :load), pd)
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :qd, ids(pm, nw, :load), qd)
end


"Create variables for shunt status"
function variable_mc_shunt_indicator(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax=false, report::Bool=true)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :shunt, :status, ids(pm, nw, :shunt), z_shunt)
end


"Create variables for bus status"
function variable_mc_bus_voltage_indicator(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :status, ids(pm, nw, :bus), z_voltage)
end


"Create variables for generator status"
function variable_mc_gen_indicator(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :gen_status, ids(pm, nw, :gen), z_gen)
end


"Create variables for storage status"
function variable_mc_storage_indicator(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :status, ids(pm, nw, :storage), z_storage)
end


# slack power variables

"generates variables for both `active` and `reactive` slack at each bus"
function variable_mc_slack_bus_power(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, kwargs...)
    variable_mc_slack_bus_power_real(pm; nw=nw, kwargs...)
    variable_mc_slack_bus_power_imaginary(pm; nw=nw, kwargs...)
end


""
function variable_mc_slack_bus_power_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, report::Bool=true)
    terminals = Dict(i => ref(pm, nw, :bus, i)["terminals"] for i in ids(pm, nw, :bus))
    p_slack = var(pm, nw)[:p_slack] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_p_slack_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "p_slack_start", t, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :p_slack, ids(pm, nw, :bus), p_slack)
end


""
function variable_mc_slack_bus_power_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, report::Bool=true)
    terminals = Dict(i => ref(pm, nw, :bus, i)["terminals"] for i in ids(pm, nw, :bus))
    q_slack = var(pm, nw)[:q_slack] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_q_slack_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "q_slack_start", t, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :q_slack, ids(pm, nw, :bus), q_slack)
end
