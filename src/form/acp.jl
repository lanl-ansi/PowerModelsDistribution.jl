""
function variable_mc_bus_voltage(pm::AbstractUnbalancedACPModel; nw=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_bus_voltage_angle(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_bus_voltage_magnitude_only(pm; nw=nw, bounded=bounded, report=report)

    # This is needed for delta loads, where division occurs by the difference
    # of voltage phasors. If the voltage phasors at one bus are initialized
    # in the same point, this would lead to division by zero.

    for id in ids(pm, nw, :bus)
        busref = ref(pm, nw, :bus, id)
        terminals = busref["terminals"]

        ncnd = length(terminals)

        vm_start = fill(1.0, 3)
        for t in 1:3
            if t in terminals
                vmax = busref["vmax"][findfirst(isequal(t), terminals)]
                vm_start[t] = min(vm_start[t], vmax)

                vmin = busref["vmin"][findfirst(isequal(t), terminals)]
                vm_start[t] = max(vm_start[t], vmin)
            end
        end

        vm = haskey(busref, "vm_start") ? busref["vm_start"] : haskey(busref, "vm") ? busref["vm"] : [vm_start..., fill(0.0, ncnd)...][terminals]
        va = haskey(busref, "va_start") ? busref["va_start"] : haskey(busref, "va") ? busref["va"] : [deg2rad.([0, -120, 120])..., zeros(length(terminals))...][terminals]

        for (idx,t) in enumerate(terminals)
            JuMP.set_start_value(var(pm, nw, :vm, id)[t], vm[idx])
            JuMP.set_start_value(var(pm, nw, :va, id)[t], va[idx])
        end
    end
end


""
function variable_mc_bus_voltage_on_off(pm::AbstractUnbalancedACPModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_bus_voltage_angle(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_bus_voltage_magnitude_on_off(pm; nw=nw, bounded=bounded, report=report)

    for id in ids(pm, nw, :bus)
        busref = ref(pm, nw, :bus, id)
        terminals = busref["terminals"]

        ncnd = length(terminals)

        vm_start = fill(1.0, 3)
        for t in 1:3
            if t in terminals
                vmax = busref["vmax"][findfirst(isequal(t), terminals)]
                vm_start[t] = min(vm_start[t], vmax)

                vmin = busref["vmin"][findfirst(isequal(t), terminals)]
                vm_start[t] = max(vm_start[t], vmin)
            end
        end

        vm = haskey(busref, "vm_start") ? busref["vm_start"] : haskey(busref, "vm") ? busref["vm"] : [vm_start..., fill(0.0, ncnd)...][terminals]
        va = haskey(busref, "va_start") ? busref["va_start"] : haskey(busref, "va") ? busref["va"] : [[_wrap_to_pi(2 * pi / 3 * (1-t)) for t in 1:3]..., zeros(length(terminals))...][terminals]

        for (idx,t) in enumerate(terminals)
            JuMP.set_start_value(var(pm, nw, :vm, id)[t], vm[idx])
            JuMP.set_start_value(var(pm, nw, :va, id)[t], va[idx])
        end
    end
end


""
function constraint_mc_switch_state_closed(pm::AbstractUnbalancedACPModel, nw::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int})
    vm_fr = var(pm, nw, :vm, f_bus)
    vm_to = var(pm, nw, :vm, t_bus)

    va_fr = var(pm, nw, :va, f_bus)
    va_to = var(pm, nw, :va, t_bus)

    for (idx,(fc,tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, vm_fr[fc] == vm_to[tc])
        JuMP.@constraint(pm.model, va_fr[fc] == va_to[fc])
    end
end


""
function constraint_mc_switch_state_on_off(pm::AbstractUnbalancedACPModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int}; relax::Bool=false)
    vm_fr = var(pm, nw, :vm, f_bus)
    vm_to = var(pm, nw, :vm, t_bus)

    va_fr = var(pm, nw, :va, f_bus)
    va_to = var(pm, nw, :va, t_bus)

    z = var(pm, nw, :switch_state, i)

    for (idx,(fc,tc)) in enumerate(zip(f_connections, t_connections))
        if relax
            vm_fr_ub = JuMP.has_upper_bound(vm_fr[fc]) ? JuMP.upper_bound(vm_fr[fc]) :  1e20
            vm_to_lb = JuMP.has_lower_bound(vm_to[tc]) ? JuMP.lower_bound(vm_to[tc]) : -1e20
            va_fr_ub = JuMP.has_upper_bound(va_fr[tc]) ? JuMP.upper_bound(va_fr[fc]) :  1e20
            va_to_lb = JuMP.has_lower_bound(va_to[tc]) ? JuMP.lower_bound(va_to[tc]) : -1e20

            M = 1e20
            JuMP.@constraint(pm.model, vm_fr[fc] - vm_to[tc] <=  M * (1-z))
            JuMP.@constraint(pm.model, vm_fr[fc] - vm_to[tc] >= -M * (1-z))
            JuMP.@constraint(pm.model, va_fr[fc] - va_to[tc] <=  M * (1-z))
            JuMP.@constraint(pm.model, va_fr[fc] - va_to[tc] >= -M * (1-z))
        else
            JuMP.@constraint(pm.model, z => {vm_fr[fc] == vm_to[tc]})
            JuMP.@constraint(pm.model, z => {va_fr[fc] == va_to[tc]})
        end
    end
end


""
function constraint_mc_power_balance_slack(pm::AbstractUnbalancedACPModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vm   = var(pm, nw, :vm, i)
    va   = var(pm, nw, :va, i)
    p    = get(var(pm, nw),    :p, Dict()); _check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(var(pm, nw),    :q, Dict()); _check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, nw),   :pg, Dict()); _check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, nw),   :qg, Dict()); _check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, nw),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, nw),   :qs, Dict()); _check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, nw),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, nw),  :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(var(pm, nw),   :pt, Dict()); _check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(var(pm, nw),   :qt, Dict()); _check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    p_slack = var(pm, nw, :p_slack, i)
    q_slack = var(pm, nw, :q_slack, i)

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    ncnds = length(terminals)
    Pd = fill(0.0, ncnds)
    Qd = fill(0.0, ncnds)
    for (ld_i, connections) in bus_loads
        load = ref(pm, nw, :load, ld_i)
        for (idx, c) in enumerate(connections)
            Pd[findfirst(isequal(c), terminals)] += load["pd"][idx]
            Qd[findfirst(isequal(c), terminals)] += load["qd"][idx]
        end
    end

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        cp = JuMP.@NLconstraint(pm.model,
              sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
            + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
            ==
              sum(pg[g][t] for (g, conns) in bus_gens if t in conns)
            - sum(ps[s][t] for (s, conns) in bus_storage if t in conns)
            - Pd[idx]
            - ( # shunt
                Gt[idx,idx] * vm[t]^2
                +sum( Gt[idx,jdx] * vm[t]*vm[u] * cos(va[t]-va[u])
                     +Bt[idx,jdx] * vm[t]*vm[u] * sin(va[t]-va[u])
                     for (jdx,u) in ungrounded_terminals if idx != jdx)
            )
            + p_slack[t]
        )
        push!(cstr_p, cp)

        cq = JuMP.@NLconstraint(pm.model,
              sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
            + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
            ==
              sum(qg[g][t] for (g, conns) in bus_gens if t in conns)
            - sum(qs[s][t] for (s, conns) in bus_storage if t in conns)
            - Qd[idx]
            - ( # shunt
                -Bt[idx,idx] * vm[t]^2
                -sum( Bt[idx,jdx] * vm[t]*vm[u] * cos(va[t]-va[u])
                     -Gt[idx,jdx] * vm[t]*vm[u] * sin(va[t]-va[u])
                     for (jdx,u) in ungrounded_terminals if idx != jdx)
            )
            + q_slack[t]
        )
        push!(cstr_q, cq)

    end

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = cstr_q

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_mc_power_balance_shed(pm::AbstractUnbalancedACPModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vm       = var(pm, nw, :vm, i)
    va       = var(pm, nw, :va, i)
    p        = get(var(pm, nw),    :p, Dict()); _check_var_keys(p, bus_arcs, "active power", "branch")
    q        = get(var(pm, nw),    :q, Dict()); _check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg       = get(var(pm, nw),   :pg, Dict()); _check_var_keys(pg, bus_gens, "active power", "generator")
    qg       = get(var(pm, nw),   :qg, Dict()); _check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps       = get(var(pm, nw),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    qs       = get(var(pm, nw),   :qs, Dict()); _check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw      = get(var(pm, nw),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw      = get(var(pm, nw),  :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt       = get(var(pm, nw),   :pt, Dict()); _check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt       = get(var(pm, nw),   :qt, Dict()); _check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")

    z_demand = var(pm, nw, :z_demand)
    z_gen = haskey(var(pm, nw), :z_gen) ? var(pm, nw, :z_gen) : Dict(i => 1.0 for i in ids(pm, nw, :gen))
    z_storage = haskey(var(pm, nw), :z_storage) ? var(pm, nw, :z_storage) : Dict(i => 1.0 for i in ids(pm, nw, :storage))
    z_shunt  = haskey(var(pm, nw), :z_shunt) ? var(pm, nw, :z_shunt) : Dict(i => 1.0 for i in ids(pm, nw, :shunt))

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        cp = JuMP.@NLconstraint(pm.model,
              sum(     p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(  pt[a_t][t] for (a_t, conns) in bus_arcs_trans if t in conns)
            - sum(pg[g][t]*z_gen[g] for (g, conns) in bus_gens if t in conns)
            + sum(ps[s][t]*z_storage[s] for (s, conns) in bus_storage if t in conns)
            + sum(ref(pm, nw, :load, d, "pd")[findfirst(isequal(t), conns)]*z_demand[d] for (d, conns) in bus_loads if t in conns)
            + sum(z_shunt[s] *
                (ref(pm, nw, :shunt, s)["gs"][findfirst(isequal(t), conns), findfirst(isequal(t), conns)] * vm[t]^2
                +sum( ref(pm, nw, :shunt, s)["gs"][findfirst(isequal(t), conns), findfirst(isequal(u), conns)] * vm[t]*vm[u] * cos(va[t]-va[u])
                     +ref(pm, nw, :shunt, s)["bs"][findfirst(isequal(t), conns), findfirst(isequal(u), conns)] * vm[t]*vm[u] * sin(va[t]-va[u])
                for (jdx, u) in ungrounded_terminals if idx != jdx ) )
            for (s, conns) in bus_shunts if t in conns )
            ==
            0.0
        )
        push!(cstr_p, cp)

        cq = JuMP.@NLconstraint(pm.model,
              sum(     q[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(qsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(  qt[a_t][t] for (a_t, conns) in bus_arcs_trans if t in conns)
            - sum(    qg[g][t]*z_gen[g] for (g, conns) in bus_gens if t in conns)
            + sum(    qs[s][t]*z_storage[s] for (s, conns) in bus_storage if t in conns)
            + sum(ref(pm, nw, :load, l, "qd")[findfirst(isequal(t), conns)]*z_demand[l] for (l, conns) in bus_loads if t in conns)
            + sum(z_shunt[sh] *
                (-ref(pm, nw, :shunt, sh)["bs"][findfirst(isequal(t), conns), findfirst(isequal(t), conns)] * vm[t]^2
                 -sum( ref(pm, nw, :shunt, sh)["bs"][findfirst(isequal(t), conns), findfirst(isequal(u), conns)] * vm[t]*vm[u] * cos(va[t]-va[u])
                      -ref(pm, nw, :shunt, sh)["gs"][findfirst(isequal(t), conns), findfirst(isequal(u), conns)] * vm[t]*vm[u] * sin(va[t]-va[u])
                for (jdx, u) in ungrounded_terminals if idx != jdx ) )
            for (sh, conns) in bus_shunts if t in conns )
            ==
            0.0
        )
        push!(cstr_q, cq)
    end

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = cstr_q

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_mc_power_balance_simple(pm::AbstractUnbalancedACPModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vm   = var(pm, nw, :vm, i)
    va   = var(pm, nw, :va, i)
    p    = get(var(pm, nw),    :p, Dict()); _check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(var(pm, nw),    :q, Dict()); _check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, nw),   :pg, Dict()); _check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, nw),   :qg, Dict()); _check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, nw),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, nw),   :qs, Dict()); _check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, nw),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, nw),  :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(var(pm, nw),   :pt, Dict()); _check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(var(pm, nw),   :qt, Dict()); _check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    ncnds = length(terminals)
    Pd = fill(0.0, ncnds)
    Qd = fill(0.0, ncnds)
    for (ld_i, connections) in bus_loads
        load = ref(pm, nw, :load, ld_i)
        for (idx, c) in enumerate(connections)
            Pd[findfirst(isequal(c), terminals)] += load["pd"][idx]
            Qd[findfirst(isequal(c), terminals)] += load["qd"][idx]
        end
    end

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        cp = JuMP.@NLconstraint(pm.model,
              sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
            + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
            ==
              sum(pg[g][t] for (g, conns) in bus_gens)
            - sum(ps[s][t] for (s, conns) in bus_storage)
            - Pd[idx]
            - ( # shunt
                Gt[idx,idx] * vm[t]^2
                +sum( Gt[idx,jdx] * vm[t]*vm[u] * cos(va[t]-va[u])
                     +Bt[idx,jdx] * vm[t]*vm[u] * sin(va[t]-va[u])
                     for (jdx,u) in ungrounded_terminals if idx != jdx)
            )
        )
        push!(cstr_p, cp)

        cq = JuMP.@NLconstraint(pm.model,
            sum(q[a][c] for a in bus_arcs)
            + sum(qsw[a_sw][c] for a_sw in bus_arcs_sw)
            + sum(qt[a_trans][c] for a_trans in bus_arcs_trans)
            ==
            sum(qg[g][c] for g in bus_gens)
            - sum(qs[s][c] for s in bus_storage)
            - Qd[idx]
            - ( # shunt
                -Bt[idx,idx] * vm[t]^2
                -sum( Bt[idx,jdx] * vm[t]*vm[u] * cos(va[t]-va[u])
                     -Gt[idx,jdx] * vm[t]*vm[u] * sin(va[t]-va[u])
                     for (jdx,u) in ungrounded_terminals if idx != jdx)
            )
        )
        push!(cstr_q, cq)
    end

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = cstr_q

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_mc_power_balance(pm::AbstractUnbalancedACPModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vm   = var(pm, nw, :vm, i)
    va   = var(pm, nw, :va, i)
    p    = get(var(pm, nw),      :p, Dict()); _check_var_keys(  p, bus_arcs, "active power", "branch")
    q    = get(var(pm, nw),      :q, Dict()); _check_var_keys(  q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, nw), :pg_bus, Dict()); _check_var_keys( pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, nw), :qg_bus, Dict()); _check_var_keys( qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, nw),     :ps, Dict()); _check_var_keys( ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, nw),     :qs, Dict()); _check_var_keys( qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, nw),    :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, nw),    :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(var(pm, nw),     :pt, Dict()); _check_var_keys( pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(var(pm, nw),     :qt, Dict()); _check_var_keys( qt, bus_arcs_trans, "reactive power", "transformer")
    pd   = get(var(pm, nw), :pd_bus, Dict()); _check_var_keys( pd, bus_loads, "active power", "load")
    qd   = get(var(pm, nw), :qd_bus, Dict()); _check_var_keys( pd, bus_loads, "reactive power", "load")

    Gs, Bs = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []
    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        if any(Bs[idx,jdx] != 0 for (jdx, u) in ungrounded_terminals if idx != jdx) || any(Gs[idx,jdx] != 0 for (jdx, u) in ungrounded_terminals if idx != jdx)
            cp = JuMP.@NLconstraint(pm.model,
                  sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( pd[l][t] for (l, conns) in bus_loads if t in conns)
                + ( # shunt
                    +Gs[idx,idx] * vm[t]^2
                    +sum( Gs[idx,jdx] * vm[t]*vm[u] * cos(va[t]-va[u])
                         +Bs[idx,jdx] * vm[t]*vm[u] * sin(va[t]-va[u])
                        for (jdx,u) in ungrounded_terminals if idx != jdx)
                )
                ==
                0.0
            )
            push!(cstr_p, cp)

            cq = JuMP.@NLconstraint(pm.model,
                  sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( qd[l][t] for (l, conns) in bus_loads if t in conns)
                + ( # shunt
                    -Bs[idx,idx] * vm[t]^2
                    -sum( Bs[idx,jdx] * vm[t]*vm[u] * cos(va[t]-va[u])
                         -Gs[idx,jdx] * vm[t]*vm[u] * sin(va[t]-va[u])
                         for (jdx,u) in ungrounded_terminals if idx != jdx)
                )
                ==
                0.0
            )
            push!(cstr_q, cq)
        else
            cp = @smart_constraint(pm.model, [p, pg, ps, psw, pt, pd, vm],
                  sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( pd[l][t] for (l, conns) in bus_loads if t in conns)
                + Gs[idx,idx] * vm[t]^2
                ==
                0.0
            )
            push!(cstr_p, cp)

            cq = @smart_constraint(pm.model, [q, qg, qs, qsw, qt, qd, vm],
                  sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( qd[l][t] for (l, conns) in bus_loads if t in conns)
                - Bs[idx,idx] * vm[t]^2
                ==
                0.0
            )
            push!(cstr_q, cq)
        end
    end

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = cstr_q

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


@doc raw"""
    constraint_mc_power_balance_capc(pm::AbstractUnbalancedACPModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})

Power balance constraints with capacitor control.

```math
\begin{align}
    & Bs = z ⋅ bs, \\
    &\text{capacitor ON: }  z = 1, \\
    &\text{capacitor OFF: } z = 0.
\end{align}
```
"""
function constraint_mc_power_balance_capc(pm::AbstractUnbalancedACPModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vm   = var(pm, nw, :vm, i)
    va   = var(pm, nw, :va, i)
    p    = get(var(pm, nw),      :p, Dict()); _check_var_keys(  p, bus_arcs, "active power", "branch")
    q    = get(var(pm, nw),      :q, Dict()); _check_var_keys(  q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, nw), :pg_bus, Dict()); _check_var_keys( pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, nw), :qg_bus, Dict()); _check_var_keys( qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, nw),     :ps, Dict()); _check_var_keys( ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, nw),     :qs, Dict()); _check_var_keys( qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, nw),    :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, nw),    :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(var(pm, nw),     :pt, Dict()); _check_var_keys( pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(var(pm, nw),     :qt, Dict()); _check_var_keys( qt, bus_arcs_trans, "reactive power", "transformer")
    pd   = get(var(pm, nw), :pd_bus, Dict()); _check_var_keys( pd, bus_loads, "active power", "load")
    qd   = get(var(pm, nw), :qd_bus, Dict()); _check_var_keys( pd, bus_loads, "reactive power", "load")

    # add constraints to model capacitor switching
    if !isempty(bus_shunts) && haskey(ref(pm, nw, :shunt, bus_shunts[1][1]), "controls")
        constraint_capacitor_on_off(pm, nw, i, bus_shunts)
    end

    # calculate Gs, Bs
    cap_state = 1.0
    ncnds = length(terminals)
    Gs = fill(0.0, ncnds, ncnds)
    Bs = convert(Matrix{JuMP.NonlinearExpression}, JuMP.@NLexpression(pm.model, [idx=1:ncnds, jdx=1:ncnds], 0.0))
    for (val, connections) in bus_shunts
        shunt = ref(pm,nw,:shunt,val)
        for (idx,c) in enumerate(connections)
            if haskey(shunt, "controls")
                cap_state = var(pm, nw, :capacitor_state, val)[c]
            end
            for (jdx,d) in enumerate(connections)
                Gs[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] += shunt["gs"][idx,jdx]
                Bs[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] = JuMP.@NLexpression(pm.model, Bs[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] + shunt["bs"][idx,jdx]*cap_state)
            end
        end
    end

    cstr_p = []
    cstr_q = []
    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        if any(Bs[idx,jdx] != 0 for (jdx, u) in ungrounded_terminals if idx != jdx) || any(Gs[idx,jdx] != 0 for (jdx, u) in ungrounded_terminals if idx != jdx)
            cp = JuMP.@NLconstraint(pm.model,
                  sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( pd[l][t] for (l, conns) in bus_loads if t in conns)
                + ( # shunt
                    +Gs[idx,idx] * vm[t]^2
                    +sum( Gs[idx,jdx] * vm[t]*vm[u] * cos(va[t]-va[u])
                         +Bs[idx,jdx] * vm[t]*vm[u] * sin(va[t]-va[u])
                        for (jdx,u) in ungrounded_terminals if idx != jdx)
                )
                ==
                0.0
            )
            push!(cstr_p, cp)

            cq = JuMP.@NLconstraint(pm.model,
                  sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( qd[l][t] for (l, conns) in bus_loads if t in conns)
                + ( # shunt
                    -Bs[idx,idx] * vm[t]^2
                    -sum( Bs[idx,jdx] * vm[t]*vm[u] * cos(va[t]-va[u])
                         -Gs[idx,jdx] * vm[t]*vm[u] * sin(va[t]-va[u])
                         for (jdx,u) in ungrounded_terminals if idx != jdx)
                )
                ==
                0.0
            )
            push!(cstr_q, cq)
        else
            cp = @smart_constraint(pm.model, [p, pg, ps, psw, pt, pd, vm],
                  sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( pd[l][t] for (l, conns) in bus_loads if t in conns)
                + Gs[idx,idx] * vm[t]^2
                ==
                0.0
            )
            push!(cstr_p, cp)

            cq = @smart_constraint(pm.model, [q, qg, qs, qsw, qt, qd, vm, Bs],
                  sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( qd[l][t] for (l, conns) in bus_loads if t in conns)
                - Bs[idx,idx] * vm[t]^2
                ==
                0.0
            )
            push!(cstr_q, cq)
        end
    end

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = cstr_q

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


@doc raw"""
    constraint_capacitor_on_off(pm::AbstractUnbalancedACPModel, i::Int; nw::Int=nw_id_default)

Add constraints to model capacitor switching

```math
\begin{align}
&\text{kvar control (ON): }  q-q_\text{on} ≤ M_q ⋅ z - ϵ ⋅ (1-z), \\
&\text{kvar control (OFF): } q-q_\text{off} ≥ -M_q ⋅ (1-z) - ϵ ⋅ z, \\
&\text{voltage control (ON): }  v-v_\text{min} ≥ -M_v ⋅ z + ϵ ⋅ (1-z), \\
&\text{voltage control (OFF): } v-v_\text{max} ≤ M_v ⋅ (1-z) - ϵ ⋅ z.
\end{align}
```
"""
function constraint_capacitor_on_off(pm::AbstractUnbalancedACPModel, nw::Int, i::Int, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    cap_state = var(pm, nw, :capacitor_state, bus_shunts[1][1])
    shunt = ref(pm, nw, :shunt, bus_shunts[1][1])
    ϵ = 1e-5
    M_q = 1e5
    M_v = 2
    elem_type = shunt["controls"]["element"]["type"]
    if shunt["controls"]["type"] == CAP_REACTIVE_POWER
        bus_idx = shunt["controls"]["terminal"] == 1 ? (shunt["controls"]["element"]["index"], shunt["controls"]["element"]["f_bus"], shunt["controls"]["element"]["t_bus"]) : (shunt["controls"]["element"]["index"], shunt["controls"]["element"]["t_bus"], shunt["controls"]["element"]["f_bus"])
        q_fr = elem_type == "branch" ? var(pm, nw, :q)[bus_idx] : elem_type == "switch" ? var(pm, nw, :qsw) : var(pm, nw, :qt, bus_idx)
        JuMP.@constraint(pm.model, sum(q_fr) - shunt["controls"]["onsetting"] ≤ M_q*cap_state[shunt["connections"][1]] - ϵ*(1-cap_state[shunt["connections"][1]]))
        JuMP.@constraint(pm.model, sum(q_fr) - shunt["controls"]["offsetting"] ≥ -M_q*(1-cap_state[shunt["connections"][1]]) - ϵ*cap_state[shunt["connections"][1]])
        JuMP.@constraint(pm.model, cap_state .== cap_state[shunt["connections"][1]])
        if shunt["controls"]["voltoverride"]
            for (idx,val) in enumerate(shunt["connections"])
                vm_cap = var(pm, nw, :vm, i)[val]
                JuMP.@constraint(pm.model, vm_cap - shunt["controls"]["vmin"] ≥ -M_v*cap_state[val] + ϵ*(1-cap_state[val]))
                JuMP.@constraint(pm.model, vm_cap - shunt["controls"]["vmax"] ≤ M_v*(1-cap_state[val]) - ϵ*cap_state[val])
            end
        end
    else
        for (idx,val) in enumerate(shunt["connections"])
            if shunt["controls"]["type"][idx] == CAP_VOLTAGE
                bus_idx = shunt["controls"]["terminal"][idx] == 1 ? shunt["controls"]["element"]["f_bus"] : shunt["controls"]["element"]["t_bus"]
                vm_cap = var(pm, nw, :vm, bus_idx)[val]
                JuMP.@constraint(pm.model, vm_cap - shunt["controls"]["onsetting"][idx] ≤ M_v*cap_state[val] - ϵ*(1-cap_state[val]))
                JuMP.@constraint(pm.model, vm_cap - shunt["controls"]["offsetting"][idx] ≥ -M_v*(1-cap_state[val]) - ϵ*cap_state[val])
            end
            if shunt["controls"]["voltoverride"][idx]
                vm_cap = var(pm, nw, :vm, i)[val]
                JuMP.@constraint(pm.model, vm_cap - shunt["controls"]["vmin"][idx] ≥ -M_v*cap_state[val] + ϵ*(1-cap_state[val]))
                JuMP.@constraint(pm.model, vm_cap - shunt["controls"]["vmax"][idx] ≤ M_v*(1-cap_state[val]) - ϵ*cap_state[val])
            end
            if shunt["controls"]["type"][idx] == CAP_DISABLED
                JuMP.@constraint(pm.model, cap_state[val] == 1 )
            end
        end
    end
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p_fr ==     g[c,c] * vm_fr[c]^2 +
            sum( g[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                 b[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in conductor_ids(pm) if d != c) +
            sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in conductor_ids(pm))
            + g_fr[c,c] * vm_fr[c]^2 +
            sum( g_fr[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                 b_fr[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in conductor_ids(pm) if d != c)
            )
q_fr == -b[c,c] *vm_fr[c]^2 -
            sum( b[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                 g[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in conductor_ids(pm) if d != c) -
            sum(-b[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                 g[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in conductor_ids(pm))
            -b_fr[c,c] *vm_fr[c]^2 -
            sum( b_fr[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                 g_fr[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in conductor_ids(pm) if d != c)
            )
```
"""
function constraint_mc_ohms_yt_from(pm::AbstractUnbalancedACPModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})
    p_fr  = var(pm, nw,  :p, f_idx)
    q_fr  = var(pm, nw,  :q, f_idx)
    vm_fr = var(pm, nw, :vm, f_bus)
    vm_to = var(pm, nw, :vm, t_bus)
    va_fr = var(pm, nw, :va, f_bus)
    va_to = var(pm, nw, :va, t_bus)

    ohms_yt_p = JuMP.ConstraintRef[]
    ohms_yt_q = JuMP.ConstraintRef[]
    for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))
        push!(ohms_yt_p, JuMP.@NLconstraint(pm.model, p_fr[fc] == (G[idx,idx]+G_fr[idx,idx])*vm_fr[fc]^2
            +sum( (G[idx,jdx]+G_fr[idx,jdx]) * vm_fr[fc]*vm_fr[fd]*cos(va_fr[fc]-va_fr[fd])
                 +(B[idx,jdx]+B_fr[idx,jdx]) * vm_fr[fc]*vm_fr[fd]*sin(va_fr[fc]-va_fr[fd])
                for (jdx, (fd,td)) in enumerate(zip(f_connections,t_connections)) if idx != jdx)
            +sum( -G[idx,jdx]*vm_fr[fc]*vm_to[td]*cos(va_fr[fc]-va_to[td])
                  -B[idx,jdx]*vm_fr[fc]*vm_to[td]*sin(va_fr[fc]-va_to[td])
                for (jdx, (fd,td)) in enumerate(zip(f_connections,t_connections)))
            )
        )

        push!(ohms_yt_q, JuMP.@NLconstraint(pm.model, q_fr[fc] == -(B[idx,idx]+B_fr[idx,idx])*vm_fr[fc]^2
            -sum( (B[idx,jdx]+B_fr[idx,jdx])*vm_fr[fc]*vm_fr[fd]*cos(va_fr[fc]-va_fr[fd])
                 -(G[idx,jdx]+G_fr[idx,jdx])*vm_fr[fc]*vm_fr[fd]*sin(va_fr[fc]-va_fr[fd])
                for (jdx, (fd,td)) in enumerate(zip(f_connections,t_connections)) if idx != jdx)
            -sum(-B[idx,jdx]*vm_fr[fc]*vm_to[td]*cos(va_fr[fc]-va_to[td])
                 +G[idx,jdx]*vm_fr[fc]*vm_to[td]*sin(va_fr[fc]-va_to[td])
                for (jdx, (fd,td)) in enumerate(zip(f_connections,t_connections)))
            )
        )
    end
    con(pm, nw, :ohms_yt)[f_idx] = [ohms_yt_p, ohms_yt_q]
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_mc_ohms_yt_to(pm::AbstractUnbalancedACPModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_to::Matrix{<:Real}, B_to::Matrix{<:Real})
    constraint_mc_ohms_yt_from(pm, nw, t_bus, f_bus, t_idx, f_idx, t_connections, f_connections, G, B, G_to, B_to)
end


""
function constraint_mc_transformer_power_yy(pm::AbstractUnbalancedACPModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    transformer = ref(pm, nw, :transformer, trans_id)

    vm_fr = var(pm, nw, :vm, f_bus)
    vm_to = var(pm, nw, :vm, t_bus)
    va_fr = var(pm, nw, :va, f_bus)
    va_to = var(pm, nw, :va, t_bus)

    p_fr = [var(pm, nw, :pt, f_idx)[c] for c in f_connections]
    p_to = [var(pm, nw, :pt, t_idx)[c] for c in t_connections]
    q_fr = [var(pm, nw, :qt, f_idx)[c] for c in f_connections]
    q_to = [var(pm, nw, :qt, t_idx)[c] for c in t_connections]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))]

    for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))
        if tm_fixed[idx]
            JuMP.@constraint(pm.model, vm_fr[fc] == tm_scale*tm[idx]*vm_to[tc])
        else
            # transformer taps without regcontrol, tap variable not required in regcontrol formulation
            JuMP.@constraint(pm.model, vm_fr[fc] == tm_scale*tm[idx]*vm_to[tc])

            # with regcontrol
            if haskey(transformer,"controls")
                v_ref = transformer["controls"]["vreg"][idx]
                δ = transformer["controls"]["band"][idx]
                r = transformer["controls"]["r"][idx]
                x = transformer["controls"]["x"][idx]

                # (cr+jci) = (p-jq)/(vm⋅cos(va)-jvm⋅sin(va))
                cr = JuMP.@NLexpression(pm.model, ( p_to[idx]*vm_to[tc]*cos(va_to[tc]) + q_to[idx]*vm_to[tc]*sin(va_to[tc]))/vm_to[tc]^2)
                ci = JuMP.@NLexpression(pm.model, (-q_to[idx]*vm_to[tc]*cos(va_to[tc]) + p_to[idx]*vm_to[tc]*sin(va_to[tc]))/vm_to[tc]^2)
                # v_drop = (cr+jci)⋅(r+jx)
                vr_drop = JuMP.@NLexpression(pm.model, r*cr-x*ci)
                vi_drop = JuMP.@NLexpression(pm.model, r*ci+x*cr)

                # v_ref-δ ≤ vm_fr-(cr+jci)⋅(r+jx)≤ v_ref+δ
                # vm_fr/1.1 ≤ vm_to ≤ vm_fr/0.9
                JuMP.@NLconstraint(pm.model, (vm_fr[fc]*cos(va_fr[fc])-vr_drop)^2 + (vm_fr[fc]*sin(va_fr[fc])-vi_drop)^2 ≥ (v_ref - δ)^2)
                JuMP.@NLconstraint(pm.model, (vm_fr[fc]*cos(va_fr[fc])-vr_drop)^2 + (vm_fr[fc]*sin(va_fr[fc])-vi_drop)^2 ≤ (v_ref + δ)^2)
                JuMP.@constraint(pm.model, vm_fr[fc]/1.1 ≤ vm_to[tc])
                JuMP.@constraint(pm.model, vm_fr[fc]/0.9 ≥ vm_to[tc])
            end
        end
        pol_angle = pol == 1 ? 0 : pi
        JuMP.@constraint(pm.model, va_fr[fc] == va_to[tc] + pol_angle)
    end

    JuMP.@constraint(pm.model, p_fr + p_to .== 0)
    JuMP.@constraint(pm.model, q_fr + q_to .== 0)
end


""
function constraint_mc_transformer_power_dy(pm::AbstractUnbalancedACPModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    vm_fr = var(pm, nw, :vm, f_bus)
    vm_to = var(pm, nw, :vm, t_bus)
    va_fr = var(pm, nw, :va, f_bus)
    va_to = var(pm, nw, :va, t_bus)

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[fc] for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))]

    nph = length(tm_set)

    # introduce auxialiary variable vd = Md*v_fr
    vd_re = Array{Any,1}(undef, nph)
    vd_im = Array{Any,1}(undef, nph)
    for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))
        # rotate by 1 to get 'previous' phase
        # e.g., for nph=3: 1->3, 2->1, 3->2
        jdx = (idx-1+1)%nph+1
        fd = f_connections[jdx]
        vd_re[idx] = JuMP.@NLexpression(pm.model, vm_fr[fc]*cos(va_fr[fc])-vm_fr[fd]*cos(va_fr[fd]))
        vd_im[idx] = JuMP.@NLexpression(pm.model, vm_fr[fc]*sin(va_fr[fc])-vm_fr[fd]*sin(va_fr[fd]))
        JuMP.@NLconstraint(pm.model, vd_re[idx] == pol*tm_scale*tm[idx]*vm_to[tc]*cos(va_to[tc]))
        JuMP.@NLconstraint(pm.model, vd_im[idx] == pol*tm_scale*tm[idx]*vm_to[tc]*sin(va_to[tc]))
    end

    p_fr = var(pm, nw, :pt, f_idx)
    p_to = var(pm, nw, :pt, t_idx)
    q_fr = var(pm, nw, :qt, f_idx)
    q_to = var(pm, nw, :qt, t_idx)

    id_re = Array{Any,1}(undef, nph)
    id_im = Array{Any,1}(undef, nph)
    # s/v      = (p+jq)/|v|^2*conj(v)
    #          = (p+jq)/|v|*(cos(va)-j*sin(va))
    # Re(s/v)  = (p*cos(va)+q*sin(va))/|v|
    # -Im(s/v) = -(q*cos(va)-p*sin(va))/|v|
    for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))
        # id = conj(s_to/v_to)./tm
        id_re[idx] = JuMP.@NLexpression(pm.model,  (p_to[tc]*cos(va_to[tc])+q_to[tc]*sin(va_to[tc]))/vm_to[tc]/(tm_scale*tm[idx])/pol)
        id_im[idx] = JuMP.@NLexpression(pm.model, -(q_to[tc]*cos(va_to[tc])-p_to[tc]*sin(va_to[tc]))/vm_to[tc]/(tm_scale*tm[idx])/pol)
    end
    for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))
        # rotate by nph-1 to get 'previous' phase
        # e.g., for nph=3: 1->3, 2->1, 3->2
        jdx = (idx-1+nph-1)%nph+1
        # s_fr  = v_fr*conj(i_fr)
        #       = v_fr*conj(id[q]-id[p])
        #       = v_fr*(id_re[q]-j*id_im[q]-id_re[p]+j*id_im[p])
        JuMP.@NLconstraint(pm.model, p_fr[fc] ==
             vm_fr[fc]*cos(va_fr[fc])*( id_re[jdx]-id_re[idx])
            -vm_fr[fc]*sin(va_fr[fc])*(-id_im[jdx]+id_im[idx])
        )
        JuMP.@NLconstraint(pm.model, q_fr[fc] ==
             vm_fr[fc]*cos(va_fr[fc])*(-id_im[jdx]+id_im[idx])
            +vm_fr[fc]*sin(va_fr[fc])*( id_re[jdx]-id_re[idx])
        )
    end
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_vuf(pm::AbstractUnbalancedACPModel, nw::Int, bus_id::Int, vufmax::Real)
    if !haskey(var(pm, nw_id_default), :vmpossqr)
        var(pm, nw_id_default)[:vmpossqr] = Dict{Int, Any}()
        var(pm, nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [var(pm, nw, :vm, bus_id)[i] for i in 1:3]
    (va_a, va_b, va_c) = [var(pm, nw, :va, bus_id)[i] for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)
    # real and imaginary components of U+
    vrepos = JuMP.@NLexpression(pm.model,
        (vm_a*cos(va_a) + are*vm_b*cos(va_b) - aim*vm_b*sin(va_b) + a2re*vm_c*cos(va_c) - a2im*vm_c*sin(va_c))/3
    )
    vimpos = JuMP.@NLexpression(pm.model,
        (vm_a*sin(va_a) + are*vm_b*sin(va_b) + aim*vm_b*cos(va_b) + a2re*vm_c*sin(va_c) + a2im*vm_c*cos(va_c))/3
    )
    # square of magnitude of U+, |U+|^2
    vmpossqr = JuMP.@NLexpression(pm.model, vrepos^2+vimpos^2)
    # real and imaginary components of U-
    vreneg = JuMP.@NLexpression(pm.model,
        (vm_a*cos(va_a) + a2re*vm_b*cos(va_b) - a2im*vm_b*sin(va_b) + are*vm_c*cos(va_c) - aim*vm_c*sin(va_c))/3
    )
    vimneg = JuMP.@NLexpression(pm.model,
        (vm_a*sin(va_a) + a2re*vm_b*sin(va_b) + a2im*vm_b*cos(va_b) + are*vm_c*sin(va_c) + aim*vm_c*cos(va_c))/3
    )
    # square of magnitude of U-, |U-|^2
    vmnegsqr = JuMP.@NLexpression(pm.model, vreneg^2+vimneg^2)
    # finally, apply constraint
    JuMP.@NLconstraint(pm.model, vmnegsqr <= vufmax^2*vmpossqr)
    # DEBUGGING: save references for post check
    #var(pm, nw_id_default, :vmpossqr)[bus_id] = vmpossqr
    #var(pm, nw_id_default, :vmnegsqr)[bus_id] = vmnegsqr
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_negative_sequence(pm::AbstractUnbalancedACPModel, nw::Int, bus_id::Int, vmnegmax::Real)
    if !haskey(var(pm, nw_id_default), :vmpossqr)
        var(pm, nw_id_default)[:vmpossqr] = Dict{Int, Any}()
        var(pm, nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [var(pm, nw, :vm, bus_id)[i] for i in 1:3]
    (va_a, va_b, va_c) = [var(pm, nw, :va, bus_id)[i] for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)
    # real and imaginary components of U-
    vreneg = JuMP.@NLexpression(pm.model,
        (vm_a*cos(va_a) + a2re*vm_b*cos(va_b) - a2im*vm_b*sin(va_b) + are*vm_c*cos(va_c) - aim*vm_c*sin(va_c))/3
    )
    vimneg = JuMP.@NLexpression(pm.model,
        (vm_a*sin(va_a) + a2re*vm_b*sin(va_b) + a2im*vm_b*cos(va_b) + are*vm_c*sin(va_c) + aim*vm_c*cos(va_c))/3
    )
    # square of magnitude of U-, |U-|^2
    vmnegsqr = JuMP.@NLexpression(pm.model, vreneg^2+vimneg^2)
    # finally, apply constraint
    JuMP.@NLconstraint(pm.model, vmnegsqr <= vmnegmax^2)
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_positive_sequence(pm::AbstractUnbalancedACPModel, nw::Int, bus_id::Int, vmposmax::Real)
    if !haskey(var(pm, nw_id_default), :vmpossqr)
        var(pm, nw_id_default)[:vmpossqr] = Dict{Int, Any}()
        var(pm, nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [var(pm, nw, :vm, bus_id)[i] for i in 1:3]
    (va_a, va_b, va_c) = [var(pm, nw, :va, bus_id)[i] for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)
    # real and imaginary components of U+
    vrepos = JuMP.@NLexpression(pm.model,
        (vm_a*cos(va_a) + are*vm_b*cos(va_b) - aim*vm_b*sin(va_b) + a2re*vm_c*cos(va_c) - a2im*vm_c*sin(va_c))/3
    )
    vimpos = JuMP.@NLexpression(pm.model,
        (vm_a*sin(va_a) + are*vm_b*sin(va_b) + aim*vm_b*cos(va_b) + a2re*vm_c*sin(va_c) + a2im*vm_c*cos(va_c))/3
    )
    # square of magnitude of U+, |U+|^2
    vmpossqr = JuMP.@NLexpression(pm.model, vrepos^2+vimpos^2)
    # finally, apply constraint
    JuMP.@NLconstraint(pm.model, vmpossqr <= vmposmax^2)
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_zero_sequence(pm::AbstractUnbalancedACPModel, nw::Int, bus_id::Int, vmzeromax::Real)
    if !haskey(var(pm, nw_id_default), :vmpossqr)
        var(pm, nw_id_default)[:vmpossqr] = Dict{Int, Any}()
        var(pm, nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [var(pm, nw, :vm, bus_id)[i] for i in 1:3]
    (va_a, va_b, va_c) = [var(pm, nw, :va, bus_id)[i] for i in 1:3]
    # real and imaginary components of U+
    vrezero = JuMP.@NLexpression(pm.model,
        (vm_a*cos(va_a) + vm_b*cos(va_b) + vm_c*cos(va_c))/3
    )
    vimzero = JuMP.@NLexpression(pm.model,
        (vm_a*sin(va_a) + vm_b*sin(va_b) + vm_c*sin(va_c))/3
    )
    # square of magnitude of U+, |U+|^2
    vmzerosqr = JuMP.@NLexpression(pm.model, vrezero^2+vimzero^2)
    # finally, apply constraint
    JuMP.@NLconstraint(pm.model, vmzerosqr <= vmzeromax^2)
end


"""
We want to express
s_ab = cp.|v_ab|+im.cq.|v_ab|
i_ab = conj(s_ab/v_ab) = |v_ab|.(cq-im.cq)/conj(v_ab) = (1/|v_ab|).(cp-im.cq)*v_ab
idem for i_bc and i_ca
And then
s_a = v_a.conj(i_a) = v_a.conj(i_ab-i_ca)
idem for s_b and s_c
"""
function constraint_mc_load_current_delta(pm::AbstractUnbalancedACPModel, nw::Int, load_id::Int, load_bus_id::Int, cp::Vector, cq::Vector)
    cp_ab, cp_bc, cp_ca = cp
    cq_ab, cq_bc, cq_ca = cq
    vm_a, vm_b, vm_c = var(pm, nw, :vm, load_bus_id)
    va_a, va_b, va_c = var(pm, nw, :va, load_bus_id)
    # v_xy = v_x - v_y
    vre_xy(vm_x, va_x, vm_y, va_y) = JuMP.@NLexpression(pm.model, vm_x*cos(va_x)-vm_y*cos(va_y))
    vim_xy(vm_x, va_x, vm_y, va_y) = JuMP.@NLexpression(pm.model, vm_x*sin(va_x)-vm_y*sin(va_y))
    vre_ab = vre_xy(vm_a, va_a, vm_b, va_b)
    vim_ab = vim_xy(vm_a, va_a, vm_b, va_b)
    vre_bc = vre_xy(vm_b, va_b, vm_c, va_c)
    vim_bc = vim_xy(vm_b, va_b, vm_c, va_c)
    vre_ca = vre_xy(vm_c, va_c, vm_a, va_a)
    vim_ca = vim_xy(vm_c, va_c, vm_a, va_a)
    # i_xy = conj(s_xy/v_xy)
    ire_xy(cp_xy, cq_xy, vre_xy, vim_xy) = JuMP.@NLexpression(pm.model, 1/sqrt(vre_xy^2+vim_xy^2)*(cp_xy*vre_xy+cq_xy*vim_xy))
    iim_xy(cp_xy, cq_xy, vre_xy, vim_xy) = JuMP.@NLexpression(pm.model, 1/sqrt(vre_xy^2+vim_xy^2)*(cp_xy*vim_xy-cq_xy*vre_xy))
    ire_ab = ire_xy(cp_ab, cq_ab, vre_ab, vim_ab)
    iim_ab = iim_xy(cp_ab, cq_ab, vre_ab, vim_ab)
    ire_bc = ire_xy(cp_bc, cq_bc, vre_bc, vim_bc)
    iim_bc = iim_xy(cp_bc, cq_bc, vre_bc, vim_bc)
    ire_ca = ire_xy(cp_ca, cq_ca, vre_ca, vim_ca)
    iim_ca = iim_xy(cp_ca, cq_ca, vre_ca, vim_ca)
    # s_x = v_x*conj(i_xy-i_zx)
    # p_x = vm_x*cos(va_x)*(ire_xy-ire_zx) + vm_x*sin(va_x)*(iim_xy-iim_zx)
    # q_x = vm_x*sin(va_x)*(ire_xy-ire_zx) - vm_x*cos(va_x)*(iim_xy-iim_zx)
    p_x(vm_x, va_x, ire_xy, iim_xy, ire_zx, iim_zx) = JuMP.@NLexpression(pm.model, vm_x*cos(va_x)*(ire_xy-ire_zx) + vm_x*sin(va_x)*(iim_xy-iim_zx))
    q_x(vm_x, va_x, ire_xy, iim_xy, ire_zx, iim_zx) = JuMP.@NLexpression(pm.model, vm_x*sin(va_x)*(ire_xy-ire_zx) - vm_x*cos(va_x)*(iim_xy-iim_zx))
    # s_x = s_x,ref
    var(pm, nw, :pd_bus)[load_id] = [
        p_x(vm_a, va_a, ire_ab, iim_ab, ire_ca, iim_ca),
        p_x(vm_b, va_b, ire_bc, iim_bc, ire_ab, iim_ab),
        p_x(vm_c, va_c, ire_ca, iim_ca, ire_bc, iim_bc)]
    var(pm, nw, :qd_bus)[load_id] = [
        q_x(vm_a, va_a, ire_ab, iim_ab, ire_ca, iim_ca),
        q_x(vm_b, va_b, ire_bc, iim_bc, ire_ab, iim_ab),
        q_x(vm_c, va_c, ire_ca, iim_ca, ire_bc, iim_bc)]
end


""
function constraint_mc_bus_voltage_magnitude_ll(pm::AbstractUnbalancedACPModel, nw::Int, bus_id::Int, vm_ll_min::Vector, vm_ll_max::Vector)
    # 3 conductors asserted in template already
    vm_ln = [var(pm, nw, :vm, bus_id)[i] for i in 1:3]
    va_ln = [var(pm, nw, :va, bus_id)[i] for i in 1:3]
    vr_ll = JuMP.@NLexpression(pm.model, [i in 1:3],
        vm_ln[i]*cos(va_ln[i]) - vm_ln[i%3+1]*cos(va_ln[i%3+1])
    )
    vi_ll = JuMP.@NLexpression(pm.model, [i in 1:3],
        vm_ln[i]*sin(va_ln[i]) - vm_ln[i%3+1]*sin(va_ln[i%3+1])
    )
    for c in 1:3
        # factor of 3 is needed because vm_ll bounds are with respect to the
        # LL base, not the LN base
        if vm_ll_min[c] > 0
            JuMP.@NLconstraint(pm.model, vr_ll[c]^2+vi_ll[c]^2 >= vm_ll_min[c]^2*3)
        end
        if vm_ll_max[c] < Inf
            JuMP.@NLconstraint(pm.model, vr_ll[c]^2+vi_ll[c]^2 <= vm_ll_max[c]^2*3)
        end
    end
end


"bus voltage on/off constraint for load shed problem"
function constraint_mc_bus_voltage_on_off(pm::AbstractUnbalancedACPModel; nw::Int=nw_id_default)
    for (i,bus) in ref(pm, nw, :bus)
        constraint_mc_bus_voltage_magnitude_on_off(pm, i; nw=nw)
    end
end


"`vm[i] == vmref`"
function constraint_mc_voltage_magnitude_only(pm::AbstractUnbalancedACPModel, nw::Int, i::Int, vm_ref::Vector{<:Real})
    bus = ref(pm, nw, :bus, i)
    terminals = bus["terminals"]
    grounded = bus["grounded"]
    vm = var(pm, nw, :vm, i)

    for (idx,t) in [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]
        JuMP.@constraint(pm.model, vm[t] == vm_ref[idx])
    end
end


""
function constraint_mc_storage_current_limit(pm::AbstractUnbalancedACPModel, nw::Int, i::Int, bus_id::Int, connections::Vector{Int}, rating::Vector{<:Real})
    vm = [var(pm, nw, :vm, bus_id)[c] for c in connections]
    ps = [var(pm, nw, :ps, i)[c] for c in connections]
    qs = [var(pm, nw, :qs, i)[c] for c in connections]

    JuMP.@constraint(pm.model, ps.^2 + qs.^2 <= rating.^2 .* vm.^2)
end


""
function constraint_mc_load_power_wye(pm::AbstractUnbalancedACPModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vm = var(pm, nw, :vm, bus_id)
    va = var(pm, nw, :va, bus_id)

    # if constant power load
    if ref(pm, nw, :load, id, "model") == POWER
        pd_bus = a
        qd_bus = b
    else
        pd_bus = Vector{JuMP.NonlinearExpression}([])
        qd_bus = Vector{JuMP.NonlinearExpression}([])

        for (idx, c) in enumerate(connections)
            crd = JuMP.@NLexpression(pm.model,
                 a[idx]*vm[c]*cos(va[c])*(vm[c]^2)^(alpha[idx]/2-1)
                +b[idx]*vm[c]*sin(va[c])*(vm[c]^2)^( beta[idx]/2-1)
            )
            cid = JuMP.@NLexpression(pm.model,
                 a[idx]*vm[c]*sin(va[c])*(vm[c]^2)^(alpha[idx]/2-1)
                -b[idx]*vm[c]*cos(va[c])*(vm[c]^2)^( beta[idx]/2-1)
            )

            push!(pd_bus, JuMP.@NLexpression(pm.model,  vm[c]*cos(va[c])*crd+vm[c]*sin(va[c])*cid))
            push!(qd_bus, JuMP.@NLexpression(pm.model, -vm[c]*cos(va[c])*cid+vm[c]*sin(va[c])*crd))
        end
    end

    pd_bus = JuMP.Containers.DenseAxisArray(pd_bus, connections)
    qd_bus = JuMP.Containers.DenseAxisArray(qd_bus, connections)

    var(pm, nw, :pd_bus)[id] = pd_bus
    var(pm, nw, :qd_bus)[id] = qd_bus

    if report
        sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        pd = Vector{JuMP.NonlinearExpression}([])
        qd = Vector{JuMP.NonlinearExpression}([])

        for (idx,c) in enumerate(connections)
            push!(pd, JuMP.@NLexpression(pm.model, a[idx]*vm[c]^alpha[idx] ))
            push!(qd, JuMP.@NLexpression(pm.model, b[idx]*vm[c]^beta[idx]  ))
        end
        sol(pm, nw, :load, id)[:pd] = JuMP.Containers.DenseAxisArray(pd, connections)
        sol(pm, nw, :load, id)[:qd] = JuMP.Containers.DenseAxisArray(qd, connections)
    end
end


""
function constraint_mc_load_power_delta(pm::AbstractUnbalancedACPModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vm = var(pm, nw, :vm, bus_id)
    va = var(pm, nw, :va, bus_id)

    nph = length(a)

    prev = Dict(c=>connections[(idx+nph-2)%nph+1] for (idx,c) in enumerate(connections))
    next = Dict(c=>connections[idx%nph+1] for (idx,c) in enumerate(connections))

    vrd = Dict()
    vid = Dict()
    for (idx, c) in enumerate(connections)
        vrd[c] = JuMP.@NLexpression(pm.model, vm[c]*cos(va[c])-vm[next[c]]*cos(va[next[c]]))
        vid[c] = JuMP.@NLexpression(pm.model, vm[c]*sin(va[c])-vm[next[c]]*sin(va[next[c]]))
    end

    crd = Dict()
    cid = Dict()
    for (idx, c) in enumerate(connections)
        crd[c] = JuMP.@NLexpression(pm.model,
             a[idx]*vrd[c]*(vrd[c]^2+vid[c]^2)^(alpha[idx]/2-1)
            +b[idx]*vid[c]*(vrd[c]^2+vid[c]^2)^(beta[idx]/2 -1)
        )
        cid[c] = JuMP.@NLexpression(pm.model,
             a[idx]*vid[c]*(vrd[c]^2+vid[c]^2)^(alpha[idx]/2-1)
            -b[idx]*vrd[c]*(vrd[c]^2+vid[c]^2)^(beta[idx]/2 -1)
        )
    end

    crd_bus = Dict()
    cid_bus = Dict()
    for (idx, c) in enumerate(connections)
        crd_bus[c] = JuMP.@NLexpression(pm.model, crd[c]-crd[prev[c]])
        cid_bus[c] = JuMP.@NLexpression(pm.model, cid[c]-cid[prev[c]])
    end

    pd_bus = Vector{JuMP.NonlinearExpression}([])
    qd_bus = Vector{JuMP.NonlinearExpression}([])
    for (idx,c) in enumerate(connections)
        push!(pd_bus, JuMP.@NLexpression(pm.model,  vm[c]*cos(va[c])*crd_bus[c]+vm[c]*sin(va[c])*cid_bus[c]))
        push!(qd_bus, JuMP.@NLexpression(pm.model, -vm[c]*cos(va[c])*cid_bus[c]+vm[c]*sin(va[c])*crd_bus[c]))
    end

    pd_bus = JuMP.Containers.DenseAxisArray(pd_bus, connections)
    qd_bus = JuMP.Containers.DenseAxisArray(qd_bus, connections)

    var(pm, nw, :pd_bus)[id] = pd_bus
    var(pm, nw, :qd_bus)[id] = qd_bus

    if report
        sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        pd = []
        qd = []
        for (idx,c) in enumerate(connections)
            push!(pd, JuMP.@NLexpression(pm.model, a[idx]*(vrd[c]^2+vid[c]^2)^(alpha[idx]/2) ))
            push!(qd, JuMP.@NLexpression(pm.model, b[idx]*(vrd[c]^2+vid[c]^2)^(beta[idx]/2)  ))
        end
        sol(pm, nw, :load, id)[:pd] = pd
        sol(pm, nw, :load, id)[:qd] = qd
    end
end


""
function constraint_mc_generator_power_delta(pm::AbstractUnbalancedACPModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vm = var(pm, nw, :vm, bus_id)
    va = var(pm, nw, :va, bus_id)
    pg = var(pm, nw, :pg, id)
    qg = var(pm, nw, :qg, id)

    nph = length(connections)
    prev = Dict(c=>connections[(idx+nph-2)%nph+1] for (idx,c) in enumerate(connections))
    next = Dict(c=>connections[idx%nph+1] for (idx,c) in enumerate(connections))

    vrg = Dict()
    vig = Dict()
    for c in connections
        vrg[c] = JuMP.@NLexpression(pm.model, vm[c]*cos(va[c])-vm[next[c]]*cos(va[next[c]]))
        vig[c] = JuMP.@NLexpression(pm.model, vm[c]*sin(va[c])-vm[next[c]]*sin(va[next[c]]))
    end

    crg = Dict()
    cig = Dict()
    for c in connections
        crg[c] = JuMP.@NLexpression(pm.model, (pg[c]*vrg[c]+qg[c]*vig[c])/(vrg[c]^2+vig[c]^2) )
        cig[c] = JuMP.@NLexpression(pm.model, (pg[c]*vig[c]-qg[c]*vrg[c])/(vrg[c]^2+vig[c]^2) )
    end

    crg_bus = Dict()
    cig_bus = Dict()
    for c in connections
        crg_bus[c] = JuMP.@NLexpression(pm.model, crg[c]-crg[prev[c]])
        cig_bus[c] = JuMP.@NLexpression(pm.model, cig[c]-cig[prev[c]])
    end

    pg_bus = Vector{JuMP.NonlinearExpression}([])
    qg_bus = Vector{JuMP.NonlinearExpression}([])
    for c in connections
        push!(pg_bus, JuMP.@NLexpression(pm.model,  vm[c]*cos(va[c])*crg_bus[c]+vm[c]*sin(va[c])*cig_bus[c]))
        push!(qg_bus, JuMP.@NLexpression(pm.model, -vm[c]*cos(va[c])*cig_bus[c]+vm[c]*sin(va[c])*crg_bus[c]))
    end
    pd_bus = JuMP.Containers.DenseAxisArray(pg_bus, connections)
    qd_bus = JuMP.Containers.DenseAxisArray(qg_bus, connections)

    var(pm, nw, :pg_bus)[id] = pg_bus
    var(pm, nw, :qg_bus)[id] = qg_bus

    if report
        sol(pm, nw, :gen, id)[:pg_bus] = pg_bus
        sol(pm, nw, :gen, id)[:qg_bus] = qg_bus
    end
end


""
function constraint_mc_storage_losses(pm::AbstractUnbalancedACPModel, n::Int, i::Int, bus::Int, connections::Vector{Int}, r::Real, x::Real, p_loss::Real, q_loss::Real)
    vm = var(pm, n, :vm, bus)
    ps = var(pm, n, :ps, i)
    qs = var(pm, n, :qs, i)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    qsc = var(pm, n, :qsc, i)

    JuMP.@NLconstraint(pm.model,
        sum(ps[c] for c in connections) + (sd - sc)
        ==
        p_loss + r * sum((ps[c]^2 + qs[c]^2)/vm[c]^2 for c in connections)
    )

    JuMP.@NLconstraint(pm.model,
        sum(qs[c] for c in connections)
        ==
        qsc + q_loss + x * sum((ps[c]^2 + qs[c]^2)/vm[c]^2 for c in connections)
    )
end


@doc raw"""
    constraint_mc_ampacity_from(pm::AbstractUnbalancedACPModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches from-side

math```
p_{fr}^2 + q_{fr}^2 \leq vm_{fr}^2 i_{max}^2
```
"""
function constraint_mc_ampacity_from(pm::AbstractUnbalancedACPModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_fr = [var(pm, nw, :p, f_idx)[c] for c in f_connections]
    q_fr = [var(pm, nw, :q, f_idx)[c] for c in f_connections]
    vm_fr = [var(pm, nw, :vm, f_idx[2])[c] for c in f_connections]

    con(pm, nw, :mu_cm_branch)[f_idx] = mu_cm_fr = [JuMP.@constraint(pm.model, p_fr[idx]^2 + q_fr[idx]^2 <= vm_fr[idx]^2 * c_rating[idx]^2) for idx in findall(c_rating .< Inf)]

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, f_idx[1])[:mu_cm_fr] = mu_cm_fr
    end

    nothing
end


@doc raw"""
    constraint_mc_ampacity_to(pm::AbstractUnbalancedACPModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches to-side

math```
p_{to}^2 + q_{to}^2 \leq vm_{to}^2 i_{max}^2
```
"""
function constraint_mc_ampacity_to(pm::AbstractUnbalancedACPModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_to = [var(pm, nw, :p, t_idx)[c] for c in t_connections]
    q_to = [var(pm, nw, :q, t_idx)[c] for c in t_connections]
    vm_to = [var(pm, nw, :vm, t_idx[2])[c] for c in t_connections]

    con(pm, nw, :mu_cm_branch)[t_idx] = mu_cm_to = [JuMP.@constraint(pm.model, p_to[idx]^2 + q_to[idx]^2 <= vm_to[idx]^2 * c_rating[idx]^2) for idx in findall(c_rating .< Inf)]

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, t_idx[1])[:mu_cm_to] = mu_cm_to
    end

    nothing
end


@doc raw"""
    constraint_mc_switch_ampacity(pm::AbstractUnbalancedACPModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on switches

math```
p_{fr}^2 + q_{fr}^2 \leq vm_{fr}^2 i_{max}^2
```
"""
function constraint_mc_switch_ampacity(pm::AbstractUnbalancedACPModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    psw_fr = [var(pm, nw, :psw, f_idx)[c] for c in f_connections]
    qsw_fr = [var(pm, nw, :qsw, f_idx)[c] for c in f_connections]
    vm_fr = [var(pm, nw, :vm, f_idx[2])[c] for c in f_connections]

    con(pm, nw, :mu_cm_switch)[f_idx] = mu_cm_fr = [JuMP.@constraint(pm.model, psw_fr[idx]^2 + qsw_fr[idx]^2 .<= vm_fr[idx]^2 * c_rating[idx]^2) for idx in findall(c_rating .< Inf)]

    if _IM.report_duals(pm)
        sol(pm, nw, :switch, f_idx[1])[:mu_cm_fr] = mu_cm_fr
    end

    nothing
end
