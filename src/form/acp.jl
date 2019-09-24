""
function variable_mc_voltage(pm::_PMs.AbstractACPModel; nw=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_voltage(pm, cnd=c; nw=nw, kwargs...)
    end
    # This is needed for delta loads, where division occurs by the difference
    # of voltage phasors. If the voltage phasors at one bus are initialized
    # in the same point, this would lead to division by zero.
    ncnd = length(_PMs.conductor_ids(pm))
    theta = [_wrap_to_pi(2 * pi / ncnd * (1-c)) for c in 1:ncnd]
    vm = 1
    for id in _PMs.ids(pm, :bus)
        busref = _PMs.ref(pm, nw, :bus, id)
        if !haskey(busref, "va_start")
        # if it has this key, it was set at PM level
            for c in 1:ncnd
                JuMP.set_start_value(_PMs.var(pm, nw, c, :va, id), theta[c])
            end
        end
    end
end


""
function variable_mc_bus_voltage_on_off(pm::_PMs.AbstractACPModel; kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_voltage_angle(pm; cnd=c, kwargs...)
    end
    variable_mc_voltage_magnitude_on_off(pm; kwargs...)

    nw = get(kwargs, :nw, pm.cnw)

    ncnd = length(_PMs.conductor_ids(pm))
    theta = [_wrap_to_pi(2 * pi / ncnd * (1-c)) for c in 1:ncnd]
    vm = 1
    for id in _PMs.ids(pm, :bus)
        busref = _PMs.ref(pm, nw, :bus, id)
        if !haskey(busref, "va_start")
        # if it has this key, it was set at PM level
            for c in 1:ncnd
                JuMP.set_start_value(_PMs.var(pm, nw, c, :va, id), theta[c])
            end
        end
    end
end


""
function constraint_mc_power_balance_slack(pm::_PMs.AbstractACPModel, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm   = _PMs.var(pm, nw, c, :vm, i)
    p    = get(_PMs.var(pm, nw, c),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(_PMs.var(pm, nw, c),    :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(_PMs.var(pm, nw, c),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(_PMs.var(pm, nw, c),   :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PMs.var(pm, nw, c),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PMs.var(pm, nw, c),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(_PMs.var(pm, nw, c),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(_PMs.var(pm, nw, c),  :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(_PMs.var(pm, nw, c),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(_PMs.var(pm, nw, c),   :qt, Dict()); _PMs._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    p_slack = _PMs.var(pm, nw, c, :p_slack, i)
    q_slack = _PMs.var(pm, nw, c, :q_slack, i)

    _PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        - sum(gs for gs in values(bus_gs))*vm^2
        + p_slack
    )
    _PMs.con(pm, nw, c, :kcl_q)[i] = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(qt[a_trans] for a_trans in bus_arcs_trans)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd for qd in values(bus_qd))
        + sum(bs for bs in values(bus_bs))*vm^2
        + q_slack
    )
end


""
function constraint_mc_power_balance_shed(pm::_PMs.AbstractACPModel, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm       = _PMs.var(pm, nw, c, :vm, i)
    p        = get(_PMs.var(pm, nw, c),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q        = get(_PMs.var(pm, nw, c),    :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg       = get(_PMs.var(pm, nw, c),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg       = get(_PMs.var(pm, nw, c),   :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps       = get(_PMs.var(pm, nw, c),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs       = get(_PMs.var(pm, nw, c),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw      = get(_PMs.var(pm, nw, c),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw      = get(_PMs.var(pm, nw, c),  :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt       = get(_PMs.var(pm, nw, c),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt       = get(_PMs.var(pm, nw, c),   :qt, Dict()); _PMs._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    z_demand = _PMs.var(pm, nw, :z_demand)
    z_shunt  = _PMs.var(pm, nw, :z_shunt)

    _PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd*z_demand[n] for (n,pd) in bus_pd)
        - sum(gs*z_shunt[n] for (n,gs) in bus_gs)*vm^2
    )
    _PMs.con(pm, nw, c, :kcl_q)[i] = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(qt[a_trans] for a_trans in bus_arcs_trans)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd*z_demand[n] for (n,qd) in bus_qd)
        + sum(bs*z_shunt[n] for (n,bs) in bus_bs)*vm^2
    )
end


""
function constraint_mc_power_balance(pm::_PMs.AbstractACPModel, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm   = _PMs.var(pm, nw, c, :vm, i)
    p    = get(_PMs.var(pm, nw, c),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(_PMs.var(pm, nw, c),    :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(_PMs.var(pm, nw, c),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(_PMs.var(pm, nw, c),   :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PMs.var(pm, nw, c),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PMs.var(pm, nw, c),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(_PMs.var(pm, nw, c),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(_PMs.var(pm, nw, c),  :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(_PMs.var(pm, nw, c),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(_PMs.var(pm, nw, c),   :qt, Dict()); _PMs._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")

    _PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        - sum(gs for gs in values(bus_gs))*vm^2
    )
    _PMs.con(pm, nw, c, :kcl_q)[i] = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(qt[a_trans] for a_trans in bus_arcs_trans)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd for qd in values(bus_qd))
        + sum(bs for bs in values(bus_bs))*vm^2
    )
end


""
function constraint_mc_power_balance_load(pm::_PMs.AbstractACPModel, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_gs, bus_bs)
    vm   = _PMs.var(pm, nw, c, :vm, i)
    p    = get(_PMs.var(pm, nw, c),   :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(_PMs.var(pm, nw, c),   :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(_PMs.var(pm, nw, c),  :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(_PMs.var(pm, nw, c),  :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PMs.var(pm, nw, c),  :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PMs.var(pm, nw, c),  :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(_PMs.var(pm, nw, c), :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(_PMs.var(pm, nw, c), :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(_PMs.var(pm, nw, c),  :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(_PMs.var(pm, nw, c),  :qt, Dict()); _PMs._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    pd   = get(_PMs.var(pm, nw, c),  :pd, Dict()); _PMs._check_var_keys(pd, bus_loads, "active power", "load")
    qd   = get(_PMs.var(pm, nw, c),  :qd, Dict()); _PMs._check_var_keys(pd, bus_loads, "reactive power", "load")

    _PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@NLconstraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd[l] for l in bus_loads)
        - sum(gs for gs in values(bus_gs))*vm^2
    )
    _PMs.con(pm, nw, c, :kcl_q)[i] = JuMP.@NLconstraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(qt[a_trans] for a_trans in bus_arcs_trans)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd[l] for l in bus_loads)
        + sum(bs for bs in values(bus_bs))*vm^2
    )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p_fr ==     g[c,c] * vm_fr[c]^2 +
            sum( g[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                 b[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in _PMs.conductor_ids(pm) if d != c) +
            sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in _PMs.conductor_ids(pm))
            + g_fr[c,c] * vm_fr[c]^2 +
            sum( g_fr[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                 b_fr[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in _PMs.conductor_ids(pm) if d != c)
            )
q_fr == -b[c,c] *vm_fr[c]^2 -
            sum( b[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                 g[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in _PMs.conductor_ids(pm) if d != c) -
            sum(-b[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                 g[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in _PMs.conductor_ids(pm))
            -b_fr[c,c] *vm_fr[c]^2 -
            sum( b_fr[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                 g_fr[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in _PMs.conductor_ids(pm) if d != c)
            )
```
"""
function constraint_mc_ohms_yt_from(pm::_PMs.AbstractACPModel, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr  = _PMs.var(pm, n, c,  :p, f_idx)
    q_fr  = _PMs.var(pm, n, c,  :q, f_idx)
    vm_fr = [_PMs.var(pm, n, d, :vm, f_bus) for d in _PMs.conductor_ids(pm)]
    vm_to = [_PMs.var(pm, n, d, :vm, t_bus) for d in _PMs.conductor_ids(pm)]
    va_fr = [_PMs.var(pm, n, d, :va, f_bus) for d in _PMs.conductor_ids(pm)]
    va_to = [_PMs.var(pm, n, d, :va, t_bus) for d in _PMs.conductor_ids(pm)]

    JuMP.@NLconstraint(pm.model, p_fr ==
                                    (g[c,c]+g_fr[c,c])*vm_fr[c]^2
                                    +sum( (g[c,d]+g_fr[c,d]) * vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d])
                                         +(b[c,d]+b_fr[c,d]) * vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d])
                                         for d in _PMs.conductor_ids(pm) if d != c)
                                    +sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d])
                                         -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d])
                                         for d in _PMs.conductor_ids(pm))
                                )
    JuMP.@NLconstraint(pm.model, q_fr == -(b[c,c]+b_fr[c,c])*vm_fr[c]^2
                                    -sum( (b[c,d]+b_fr[c,d])*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d])
                                         -(g[c,d]+g_fr[c,d])*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d])
                                         for d in _PMs.conductor_ids(pm) if d != c)
                                    -sum(-b[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d])
                                         +g[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d])
                                         for d in _PMs.conductor_ids(pm))
                                )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_mc_ohms_yt_to(pm::_PMs.AbstractACPModel, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    constraint_mc_ohms_yt_from(pm, n, c, t_bus, f_bus, t_idx, f_idx, g, b, g_to, b_to, tr, ti, tm)
end


#TODO: extend for matrix shunts
"""
```
p[f_idx] == z*(g/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
q[f_idx] == z*(-(b+c/2)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
```
"""
function constraint_mc_ohms_yt_from_on_off(pm::_PMs.AbstractACPModel, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
    p_fr  = _PMs.var(pm, n, c,  :p, f_idx)
    q_fr  = _PMs.var(pm, n, c,  :q, f_idx)
    vm_fr = [_PMs.var(pm, n, d, :vm, f_bus) for d in _PMs.conductor_ids(pm)]
    vm_to = [_PMs.var(pm, n, d, :vm, t_bus) for d in _PMs.conductor_ids(pm)]
    va_fr = [_PMs.var(pm, n, d, :va, f_bus) for d in _PMs.conductor_ids(pm)]
    va_to = [_PMs.var(pm, n, d, :va, t_bus) for d in _PMs.conductor_ids(pm)]
    z = _PMs.var(pm, n, c, :branch_z, i)

    JuMP.@NLconstraint(pm.model, p_fr == z*((g_fr[c]+g[c,c]) * vm_fr[c]^2 +
                                        sum( g[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                                             b[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in _PMs.conductor_ids(pm) if d != c) +
                                        sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                            -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in _PMs.conductor_ids(pm))) )
    JuMP.@NLconstraint(pm.model, q_fr == z*(-(b_fr[c]+b[c,c]) *vm_fr[c]^2 -
                                        sum( b[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                                             g[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in _PMs.conductor_ids(pm) if d != c) -
                                        sum(-b[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                             g[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in _PMs.conductor_ids(pm))) )
end


#TODO: extend for matrix shunts
"""
```
p[t_idx] == z*(g*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
q[t_idx] == z*(-(b+c/2)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
```
"""
function constraint_mc_ohms_yt_to_on_off(pm::_PMs.AbstractACPModel, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
    p_to  = _PMs.var(pm, n, c,  :p, t_idx)
    q_to  = _PMs.var(pm, n, c,  :q, t_idx)
    vm_fr = [_PMs.var(pm, n, d, :vm, f_bus) for d in _PMs.conductor_ids(pm)]
    vm_to = [_PMs.var(pm, n, d, :vm, t_bus) for d in _PMs.conductor_ids(pm)]
    va_fr = [_PMs.var(pm, n, d, :va, f_bus) for d in _PMs.conductor_ids(pm)]
    va_to = [_PMs.var(pm, n, d, :va, t_bus) for d in _PMs.conductor_ids(pm)]
    z = _PMs.var(pm, n, c, :branch_z, i)

    JuMP.@NLconstraint(pm.model, p_to == z*((g_to[c]+g[c,c])*vm_to[c]^2 +
                                        sum( g[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) +
                                             b[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in _PMs.conductor_ids(pm) if d != c) +
                                        sum(-g[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                            -b[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in _PMs.conductor_ids(pm))) )
    JuMP.@NLconstraint(pm.model, q_to == z*(-(b_to[c]+b[c,c])*vm_to[c]^2 -
                                        sum( b[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) -
                                             g[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in _PMs.conductor_ids(pm) if d != c) -
                                        sum(-b[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                             g[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in _PMs.conductor_ids(pm))) )
end


"Links the voltage at both windings of a fixed tap transformer"
function constraint_mc_transformer_voltage(pm::_PMs.AbstractACPModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, tm::_PMs.MultiConductorVector, Tv_fr, Tv_im, Cv_to)
    ncnd  = 3
    # from side
    vm_fr = [_PMs.var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [_PMs.var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    # to side
    vm_to = [_PMs.var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
    va_to = [_PMs.var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
    # the intermediate bus voltage is saved as an expression
    # vm_im[c] = vm_to[c]*tm[c]*Cv_to
    # va_im = va_to
    for n in 1:size(Tv_fr)[1]
        JuMP.@NLconstraint(pm.model,
              sum(Tv_fr[n,c]*vm_fr[c]*cos(va_fr[c]) for c in 1:ncnd)
            ==sum(Tv_im[n,c]*(vm_to[c]*tm[c]*Cv_to)*cos(va_to[c]) for c in 1:ncnd)
        )
        JuMP.@NLconstraint(pm.model,
              sum(Tv_fr[n,c]*vm_fr[c]*sin(va_fr[c]) for c in 1:ncnd)
            ==sum(Tv_im[n,c]*(vm_to[c]*tm[c]*Cv_to)*sin(va_to[c]) for c in 1:ncnd)
        )
    end
end


"Links the power flowing into both windings of a fixed tap transformer"
function constraint_mc_transformer_flow(pm::_PMs.AbstractACPModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, tm::_PMs.MultiConductorVector, Ti_fr, Ti_im, Cv_to)
    ncnd  = 3
    # from side variables
    vm_fr = [_PMs.var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [_PMs.var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    p_fr = [_PMs.var(pm, nw, c, :pt, f_idx) for c in 1:ncnd]
    q_fr = [_PMs.var(pm, nw, c, :qt, f_idx) for c in 1:ncnd]
    # to side
    vm_to = [_PMs.var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
    va_to = [_PMs.var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
    p_to = [_PMs.var(pm, nw, c, :pt, t_idx) for c in 1:ncnd]
    q_to = [_PMs.var(pm, nw, c, :qt, t_idx) for c in 1:ncnd]
    # the intermediate bus voltage is saved as an expression
    # vm_im[c] = vm_to[c]*tm[c]*Cv_to
    # va_im = va_to
    for n in 1:size(Ti_fr)[1]
        # i_fr_re[c] = 1/vm_fr[c]*(p_fr[c]*cos(va_fr[c])+q_fr[c]*sin(va_fr[c]))
        # i_fr_im[c] = 1/vm_fr[c]*(p_fr[c]*sin(va_fr[c])-q_fr[c]*cos(va_fr[c]))
        # i_to_re[c] = 1/vm_to[c]*(p_to[c]*cos(va_im[c])+q_to[c]*sin(va_im[c]))
        # i_to_im[c] = 1/vm_to[c]*(p_to[c]*sin(va_im[c])-q_to[c]*cos(va_im[c]))
        JuMP.@NLconstraint(pm.model,
              sum(Ti_fr[n,c]*
                    1/vm_fr[c]*(p_fr[c]*cos(va_fr[c])+q_fr[c]*sin(va_fr[c])) # i_fr_re[c]
              for c in 1:ncnd)
            + sum(Ti_im[n,c]*
                    1/(vm_to[c]*tm[c]*Cv_to)*(p_to[c]*cos(va_to[c])+q_to[c]*sin(va_to[c])) # i_to_re[c]
              for c in 1:ncnd)
            == 0
        )
        JuMP.@NLconstraint(pm.model,
              sum(Ti_fr[n,c]*
                    1/vm_fr[c]*(p_fr[c]*sin(va_fr[c])-q_fr[c]*cos(va_fr[c])) # i_fr_im[c]
              for c in 1:ncnd)
            + sum(Ti_im[n,c]*
                    1/(vm_to[c]*tm[c]*Cv_to)*(p_to[c]*sin(va_to[c])-q_to[c]*cos(va_to[c])) # i_to_im[c]
              for c in 1:ncnd)
            == 0
        )
    end
end


"Links the voltage at both windings of a variable tap transformer."
function constraint_mc_oltc_voltage(pm::_PMs.AbstractACPModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, Tv_fr, Tv_im, Cv_to)
    ncnd  = 3
    # from side
    vm_fr = [_PMs.var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [_PMs.var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    # tap
    tap = [_PMs.var(pm, nw, c, :tap)[i] for c in 1:ncnd]
    # to side
    vm_to = [_PMs.var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
    va_to = [_PMs.var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
    # the intermediate bus voltage is saved as an expression
    # vm_im[c] = vm_to[c]*tap[c]*Cv_to
    # va_im = va_to
    for n in 1:size(Tv_fr)[1]
        JuMP.@NLconstraint(pm.model,
              sum(Tv_fr[n,c]*vm_fr[c]*cos(va_fr[c]) for c in 1:ncnd)
            ==sum(Tv_im[n,c]*(vm_to[c]*tap[c]*Cv_to)*cos(va_to[c]) for c in 1:ncnd)
        )
        JuMP.@NLconstraint(pm.model,
              sum(Tv_fr[n,c]*vm_fr[c]*sin(va_fr[c]) for c in 1:ncnd)
            ==sum(Tv_im[n,c]*(vm_to[c]*tap[c]*Cv_to)*sin(va_to[c]) for c in 1:ncnd)
        )
    end
end


"Links the power flowing into both windings of a variable tap transformer."
function constraint_mc_oltc_flow(pm::_PMs.AbstractACPModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, Ti_fr, Ti_im, Cv_to)
    ncnd  = 3
    # from side variables
    vm_fr = [_PMs.var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [_PMs.var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    p_fr = [_PMs.var(pm, nw, c, :pt, f_idx) for c in 1:ncnd]
    q_fr = [_PMs.var(pm, nw, c, :qt, f_idx) for c in 1:ncnd]
    # tap
    tap = [_PMs.var(pm, nw, c, :tap, i) for c in 1:ncnd]
    # to side
    vm_to = [_PMs.var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
    va_to = [_PMs.var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
    p_to = [_PMs.var(pm, nw, c, :pt, t_idx) for c in 1:ncnd]
    q_to = [_PMs.var(pm, nw, c, :qt, t_idx) for c in 1:ncnd]
    # the intermediate bus voltage is saved as an expression
    # vm_im[c] = vm_to[c]*tap[c]*Cv_to
    # va_im = va_to
    for n in 1:size(Ti_fr)[1]
        # i_fr_re[c] = 1/vm_fr[c]*(p_fr[c]*cos(va_fr[c])+q_fr[c]*sin(va_fr[c]))
        # i_fr_im[c] = 1/vm_fr[c]*(p_fr[c]*sin(va_fr[c])-q_fr[c]*cos(va_fr[c]))
        # i_to_re[c] = 1/vm_to[c]*(p_to[c]*cos(va_im[c])+q_to[c]*sin(va_im[c]))
        # i_to_im[c] = 1/vm_to[c]*(p_to[c]*sin(va_im[c])-q_to[c]*cos(va_im[c]))
        JuMP.@NLconstraint(pm.model,
              sum(Ti_fr[n,c]*
                    1/vm_fr[c]*(p_fr[c]*cos(va_fr[c])+q_fr[c]*sin(va_fr[c])) # i_fr_re[c]
              for c in 1:ncnd)
            + sum(Ti_im[n,c]*
                    1/(vm_to[c]*tap[c]*Cv_to)*(p_to[c]*cos(va_to[c])+q_to[c]*sin(va_to[c])) # i_to_re[c]
              for c in 1:ncnd)
            == 0
        )
        JuMP.@NLconstraint(pm.model,
              sum(Ti_fr[n,c]*
                    1/vm_fr[c]*(p_fr[c]*sin(va_fr[c])-q_fr[c]*cos(va_fr[c])) # i_fr_im[c]
              for c in 1:ncnd)
            + sum(Ti_im[n,c]*
                    1/(vm_to[c]*tap[c]*Cv_to)*(p_to[c]*sin(va_to[c])-q_to[c]*cos(va_to[c])) # i_to_im[c]
              for c in 1:ncnd)
            == 0
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
function constraint_mc_vm_vuf(pm::_PMs.AbstractACPModel, nw::Int, bus_id::Int, vufmax::Float64)
    if !haskey(_PMs.var(pm, pm.cnw), :vmpossqr)
        _PMs.var(pm, pm.cnw)[:vmpossqr] = Dict{Int, Any}()
        _PMs.var(pm, pm.cnw)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [_PMs.var(pm, nw, i, :vm, bus_id) for i in 1:3]
    (va_a, va_b, va_c) = [_PMs.var(pm, nw, i, :va, bus_id) for i in 1:3]
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
    #var(pm, pm.cnw, :vmpossqr)[bus_id] = vmpossqr
    #var(pm, pm.cnw, :vmnegsqr)[bus_id] = vmnegsqr
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_vm_neg_seq(pm::_PMs.AbstractACPModel, nw::Int, bus_id::Int, vmnegmax::Float64)
    if !haskey(_PMs.var(pm, pm.cnw), :vmpossqr)
        _PMs.var(pm, pm.cnw)[:vmpossqr] = Dict{Int, Any}()
        _PMs.var(pm, pm.cnw)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [_PMs.var(pm, nw, i, :vm, bus_id) for i in 1:3]
    (va_a, va_b, va_c) = [_PMs.var(pm, nw, i, :va, bus_id) for i in 1:3]
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


"Links the power flowing into both windings of a variable tap transformer."
function constraint_mc_transformer_flow_var(pm::_PMs.AbstractPowerModel, i::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, Ti_fr, Ti_im; nw::Int=pm.cnw)
    # for ac formulation, indentical to fixed tap
    constraint_mc_transformer_flow(pm, i, f_bus, t_bus, f_idx, t_idx, Ti_fr, Ti_im)
end


""
function constraint_load_power_wye(pm::_PMs.AbstractACPModel, nw::Int, cnd::Int, load_id::Int, pd::Real, qd::Real)
    _PMs.var(pm, nw, cnd, :pd)[load_id] = pd
    _PMs.var(pm, nw, cnd, :qd)[load_id] = qd
end


""
function constraint_load_current_wye(pm::_PMs.AbstractACPModel, nw::Int, cnd::Int, load_id::Int, load_bus_id::Int, scale_p::Real, scale_q::Real)
    vm = _PMs.var(pm, nw, cnd, :vm, load_bus_id)
    # this has to be a NLexpression and not an expression;
    # might be mixed in KCL with NLexpression, so has to also be NLexpression
    _PMs.var(pm, nw, cnd, :pd)[load_id] = JuMP.@NLexpression(pm.model, scale_p*vm)
    _PMs.var(pm, nw, cnd, :qd)[load_id] = JuMP.@NLexpression(pm.model, scale_q*vm)
end


""
function constraint_load_impedance_wye(pm::_PMs.AbstractACPModel, nw::Int, cnd::Int, load_id::Int, load_bus_id::Int, cp::Real, cq::Real)
    vm = _PMs.var(pm, nw, cnd, :vm, load_bus_id)
    _PMs.var(pm, nw, cnd, :pd)[load_id] = JuMP.@NLexpression(pm.model, cp*vm^2)
    _PMs.var(pm, nw, cnd, :qd)[load_id] = JuMP.@NLexpression(pm.model, cq*vm^2)
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_vm_pos_seq(pm::_PMs.AbstractACPModel, nw::Int, bus_id::Int, vmposmax::Float64)
    if !haskey(_PMs.var(pm, pm.cnw), :vmpossqr)
        _PMs.var(pm, pm.cnw)[:vmpossqr] = Dict{Int, Any}()
        _PMs.var(pm, pm.cnw)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [_PMs.var(pm, nw, i, :vm, bus_id) for i in 1:3]
    (va_a, va_b, va_c) = [_PMs.var(pm, nw, i, :va, bus_id) for i in 1:3]
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
For a delta load, sd = (s_ab, s_bc, s_ca), but we want to fix s = (s_a, s_b, s_c)
s is a non-linear transform of v and sd, s=f(v,sd)
s_a = v_a*conj(s_ab/(v_a-v_b) - s_ca/(v_c-v_a))
s_b = v_b*conj(s_ab/(v_a-v_b) - s_ca/(v_c-v_a))
s_c = v_c*conj(s_ab/(v_a-v_b) - s_ca/(v_c-v_a))
"""
function constraint_mc_load_power_delta(pm::_PMs.AbstractACPModel, nw::Int, load_id::Int, load_bus_id::Int, pd::_PMs.MultiConductorVector, qd::_PMs.MultiConductorVector)
    p_ab, p_bc, p_ca = pd
    q_ab, q_bc, q_ca = qd
    vm_a, vm_b, vm_c = [_PMs.var(pm, nw, c, :vm, load_bus_id) for c in 1:3]
    va_a, va_b, va_c = [_PMs.var(pm, nw, c, :va, load_bus_id) for c in 1:3]
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
    ire_xy(p_xy, q_xy, vre_xy, vim_xy) = JuMP.@NLexpression(pm.model, (p_xy*vre_xy+q_xy*vim_xy)/(vre_xy^2+vim_xy^2))
    iim_xy(p_xy, q_xy, vre_xy, vim_xy) = JuMP.@NLexpression(pm.model, (p_xy*vim_xy-q_xy*vre_xy)/(vre_xy^2+vim_xy^2))
    ire_ab = ire_xy(p_ab, q_ab, vre_ab, vim_ab)
    iim_ab = iim_xy(p_ab, q_ab, vre_ab, vim_ab)
    ire_bc = ire_xy(p_bc, q_bc, vre_bc, vim_bc)
    iim_bc = iim_xy(p_bc, q_bc, vre_bc, vim_bc)
    ire_ca = ire_xy(p_ca, q_ca, vre_ca, vim_ca)
    iim_ca = iim_xy(p_ca, q_ca, vre_ca, vim_ca)
    # s_x = v_x*conj(i_xy-i_zx)
    # p_x = vm_x*cos(va_x)*(ire_xy-ire_zx) + vm_x*sin(va_x)*(iim_xy-iim_zx)
    # q_x = vm_x*sin(va_x)*(ire_xy-ire_zx) - vm_x*cos(va_x)*(iim_xy-iim_zx)
    p_x(vm_x, va_x, ire_xy, iim_xy, ire_zx, iim_zx) = JuMP.@NLexpression(pm.model, vm_x*cos(va_x)*(ire_xy-ire_zx) + vm_x*sin(va_x)*(iim_xy-iim_zx))
    q_x(vm_x, va_x, ire_xy, iim_xy, ire_zx, iim_zx) = JuMP.@NLexpression(pm.model, vm_x*sin(va_x)*(ire_xy-ire_zx) - vm_x*cos(va_x)*(iim_xy-iim_zx))
    # s_x = s_x,ref
    _PMs.var(pm, nw, 1, :pd)[load_id] = p_x(vm_a, va_a, ire_ab, iim_ab, ire_ca, iim_ca)
    _PMs.var(pm, nw, 2, :pd)[load_id] = p_x(vm_b, va_b, ire_bc, iim_bc, ire_ab, iim_ab)
    _PMs.var(pm, nw, 3, :pd)[load_id] = p_x(vm_c, va_c, ire_ca, iim_ca, ire_bc, iim_bc)
    _PMs.var(pm, nw, 1, :qd)[load_id] = q_x(vm_a, va_a, ire_ab, iim_ab, ire_ca, iim_ca)
    _PMs.var(pm, nw, 2, :qd)[load_id] = q_x(vm_b, va_b, ire_bc, iim_bc, ire_ab, iim_ab)
    _PMs.var(pm, nw, 3, :qd)[load_id] = q_x(vm_c, va_c, ire_ca, iim_ca, ire_bc, iim_bc)
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_vm_zero_seq(pm::_PMs.AbstractACPModel, nw::Int, bus_id::Int, vmzeromax::Float64)
    if !haskey(_PMs.var(pm, pm.cnw), :vmpossqr)
        _PMs.var(pm, pm.cnw)[:vmpossqr] = Dict{Int, Any}()
        _PMs.var(pm, pm.cnw)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [_PMs.var(pm, nw, i, :vm, bus_id) for i in 1:3]
    (va_a, va_b, va_c) = [_PMs.var(pm, nw, i, :va, bus_id) for i in 1:3]
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
function constraint_mc_load_current_delta(pm::_PMs.AbstractACPModel, nw::Int, load_id::Int, load_bus_id::Int, cp::_PMs.MultiConductorVector, cq::_PMs.MultiConductorVector)
    cp_ab, cp_bc, cp_ca = cp
    cq_ab, cq_bc, cq_ca = cq
    vm_a, vm_b, vm_c = [_PMs.var(pm, nw, c, :vm, load_bus_id) for c in 1:3]
    va_a, va_b, va_c = [_PMs.var(pm, nw, c, :va, load_bus_id) for c in 1:3]
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
    _PMs.var(pm, nw, 1, :pd)[load_id] = p_x(vm_a, va_a, ire_ab, iim_ab, ire_ca, iim_ca)
    _PMs.var(pm, nw, 2, :pd)[load_id] = p_x(vm_b, va_b, ire_bc, iim_bc, ire_ab, iim_ab)
    _PMs.var(pm, nw, 3, :pd)[load_id] = p_x(vm_c, va_c, ire_ca, iim_ca, ire_bc, iim_bc)
    _PMs.var(pm, nw, 1, :qd)[load_id] = q_x(vm_a, va_a, ire_ab, iim_ab, ire_ca, iim_ca)
    _PMs.var(pm, nw, 2, :qd)[load_id] = q_x(vm_b, va_b, ire_bc, iim_bc, ire_ab, iim_ab)
    _PMs.var(pm, nw, 3, :qd)[load_id] = q_x(vm_c, va_c, ire_ca, iim_ca, ire_bc, iim_bc)
end


""
function constraint_mc_vm_ll(pm::_PMs.AbstractACPModel, nw::Int, bus_id::Int, vm_ll_min::_PMs.MultiConductorVector, vm_ll_max::_PMs.MultiConductorVector)
    # 3 conductors asserted in template already
    vm_ln = [_PMs.var(pm, nw, i, :vm, bus_id) for i in 1:3]
    va_ln = [_PMs.var(pm, nw, i, :va, bus_id) for i in 1:3]
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


"""
We want to express
s_ab = cp.|v_ab|^2+im.cq.|v_ab|^2
i_ab = conj(s_ab/v_ab) = |v_ab|^2.(cq-im.cq)/conj(v_ab) = (cp-im.cq)*v_ab
idem for i_bc and i_ca
And then
s_a = v_a.conj(i_a) = v_a.conj(i_ab-i_ca)
idem for s_b and s_c
"""
function constraint_mc_load_impedance_delta(pm::_PMs.AbstractACPModel, nw::Int, load_id::Int, load_bus_id::Int, cp::_PMs.MultiConductorVector, cq::_PMs.MultiConductorVector)
    cp_ab, cp_bc, cp_ca = cp
    cq_ab, cq_bc, cq_ca = cq
    vm_a, vm_b, vm_c = [_PMs.var(pm, nw, c, :vm, load_bus_id) for c in 1:3]
    va_a, va_b, va_c = [_PMs.var(pm, nw, c, :va, load_bus_id) for c in 1:3]
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
    ire_xy(cp_xy, cq_xy, vre_xy, vim_xy) = JuMP.@NLexpression(pm.model, cp_xy*vre_xy+cq_xy*vim_xy)
    iim_xy(cp_xy, cq_xy, vre_xy, vim_xy) = JuMP.@NLexpression(pm.model, cp_xy*vim_xy-cq_xy*vre_xy)
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
    _PMs.var(pm, nw, 1, :pd)[load_id] = p_x(vm_a, va_a, ire_ab, iim_ab, ire_ca, iim_ca)
    _PMs.var(pm, nw, 2, :pd)[load_id] = p_x(vm_b, va_b, ire_bc, iim_bc, ire_ab, iim_ab)
    _PMs.var(pm, nw, 3, :pd)[load_id] = p_x(vm_c, va_c, ire_ca, iim_ca, ire_bc, iim_bc)
    _PMs.var(pm, nw, 1, :qd)[load_id] = q_x(vm_a, va_a, ire_ab, iim_ab, ire_ca, iim_ca)
    _PMs.var(pm, nw, 2, :qd)[load_id] = q_x(vm_b, va_b, ire_bc, iim_bc, ire_ab, iim_ab)
    _PMs.var(pm, nw, 3, :qd)[load_id] = q_x(vm_c, va_c, ire_ca, iim_ca, ire_bc, iim_bc)
end


"bus voltage on/off constraint for load shed problem"
function constraint_mc_bus_voltage_on_off(pm::_PMs.AbstractACPModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, kwargs...)
    for (i,bus) in _PMs.ref(pm, nw, :bus)
        constraint_mc_voltage_magnitude_on_off(pm, i; nw=nw, cnd=cnd)
    end
end
