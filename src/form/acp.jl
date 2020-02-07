""
function variable_mc_voltage(pm::_PMs.AbstractACPModel; nw=pm.cnw, kwargs...)
    variable_mc_voltage_angle(pm; nw=nw, kwargs...)
    variable_mc_voltage_magnitude(pm; nw=nw, kwargs...)

    # This is needed for delta loads, where division occurs by the difference
    # of voltage phasors. If the voltage phasors at one bus are initialized
    # in the same point, this would lead to division by zero.

    ncnds = length(_PMs.conductor_ids(pm, nw))

    bus_t1 = [bus for (_, bus) in _PMs.ref(pm, nw, :bus) if bus["bus_type"]==3]
    if length(bus_t1)>0
        theta = bus_t1[1]["va"]
    else
        theta = [_wrap_to_pi(2 * pi / ncnds * (1-c)) for c in 1:ncnds]
    end
    vm = 1
    for id in _PMs.ids(pm, nw, :bus)
        busref = _PMs.ref(pm, nw, :bus, id)
        if !haskey(busref, "va_start")
        # if it has this key, it was set at PM level
            for c in 1:ncnds
                JuMP.set_start_value(_PMs.var(pm, nw, :va, id)[c], theta[c])
            end
        end
    end
end


""
function variable_mc_bus_voltage_on_off(pm::_PMs.AbstractACPModel; kwargs...)
    variable_mc_voltage_angle(pm; kwargs...)
    variable_mc_voltage_magnitude_on_off(pm; kwargs...)

    nw = get(kwargs, :nw, pm.cnw)


    vm = 1
    for id in _PMs.ids(pm, nw, :bus)
        busref = _PMs.ref(pm, nw, :bus, id)
        ncnd = length(_PMs.conductor_ids(pm, nw))
        theta = [_wrap_to_pi(2 * pi / ncnd * (1-c)) for c in 1:ncnd]
        if !haskey(busref, "va_start")
        # if it has this key, it was set at PM level
            for c in 1:ncnd
                JuMP.set_start_value(_PMs.var(pm, nw, :va, id)[c], theta[c])
            end
        end
    end
end


""
function constraint_mc_power_balance_slack(pm::_PMs.AbstractACPModel, nw::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm   = _PMs.var(pm, nw, :vm, i)
    p    = get(_PMs.var(pm, nw),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(_PMs.var(pm, nw),    :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(_PMs.var(pm, nw),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(_PMs.var(pm, nw),   :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PMs.var(pm, nw),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PMs.var(pm, nw),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(_PMs.var(pm, nw),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(_PMs.var(pm, nw),  :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(_PMs.var(pm, nw),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(_PMs.var(pm, nw),   :qt, Dict()); _PMs._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    p_slack = _PMs.var(pm, nw, :p_slack, i)
    q_slack = _PMs.var(pm, nw, :q_slack, i)

    cstr_p = []
    cstr_q = []

    for c in _PMs.conductor_ids(pm; nw=nw)
        cp = JuMP.@constraint(pm.model,
            sum(p[a][c] for a in bus_arcs)
            + sum(psw[a_sw][c] for a_sw in bus_arcs_sw)
            + sum(pt[a_trans][c] for a_trans in bus_arcs_trans)
            ==
            sum(pg[g][c] for g in bus_gens)
            - sum(ps[s][c] for s in bus_storage)
            - sum(pd[c] for pd in values(bus_pd))
            - sum(gs[c] for gs in values(bus_gs))*vm[c]^2
            + p_slack[c]
        )
        push!(cstr_p, cp)

        cq = JuMP.@constraint(pm.model,
            sum(q[a][c] for a in bus_arcs)
            + sum(qsw[a_sw][c] for a_sw in bus_arcs_sw)
            + sum(qt[a_trans][c] for a_trans in bus_arcs_trans)
            ==
            sum(qg[g][c] for g in bus_gens)
            - sum(qs[s][c] for s in bus_storage)
            - sum(qd[c] for qd in values(bus_qd))
            + sum(bs[c] for bs in values(bus_bs))*vm[c]^2
            + q_slack[c]
        )
        push!(cstr_q, cq)

    end

    if _PMs.report_duals(pm)
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_mc_power_balance_shed(pm::_PMs.AbstractACPModel, nw::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm       = _PMs.var(pm, nw, :vm, i)
    p        = get(_PMs.var(pm, nw),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q        = get(_PMs.var(pm, nw),    :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg       = get(_PMs.var(pm, nw),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg       = get(_PMs.var(pm, nw),   :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps       = get(_PMs.var(pm, nw),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs       = get(_PMs.var(pm, nw),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw      = get(_PMs.var(pm, nw),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw      = get(_PMs.var(pm, nw),  :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt       = get(_PMs.var(pm, nw),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt       = get(_PMs.var(pm, nw),   :qt, Dict()); _PMs._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    z_demand = _PMs.var(pm, nw, :z_demand)
    z_shunt  = _PMs.var(pm, nw, :z_shunt)

    cstr_p = []
    cstr_q = []

    for c in _PMs.conductor_ids(pm; nw=nw)
        cp = JuMP.@constraint(pm.model,
            sum(p[a][c] for a in bus_arcs)
            + sum(psw[a_sw][c] for a_sw in bus_arcs_sw)
            + sum(pt[a_trans][c] for a_trans in bus_arcs_trans)
            ==
            sum(pg[g][c] for g in bus_gens)
            - sum(ps[s][c] for s in bus_storage)
            - sum(pd[c]*z_demand[n] for (n,pd) in bus_pd)
            - sum(gs[c]*z_shunt[n] for (n,gs) in bus_gs)*vm[c]^2
        )
        push!(cstr_p, cp)

        cq = JuMP.@constraint(pm.model,
            sum(q[a][c] for a in bus_arcs)
            + sum(qsw[a_sw][c] for a_sw in bus_arcs_sw)
            + sum(qt[a_trans][c] for a_trans in bus_arcs_trans)
            ==
            sum(qg[g][c] for g in bus_gens)
            - sum(qs[s][c] for s in bus_storage)
            - sum(qd[c]*z_demand[n] for (n,qd) in bus_qd)
            + sum(bs[c]*z_shunt[n] for (n,bs) in bus_bs)*vm[c]^2
        )
        push!(cstr_q, cq)

    end

    if _PMs.report_duals(pm)
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_mc_power_balance(pm::_PMs.AbstractACPModel, nw::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm   = _PMs.var(pm, nw, :vm, i)
    p    = get(_PMs.var(pm, nw),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(_PMs.var(pm, nw),    :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(_PMs.var(pm, nw),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(_PMs.var(pm, nw),   :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PMs.var(pm, nw),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PMs.var(pm, nw),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(_PMs.var(pm, nw),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(_PMs.var(pm, nw),  :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(_PMs.var(pm, nw),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(_PMs.var(pm, nw),   :qt, Dict()); _PMs._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")

    cstr_p = []
    cstr_q = []

    for c in _PMs.conductor_ids(pm; nw=nw)
        cp = JuMP.@constraint(pm.model,
            sum(p[a][c] for a in bus_arcs)
            + sum(psw[a_sw][c] for a_sw in bus_arcs_sw)
            + sum(pt[a_trans][c] for a_trans in bus_arcs_trans)
            ==
            sum(pg[g][c] for g in bus_gens)
            - sum(ps[s][c] for s in bus_storage)
            - sum(pd[c] for pd in values(bus_pd))
            - sum(gs[c] for gs in values(bus_gs))*vm[c]^2
        )
        push!(cstr_p, cp)

        cq = JuMP.@constraint(pm.model,
            sum(q[a][c] for a in bus_arcs)
            + sum(qsw[a_sw][c] for a_sw in bus_arcs_sw)
            + sum(qt[a_trans][c] for a_trans in bus_arcs_trans)
            ==
            sum(qg[g][c] for g in bus_gens)
            - sum(qs[s][c] for s in bus_storage)
            - sum(qd[c] for qd in values(bus_qd))
            + sum(bs[c] for bs in values(bus_bs))*vm[c]^2
        )
        push!(cstr_q, cq)
    end

    if _PMs.report_duals(pm)
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_mc_power_balance_load(pm::_PMs.AbstractACPModel, nw::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_gs, bus_bs)
    vm   = _PMs.var(pm, nw, :vm, i)
    p    = get(_PMs.var(pm, nw),   :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(_PMs.var(pm, nw),   :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(_PMs.var(pm, nw),  :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(_PMs.var(pm, nw),  :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PMs.var(pm, nw),  :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PMs.var(pm, nw),  :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(_PMs.var(pm, nw), :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(_PMs.var(pm, nw), :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(_PMs.var(pm, nw),  :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(_PMs.var(pm, nw),  :qt, Dict()); _PMs._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    pd   = get(_PMs.var(pm, nw),  :pd, Dict()); _PMs._check_var_keys(pd, bus_loads, "active power", "load")
    qd   = get(_PMs.var(pm, nw),  :qd, Dict()); _PMs._check_var_keys(pd, bus_loads, "reactive power", "load")

    cstr_p = []
    cstr_q = []

    for c in _PMs.conductor_ids(pm; nw=nw)
        cp = JuMP.@NLconstraint(pm.model,
            sum(p[a][c] for a in bus_arcs)
            + sum(psw[a_sw][c] for a_sw in bus_arcs_sw)
            + sum(pt[a_trans][c] for a_trans in bus_arcs_trans)
            ==
            sum(pg[g][c] for g in bus_gens)
            - sum(ps[s][c] for s in bus_storage)
            - sum(pd[l][c] for l in bus_loads)
            - sum(gs[c] for gs in values(bus_gs))*vm[c]^2
        )
        push!(cstr_p, cp)

        cq = JuMP.@NLconstraint(pm.model,
            sum(q[a][c] for a in bus_arcs)
            + sum(qsw[a_sw][c] for a_sw in bus_arcs_sw)
            + sum(qt[a_trans][c] for a_trans in bus_arcs_trans)
            ==
            sum(qg[g][c] for g in bus_gens)
            - sum(qs[s][c] for s in bus_storage)
            - sum(qd[l][c] for l in bus_loads)
            + sum(bs[c] for bs in values(bus_bs))*vm[c]^2
        )
        push!(cstr_q, cq)
    end

    if _PMs.report_duals(pm)
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
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
function constraint_mc_ohms_yt_from(pm::_PMs.AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr  = _PMs.var(pm, n,  :p, f_idx)
    q_fr  = _PMs.var(pm, n,  :q, f_idx)
    vm_fr = _PMs.var(pm, n, :vm, f_bus)
    vm_to = _PMs.var(pm, n, :vm, t_bus)
    va_fr = _PMs.var(pm, n, :va, f_bus)
    va_to = _PMs.var(pm, n, :va, t_bus)

    cnds = _PMs.conductor_ids(pm; nw=n)
    for c in cnds
        JuMP.@NLconstraint(pm.model, p_fr[c] ==
                                        (g[c,c]+g_fr[c,c])*vm_fr[c]^2
                                        +sum( (g[c,d]+g_fr[c,d]) * vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d])
                                             +(b[c,d]+b_fr[c,d]) * vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d])
                                             for d in cnds if d != c)
                                        +sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d])
                                             -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d])
                                             for d in cnds)
                                    )
        JuMP.@NLconstraint(pm.model, q_fr[c] == -(b[c,c]+b_fr[c,c])*vm_fr[c]^2
                                        -sum( (b[c,d]+b_fr[c,d])*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d])
                                             -(g[c,d]+g_fr[c,d])*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d])
                                             for d in cnds if d != c)
                                        -sum(-b[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d])
                                             +g[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d])
                                             for d in cnds)
                                    )
    end
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_mc_ohms_yt_to(pm::_PMs.AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    constraint_mc_ohms_yt_from(pm, n, t_bus, f_bus, t_idx, f_idx, g, b, g_to, b_to, tr, ti, tm)
end


function constraint_mc_trans_yy(pm::_PMs.AbstractACPModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
    vm_fr = _PMs.var(pm, nw, :vm, f_bus)[f_cnd]
    vm_to = _PMs.var(pm, nw, :vm, t_bus)[t_cnd]
    va_fr = _PMs.var(pm, nw, :va, f_bus)[f_cnd]
    va_to = _PMs.var(pm, nw, :va, t_bus)[t_cnd]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[p] ? tm_set[p] : _PMs.var(pm, nw, :tap, trans_id)[p] for p in _PMs.conductor_ids(pm)]


    for p in _PMs.conductor_ids(pm)
        if tm_fixed[p]
            JuMP.@constraint(pm.model, vm_fr[p] == tm_scale*tm[p]*vm_to[p])
        else
            JuMP.@NLconstraint(pm.model, vm_fr[p] == tm_scale*tm[p]*vm_to[p])
        end
        pol_angle = pol==1 ? 0 : pi
        JuMP.@constraint(pm.model, va_fr[p] == va_to[p] + pol_angle)
    end

    p_fr = _PMs.var(pm, nw, :pt, f_idx)[f_cnd]
    p_to = _PMs.var(pm, nw, :pt, t_idx)[t_cnd]
    q_fr = _PMs.var(pm, nw, :qt, f_idx)[f_cnd]
    q_to = _PMs.var(pm, nw, :qt, t_idx)[t_cnd]

    JuMP.@constraint(pm.model, p_fr + p_to .== 0)
    JuMP.@constraint(pm.model, q_fr + q_to .== 0)
end


function constraint_mc_trans_dy(pm::_PMs.AbstractACPModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
    vm_fr = [_PMs.var(pm, nw, :vm, f_bus)[p] for p in f_cnd]
    vm_to = [_PMs.var(pm, nw, :vm, t_bus)[p] for p in t_cnd]
    va_fr = [_PMs.var(pm, nw, :va, f_bus)[p] for p in f_cnd]
    va_to = [_PMs.var(pm, nw, :va, t_bus)[p] for p in t_cnd]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[p] ? tm_set[p] : _PMs.var(pm, nw, p, :tap, trans_id) for p in _PMs.conductor_ids(pm)]

    # introduce auxialiary variable vd = Md*v_fr
    nph = length(_PMs.conductor_ids(pm))
    vd_re = Array{Any,1}(undef, nph)
    vd_im = Array{Any,1}(undef, nph)
    for p in 1:nph
        # rotate by 1 to get 'previous' phase
        # e.g., for nph=3: 1->3, 2->1, 3->2
        q = (p-1+1)%nph+1
        vd_re[p] = JuMP.@NLexpression(pm.model, vm_fr[p]*cos(va_fr[p])-vm_fr[q]*cos(va_fr[q]))
        vd_im[p] = JuMP.@NLexpression(pm.model, vm_fr[p]*sin(va_fr[p])-vm_fr[q]*sin(va_fr[q]))
        JuMP.@NLconstraint(pm.model, vd_re[p] == pol*tm_scale*tm[p]*vm_to[p]*cos(va_to[p]))
        JuMP.@NLconstraint(pm.model, vd_im[p] == pol*tm_scale*tm[p]*vm_to[p]*sin(va_to[p]))
    end

    p_fr = [_PMs.var(pm, nw, :pt, f_idx)[p] for p in f_cnd]
    p_to = [_PMs.var(pm, nw, :pt, t_idx)[p] for p in t_cnd]
    q_fr = [_PMs.var(pm, nw, :qt, f_idx)[p] for p in f_cnd]
    q_to = [_PMs.var(pm, nw, :qt, t_idx)[p] for p in t_cnd]

    id_re = Array{Any,1}(undef, nph)
    id_im = Array{Any,1}(undef, nph)
    # s/v      = (p+jq)/|v|^2*conj(v)
    #          = (p+jq)/|v|*(cos(va)-j*sin(va))
    # Re(s/v)  = (p*cos(va)+q*sin(va))/|v|
    # -Im(s/v) = -(q*cos(va)-p*sin(va))/|v|
    for p in _PMs.conductor_ids(pm)
        # id = conj(s_to/v_to)./tm
        id_re[p] = JuMP.@NLexpression(pm.model,  (p_to[p]*cos(va_to[p])+q_to[p]*sin(va_to[p]))/vm_to[p]/(tm_scale*tm[p])/pol)
        id_im[p] = JuMP.@NLexpression(pm.model, -(q_to[p]*cos(va_to[p])-p_to[p]*sin(va_to[p]))/vm_to[p]/(tm_scale*tm[p])/pol)
    end
    for p in _PMs.conductor_ids(pm)
        # rotate by nph-1 to get 'previous' phase
        # e.g., for nph=3: 1->3, 2->1, 3->2
        q = (p-1+nph-1)%nph+1
        # s_fr  = v_fr*conj(i_fr)
        #       = v_fr*conj(id[q]-id[p])
        #       = v_fr*(id_re[q]-j*id_im[q]-id_re[p]+j*id_im[p])
        JuMP.@NLconstraint(pm.model, p_fr[p] ==
             vm_fr[p]*cos(va_fr[p])*(id_re[q]-id_re[p])
            -vm_fr[p]*sin(va_fr[p])*(-id_im[q]+id_im[p])
        )
        JuMP.@NLconstraint(pm.model, q_fr[p] ==
             vm_fr[p]*cos(va_fr[p])*(-id_im[q]+id_im[p])
            +vm_fr[p]*sin(va_fr[p])*(id_re[q]-id_re[p])
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


#= TODO unused function, remove?
"Links the power flowing into both windings of a variable tap transformer."
function constraint_mc_transformer_flow_var(pm::_PMs.AbstractPowerModel, i::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, Ti_fr, Ti_im; nw::Int=pm.cnw)
    # for ac formulation, indentical to fixed tap
    constraint_mc_transformer_flow(pm, i, f_bus, t_bus, f_idx, t_idx, Ti_fr, Ti_im)
end
=#


""
function constraint_mc_load_power_wye(pm::_PMs.AbstractACPModel, nw::Int, load_id::Int, pd, qd)
    _PMs.var(pm, nw, :pd)[load_id] = pd
    _PMs.var(pm, nw, :qd)[load_id] = qd
end


""
function constraint_mc_load_current_wye(pm::_PMs.AbstractACPModel, nw::Int, load_id::Int, load_bus_id::Int, scale_p, scale_q)
    vm = _PMs.var(pm, nw, :vm, load_bus_id)
    # this has to be a NLexpression and not an expression;
    # might be mixed in KCL with NLexpression, so has to also be NLexpression

    ncnds = length(_PMs.conductor_ids(pm, nw))

    _PMs.var(pm, nw, :pd)[load_id] = JuMP.@NLexpression(pm.model, [c in 1:ncnds], scale_p[c]*vm[c])
    _PMs.var(pm, nw, :qd)[load_id] = JuMP.@NLexpression(pm.model, [c in 1:ncnds], scale_q[c]*vm[c])
end


""
function constraint_mc_load_impedance_wye(pm::_PMs.AbstractACPModel, nw::Int, load_id::Int, load_bus_id::Int, cp, cq)
    vm = _PMs.var(pm, nw, :vm, load_bus_id)

    ncnds = length(_PMs.conductor_ids(pm, nw))

    _PMs.var(pm, nw, :pd)[load_id] = JuMP.@NLexpression(pm.model, [c in 1:ncnds], cp[c]*vm[c]^2)
    _PMs.var(pm, nw, :qd)[load_id] = JuMP.@NLexpression(pm.model, [c in 1:ncnds], cq[c]*vm[c]^2)
end


""
function constraint_mc_load_exponential_wye(pm::_PMs.AbstractACPModel, nw::Int, load_id::Int, load_bus_id::Int, a, alpha, b, beta)
    vm = _PMs.var(pm, nw, :vm, load_bus_id)

    cnds = length(_PMs.conductor_ids(pm, nw))

    _PMs.var(pm, nw, cnd, :pd)[load_id] = JuMP.@NLexpression(pm.model, [c in 1:ncnds], a[c]*vm[c]^alpha[c])
    _PMs.var(pm, nw, cnd, :qd)[load_id] = JuMP.@NLexpression(pm.model, [c in 1:ncnds],b[c]*vm[c]^beta[c])
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
function constraint_mc_load_power_delta(pm::_PMs.AbstractACPModel, nw::Int, load_id::Int, load_bus_id::Int, pd::Vector, qd::Vector)
    p_ab, p_bc, p_ca = pd
    q_ab, q_bc, q_ca = qd
    vm_a, vm_b, vm_c = _PMs.var(pm, nw, :vm, load_bus_id)
    va_a, va_b, va_c = _PMs.var(pm, nw, :va, load_bus_id)
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
    _PMs.var(pm, nw, :pd)[load_id] = [
        p_x(vm_a, va_a, ire_ab, iim_ab, ire_ca, iim_ca),
        p_x(vm_b, va_b, ire_bc, iim_bc, ire_ab, iim_ab),
        p_x(vm_c, va_c, ire_ca, iim_ca, ire_bc, iim_bc)]
    _PMs.var(pm, nw, :qd)[load_id] = [
        q_x(vm_a, va_a, ire_ab, iim_ab, ire_ca, iim_ca),
        q_x(vm_b, va_b, ire_bc, iim_bc, ire_ab, iim_ab),
        q_x(vm_c, va_c, ire_ca, iim_ca, ire_bc, iim_bc)]
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
function constraint_mc_load_current_delta(pm::_PMs.AbstractACPModel, nw::Int, load_id::Int, load_bus_id::Int, cp::Vector, cq::Vector)
    cp_ab, cp_bc, cp_ca = cp
    cq_ab, cq_bc, cq_ca = cq
    vm_a, vm_b, vm_c = _PMs.var(pm, nw, :vm, load_bus_id)
    va_a, va_b, va_c = _PMs.var(pm, nw, :va, load_bus_id)
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
    _PMs.var(pm, nw, :pd)[load_id] = [
        p_x(vm_a, va_a, ire_ab, iim_ab, ire_ca, iim_ca),
        p_x(vm_b, va_b, ire_bc, iim_bc, ire_ab, iim_ab),
        p_x(vm_c, va_c, ire_ca, iim_ca, ire_bc, iim_bc)]
    _PMs.var(pm, nw, :qd)[load_id] = [
        q_x(vm_a, va_a, ire_ab, iim_ab, ire_ca, iim_ca),
        q_x(vm_b, va_b, ire_bc, iim_bc, ire_ab, iim_ab),
        q_x(vm_c, va_c, ire_ca, iim_ca, ire_bc, iim_bc)]
end


""
function constraint_mc_vm_ll(pm::_PMs.AbstractACPModel, nw::Int, bus_id::Int, vm_ll_min::Vector, vm_ll_max::Vector)
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
function constraint_mc_load_impedance_delta(pm::_PMs.AbstractACPModel, nw::Int, load_id::Int, load_bus_id::Int, cp::Vector, cq::Vector)
    cp_ab, cp_bc, cp_ca = cp
    cq_ab, cq_bc, cq_ca = cq
    vm_a, vm_b, vm_c = _PMs.var(pm, nw, :vm, load_bus_id)
    va_a, va_b, va_c = _PMs.var(pm, nw, :va, load_bus_id)
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
    _PMs.var(pm, nw, :pd)[load_id] = [
        p_x(vm_a, va_a, ire_ab, iim_ab, ire_ca, iim_ca),
        p_x(vm_b, va_b, ire_bc, iim_bc, ire_ab, iim_ab),
        p_x(vm_c, va_c, ire_ca, iim_ca, ire_bc, iim_bc)]
    _PMs.var(pm, nw, :qd)[load_id] = [
        q_x(vm_a, va_a, ire_ab, iim_ab, ire_ca, iim_ca),
        q_x(vm_b, va_b, ire_bc, iim_bc, ire_ab, iim_ab),
        q_x(vm_c, va_c, ire_ca, iim_ca, ire_bc, iim_bc)]
end


"""
exponential model in acp voltage coordinates
"""
function constraint_mc_load_exponential_delta(pm::_PMs.AbstractACPModel, nw::Int, load_id::Int, load_bus_id::Int, a, alpha, b, beta)
    nph = length(_PMs.conductor_ids(pm))

    vm = _PMs.var(pm, nw, c, :vm, load_bus_id)
    va = _PMs.var(pm, nw, c, :va, load_bus_id)

    # vd = Td*v
    vdr = JuMP.@NLexpression(pm.model, [i in 1:nph], vm[i]*cos(va[i])-vm[mod(i,nph)+1]*cos(va[mod(i,nph)+1]))
    vdi = JuMP.@NLexpression(pm.model, [i in 1:nph], vm[i]*sin(va[i])-vm[mod(i,nph)+1]*sin(va[mod(i,nph)+1]))
    pl = JuMP.@NLexpression(pm.model, [i in 1:nph], a[i]*(vdr[i]^2+vdi[i]^2)^(alpha[i]/2))
    ql = JuMP.@NLexpression(pm.model, [i in 1:nph], b[i]*(vdr[i]^2+vdi[i]^2)^(beta[i]/2))

    # idc: conjugate of delta current
    # idc = sl./vd
    #     = (pl + j.ql)/vd = (a.|vd|^alpha + j.b.|vd|^beta)/vd, now use 1/c = conj(c)/|c|^2
    #     = (vdr-j.vdi).(a.|vd|^alpha + j.b.|vd|^beta)/|vd|^2
    #     = (vdr-j.vdi).(a.|vd|^(alpha-2) + j.b.|vd|^(beta-2))
    #  now also use that |vd| = (vdr^2 + wdi^2)^(1/2)
    idcr = JuMP.@NLexpression(pm.model, [i in 1:nph], vdr[i]*a[i]*(vdr[i]^2+vdi[i]^2)^(alpha[i]/2-1) + vdi[i]*b[i]*(vdr[i]^2+vdi[i]^2)^(beta[i]/2-1))
    idci = JuMP.@NLexpression(pm.model, [i in 1:nph], vdr[i]*b[i]*(vdr[i]^2+vdi[i]^2)^(beta[i]/2-1) - vdi[i]*a[i]*(vdr[i]^2+vdi[i]^2)^(alpha[i]/2-1))

    # icr = Td'*idcr
    icr = JuMP.@NLexpression(pm.model, [i in 1:nph], idcr[i]-idcr[mod(i+nph-2,nph)+1])
    ici = JuMP.@NLexpression(pm.model, [i in 1:nph], idci[i]-idci[mod(i+nph-2,nph)+1])

    # s = v.*conj(ic)
    pd = JuMP.@NLexpression(pm.model, [i in 1:nph], vm[i]*cos(va[i])*icr[i] - vm[i]*sin(va[i])*ici[i])
    qd = JuMP.@NLexpression(pm.model, [i in 1:nph], vm[i]*cos(va[i])*ici[i] + vm[i]*sin(va[i])*icr[i])

    _PMs.var(pm, nw, :pd)[load_id] = pd
    _PMs.var(pm, nw, :qd)[load_id] = qd

    # also save pl, ql for inspection
    for c in _PMs.conductor_ids(pm)
        if !haskey(_PMs.var(pm, nw, c), :pl)
            _PMs.var(pm, nw, c)[:pl] = Dict()
            _PMs.var(pm, nw, c)[:ql] = Dict()
        end
        _PMs.var(pm, nw, c, :pl)[load_id] = pl[c]
        _PMs.var(pm, nw, c, :ql)[load_id] = ql[c]
    end
end


"bus voltage on/off constraint for load shed problem"
function constraint_mc_bus_voltage_on_off(pm::_PMs.AbstractACPModel; nw::Int=pm.cnw, kwargs...)
    for (i,bus) in _PMs.ref(pm, nw, :bus)
        constraint_mc_voltage_magnitude_on_off(pm, i; nw=nw)
    end
end

"`vm[i] == vmref`"
function constraint_mc_voltage_magnitude_setpoint(pm::_PMs.AbstractACPModel, n::Int, i::Int, vmref)
    vm = _PMs.var(pm, n, :vm, i)

    JuMP.@constraint(pm.model, vm .== vmref)
end

""
function constraint_mc_storage_current_limit(pm::_PMs.AbstractACPModel, n::Int, i, bus, rating)
    vm = var(pm, n, :vm, bus)
    ps = var(pm, n, :ps, i)
    qs = var(pm, n, :qs, i)

    JuMP.@constraint(pm.model, ps.^2 + qs.^2 <= rating.^2 .* vm.^2)
end
