""
function constraint_mc_power_balance_slack(pm::_PMs.AbstractWModels, nw::Int, c::Int, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    w    = _PMs.var(pm, nw, c, :w, i)
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
        - sum(gs for gs in values(bus_gs))*w
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
        + sum(bs for bs in values(bus_bs))*w
        + q_slack
    )
end


"Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)"
function constraint_mc_ohms_yt_from(pm::_PMs.AbstractWRModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr = _PMs.var(pm, n, c, :p, f_idx)
    q_fr = _PMs.var(pm, n, c, :q, f_idx)
    w    = _PMs.var(pm, n, :w)
    wr   = _PMs.var(pm, n, :wr)
    wi   = _PMs.var(pm, n, :wi)

    #TODO extend to shunt matrices; this ignores the off-diagonals
    JuMP.@constraint(pm.model, p_fr ==  ( g_fr[c,c]+g[c,c]) * w[(f_bus, c)] +
                                sum( g[c,d] * wr[(f_bus, f_bus, c, d)] +
                                     b[c,d] * wi[(f_bus, f_bus, c, d)] for d in _PMs.conductor_ids(pm) if d != c) +
                                sum(-g[c,d] * wr[(f_bus, t_bus, c, d)] +
                                    -b[c,d] * wi[(f_bus, t_bus, c, d)] for d in _PMs.conductor_ids(pm)) )
    JuMP.@constraint(pm.model, q_fr == -( b_fr[c,c]+b[c,c]) * w[(f_bus, c)] -
                                sum( b[c,d] * wr[(f_bus, f_bus, c, d)] -
                                     g[c,d] * wi[(f_bus, f_bus, c, d)] for d in _PMs.conductor_ids(pm) if d != c) -
                                sum(-b[c,d] * wr[(f_bus, t_bus, c, d)] +
                                     g[c,d] * wi[(f_bus, t_bus, c, d)] for d in _PMs.conductor_ids(pm)) )
end


"Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)"
function constraint_mc_ohms_yt_to(pm::_PMs.AbstractWRModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    q_to = _PMs.var(pm, n, c, :q, t_idx)
    p_to = _PMs.var(pm, n, c, :p, t_idx)
    w    = _PMs.var(pm, n, :w)
    wr   = _PMs.var(pm, n, :wr)
    wi   = _PMs.var(pm, n, :wi)

    #TODO extend to shunt matrices; this ignores the off-diagonals
    JuMP.@constraint(pm.model, p_to ==  ( g_to[c,c]+g[c,c]) * w[(t_bus, c)] +
                                sum( g[c,d] * wr[(t_bus, t_bus, c, d)] +
                                     b[c,d] *-wi[(t_bus, t_bus, c, d)] for d in _PMs.conductor_ids(pm) if d != c) +
                                sum(-g[c,d] * wr[(f_bus, t_bus, c, d)] +
                                    -b[c,d] *-wi[(f_bus, t_bus, c, d)] for d in _PMs.conductor_ids(pm)) )
    JuMP.@constraint(pm.model, q_to == -( b_to[c,c]+b[c,c]) * w[(t_bus, c)] -
                                sum( b[c,d] * wr[(t_bus, t_bus, c, d)] -
                                     g[c,d] *-wi[(t_bus, t_bus, c, d)] for d in _PMs.conductor_ids(pm) if d != c) -
                                sum(-b[c,d] * wr[(f_bus, t_bus, c, d)] +
                                     g[c,d] *-wi[(f_bus, t_bus, c, d)] for d in _PMs.conductor_ids(pm)) )
end


"do nothing, no way to represent this in these variables"
function constraint_mc_theta_ref(pm::_PMs.AbstractWModels, n::Int, c::Int, d)
end


"Creates phase angle constraints at reference buses"
function constraint_mc_theta_ref(pm::_PMs.AbstractPolarModels, n::Int, c::Int, d)
    va = _PMs.var(pm, n, c, :va, d)
    nconductors = length(_PMs.conductor_ids(pm))

    JuMP.@constraint(pm.model, va == _wrap_to_pi(2 * pi / nconductors * (1-c)))
end


"""
For a variable tap transformer, fix the tap variables which are fixed. For
example, an OLTC where the third phase is fixed, will have tap variables for
all phases, but the third tap variable should be fixed.
"""
function constraint_mc_oltc_tap_fix(pm::_PMs.AbstractPowerModel, i::Int, fixed::MultiConductorVector, tm::MultiConductorVector; nw=pm.cnw)
    for (c,fixed) in enumerate(fixed)
        if fixed
            JuMP.@constraint(pm.model, _PMs.var(pm, nw, c, :tap)[i]==tm[c])
        end
    end
end


"KCL for load shed problem with transformers (AbstractWForms)"
function constraint_mc_power_balance_shed(pm::_PMs.AbstractWModels, nw::Int, c::Int, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    w    = _PMs.var(pm, nw, c, :w, i)
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
        - sum(gs*1.0^2*z_shunt[n] for (n,gs) in bus_gs)*w
    )
    _PMs.con(pm, nw, c, :kcl_q)[i] = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(qt[a_trans] for a_trans in bus_arcs_trans)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd*z_demand[n] for (n,qd) in bus_qd)
        + sum(bs*1.0^2*z_shunt[n] for (n,bs) in bus_bs)*w
    )
end


""
function constraint_mc_power_balance(pm::_PMs.AbstractWModels, nw::Int, c::Int, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    w    = _PMs.var(pm, nw, c, :w, i)
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
        - sum(gs for gs in values(bus_gs))*w
    )
    _PMs.con(pm, nw, c, :kcl_q)[i] = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(qt[a_trans] for a_trans in bus_arcs_trans)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd for qd in values(bus_qd))
        + sum(bs for bs in values(bus_bs))*w
    )
end


"delegate back to PowerModels"
function constraint_mc_ohms_yt_from(pm::_PMs.AbstractWModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    _PMs.constraint_ohms_yt_from(pm, n, c, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


"delegate back to PowerModels"
function constraint_mc_ohms_yt_to(pm::_PMs.AbstractWModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    _PMs.constraint_ohms_yt_to(pm, n, c, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end


"on/off bus voltage constraint for relaxed forms"
function constraint_mc_bus_voltage_on_off(pm::_PMs.AbstractWModels, n::Int, c::Int; kwargs...)
    for (i, bus) in _PMs.ref(pm, n, :bus)
        constraint_mc_voltage_magnitude_sqr_on_off(pm, i; nw=n, cnd=c)
    end
end


"By default, delegate back to PM; only certain formulations differ between PMD and PMs."
function constraint_mc_voltage_angle_difference(pm::_PMs.AbstractPowerModel, n::Int, f_idx, angmin, angmax)
    for c in _PMs.conductor_ids(pm; nw=n)
        _PMs.constraint_voltage_angle_difference(pm, n, c, f_idx, angmin[c], angmax[c])
    end
end
