
"`vm[i] == vmref`"
function constraint_mc_voltage_magnitude_setpoint(pm::_PMs.AbstractWModels, n::Int, i::Int, vmref)
    w = _PMs.var(pm, n, :w, i)
    JuMP.@constraint(pm.model, w .== vmref.^2)
end

""
function constraint_mc_power_balance_slack(pm::_PMs.AbstractWModels, nw::Int, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    w    = _PMs.var(pm, nw, :w, i)
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
            - sum(gs[c] for gs in values(bus_gs))*w[c]
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
            + sum(bs[c] for bs in values(bus_bs))*w
            + q_slack[c]
        )
        push!(cstr_q, cq)

    end

    if _PMs.report_duals(pm)
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end




"Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)"
function constraint_mc_ohms_yt_from(pm::_PMs.AbstractWRModels, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    for c in _PMs.conductor_ids(pm, n)
        p_fr = _PMs.var(pm, n, :p, f_idx)[c]
        q_fr = _PMs.var(pm, n, :q, f_idx)[c]
        w    = _PMs.var(pm, n, :w)
        wr   = _PMs.var(pm, n, :wr)
        wi   = _PMs.var(pm, n, :wi)

        #TODO extend to shunt matrices; this ignores the off-diagonals
        JuMP.@constraint(pm.model, p_fr ==  ( g_fr[c,c]+g[c,c]) * w[f_bus][c] +
                                    sum( g[c,d] * wr[(f_bus, f_bus, c, d)] +
                                         b[c,d] * wi[(f_bus, f_bus, c, d)] for d in _PMs.conductor_ids(pm) if d != c) +
                                    sum(-g[c,d] * wr[(f_bus, t_bus, c, d)] +
                                        -b[c,d] * wi[(f_bus, t_bus, c, d)] for d in _PMs.conductor_ids(pm)) )
        JuMP.@constraint(pm.model, q_fr == -( b_fr[c,c]+b[c,c]) * w[f_bus][c] -
                                    sum( b[c,d] * wr[(f_bus, f_bus, c, d)] -
                                         g[c,d] * wi[(f_bus, f_bus, c, d)] for d in _PMs.conductor_ids(pm) if d != c) -
                                    sum(-b[c,d] * wr[(f_bus, t_bus, c, d)] +
                                         g[c,d] * wi[(f_bus, t_bus, c, d)] for d in _PMs.conductor_ids(pm)) )
        end
end


"Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)"
function constraint_mc_ohms_yt_to(pm::_PMs.AbstractWRModels, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    for c in _PMs.conductor_ids(pm, n)
        q_to = _PMs.var(pm, n, :q, t_idx)[c]
        p_to = _PMs.var(pm, n, :p, t_idx)[c]
        w    = _PMs.var(pm, n, :w)
        wr   = _PMs.var(pm, n, :wr)
        wi   = _PMs.var(pm, n, :wi)

        #TODO extend to shunt matrices; this ignores the off-diagonals
        JuMP.@constraint(pm.model, p_to ==  ( g_to[c,c]+g[c,c]) * w[t_bus][c] +
                                    sum( g[c,d] * wr[(t_bus, t_bus, c, d)] +
                                         b[c,d] *-wi[(t_bus, t_bus, c, d)] for d in _PMs.conductor_ids(pm) if d != c) +
                                    sum(-g[c,d] * wr[(f_bus, t_bus, c, d)] +
                                        -b[c,d] *-wi[(f_bus, t_bus, c, d)] for d in _PMs.conductor_ids(pm)) )
        JuMP.@constraint(pm.model, q_to == -( b_to[c,c]+b[c,c]) * w[t_bus][c] -
                                    sum( b[c,d] * wr[(t_bus, t_bus, c, d)] -
                                         g[c,d] *-wi[(t_bus, t_bus, c, d)] for d in _PMs.conductor_ids(pm) if d != c) -
                                    sum(-b[c,d] * wr[(f_bus, t_bus, c, d)] +
                                         g[c,d] *-wi[(f_bus, t_bus, c, d)] for d in _PMs.conductor_ids(pm)) )
    end
end


"do nothing, no way to represent this in these variables"
function constraint_mc_theta_ref(pm::_PMs.AbstractWModels, n::Int, d)
end


"Creates phase angle constraints at reference buses"
function constraint_mc_theta_ref(pm::_PMs.AbstractPolarModels, n::Int, d)
    cnds = _PMs.conductor_ids(pm; nw=n)
    nconductors = length(cnds)

    va = _PMs.var(pm, n, :va, d)

    for c in cnds
        JuMP.@constraint(pm.model, va[c] == _wrap_to_pi(2 * pi / nconductors * (1-c)))
    end
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
function constraint_mc_power_balance_shed(pm::_PMs.AbstractWModels, nw::Int, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    w    = _PMs.var(pm, nw, :w, i)
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

    _PMs.con(pm, nw, :kcl_p)[i] = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        .==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd .*z_demand[n] for (n,pd) in bus_pd)
        - sum(gs*1.0^2 .*z_shunt[n] for (n,gs) in bus_gs)*w
    )
    _PMs.con(pm, nw, :kcl_q)[i] = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(qt[a_trans] for a_trans in bus_arcs_trans)
        .==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd.*z_demand[n] for (n,qd) in bus_qd)
        + sum(bs*1.0^2  .*z_shunt[n] for (n,bs) in bus_bs)*w
    )
end


""
function constraint_mc_power_balance(pm::_PMs.AbstractWRModels, nw::Int, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    w = _PMs.var(pm, nw, :w, i)

    p = get(_PMs.var(pm, nw), :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q = get(_PMs.var(pm, nw), :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")

    psw  = get(_PMs.var(pm, nw),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(_PMs.var(pm, nw),  :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(_PMs.var(pm, nw),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(_PMs.var(pm, nw),   :qt, Dict()); _PMs._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")

    pg = get(_PMs.var(pm, nw), :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg = get(_PMs.var(pm, nw), :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PMs.var(pm, nw),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PMs.var(pm, nw),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")

    cstr_p = []
    cstr_q = []

    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        .==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        - sum(gs for gs in values(bus_gs))*w
    )

    cstr_q = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(qt[a_trans] for a_trans in bus_arcs_trans)
        .==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd for qd in values(bus_qd))
        + sum(bs for bs in values(bus_bs))*w
    )

    if _PMs.report_duals(pm)
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_mc_power_balance(pm::_PMs.AbstractWModels, nw::Int, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    Wr = _PMs.var(pm, nw, :Wr, i)
    # Wi = _PMs.var(pm, nw, :Wi, i)
    P = get(_PMs.var(pm, nw), :P, Dict()); _PMs._check_var_keys(P, bus_arcs, "active power", "branch")
    Q = get(_PMs.var(pm, nw), :Q, Dict()); _PMs._check_var_keys(Q, bus_arcs, "reactive power", "branch")
    Psw  = get(_PMs.var(pm, nw),  :Psw, Dict()); _PMs._check_var_keys(Psw, bus_arcs_sw, "active power", "switch")
    Qsw  = get(_PMs.var(pm, nw),  :Qsw, Dict()); _PMs._check_var_keys(Qsw, bus_arcs_sw, "reactive power", "switch")
    Pt   = get(_PMs.var(pm, nw),   :Pt, Dict()); _PMs._check_var_keys(Pt, bus_arcs_trans, "active power", "transformer")
    Qt   = get(_PMs.var(pm, nw),   :Qt, Dict()); _PMs._check_var_keys(Qt, bus_arcs_trans, "reactive power", "transformer")

    pg = get(_PMs.var(pm, nw), :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg = get(_PMs.var(pm, nw), :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PMs.var(pm, nw),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PMs.var(pm, nw),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")

    cstr_p = []
    cstr_q = []

    for c in _PMs.conductor_ids(pm; nw=nw)
        cp = JuMP.@constraint(pm.model,
            sum(P[a][c,c] for a in bus_arcs)
            + sum(Psw[a_sw][c,c] for a_sw in bus_arcs_sw)
            + sum(Pt[a_trans][c,c] for a_trans in bus_arcs_trans)
            ==
            sum(pg[g][c] for g in bus_gens)
            - sum(ps[s][c] for s in bus_storage)
            - sum(pd[c] for pd in values(bus_pd))
            - sum(gs[c] for gs in values(bus_gs))*Wr[c,c]
        )
        push!(cstr_p, cp)

        cq = JuMP.@constraint(pm.model,
            sum(Q[a][c,c] for a in bus_arcs)
            + sum(Qsw[a_sw][c,c] for a_sw in bus_arcs_sw)
            + sum(Qt[a_trans][c,c] for a_trans in bus_arcs_trans)
            ==
            sum(qg[g][c] for g in bus_gens)
            - sum(qs[s][c] for s in bus_storage)
            - sum(qd[c] for qd in values(bus_qd))
            + sum(bs[c] for bs in values(bus_bs))*Wr[c,c]
        )
        push!(cstr_q, cq)
    end

    if _PMs.report_duals(pm)
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
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
function constraint_mc_bus_voltage_on_off(pm::_PMs.AbstractWModels, n::Int; kwargs...)
    for (i, bus) in _PMs.ref(pm, n, :bus)
        constraint_mc_voltage_magnitude_sqr_on_off(pm, i, nw=n)
    end
end


function constraint_mc_voltage_angle_difference(pm::_PMs.AbstractPolarModels, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx

    va_fr = _PMs.var(pm, n, :va, f_bus)
    va_to = _PMs.var(pm, n, :va, t_bus)

    for c in _PMs.conductor_ids(pm; nw=n)
        JuMP.@constraint(pm.model, va_fr[c] - va_to[c] <= angmax[c])
        JuMP.@constraint(pm.model, va_fr[c] - va_to[c] >= angmin[c])
    end
end


""
function constraint_mc_voltage_angle_difference(pm::_PMs.AbstractWModels, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx

    ncnds = length(_PMs.conductor_ids(pm, n))

    w_fr = _PMs.var(pm, n, :w, f_bus)
    w_to = _PMs.var(pm, n, :w, t_bus)
    wr   = [_PMs.var(pm, n, :wr)[(f_bus, t_bus, c, c)] for c in 1:ncnds]
    wi   = [_PMs.var(pm, n, :wi)[(f_bus, t_bus, c, c)] for c in 1:ncnds]

    JuMP.@constraint(pm.model, wi .<= tan.(angmax).*wr)
    JuMP.@constraint(pm.model, wi .>= tan.(angmin).*wr)

    for c in 1:ncnds
        _PMs.cut_complex_product_and_angle_difference(pm.model, w_fr[c], w_to[c], wr[c], wi[c], angmin[c], angmax[c])
    end
end
