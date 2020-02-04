import LinearAlgebra: diag

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

    Gt = isempty(bus_gs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_gs))
    Bt = isempty(bus_bs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_bs))

    cstr_p = JuMP.@constraint(pm.model,
        sum(diag(P[a]) for a in bus_arcs)
        + sum(diag(Psw[a_sw]) for a_sw in bus_arcs_sw)
        + sum(diag(Pt[a_trans]) for a_trans in bus_arcs_trans)
        .==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        - diag(Wr*G'+Wi*B')
        + p_slack
    )

    cstr_q = JuMP.@constraint(pm.model,
        sum(diag(Q[a]) for a in bus_arcs)
        + sum(diag(Qsw[a_sw]) for a_sw in bus_arcs_sw)
        + sum(diag(Qt[a_trans]) for a_trans in bus_arcs_trans)
        .==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd for qd in values(bus_qd))
        - diag(-Wr*B'+Wi*G')
        + q_slack
    )

    if _PMs.report_duals(pm)
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


"do nothing, no way to represent this in these variables"
function constraint_mc_theta_ref(pm::_PMs.AbstractWModels, n::Int, d::Int, va_ref)
end


"Creates phase angle constraints at reference buses"
function constraint_mc_theta_ref(pm::_PMs.AbstractPolarModels, n::Int, d::Int, va_ref)
    cnds = _PMs.conductor_ids(pm; nw=n)
    nconductors = length(cnds)

    va = _PMs.var(pm, n, :va, d)

    JuMP.@constraint(pm.model, va .== va_ref)
end


"""
For a variable tap transformer, fix the tap variables which are fixed. For
example, an OLTC where the third phase is fixed, will have tap variables for
all phases, but the third tap variable should be fixed.
"""
function constraint_mc_oltc_tap_fix(pm::_PMs.AbstractPowerModel, i::Int, fixed::Vector, tm::Vector; nw=pm.cnw)
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

    bus_GsBs = [(n,bus_gs[n], bus_bs[n]) for n in keys(bus_gs)]

    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        .==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd .*z_demand[n] for (n,pd) in bus_pd)
        - sum(z_shunt[n].*diag(Wr*Gt'+Wi*Bt') for (n,Gs,Bs) in bus_GsBs)
    )
    cstr_q = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(qt[a_trans] for a_trans in bus_arcs_trans)
        .==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd.*z_demand[n] for (n,qd) in bus_qd)
        - diag(-Wr*Bt'+Wi*Gt')
        - sum(z_shunt[n].*diag(-Wr*Bt'+Wi*Gt') for (n,Gs,Bs) in bus_GsBs)
    )

    if _PMs.report_duals(pm)
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_mc_power_balance(pm::_PMs.AbstractWModels, nw::Int, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    Wr = _PMs.var(pm, nw, :Wr, i)
    Wi = _PMs.var(pm, nw, :Wi, i)
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

    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    Gt = isempty(bus_gs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_gs))
    Bt = isempty(bus_bs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_bs))

    cstr_p = JuMP.@constraint(pm.model,
        sum(diag(P[a]) for a in bus_arcs)
        + sum(diag(Psw[a_sw]) for a_sw in bus_arcs_sw)
        + sum(diag(Pt[a_trans]) for a_trans in bus_arcs_trans)
        .==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        - diag(Wr*Gt'+Wi*Bt')
    )

    cstr_q = JuMP.@constraint(pm.model,
        sum(diag(Q[a]) for a in bus_arcs)
        + sum(diag(Qsw[a_sw]) for a_sw in bus_arcs_sw)
        + sum(diag(Qt[a_trans]) for a_trans in bus_arcs_trans)
        .==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd for qd in values(bus_qd))
        - diag(-Wr*Bt'+Wi*Gt')
    )

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


function constraint_mc_storage_on_off(pm::_PMs.AbstractPowerModel, n::Int, i, pmin, pmax, qmin, qmax, charge_ub, discharge_ub)
    z_storage =_PMs.var(pm, n, :z_storage, i)
    ps =_PMs.var(pm, n, :ps, i)
    qs =_PMs.var(pm, n, :qs, i)

    JuMP.@constraint(pm.model, ps .<= z_storage.*pmax)
    JuMP.@constraint(pm.model, ps .>= z_storage.*pmin)

    JuMP.@constraint(pm.model, qs .<= z_storage.*qmax)
    JuMP.@constraint(pm.model, qs .>= z_storage.*qmin)
end
