""
function constraint_tp_power_balance_shunt_slack(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, i, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: _PMs.AbstractWForms
    w    = _PMs.var(pm, n, c, :w, i)
    p_slack = _PMs.var(pm, n, c, :p_slack, i)
    q_slack = _PMs.var(pm, n, c, :q_slack, i)
    pg   = _PMs.var(pm, n, c, :pg)
    qg   = _PMs.var(pm, n, c, :qg)
    p    = _PMs.var(pm, n, c, :p)
    q    = _PMs.var(pm, n, c, :q)
    p_dc = _PMs.var(pm, n, c, :p_dc)
    q_dc = _PMs.var(pm, n, c, :q_dc)

    JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*w + p_slack)
    JuMP.@constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*w + q_slack)
end


"Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)"
function constraint_tp_ohms_yt_from(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: _PMs.AbstractWRForms
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
function constraint_tp_ohms_yt_to(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: _PMs.AbstractWRForms
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
function constraint_tp_theta_ref(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, d) where T <: _PMs.AbstractWForms
end


"Creates phase angle constraints at reference buses"
function constraint_tp_theta_ref(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, d) where T <: _PMs.AbstractPForms
    va = _PMs.var(pm, n, c, :va, d)
    nconductors = length(_PMs.conductor_ids(pm))

    JuMP.@constraint(pm.model, va == _wrap_to_pi(2 * pi / nconductors * (1-c)))
end


"""
For a variable tap transformer, fix the tap variables which are fixed. For
example, an OLTC where the third phase is fixed, will have tap variables for
all phases, but the third tap variable should be fixed.
"""
function constraint_tp_oltc_tap_fix(pm::_PMs.GenericPowerModel, i::Int, fixed::_PMs.MultiConductorVector, tm::_PMs.MultiConductorVector; nw=pm.cnw)
    for (c,fixed) in enumerate(fixed)
        if fixed
            JuMP.@constraint(pm.model, _PMs.var(pm, nw, c, :tap)[i]==tm[c])
        end
    end
end


"KCL for load shed problem with transformers (AbstractWForms)"
function constraint_tp_power_balance_shunt_trans_shed(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, i, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: _PMs.AbstractWForms
    w    = _PMs.var(pm, n, c, :w, i)
    pg   = _PMs.var(pm, n, c, :pg)
    qg   = _PMs.var(pm, n, c, :qg)
    p    = _PMs.var(pm, n, c, :p)
    q    = _PMs.var(pm, n, c, :q)
    p_dc = _PMs.var(pm, n, c, :p_dc)
    q_dc = _PMs.var(pm, n, c, :q_dc)
    pt = _PMs.var(pm, n, c, :pt)
    qt = _PMs.var(pm, n, c, :qt)
    z_demand = _PMs.var(pm, n, :z_demand)
    z_shunt = _PMs.var(pm, n, :z_shunt)

    JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) + sum(pt[a_trans] for a_trans in bus_arcs_trans) == sum(pg[g] for g in bus_gens) - sum(pd*z_demand[i] for (i,pd) in bus_pd) - sum(gs*1.0^2*z_shunt[i] for (i,gs) in bus_gs)*w)
    JuMP.@constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) + sum(qt[a_trans] for a_trans in bus_arcs_trans) == sum(qg[g] for g in bus_gens) - sum(qd*z_demand[i] for (i,qd) in bus_qd) + sum(bs*1.0^2*z_shunt[i] for (i,bs) in bus_bs)*w)
end


"delegate back to PowerModels"
function constraint_tp_ohms_yt_from(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: _PMs.AbstractWForms
    _PMs.constraint_tp_ohms_yt_from(pm, n, c, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


"delegate back to PowerModels"
function constraint_tp_ohms_yt_to(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: _PMs.AbstractWForms
    _PMs.constraint_ohms_yt_to(pm, n, c, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end


"on/off bus voltage constraint for relaxed forms"
function constraint_tp_bus_voltage_on_off(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int; kwargs...) where T <: _PMs.AbstractWForms
    for (i, bus) in _PMs.ref(pm, n, :bus)
        constraint_tp_voltage_magnitude_sqr_on_off(pm, i; nw=n, cnd=c)
    end
end
