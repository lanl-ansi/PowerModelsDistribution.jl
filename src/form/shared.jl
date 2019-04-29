# Three-phase specific constraints


""
function constraint_kcl_shunt_slack(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, i, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: PMs.AbstractWForms
    w    = PMs.var(pm, n, c, :w, i)
    p_slack = PMs.var(pm, n, c, :p_slack, i)
    q_slack = PMs.var(pm, n, c, :q_slack, i)
    pg   = PMs.var(pm, n, c, :pg)
    qg   = PMs.var(pm, n, c, :qg)
    p    = PMs.var(pm, n, c, :p)
    q    = PMs.var(pm, n, c, :q)
    p_dc = PMs.var(pm, n, c, :p_dc)
    q_dc = PMs.var(pm, n, c, :q_dc)

    JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*w + p_slack)
    JuMP.@constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*w + q_slack)
end



"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_tp_yt_from(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractWRForms
    p_fr = PMs.var(pm, n, c, :p, f_idx)
    q_fr = PMs.var(pm, n, c, :q, f_idx)
    w    = PMs.var(pm, n, :w)
    wr   = PMs.var(pm, n, :wr)
    wi   = PMs.var(pm, n, :wi)

    JuMP.@constraint(pm.model, p_fr ==  ( g_fr[c]+g[c,c]) * w[(f_bus, c)] +
                                sum( g[c,d] * wr[(f_bus, f_bus, c, d)] +
                                     b[c,d] * wi[(f_bus, f_bus, c, d)] for d in PMs.conductor_ids(pm) if d != c) +
                                sum(-g[c,d] * wr[(f_bus, t_bus, c, d)] +
                                    -b[c,d] * wi[(f_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm)) )
    JuMP.@constraint(pm.model, q_fr == -( b_fr[c]+b[c,c]) * w[(f_bus, c)] -
                                sum( b[c,d] * wr[(f_bus, f_bus, c, d)] -
                                     g[c,d] * wi[(f_bus, f_bus, c, d)] for d in PMs.conductor_ids(pm) if d != c) -
                                sum(-b[c,d] * wr[(f_bus, t_bus, c, d)] +
                                     g[c,d] * wi[(f_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm)) )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_tp_yt_to(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractWRForms
    q_to = PMs.var(pm, n, c, :q, t_idx)
    p_to = PMs.var(pm, n, c, :p, t_idx)
    w    = PMs.var(pm, n, :w)
    wr   = PMs.var(pm, n, :wr)
    wi   = PMs.var(pm, n, :wi)

    JuMP.@constraint(pm.model, p_to ==  ( g_to[c]+g[c,c]) * w[(t_bus, c)] +
                                sum( g[c,d] * wr[(t_bus, t_bus, c, d)] +
                                     b[c,d] *-wi[(t_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm) if d != c) +
                                sum(-g[c,d] * wr[(f_bus, t_bus, c, d)] +
                                    -b[c,d] *-wi[(f_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm)) )
    JuMP.@constraint(pm.model, q_to == -( b_to[c]+b[c,c]) * w[(t_bus, c)] -
                                sum( b[c,d] * wr[(t_bus, t_bus, c, d)] -
                                     g[c,d] *-wi[(t_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm) if d != c) -
                                sum(-b[c,d] * wr[(f_bus, t_bus, c, d)] +
                                     g[c,d] *-wi[(f_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm)) )
end


"do nothing, no way to represent this in these variables"
function constraint_tp_theta_ref(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, d) where T <: PMs.AbstractWForms
end


"Creates phase angle constraints at reference buses"
function constraint_tp_theta_ref(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, d) where T <: PMs.AbstractPForms
    va = PMs.var(pm, n, c, :va, d)
    nconductors = length(PMs.conductor_ids(pm))

    JuMP.@constraint(pm.model, va == wraptopi(2 * pi / nconductors * (1-c)))
end


"""
For a variable tap transformer, fix the tap variables which are fixed. For
example, an OLTC where the third phase is fixed, will have tap variables for
all phases, but the third tap variable should be fixed.
"""
function constraint_tp_oltc_tap_fix(pm::PMs.GenericPowerModel, i::Int, fixed::PMs.MultiConductorVector, tm::PMs.MultiConductorVector; nw=pm.cnw)
    for (c,fixed) in enumerate(fixed)
        if fixed
            JuMP.@constraint(pm.model, PMs.var(pm, nw, c, :tap)[i]==tm[c])
        end
    end
end
