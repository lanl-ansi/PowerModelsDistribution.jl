# Three-phase specific constraints


""
function constraint_kcl_shunt_slack(pm::GenericPowerModel{T}, n::Int, c::Int, i, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: PMs.AbstractWForms
    w    = var(pm, n, c, :w, i)
    p_slack = var(pm, n, c, :p_slack, i)
    q_slack = var(pm, n, c, :q_slack, i)
    pg   = var(pm, n, c, :pg)
    qg   = var(pm, n, c, :qg)
    p    = var(pm, n, c, :p)
    q    = var(pm, n, c, :q)
    p_dc = var(pm, n, c, :p_dc)
    q_dc = var(pm, n, c, :q_dc)

    @constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*w + p_slack)
    @constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*w + q_slack)
end



"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_tp_yt_from(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractWRForms
    p_fr = var(pm, n, c, :p, f_idx)
    q_fr = var(pm, n, c, :q, f_idx)
    w    = var(pm, n, :w)
    wr   = var(pm, n, :wr)
    wi   = var(pm, n, :wi)

    @constraint(pm.model, p_fr ==  ( g_fr[c]+g[c,c]) * w[(f_bus, c)] +
                                sum( g[c,d] * wr[(f_bus, f_bus, c, d)] +
                                     b[c,d] * wi[(f_bus, f_bus, c, d)] for d in PMs.conductor_ids(pm) if d != c) +
                                sum(-g[c,d] * wr[(f_bus, t_bus, c, d)] +
                                    -b[c,d] * wi[(f_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm)) )
    @constraint(pm.model, q_fr == -( b_fr[c]+b[c,c]) * w[(f_bus, c)] -
                                sum( b[c,d] * wr[(f_bus, f_bus, c, d)] -
                                     g[c,d] * wi[(f_bus, f_bus, c, d)] for d in PMs.conductor_ids(pm) if d != c) -
                                sum(-b[c,d] * wr[(f_bus, t_bus, c, d)] +
                                     g[c,d] * wi[(f_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm)) )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_tp_yt_to(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractWRForms
    q_to = var(pm, n, c, :q, t_idx)
    p_to = var(pm, n, c, :p, t_idx)
    w    = var(pm, n, :w)
    wr   = var(pm, n, :wr)
    wi   = var(pm, n, :wi)

    @constraint(pm.model, p_to ==  ( g_to[c]+g[c,c]) * w[(t_bus, c)] +
                                sum( g[c,d] * wr[(t_bus, t_bus, c, d)] +
                                     b[c,d] *-wi[(t_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm) if d != c) +
                                sum(-g[c,d] * wr[(f_bus, t_bus, c, d)] +
                                    -b[c,d] *-wi[(f_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm)) )
    @constraint(pm.model, q_to == -( b_to[c]+b[c,c]) * w[(t_bus, c)] -
                                sum( b[c,d] * wr[(t_bus, t_bus, c, d)] -
                                     g[c,d] *-wi[(t_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm) if d != c) -
                                sum(-b[c,d] * wr[(f_bus, t_bus, c, d)] +
                                     g[c,d] *-wi[(f_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm)) )
end


"do nothing, no way to represent this in these variables"
function constraint_tp_theta_ref(pm::GenericPowerModel{T}, n::Int, c::Int, d) where T <: PMs.AbstractWForms
end


"Creates phase angle constraints at reference buses"
function constraint_tp_theta_ref(pm::GenericPowerModel{T}, n::Int, c::Int, d) where T <: PMs.AbstractPForms
    va = var(pm, n, c, :va, d)
    nconductors = length(PMs.conductor_ids(pm))

    @constraint(pm.model, va == wraptopi(2 * pi / nconductors * (1-c)))
end
