# Three-phase specific constraints


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
                                sum( g[c,i] * wr[(f_bus, f_bus, c, i)] +
                                     b[c,i] * wi[(f_bus, f_bus, c, i)] for i in PMs.conductor_ids(pm) if i != c) +
                                sum(-g[c,i] * wr[(f_bus, t_bus, c, i)] +
                                    -b[c,i] * wi[(f_bus, t_bus, c, i)] for i in PMs.conductor_ids(pm)) )
    @constraint(pm.model, q_fr == -( b_fr[c]+b[c,c]) * w[(f_bus, c)] -
                                sum( b[c,i] * wr[(f_bus, f_bus, c, i)] -
                                     g[c,i] * wi[(f_bus, f_bus, c, i)] for i in PMs.conductor_ids(pm) if i != c) -
                                sum(-b[c,i] * wr[(f_bus, t_bus, c, i)] +
                                     g[c,i] * wi[(f_bus, t_bus, c, i)] for i in PMs.conductor_ids(pm)) )
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
                                sum( g[c,i] * wr[(t_bus, t_bus, c, i)] +
                                     b[c,i] *-wi[(t_bus, t_bus, c, i)] for i in PMs.conductor_ids(pm) if i != c) +
                                sum(-g[c,i] * wr[(f_bus, t_bus, c, i)] +
                                    -b[c,i] *-wi[(f_bus, t_bus, c, i)] for i in PMs.conductor_ids(pm)) )
    @constraint(pm.model, q_to == -( b_to[c]+b[c,c]) * w[(t_bus, c)] -
                                sum( b[c,i] * wr[(t_bus, t_bus, c, i)] -
                                     g[c,i] *-wi[(t_bus, t_bus, c, i)] for i in PMs.conductor_ids(pm) if i != c) -
                                sum(-b[c,i] * wr[(f_bus, t_bus, c, i)] +
                                     g[c,i] *-wi[(f_bus, t_bus, c, i)] for i in PMs.conductor_ids(pm)) )
end


"do nothing, no way to represent this in these variables"
function constraint_tp_theta_ref(pm::GenericPowerModel{T}, n::Int, c::Int, i) where T <: PMs.AbstractWForms
end


"Creates phase angle constraints at reference buses"
function constraint_tp_theta_ref(pm::GenericPowerModel{T}, n::Int, c::Int, i) where T <: PMs.AbstractPForms
    va = var(pm, n, c, :va, i)
    num_conductors = length(PMs.conductor_ids(pm))

    @constraint(pm.model, va == 2 * pi / num_conductors * (c - 1))
end
