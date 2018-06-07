# Three-phase specific constraints


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_tp_yt_from(pm::GenericPowerModel{T}, n::Int, h::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractWRForms
    p_fr = var(pm, n, h, :p, f_idx)
    q_fr = var(pm, n, h, :q, f_idx)
    w    = var(pm, n, :w)
    wr   = var(pm, n, :wr)
    wi   = var(pm, n, :wi)

    @constraint(pm.model, p_fr ==  ( g_fr[h]+g[h,h]) * w[(f_bus, h)] +
                                sum( g[h,i] * wr[(f_bus, f_bus, h, i)] +
                                     b[h,i] * wi[(f_bus, f_bus, h, i)] for i in PMs.phase_ids(pm) if i != h) +
                                sum(-g[h,i] * wr[(f_bus, t_bus, h, i)] +
                                    -b[h,i] * wi[(f_bus, t_bus, h, i)] for i in PMs.phase_ids(pm)) )
    @constraint(pm.model, q_fr == -( b_fr[h]+b[h,h]) * w[(f_bus, h)] -
                                sum( b[h,i] * wr[(f_bus, f_bus, h, i)] -
                                     g[h,i] * wi[(f_bus, f_bus, h, i)] for i in PMs.phase_ids(pm) if i != h) -
                                sum(-b[h,i] * wr[(f_bus, t_bus, h, i)] +
                                     g[h,i] * wi[(f_bus, t_bus, h, i)] for i in PMs.phase_ids(pm)) )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_tp_yt_to(pm::GenericPowerModel{T}, n::Int, h::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractWRForms
    q_to = var(pm, n, h, :q, t_idx)
    p_to = var(pm, n, h, :p, t_idx)
    w    = var(pm, n, :w)
    wr   = var(pm, n, :wr)
    wi   = var(pm, n, :wi)

    @constraint(pm.model, p_to ==  ( g_to[h]+g[h,h]) * w[(t_bus, h)] +
                                sum( g[h,i] * wr[(t_bus, t_bus, h, i)] +
                                     b[h,i] *-wi[(t_bus, t_bus, h, i)] for i in PMs.phase_ids(pm) if i != h) +
                                sum(-g[h,i] * wr[(f_bus, t_bus, h, i)] +
                                    -b[h,i] *-wi[(f_bus, t_bus, h, i)] for i in PMs.phase_ids(pm)) )
    @constraint(pm.model, q_to == -( b_to[h]+b[h,h]) * w[(t_bus, h)] -
                                sum( b[h,i] * wr[(t_bus, t_bus, h, i)] -
                                     g[h,i] *-wi[(t_bus, t_bus, h, i)] for i in PMs.phase_ids(pm) if i != h) -
                                sum(-b[h,i] * wr[(f_bus, t_bus, h, i)] +
                                     g[h,i] *-wi[(f_bus, t_bus, h, i)] for i in PMs.phase_ids(pm)) )
end


"do nothing, no way to represent this in these variables"
function constraint_tp_theta_ref(pm::GenericPowerModel{T}, n::Int, h::Int, i) where T <: PMs.AbstractWForms
end


"Creates phase angle constraints at reference buses"
function constraint_tp_theta_ref(pm::GenericPowerModel{T}, n::Int, h::Int, i) where T <: PMs.AbstractPForms
    va = var(pm, n, h, :va, i)
    nphases = length(PMs.phase_ids(pm))

    @constraint(pm.model, va == 2 * pi / nphases * (h - 1))
end
