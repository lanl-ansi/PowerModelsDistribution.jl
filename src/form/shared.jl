# Three-phase specific constraints


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_tp_yt_from(pm::GenericPowerModel{T}, n::Int, h::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractWRForms
    p_fr = var(pm, n, h, :p, f_idx)
    q_fr = var(pm, n, h, :q, f_idx)
    w_fr = var(pm, n, h, :w, f_bus)
    wr   = var(pm, n, h, :wr, (f_bus, t_bus))
    wi   = var(pm, n, h, :wi, (f_bus, t_bus))

    @constraint(pm.model, p_fr ==  (g[h,h]+g_fr)/tm^2*w_fr + (-g*tr+b*ti)[h]/tm^2*wr + (-b*tr-g*ti)[h]/tm^2*wi )
    @constraint(pm.model, q_fr == -(b[h,h]+b_fr)/tm^2*w_fr - (-b*tr-g*ti)[h]/tm^2*wr + (-g*tr+b*ti)[h]/tm^2*wi )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_tp_yt_to(pm::GenericPowerModel{T}, n::Int, h::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractWRForms
    q_to = var(pm, n, h, :q, t_idx)
    p_to = var(pm, n, h, :p, t_idx)
    w_to = var(pm, n, h, :w, t_bus)
    wr   = var(pm, n, h, :wr, (f_bus, t_bus))
    wi   = var(pm, n, h, :wi, (f_bus, t_bus))

    @constraint(pm.model, p_to ==  (g[h,h]+g_to)*w_to + (-g*tr-b*ti)[h]/tm^2*wr + (-b*tr+g*ti)[h]/tm^2*-wi )
    @constraint(pm.model, q_to == -(b[h,h]+b_to)*w_to - (-b*tr+g*ti)[h]/tm^2*wr + (-g*tr-b*ti)[h]/tm^2*-wi )
end