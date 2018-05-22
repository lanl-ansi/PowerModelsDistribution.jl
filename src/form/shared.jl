# Three-phase specific constraints


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_tp_yt_from(pm::GenericPowerModel{T}, n::Int, h::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractWRForms
    p_fr = var(pm, n, h, :p, f_idx)
    q_fr = var(pm, n, h, :q, f_idx)
    w_fr = [var(pm, n, j, :w, f_bus) for j in PMs.phase_ids(pm)]
    wr   = [var(pm, n, j, :wr, (f_bus, t_bus)) for j in PMs.phase_ids(pm)]
    wi   = [var(pm, n, j, :wi, (f_bus, t_bus)) for j in PMs.phase_ids(pm)]

    @constraint(pm.model, p_fr ==  g_fr[h]/tm[h]^2*w_fr[h] + sum(g[h,i]/tm[h]^2*w_fr[i] for i in PMs.phase_ids(pm)) + sum((-g[h,i]*tr[h]+b[h,i]*ti[h])/tm[h]^2*wr[i] + (-b[h,i]*tr[h]-g[h,i]*ti[h])/tm[h]^2*wi[i] for i in PMs.phase_ids(pm)) )
    @constraint(pm.model, q_fr == -b_fr[h]/tm[h]^2*w_fr[h] - sum(b[h,i]/tm[h]^2*w_fr[i] for i in PMs.phase_ids(pm)) - sum((-b[h,i]*tr[h]-g[h,i]*ti[h])/tm[h]^2*wr[i] + (-g[h,i]*tr[h]+b[h,i]*ti[h])/tm[h]^2*wi[i] for i in PMs.phase_ids(pm)) )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_tp_yt_to(pm::GenericPowerModel{T}, n::Int, h::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractWRForms
    q_to = var(pm, n, h, :q, t_idx)
    p_to = var(pm, n, h, :p, t_idx)
    w_to = [var(pm, n, j, :w, t_bus) for j in PMs.phase_ids(pm)]
    wr   = [var(pm, n, j, :wr, (f_bus, t_bus)) for j in PMs.phase_ids(pm)]
    wi   = [var(pm, n, j, :wi, (f_bus, t_bus)) for j in PMs.phase_ids(pm)]

    @constraint(pm.model, p_to ==  g_to[h]*w_to[h] + sum(g[h,i]*w_to[i] for i in PMs.phase_ids(pm)) + sum((-g[h,i]*tr[h]-b[h,i]*ti[h])/tm[h]^2*wr[i] + (-b[h,i]*tr[h]+g[h,i]*ti[h])/tm[h]^2*-wi[i] for i in PMs.phase_ids(pm)) )
    @constraint(pm.model, q_to == -b_to[h]*w_to[h] - sum(b[h,i]*w_to[i] for i in PMs.phase_ids(pm)) - sum((-b[h,i]*tr[h]+g[h,i]*ti[h])/tm[h]^2*wr[i] + (-g[h,i]*tr[h]-b[h,i]*ti[h])/tm[h]^2*-wi[i] for i in PMs.phase_ids(pm)) )
end


"do nothing, no way to represent this in these variables"
function constraint_tp_theta_ref(pm::GenericPowerModel{T}, n::Int, h::Int, i) where T <: PMs.AbstractWRForms
end