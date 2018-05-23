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

    @constraint(pm.model, p_fr ==  g_fr[h]*w_fr[h] + sum(
                                    g[h,i]*w_fr[i] -
                                    g[h,i]*wr[i] -
                                    b[h,i]*wi[i] for i in PMs.phase_ids(pm)) )
    @constraint(pm.model, q_fr == -b_fr[h]*w_fr[h] + sum(
                                   -b[h,i]*w_fr[i] +
                                    b[h,i]*wr[i] -
                                    g[h,i]*wi[i] for i in PMs.phase_ids(pm)) )
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

    @constraint(pm.model, p_to ==  g_to[h]*w_to[h] + sum(
                                    g[h,i]*w_to[i] -
                                    g[h,i]*wr[i] -
                                    b[h,i]*-wi[i] for i in PMs.phase_ids(pm)) )
    @constraint(pm.model, q_to == -b_to[h]*w_to[h] + sum(
                                   -b[h,i]*w_to[i] +
                                    b[h,i]*wr[i] -
                                    g[h,i]*-wi[i] for i in PMs.phase_ids(pm)) )
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
