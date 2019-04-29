
"nothing to do, these models do not have complex voltage constraints"
function constraint_tp_voltage(pm::PMs.GenericPowerModel{T}, n::Int, c::Int) where T <: PMs.AbstractDCPForm
end


### DC Power Flow Approximation ###

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] == -b*(t[f_bus] - t[t_bus])
```
"""
function constraint_ohms_tp_yt_from(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractDCPForm
    p_fr  = PMs.var(pm, n, c,  :p, f_idx)
    va_fr = [PMs.var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [PMs.var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]

    JuMP.@constraint(pm.model, p_fr == -sum(b[c,d]*(va_fr[c] - va_to[d]) for d in PMs.conductor_ids(pm)))
    # omit reactive constraint
end


"Do nothing, this model is symmetric"
function constraint_ohms_tp_yt_to(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractDCPForm
end


"`-b*(t[f_bus] - t[t_bus] + vad_min*(1-branch_z[i])) <= p[f_idx] <= -b*(t[f_bus] - t[t_bus] + vad_max*(1-branch_z[i]))`"
function constraint_ohms_tp_yt_from_on_off(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractDCPForm
    p_fr  = PMs.var(pm, n, c,  :p, f_idx)
    va_fr = [PMs.var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [PMs.var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]
    z = [PMs.var(pm, n, d, :branch_z, i) for d in PMs.conductor_ids(pm)]

    JuMP.@constraint(pm.model, p_fr <= sum(-b[c,d]*(va_fr[d] - va_to[d] + vad_max[d]*(1-z[d])) for d in PMs.conductor_ids(pm)) )
    JuMP.@constraint(pm.model, p_fr >= sum(-b[c,d]*(va_fr[d] - va_to[d] + vad_min[d]*(1-z[d])) for d in PMs.conductor_ids(pm)) )
end


"Do nothing, this model is symmetric"
function constraint_ohms_tp_yt_to_on_off(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractDCPForm
end



function constraint_tp_storage_loss(pm::PMs.GenericPowerModel{T}, n::Int, i, bus, r, x, standby_loss) where T <: PMs.AbstractDCPForm
    conductors = PMs.conductor_ids(pm)
    ps = [PMs.var(pm, n, c, :ps, i) for c in conductors]
    sc = PMs.var(pm, n, :sc, i)
    sd = PMs.var(pm, n, :sd, i)

    JuMP.@NLconstraint(pm.model, sum(ps[c] for c in conductors) + (sd - sc) == standby_loss + sum( r[c]*ps[c]^2 for c in conductors) )
end


### Network Flow Approximation ###

"nothing to do, no voltage angle variables"
function constraint_tp_theta_ref(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, d) where T <: PMs.NFAForm
end

"nothing to do, no voltage angle variables"
function constraint_ohms_tp_yt_from(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.NFAForm
end

"nothing to do, this model is symmetric"
function constraint_ohms_tp_yt_to(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.NFAForm
end


""
function constraint_kcl_shunt_trans(pm::PMs.GenericPowerModel{T}, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: PMs.AbstractDCPForm
    pg   = PMs.var(pm, nw, c, :pg)
    p    = PMs.var(pm, nw, c, :p)
    p_dc = PMs.var(pm, nw, c, :p_dc)
    p_trans = PMs.var(pm, nw, c, :pt)

    PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) + sum(p_trans[a_trans] for a_trans in bus_arcs_trans) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*1.0^2)
    # omit reactive constraint
end
