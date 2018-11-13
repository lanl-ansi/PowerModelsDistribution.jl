
"nothing to do, these models do not have complex voltage constraints"
function constraint_tp_voltage(pm::GenericPowerModel{T}, n::Int, c::Int) where T <: PMs.AbstractDCPForm
end


### DC Power Flow Approximation ###

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] == -b*(t[f_bus] - t[t_bus])
```
"""
function constraint_ohms_tp_yt_from(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractDCPForm
    p_fr  = var(pm, n, c,  :p, f_idx)
    va_fr = [var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]

    @constraint(pm.model, p_fr == -sum(b[c,d]*(va_fr[c] - va_to[d]) for d in PMs.conductor_ids(pm)))
    # omit reactive constraint
end


"Do nothing, this model is symmetric"
function constraint_ohms_tp_yt_to(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractDCPForm
end


"`-b*(t[f_bus] - t[t_bus] + vad_min*(1-branch_z[i])) <= p[f_idx] <= -b*(t[f_bus] - t[t_bus] + vad_max*(1-branch_z[i]))`"
function constraint_ohms_tp_yt_from_on_off(pm::GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractDCPForm
    p_fr  = var(pm, n, c,  :p, f_idx)
    va_fr = [var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]
    z = [var(pm, n, d, :branch_z, i) for d in PMs.conductor_ids(pm)]

    @constraint(pm.model, p_fr <= sum(-b[c,d]*(va_fr[d] - va_to[d] + vad_max[d]*(1-z[d])) for d in PMs.conductor_ids(pm)) )
    @constraint(pm.model, p_fr >= sum(-b[c,d]*(va_fr[d] - va_to[d] + vad_min[d]*(1-z[d])) for d in PMs.conductor_ids(pm)) )
end


"Do nothing, this model is symmetric"
function constraint_ohms_tp_yt_to_on_off(pm::GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractDCPForm
end



function constraint_tp_storage_loss(pm::GenericPowerModel{T}, n::Int, i, bus, r, x, standby_loss) where T <: PMs.AbstractDCPForm
    conductors = PMs.conductor_ids(pm)
    ps = [var(pm, n, c, :ps, i) for c in conductors]
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)

    @NLconstraint(pm.model, sum(ps[c] for c in conductors) + (sd - sc) == standby_loss + sum( r[c]*ps[c]^2 for c in conductors) )
end



### Network Flow Approximation ###

"nothing to do, no voltage angle variables"
function constraint_tp_theta_ref(pm::GenericPowerModel{T}, n::Int, c::Int, d) where T <: PMs.NFAForm
end

"nothing to do, no voltage angle variables"
function constraint_ohms_tp_yt_from(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.NFAForm
end

"nothing to do, this model is symmetric"
function constraint_ohms_tp_yt_to(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.NFAForm
end

