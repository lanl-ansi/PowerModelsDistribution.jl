"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] == -b*(t[f_bus] - t[t_bus])
```
"""
function constraint_ohms_tp_yt_from(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractDCPForm
    p_fr  = var(pm, n, c,  :p, f_idx)
    va_fr = [var(pm, n, i, :va, f_bus) for i in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, i, :va, t_bus) for i in PMs.conductor_ids(pm)]

    @constraint(pm.model, p_fr == -sum(b[c,i]*(va_fr[c] - va_to[i]) for i in PMs.conductor_ids(pm)))
    # omit reactive constraint
end


"Do nothing, this model is symmetric"
function constraint_ohms_tp_yt_to(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractDCPForm
end


"`-b*(t[f_bus] - t[t_bus] + vad_min*(1-branch_z[i])) <= p[f_idx] <= -b*(t[f_bus] - t[t_bus] + vad_max*(1-branch_z[i]))`"
function constraint_ohms_tp_yt_from_on_off(pm::GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractDCPForm
    p_fr  = var(pm, n, c,  :p, f_idx)
    va_fr = [var(pm, n, i, :va, f_bus) for i in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, i, :va, t_bus) for i in PMs.conductor_ids(pm)]
    z = [var(pm, n, d, :branch_z, i) for d in PMs.conductor_ids(pm)]

    @constraint(pm.model, p_fr <= sum(-b[c,i]*(va_fr[i] - va_to[i] + vad_max[i]*(1-z[i])) for i in PMs.conductor_ids(pm)) )
    @constraint(pm.model, p_fr >= sum(-b[c,i]*(va_fr[i] - va_to[i] + vad_min[i]*(1-z[i])) for i in PMs.conductor_ids(pm)) )
end


"Do nothing, this model is symmetric"
function constraint_ohms_tp_yt_to_on_off(pm::GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractDCPForm
end
