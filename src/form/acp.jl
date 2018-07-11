# Three-phase specific constraints
""
function variable_tp_voltage(pm::GenericPowerModel{T}; kwargs...) where T <: PMs.AbstractACPForm
    PMs.variable_voltage_angle(pm; kwargs...)
    PMs.variable_voltage_magnitude(pm; kwargs...)
end


"do nothing, this model does not have complex voltage constraints"
function constraint_tp_voltage(pm::GenericPowerModel{T}; kwargs...) where T <: PMs.AbstractACPForm
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] ==  (g+g_fr)/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
q[f_idx] == -(b+b_fr)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
```
"""
function constraint_ohms_tp_yt_from(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractACPForm
    p_fr  = var(pm, n, c,  :p, f_idx)
    q_fr  = var(pm, n, c,  :q, f_idx)
    vm_fr = [var(pm, n, d, :vm, f_bus) for d in PMs.conductor_ids(pm)]
    vm_to = [var(pm, n, d, :vm, t_bus) for d in PMs.conductor_ids(pm)]
    va_fr = [var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]

    @NLconstraint(pm.model, p_fr ==  (g_fr[c]+g[c,c]) * vm_fr[c]^2 +
                                    sum( g[c,i]*vm_fr[c]*vm_fr[i]*cos(va_fr[c]-va_fr[i]) +
                                         b[c,i]*vm_fr[c]*vm_fr[i]*sin(va_fr[c]-va_fr[i]) for i in PMs.conductor_ids(pm) if i != c) +
                                    sum(-g[c,i]*vm_fr[c]*vm_to[i]*cos(va_fr[c]-va_to[i]) +
                                        -b[c,i]*vm_fr[c]*vm_to[i]*sin(va_fr[c]-va_to[i]) for i in PMs.conductor_ids(pm)) )
    @NLconstraint(pm.model, q_fr == -(b_fr[c]+b[c,c]) *vm_fr[c]^2 -
                                    sum( b[c,i]*vm_fr[c]*vm_fr[i]*cos(va_fr[c]-va_fr[i]) -
                                         g[c,i]*vm_fr[c]*vm_fr[i]*sin(va_fr[c]-va_fr[i]) for i in PMs.conductor_ids(pm) if i != c) -
                                    sum(-b[c,i]*vm_fr[c]*vm_to[i]*cos(va_fr[c]-va_to[i]) +
                                         g[c,i]*vm_fr[c]*vm_to[i]*sin(va_fr[c]-va_to[i]) for i in PMs.conductor_ids(pm)) )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_ohms_tp_yt_to(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractACPForm
    p_to  = var(pm, n, c,  :p, t_idx)
    q_to  = var(pm, n, c,  :q, t_idx)
    vm_fr = [var(pm, n, i, :vm, f_bus) for i in PMs.conductor_ids(pm)]
    vm_to = [var(pm, n, i, :vm, t_bus) for i in PMs.conductor_ids(pm)]
    va_fr = [var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]

    @NLconstraint(pm.model, p_to == (g_to[c]+g[c,c])*vm_to[c]^2 +
                                   sum( g[c,i]*vm_to[c]*vm_to[i]*cos(va_to[c]-va_to[i]) +
                                        b[c,i]*vm_to[c]*vm_to[i]*sin(va_to[c]-va_to[i]) for i in PMs.conductor_ids(pm) if i != c) +
                                   sum(-g[c,i]*vm_to[c]*vm_fr[i]*cos(va_to[c]-va_fr[i]) +
                                       -b[c,i]*vm_to[c]*vm_fr[i]*sin(va_to[c]-va_fr[i]) for i in PMs.conductor_ids(pm)) )
    @NLconstraint(pm.model, q_to == -(b_to[c]+b[c,c])*vm_to[c]^2 -
                                   sum( b[c,i]*vm_to[c]*vm_to[i]*cos(va_to[c]-va_to[i]) -
                                        g[c,i]*vm_to[c]*vm_to[i]*sin(va_to[c]-va_to[i]) for i in PMs.conductor_ids(pm) if i != c) -
                                   sum(-b[c,i]*vm_to[c]*vm_fr[i]*cos(va_to[c]-va_fr[i]) +
                                        g[c,i]*vm_to[c]*vm_fr[i]*sin(va_to[c]-va_fr[i]) for i in PMs.conductor_ids(pm)) )
end


"""
```
p[f_idx] == z*(g/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
q[f_idx] == z*(-(b+c/2)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
```
"""
function constraint_ohms_tp_yt_from_on_off(pm::GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractACPForm
    p_fr  = var(pm, n, c,  :p, f_idx)
    q_fr  = var(pm, n, c,  :q, f_idx)
    vm_fr = [var(pm, n, i, :vm, f_bus) for i in PMs.conductor_ids(pm)]
    vm_to = [var(pm, n, i, :vm, t_bus) for i in PMs.conductor_ids(pm)]
    va_fr = [var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]
    z = var(pm, n, c, :branch_z, i)

    @NLconstraint(pm.model, p_fr == z*( g_fr[c]*vm_fr[c]^2 + sum(
                                            g[c,i]*vm_fr[c]*vm_fr[i]*cos(va_fr[c]-va_fr[i]) +
                                            b[c,i]*vm_fr[c]*vm_fr[i]*sin(va_fr[c]-va_fr[i]) -
                                            g[c,i]*vm_fr[c]*vm_to[i]*cos(va_fr[c]-va_to[i]) -
                                            b[c,i]*vm_fr[c]*vm_to[i]*sin(va_fr[c]-va_to[i]) for i in PMs.conductor_ids(pm)) ) )
    @NLconstraint(pm.model, q_fr == z*(-b_fr[c]*vm_fr[c]^2 - sum(
                                            b[c,i]*vm_fr[c]*vm_fr[i]*cos(va_fr[c]-va_fr[i]) -
                                            g[c,i]*vm_fr[c]*vm_fr[i]*sin(va_fr[c]-va_fr[i]) -
                                            b[c,i]*vm_fr[c]*vm_to[i]*cos(va_fr[c]-va_to[i]) +
                                            g[c,i]*vm_fr[c]*vm_to[i]*sin(va_fr[c]-va_to[i]) for i in PMs.conductor_ids(pm)) ) )
end

"""
```
p[t_idx] == z*(g*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
q[t_idx] == z*(-(b+c/2)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
```
"""
function constraint_ohms_tp_yt_to_on_off(pm::GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractACPForm
    p_to  = var(pm, n, c,  :p, t_idx)
    q_to  = var(pm, n, c,  :q, t_idx)
    vm_fr = [var(pm, n, i, :vm, f_bus) for i in PMs.conductor_ids(pm)]
    vm_to = [var(pm, n, i, :vm, t_bus) for i in PMs.conductor_ids(pm)]
    va_fr = [var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]
    z = var(pm, n, c, :branch_z, i)

    @NLconstraint(pm.model, p_to == z*( g_to[c]*vm_to[c]^2 + sum(
                                        g[c,i]*vm_to[c]*vm_to[i]*cos(va_to[c]-va_to[i]) +
                                        b[c,i]*vm_to[c]*vm_to[i]*sin(va_to[c]-va_to[i]) -
                                        g[c,i]*vm_to[c]*vm_fr[i]*cos(va_to[c]-va_fr[i]) -
                                        b[c,i]*vm_to[c]*vm_fr[i]*sin(va_to[c]-va_fr[i]) for i in PMs.conductor_ids(pm)) ) )
    @NLconstraint(pm.model, q_to == z*(-b_to[c]*vm_to[c]^2 - sum(
                                        b[c,i]*vm_to[c]*vm_to[i]*cos(va_to[c]-va_to[i]) -
                                        g[c,i]*vm_to[c]*vm_to[i]*sin(va_to[c]-va_to[i]) -
                                        b[c,i]*vm_to[c]*vm_fr[i]*cos(va_to[c]-va_fr[i]) +
                                        g[c,i]*vm_to[c]*vm_fr[i]*sin(va_to[c]-va_fr[i]) for i in PMs.conductor_ids(pm)) ) )
end
