# Three-phase specific constraints


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] ==  (g+g_fr)/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
q[f_idx] == -(b+b_fr)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
```
"""
function constraint_ohms_tp_yt_from(pm::GenericPowerModel{T}, n::Int, h::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractACPForm
    p_fr  = var(pm, n, h,  :p, f_idx)
    q_fr  = var(pm, n, h,  :q, f_idx)
    vm_fr = [var(pm, n, j, :vm, f_bus) for j in PMs.phase_ids(pm)]
    vm_to = [var(pm, n, j, :vm, t_bus) for j in PMs.phase_ids(pm)]
    va_fr = [var(pm, n, j, :va, f_bus) for j in PMs.phase_ids(pm)]
    va_to = [var(pm, n, j, :va, t_bus) for j in PMs.phase_ids(pm)]

    @NLconstraint(pm.model, p_fr ==  g_fr[h]*vm_fr[h]^2 + sum(
                                        g[h,i]*vm_fr[h]*vm_fr[i]*cos(va_fr[h]-va_fr[i]) +
                                        b[h,i]*vm_fr[h]*vm_fr[i]*sin(va_fr[h]-va_fr[i]) -
                                        g[h,i]*vm_fr[h]*vm_to[i]*cos(va_fr[h]-va_to[i]) -
                                        b[h,i]*vm_fr[h]*vm_to[i]*sin(va_fr[h]-va_to[i]) for i in PMs.phase_ids(pm)) )
    @NLconstraint(pm.model, q_fr == -b_fr[h]*vm_fr[h]^2 - sum(
                                        b[h,i]*vm_fr[h]*vm_fr[i]*cos(va_fr[h]-va_fr[i]) -
                                        g[h,i]*vm_fr[h]*vm_fr[i]*sin(va_fr[h]-va_fr[i]) -
                                        b[h,i]*vm_fr[h]*vm_to[i]*cos(va_fr[h]-va_to[i]) +
                                        g[h,i]*vm_fr[h]*vm_to[i]*sin(va_fr[h]-va_to[i]) for i in PMs.phase_ids(pm)) )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_ohms_tp_yt_to(pm::GenericPowerModel{T}, n::Int, h::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractACPForm
    p_to  = var(pm, n, h,  :p, t_idx)
    q_to  = var(pm, n, h,  :q, t_idx)
    vm_fr = [var(pm, n, i, :vm, f_bus) for i in PMs.phase_ids(pm)]
    vm_to = [var(pm, n, i, :vm, t_bus) for i in PMs.phase_ids(pm)]
    va_fr = [var(pm, n, j, :va, f_bus) for j in PMs.phase_ids(pm)]
    va_to = [var(pm, n, j, :va, t_bus) for j in PMs.phase_ids(pm)]

    @NLconstraint(pm.model, p_to ==  g_to[h]*vm_to[h]^2 + sum(
                                        g[h,i]*vm_to[h]*vm_to[i]*cos(va_to[h]-va_to[i]) +
                                        b[h,i]*vm_to[h]*vm_to[i]*sin(va_to[h]-va_to[i]) -
                                        g[h,i]*vm_to[h]*vm_fr[i]*cos(va_to[h]-va_fr[i]) -
                                        b[h,i]*vm_to[h]*vm_fr[i]*sin(va_to[h]-va_fr[i]) for i in PMs.phase_ids(pm)) )
    @NLconstraint(pm.model, q_to == -b_to[h]*vm_to[h]^2 - sum(
                                        b[h,i]*vm_to[h]*vm_to[i]*cos(va_to[h]-va_to[i]) -
                                        g[h,i]*vm_to[h]*vm_to[i]*sin(va_to[h]-va_to[i]) -
                                        b[h,i]*vm_to[h]*vm_fr[i]*cos(va_to[h]-va_fr[i]) +
                                        g[h,i]*vm_to[h]*vm_fr[i]*sin(va_to[h]-va_fr[i]) for i in PMs.phase_ids(pm)) )
end


"""
```
p[f_idx] == z*(g/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
q[f_idx] == z*(-(b+c/2)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
```
"""
function constraint_ohms_tp_yt_from_on_off(pm::GenericPowerModel{T}, n::Int, h::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractACPForm
    p_fr  = var(pm, n, h,  :p, f_idx)
    q_fr  = var(pm, n, h,  :q, f_idx)
    vm_fr = [var(pm, n, i, :vm, f_bus) for i in PMs.phase_ids(pm)]
    vm_to = [var(pm, n, i, :vm, t_bus) for i in PMs.phase_ids(pm)]
    va_fr = [var(pm, n, j, :va, f_bus) for j in PMs.phase_ids(pm)]
    va_to = [var(pm, n, j, :va, t_bus) for j in PMs.phase_ids(pm)]
    z = var(pm, n, h, :branch_z, i)

    @NLconstraint(pm.model, p_fr == z*( g_fr[h]*vm_fr[h]^2 + sum(
                                            g[h,i]*vm_fr[h]*vm_fr[i]*cos(va_fr[h]-va_fr[i]) +
                                            b[h,i]*vm_fr[h]*vm_fr[i]*sin(va_fr[h]-va_fr[i]) -
                                            g[h,i]*vm_fr[h]*vm_to[i]*cos(va_fr[h]-va_to[i]) -
                                            b[h,i]*vm_fr[h]*vm_to[i]*sin(va_fr[h]-va_to[i]) for i in PMs.phase_ids(pm)) ) )
    @NLconstraint(pm.model, q_fr == z*(-b_fr[h]*vm_fr[h]^2 - sum(
                                            b[h,i]*vm_fr[h]*vm_fr[i]*cos(va_fr[h]-va_fr[i]) -
                                            g[h,i]*vm_fr[h]*vm_fr[i]*sin(va_fr[h]-va_fr[i]) -
                                            b[h,i]*vm_fr[h]*vm_to[i]*cos(va_fr[h]-va_to[i]) +
                                            g[h,i]*vm_fr[h]*vm_to[i]*sin(va_fr[h]-va_to[i]) for i in PMs.phase_ids(pm)) ) )
end

"""
```
p[t_idx] == z*(g*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
q[t_idx] == z*(-(b+c/2)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
```
"""
function constraint_ohms_tp_yt_to_on_off(pm::GenericPowerModel{T}, n::Int, h::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractACPForm
    p_to  = var(pm, n, h,  :p, t_idx)
    q_to  = var(pm, n, h,  :q, t_idx)
    vm_fr = [var(pm, n, i, :vm, f_bus) for i in PMs.phase_ids(pm)]
    vm_to = [var(pm, n, i, :vm, t_bus) for i in PMs.phase_ids(pm)]
    va_fr = [var(pm, n, j, :va, f_bus) for j in PMs.phase_ids(pm)]
    va_to = [var(pm, n, j, :va, t_bus) for j in PMs.phase_ids(pm)]
    z = var(pm, n, h, :branch_z, i)

    @NLconstraint(pm.model, p_to == z*( g_to[h]*vm_to[h]^2 + sum(
                                        g[h,i]*vm_to[h]*vm_to[i]*cos(va_to[h]-va_to[i]) +
                                        b[h,i]*vm_to[h]*vm_to[i]*sin(va_to[h]-va_to[i]) -
                                        g[h,i]*vm_to[h]*vm_fr[i]*cos(va_to[h]-va_fr[i]) -
                                        b[h,i]*vm_to[h]*vm_fr[i]*sin(va_to[h]-va_fr[i]) for i in PMs.phase_ids(pm)) ) )
    @NLconstraint(pm.model, q_to == z*(-b_to[h]*vm_to[h]^2 - sum(
                                        b[h,i]*vm_to[h]*vm_to[i]*cos(va_to[h]-va_to[i]) -
                                        g[h,i]*vm_to[h]*vm_to[i]*sin(va_to[h]-va_to[i]) -
                                        b[h,i]*vm_to[h]*vm_fr[i]*cos(va_to[h]-va_fr[i]) +
                                        g[h,i]*vm_to[h]*vm_fr[i]*sin(va_to[h]-va_fr[i]) for i in PMs.phase_ids(pm)) ) )
end
