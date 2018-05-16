# Three-phase specific constraints


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] ==  (g+g_fr)/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
q[f_idx] == -(b+b_fr)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
```
"""
function constraint_ohms_yt_from(pm::GenericPowerModel{T}, n::Int, h::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractACPForm
    p_fr  = var(pm, n, h,  :p, f_idx)
    q_fr  = var(pm, n, h,  :q, f_idx)
    vm_fr = var(pm, n, h, :vm, f_bus)
    vm_to = var(pm, n, h, :vm, t_bus)
    va_fr = var(pm, n, h, :va, f_bus)
    va_to = var(pm, n, h, :va, t_bus)

    @NLconstraint(pm.model, p_fr ==  (g[h,h]+g_fr)/tm^2*vm_fr^2 + (-g*tr+b*ti)[h]/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)[h]/tm^2*(vm_fr*vm_to*sin(va_fr-va_to)) )
    @NLconstraint(pm.model, q_fr == -(b[h,h]+b_fr)/tm^2*vm_fr^2 - (-b*tr-g*ti)[h]/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)[h]/tm^2*(vm_fr*vm_to*sin(va_fr-va_to)) )
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_ohms_yt_to(pm::GenericPowerModel{T}, n::Int, h::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractACPForm
    p_to  = var(pm, n, h,  :p, t_idx)
    q_to  = var(pm, n, h,  :q, t_idx)
    vm_fr = var(pm, n, h, :vm, f_bus)
    vm_to = var(pm, n, h, :vm, t_bus)
    va_fr = var(pm, n, h, :va, f_bus)
    va_to = var(pm, n, h, :va, t_bus)

    @NLconstraint(pm.model, p_to ==  (g[h,h]+g_to)*vm_to^2 + (-g*tr-b*ti)[h]/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)[h]/tm^2*(vm_to*vm_fr*sin(va_to-va_fr)) )
    @NLconstraint(pm.model, q_to == -(b[h,h]+b_to)*vm_to^2 - (-b*tr+g*ti)[h]/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)[h]/tm^2*(vm_to*vm_fr*sin(va_to-va_fr)) )
end

"""
Creates Ohms constraints for AC models (y post fix indicates that Y values are in rectangular form)

```
p[f_idx] ==  (g+g_fr)*(v[f_bus]/tr)^2 + -g*v[f_bus]/tr*v[t_bus]*cos(t[f_bus]-t[t_bus]-as) + -b*v[f_bus]/tr*v[t_bus]*sin(t[f_bus]-t[t_bus]-as)
q[f_idx] == -(b+b_fr)*(v[f_bus]/tr)^2 + b*v[f_bus]/tr*v[t_bus]*cos(t[f_bus]-t[t_bus]-as) + -g*v[f_bus]/tr*v[t_bus]*sin(t[f_bus]-t[t_bus]-as)
```
"""
function constraint_ohms_y_from(pm::GenericPowerModel{T}, n::Int, h::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tm, ta) where T <: PMs.AbstractACPForm
    p_fr  = var(pm, n, h,  :p, f_idx)
    q_fr  = var(pm, n, h,  :q, f_idx)
    vm_fr = var(pm, n, h, :vm, f_bus)
    vm_to = var(pm, n, h, :vm, t_bus)
    va_fr = var(pm, n, h, :va, f_bus)
    va_to = var(pm, n, h, :va, t_bus)

    @NLconstraint(pm.model, p_fr ==  (g[h,h]+g_fr)*(vm_fr/tm)^2 - g[h,h]*vm_fr/tm*vm_to*cos(va_fr-va_to-ta) + -b[h,h]*vm_fr/tm*vm_to*sin(va_fr-va_to-ta) )
    @NLconstraint(pm.model, q_fr == -(b[h,h]+b_fr)*(vm_fr/tm)^2 + b[h,h]*vm_fr/tm*vm_to*cos(va_fr-va_to-ta) + -g[h,h]*vm_fr/tm*vm_to*sin(va_fr-va_to-ta) )
end

"""
Creates Ohms constraints for AC models (y post fix indicates that Y values are in rectangular form)

```
p[t_idx] == (g+g_to)*v[t_bus]^2 + -g*v[t_bus]*v[f_bus]/tr*cos(t[t_bus]-t[f_bus]+as) + -b*v[t_bus]*v[f_bus]/tr*sin(t[t_bus]-t[f_bus]+as)
q_to == -(b+b_to)*v[t_bus]^2 + b*v[t_bus]*v[f_bus]/tr*cos(t[f_bus]-t[t_bus]+as) + -g*v[t_bus]*v[f_bus]/tr*sin(t[t_bus]-t[f_bus]+as)
```
"""
function constraint_ohms_y_to(pm::GenericPowerModel{T}, n::Int, h::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tm, ta) where T <: PMs.AbstractACPForm
    p_to  = var(pm, n, h,  :p, t_idx)
    q_to  = var(pm, n, h,  :q, t_idx)
    vm_fr = var(pm, n, h, :vm, f_bus)
    vm_to = var(pm, n, h, :vm, t_bus)
    va_fr = var(pm, n, h, :va, f_bus)
    va_to = var(pm, n, h, :va, t_bus)

    @NLconstraint(pm.model, p_to ==  (g[h,h]+g_to)*vm_to^2 - g[h,h]*vm_to*vm_fr/tm*cos(va_to-va_fr+ta) + -b[h,h]*vm_to*vm_fr/tm*sin(va_to-va_fr+ta) )
    @NLconstraint(pm.model, q_to == -(b[h,h]+b_to)*vm_to^2 + b[h,h]*vm_to*vm_fr/tm*cos(va_to-va_fr+ta) + -g[h,h]*vm_to*vm_fr/tm*sin(va_to-va_fr+ta) )
end


"""
```
p[f_idx] == z*(g/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
q[f_idx] == z*(-(b+c/2)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
```
"""
function constraint_ohms_yt_from_on_off(pm::GenericPowerModel{T}, n::Int, h::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractACPForm
    p_fr  = var(pm, n, h,  :p, f_idx)
    q_fr  = var(pm, n, h,  :q, f_idx)
    vm_fr = var(pm, n, h, :vm, f_bus)
    vm_to = var(pm, n, h, :vm, t_bus)
    va_fr = var(pm, n, h, :va, f_bus)
    va_to = var(pm, n, h, :va, t_bus)
    z = var(pm, n, h, :branch_z, i)

    @NLconstraint(pm.model, p_fr == z*( (g[h,h]+g_fr)/tm^2*vm_fr^2 + (-g*tr+b*ti)[h]/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)[h]/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))) )
    @NLconstraint(pm.model, q_fr == z*(-(b[h,h]+b_fr)/tm^2*vm_fr^2 - (-b*tr-g*ti)[h]/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)[h]/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))) )
end

"""
```
p[t_idx] == z*(g*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
q[t_idx] == z*(-(b+c/2)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
```
"""
function constraint_ohms_yt_to_on_off(pm::GenericPowerModel{T}, n::Int, h::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractACPForm
    p_to  = var(pm, n, h,  :p, t_idx)
    q_to  = var(pm, n, h,  :q, t_idx)
    vm_fr = var(pm, n, h, :vm, f_bus)
    vm_to = var(pm, n, h, :vm, t_bus)
    va_fr = var(pm, n, h, :va, f_bus)
    va_to = var(pm, n, h, :va, t_bus)
    z = var(pm, n, h, :branch_z, i)

    @NLconstraint(pm.model, p_to == z*( (g[h,h]+g_to)*vm_to^2 + (-g*tr-b*ti)[h]/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)[h]/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))) )
    @NLconstraint(pm.model, q_to == z*(-(b[h,h]+b_to)*vm_to^2 - (-b*tr+g*ti)[h]/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)[h]/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))) )
end

"""
```
p_ne[f_idx] == z*(g/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
q_ne[f_idx] == z*(-(b+c/2)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
```
"""
function constraint_ohms_yt_from_ne(pm::GenericPowerModel{T}, n::Int, h::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractACPForm
    p_fr  = var(pm, n, h, :p_ne, f_idx)
    q_fr  = var(pm, n, h, :q_ne, f_idx)
    vm_fr = var(pm, n, h,   :vm, f_bus)
    vm_to = var(pm, n, h,   :vm, t_bus)
    va_fr = var(pm, n, h,   :va, f_bus)
    va_to = var(pm, n, h,   :va, t_bus)
    z = var(pm, n, h, :branch_ne, i)

    @NLconstraint(pm.model, p_fr == z*( (g[h,h]+g_fr)/tm^2*vm_fr^2 + (-g*tr+b*ti)[h]/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)[h]/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))) )
    @NLconstraint(pm.model, q_fr == z*(-(b[h,h]+b_fr)/tm^2*vm_fr^2 - (-b*tr-g*ti)[h]/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)[h]/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))) )
end

"""
```
p_ne[t_idx] == z*(g*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
q_ne[t_idx] == z*(-(b+c/2)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
```
"""
function constraint_ohms_yt_to_ne(pm::GenericPowerModel{T}, n::Int, h::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractACPForm
    p_to = var(pm, n, h, :p_ne, t_idx)
    q_to = var(pm, n, h, :q_ne, t_idx)
    vm_fr = var(pm, n, h, :vm, f_bus)
    vm_to = var(pm, n, h, :vm, t_bus)
    va_fr = var(pm, n, h, :va, f_bus)
    va_to = var(pm, n, h, :va, t_bus)
    z = var(pm, n, h, :branch_ne, i)

    @NLconstraint(pm.model, p_to == z*( (g[h,h]+g_to)*vm_to^2 + (-g*tr-b*ti)[h]/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)[h]/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))) )
    @NLconstraint(pm.model, q_to == z*(-(b[h,h]+b_to)*vm_to^2 - (-b*tr+g*ti)[h]/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)[h]/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))) )
end
