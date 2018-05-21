# Three-phase specific constraints


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] == g/tm*w_fr_ne[i] + (-g*tr+b*ti)/tm*(wr_ne[i]) + (-b*tr-g*ti)/tm*(wi_ne[i])
q[f_idx] == -(b+c/2)/tm*w_fr_ne[i] - (-b*tr-g*ti)/tm*(wr_ne[i]) + (-g*tr+b*ti)/tm*(wi_ne[i])
```
"""
function constraint_ohms_tp_yt_from_ne(pm::GenericPowerModel{T}, n::Int, h::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractWRForm
    p_fr = var(pm, n, h,    :p_ne, f_idx)
    q_fr = var(pm, n, h,    :q_ne, f_idx)
    w_fr = var(pm, n, h, :w_fr_ne, i)
    wr   = var(pm, n, h,   :wr_ne, i)
    wi   = var(pm, n, h,   :wi_ne, i)

    @constraint(pm.model, p_fr ==  (g[h,h]+g_fr)/tm^2*w_fr + (-g*tr+b*ti)[h]/tm^2*wr + (-b*tr-g*ti)[h]/tm^2*wi )
    @constraint(pm.model, q_fr == -(b[h,h]+b_fr)/tm^2*w_fr - (-b*tr-g*ti)[h]/tm^2*wr + (-g*tr+b*ti)[h]/tm^2*wi )
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] == g*w_to_ne[i] + (-g*tr-b*ti)/tm*(wr_ne[i]) + (-b*tr+g*ti)/tm*(-wi_ne[i])
q[t_idx] == -(b+c/2)*w_to_ne[i] - (-b*tr+g*ti)/tm*(wr_ne[i]) + (-g*tr-b*ti)/tm*(-wi_ne[i])
```
"""
function constraint_ohms_tp_yt_to_ne(pm::GenericPowerModel{T}, n::Int, h::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractWRForm
    p_to = var(pm, n, h,    :p_ne, t_idx)
    q_to = var(pm, n, h,    :q_ne, t_idx)
    w_to = var(pm, n, h, :w_to_ne, i)
    wr   = var(pm, n, h,   :wr_ne, i)
    wi   = var(pm, n, h,   :wi_ne, i)

    @constraint(pm.model, p_to ==  (g[h,h]+g_to)*w_to + (-g*tr-b*ti)[h]/tm^2*wr + (-b*tr+g*ti)[h]/tm^2*-wi )
    @constraint(pm.model, q_to == -(b[h,h]+b_to)*w_to - (-b*tr+g*ti)[h]/tm^2*wr + (-g*tr-b*ti)[h]/tm^2*-wi )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] ==        g/tm*w_fr[i] + (-g*tr+b*ti)/tm*(wr[i]) + (-b*tr-g*ti)/tm*(wi[i])
q[f_idx] == -(b+c/2)/tm*w_fr[i] - (-b*tr-g*ti)/tm*(wr[i]) + (-g*tr+b*ti)/tm*(wi[i])
```
"""
function constraint_ohms_tp_yt_from_on_off(pm::GenericPowerModel{T}, n::Int, h::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractWRForm
    p_fr = var(pm, n, h, :p, f_idx)
    q_fr = var(pm, n, h, :q, f_idx)
    w_fr = var(pm, n, h, :w_fr, i)
    wr   = var(pm, n, h, :wr, i)
    wi   = var(pm, n, h, :wi, i)

    @constraint(pm.model, p_fr ==  (g[h,h]+g_fr)/tm^2*w_fr + (-g*tr+b*ti)[h]/tm^2*wr + (-b*tr-g*ti)[h]/tm^2*wi )
    @constraint(pm.model, q_fr == -(b[h,h]+b_fr)/tm^2*w_fr - (-b*tr-g*ti)[h]/tm^2*wr + (-g*tr+b*ti)[h]/tm^2*wi )
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==        g*w_to[i] + (-g*tr-b*ti)/tm*(wr[i]) + (-b*tr+g*ti)/tm*(-wi[i])
q[t_idx] == -(b+c/2)*w_to[i] - (-b*tr+g*ti)/tm*(wr[i]) + (-g*tr-b*ti)/tm*(-wi[i])
```
"""
function constraint_ohms_tp_yt_to_on_off(pm::GenericPowerModel{T}, n::Int, h::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractWRForm
    p_to = var(pm, n, h, :p, t_idx)
    q_to = var(pm, n, h, :q, t_idx)
    w_to = var(pm, n, h, :w_to, i)
    wr = var(pm, n, h, :wr, i)
    wi = var(pm, n, h, :wi, i)

    @constraint(pm.model, p_to ==  (g[h,h]+g_to)*w_to + (-g*tr-b*ti)[h]/tm^2*wr + (-b*tr+g*ti)[h]/tm^2*-wi )
    @constraint(pm.model, q_to == -(b[h,h]+b_to)*w_to - (-b*tr+g*ti)[h]/tm^2*wr + (-g*tr-b*ti)[h]/tm^2*-wi )
end