# Three-phase specific constraints


""
function variable_tp_voltage(pm::GenericPowerModel{T}; kwargs...) where T <: PMs.AbstractACRForm
    for c in PMs.conductor_ids(pm)
        PMs.variable_voltage(pm, cnd=c; kwargs...)
    end
    # local infeasbility issues without proper initialization;
    # convergence issues start when the equivalent angles of the starting point
    # are further away than 90 degrees from the solution (as given by ACP)
    # this is the default behaviour of PMs, initialize all phases as (1,0)
    # the magnitude seems to have a little effect on the convergence (>0.05)
    # updating the starting point to a balanced phasor does the job
    ncnd = length(PMs.conductor_ids(pm))
    theta = [wraptopi(2 * pi / ncnd * (1-c)) for c in 1:ncnd]
    vm = 1
    for c in 1:ncnd
        vr = vm*cos(theta[c])
        vi = vm*sin(theta[c])
        for id in PMs.ids(pm, :bus)
            setvalue(var(pm, pm.cnw, c, :vr, id), vr)
            setvalue(var(pm, pm.cnw, c, :vi, id), vi)
        end
    end
end


"delegate back to PowerModels"
function constraint_tp_voltage(pm::GenericPowerModel{T}, n::Int, c::Int) where T <: PMs.AbstractACRForm
    PMs.constraint_voltage(pm, n, c)
end


"Creates phase angle constraints at reference buses"
function constraint_tp_theta_ref(pm::GenericPowerModel{T}, n::Int, c::Int, d) where T <: PMs.AbstractACRForm
    vr = var(pm, n, c, :vr, d)
    vi = var(pm, n, c, :vi, d)
    nconductors = length(PMs.conductor_ids(pm))
    theta = wraptopi(2 * pi / nconductors * (1-c))
    # deal with cases first where tan(theta)==Inf or tan(theta)==0
    if theta == pi/2
        @constraint(pm.model, vr == 0)
        @constraint(pm.model, vi >= 0)
    elseif theta == -pi/2
        @constraint(pm.model, vr == 0)
        @constraint(pm.model, vi <= 0)
    elseif theta == 0
        @constraint(pm.model, vr >= 0)
        @constraint(pm.model, vi == 0)
    elseif theta == pi
        @constraint(pm.model, vr >= 0)
        @constraint(pm.model, vi == 0)
    else
        tanth = tan(theta)
        @constraint(pm.model, vi == tanth*vr)
        # theta also implies a sign for vr, vi
        if 0<=theta && theta <= pi
            @constraint(pm.model, vi >= 0)
        else
            @constraint(pm.model, vi <= 0)
        end
    end
end


""
function constraint_kcl_shunt_slack(pm::GenericPowerModel{T}, n::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: PMs.AbstractACRForm
    vr = var(pm, n, c, :vr, i)
    vr = var(pm, n, c, :vi, i)
    p_slack = var(pm, n, c, :p_slack, i)
    q_slack = var(pm, n, c, :q_slack, i)
    p = var(pm, n, c, :p)
    q = var(pm, n, c, :q)
    pg = var(pm, n, c, :pg)
    qg = var(pm, n, c, :qg)
    p_dc = var(pm, n, c, :p_dc)
    q_dc = var(pm, n, c, :q_dc)

    con(pm, n, c, :kcl_p)[i] = @constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*(vr^2 + vi^2) + p_slack)
    con(pm, n, c, :kcl_q)[i] = @constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*(vr^2 + vi^2) + q_slack)
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] ==  (g+g_fr)/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
q[f_idx] == -(b+b_fr)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
```
"""
function constraint_ohms_tp_yt_from(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractACRForm
    p_fr  = var(pm, n, c,  :p, f_idx)
    q_fr  = var(pm, n, c,  :q, f_idx)
    vr_fr = [var(pm, n, d, :vr, f_bus) for d in PMs.conductor_ids(pm)]
    vr_to = [var(pm, n, d, :vr, t_bus) for d in PMs.conductor_ids(pm)]
    vi_fr = [var(pm, n, d, :vi, f_bus) for d in PMs.conductor_ids(pm)]
    vi_to = [var(pm, n, d, :vi, t_bus) for d in PMs.conductor_ids(pm)]

    @NLconstraint(pm.model, p_fr ==  g_fr[c]*(vr_fr[c]^2+vi_fr[c]^2)+ sum(
                                         vr_fr[c]*(g[c,d]*(vr_fr[d]-vr_to[d])-b[c,d]*(vi_fr[d]-vi_to[d]))
                                        -vi_fr[c]*(-b[c,d]*(vr_fr[d]-vr_to[d])-g[c,d]*(vi_fr[d]-vi_to[d]))
                                    for d in PMs.conductor_ids(pm))
    )
    @NLconstraint(pm.model, q_fr ==  -b_fr[c]*(vr_fr[c]^2+vi_fr[c]^2)+ sum(
                                            -vr_fr[c]*(b[c,d]*(vr_fr[d]-vr_to[d])+g[c,d]*(vi_fr[d]-vi_to[d]))
                                            +vi_fr[c]*(g[c,d]*(vr_fr[d]-vr_to[d])-b[c,d]*(vi_fr[d]-vi_to[d]))
                                        for d in PMs.conductor_ids(pm))
    )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_ohms_tp_yt_to(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractACRForm
    constraint_ohms_tp_yt_from(pm, n, c, t_bus, f_bus, t_idx, f_idx, g, b, g_to, b_to, tr, ti, tm)
end


"""
```
p[f_idx] == z*(g/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
q[f_idx] == z*(-(b+c/2)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
```
"""
function constraint_ohms_tp_yt_from_on_off(pm::GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractACRForm
    p_fr  = var(pm, n, c,  :p, f_idx)
    q_fr  = var(pm, n, c,  :q, f_idx)
    vm_fr = [var(pm, n, d, :vm, f_bus) for d in PMs.conductor_ids(pm)]
    vm_to = [var(pm, n, d, :vm, t_bus) for d in PMs.conductor_ids(pm)]
    va_fr = [var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]
    z = var(pm, n, c, :branch_z, i)

    @NLconstraint(pm.model, p_fr == z*((g_fr[c]+g[c,c]) * vm_fr[c]^2 +
                                        sum( g[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                                             b[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in PMs.conductor_ids(pm) if d != c) +
                                        sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                            -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in PMs.conductor_ids(pm))) )
    @NLconstraint(pm.model, q_fr == z*(-(b_fr[c]+b[c,c]) *vm_fr[c]^2 -
                                        sum( b[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                                             g[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in PMs.conductor_ids(pm) if d != c) -
                                        sum(-b[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                             g[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in PMs.conductor_ids(pm))) )
end

"""
```
p[t_idx] == z*(g*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
q[t_idx] == z*(-(b+c/2)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
```
"""
function constraint_ohms_tp_yt_to_on_off(pm::GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractACRForm
    p_to  = var(pm, n, c,  :p, t_idx)
    q_to  = var(pm, n, c,  :q, t_idx)
    vm_fr = [var(pm, n, d, :vm, f_bus) for d in PMs.conductor_ids(pm)]
    vm_to = [var(pm, n, d, :vm, t_bus) for d in PMs.conductor_ids(pm)]
    va_fr = [var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]
    z = var(pm, n, c, :branch_z, i)

    @NLconstraint(pm.model, p_to == z*((g_to[c]+g[c,c])*vm_to[c]^2 +
                                        sum( g[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) +
                                             b[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in PMs.conductor_ids(pm) if d != c) +
                                        sum(-g[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                            -b[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in PMs.conductor_ids(pm))) )
    @NLconstraint(pm.model, q_to == z*(-(b_to[c]+b[c,c])*vm_to[c]^2 -
                                        sum( b[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) -
                                             g[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in PMs.conductor_ids(pm) if d != c) -
                                        sum(-b[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                             g[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in PMs.conductor_ids(pm))) )
end
