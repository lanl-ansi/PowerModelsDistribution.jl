# Three-phase specific constraints


"do nothing, this model does not have complex voltage constraints"
function constraint_tp_voltage(pm::PMs.GenericPowerModel{T}, n::Int, c::Int) where T <: PMs.AbstractACPForm
end


""
function constraint_kcl_shunt_slack(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: PMs.AbstractACPForm
    vm = PMs.var(pm, n, c, :vm, i)
    p_slack = PMs.var(pm, n, c, :p_slack, i)
    q_slack = PMs.var(pm, n, c, :q_slack, i)
    p = PMs.var(pm, n, c, :p)
    q = PMs.var(pm, n, c, :q)
    pg = PMs.var(pm, n, c, :pg)
    qg = PMs.var(pm, n, c, :qg)
    p_dc = PMs.var(pm, n, c, :p_dc)
    q_dc = PMs.var(pm, n, c, :q_dc)

    PMs.con(pm, n, c, :kcl_p)[i] = JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*vm^2 + p_slack)
    PMs.con(pm, n, c, :kcl_q)[i] = JuMP.@constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*vm^2 + q_slack)
end


""
function constraint_kcl_shunt_trans(pm::PMs.GenericPowerModel{T}, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: PMs.AbstractACPForm
    vm = PMs.var(pm, nw, c, :vm, i)
    p = PMs.var(pm, nw, c, :p)
    q = PMs.var(pm, nw, c, :q)
    pg = PMs.var(pm, nw, c, :pg)
    qg = PMs.var(pm, nw, c, :qg)
    p_dc = PMs.var(pm, nw, c, :p_dc)
    q_dc = PMs.var(pm, nw, c, :q_dc)
    p_trans = PMs.var(pm, nw, c, :pt)
    q_trans = PMs.var(pm,  nw, c, :qt)

    PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) + sum(p_trans[a_trans] for a_trans in bus_arcs_trans) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*vm^2)
    PMs.con(pm, nw, c, :kcl_q)[i] = JuMP.@constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) + sum(q_trans[a_trans] for a_trans in bus_arcs_trans) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*vm^2)
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] ==  (g+g_fr)/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
q[f_idx] == -(b+b_fr)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
```
"""
function constraint_ohms_tp_yt_from(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractACPForm
    p_fr  = PMs.var(pm, n, c,  :p, f_idx)
    q_fr  = PMs.var(pm, n, c,  :q, f_idx)
    vm_fr = [PMs.var(pm, n, d, :vm, f_bus) for d in PMs.conductor_ids(pm)]
    vm_to = [PMs.var(pm, n, d, :vm, t_bus) for d in PMs.conductor_ids(pm)]
    va_fr = [PMs.var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [PMs.var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]

    JuMP.@NLconstraint(pm.model, p_fr ==  (g_fr[c]+g[c,c]) * vm_fr[c]^2 +
                                    sum( g[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                                         b[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in PMs.conductor_ids(pm) if d != c) +
                                    sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                        -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in PMs.conductor_ids(pm)) )
    JuMP.@NLconstraint(pm.model, q_fr == -(b_fr[c]+b[c,c]) *vm_fr[c]^2 -
                                    sum( b[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                                         g[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in PMs.conductor_ids(pm) if d != c) -
                                    sum(-b[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                         g[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in PMs.conductor_ids(pm)) )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_ohms_tp_yt_to(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractACPForm
    p_to  = PMs.var(pm, n, c,  :p, t_idx)
    q_to  = PMs.var(pm, n, c,  :q, t_idx)
    vm_fr = [PMs.var(pm, n, d, :vm, f_bus) for d in PMs.conductor_ids(pm)]
    vm_to = [PMs.var(pm, n, d, :vm, t_bus) for d in PMs.conductor_ids(pm)]
    va_fr = [PMs.var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [PMs.var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]

    JuMP.@NLconstraint(pm.model, p_to == (g_to[c]+g[c,c])*vm_to[c]^2 +
                                   sum( g[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) +
                                        b[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in PMs.conductor_ids(pm) if d != c) +
                                   sum(-g[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                       -b[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in PMs.conductor_ids(pm)) )
    JuMP.@NLconstraint(pm.model, q_to == -(b_to[c]+b[c,c])*vm_to[c]^2 -
                                   sum( b[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) -
                                        g[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in PMs.conductor_ids(pm) if d != c) -
                                   sum(-b[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                        g[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in PMs.conductor_ids(pm)) )
end


"""
```
p[f_idx] == z*(g/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
q[f_idx] == z*(-(b+c/2)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
```
"""
function constraint_ohms_tp_yt_from_on_off(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractACPForm
    p_fr  = PMs.var(pm, n, c,  :p, f_idx)
    q_fr  = PMs.var(pm, n, c,  :q, f_idx)
    vm_fr = [PMs.var(pm, n, d, :vm, f_bus) for d in PMs.conductor_ids(pm)]
    vm_to = [PMs.var(pm, n, d, :vm, t_bus) for d in PMs.conductor_ids(pm)]
    va_fr = [PMs.var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [PMs.var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]
    z = PMs.var(pm, n, c, :branch_z, i)

    JuMP.@NLconstraint(pm.model, p_fr == z*((g_fr[c]+g[c,c]) * vm_fr[c]^2 +
                                        sum( g[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                                             b[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in PMs.conductor_ids(pm) if d != c) +
                                        sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                            -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in PMs.conductor_ids(pm))) )
    JuMP.@NLconstraint(pm.model, q_fr == z*(-(b_fr[c]+b[c,c]) *vm_fr[c]^2 -
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
function constraint_ohms_tp_yt_to_on_off(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractACPForm
    p_to  = PMs.var(pm, n, c,  :p, t_idx)
    q_to  = PMs.var(pm, n, c,  :q, t_idx)
    vm_fr = [PMs.var(pm, n, d, :vm, f_bus) for d in PMs.conductor_ids(pm)]
    vm_to = [PMs.var(pm, n, d, :vm, t_bus) for d in PMs.conductor_ids(pm)]
    va_fr = [PMs.var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [PMs.var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]
    z = PMs.var(pm, n, c, :branch_z, i)

    JuMP.@NLconstraint(pm.model, p_to == z*((g_to[c]+g[c,c])*vm_to[c]^2 +
                                        sum( g[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) +
                                             b[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in PMs.conductor_ids(pm) if d != c) +
                                        sum(-g[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                            -b[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in PMs.conductor_ids(pm))) )
    JuMP.@NLconstraint(pm.model, q_to == z*(-(b_to[c]+b[c,c])*vm_to[c]^2 -
                                        sum( b[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) -
                                             g[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in PMs.conductor_ids(pm) if d != c) -
                                        sum(-b[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                             g[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in PMs.conductor_ids(pm))) )
end


"Links the voltage at both windings of a fixed tap transformer."
function constraint_tp_trans_voltage(pm::PMs.GenericPowerModel{T}, nw::Int, i::Int, f_bus::Int, t_bus::Int, tm::PMs.MultiConductorVector, Tv_fr, Tv_im, Cv_to) where T <: PMs.AbstractACPForm
    ncnd  = 3
    # from side
    vm_fr = [PMs.var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [PMs.var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    # to side
    vm_to = [PMs.var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
    va_to = [PMs.var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
    # the intermediate bus voltage is saved as an expression
    # vm_im[c] = vm_to[c]*tm[c]*Cv_to
    # va_im = va_to
    for n in 1:size(Tv_fr)[1]
        JuMP.@NLconstraint(pm.model,
              sum(Tv_fr[n,c]*vm_fr[c]*cos(va_fr[c]) for c in 1:ncnd)
            ==sum(Tv_im[n,c]*(vm_to[c]*tm[c]*Cv_to)*cos(va_to[c]) for c in 1:ncnd)
        )
        JuMP.@NLconstraint(pm.model,
              sum(Tv_fr[n,c]*vm_fr[c]*sin(va_fr[c]) for c in 1:ncnd)
            ==sum(Tv_im[n,c]*(vm_to[c]*tm[c]*Cv_to)*sin(va_to[c]) for c in 1:ncnd)
        )
    end
end


"Links the power flowing into both windings of a fixed tap transformer."
function constraint_tp_trans_flow(pm::PMs.GenericPowerModel{T}, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, tm::PMs.MultiConductorVector, Ti_fr, Ti_im, Cv_to) where T <: PMs.AbstractACPForm
    ncnd  = 3
    # from side variables
    vm_fr = [PMs.var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [PMs.var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    p_fr = [PMs.var(pm, nw, c, :pt, f_idx) for c in 1:ncnd]
    q_fr = [PMs.var(pm, nw, c, :qt, f_idx) for c in 1:ncnd]
    # to side
    vm_to = [PMs.var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
    va_to = [PMs.var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
    p_to = [PMs.var(pm, nw, c, :pt, t_idx) for c in 1:ncnd]
    q_to = [PMs.var(pm, nw, c, :qt, t_idx) for c in 1:ncnd]
    # the intermediate bus voltage is saved as an expression
    # vm_im[c] = vm_to[c]*tm[c]*Cv_to
    # va_im = va_to
    for n in 1:size(Ti_fr)[1]
        # i_fr_re[c] = 1/vm_fr[c]*(p_fr[c]*cos(va_fr[c])+q_fr[c]*sin(va_fr[c]))
        # i_fr_im[c] = 1/vm_fr[c]*(p_fr[c]*sin(va_fr[c])-q_fr[c]*cos(va_fr[c]))
        # i_to_re[c] = 1/vm_to[c]*(p_to[c]*cos(va_im[c])+q_to[c]*sin(va_im[c]))
        # i_to_im[c] = 1/vm_to[c]*(p_to[c]*sin(va_im[c])-q_to[c]*cos(va_im[c]))
        JuMP.@NLconstraint(pm.model,
              sum(Ti_fr[n,c]*
                    1/vm_fr[c]*(p_fr[c]*cos(va_fr[c])+q_fr[c]*sin(va_fr[c])) # i_fr_re[c]
              for c in 1:ncnd)
            + sum(Ti_im[n,c]*
                    1/(vm_to[c]*tm[c]*Cv_to)*(p_to[c]*cos(va_to[c])+q_to[c]*sin(va_to[c])) # i_to_re[c]
              for c in 1:ncnd)
            == 0
        )
        JuMP.@NLconstraint(pm.model,
              sum(Ti_fr[n,c]*
                    1/vm_fr[c]*(p_fr[c]*sin(va_fr[c])-q_fr[c]*cos(va_fr[c])) # i_fr_im[c]
              for c in 1:ncnd)
            + sum(Ti_im[n,c]*
                    1/(vm_to[c]*tm[c]*Cv_to)*(p_to[c]*sin(va_to[c])-q_to[c]*cos(va_to[c])) # i_to_im[c]
              for c in 1:ncnd)
            == 0
        )
    end
end


"Links the voltage at both windings of a variable tap transformer."
function constraint_tp_oltc_voltage(pm::PMs.GenericPowerModel{T}, nw::Int, i::Int, f_bus::Int, t_bus::Int, Tv_fr, Tv_im, Cv_to) where T <: PMs.AbstractACPForm
    ncnd  = 3
    # from side
    vm_fr = [PMs.var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [PMs.var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    # tap
    tap = [PMs.var(pm, nw, c, :tap)[i] for c in 1:ncnd]
    # to side
    vm_to = [PMs.var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
    va_to = [PMs.var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
    # the intermediate bus voltage is saved as an expression
    # vm_im[c] = vm_to[c]*tap[c]*Cv_to
    # va_im = va_to
    for n in 1:size(Tv_fr)[1]
        JuMP.@NLconstraint(pm.model,
              sum(Tv_fr[n,c]*vm_fr[c]*cos(va_fr[c]) for c in 1:ncnd)
            ==sum(Tv_im[n,c]*(vm_to[c]*tap[c]*Cv_to)*cos(va_to[c]) for c in 1:ncnd)
        )
        JuMP.@NLconstraint(pm.model,
              sum(Tv_fr[n,c]*vm_fr[c]*sin(va_fr[c]) for c in 1:ncnd)
            ==sum(Tv_im[n,c]*(vm_to[c]*tap[c]*Cv_to)*sin(va_to[c]) for c in 1:ncnd)
        )
    end
end


"Links the power flowing into both windings of a variable tap transformer."
function constraint_tp_oltc_flow(pm::PMs.GenericPowerModel{T}, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, Ti_fr, Ti_im, Cv_to) where T <: PMs.AbstractACPForm
    ncnd  = 3
    # from side variables
    vm_fr = [PMs.var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [PMs.var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    p_fr = [PMs.var(pm, nw, c, :pt, f_idx) for c in 1:ncnd]
    q_fr = [PMs.var(pm, nw, c, :qt, f_idx) for c in 1:ncnd]
    # tap
    tap = [PMs.var(pm, nw, c, :tap, i) for c in 1:ncnd]
    # to side
    vm_to = [PMs.var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
    va_to = [PMs.var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
    p_to = [PMs.var(pm, nw, c, :pt, t_idx) for c in 1:ncnd]
    q_to = [PMs.var(pm, nw, c, :qt, t_idx) for c in 1:ncnd]
    # the intermediate bus voltage is saved as an expression
    # vm_im[c] = vm_to[c]*tap[c]*Cv_to
    # va_im = va_to
    for n in 1:size(Ti_fr)[1]
        # i_fr_re[c] = 1/vm_fr[c]*(p_fr[c]*cos(va_fr[c])+q_fr[c]*sin(va_fr[c]))
        # i_fr_im[c] = 1/vm_fr[c]*(p_fr[c]*sin(va_fr[c])-q_fr[c]*cos(va_fr[c]))
        # i_to_re[c] = 1/vm_to[c]*(p_to[c]*cos(va_im[c])+q_to[c]*sin(va_im[c]))
        # i_to_im[c] = 1/vm_to[c]*(p_to[c]*sin(va_im[c])-q_to[c]*cos(va_im[c]))
        JuMP.@NLconstraint(pm.model,
              sum(Ti_fr[n,c]*
                    1/vm_fr[c]*(p_fr[c]*cos(va_fr[c])+q_fr[c]*sin(va_fr[c])) # i_fr_re[c]
              for c in 1:ncnd)
            + sum(Ti_im[n,c]*
                    1/(vm_to[c]*tap[c]*Cv_to)*(p_to[c]*cos(va_to[c])+q_to[c]*sin(va_to[c])) # i_to_re[c]
              for c in 1:ncnd)
            == 0
        )
        JuMP.@NLconstraint(pm.model,
              sum(Ti_fr[n,c]*
                    1/vm_fr[c]*(p_fr[c]*sin(va_fr[c])-q_fr[c]*cos(va_fr[c])) # i_fr_im[c]
              for c in 1:ncnd)
            + sum(Ti_im[n,c]*
                    1/(vm_to[c]*tap[c]*Cv_to)*(p_to[c]*sin(va_to[c])-q_to[c]*cos(va_to[c])) # i_to_im[c]
              for c in 1:ncnd)
            == 0
        )
    end
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_tp_vm_vuf(pm::PMs.GenericPowerModel{T}, nw::Int, bus_id::Int, vufmax::Float64) where T <: PMs.AbstractACPForm
    if !haskey(PMs.var(pm, pm.cnw), :vmpossqr)
        PMs.var(pm, pm.cnw)[:vmpossqr] = Dict{Int, Any}()
        PMs.var(pm, pm.cnw)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [PMs.var(pm, nw, i, :vm, bus_id) for i in 1:3]
    (va_a, va_b, va_c) = [PMs.var(pm, nw, i, :va, bus_id) for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)
    # real and imaginary components of U+
    vrepos = JuMP.@NLexpression(pm.model,
        (vm_a*cos(va_a) + are*vm_b*cos(va_b) - aim*vm_b*sin(va_b) + a2re*vm_c*cos(va_c) - a2im*vm_c*sin(va_c))/3
    )
    vimpos = JuMP.@NLexpression(pm.model,
        (vm_a*sin(va_a) + are*vm_b*sin(va_b) + aim*vm_b*cos(va_b) + a2re*vm_c*sin(va_c) + a2im*vm_c*cos(va_c))/3
    )
    # square of magnitude of U+, |U+|^2
    vmpossqr = JuMP.@NLexpression(pm.model, vrepos^2+vimpos^2)
    # real and imaginary components of U-
    vreneg = JuMP.@NLexpression(pm.model,
        (vm_a*cos(va_a) + a2re*vm_b*cos(va_b) - a2im*vm_b*sin(va_b) + are*vm_c*cos(va_c) - aim*vm_c*sin(va_c))/3
    )
    vimneg = JuMP.@NLexpression(pm.model,
        (vm_a*sin(va_a) + a2re*vm_b*sin(va_b) + a2im*vm_b*cos(va_b) + are*vm_c*sin(va_c) + aim*vm_c*cos(va_c))/3
    )
    # square of magnitude of U-, |U-|^2
    vmnegsqr = JuMP.@NLexpression(pm.model, vreneg^2+vimneg^2)
    # finally, apply constraint
    JuMP.@NLconstraint(pm.model, vmnegsqr <= vufmax^2*vmpossqr)
    # DEBUGGING: save references for post check
    #var(pm, pm.cnw, :vmpossqr)[bus_id] = vmpossqr
    #var(pm, pm.cnw, :vmnegsqr)[bus_id] = vmnegsqr
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_tp_vm_neg_seq(pm::PMs.GenericPowerModel{T}, nw::Int, bus_id::Int, vmnegmax::Float64) where T <: PMs.AbstractACPForm
    if !haskey(PMs.var(pm, pm.cnw), :vmpossqr)
        PMs.var(pm, pm.cnw)[:vmpossqr] = Dict{Int, Any}()
        PMs.var(pm, pm.cnw)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [PMs.var(pm, nw, i, :vm, bus_id) for i in 1:3]
    (va_a, va_b, va_c) = [PMs.var(pm, nw, i, :va, bus_id) for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)
    # real and imaginary components of U-
    vreneg = JuMP.@NLexpression(pm.model,
        (vm_a*cos(va_a) + a2re*vm_b*cos(va_b) - a2im*vm_b*sin(va_b) + are*vm_c*cos(va_c) - aim*vm_c*sin(va_c))/3
    )
    vimneg = JuMP.@NLexpression(pm.model,
        (vm_a*sin(va_a) + a2re*vm_b*sin(va_b) + a2im*vm_b*cos(va_b) + are*vm_c*sin(va_c) + aim*vm_c*cos(va_c))/3
    )
    # square of magnitude of U-, |U-|^2
    vmnegsqr = JuMP.@NLexpression(pm.model, vreneg^2+vimneg^2)
    # finally, apply constraint
    JuMP.@NLconstraint(pm.model, vmnegsqr <= vmnegmax^2)
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_tp_vm_pos_seq(pm::PMs.GenericPowerModel{T}, nw::Int, bus_id::Int, vmposmax::Float64) where T <: PMs.AbstractACPForm
    if !haskey(PMs.var(pm, pm.cnw), :vmpossqr)
        PMs.var(pm, pm.cnw)[:vmpossqr] = Dict{Int, Any}()
        PMs.var(pm, pm.cnw)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [PMs.var(pm, nw, i, :vm, bus_id) for i in 1:3]
    (va_a, va_b, va_c) = [PMs.var(pm, nw, i, :va, bus_id) for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)
    # real and imaginary components of U+
    vrepos = JuMP.@NLexpression(pm.model,
        (vm_a*cos(va_a) + are*vm_b*cos(va_b) - aim*vm_b*sin(va_b) + a2re*vm_c*cos(va_c) - a2im*vm_c*sin(va_c))/3
    )
    vimpos = JuMP.@NLexpression(pm.model,
        (vm_a*sin(va_a) + are*vm_b*sin(va_b) + aim*vm_b*cos(va_b) + a2re*vm_c*sin(va_c) + a2im*vm_c*cos(va_c))/3
    )
    # square of magnitude of U+, |U+|^2
    vmpossqr = JuMP.@NLexpression(pm.model, vrepos^2+vimpos^2)
    # finally, apply constraint
    JuMP.@NLconstraint(pm.model, vmpossqr <= vmposmax^2)
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_tp_vm_zero_seq(pm::PMs.GenericPowerModel{T}, nw::Int, bus_id::Int, vmzeromax::Float64) where T <: PMs.AbstractACPForm
    if !haskey(PMs.var(pm, pm.cnw), :vmpossqr)
        PMs.var(pm, pm.cnw)[:vmpossqr] = Dict{Int, Any}()
        PMs.var(pm, pm.cnw)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [PMs.var(pm, nw, i, :vm, bus_id) for i in 1:3]
    (va_a, va_b, va_c) = [PMs.var(pm, nw, i, :va, bus_id) for i in 1:3]
    # real and imaginary components of U+
    vrezero = JuMP.@NLexpression(pm.model,
        (vm_a*cos(va_a) + vm_b*cos(va_b) + vm_c*cos(va_c))/3
    )
    vimzero = JuMP.@NLexpression(pm.model,
        (vm_a*sin(va_a) + vm_b*sin(va_b) + vm_c*sin(va_c))/3
    )
    # square of magnitude of U+, |U+|^2
    vmzerosqr = JuMP.@NLexpression(pm.model, vrezero^2+vimzero^2)
    # finally, apply constraint
    JuMP.@NLconstraint(pm.model, vmzerosqr <= vmzeromax^2)
end


"""
"""
function constraint_tp_vm_ll(pm::PMs.GenericPowerModel{T}, nw::Int, bus_id::Int, vm_ll_min::PMs.MultiConductorVector, vm_ll_max::PMs.MultiConductorVector) where T <: PMs.AbstractACPForm
    # 3 conductors asserted in template already
    vm_ln = [PMs.var(pm, nw, i, :vm, bus_id) for i in 1:3]
    va_ln = [PMs.var(pm, nw, i, :va, bus_id) for i in 1:3]
    vr_ll = JuMP.@NLexpression(pm.model, [i in 1:3],
        vm_ln[i]*cos(va_ln[i]) - vm_ln[i%3+1]*cos(va_ln[i%3+1])
    )
    vi_ll = JuMP.@NLexpression(pm.model, [i in 1:3],
        vm_ln[i]*sin(va_ln[i]) - vm_ln[i%3+1]*sin(va_ln[i%3+1])
    )
    for c in 1:3
        # factor of 3 is needed because vm_ll bounds are with respect to the
        # LL base, not the LN base
        if vm_ll_min[c] > 0
            JuMP.@NLconstraint(pm.model, vr_ll[c]^2+vi_ll[c]^2 >= vm_ll_min[c]^2*3)
        end
        if vm_ll_max[c] < Inf
            JuMP.@NLconstraint(pm.model, vr_ll[c]^2+vi_ll[c]^2 <= vm_ll_max[c]^2*3)
        end
    end
end
