# Three-phase specific constraints


"do nothing, this model does not have complex voltage constraints"
function constraint_tp_model_voltage(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int) where T <: _PMs.AbstractACPForm
end


""
function constraint_tp_power_balance_shunt_slack(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: _PMs.AbstractACPForm
    vm = _PMs.var(pm, n, c, :vm, i)
    p_slack = _PMs.var(pm, n, c, :p_slack, i)
    q_slack = _PMs.var(pm, n, c, :q_slack, i)
    p = _PMs.var(pm, n, c, :p)
    q = _PMs.var(pm, n, c, :q)
    pg = _PMs.var(pm, n, c, :pg)
    qg = _PMs.var(pm, n, c, :qg)
    p_dc = _PMs.var(pm, n, c, :p_dc)
    q_dc = _PMs.var(pm, n, c, :q_dc)

    _PMs.con(pm, n, c, :kcl_p)[i] = JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*vm^2 + p_slack)
    _PMs.con(pm, n, c, :kcl_q)[i] = JuMP.@constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*vm^2 + q_slack)
end


""
function constraint_tp_power_balance_shunt_trans(pm::_PMs.GenericPowerModel{T}, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: _PMs.AbstractACPForm
    vm = _PMs.var(pm, nw, c, :vm, i)
    p = _PMs.var(pm, nw, c, :p)
    q = _PMs.var(pm, nw, c, :q)
    pg = _PMs.var(pm, nw, c, :pg)
    qg = _PMs.var(pm, nw, c, :qg)
    p_dc = _PMs.var(pm, nw, c, :p_dc)
    q_dc = _PMs.var(pm, nw, c, :q_dc)
    p_trans = _PMs.var(pm, nw, c, :pt)
    q_trans = _PMs.var(pm,  nw, c, :qt)

    _PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) + sum(p_trans[a_trans] for a_trans in bus_arcs_trans) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*vm^2)
    _PMs.con(pm, nw, c, :kcl_q)[i] = JuMP.@constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) + sum(q_trans[a_trans] for a_trans in bus_arcs_trans) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*vm^2)
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] ==  (g+g_fr)/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
q[f_idx] == -(b+b_fr)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
```
"""
function constraint_tp_ohms_yt_from(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: _PMs.AbstractACPForm
    p_fr  = _PMs.var(pm, n, c,  :p, f_idx)
    q_fr  = _PMs.var(pm, n, c,  :q, f_idx)
    vm_fr = [_PMs.var(pm, n, d, :vm, f_bus) for d in _PMs.conductor_ids(pm)]
    vm_to = [_PMs.var(pm, n, d, :vm, t_bus) for d in _PMs.conductor_ids(pm)]
    va_fr = [_PMs.var(pm, n, d, :va, f_bus) for d in _PMs.conductor_ids(pm)]
    va_to = [_PMs.var(pm, n, d, :va, t_bus) for d in _PMs.conductor_ids(pm)]

    JuMP.@NLconstraint(pm.model, p_fr ==  (g_fr[c]+g[c,c]) * vm_fr[c]^2 +
                                    sum( g[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                                         b[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in _PMs.conductor_ids(pm) if d != c) +
                                    sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                        -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in _PMs.conductor_ids(pm)) )
    JuMP.@NLconstraint(pm.model, q_fr == -(b_fr[c]+b[c,c]) *vm_fr[c]^2 -
                                    sum( b[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                                         g[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in _PMs.conductor_ids(pm) if d != c) -
                                    sum(-b[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                         g[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in _PMs.conductor_ids(pm)) )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_tp_ohms_yt_to(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: _PMs.AbstractACPForm
    p_to  = _PMs.var(pm, n, c,  :p, t_idx)
    q_to  = _PMs.var(pm, n, c,  :q, t_idx)
    vm_fr = [_PMs.var(pm, n, d, :vm, f_bus) for d in _PMs.conductor_ids(pm)]
    vm_to = [_PMs.var(pm, n, d, :vm, t_bus) for d in _PMs.conductor_ids(pm)]
    va_fr = [_PMs.var(pm, n, d, :va, f_bus) for d in _PMs.conductor_ids(pm)]
    va_to = [_PMs.var(pm, n, d, :va, t_bus) for d in _PMs.conductor_ids(pm)]

    JuMP.@NLconstraint(pm.model, p_to == (g_to[c]+g[c,c])*vm_to[c]^2 +
                                   sum( g[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) +
                                        b[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in _PMs.conductor_ids(pm) if d != c) +
                                   sum(-g[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                       -b[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in _PMs.conductor_ids(pm)) )
    JuMP.@NLconstraint(pm.model, q_to == -(b_to[c]+b[c,c])*vm_to[c]^2 -
                                   sum( b[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) -
                                        g[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in _PMs.conductor_ids(pm) if d != c) -
                                   sum(-b[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                        g[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in _PMs.conductor_ids(pm)) )
end


"""
```
p[f_idx] == z*(g/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
q[f_idx] == z*(-(b+c/2)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
```
"""
function constraint_tp_ohms_yt_from_on_off(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: _PMs.AbstractACPForm
    p_fr  = _PMs.var(pm, n, c,  :p, f_idx)
    q_fr  = _PMs.var(pm, n, c,  :q, f_idx)
    vm_fr = [_PMs.var(pm, n, d, :vm, f_bus) for d in _PMs.conductor_ids(pm)]
    vm_to = [_PMs.var(pm, n, d, :vm, t_bus) for d in _PMs.conductor_ids(pm)]
    va_fr = [_PMs.var(pm, n, d, :va, f_bus) for d in _PMs.conductor_ids(pm)]
    va_to = [_PMs.var(pm, n, d, :va, t_bus) for d in _PMs.conductor_ids(pm)]
    z = _PMs.var(pm, n, c, :branch_z, i)

    JuMP.@NLconstraint(pm.model, p_fr == z*((g_fr[c]+g[c,c]) * vm_fr[c]^2 +
                                        sum( g[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                                             b[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in _PMs.conductor_ids(pm) if d != c) +
                                        sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                            -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in _PMs.conductor_ids(pm))) )
    JuMP.@NLconstraint(pm.model, q_fr == z*(-(b_fr[c]+b[c,c]) *vm_fr[c]^2 -
                                        sum( b[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                                             g[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in _PMs.conductor_ids(pm) if d != c) -
                                        sum(-b[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                             g[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in _PMs.conductor_ids(pm))) )
end

"""
```
p[t_idx] == z*(g*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
q[t_idx] == z*(-(b+c/2)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
```
"""
function constraint_tp_ohms_yt_to_on_off(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: _PMs.AbstractACPForm
    p_to  = _PMs.var(pm, n, c,  :p, t_idx)
    q_to  = _PMs.var(pm, n, c,  :q, t_idx)
    vm_fr = [_PMs.var(pm, n, d, :vm, f_bus) for d in _PMs.conductor_ids(pm)]
    vm_to = [_PMs.var(pm, n, d, :vm, t_bus) for d in _PMs.conductor_ids(pm)]
    va_fr = [_PMs.var(pm, n, d, :va, f_bus) for d in _PMs.conductor_ids(pm)]
    va_to = [_PMs.var(pm, n, d, :va, t_bus) for d in _PMs.conductor_ids(pm)]
    z = _PMs.var(pm, n, c, :branch_z, i)

    JuMP.@NLconstraint(pm.model, p_to == z*((g_to[c]+g[c,c])*vm_to[c]^2 +
                                        sum( g[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) +
                                             b[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in _PMs.conductor_ids(pm) if d != c) +
                                        sum(-g[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                            -b[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in _PMs.conductor_ids(pm))) )
    JuMP.@NLconstraint(pm.model, q_to == z*(-(b_to[c]+b[c,c])*vm_to[c]^2 -
                                        sum( b[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) -
                                             g[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in _PMs.conductor_ids(pm) if d != c) -
                                        sum(-b[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                             g[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in _PMs.conductor_ids(pm))) )
end


"Links the voltage at both windings of a fixed tap transformer."
function constraint_tp_trans_voltage(pm::_PMs.GenericPowerModel{T}, nw::Int, i::Int, f_bus::Int, t_bus::Int, tm::_PMs.MultiConductorVector, Tv_fr, Tv_im, Cv_to) where T <: _PMs.AbstractACPForm
    ncnd  = 3
    # from side
    vm_fr = [_PMs.var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [_PMs.var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    # to side
    vm_to = [_PMs.var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
    va_to = [_PMs.var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
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
function constraint_tp_trans_flow(pm::_PMs.GenericPowerModel{T}, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, tm::_PMs.MultiConductorVector, Ti_fr, Ti_im, Cv_to) where T <: _PMs.AbstractACPForm
    ncnd  = 3
    # from side variables
    vm_fr = [_PMs.var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [_PMs.var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    p_fr = [_PMs.var(pm, nw, c, :pt, f_idx) for c in 1:ncnd]
    q_fr = [_PMs.var(pm, nw, c, :qt, f_idx) for c in 1:ncnd]
    # to side
    vm_to = [_PMs.var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
    va_to = [_PMs.var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
    p_to = [_PMs.var(pm, nw, c, :pt, t_idx) for c in 1:ncnd]
    q_to = [_PMs.var(pm, nw, c, :qt, t_idx) for c in 1:ncnd]
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
function constraint_tp_oltc_voltage(pm::_PMs.GenericPowerModel{T}, nw::Int, i::Int, f_bus::Int, t_bus::Int, Tv_fr, Tv_im, Cv_to) where T <: _PMs.AbstractACPForm
    ncnd  = 3
    # from side
    vm_fr = [_PMs.var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [_PMs.var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    # tap
    tap = [_PMs.var(pm, nw, c, :tap)[i] for c in 1:ncnd]
    # to side
    vm_to = [_PMs.var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
    va_to = [_PMs.var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
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
function constraint_tp_oltc_flow(pm::_PMs.GenericPowerModel{T}, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, Ti_fr, Ti_im, Cv_to) where T <: _PMs.AbstractACPForm
    ncnd  = 3
    # from side variables
    vm_fr = [_PMs.var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [_PMs.var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    p_fr = [_PMs.var(pm, nw, c, :pt, f_idx) for c in 1:ncnd]
    q_fr = [_PMs.var(pm, nw, c, :qt, f_idx) for c in 1:ncnd]
    # tap
    tap = [_PMs.var(pm, nw, c, :tap, i) for c in 1:ncnd]
    # to side
    vm_to = [_PMs.var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
    va_to = [_PMs.var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
    p_to = [_PMs.var(pm, nw, c, :pt, t_idx) for c in 1:ncnd]
    q_to = [_PMs.var(pm, nw, c, :qt, t_idx) for c in 1:ncnd]
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
