# Three-phase specific constraints


""
function variable_tp_voltage(pm::_PMs.GenericPowerModel{T}; kwargs...) where T <: _PMs.AbstractACRForm
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_voltage(pm, cnd=c; kwargs...)
    end
    # local infeasbility issues without proper initialization;
    # convergence issues start when the equivalent angles of the starting point
    # are further away than 90 degrees from the solution (as given by ACP)
    # this is the default behaviour of PMs, initialize all phases as (1,0)
    # the magnitude seems to have little effect on the convergence (>0.05)
    # updating the starting point to a balanced phasor does the job
    ncnd = length(_PMs.conductor_ids(pm))
    theta = [_wrap_to_pi(2 * pi / ncnd * (1-c)) for c in 1:ncnd]
    vm = 1
    for c in 1:ncnd
        vr = vm*cos(theta[c])
        vi = vm*sin(theta[c])
        for id in _PMs.ids(pm, :bus)
            JuMP.set_start_value(_PMs.var(pm, pm.cnw, c, :vr, id), vr)
            JuMP.set_start_value(_PMs.var(pm, pm.cnw, c, :vi, id), vi)
        end
    end
end


"delegate back to PowerModels"
function constraint_tp_model_voltage(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int) where T <: _PMs.AbstractACRForm
    _PMs.constraint_model_voltage(pm, n, c)
end


"Creates phase angle constraints at reference buses"
function constraint_tp_theta_ref(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, d) where T <: _PMs.AbstractACRForm
    vr = _PMs.var(pm, n, c, :vr, d)
    vi = _PMs.var(pm, n, c, :vi, d)
    nconductors = length(_PMs.conductor_ids(pm))
    theta = _wrap_to_pi(2 * pi / nconductors * (1-c))
    # deal with cases first where tan(theta)==Inf or tan(theta)==0
    if theta == pi/2
        JuMP.@constraint(pm.model, vr == 0)
        JuMP.@constraint(pm.model, vi >= 0)
    elseif theta == -pi/2
        JuMP.@constraint(pm.model, vr == 0)
        JuMP.@constraint(pm.model, vi <= 0)
    elseif theta == 0
        JuMP.@constraint(pm.model, vr >= 0)
        JuMP.@constraint(pm.model, vi == 0)
    elseif theta == pi
        JuMP.@constraint(pm.model, vr >= 0)
        JuMP.@constraint(pm.model, vi == 0)
    else
        JuMP.@constraint(pm.model, vi == tan(theta)*vr)
        # theta also implies a sign for vr, vi
        if 0<=theta && theta <= pi
            JuMP.@constraint(pm.model, vi >= 0)
        else
            JuMP.@constraint(pm.model, vi <= 0)
        end
    end
end


""
function constraint_tp_power_balance_shunt_slack(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: _PMs.AbstractACRForm
    vr = _PMs.var(pm, n, c, :vr, i)
    vi = _PMs.var(pm, n, c, :vi, i)
    p_slack = _PMs.var(pm, n, c, :p_slack, i)
    q_slack = _PMs.var(pm, n, c, :q_slack, i)
    p = _PMs.var(pm, n, c, :p)
    q = _PMs.var(pm, n, c, :q)
    pg = _PMs.var(pm, n, c, :pg)
    qg = _PMs.var(pm, n, c, :qg)
    p_dc = _PMs.var(pm, n, c, :p_dc)
    q_dc = _PMs.var(pm, n, c, :q_dc)

    _PMs.con(pm, n, c, :kcl_p)[i] = JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*(vr^2 + vi^2) + p_slack)
    _PMs.con(pm, n, c, :kcl_q)[i] = JuMP.@constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*(vr^2 + vi^2) + q_slack)
end


""
function constraint_tp_power_balance_shunt_trans(pm::_PMs.GenericPowerModel{T}, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: _PMs.AbstractACRForm
    vr = _PMs.var(pm, nw, c, :vr, i)
    vi = _PMs.var(pm, nw, c, :vi, i)
    p = _PMs.var(pm, nw, c, :p)
    q = _PMs.var(pm, nw, c, :q)
    pg = _PMs.var(pm, nw, c, :pg)
    qg = _PMs.var(pm, nw, c, :qg)
    p_dc = _PMs.var(pm, nw, c, :p_dc)
    q_dc = _PMs.var(pm, nw, c, :q_dc)
    p_trans = _PMs.var(pm, nw, c, :pt)
    q_trans = _PMs.var(pm,  nw, c, :qt)

    _PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) + sum(p_trans[a_trans] for a_trans in bus_arcs_trans) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*(vr^2 + vi^2))
    _PMs.con(pm, nw, c, :kcl_q)[i] = JuMP.@constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) + sum(q_trans[a_trans] for a_trans in bus_arcs_trans) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*(vr^2 + vi^2))
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p_fr =  g_fr[c]*(vr_fr[c]^2+vi_fr[c]^2)+ sum(
                                     vr_fr[c]*(g[c,d]*(vr_fr[d]-vr_to[d])-b[c,d]*(vi_fr[d]-vi_to[d]))
                                    -vi_fr[c]*(-b[c,d]*(vr_fr[d]-vr_to[d])-g[c,d]*(vi_fr[d]-vi_to[d]))
                                for d in _PMs.conductor_ids(pm))
q_fr =  -b_fr[c]*(vr_fr[c]^2+vi_fr[c]^2)+ sum(
                                        -vr_fr[c]*(b[c,d]*(vr_fr[d]-vr_to[d])+g[c,d]*(vi_fr[d]-vi_to[d]))
                                        +vi_fr[c]*(g[c,d]*(vr_fr[d]-vr_to[d])-b[c,d]*(vi_fr[d]-vi_to[d]))
                                    for d in _PMs.conductor_ids(pm))
```
"""
function constraint_tp_ohms_yt_from(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: _PMs.AbstractACRForm
    p_fr  = _PMs.var(pm, n, c,  :p, f_idx)
    q_fr  = _PMs.var(pm, n, c,  :q, f_idx)
    vr_fr = [_PMs.var(pm, n, d, :vr, f_bus) for d in _PMs.conductor_ids(pm)]
    vr_to = [_PMs.var(pm, n, d, :vr, t_bus) for d in _PMs.conductor_ids(pm)]
    vi_fr = [_PMs.var(pm, n, d, :vi, f_bus) for d in _PMs.conductor_ids(pm)]
    vi_to = [_PMs.var(pm, n, d, :vi, t_bus) for d in _PMs.conductor_ids(pm)]

    JuMP.@NLconstraint(pm.model, p_fr ==  g_fr[c]*(vr_fr[c]^2+vi_fr[c]^2)+ sum(
                                             vr_fr[c]*(g[c,d]*(vr_fr[d]-vr_to[d])-b[c,d]*(vi_fr[d]-vi_to[d]))
                                            -vi_fr[c]*(-b[c,d]*(vr_fr[d]-vr_to[d])-g[c,d]*(vi_fr[d]-vi_to[d]))
                                        for d in _PMs.conductor_ids(pm))
    )
    JuMP.@NLconstraint(pm.model, q_fr ==  -b_fr[c]*(vr_fr[c]^2+vi_fr[c]^2)+ sum(
                                            -vr_fr[c]*(b[c,d]*(vr_fr[d]-vr_to[d])+g[c,d]*(vi_fr[d]-vi_to[d]))
                                            +vi_fr[c]*(g[c,d]*(vr_fr[d]-vr_to[d])-b[c,d]*(vi_fr[d]-vi_to[d]))
                                        for d in _PMs.conductor_ids(pm))
    )
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_tp_ohms_yt_to(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: _PMs.AbstractACRForm
    constraint_tp_ohms_yt_from(pm, n, c, t_bus, f_bus, t_idx, f_idx, g, b, g_to, b_to, tr, ti, tm)
end
