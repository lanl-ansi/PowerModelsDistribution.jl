# Three-phase specific constraints


""
function variable_tp_voltage(pm::PMs.GenericPowerModel{T}; kwargs...) where T <: PMs.AbstractACRForm
    for c in PMs.conductor_ids(pm)
        PMs.variable_voltage(pm, cnd=c; kwargs...)
    end
    # local infeasbility issues without proper initialization;
    # convergence issues start when the equivalent angles of the starting point
    # are further away than 90 degrees from the solution (as given by ACP)
    # this is the default behaviour of PMs, initialize all phases as (1,0)
    # the magnitude seems to have little effect on the convergence (>0.05)
    # updating the starting point to a balanced phasor does the job
    ncnd = length(PMs.conductor_ids(pm))
    theta = [wraptopi(2 * pi / ncnd * (1-c)) for c in 1:ncnd]
    vm = 1
    for c in 1:ncnd
        vr = vm*cos(theta[c])
        vi = vm*sin(theta[c])
        for id in PMs.ids(pm, :bus)
            JuMP.setvalue(PMs.var(pm, pm.cnw, c, :vr, id), vr)
            JuMP.setvalue(PMs.var(pm, pm.cnw, c, :vi, id), vi)
        end
    end
end


"delegate back to PowerModels"
function constraint_tp_voltage(pm::PMs.GenericPowerModel{T}, n::Int, c::Int) where T <: PMs.AbstractACRForm
    PMs.constraint_voltage(pm, n, c)
end


"Creates phase angle constraints at reference buses"
function constraint_tp_theta_ref(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, d) where T <: PMs.AbstractACRForm
    vr = PMs.var(pm, n, c, :vr, d)
    vi = PMs.var(pm, n, c, :vi, d)
    nconductors = length(PMs.conductor_ids(pm))
    theta = wraptopi(2 * pi / nconductors * (1-c))
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
function constraint_kcl_shunt_slack(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: PMs.AbstractACRForm
    vr = PMs.var(pm, n, c, :vr, i)
    vi = PMs.var(pm, n, c, :vi, i)
    p_slack = PMs.var(pm, n, c, :p_slack, i)
    q_slack = PMs.var(pm, n, c, :q_slack, i)
    p = PMs.var(pm, n, c, :p)
    q = PMs.var(pm, n, c, :q)
    pg = PMs.var(pm, n, c, :pg)
    qg = PMs.var(pm, n, c, :qg)
    p_dc = PMs.var(pm, n, c, :p_dc)
    q_dc = PMs.var(pm, n, c, :q_dc)

    PMs.con(pm, n, c, :kcl_p)[i] = JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*(vr^2 + vi^2) + p_slack)
    PMs.con(pm, n, c, :kcl_q)[i] = JuMP.@constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*(vr^2 + vi^2) + q_slack)
end


""
function constraint_kcl_shunt_trans(pm::PMs.GenericPowerModel{T}, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: PMs.AbstractACRForm
    vr = PMs.var(pm, nw, c, :vr, i)
    vi = PMs.var(pm, nw, c, :vi, i)
    p = PMs.var(pm, nw, c, :p)
    q = PMs.var(pm, nw, c, :q)
    pg = PMs.var(pm, nw, c, :pg)
    qg = PMs.var(pm, nw, c, :qg)
    p_dc = PMs.var(pm, nw, c, :p_dc)
    q_dc = PMs.var(pm, nw, c, :q_dc)
    p_trans = PMs.var(pm, nw, c, :pt)
    q_trans = PMs.var(pm,  nw, c, :qt)

    PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) + sum(p_trans[a_trans] for a_trans in bus_arcs_trans) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*(vr^2 + vi^2))
    PMs.con(pm, nw, c, :kcl_q)[i] = JuMP.@constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) + sum(q_trans[a_trans] for a_trans in bus_arcs_trans) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*(vr^2 + vi^2))
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p_fr =  g_fr[c]*(vr_fr[c]^2+vi_fr[c]^2)+ sum(
                                     vr_fr[c]*(g[c,d]*(vr_fr[d]-vr_to[d])-b[c,d]*(vi_fr[d]-vi_to[d]))
                                    -vi_fr[c]*(-b[c,d]*(vr_fr[d]-vr_to[d])-g[c,d]*(vi_fr[d]-vi_to[d]))
                                for d in PMs.conductor_ids(pm))
q_fr =  -b_fr[c]*(vr_fr[c]^2+vi_fr[c]^2)+ sum(
                                        -vr_fr[c]*(b[c,d]*(vr_fr[d]-vr_to[d])+g[c,d]*(vi_fr[d]-vi_to[d]))
                                        +vi_fr[c]*(g[c,d]*(vr_fr[d]-vr_to[d])-b[c,d]*(vi_fr[d]-vi_to[d]))
                                    for d in PMs.conductor_ids(pm))
```
"""
function constraint_ohms_tp_yt_from(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractACRForm
    p_fr  = PMs.var(pm, n, c,  :p, f_idx)
    q_fr  = PMs.var(pm, n, c,  :q, f_idx)
    vr_fr = [PMs.var(pm, n, d, :vr, f_bus) for d in PMs.conductor_ids(pm)]
    vr_to = [PMs.var(pm, n, d, :vr, t_bus) for d in PMs.conductor_ids(pm)]
    vi_fr = [PMs.var(pm, n, d, :vi, f_bus) for d in PMs.conductor_ids(pm)]
    vi_to = [PMs.var(pm, n, d, :vi, t_bus) for d in PMs.conductor_ids(pm)]

    JuMP.@NLconstraint(pm.model, p_fr ==  g_fr[c]*(vr_fr[c]^2+vi_fr[c]^2)+ sum(
                                             vr_fr[c]*(g[c,d]*(vr_fr[d]-vr_to[d])-b[c,d]*(vi_fr[d]-vi_to[d]))
                                            -vi_fr[c]*(-b[c,d]*(vr_fr[d]-vr_to[d])-g[c,d]*(vi_fr[d]-vi_to[d]))
                                        for d in PMs.conductor_ids(pm))
    )
    JuMP.@NLconstraint(pm.model, q_fr ==  -b_fr[c]*(vr_fr[c]^2+vi_fr[c]^2)+ sum(
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
function constraint_ohms_tp_yt_to(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractACRForm
    constraint_ohms_tp_yt_from(pm, n, c, t_bus, f_bus, t_idx, f_idx, g, b, g_to, b_to, tr, ti, tm)
end
