# Three-phase specific constraints


""
function variable_tp_voltage(pm::PMs.GenericPowerModel{T}; nw=pm.cnw, kwargs...) where T <: PMs.AbstractACRForm
    for c in PMs.conductor_ids(pm)
        PMs.variable_voltage(pm, cnd=c; nw=nw, kwargs...)
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
    for id in PMs.ids(pm, nw, :bus)
        busref = PMs.ref(pm, nw, :bus, id)
        if !haskey(busref, "va_start")
            for c in 1:ncnd
                vr = vm*cos(theta[c])
                vi = vm*sin(theta[c])
                JuMP.set_start_value(PMs.var(pm, nw, c, :vr, id), vr)
                JuMP.set_start_value(PMs.var(pm, nw, c, :vi, id), vi)
            end
        end
    end
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
    pt = PMs.var(pm, nw, c, :pt)
    qt = PMs.var(pm,  nw, c, :qt)

    PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) + sum(pt[a_trans] for a_trans in bus_arcs_trans) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*(vr^2 + vi^2))
    PMs.con(pm, nw, c, :kcl_q)[i] = JuMP.@constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) + sum(qt[a_trans] for a_trans in bus_arcs_trans) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*(vr^2 + vi^2))
end


""
function constraint_kcl_shunt_trans_load(pm::PMs.GenericPowerModel{T}, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_loads, bus_gs, bus_bs) where T <: PMs.AbstractACRForm
    vr = PMs.var(pm, nw, c, :vr, i)
    vi = PMs.var(pm, nw, c, :vi, i)
    p = PMs.var(pm, nw, c, :p)
    q = PMs.var(pm, nw, c, :q)
    pg = PMs.var(pm, nw, c, :pg)
    qg = PMs.var(pm, nw, c, :qg)
    pd = PMs.var(pm, nw, c, :pd)
    qd = PMs.var(pm, nw, c, :qd)
    p_dc = PMs.var(pm, nw, c, :p_dc)
    q_dc = PMs.var(pm, nw, c, :q_dc)
    pt = PMs.var(pm, nw, c, :pt)
    qt = PMs.var(pm,  nw, c, :qt)

    PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@NLconstraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) + sum(pt[a_trans] for a_trans in bus_arcs_trans) == sum(pg[g] for g in bus_gens) - sum(pd[l] for l in bus_loads) - sum(gs for gs in values(bus_gs))*(vr^2 + vi^2))
    PMs.con(pm, nw, c, :kcl_q)[i] = JuMP.@NLconstraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) + sum(qt[a_trans] for a_trans in bus_arcs_trans) == sum(qg[g] for g in bus_gens) - sum(qd[l] for l in bus_loads) + sum(bs for bs in values(bus_bs))*(vr^2 + vi^2))
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


function constraint_load_power_wye(pm::PMs.GenericPowerModel{T}, nw::Int, cnd::Int, load_id::Int, pd::Real, qd::Real) where T <: PMs.AbstractACRForm
    PMs.var(pm, nw, cnd, :pd)[load_id] = pd
    PMs.var(pm, nw, cnd, :qd)[load_id] = qd
end


function constraint_load_current_wye(pm::PMs.GenericPowerModel{T}, nw::Int, cnd::Int, load_id::Int, load_bus_id::Int, cp::Real, cq::Real) where T <: PMs.AbstractACRForm
    vr = PMs.var(pm, nw, cnd, :vr, load_bus_id)
    vi = PMs.var(pm, nw, cnd, :vi, load_bus_id)
    PMs.var(pm, nw, cnd, :pd)[load_id] = JuMP.@NLexpression(pm.model, cp*sqrt(vr^2+vi^2))
    PMs.var(pm, nw, cnd, :qd)[load_id] = JuMP.@NLexpression(pm.model, cq*sqrt(vr^2+vi^2))
end


function constraint_load_impedance_wye(pm::PMs.GenericPowerModel{T}, nw::Int, cnd::Int, load_id::Int, load_bus_id::Int, cp::Real, cq::Real) where T <: PMs.AbstractACRForm
    vr = PMs.var(pm, nw, cnd, :vr, load_bus_id)
    vi = PMs.var(pm, nw, cnd, :vi, load_bus_id)
    PMs.var(pm, nw, cnd, :pd)[load_id] = JuMP.@NLexpression(pm.model, cp*(vr^2+vi^2))
    PMs.var(pm, nw, cnd, :qd)[load_id] = JuMP.@NLexpression(pm.model, cq*(vr^2+vi^2))
end


"""
For a delta load, sd = (s_ab, s_bc, s_ca), but we want to fix s = (s_a, s_b, s_c)
s is a non-linear transform of v and sd, s=f(v,sd)
s_a = v_a*conj(s_ab/(v_a-v_b) - s_ca/(v_c-v_a))
s_b = v_b*conj(s_ab/(v_a-v_b) - s_ca/(v_c-v_a))
s_c = v_c*conj(s_ab/(v_a-v_b) - s_ca/(v_c-v_a))
"""
function constraint_tp_load_power_delta(pm::PMs.GenericPowerModel{T}, nw::Int, load_id::Int, load_bus_id::Int, pd::PMs.MultiConductorVector, qd::PMs.MultiConductorVector) where T <: PMs.AbstractACRForm
    p_ab, p_bc, p_ca = pd
    q_ab, q_bc, q_ca = qd
    vre_a, vre_b, vre_c = [PMs.var(pm, nw, c, :vr, load_bus_id) for c in 1:3]
    vim_a, vim_b, vim_c = [PMs.var(pm, nw, c, :vi, load_bus_id) for c in 1:3]
    # v_xy = v_x - v_y
    vre_xy(vre_x, vre_y) = JuMP.@NLexpression(pm.model, vre_x-vre_y)
    vim_xy(vim_x, vim_y) = JuMP.@NLexpression(pm.model, vim_x-vim_y)
    vre_ab = vre_xy(vre_a, vre_b)
    vim_ab = vim_xy(vim_a, vim_b)
    vre_bc = vre_xy(vre_b, vre_c)
    vim_bc = vim_xy(vim_b, vim_c)
    vre_ca = vre_xy(vre_c, vre_a)
    vim_ca = vim_xy(vim_c, vim_a)
    # i_xy = conj(s_xy/v_xy)
    ire_xy(p_xy, q_xy, vre_xy, vim_xy) = JuMP.@NLexpression(pm.model, (p_xy*vre_xy+q_xy*vim_xy)/(vre_xy^2+vim_xy^2))
    iim_xy(p_xy, q_xy, vre_xy, vim_xy) = JuMP.@NLexpression(pm.model, (p_xy*vim_xy-q_xy*vre_xy)/(vre_xy^2+vim_xy^2))
    ire_ab = ire_xy(p_ab, q_ab, vre_ab, vim_ab)
    iim_ab = iim_xy(p_ab, q_ab, vre_ab, vim_ab)
    ire_bc = ire_xy(p_bc, q_bc, vre_bc, vim_bc)
    iim_bc = iim_xy(p_bc, q_bc, vre_bc, vim_bc)
    ire_ca = ire_xy(p_ca, q_ca, vre_ca, vim_ca)
    iim_ca = iim_xy(p_ca, q_ca, vre_ca, vim_ca)
    # s_x = v_x*conj(i_xy-i_zx)
    # p_x = vm_x*cos(va_x)*(ire_xy-ire_zx) + vm_x*sin(va_x)*(iim_xy-iim_zx)
    # q_x = vm_x*sin(va_x)*(ire_xy-ire_zx) - vm_x*cos(va_x)*(iim_xy-iim_zx)
    p_x(vre_x, vim_x, ire_xy, iim_xy, ire_zx, iim_zx) = JuMP.@NLexpression(pm.model, vre_x*(ire_xy-ire_zx) + vim_x*(iim_xy-iim_zx))
    q_x(vre_x, vim_x, ire_xy, iim_xy, ire_zx, iim_zx) = JuMP.@NLexpression(pm.model, vim_x*(ire_xy-ire_zx) - vre_x*(iim_xy-iim_zx))
    # s_x = s_x,ref
    PMs.var(pm, nw, 1, :pd)[load_id] = p_x(vre_a, vim_a, ire_ab, iim_ab, ire_ca, iim_ca)
    PMs.var(pm, nw, 2, :pd)[load_id] = p_x(vre_b, vim_b, ire_bc, iim_bc, ire_ab, iim_ab)
    PMs.var(pm, nw, 3, :pd)[load_id] = p_x(vre_c, vim_c, ire_ca, iim_ca, ire_bc, iim_bc)
    PMs.var(pm, nw, 1, :qd)[load_id] = q_x(vre_a, vim_a, ire_ab, iim_ab, ire_ca, iim_ca)
    PMs.var(pm, nw, 2, :qd)[load_id] = q_x(vre_b, vim_b, ire_bc, iim_bc, ire_ab, iim_ab)
    PMs.var(pm, nw, 3, :qd)[load_id] = q_x(vre_c, vim_c, ire_ca, iim_ca, ire_bc, iim_bc)

end


"""
We want to express
s_ab = cp.|v_ab|+im.cq.|v_ab|
i_ab = conj(s_ab/v_ab) = |v_ab|.(cq-im.cq)/conj(v_ab) = (1/|v_ab|).(cp-im.cq)*v_ab
idem for i_bc and i_ca
And then
s_a = v_a.conj(i_a) = v_a.conj(i_ab-i_ca)
idem for s_b and s_c
"""
function constraint_tp_load_current_delta(pm::PMs.GenericPowerModel{T}, nw::Int, load_id::Int, load_bus_id::Int, cp::PMs.MultiConductorVector, cq::PMs.MultiConductorVector) where T <: PMs.AbstractACRForm
    cp_ab, cp_bc, cp_ca = cp
    cq_ab, cq_bc, cq_ca = cq
    vre_a, vre_b, vre_c = [PMs.var(pm, nw, c, :vr, load_bus_id) for c in 1:3]
    vim_a, vim_b, vim_c = [PMs.var(pm, nw, c, :vi, load_bus_id) for c in 1:3]
    # v_xy = v_x - v_y
    vre_xy(vre_x, vre_y) = JuMP.@NLexpression(pm.model, vre_x-vre_y)
    vim_xy(vim_x, vim_y) = JuMP.@NLexpression(pm.model, vim_x-vim_y)
    vre_ab = vre_xy(vre_a, vre_b)
    vim_ab = vim_xy(vim_a, vim_b)
    vre_bc = vre_xy(vre_b, vre_c)
    vim_bc = vim_xy(vim_b, vim_c)
    vre_ca = vre_xy(vre_c, vre_a)
    vim_ca = vim_xy(vim_c, vim_a)
    # i_xy = conj(s_xy/v_xy)
    ire_xy(cp_xy, cq_xy, vre_xy, vim_xy) = JuMP.@NLexpression(pm.model, 1/sqrt(vre_xy^2+vim_xy^2)*(cp_xy*vre_xy+cq_xy*vim_xy))
    iim_xy(cp_xy, cq_xy, vre_xy, vim_xy) = JuMP.@NLexpression(pm.model, 1/sqrt(vre_xy^2+vim_xy^2)*(cp_xy*vim_xy-cq_xy*vre_xy))
    ire_ab = ire_xy(cp_ab, cq_ab, vre_ab, vim_ab)
    iim_ab = iim_xy(cp_ab, cq_ab, vre_ab, vim_ab)
    ire_bc = ire_xy(cp_bc, cq_bc, vre_bc, vim_bc)
    iim_bc = iim_xy(cp_bc, cq_bc, vre_bc, vim_bc)
    ire_ca = ire_xy(cp_ca, cq_ca, vre_ca, vim_ca)
    iim_ca = iim_xy(cp_ca, cq_ca, vre_ca, vim_ca)
    # s_x = v_x*conj(i_xy-i_zx)
    # p_x = vm_x*cos(va_x)*(ire_xy-ire_zx) + vm_x*sin(va_x)*(iim_xy-iim_zx)
    # q_x = vm_x*sin(va_x)*(ire_xy-ire_zx) - vm_x*cos(va_x)*(iim_xy-iim_zx)
    p_x(vre_x, vim_x, ire_xy, iim_xy, ire_zx, iim_zx) = JuMP.@NLexpression(pm.model, vre_x*(ire_xy-ire_zx) + vim_x*(iim_xy-iim_zx))
    q_x(vre_x, vim_x, ire_xy, iim_xy, ire_zx, iim_zx) = JuMP.@NLexpression(pm.model, vim_x*(ire_xy-ire_zx) - vre_x*(iim_xy-iim_zx))
    # s_x = s_x,ref
    PMs.var(pm, nw, 1, :pd)[load_id] = p_x(vre_a, vim_a, ire_ab, iim_ab, ire_ca, iim_ca)
    PMs.var(pm, nw, 2, :pd)[load_id] = p_x(vre_b, vim_b, ire_bc, iim_bc, ire_ab, iim_ab)
    PMs.var(pm, nw, 3, :pd)[load_id] = p_x(vre_c, vim_c, ire_ca, iim_ca, ire_bc, iim_bc)
    PMs.var(pm, nw, 1, :qd)[load_id] = q_x(vre_a, vim_a, ire_ab, iim_ab, ire_ca, iim_ca)
    PMs.var(pm, nw, 2, :qd)[load_id] = q_x(vre_b, vim_b, ire_bc, iim_bc, ire_ab, iim_ab)
    PMs.var(pm, nw, 3, :qd)[load_id] = q_x(vre_c, vim_c, ire_ca, iim_ca, ire_bc, iim_bc)
end


"""
We want to express
s_ab = cp.|v_ab|^2+im.cq.|v_ab|^2
i_ab = conj(s_ab/v_ab) = |v_ab|^2.(cq-im.cq)/conj(v_ab) = (cp-im.cq)*v_ab
idem for i_bc and i_ca
And then
s_a = v_a.conj(i_a) = v_a.conj(i_ab-i_ca)
idem for s_b and s_c
"""
function constraint_tp_load_impedance_delta(pm::PMs.GenericPowerModel{T}, nw::Int, load_id::Int, load_bus_id::Int, cp::PMs.MultiConductorVector, cq::PMs.MultiConductorVector) where T <: PMs.AbstractACRForm
    cp_ab, cp_bc, cp_ca = cp
    cq_ab, cq_bc, cq_ca = cq
    vre_a, vre_b, vre_c = [PMs.var(pm, nw, c, :vr, load_bus_id) for c in 1:3]
    vim_a, vim_b, vim_c = [PMs.var(pm, nw, c, :vi, load_bus_id) for c in 1:3]
    # v_xy = v_x - v_y
    vre_xy(vre_x, vre_y) = JuMP.@NLexpression(pm.model, vre_x-vre_y)
    vim_xy(vim_x, vim_y) = JuMP.@NLexpression(pm.model, vim_x-vim_y)
    vre_ab = vre_xy(vre_a, vre_b)
    vim_ab = vim_xy(vim_a, vim_b)
    vre_bc = vre_xy(vre_b, vre_c)
    vim_bc = vim_xy(vim_b, vim_c)
    vre_ca = vre_xy(vre_c, vre_a)
    vim_ca = vim_xy(vim_c, vim_a)
    # i_xy = conj(s_xy/v_xy)
    ire_xy(cp_xy, cq_xy, vre_xy, vim_xy) = JuMP.@NLexpression(pm.model, cp_xy*vre_xy+cq_xy*vim_xy)
    iim_xy(cp_xy, cq_xy, vre_xy, vim_xy) = JuMP.@NLexpression(pm.model, cp_xy*vim_xy-cq_xy*vre_xy)
    ire_ab = ire_xy(cp_ab, cq_ab, vre_ab, vim_ab)
    iim_ab = iim_xy(cp_ab, cq_ab, vre_ab, vim_ab)
    ire_bc = ire_xy(cp_bc, cq_bc, vre_bc, vim_bc)
    iim_bc = iim_xy(cp_bc, cq_bc, vre_bc, vim_bc)
    ire_ca = ire_xy(cp_ca, cq_ca, vre_ca, vim_ca)
    iim_ca = iim_xy(cp_ca, cq_ca, vre_ca, vim_ca)
    # s_x = v_x*conj(i_xy-i_zx)
    # p_x = vm_x*cos(va_x)*(ire_xy-ire_zx) + vm_x*sin(va_x)*(iim_xy-iim_zx)
    # q_x = vm_x*sin(va_x)*(ire_xy-ire_zx) - vm_x*cos(va_x)*(iim_xy-iim_zx)
    p_x(vre_x, vim_x, ire_xy, iim_xy, ire_zx, iim_zx) = JuMP.@NLexpression(pm.model, vre_x*(ire_xy-ire_zx) + vim_x*(iim_xy-iim_zx))
    q_x(vre_x, vim_x, ire_xy, iim_xy, ire_zx, iim_zx) = JuMP.@NLexpression(pm.model, vim_x*(ire_xy-ire_zx) - vre_x*(iim_xy-iim_zx))
    # s_x = s_x,ref
    PMs.var(pm, nw, 1, :pd)[load_id] = p_x(vre_a, vim_a, ire_ab, iim_ab, ire_ca, iim_ca)
    PMs.var(pm, nw, 2, :pd)[load_id] = p_x(vre_b, vim_b, ire_bc, iim_bc, ire_ab, iim_ab)
    PMs.var(pm, nw, 3, :pd)[load_id] = p_x(vre_c, vim_c, ire_ca, iim_ca, ire_bc, iim_bc)
    PMs.var(pm, nw, 1, :qd)[load_id] = q_x(vre_a, vim_a, ire_ab, iim_ab, ire_ca, iim_ca)
    PMs.var(pm, nw, 2, :qd)[load_id] = q_x(vre_b, vim_b, ire_bc, iim_bc, ire_ab, iim_ab)
    PMs.var(pm, nw, 3, :qd)[load_id] = q_x(vre_c, vim_c, ire_ca, iim_ca, ire_bc, iim_bc)
end
