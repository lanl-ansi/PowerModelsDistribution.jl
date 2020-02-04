""
function variable_mc_voltage(pm::_PMs.AbstractACRModel; nw=pm.cnw, bounded::Bool=true, kwargs...)
    variable_mc_voltage_real(pm; nw=nw, bounded=bounded, kwargs...)
    variable_mc_voltage_imaginary(pm; nw=nw, bounded=bounded, kwargs...)

    # local infeasbility issues without proper initialization;
    # convergence issues start when the equivalent angles of the starting point
    # are further away than 90 degrees from the solution (as given by ACP)
    # this is the default behaviour of _PMs, initialize all phases as (1,0)
    # the magnitude seems to have little effect on the convergence (>0.05)
    # updating the starting point to a balanced phasor does the job
    ncnd = length(_PMs.conductor_ids(pm))
    theta = [_wrap_to_pi(2 * pi / ncnd * (1-c)) for c in 1:ncnd]
    vm = 1
    for id in _PMs.ids(pm, nw, :bus)
        busref = _PMs.ref(pm, nw, :bus, id)
        if !haskey(busref, "va_start")
            for c in 1:ncnd
                vr = vm*cos(theta[c])
                vi = vm*sin(theta[c])
                JuMP.set_start_value(_PMs.var(pm, nw, :vr, id)[c], vr)
                JuMP.set_start_value(_PMs.var(pm, nw, :vi, id)[c], vi)
            end
        end
    end

    # apply bounds if bounded
    if bounded
        for i in _PMs.ids(pm, nw, :bus)
            constraint_mc_voltage_magnitude_bounds(pm, i, nw=nw)
        end
    end
end


"`vmin <= vm[i] <= vmax`"
function constraint_mc_voltage_magnitude_bounds(pm::_PMs.AbstractACRModel, n::Int, i, vmin, vmax)
    @assert all(vmin .<= vmax)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)

    JuMP.@constraint(pm.model, vmin.^2 .<= vr.^2 + vi.^2)
    JuMP.@constraint(pm.model, vmax.^2 .>= vr.^2 + vi.^2)
end


"Creates phase angle constraints at reference buses"
function constraint_mc_theta_ref(pm::_PMs.AbstractACRModel, n::Int, d::Int, va_ref)
    vr = _PMs.var(pm, n, :vr, d)
    vi = _PMs.var(pm, n, :vi, d)
    cnds = _PMs.conductor_ids(pm; nw=n)
    # deal with cases first where tan(theta)==Inf or tan(theta)==0

    for c in cnds
        if va_ref[c] == pi/2
            JuMP.@constraint(pm.model, vr[c] == 0)
            JuMP.@constraint(pm.model, vi[c] >= 0)
        elseif va_ref[c] == -pi/2
            JuMP.@constraint(pm.model, vr[c] == 0)
            JuMP.@constraint(pm.model, vi[c] <= 0)
        elseif va_ref[c] == 0
            JuMP.@constraint(pm.model, vr[c] >= 0)
            JuMP.@constraint(pm.model, vi[c] == 0)
        elseif va_ref[c] == pi
            JuMP.@constraint(pm.model, vr[c] >= 0)
            JuMP.@constraint(pm.model, vi[c] == 0)
        else
            JuMP.@constraint(pm.model, vi[c] == tan(va_ref[c])*vr[c])
            # va_ref also implies a sign for vr, vi
            if 0<=va_ref[c] && va_ref[c] <= pi
                JuMP.@constraint(pm.model, vi[c] >= 0)
            else
                JuMP.@constraint(pm.model, vi[c] <= 0)
            end
        end
    end
end

function constraint_mc_voltage_angle_difference(pm::_PMs.AbstractACRModel, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx

    vr_fr = _PMs.var(pm, n, :vr, f_bus)
    vi_fr = _PMs.var(pm, n, :vi, f_bus)
    vr_to = _PMs.var(pm, n, :vr, t_bus)
    vi_to = _PMs.var(pm, n, :vi, t_bus)

    for c in _PMs.conductor_ids(pm; nw=n)
        JuMP.@constraint(pm.model, (vi_fr[c]*vr_to[c] - vr_fr[c]*vi_to[c]) <= tan(angmax[c])*(vr_fr[c]*vr_to[c] + vi_fr[c]*vi_to[c]))
        JuMP.@constraint(pm.model, (vi_fr[c]*vr_to[c] - vr_fr[c]*vi_to[c]) >= tan(angmin[c])*(vr_fr[c]*vr_to[c] + vi_fr[c]*vi_to[c]))
    end
end



""
function constraint_mc_power_balance_slack(pm::_PMs.AbstractACRModel, nw::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vr = _PMs.var(pm, nw, :vr, i)
    vi = _PMs.var(pm, nw, :vi, i)
    p    = get(_PMs.var(pm, nw),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(_PMs.var(pm, nw),    :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(_PMs.var(pm, nw),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(_PMs.var(pm, nw),   :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PMs.var(pm, nw),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PMs.var(pm, nw),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(_PMs.var(pm, nw),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(_PMs.var(pm, nw),  :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(_PMs.var(pm, nw),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(_PMs.var(pm, nw),   :qt, Dict()); _PMs._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    p_slack = _PMs.var(pm, nw, :p_slack, i)
    q_slack = _PMs.var(pm, nw, :q_slack, i)

    cstr_p = []
    cstr_q = []

    for c in _PMs.conductor_ids(pm; nw=nw)
        cp = JuMP.@constraint(pm.model,
            sum(p[a][c] for a in bus_arcs)
            + sum(psw[a_sw][c] for a_sw in bus_arcs_sw)
            + sum(pt[a_trans][c] for a_trans in bus_arcs_trans)
            ==
            sum(pg[g][c] for g in bus_gens)
            - sum(ps[s][c] for s in bus_storage)
            - sum(pd[c] for pd in values(bus_pd))
            - sum(gs[c] for gs in values(bus_gs))*(vr[c]^2 + vi[c]^2)
            + p_slack[c]
        )
        push!(cstr_p, cp)

        cq = JuMP.@constraint(pm.model,
            sum(q[a][c] for a in bus_arcs)
            + sum(qsw[a_sw][c] for a_sw in bus_arcs_sw)
            + sum(qt[a_trans][c] for a_trans in bus_arcs_trans)
            ==
            sum(qg[g][c] for g in bus_gens)
            - sum(qs[s][c] for s in bus_storage)
            - sum(qd[c] for qd in values(bus_qd))
            + sum(bs[c] for bs in values(bus_bs))*(vr[c]^2 + vi[c]^2)
            + q_slack[c]
        )
        push!(cstr_q, cq)
    end

    if _PMs.report_duals(pm)
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_mc_power_balance(pm::_PMs.AbstractACRModel, nw::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vr = _PMs.var(pm, nw, :vr, i)
    vi = _PMs.var(pm, nw, :vi, i)
    p    = get(_PMs.var(pm, nw),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(_PMs.var(pm, nw),    :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(_PMs.var(pm, nw),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(_PMs.var(pm, nw),   :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PMs.var(pm, nw),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PMs.var(pm, nw),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(_PMs.var(pm, nw),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(_PMs.var(pm, nw),  :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(_PMs.var(pm, nw),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(_PMs.var(pm, nw),   :qt, Dict()); _PMs._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")

    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    Gt = isempty(bus_gs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_gs))
    Bt = isempty(bus_bs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_bs))

    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        .==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        # shunt
        - (vr.*(Gt*vr-Bt*vi) + vi.*(Gt*vi+Bt*vr))
    )

    cstr_q = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(qt[a_trans] for a_trans in bus_arcs_trans)
        .==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd for qd in values(bus_qd))
        # shunt
        - (-vr.*(Gt*vi+Bt*vr) + vi.*(Gt*vr-Bt*vi))
    )

    if _PMs.report_duals(pm)
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_mc_power_balance_load(pm::_PMs.AbstractACRModel, nw::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_gs, bus_bs)
    vr = _PMs.var(pm, nw, :vr, i)
    vi = _PMs.var(pm, nw, :vi, i)
    p    = get(_PMs.var(pm, nw),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(_PMs.var(pm, nw),    :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(_PMs.var(pm, nw),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(_PMs.var(pm, nw),   :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PMs.var(pm, nw),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PMs.var(pm, nw),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(_PMs.var(pm, nw),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(_PMs.var(pm, nw),  :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(_PMs.var(pm, nw),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(_PMs.var(pm, nw),   :qt, Dict()); _PMs._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    pd   = get(_PMs.var(pm, nw),  :pd, Dict()); _PMs._check_var_keys(pd, bus_loads, "active power", "load")
    qd   = get(_PMs.var(pm, nw),  :qd, Dict()); _PMs._check_var_keys(pd, bus_loads, "reactive power", "load")

    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    Gt = isempty(bus_gs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_gs))
    Bt = isempty(bus_bs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_bs))

    cstr_p = []
    cstr_q = []

    # pd/qd can be NLexpressions, so cannot be vectorized
    for c in _PMs.conductor_ids(pm; nw=nw)
        cp = JuMP.@NLconstraint(pm.model,
            sum(p[a][c] for a in bus_arcs)
            + sum(psw[a_sw][c] for a_sw in bus_arcs_sw)
            + sum(pt[a_trans][c] for a_trans in bus_arcs_trans)
            ==
            sum(pg[g][c] for g in bus_gens)
            - sum(ps[s][c] for s in bus_storage)
            - sum(pd[l][c] for l in bus_loads)
            - sum( # shunt
                   vr[c] * ( Gt[c,d]*vr[d] - Bt[c,d]*vi[d])
                  -vi[c] * (-Bt[c,d]*vr[d] - Gt[c,d]*vi[d])
              for d in cnds)
        )
        push!(cstr_p, cp)

        cq = JuMP.@NLconstraint(pm.model,
            sum(q[a][c] for a in bus_arcs)
            + sum(qsw[a_sw][c] for a_sw in bus_arcs_sw)
            + sum(qt[a_trans][c] for a_trans in bus_arcs_trans)
            ==
            sum(qg[g][c] for g in bus_gens)
            - sum(qs[s][c] for s in bus_storage)
            - sum(qd[l][c] for l in bus_loads)
            - sum( # shunt
                  -vr[c] * (Bt[c,d]*vr[d] + Gt[c,d]*vi[d])
                  +vi[c] * (Gt[c,d]*vr[d] - Bt[c,d]*vi[d])
              for d in cnds)
        )
        push!(cstr_q, cq)
    end

    if _PMs.report_duals(pm)
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


"""
Creates Ohms constraints

s_fr = v_fr.*conj(Y*(v_fr-v_to))
s_fr = (vr_fr+im*vi_fr).*(G-im*B)*([vr_fr-vr_to]-im*[vi_fr-vi_to])
s_fr = (vr_fr+im*vi_fr).*([G*vr_fr-G*vr_to-B*vi_fr+B*vi_to]-im*[G*vi_fr-G*vi_to+B*vr_fr-B*vr_to])
"""
function constraint_mc_ohms_yt_from(pm::_PMs.AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, G, B, G_fr, B_fr, tr, ti, tm)
    p_fr  = _PMs.var(pm, n, :p, f_idx)
    q_fr  = _PMs.var(pm, n, :q, f_idx)
    vr_fr = _PMs.var(pm, n, :vr, f_bus)
    vr_to = _PMs.var(pm, n, :vr, t_bus)
    vi_fr = _PMs.var(pm, n, :vi, f_bus)
    vi_to = _PMs.var(pm, n, :vi, t_bus)

    cnds = _PMs.conductor_ids(pm; nw=n)

    JuMP.@constraint(pm.model,
            p_fr .==  vr_fr.*(G*vr_fr-G*vr_to-B*vi_fr+B*vi_to)
                     +vi_fr.*(G*vi_fr-G*vi_to+B*vr_fr-B*vr_to)
                     # shunt
                     +vr_fr.*(G_fr*vr_fr-B_fr*vi_fr)
                     +vi_fr.*(G_fr*vi_fr+B_fr*vr_fr)
    )
    JuMP.@constraint(pm.model,
            q_fr .== -vr_fr.*(G*vi_fr-G*vi_to+B*vr_fr-B*vr_to)
                     +vi_fr.*(G*vr_fr-G*vr_to-B*vi_fr+B*vi_to)
                     # shunt
                     -vr_fr.*(G_fr*vi_fr+B_fr*vr_fr)
                     +vi_fr.*(G_fr*vr_fr-B_fr*vi_fr)
    )



    # for c in cnds
    #     JuMP.@NLconstraint(pm.model,
    #             p_fr[c] ==  sum(
    #                              vr_fr[c]*(g[c,d]*(vr_fr[d]-vr_to[d])-b[c,d]*(vi_fr[d]-vi_to[d]))
    #                             -vi_fr[c]*(-b[c,d]*(vr_fr[d]-vr_to[d])-g[c,d]*(vi_fr[d]-vi_to[d]))
    #                         for d in cnds)
    #                       # the shunt element is identical but vr_to=vi_to=0
    #                       + sum(
    #                              vr_fr[c]*(g_fr[c,d]*vr_fr[d]-b_fr[c,d]*vi_fr[d])
    #                             -vi_fr[c]*(-b_fr[c,d]*vr_fr[d]-g_fr[c,d]*vi_fr[d])
    #                         for d in cnds)
    #     )
    #     JuMP.@NLconstraint(pm.model,
    #             q_fr[c] ==  sum(
    #                             -vr_fr[c]*(b[c,d]*(vr_fr[d]-vr_to[d])+g[c,d]*(vi_fr[d]-vi_to[d]))
    #                             +vi_fr[c]*(g[c,d]*(vr_fr[d]-vr_to[d])-b[c,d]*(vi_fr[d]-vi_to[d]))
    #                         for d in cnds)
    #                       # the shunt element is identical but vr_to=vi_to=0
    #                       + sum(
    #                             -vr_fr[c]*(b_fr[c,d]*vr_fr[d]+g_fr[c,d]*vi_fr[d])
    #                             +vi_fr[c]*(g_fr[c,d]*vr_fr[d]-b_fr[c,d]*vi_fr[d])
    #                         for d in cnds)
    #     )
    # end
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_mc_ohms_yt_to(pm::_PMs.AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    constraint_mc_ohms_yt_from(pm, n, t_bus, f_bus, t_idx, f_idx, g, b, g_to, b_to, tr, ti, tm)
end


""
function constraint_mc_load_power_wye(pm::_PMs.AbstractACRModel, nw::Int, load_id::Int, pd, qd)
    _PMs.var(pm, nw, :pd)[load_id] = pd
    _PMs.var(pm, nw, :qd)[load_id] = qd
end


""
function constraint_mc_load_current_wye(pm::_PMs.AbstractACRModel, nw::Int, load_id::Int, load_bus_id::Int, cp, cq)
    vr = _PMs.var(pm, nw, :vr, load_bus_id)
    vi = _PMs.var(pm, nw, :vi, load_bus_id)

    ncnds = length(_PMs.conductor_ids(pm, nw))

    _PMs.var(pm, nw, :pd)[load_id] = JuMP.@NLexpression(pm.model, [c in 1:ncnds], cp[c]*sqrt(vr[c]^2+vi[c]^2))
    _PMs.var(pm, nw, :qd)[load_id] = JuMP.@NLexpression(pm.model, [c in 1:ncnds], cq[c]*sqrt(vr[c]^2+vi[c]^2))
end


""
function constraint_mc_load_impedance_wye(pm::_PMs.AbstractACRModel, nw::Int, load_id::Int, load_bus_id::Int, cp, cq)
    vr = _PMs.var(pm, nw, :vr, load_bus_id)
    vi = _PMs.var(pm, nw, :vi, load_bus_id)

    ncnds = length(_PMs.conductor_ids(pm, nw))

    _PMs.var(pm, nw, :pd)[load_id] = JuMP.@NLexpression(pm.model, [c in 1:ncnds], cp[c]*(vr[c]^2+vi[c]^2))
    _PMs.var(pm, nw, :qd)[load_id] = JuMP.@NLexpression(pm.model, [c in 1:ncnds], cq[c]*(vr[c]^2+vi[c]^2))
end


"""
For a delta load, sd = (s_ab, s_bc, s_ca), but we want to fix s = (s_a, s_b, s_c)
s is a non-linear transform of v and sd, s=f(v,sd)
s_a = v_a*conj(s_ab/(v_a-v_b) - s_ca/(v_c-v_a))
s_b = v_b*conj(s_ab/(v_a-v_b) - s_ca/(v_c-v_a))
s_c = v_c*conj(s_ab/(v_a-v_b) - s_ca/(v_c-v_a))
"""
function constraint_mc_load_power_delta(pm::_PMs.AbstractACRModel, nw::Int, load_id::Int, load_bus_id::Int, pd::Vector, qd::Vector)
    p_ab, p_bc, p_ca = pd
    q_ab, q_bc, q_ca = qd
    vre_a, vre_b, vre_c = _PMs.var(pm, nw, :vr, load_bus_id)
    vim_a, vim_b, vim_c = _PMs.var(pm, nw, :vi, load_bus_id)
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
    _PMs.var(pm, nw, :pd)[load_id] = [
        p_x(vre_a, vim_a, ire_ab, iim_ab, ire_ca, iim_ca),
        p_x(vre_b, vim_b, ire_bc, iim_bc, ire_ab, iim_ab),
        p_x(vre_c, vim_c, ire_ca, iim_ca, ire_bc, iim_bc)]
    _PMs.var(pm, nw, :qd)[load_id] = [
        q_x(vre_a, vim_a, ire_ab, iim_ab, ire_ca, iim_ca),
        q_x(vre_b, vim_b, ire_bc, iim_bc, ire_ab, iim_ab),
        q_x(vre_c, vim_c, ire_ca, iim_ca, ire_bc, iim_bc)]

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
function constraint_mc_load_current_delta(pm::_PMs.AbstractACRModel, nw::Int, load_id::Int, load_bus_id::Int, cp::Vector, cq::Vector)
    cp_ab, cp_bc, cp_ca = cp
    cq_ab, cq_bc, cq_ca = cq
    vre_a, vre_b, vre_c = _PMs.var(pm, nw, :vr, load_bus_id)
    vim_a, vim_b, vim_c = _PMs.var(pm, nw, :vi, load_bus_id)
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
    _PMs.var(pm, nw, :pd)[load_id] = [
        p_x(vre_a, vim_a, ire_ab, iim_ab, ire_ca, iim_ca),
        p_x(vre_b, vim_b, ire_bc, iim_bc, ire_ab, iim_ab),
        p_x(vre_c, vim_c, ire_ca, iim_ca, ire_bc, iim_bc)]
    _PMs.var(pm, nw, :qd)[load_id] = [
        q_x(vre_a, vim_a, ire_ab, iim_ab, ire_ca, iim_ca),
        q_x(vre_b, vim_b, ire_bc, iim_bc, ire_ab, iim_ab),
        q_x(vre_c, vim_c, ire_ca, iim_ca, ire_bc, iim_bc)]
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
function constraint_mc_load_impedance_delta(pm::_PMs.AbstractACRModel, nw::Int, load_id::Int, load_bus_id::Int, cp::Vector, cq::Vector)
    cp_ab, cp_bc, cp_ca = cp
    cq_ab, cq_bc, cq_ca = cq
    vre_a, vre_b, vre_c = _PMs.var(pm, nw, :vr, load_bus_id)
    vim_a, vim_b, vim_c = _PMs.var(pm, nw, :vi, load_bus_id)
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
    _PMs.var(pm, nw, :pd)[load_id] = [
        p_x(vre_a, vim_a, ire_ab, iim_ab, ire_ca, iim_ca),
        p_x(vre_b, vim_b, ire_bc, iim_bc, ire_ab, iim_ab),
        p_x(vre_c, vim_c, ire_ca, iim_ca, ire_bc, iim_bc)]
    _PMs.var(pm, nw, :qd)[load_id] = [
        q_x(vre_a, vim_a, ire_ab, iim_ab, ire_ca, iim_ca),
        q_x(vre_b, vim_b, ire_bc, iim_bc, ire_ab, iim_ab),
        q_x(vre_c, vim_c, ire_ca, iim_ca, ire_bc, iim_bc)]
end


"`vm[i] == vmref`"
function constraint_mc_voltage_magnitude_setpoint(pm::_PMs.AbstractACRModel, n::Int, i::Int, vmref)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)

    JuMP.@constraint(pm.model, vr.^2 + vi.^2  .== vmref.^2)
end
