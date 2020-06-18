""
function variable_mc_bus_voltage(pm::_PM.AbstractACRModel; nw=pm.cnw, bounded::Bool=true, kwargs...)
    variable_mc_bus_voltage_real(pm; nw=nw, bounded=bounded, kwargs...)
    variable_mc_bus_voltage_imaginary(pm; nw=nw, bounded=bounded, kwargs...)

    # local infeasbility issues without proper initialization;
    # convergence issues start when the equivalent angles of the starting point
    # are further away than 90 degrees from the solution (as given by ACP)
    # this is the default behaviour of _PM, initialize all phases as (1,0)
    # the magnitude seems to have little effect on the convergence (>0.05)
    # updating the starting point to a balanced phasor does the job
    ncnd = length(conductor_ids(pm))
    theta = [_wrap_to_pi(2 * pi / ncnd * (1-c)) for c in 1:ncnd]
    vm = haskey(ref(pm, nw, :bus, id), "vm_start") ? busref["vm_start"] : 1.0
    for id in ids(pm, nw, :bus)
        busref = ref(pm, nw, :bus, id)
        if !haskey(busref, "va_start")
            for c in 1:ncnd
                vr = vm*cos(theta[c])
                vi = vm*sin(theta[c])
                JuMP.set_start_value(var(pm, nw, :vr, id)[c], vr)
                JuMP.set_start_value(var(pm, nw, :vi, id)[c], vi)
            end
        end
    end

    # apply bounds if bounded
    if bounded
        for i in ids(pm, nw, :bus)
            constraint_mc_voltage_magnitude_bounds(pm, i, nw=nw)
        end
    end
end


"`vmin <= vm[i] <= vmax`"
function constraint_mc_voltage_magnitude_bounds(pm::_PM.AbstractACRModel, n::Int, i, vmin, vmax)
    @assert all(vmin .<= vmax)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    for t in 1:length(vr)
        JuMP.@constraint(pm.model, vmin[t]^2 .<= vr[t]^2 + vi[t]^2)
        if vmax[t] < Inf
            JuMP.@constraint(pm.model, vmax[t]^2 .>= vr[t]^2 + vi[t]^2)
        end
    end
end


"Creates phase angle constraints at reference buses"
function constraint_mc_theta_ref(pm::_PM.AbstractACRModel, n::Int, d::Int, va_ref)
    vr = var(pm, n, :vr, d)
    vi = var(pm, n, :vi, d)
    cnds = conductor_ids(pm; nw=n)
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


""
function constraint_mc_voltage_angle_difference(pm::_PM.AbstractACRModel, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx

    vr_fr = var(pm, n, :vr, f_bus)
    vi_fr = var(pm, n, :vi, f_bus)
    vr_to = var(pm, n, :vr, t_bus)
    vi_to = var(pm, n, :vi, t_bus)

    for c in conductor_ids(pm; nw=n)
        JuMP.@constraint(pm.model, (vi_fr[c]*vr_to[c] - vr_fr[c]*vi_to[c]) <= tan(angmax[c])*(vr_fr[c]*vr_to[c] + vi_fr[c]*vi_to[c]))
        JuMP.@constraint(pm.model, (vi_fr[c]*vr_to[c] - vr_fr[c]*vi_to[c]) >= tan(angmin[c])*(vr_fr[c]*vr_to[c] + vi_fr[c]*vi_to[c]))
    end
end


""
function constraint_mc_slack_power_balance(pm::_PM.AbstractACRModel, nw::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)
    p    = get(var(pm, nw),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(var(pm, nw),    :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, nw),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, nw),   :qg, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, nw),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, nw),   :qs, Dict()); _PM._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, nw),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, nw),  :qsw, Dict()); _PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(var(pm, nw),   :pt, Dict()); _PM._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(var(pm, nw),   :qt, Dict()); _PM._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    p_slack = var(pm, nw, :p_slack, i)
    q_slack = var(pm, nw, :q_slack, i)

    cstr_p = []
    cstr_q = []

    for c in conductor_ids(pm; nw=nw)
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

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = cstr_q

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_mc_power_balance(pm::_PM.AbstractACRModel, nw::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)
    p    = get(var(pm, nw),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(var(pm, nw),    :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, nw),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, nw),   :qg, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, nw),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, nw),   :qs, Dict()); _PM._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, nw),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, nw),  :qsw, Dict()); _PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(var(pm, nw),   :pt, Dict()); _PM._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(var(pm, nw),   :qt, Dict()); _PM._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")

    cnds = conductor_ids(pm; nw=nw)
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

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = cstr_q

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_mc_load_power_balance(pm::_PM.AbstractACRModel, nw::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_gs, bus_bs)
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)
    p    = get(var(pm, nw),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(var(pm, nw),    :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, nw),   :pg_bus, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, nw),   :qg_bus, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, nw),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, nw),   :qs, Dict()); _PM._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, nw),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, nw),  :qsw, Dict()); _PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(var(pm, nw),   :pt, Dict()); _PM._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(var(pm, nw),   :qt, Dict()); _PM._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    pd   = get(var(pm, nw),  :pd_bus, Dict()); _PM._check_var_keys(pd, bus_loads, "active power", "load")
    qd   = get(var(pm, nw),  :qd_bus, Dict()); _PM._check_var_keys(pd, bus_loads, "reactive power", "load")

    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    Gt = isempty(bus_gs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_gs))
    Bt = isempty(bus_bs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_bs))

    cstr_p = []
    cstr_q = []

    # pd/qd can be NLexpressions, so cannot be vectorized
    for c in conductor_ids(pm; nw=nw)
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

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = cstr_q

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


"""
Creates Ohms constraints

s_fr = v_fr.*conj(Y*(v_fr-v_to))
s_fr = (vr_fr+im*vi_fr).*(G-im*B)*([vr_fr-vr_to]-im*[vi_fr-vi_to])
s_fr = (vr_fr+im*vi_fr).*([G*vr_fr-G*vr_to-B*vi_fr+B*vi_to]-im*[G*vi_fr-G*vi_to+B*vr_fr-B*vr_to])
"""
function constraint_mc_ohms_yt_from(pm::_PM.AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, G, B, G_fr, B_fr, tr, ti, tm)
    p_fr  = var(pm, n, :p, f_idx)
    q_fr  = var(pm, n, :q, f_idx)
    vr_fr = var(pm, n, :vr, f_bus)
    vr_to = var(pm, n, :vr, t_bus)
    vi_fr = var(pm, n, :vi, f_bus)
    vi_to = var(pm, n, :vi, t_bus)

    cnds = conductor_ids(pm; nw=n)

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
function constraint_mc_ohms_yt_to(pm::_PM.AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    constraint_mc_ohms_yt_from(pm, n, t_bus, f_bus, t_idx, f_idx, g, b, g_to, b_to, tr, ti, tm)
end


""
function constraint_mc_load_setpoint_wye(pm::_PM.AbstractACRModel, nw::Int, id::Int, bus_id::Int, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    nph = 3

    # if constant power load
    if all(alpha.==0) && all(beta.==0)
        pd_bus = a
        qd_bus = b
    else
        crd = JuMP.@NLexpression(pm.model, [i in 1:nph],
            a[i]*vr[i]*(vr[i]^2+vi[i]^2)^(alpha[i]/2-1)
           +b[i]*vi[i]*(vr[i]^2+vi[i]^2)^(beta[i]/2 -1)
        )
        cid = JuMP.@NLexpression(pm.model, [i in 1:nph],
            a[i]*vi[i]*(vr[i]^2+vi[i]^2)^(alpha[i]/2-1)
           -b[i]*vr[i]*(vr[i]^2+vi[i]^2)^(beta[i]/2 -1)
        )

        pd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph],  vr[i]*crd[i]+vi[i]*cid[i])
        qd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], -vr[i]*cid[i]+vi[i]*crd[i])
    end

    var(pm, nw, :pd_bus)[id] = pd_bus
    var(pm, nw, :qd_bus)[id] = qd_bus

    if report
        sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        pd = JuMP.@NLexpression(pm.model, [i in 1:nph], a[i]*(vr[i]^2+vi[i]^2)^(alpha[i]/2) )
        qd = JuMP.@NLexpression(pm.model, [i in 1:nph], b[i]*(vr[i]^2+vi[i]^2)^(beta[i]/2)  )
        sol(pm, nw, :load, id)[:pd] = pd
        sol(pm, nw, :load, id)[:qd] = qd
    end
end


""
function constraint_mc_load_setpoint_delta(pm::_PM.AbstractACRModel, nw::Int, id::Int, bus_id::Int, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    nph = 3
    prev = Dict(i=>(i+nph-2)%nph+1 for i in 1:nph)
    next = Dict(i=>i%nph+1 for i in 1:nph)

    vrd = JuMP.@NLexpression(pm.model, [i in 1:nph], vr[i]-vr[next[i]])
    vid = JuMP.@NLexpression(pm.model, [i in 1:nph], vi[i]-vi[next[i]])

    crd = JuMP.@NLexpression(pm.model, [i in 1:nph],
        a[i]*vrd[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2-1)
       +b[i]*vid[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2 -1)
    )
    cid = JuMP.@NLexpression(pm.model, [i in 1:nph],
        a[i]*vid[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2-1)
       -b[i]*vrd[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2 -1)
    )

    crd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], crd[i]-crd[prev[i]])
    cid_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], cid[i]-cid[prev[i]])

    pd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph],  vr[i]*crd_bus[i]+vi[i]*cid_bus[i])
    qd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], -vr[i]*cid_bus[i]+vi[i]*crd_bus[i])

    var(pm, nw, :pd_bus)[id] = pd_bus
    var(pm, nw, :qd_bus)[id] = qd_bus

    if report
        sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        pd = JuMP.@NLexpression(pm.model, [i in 1:nph], a[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2) )
        qd = JuMP.@NLexpression(pm.model, [i in 1:nph], b[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2)  )
        sol(pm, nw, :load, id)[:pd] = pd
        sol(pm, nw, :load, id)[:qd] = qd
    end
end


"`vm[i] == vmref`"
function constraint_mc_voltage_magnitude_only(pm::_PM.AbstractACRModel, n::Int, i::Int, vmref)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    JuMP.@constraint(pm.model, vr.^2 + vi.^2  .== vmref.^2)
end


""
function constraint_mc_gen_setpoint_delta(pm::_PM.AbstractACRModel, nw::Int, id::Int, bus_id::Int, pmin::Vector, pmax::Vector, qmin::Vector, qmax::Vector; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    pg = var(pm, nw, :pg, id)
    qg = var(pm, nw, :qg, id)

    crg = []
    cig = []

    nph = 3
    prev = Dict(i=>(i+nph-2)%nph+1 for i in 1:nph)
    next = Dict(i=>i%nph+1 for i in 1:nph)

    vrg = JuMP.@NLexpression(pm.model, [i in 1:nph], vr[i]-vr[next[i]])
    vig = JuMP.@NLexpression(pm.model, [i in 1:nph], vi[i]-vi[next[i]])

    crg = JuMP.@NLexpression(pm.model, [i in 1:nph], (pg[i]*vrg[i]+qg[i]*vig[i])/(vrg[i]^2+vig[i]^2) )
    cig = JuMP.@NLexpression(pm.model, [i in 1:nph], (pg[i]*vig[i]-qg[i]*vrg[i])/(vrg[i]^2+vig[i]^2) )

    crg_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], crg[i]-crg[prev[i]])
    cig_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], cig[i]-cig[prev[i]])

    pg_bus = JuMP.@NLexpression(pm.model, [i in 1:nph],  vr[i]*crg_bus[i]+vi[i]*cig_bus[i])
    qg_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], -vr[i]*cig_bus[i]+vi[i]*crg_bus[i])

    var(pm, nw, :pg_bus)[id] = pg_bus
    var(pm, nw, :qg_bus)[id] = qg_bus

    if report
        sol(pm, nw, :gen, id)[:pg_bus] = pg_bus
        sol(pm, nw, :gen, id)[:qg_bus] = qg_bus
    end
end


"This function adds all constraints required to model a two-winding, wye-wye connected transformer."
function constraint_mc_transformer_power_yy(pm::_PM.AbstractACRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
    vr_fr = [var(pm, nw, :vr, f_bus)[p] for p in f_cnd]
    vr_to = [var(pm, nw, :vr, t_bus)[p] for p in t_cnd]
    vi_fr = [var(pm, nw, :vi, f_bus)[p] for p in f_cnd]
    vi_to = [var(pm, nw, :vi, t_bus)[p] for p in t_cnd]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[p] ? tm_set[p] : var(pm, nw, :tap, trans_id)[p] for p in conductor_ids(pm)]


    for p in conductor_ids(pm)
        if tm_fixed[p]
            JuMP.@constraint(pm.model, vr_fr[p] == pol*tm_scale*tm[p]*vr_to[p])
            JuMP.@constraint(pm.model, vi_fr[p] == pol*tm_scale*tm[p]*vi_to[p])
        else
            JuMP.@constraint(pm.model, vr_fr[p] == pol*tm_scale*tm[p]*vr_to[p])
            JuMP.@constraint(pm.model, vi_fr[p] == pol*tm_scale*tm[p]*vi_to[p])
        end
    end

    p_fr = var(pm, nw, :pt, f_idx)[f_cnd]
    p_to = var(pm, nw, :pt, t_idx)[t_cnd]
    q_fr = var(pm, nw, :qt, f_idx)[f_cnd]
    q_to = var(pm, nw, :qt, t_idx)[t_cnd]

    JuMP.@constraint(pm.model, p_fr + p_to .== 0)
    JuMP.@constraint(pm.model, q_fr + q_to .== 0)
end


"This function adds all constraints required to model a two-winding, delta-wye connected transformer."
function constraint_mc_transformer_power_dy(pm::_PM.AbstractACRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
    vr_p_fr = [var(pm, nw, :vr, f_bus)[p] for p in f_cnd]
    vr_p_to = [var(pm, nw, :vr, t_bus)[p] for p in t_cnd]
    vi_p_fr = [var(pm, nw, :vi, f_bus)[p] for p in f_cnd]
    vi_p_to = [var(pm, nw, :vi, t_bus)[p] for p in t_cnd]

    @assert length(f_cnd)==3 && length(t_cnd)==3

    M = _get_delta_transformation_matrix(3)

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[p] ? tm_set[p] : var(pm, nw, p, :tap, trans_id) for p in conductor_ids(pm)]

    # introduce auxialiary variable vd = Md*v_fr
    vrd = M*vr_p_fr
    vid = M*vi_p_fr

    JuMP.@constraint(pm.model, vrd .== (pol*tm_scale)*tm.*vr_p_to)
    JuMP.@constraint(pm.model, vid .== (pol*tm_scale)*tm.*vi_p_to)

    p_fr = [var(pm, nw, :pt, f_idx)[p] for p in f_cnd]
    p_to = [var(pm, nw, :pt, t_idx)[p] for p in t_cnd]
    q_fr = [var(pm, nw, :qt, f_idx)[p] for p in f_cnd]
    q_to = [var(pm, nw, :qt, t_idx)[p] for p in t_cnd]

    id_re = Array{Any,1}(undef, 3)
    id_im = Array{Any,1}(undef, 3)
    # s/v      = (p+jq)/|v|^2*conj(v)
    #          = (p+jq)/|v|*(cos(va)-j*sin(va))
    # Re(s/v)  = (p*cos(va)+q*sin(va))/|v|
    # -Im(s/v) = -(q*cos(va)-p*sin(va))/|v|
    for p in conductor_ids(pm)
        # id = conj(s_to/v_to)./tm
        id_re[p] = JuMP.@NLexpression(pm.model, (p_to[p]*vr_p_to[p]+q_to[p]*vi_p_to[p])/(tm_scale*tm[p]*pol)/(vr_p_to[p]^2+vi_p_to[p]^2))
        id_im[p] = JuMP.@NLexpression(pm.model, (p_to[p]*vi_p_to[p]-q_to[p]*vr_p_to[p])/(tm_scale*tm[p]*pol)/(vr_p_to[p]^2+vi_p_to[p]^2))
    end
    for (p,q) in zip([1:3...], [3,1,2])
        # s_fr  = v_fr*conj(i_fr)
        #       = v_fr*conj(id[q]-id[p])
        #       = v_fr*(id_re[q]-j*id_im[q]-id_re[p]+j*id_im[p])
        JuMP.@NLconstraint(pm.model, p_fr[p] ==
             vr_p_fr[p]*(id_re[q]-id_re[p])
            -vi_p_fr[p]*(-id_im[q]+id_im[p])
        )
        JuMP.@NLconstraint(pm.model, q_fr[p] ==
             vr_p_fr[p]*(-id_im[q]+id_im[p])
            +vi_p_fr[p]*(id_re[q]-id_re[p])
        )
    end
end
