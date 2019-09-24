### simple active power only approximations (e.g. DC Power Flow)


######## AbstractDCPForm Models (has va but assumes vm is 1.0) ########
"nothing to do, these models do not have complex voltage constraints"
function constraint_mc_model_voltage(pm::_PMs.AbstractDCPModel, n::Int, c::Int)
end


""
function variable_mc_bus_voltage_on_off(pm::_PMs.AbstractDCPModel; kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_voltage_angle(pm; cnd=c, kwargs...)
    end
end


### DC Power Flow Approximation ###
"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] == -b*(t[f_bus] - t[t_bus])
```
"""
function constraint_mc_ohms_yt_from(pm::_PMs.AbstractDCPModel, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr  = _PMs.var(pm, n, c,  :p, f_idx)
    va_fr = [_PMs.var(pm, n, d, :va, f_bus) for d in _PMs.conductor_ids(pm)]
    va_to = [_PMs.var(pm, n, d, :va, t_bus) for d in _PMs.conductor_ids(pm)]

    JuMP.@constraint(pm.model, p_fr == -sum(b[c,d]*(va_fr[c] - va_to[d]) for d in _PMs.conductor_ids(pm)))
end


"power balance constraint with line shunts and transformers for load shed problem, DCP formulation"
function constraint_mc_power_balance_shed(pm::_PMs.AbstractDCPModel, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    p        = get(_PMs.var(pm, nw, c),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    pg       = get(_PMs.var(pm, nw, c),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    ps       = get(_PMs.var(pm, nw, c),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    psw      = get(_PMs.var(pm, nw, c),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    pt       = get(_PMs.var(pm, nw, c),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    z_demand = _PMs.var(pm, nw, :z_demand)
    z_shunt  = _PMs.var(pm, nw, :z_shunt)

    _PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd*z_demand[n] for (n,pd) in bus_pd)
        - sum(gs*1.0^2*z_shunt[n] for (n,gs) in bus_gs)
    )
end


"on/off bus voltage constraint for DCP formulation, nothing to do"
function constraint_mc_bus_voltage_on_off(pm::_PMs.AbstractDCPModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, kwargs...)
end
