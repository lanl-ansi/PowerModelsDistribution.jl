""
function variable_mc_bus_voltage(pm::_PM.AbstractDCPModel; nw=pm.cnw, kwargs...)
    variable_mc_bus_voltage_angle(pm; nw=nw, kwargs...)
end


######## AbstractDCPForm Models (has va but assumes vm is 1.0) ########

"nothing to do, these models do not have complex voltage constraints"
function constraint_mc_model_voltage(pm::_PM.AbstractDCPModel, n::Int, c::Int)
end


""
function variable_mc_bus_voltage_on_off(pm::_PM.AbstractDCPModel; kwargs...)
    variable_mc_bus_voltage_angle(pm; kwargs...)
end


### DC Power Flow Approximation ###

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] == -b*(t[f_bus] - t[t_bus])
```
"""
function constraint_mc_ohms_yt_from(pm::_PM.AbstractDCPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr  = var(pm, n,  :p, f_idx)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    for c in conductor_ids(pm, n)
        JuMP.@constraint(pm.model, p_fr[c] == -sum(b[c,d]*(va_fr[c] - va_to[d]) for d in conductor_ids(pm)))
    end
end


"power balance constraint with line shunts and transformers for load shed problem, DCP formulation"
function constraint_mc_shed_power_balance(pm::_PM.AbstractDCPModel, nw::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    p        = get(var(pm, nw),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    pg       = get(var(pm, nw),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    ps       = get(var(pm, nw),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    psw      = get(var(pm, nw),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    pt       = get(var(pm, nw),   :pt, Dict()); _PM._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    z_demand = var(pm, nw, :z_demand)
    z_shunt  = var(pm, nw, :z_shunt)

    cp = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        .==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd.*z_demand[n] for (n,pd) in bus_pd)
        - sum(diag(gs)*1.0^2 .*z_shunt[n] for (n,gs) in bus_gs)
    )

    con(pm, nw, :lam_kcl_r)[i] = isa(cp, JuMP.ConstraintRef) ? [cp] : cp
    con(pm, nw, :lam_kcl_i)[i] = []

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cp
    end
end


"on/off bus voltage constraint for DCP formulation, nothing to do"
function constraint_mc_bus_voltage_on_off(pm::_PM.AbstractDCPModel; nw::Int=pm.cnw, kwargs...)
end


""
function variable_mc_branch_power_real(pm::_PM.AbstractAPLossLessModels; nw::Int=pm.cnw, bounded::Bool = true, report::Bool = true)
    cnds = conductor_ids(pm)
    ncnds = length(cnds)

    p = Dict((l,i,j) => JuMP.@variable(pm.model,
        [c in 1:ncnds], base_name="$(nw)_($l,$i,$j)_p",
        start = comp_start_value(ref(pm, nw, :branch, l), "p_start", c, 0.0)
    ) for (l,i,j) in ref(pm, nw, :arcs_from))

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_from)
            smax = _calc_branch_power_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
            if !ismissing(smax)
                set_upper_bound.(p[(l,i,j)],  smax)
                set_lower_bound.(p[(l,i,j)], -smax)
            end
        end
    end

    for (l,branch) in ref(pm, nw, :branch)
        if haskey(branch, "pf_start")
            f_idx = (l, branch["f_bus"], branch["t_bus"])
            JuMP.set_start_value(p[f_idx], branch["pf_start"])
        end
    end

    # this explicit type erasure is necessary
    p_expr = Dict{Any,Any}( ((l,i,j), p[(l,i,j)]) for (l,i,j) in ref(pm, nw, :arcs_from) )
    p_expr = merge(p_expr, Dict( ((l,j,i), -1.0*p[(l,i,j)]) for (l,i,j) in ref(pm, nw, :arcs_from)))
    var(pm, nw)[:p] = p_expr
end
