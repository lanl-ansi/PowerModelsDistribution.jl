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


######## Lossless Models ########
"Create variables for the active power flowing into all transformer windings"
function variable_mc_transformer_active_flow(pm::_PMs.AbstractAPLossLessModels; nw::Int=pm.cnw, bounded=true)
    for cnd in _PMs.conductor_ids(pm)
        pt = _PMs.var(pm, nw, cnd)[:pt] = JuMP.@variable(pm.model,
            [(l,i,j) in _PMs.ref(pm, nw, :arcs_from_trans)],
            base_name="$(nw)_$(cnd)_p_trans",
            start=0
        )
        if bounded
            for arc in _PMs.ref(pm, nw, :arcs_from_trans)
                tr_id = arc[1]
                flow_lb  = -_PMs.ref(pm, nw, :transformer, tr_id, "rate_a")[cnd]
                flow_ub  =  _PMs.ref(pm, nw, :transformer, tr_id, "rate_a")[cnd]
                JuMP.set_lower_bound(pt[arc], flow_lb)
                JuMP.set_upper_bound(pt[arc], flow_ub)
            end
        end

        for (l,branch) in _PMs.ref(pm, nw, :branch)
            if haskey(branch, "pf_start")
                f_idx = (l, branch["f_bus"], branch["t_bus"])
                JuMP.set_value(p[f_idx], branch["pf_start"])
            end
        end

        # this explicit type erasure is necessary
        p_expr = Dict{Any,Any}( ((l,i,j), pt[(l,i,j)]) for (l,i,j) in _PMs.ref(pm, nw, :arcs_from_trans) )
        p_expr = merge(p_expr, Dict( ((l,j,i), -1.0*pt[(l,i,j)]) for (l,i,j) in _PMs.ref(pm, nw, :arcs_from_trans)))
        _PMs.var(pm, nw, cnd)[:pt] = p_expr
    end
end


"Do nothing, this model is symmetric"
function constraint_mc_ohms_yt_to(pm::_PMs.AbstractAPLossLessModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end


"Do nothing, this model is symmetric"
function constraint_mc_ohms_yt_to_on_off(pm::_PMs.AbstractAPLossLessModels, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
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
    # omit reactive constraint
end


"`-b*(t[f_bus] - t[t_bus] + vad_min*(1-branch_z[i])) <= p[f_idx] <= -b*(t[f_bus] - t[t_bus] + vad_max*(1-branch_z[i]))`"
function constraint_mc_ohms_yt_from_on_off(pm::_PMs.AbstractDCPModel, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
    p_fr  = _PMs.var(pm, n, c,  :p, f_idx)
    va_fr = [_PMs.var(pm, n, d, :va, f_bus) for d in _PMs.conductor_ids(pm)]
    va_to = [_PMs.var(pm, n, d, :va, t_bus) for d in _PMs.conductor_ids(pm)]
    z = [_PMs.var(pm, n, d, :branch_z, i) for d in _PMs.conductor_ids(pm)]

    JuMP.@constraint(pm.model, p_fr <= sum(-b[c,d]*(va_fr[d] - va_to[d] + vad_max[d]*(1-z[d])) for d in _PMs.conductor_ids(pm)) )
    JuMP.@constraint(pm.model, p_fr >= sum(-b[c,d]*(va_fr[d] - va_to[d] + vad_min[d]*(1-z[d])) for d in _PMs.conductor_ids(pm)) )
end


### Network Flow Approximation ###
"nothing to do, no voltage angle variables"
function constraint_mc_theta_ref(pm::_PMs.AbstractNFAModel, n::Int, c::Int, d)
end


"nothing to do, no voltage angle variables"
function constraint_mc_ohms_yt_from(pm::_PMs.AbstractNFAModel, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


"nothing to do, this model is symmetric"
function constraint_mc_ohms_yt_to(pm::_PMs.AbstractNFAModel, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end


"nothing to do, no voltage variables"
function constraint_mc_transformer_voltage(pm::_PMs.AbstractNFAModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, tm::_PMs.MultiConductorVector, Tv_fr, Tv_im, Cv_to)
end


"nothing to do, this model is symmetric"
function constraint_mc_transformer_flow(pm::_PMs.AbstractNFAModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, tm::_PMs.MultiConductorVector, Ti_fr, Ti_im, Cv_to)
end


"power balance constraint with line shunts and transformers for load shed problem, DCP formulation"
function constraint_mc_power_balance_shed(pm::_PMs.AbstractDCPModel, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
    p = _PMs.var(pm, nw, c, :p)
    pg = _PMs.var(pm, nw, c, :pg)
    p_dc = _PMs.var(pm, nw, c, :p_dc)
    p_trans = _PMs.var(pm, nw, c, :pt)
    z_demand = _PMs.var(pm, nw, :z_demand)
    z_shunt = _PMs.var(pm, nw, :z_shunt)

    _PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) + sum(p_trans[a_trans] for a_trans in bus_arcs_trans) == sum(pg[g] for g in bus_gens) - sum(pd*z_demand[n] for (n,pd) in bus_pd) - sum(gs*1.0^2*z_shunt[n] for (n,gs) in bus_gs))
end


"on/off bus voltage constraint for DCP formulation, nothing to do"
function constraint_mc_bus_voltage_on_off(pm::_PMs.AbstractDCPModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, kwargs...)
end
