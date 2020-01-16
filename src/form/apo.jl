### generic features that apply to all active-power-only (apo) approximations

"do nothing, no reactive power in this model"
function variable_mc_transformer_reactive_flow(pm::_PMs.AbstractActivePowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded=true)
end


"power balanace constraint with line shunts and transformers, active power only"
function constraint_mc_power_balance(pm::_PMs.AbstractActivePowerModel, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    p    = get(_PMs.var(pm, nw, c),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    pg   = get(_PMs.var(pm, nw, c),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    ps   = get(_PMs.var(pm, nw, c),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    psw  = get(_PMs.var(pm, nw, c),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    pt   = get(_PMs.var(pm, nw, c),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")

    _PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        - sum(gs for gs in values(bus_gs))*1.0^2
    )
    # omit reactive constraint
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


"nothing to do, this model is symmetric"
function constraint_mc_trans(pm::_PMs.AbstractNFAModel, i::Int; nw::Int=pm.cnw)
end

