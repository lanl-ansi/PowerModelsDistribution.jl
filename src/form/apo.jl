### generic features that apply to all active-power-only (apo) approximations

"apo models ignore reactive power flows"
function variable_mc_generation_reactive(pm::_PMs.AbstractActivePowerModel; kwargs...)
end


"apo models ignore reactive power flows"
function variable_mc_reactive_generation_on_off(pm::_PMs.AbstractActivePowerModel; kwargs...)
end


"apo models ignore reactive power flows"
function variable_mc_storage_reactive(pm::_PMs.AbstractActivePowerModel; kwargs...)
end


"apo models ignore reactive power flows"
function variable_mc_on_off_storage_reactive(pm::_PMs.AbstractActivePowerModel; kwargs...)
end


"apo models ignore reactive power flows"
function variable_mc_branch_flow_reactive(pm::_PMs.AbstractActivePowerModel; kwargs...)
end


"apo models ignore reactive power flows"
function variable_mc_branch_flow_ne_reactive(pm::_PMs.AbstractActivePowerModel; kwargs...)
end


# "do nothing, apo models do not have reactive variables"
# function constraint_mc_gen_setpoint_reactive(pm::_PMs.AbstractActivePowerModel, n::Int, c::Int, i, qg)
# end


"nothing to do, these models do not have complex voltage variables"
function variable_mc_voltage(pm::_PMs.AbstractNFAModel; nw=pm.cnw, kwargs...)
end

"nothing to do, these models do not have angle difference constraints"
function constraint_mc_voltage_angle_difference(pm::_PMs.AbstractNFAModel, n::Int, f_idx, angmin, angmax)
end




"apo models ignore reactive power flows"
function variable_mc_transformer_flow_reactive(pm::_PMs.AbstractActivePowerModel; nw::Int=pm.cnw, bounded=true)
end


"power balanace constraint with line shunts and transformers, active power only"
function constraint_mc_power_balance(pm::_PMs.AbstractActivePowerModel, nw::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    p    = get(_PMs.var(pm, nw),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    pg   = get(_PMs.var(pm, nw),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    ps   = get(_PMs.var(pm, nw),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    psw  = get(_PMs.var(pm, nw),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    pt   = get(_PMs.var(pm, nw),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")

    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        .==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd.values for pd in values(bus_pd))
        - sum(gs.values for gs in values(bus_gs))*1.0^2
    )
    # omit reactive constraint
    cnds = _PMs.conductor_ids(pm, nw)
    ncnds = length(cnds)

    if _PMs.report_duals(pm)
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_i] = [NaN for i in 1:ncnds]
    end
end


######## Lossless Models ########
"Create variables for the active power flowing into all transformer windings"
function variable_mc_transformer_flow_active(pm::_PMs.AbstractAPLossLessModels; nw::Int=pm.cnw, bounded=true)
    ncnds = length(_PMs.conductor_ids(pm))

    pt = _PMs.var(pm, nw)[:pt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_pt_$((l,i,j))", start=0
        ) for (l,i,j) in _PMs.ref(pm, nw, :arcs_from_trans)
    )

    for cnd in _PMs.conductor_ids(pm)

        if bounded
            for arc in _PMs.ref(pm, nw, :arcs_from_trans)
                tr_id = arc[1]
                flow_lb  = -_PMs.ref(pm, nw, :transformer, tr_id, "rate_a")[cnd]
                flow_ub  =  _PMs.ref(pm, nw, :transformer, tr_id, "rate_a")[cnd]
                JuMP.set_lower_bound(pt[arc][cnd], flow_lb)
                JuMP.set_upper_bound(pt[arc][cnd], flow_ub)
            end
        end

        #TODO what does this dfo, and p does not seem to be defined!
        for (l,branch) in _PMs.ref(pm, nw, :branch)
            if haskey(branch, "pf_start")
                f_idx = (l, branch["f_bus"], branch["t_bus"])
                JuMP.set_value(p[f_idx], branch["pf_start"])
            end
        end

    end

    # this explicit type erasure is necessary
    p_expr = Dict{Any,Any}( ((l,i,j), pt[(l,i,j)]) for (l,i,j) in _PMs.ref(pm, nw, :arcs_from_trans) )
    p_expr = merge(p_expr, Dict( ((l,j,i), -1.0*pt[(l,i,j)]) for (l,i,j) in _PMs.ref(pm, nw, :arcs_from_trans)))
    _PMs.var(pm, nw)[:pt] = p_expr
end


"Do nothing, this model is symmetric"
function constraint_mc_ohms_yt_to(pm::_PMs.AbstractAPLossLessModels, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end

### Network Flow Approximation ###
"nothing to do, no voltage angle variables"
function constraint_mc_theta_ref(pm::_PMs.AbstractNFAModel, n::Int, d)
end


"nothing to do, no voltage angle variables"
function constraint_mc_ohms_yt_from(pm::_PMs.AbstractNFAModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


"nothing to do, this model is symmetric"
function constraint_mc_ohms_yt_to(pm::_PMs.AbstractNFAModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end


"nothing to do, this model is symmetric"
function constraint_mc_trans(pm::_PMs.AbstractNFAModel, i::Int; nw::Int=pm.cnw)
end

## From PowerModels

"`-rate_a <= p[f_idx] <= rate_a`"
function constraint_mc_thermal_limit_from(pm::_PMs.AbstractActivePowerModel, n::Int, f_idx, rate_a)
    cnds = _PMs.conductor_ids(pm, n)
    ncnds = length(cnds)
    mu_sm_fr = []

    for c in 1:ncnds
        p_fr =_PMs.var(pm, n, :p, f_idx)[c]
        if isa(p_fr, JuMP.VariableRef) && JuMP.has_lower_bound(p_fr)
           push!(mu_sm_fr,JuMP.LowerBoundRef(p_fr))
            JuMP.lower_bound(p_fr) < -rate_a[c] && JuMP.set_lower_bound(p_fr, -rate_a[c])
            if JuMP.has_upper_bound(p_fr)
                JuMP.upper_bound(p_fr) > rate_a[c] && JuMP.set_upper_bound(p_fr, rate_a[c])
            end
        else
           push!(mu_sm_fr, JuMP.@constraint(pm.model, p_fr <= rate_a[c]))
        end
    end

    if _PMs.report_duals(pm)
        _PMs.sol(pm, n, :branch, f_idx[1])[:mu_sm_fr] = mu_sm_fr
    end
end

""
function constraint_mc_thermal_limit_to(pm::_PMs.AbstractActivePowerModel, n::Int, t_idx, rate_a)
    cnds = _PMs.conductor_ids(pm, n)
    ncnds = length(cnds)
    mu_sm_to = []

    for c in 1:ncnds
        p_to =_PMs.var(pm, n, :p, t_idx)[c]
        if isa(p_to, JuMP.VariableRef) && JuMP.has_lower_bound(p_to)
           push!(mu_sm_to, JuMP.LowerBoundRef(p_to))
            JuMP.lower_bound(p_to) < -rate_a[c] && JuMP.set_lower_bound(p_to, -rate_a[c])
            if JuMP.has_upper_bound(p_to)
                JuMP.upper_bound(p_to) >  rate_a[c] && JuMP.set_upper_bound(p_to,  rate_a[c])
            end
        else
           push!(mu_sm_to, JuMP.@constraint(pm.model, p_to <= rate_a[c]))
        end
    end

    if _PMs.report_duals(pm)
        _PMs.sol(pm, n, :branch, t_idx[1])[:mu_sm_to] = mu_sm_to
    end
end

""
function constraint_mc_current_limit(pm::_PMs.AbstractActivePowerModel, n::Int, f_idx, c_rating_a)
    p_fr =_PMs.var(pm, n, :p, f_idx)

    cnds = _PMs.conductor_ids(pm, n)
    ncnds = length(cnds)

    for c in 1:ncnds
        JuMP.lower_bound(p_fr[c]) < -c_rating_a[c] && JuMP.set_lower_bound(p_fr[c], -c_rating_a[c])
        JuMP.upper_bound(p_fr[c]) >  c_rating_a[c] && JuMP.set_upper_bound(p_fr[c],  c_rating_a[c])
    end
end


""
function constraint_mc_thermal_limit_from_on_off(pm::_PMs.AbstractActivePowerModel, n::Int, i, f_idx, rate_a)
    p_fr =_PMs.var(pm, n, :p, f_idx)
    z =_PMs.var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, p_fr .<=  rate_a.*z)
    JuMP.@constraint(pm.model, p_fr .>= -rate_a.*z)
end

""
function constraint_mc_thermal_limit_to_on_off(pm::_PMs.AbstractActivePowerModel, n::Int, i, t_idx, rate_a)
    p_to =_PMs.var(pm, n, :p, t_idx)
    z =_PMs.var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, p_to .<=  rate_a.*z)
    JuMP.@constraint(pm.model, p_to .>= -rate_a.*z)
end

""
function constraint_mc_thermal_limit_from_ne(pm::_PMs.AbstractActivePowerModel, n::Int, i, f_idx, rate_a)
    p_fr =_PMs.var(pm, n, :p_ne, f_idx)
    z =_PMs.var(pm, n, :branch_ne, i)

    JuMP.@constraint(pm.model, p_fr .<=  rate_a.*z)
    JuMP.@constraint(pm.model, p_fr .>= -rate_a.*z)
end

""
function constraint_mc_thermal_limit_to_ne(pm::_PMs.AbstractActivePowerModel, n::Int, i, t_idx, rate_a)
    p_to =_PMs.var(pm, n, :p_ne, t_idx)
    z =_PMs.var(pm, n, :branch_ne, i)

    JuMP.@constraint(pm.model, p_to .<=  rate_a.*z)
    JuMP.@constraint(pm.model, p_to .>= -rate_a.*z)
end


""
function constraint_mc_switch_thermal_limit(pm::_PMs.AbstractActivePowerModel, n::Int, f_idx, rating)
    psw =_PMs.var(pm, n, :psw, f_idx)

    cnds = _PMs.conductor_ids(pm, n)
    ncnds = length(cnds)

    for c in 1:ncnds
        JuMP.lower_bound(psw[c]) < -rating[c] && JuMP.set_lower_bound(psw[c], -rating[c])
        JuMP.upper_bound(psw[c]) >  rating[c] && JuMP.set_upper_bound(psw[c],  rating[c])
    end
end


""
function constraint_mc_storage_thermal_limit(pm::_PMs.AbstractActivePowerModel, n::Int, i, rating)
    ps =_PMs.var(pm, n, :ps, i)

    cnds = _PMs.conductor_ids(pm, n)
    ncnds = length(cnds)

    for c in 1:ncnds
        JuMP.lower_bound(ps[c]) < -rating[c] && JuMP.set_lower_bound(ps[c], -rating[c])
        JuMP.upper_bound(ps[c]) >  rating[c] && JuMP.set_upper_bound(ps[c],  rating[c])
    end
end

""
function constraint_mc_storage_current_limit(pm::_PMs.AbstractActivePowerModel, n::Int, i, bus, rating)
    ps =_PMs.var(pm, n, :ps, i)

    cnds = _PMs.conductor_ids(pm, n)
    ncnds = length(cnds)

    for c in 1:ncnds
        JuMP.lower_bound(ps[c]) < -rating[c] && JuMP.set_lower_bound(ps[c], -rating[c])
        JuMP.upper_bound(ps[c]) >  rating[c] && JuMP.set_upper_bound(ps[c],  rating[c])
    end
end

""
function constraint_mc_storage_loss(pm::_PMs.AbstractActivePowerModel, n::Int, i, bus, conductors, r, x, p_loss, q_loss)
    ps = _PMs.var(pm, n, :ps, i)
    sc = _PMs.var(pm, n, :sc, i)
    sd = _PMs.var(pm, n, :sd, i)

    JuMP.@constraint(pm.model,
        sum(ps[c] for c in conductors) + (sd - sc)
        ==
        p_loss + sum(r[c]*ps[c]^2 for c in conductors)
    )
end

function constraint_mc_storage_on_off(pm::_PMs.AbstractActivePowerModel, n::Int, i, pmin, pmax, qmin, qmax, charge_ub, discharge_ub)

    z_storage =_PMs.var(pm, n, :z_storage, i)
    ps =_PMs.var(pm, n, :ps, i)

    JuMP.@constraint(pm.model, ps .<= pmax.*z_storage)
    JuMP.@constraint(pm.model, ps .>= pmin.*z_storage)

end

#
# ""
# function add_setpoint_switch_flow!(sol, pm::_PMs.AbstractActivePowerModel)
#     add_setpoint!(sol, pm, "switch", "psw", :psw, var_key = (idx,item) -> (idx, item["f_bus"], item["t_bus"]))
#     add_setpoint_fixed!(sol, pm, "switch", "qsw")
# end
