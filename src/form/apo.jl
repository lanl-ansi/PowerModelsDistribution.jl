import LinearAlgebra: diag


"apo models ignore reactive power flows"
function variable_mc_gen_power_setpoint_imaginary(pm::_PM.AbstractActivePowerModel; kwargs...)
end


"apo models ignore reactive power flows"
function variable_mc_gen_power_setpoint_imaginary_on_off(pm::_PM.AbstractActivePowerModel; kwargs...)
end


"on/off constraint for generators"
function constraint_mc_gen_power_on_off(pm::_PM.AbstractActivePowerModel, n::Int, i::Int, pmin, pmax, qmin, qmax)
    pg = var(pm, n, :pg, i)
    z = var(pm, n, :z_gen, i)

    for c in conductor_ids(pm, n)
        if isfinite(pmax[c])
            JuMP.@constraint(pm.model, pg[c] .<= pmax[c].*z)
        end

        if isfinite(pmin[c])
            JuMP.@constraint(pm.model, pg[c] .>= pmin[c].*z)
        end
    end
end


"nothing to do"
function constraint_mc_voltage_magnitude_only(pm::_PM.AbstractActivePowerModel, n::Int, i::Int, vmref)
end


"apo models ignore reactive power flows"
function variable_mc_storage_power_imaginary(pm::_PM.AbstractActivePowerModel; kwargs...)
end


"apo models ignore reactive power flows"
function variable_mc_storage_power_imaginary_on_off(pm::_PM.AbstractActivePowerModel; kwargs...)
end


"apo models ignore reactive power flows"
function variable_mc_branch_power_imaginary(pm::_PM.AbstractActivePowerModel; kwargs...)
end


"apo models ignore reactive power flows"
function variable_mc_branch_flow_ne_reactive(pm::_PM.AbstractActivePowerModel; kwargs...)
end


"nothing to do, these models do not have complex voltage variables"
function variable_mc_bus_voltage(pm::_PM.AbstractNFAModel; nw=pm.cnw, kwargs...)
end

"nothing to do, these models do not have angle difference constraints"
function constraint_mc_voltage_angle_difference(pm::_PM.AbstractNFAModel, n::Int, f_idx, angmin, angmax)
end


"apo models ignore reactive power flows"
function variable_mc_transformer_power_imaginary(pm::_PM.AbstractActivePowerModel; nw::Int=pm.cnw, bounded=true)
end


"power balanace constraint with line shunts and transformers, active power only"
function constraint_mc_load_power_balance(pm::_PM.AbstractActivePowerModel, nw::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_gs, bus_bs)
    p    = get(var(pm, nw),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    pg   = get(var(pm, nw),   :pg_bus, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    ps   = get(var(pm, nw),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    psw  = get(var(pm, nw),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    pt   = get(var(pm, nw),   :pt, Dict()); _PM._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    pd   = get(var(pm, nw),   :pd_bus, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")

    cstr_p = []

    for c in conductor_ids(pm; nw=nw)
        cp = JuMP.@constraint(pm.model,
            sum(p[a][c] for a in bus_arcs)
            + sum(psw[a_sw][c] for a_sw in bus_arcs_sw)
            + sum(pt[a_trans][c] for a_trans in bus_arcs_trans)
            ==
            sum(pg[g][c] for g in bus_gens)
            - sum(ps[s][c] for s in bus_storage)
            - sum(pd[d][c] for d in bus_loads)
            - sum(diag(gs)[c] for gs in values(bus_gs))*1.0^2
        )
        push!(cstr_p, cp)
    end
    # omit reactive constraint
    cnds = conductor_ids(pm, nw)
    ncnds = length(cnds)

    con(pm, nw, :lam_kcl_r)[i] = isa(cstr_p, Array) ? cstr_p : [cstr_p]
    con(pm, nw, :lam_kcl_i)[i] = []

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = [NaN for i in 1:ncnds]
    end
end


######## Lossless Models ########

"Create variables for the active power flowing into all transformer windings"
function variable_mc_transformer_power_real(pm::_PM.AbstractAPLossLessModels; nw::Int=pm.cnw, bounded=true)
    ncnds = length(conductor_ids(pm))

    pt = var(pm, nw)[:pt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_pt_$((l,i,j))", start=0
        ) for (l,i,j) in ref(pm, nw, :arcs_from_trans)
    )

    if bounded
        for arc in ref(pm, nw, :arcs_from_trans)
            (t,i,j) = arc
            rate_a_fr, rate_a_to = _calc_transformer_power_ub_frto(ref(pm, nw, :transformer, t), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))

            set_lower_bound.(pt[(t,i,j)], -min.(rate_a_fr, rate_a_to))
            set_upper_bound.(pt[(t,i,j)],  min.(rate_a_fr, rate_a_to))
        end
    end

    for cnd in conductor_ids(pm)

        #TODO what does this dfo, and p does not seem to be defined!
        for (l,branch) in ref(pm, nw, :branch)
            if haskey(branch, "pf_start")
                f_idx = (l, branch["f_bus"], branch["t_bus"])
                JuMP.set_value(p[f_idx], branch["pf_start"])
            end
        end

    end

    # this explicit type erasure is necessary
    p_expr = Dict{Any,Any}( ((l,i,j), pt[(l,i,j)]) for (l,i,j) in ref(pm, nw, :arcs_from_trans) )
    p_expr = merge(p_expr, Dict( ((l,j,i), -1.0*pt[(l,i,j)]) for (l,i,j) in ref(pm, nw, :arcs_from_trans)))
    var(pm, nw)[:pt] = p_expr
end


""
function constraint_mc_network_power_balance(pm::_PM.AbstractAPLossLessModels, n::Int, i, comp_gen_ids, comp_pd, comp_qd, comp_gs, comp_bs, comp_branch_g, comp_branch_b)
    pg = var(pm, n, :pg)

    for c in conductor_ids(pm)
        JuMP.@constraint(pm.model, sum(pg[g][c] for g in comp_gen_ids) == sum(pd[c] for (i,pd) in values(comp_pd)) + sum(gs[c]*1.0^2 for (i,gs) in values(comp_gs)))
        # omit reactive constraint
    end
end


"Do nothing, this model is symmetric"
function constraint_mc_ohms_yt_to(pm::_PM.AbstractAPLossLessModels, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end

### Network Flow Approximation ###

"nothing to do, no voltage angle variables"
function constraint_mc_theta_ref(pm::_PM.AbstractNFAModel, n::Int, d::Int, va_ref)
end


"nothing to do, no voltage angle variables"
function constraint_mc_ohms_yt_from(pm::_PM.AbstractNFAModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


"nothing to do, this model is symmetric"
function constraint_mc_ohms_yt_to(pm::_PM.AbstractNFAModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end


"nothing to do, this model is symmetric"
function constraint_mc_transformer_power(pm::_PM.AbstractNFAModel, i::Int; nw::Int=pm.cnw)
end


"`-rate_a <= p[f_idx] <= rate_a`"
function constraint_mc_thermal_limit_from(pm::_PM.AbstractActivePowerModel, n::Int, f_idx, rate_a)
    cnds = conductor_ids(pm, n)
    ncnds = length(cnds)
    mu_sm_fr = []

    for c in 1:ncnds
        p_fr =var(pm, n, :p, f_idx)[c]
        if isa(p_fr, JuMP.VariableRef) && JuMP.has_lower_bound(p_fr)
           push!(mu_sm_fr,JuMP.LowerBoundRef(p_fr))
            JuMP.lower_bound(p_fr) < -rate_a[c] && set_lower_bound(p_fr, -rate_a[c])
            if JuMP.has_upper_bound(p_fr)
                JuMP.upper_bound(p_fr) > rate_a[c] && set_upper_bound(p_fr, rate_a[c])
            end
        else
           push!(mu_sm_fr, JuMP.@constraint(pm.model, p_fr <= rate_a[c]))
        end
    end

    if _IM.report_duals(pm)
        sol(pm, n, :branch, f_idx[1])[:mu_sm_fr] = mu_sm_fr
    end
end


""
function constraint_mc_thermal_limit_to(pm::_PM.AbstractActivePowerModel, n::Int, t_idx, rate_a)
    cnds = conductor_ids(pm, n)
    ncnds = length(cnds)
    mu_sm_to = []

    for c in 1:ncnds
        p_to =var(pm, n, :p, t_idx)[c]
        if isa(p_to, JuMP.VariableRef) && JuMP.has_lower_bound(p_to)
           push!(mu_sm_to, JuMP.LowerBoundRef(p_to))
            JuMP.lower_bound(p_to) < -rate_a[c] && set_lower_bound(p_to, -rate_a[c])
            if JuMP.has_upper_bound(p_to)
                JuMP.upper_bound(p_to) >  rate_a[c] && set_upper_bound(p_to,  rate_a[c])
            end
        else
           push!(mu_sm_to, JuMP.@constraint(pm.model, p_to <= rate_a[c]))
        end
    end

    if _IM.report_duals(pm)
        sol(pm, n, :branch, t_idx[1])[:mu_sm_to] = mu_sm_to
    end
end

""
function constraint_mc_current_limit(pm::_PM.AbstractActivePowerModel, n::Int, f_idx, c_rating_a)
    p_fr =var(pm, n, :p, f_idx)

    cnds = conductor_ids(pm, n)
    ncnds = length(cnds)

    for c in 1:ncnds
        JuMP.lower_bound(p_fr[c]) < -c_rating_a[c] && set_lower_bound(p_fr[c], -c_rating_a[c])
        JuMP.upper_bound(p_fr[c]) >  c_rating_a[c] && set_upper_bound(p_fr[c],  c_rating_a[c])
    end
end


""
function constraint_mc_thermal_limit_from_on_off(pm::_PM.AbstractActivePowerModel, n::Int, i, f_idx, rate_a)
    p_fr =var(pm, n, :p, f_idx)
    z =var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, p_fr .<=  rate_a.*z)
    JuMP.@constraint(pm.model, p_fr .>= -rate_a.*z)
end

""
function constraint_mc_thermal_limit_to_on_off(pm::_PM.AbstractActivePowerModel, n::Int, i, t_idx, rate_a)
    p_to =var(pm, n, :p, t_idx)
    z =var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, p_to .<=  rate_a.*z)
    JuMP.@constraint(pm.model, p_to .>= -rate_a.*z)
end


""
function constraint_mc_thermal_limit_from_ne(pm::_PM.AbstractActivePowerModel, n::Int, i, f_idx, rate_a)
    p_fr =var(pm, n, :p_ne, f_idx)
    z =var(pm, n, :branch_ne, i)

    JuMP.@constraint(pm.model, p_fr .<=  rate_a.*z)
    JuMP.@constraint(pm.model, p_fr .>= -rate_a.*z)
end


""
function constraint_mc_thermal_limit_to_ne(pm::_PM.AbstractActivePowerModel, n::Int, i, t_idx, rate_a)
    p_to =var(pm, n, :p_ne, t_idx)
    z =var(pm, n, :branch_ne, i)

    JuMP.@constraint(pm.model, p_to .<=  rate_a.*z)
    JuMP.@constraint(pm.model, p_to .>= -rate_a.*z)
end


""
function constraint_mc_switch_thermal_limit(pm::_PM.AbstractActivePowerModel, n::Int, f_idx, rating)
    psw =var(pm, n, :psw, f_idx)

    cnds = conductor_ids(pm, n)
    ncnds = length(cnds)

    for c in 1:ncnds
        JuMP.lower_bound(psw[c]) < -rating[c] && set_lower_bound(psw[c], -rating[c])
        JuMP.upper_bound(psw[c]) >  rating[c] && set_upper_bound(psw[c],  rating[c])
    end
end


""
function constraint_mc_storage_thermal_limit(pm::_PM.AbstractActivePowerModel, n::Int, i, rating)
    ps =var(pm, n, :ps, i)

    cnds = conductor_ids(pm, n)
    ncnds = length(cnds)

    for c in 1:ncnds
        JuMP.has_lower_bound(ps[c]) && JuMP.lower_bound(ps[c]) < -rating[c] && set_lower_bound(ps[c], -rating[c])
        JuMP.has_upper_bound(ps[c]) && JuMP.upper_bound(ps[c]) >  rating[c] && set_upper_bound(ps[c],  rating[c])
    end
end


""
function constraint_mc_storage_current_limit(pm::_PM.AbstractActivePowerModel, n::Int, i, bus, rating)
    ps =var(pm, n, :ps, i)

    cnds = conductor_ids(pm, n)
    ncnds = length(cnds)

    for c in 1:ncnds
        JuMP.lower_bound(ps[c]) < -rating[c] && set_lower_bound(ps[c], -rating[c])
        JuMP.upper_bound(ps[c]) >  rating[c] && set_upper_bound(ps[c],  rating[c])
    end
end


""
function constraint_mc_storage_losses(pm::_PM.AbstractActivePowerModel, n::Int, i, bus, conductors, r, x, p_loss, q_loss)
    ps = var(pm, n, :ps, i)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)

    JuMP.@constraint(pm.model,
        sum(ps[c] for c in conductors) + (sd - sc)
        ==
        p_loss + sum(r[c]*ps[c]^2 for c in conductors)
    )
end


""
function constraint_mc_storage_on_off(pm::_PM.AbstractActivePowerModel, n::Int, i, pmin, pmax, qmin, qmax, charge_ub, discharge_ub)

    z_storage =var(pm, n, :z_storage, i)
    ps =var(pm, n, :ps, i)

    JuMP.@constraint(pm.model, ps .<= pmax.*z_storage)
    JuMP.@constraint(pm.model, ps .>= pmin.*z_storage)

end


"Only support wye-connected generators."
function constraint_mc_gen_setpoint(pm::_PM.AbstractActivePowerModel, id::Int; nw::Int=pm.cnw, report::Bool=true)
    var(pm, nw, :pg_bus)[id] = var(pm, nw, :pg, id)
end


"Only support wye-connected, constant-power loads."
function constraint_mc_load_setpoint(pm::_PM.AbstractActivePowerModel, id::Int; nw::Int=pm.cnw, report::Bool=true)
    load = ref(pm, nw, :load, id)

    pd = load["pd"]

    var(pm, nw, :pd)[id] = pd
    var(pm, nw, :pd_bus)[id] = var(pm, nw, :pd, id)


    if report
        sol(pm, nw, :load, id)[:pd] = var(pm, nw, :pd, id)
        sol(pm, nw, :load, id)[:pd_bus] = var(pm, nw, :pd_bus, id)
    end
end
