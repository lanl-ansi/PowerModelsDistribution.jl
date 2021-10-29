import LinearAlgebra: diag


"apo models ignore reactive power flows"
function variable_mc_generator_power_imaginary(pm::AbstractUnbalancedActivePowerModel; kwargs...)
end


"apo models ignore reactive power flows"
function variable_mc_generator_power_imaginary_on_off(pm::AbstractUnbalancedActivePowerModel; kwargs...)
end


"on/off constraint for generators"
function constraint_mc_gen_power_on_off(pm::AbstractUnbalancedActivePowerModel, nw::Int, i::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real})
    z = var(pm, nw, :z_gen, i)
    pg = var(pm, nw, :pg_bus)[i] = var(pm, nw, :pg, i)

    for (idx, c) in enumerate(connections)
        if isfinite(pmax[idx])
            JuMP.@constraint(pm.model, pg[c] .<= pmax[idx].*z)
        end

        if isfinite(pmin[idx])
            JuMP.@constraint(pm.model, pg[c] .>= pmin[idx].*z)
        end
    end
end


"nothing to do"
function constraint_mc_voltage_magnitude_only(pm::AbstractUnbalancedActivePowerModel, nw::Int, i::Int, vm_ref::Vector{<:Real})
end


"apo models ignore reactive power flows"
function variable_mc_storage_power_imaginary(pm::AbstractUnbalancedActivePowerModel; kwargs...)
end


"apo models ignore reactive power flows"
function variable_mc_storage_power_imaginary_on_off(pm::AbstractUnbalancedActivePowerModel; kwargs...)
end


"apo models ignore reactive power flows"
function variable_mc_branch_power_imaginary(pm::AbstractUnbalancedActivePowerModel; kwargs...)
end


"apo models ignore reactive power flows"
function variable_mc_branch_flow_ne_reactive(pm::AbstractUnbalancedActivePowerModel; kwargs...)
end


"nothing to do, these models do not have complex voltage variables"
function variable_mc_bus_voltage(pm::AbstractUnbalancedNFAModel; nw::Int=nw_id_default, kwargs...)
end


"nothing to do, these models do not have complex voltage variables"
function variable_mc_capcontrol(pm::AbstractUnbalancedNFAModel; nw::Int=nw_id_default, kwargs...)
end


"nothing to do"
function constraint_mc_switch_state_closed(pm::AbstractUnbalancedNFAModel, nw::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int})
end


"nothing to do, these models do not have angle difference constraints"
function constraint_mc_voltage_angle_difference(pm::AbstractUnbalancedNFAModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})
end


"apo models ignore reactive power flows"
function variable_mc_transformer_power_imaginary(pm::AbstractUnbalancedActivePowerModel; nw::Int=nw_id_default, bounded=true)
end


@doc raw"""
    constraint_mc_power_balance(pm::AbstractUnbalancedActivePowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})

power balanace constraint with line shunts and transformers, active power only

```math
p_{br} + p_{tr} + p_{sw} == p_{g} - p_{s} - p_{d} - G
```
"""
function constraint_mc_power_balance(pm::AbstractUnbalancedActivePowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    p   = get(var(pm, nw),      :p,  Dict()); _check_var_keys(p,   bus_arcs, "active power", "branch")
    psw = get(var(pm, nw),    :psw,  Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    pt  = get(var(pm, nw),     :pt,  Dict()); _check_var_keys(pt,  bus_arcs_trans, "active power", "transformer")
    pg  = get(var(pm, nw),     :pg,  Dict()); _check_var_keys(pg,  bus_gens, "active power", "generator")
    ps  = get(var(pm, nw),     :ps,  Dict()); _check_var_keys(ps,  bus_storage, "active power", "storage")
    pd  = get(var(pm, nw), :pd_bus,  Dict()); _check_var_keys(pd,  bus_loads, "active power", "load")

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
              sum(p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(pt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
            ==
              sum(pg[g][t] for (g, conns) in bus_gens if t in conns)
            - sum(ps[s][t] for (s, conns) in bus_storage if t in conns)
            - sum(pd[d][t] for (d, conns) in bus_loads if t in conns)
            - diag(Gt)[idx]
        )
        push!(cstr_p, cp)
    end
    # omit reactive constraint

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = []

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = []
    end
end


"nothing to do in capc, no complex voltage variables in these models"
constraint_mc_power_balance_capc(pm::AbstractUnbalancedActivePowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}}) = constraint_mc_power_balance(pm, nw, i, terminals, grounded, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_shunts)


######## Lossless Models ########

""
function constraint_mc_network_power_balance(pm::AbstractUnbalancedAPLossLessModels, nw::Int, i, comp_gen_ids, comp_pd, comp_qd, comp_gs, comp_bs, comp_branch_g, comp_branch_b)
    pg = var(pm, nw, :pg)

    for (idx, t) in enumerate(ref(pm, nw, :bus, i, "terminals"))
        JuMP.@constraint(pm.model, sum(pg[g][t] for (g, conns) in comp_gen_ids if t in conns) == sum(pd[findfirst(isequal(t), conns)] for (i,conns,pd) in values(comp_pd) if t in conns) + sum(gs[findfirst(isequal(t), conns)]*1.0^2 for (i,conns,gs) in values(comp_gs) if t in conns))
        # omit reactive constraint
    end
end


"Do nothing, this model is symmetric"
function constraint_mc_ohms_yt_to(pm::AbstractUnbalancedAPLossLessModels, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})
end

### Network Flow Approximation ###

"nothing to do, no voltage angle variables"
function constraint_mc_theta_ref(pm::NFAUPowerModel, i::Int; nw::Int=nw_id_default)
end


"nothing to do, no voltage angle variables"
function constraint_mc_ohms_yt_from(pm::NFAUPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})
end


"nothing to do, this model is symmetric"
function constraint_mc_ohms_yt_to(pm::NFAUPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})
end


@doc raw"""
    constraint_mc_transformer_power(pm::NFAUPowerModel, i::Int; nw::Int=nw_id_default)

transformer active power only constraint pf=-pt

```math
p_f[fc] == -pt[tc]
```
"""
function constraint_mc_transformer_power(pm::NFAUPowerModel, i::Int; nw::Int=nw_id_default, fix_taps::Bool=false)
    transformer = ref(pm, nw, :transformer, i)

    pf = var(pm, nw, :pt, (i, transformer["f_bus"], transformer["t_bus"]))
    pt = var(pm, nw, :pt, (i, transformer["t_bus"], transformer["f_bus"]))

    for (fc, tc) in zip(transformer["f_connections"],transformer["t_connections"])
        JuMP.@constraint(pm.model, pf[fc] == -pt[tc])
    end
end


"`-rate_a <= p[f_idx] <= rate_a`"
function constraint_mc_thermal_limit_from(pm::AbstractUnbalancedActivePowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})
    mu_sm_fr = JuMP.ConstraintRef[]

    for (idx, c) in enumerate(f_connections)
        p_fr = var(pm, nw, :p, f_idx)[c]
        if isa(p_fr, JuMP.VariableRef) && JuMP.has_lower_bound(p_fr)
           push!(mu_sm_fr,JuMP.LowerBoundRef(p_fr))
            JuMP.lower_bound(p_fr) < -rate_a[idx] && set_lower_bound(p_fr, -rate_a[idx])
            if JuMP.has_upper_bound(p_fr)
                JuMP.upper_bound(p_fr) > rate_a[idx] && set_upper_bound(p_fr, rate_a[idx])
            end
        else
            if rate_a[idx] < Inf
                push!(mu_sm_fr, JuMP.@constraint(pm.model, p_fr <= rate_a[idx]))
            end
        end
    end

    con(pm, nw, :mu_sm_branch)[f_idx] = mu_sm_fr

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, f_idx[1])[:mu_sm_fr] = mu_sm_fr
    end
end


""
function constraint_mc_thermal_limit_to(pm::AbstractUnbalancedActivePowerModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})
    mu_sm_to = JuMP.ConstraintRef[]

    for (idx,c) in enumerate(t_connections)
        p_to = var(pm, nw, :p, t_idx)[c]
        if isa(p_to, JuMP.VariableRef) && JuMP.has_lower_bound(p_to)
           push!(mu_sm_to, JuMP.LowerBoundRef(p_to))
            JuMP.lower_bound(p_to) < -rate_a[idx] && set_lower_bound(p_to, -rate_a[idx])
            if JuMP.has_upper_bound(p_to)
                JuMP.upper_bound(p_to) >  rate_a[idx] && set_upper_bound(p_to,  rate_a[idx])
            end
        else
            if rate_a[idx] < Inf
                push!(mu_sm_to, JuMP.@constraint(pm.model, p_to <= rate_a[idx]))
            end
        end
    end

    con(pm, nw, :mu_sm_branch)[t_idx] = mu_sm_to

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, t_idx[1])[:mu_sm_to] = mu_sm_to
    end
end

@doc raw"""
    constraint_mc_switch_thermal_limit(pm::AbstractUnbalancedActivePowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{<:Real})::Nothing

Active power only switch thermal limit constraint

math```
-S_{max} \leq p_{fr} \leq S_{max}
```
"""
function constraint_mc_switch_thermal_limit(pm::AbstractUnbalancedActivePowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{<:Real})::Nothing
    psw = var(pm, nw, :psw, f_idx)

    mu_sm_fr = JuMP.ConstraintRef[]
    for (idx, c) in enumerate(f_connections)
        if rating[idx] < Inf
            JuMP.lower_bound(psw[c]) < -rating[idx] && set_lower_bound(psw[c], -rating[idx])
            JuMP.upper_bound(psw[c]) >  rating[idx] && set_upper_bound(psw[c],  rating[idx])
        end
    end

    con(pm, nw, :mu_sm_switch)[f_idx] = mu_sm_fr

    nothing
end


""
function constraint_mc_storage_thermal_limit(pm::AbstractUnbalancedActivePowerModel, nw::Int, i::Int, connections::Vector{Int}, rating::Vector{<:Real})
    ps =var(pm, nw, :ps, i)

    for (idx,c) in enumerate(connections)
        JuMP.has_lower_bound(ps[c]) && JuMP.lower_bound(ps[c]) < -rating[idx] && set_lower_bound(ps[c], -rating[idx])
        JuMP.has_upper_bound(ps[c]) && JuMP.upper_bound(ps[c]) >  rating[idx] && set_upper_bound(ps[c],  rating[idx])
    end
end


""
function constraint_mc_storage_current_limit(pm::AbstractUnbalancedActivePowerModel, nw::Int, i::Int, bus::Int, connections::Vector{Int}, rating::Vector{<:Real})
    ps = var(pm, nw, :ps, i)

    for (idx,c) in enumerate(connections)
        JuMP.lower_bound(ps[c]) < -rating[idx] && set_lower_bound(ps[c], -rating[idx])
        JuMP.upper_bound(ps[c]) >  rating[idx] && set_upper_bound(ps[c],  rating[idx])
    end
end


""
function constraint_mc_storage_losses(pm::AbstractUnbalancedActivePowerModel, nw::Int, i::Int, bus, connections::Vector{Int}, r::Vector{<:Real}, x::Vector{<:Real}, p_loss::Real, q_loss::Real)
    ps = var(pm, nw, :ps, i)
    sc = var(pm, nw, :sc, i)
    sd = var(pm, nw, :sd, i)

    JuMP.@constraint(pm.model,
        sum(ps[c] for c in connections) + (sd - sc)
        ==
        p_loss + sum(r[idx]*ps[c]^2 for (idx,c) in enumerate(connections))
    )
end


""
function constraint_mc_storage_on_off(pm::AbstractUnbalancedActivePowerModel, nw::Int, i::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}, charge_ub, discharge_ub)
    z_storage = var(pm, nw, :z_storage, i)
    ps = [var(pm, nw, :ps, i)[c] for c in connections]

    JuMP.@constraint(pm.model, ps .<= pmax.*z_storage)
    JuMP.@constraint(pm.model, ps .>= pmin.*z_storage)

end


"Only support wye-connected generators."
function constraint_mc_generator_power(pm::AbstractUnbalancedActivePowerModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
    var(pm, nw, :pg_bus)[id] = var(pm, nw, :pg, id)
end


"Only support wye-connected, constant-power loads."
function constraint_mc_load_power(pm::AbstractUnbalancedActivePowerModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
    load = ref(pm, nw, :load, id)
    connections = load["connections"]

    pd = load["pd"]

    var(pm, nw, :pd)[id] = JuMP.Containers.DenseAxisArray(pd, connections)
    var(pm, nw, :pd_bus)[id] = var(pm, nw, :pd, id)


    if report
        sol(pm, nw, :load, id)[:pd] = var(pm, nw, :pd, id)
        sol(pm, nw, :load, id)[:pd_bus] = var(pm, nw, :pd_bus, id)
    end
end


""
function constraint_storage_losses(pm::AbstractUnbalancedActivePowerModel, n::Int, i, bus, r, x, p_loss, q_loss; conductors=[1])
    ps = var(pm, n, :ps, i)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)

    JuMP.@constraint(pm.model,
        sum(ps[c] for c in conductors) + (sd - sc)
        ==
        p_loss + sum(r[c]*ps[c]^2 for c in conductors)
    )
end


"nothing to do, no voltage variables"
function constraint_mc_ampacity_from(pm::AbstractUnbalancedActivePowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{Float64})::Nothing
end


"nothing to do, no voltage variables"
function constraint_mc_ampacity_to(pm::AbstractUnbalancedActivePowerModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rating::Vector{Float64})::Nothing
end


"nothing to do, no voltage variables"
function constraint_mc_switch_ampacity(pm::AbstractUnbalancedActivePowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{Float64})::Nothing
end
