import LinearAlgebra: diag


"apo models ignore reactive power flows"
function variable_mc_generator_power_imaginary(pm::_PM.AbstractActivePowerModel; kwargs...)
end


"apo models ignore reactive power flows"
function variable_mc_generator_power_imaginary_on_off(pm::_PM.AbstractActivePowerModel; kwargs...)
end


"on/off constraint for generators"
function constraint_mc_gen_power_on_off(pm::_PM.AbstractActivePowerModel, nw::Int, i::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real})
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
function constraint_mc_voltage_magnitude_only(pm::_PM.AbstractActivePowerModel, nw::Int, i::Int, vm_ref::Vector{<:Real})
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
function variable_mc_bus_voltage(pm::_PM.AbstractNFAModel; nw::Int=nw_id_default, kwargs...)
end


"nothing to do"
function constraint_mc_switch_state_closed(pm::_PM.AbstractNFAModel, nw::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int})
end


"nothing to do, these models do not have angle difference constraints"
function constraint_mc_voltage_angle_difference(pm::_PM.AbstractNFAModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})
end


"apo models ignore reactive power flows"
function variable_mc_transformer_power_imaginary(pm::_PM.AbstractActivePowerModel; nw::Int=nw_id_default, bounded=true)
end


"power balanace constraint with line shunts and transformers, active power only"
function constraint_mc_power_balance(pm::_PM.AbstractActivePowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    p    = get(var(pm, nw),    :p, Dict())#; _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    pg   = get(var(pm, nw),   :pg_bus, Dict())#; _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    ps   = get(var(pm, nw),   :ps, Dict())#; _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    psw  = get(var(pm, nw),  :psw, Dict())#; _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    pt   = get(var(pm, nw),   :pt, Dict())#; _PM._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    pd   = get(var(pm, nw),   :pd_bus, Dict())#; _PM._check_var_keys(pg, bus_gens, "active power", "generator")

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


######## Lossless Models ########

"Create variables for the active power flowing into all transformer windings"
function variable_mc_transformer_power_real(pm::_PM.AbstractAPLossLessModels; nw::Int=nw_id_default, bounded::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_transformer) for ((l,i,j), connections) in entry)
    pt = var(pm, nw)[:pt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_pt_$((l,i,j))", start=0
        ) for (l,i,j) in ref(pm, nw, :arcs_from_trans)
    )

    if bounded
        for arc in ref(pm, nw, :arcs_from_trans)
            (t,i,j) = arc
            rate_a_fr, rate_a_to = _calc_transformer_power_ub_frto(ref(pm, nw, :transformer, t), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))

            for (idx, c) in enumerate(connections[arc])
                set_lower_bound(pt[(t,i,j)][c], -min.(rate_a_fr, rate_a_to)[idx])
                set_upper_bound(pt[(t,i,j)][c],  min.(rate_a_fr, rate_a_to)[idx])
            end
        end
    end

    for (l,transformer) in ref(pm, nw, :transformer)
        if haskey(transformer, "pf_start")
            f_idx = (l, transformer["f_bus"], transformer["t_bus"])
            for (idx, c) in enumerate(connections[f_idx])
                JuMP.set_value(pt[f_idx][c], transformer["pf_start"][idx])
            end
        end
        if haskey(transformer, "pt_start")
            t_idx = (l, transformer["t_bus"], transformer["f_bus"])
            for (idx, c) in enumerate(connections[t_idx])
                JuMP.set_start_value(pt[t_idx][c], transformer["pt_start"][idx])
            end
        end
    end

    # this explicit type erasure is necessary
    p_expr = Dict{Any,Any}( ((l,i,j), pt[(l,i,j)]) for (l,i,j) in ref(pm, nw, :arcs_from_trans) )
    p_expr = merge(p_expr, Dict( ((l,j,i), -1.0*pt[(l,i,j)]) for (l,i,j) in ref(pm, nw, :arcs_from_trans)))
    var(pm, nw)[:pt] = p_expr
end


""
function constraint_mc_network_power_balance(pm::_PM.AbstractAPLossLessModels, nw::Int, i, comp_gen_ids, comp_pd, comp_qd, comp_gs, comp_bs, comp_branch_g, comp_branch_b)
    pg = var(pm, nw, :pg)

    for (idx, t) in enumerate(ref(pm, nw, :bus, i, "terminals"))
        JuMP.@constraint(pm.model, sum(pg[g][t] for (g, conns) in comp_gen_ids if t in conns) == sum(pd[findfirst(isequal(t), conns)] for (i,conns,pd) in values(comp_pd) if t in conns) + sum(gs[findfirst(isequal(t), conns)]*1.0^2 for (i,conns,gs) in values(comp_gs) if t in conns))
        # omit reactive constraint
    end
end


"Do nothing, this model is symmetric"
function constraint_mc_ohms_yt_to(pm::_PM.AbstractAPLossLessModels, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})
end

### Network Flow Approximation ###

"nothing to do, no voltage angle variables"
function constraint_mc_theta_ref(pm::_PM.NFAPowerModel, i::Int; nw::Int=nw_id_default)
end


"nothing to do, no voltage angle variables"
function constraint_mc_ohms_yt_from(pm::_PM.NFAPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})
end


"nothing to do, this model is symmetric"
function constraint_mc_ohms_yt_to(pm::_PM.NFAPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})
end


"nothing to do, this model is symmetric"
function constraint_mc_transformer_power(pm::_PM.NFAPowerModel, i::Int; nw::Int=nw_id_default)
end


"`-rate_a <= p[f_idx] <= rate_a`"
function constraint_mc_thermal_limit_from(pm::_PM.AbstractActivePowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})
    mu_sm_fr = []

    for (idx, c) in enumerate(f_connections)
        p_fr = var(pm, nw, :p, f_idx)[c]
        if isa(p_fr, JuMP.VariableRef) && JuMP.has_lower_bound(p_fr)
           push!(mu_sm_fr,JuMP.LowerBoundRef(p_fr))
            JuMP.lower_bound(p_fr) < -rate_a[idx] && set_lower_bound(p_fr, -rate_a[idx])
            if JuMP.has_upper_bound(p_fr)
                JuMP.upper_bound(p_fr) > rate_a[idx] && set_upper_bound(p_fr, rate_a[idx])
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
function constraint_mc_thermal_limit_to(pm::_PM.AbstractActivePowerModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})
    mu_sm_to = []

    for (idx,c) in enumerate(t_connections)
        p_to = var(pm, nw, :p, t_idx)[c]
        if isa(p_to, JuMP.VariableRef) && JuMP.has_lower_bound(p_to)
           push!(mu_sm_to, JuMP.LowerBoundRef(p_to))
            JuMP.lower_bound(p_to) < -rate_a[idx] && set_lower_bound(p_to, -rate_a[idx])
            if JuMP.has_upper_bound(p_to)
                JuMP.upper_bound(p_to) >  rate_a[idx] && set_upper_bound(p_to,  rate_a[idx])
            end
        else
           push!(mu_sm_to, JuMP.@constraint(pm.model, p_to <= rate_a[idx]))
        end
    end

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, t_idx[1])[:mu_sm_to] = mu_sm_to
    end
end

""
function constraint_mc_current_limit(pm::_PM.AbstractActivePowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, c_rating_a::Vector{<:Real})
    p_fr = var(pm, nw, :p, f_idx)

    for (idx,c) in enumerate(f_connections)
        JuMP.lower_bound(p_fr[c]) < -c_rating_a[idx] && set_lower_bound(p_fr[c], -c_rating_a[idx])
        JuMP.upper_bound(p_fr[c]) >  c_rating_a[idx] && set_upper_bound(p_fr[c],  c_rating_a[idx])
    end
end


""
function constraint_mc_thermal_limit_from_on_off(pm::_PM.AbstractActivePowerModel, nw::Int, i::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})
    p_fr = [var(pm, nw, :p, f_idx)[c] for c in f_connections]
    z = var(pm, nw, :z_branch, i)

    JuMP.@constraint(pm.model, p_fr .<=  rate_a.*z)
    JuMP.@constraint(pm.model, p_fr .>= -rate_a.*z)
end

""
function constraint_mc_thermal_limit_to_on_off(pm::_PM.AbstractActivePowerModel, nw::Int, i::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})
    p_to = [var(pm, nw, :p, t_idx)[c] for c in t_connections]
    z = var(pm, nw, :z_branch, i)

    JuMP.@constraint(pm.model, p_to .<=  rate_a.*z)
    JuMP.@constraint(pm.model, p_to .>= -rate_a.*z)
end


""
function constraint_mc_thermal_limit_from_ne(pm::_PM.AbstractActivePowerModel, nw::Int, i::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})
    p_fr = [var(pm, nw, :p_ne, f_idx)[c] for c in f_connections]
    z =var(pm, nw, :branch_ne, i)

    JuMP.@constraint(pm.model, p_fr .<=  rate_a.*z)
    JuMP.@constraint(pm.model, p_fr .>= -rate_a.*z)
end


""
function constraint_mc_thermal_limit_to_ne(pm::_PM.AbstractActivePowerModel, nw::Int, i::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})
    p_to = [var(pm, nw, :p_ne, t_idx)[c] for c in t_connections]
    z =var(pm, nw, :branch_ne, i)

    JuMP.@constraint(pm.model, p_to .<=  rate_a.*z)
    JuMP.@constraint(pm.model, p_to .>= -rate_a.*z)
end


""
function constraint_mc_switch_thermal_limit(pm::_PM.AbstractActivePowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{<:Real})
    psw = var(pm, nw, :psw, f_idx)

    for (idx, c) in enumerate(f_connections)
        JuMP.lower_bound(psw[c]) < -rating[idx] && set_lower_bound(psw[c], -rating[idx])
        JuMP.upper_bound(psw[c]) >  rating[idx] && set_upper_bound(psw[c],  rating[idx])
    end
end


""
function constraint_mc_storage_thermal_limit(pm::_PM.AbstractActivePowerModel, nw::Int, i::Int, connections::Vector{Int}, rating::Vector{<:Real})
    ps =var(pm, nw, :ps, i)

    for (idx,c) in enumerate(connections)
        JuMP.has_lower_bound(ps[c]) && JuMP.lower_bound(ps[c]) < -rating[idx] && set_lower_bound(ps[c], -rating[idx])
        JuMP.has_upper_bound(ps[c]) && JuMP.upper_bound(ps[c]) >  rating[idx] && set_upper_bound(ps[c],  rating[idx])
    end
end


""
function constraint_mc_storage_current_limit(pm::_PM.AbstractActivePowerModel, nw::Int, i::Int, bus::Int, connections::Vector{Int}, rating::Vector{<:Real})
    ps = var(pm, nw, :ps, i)

    for (idx,c) in enumerate(connections)
        JuMP.lower_bound(ps[c]) < -rating[idx] && set_lower_bound(ps[c], -rating[idx])
        JuMP.upper_bound(ps[c]) >  rating[idx] && set_upper_bound(ps[c],  rating[idx])
    end
end


""
function constraint_mc_storage_losses(pm::_PM.AbstractActivePowerModel, nw::Int, i::Int, bus, connections::Vector{Int}, r::Vector{<:Real}, x::Vector{<:Real}, p_loss::Real, q_loss::Real)
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
function constraint_mc_storage_on_off(pm::_PM.AbstractActivePowerModel, nw::Int, i::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}, charge_ub, discharge_ub)
    z_storage = var(pm, nw, :z_storage, i)
    ps = [var(pm, nw, :ps, i)[c] for c in connections]

    JuMP.@constraint(pm.model, ps .<= pmax.*z_storage)
    JuMP.@constraint(pm.model, ps .>= pmin.*z_storage)

end


"Only support wye-connected generators."
function constraint_mc_generator_power(pm::_PM.AbstractActivePowerModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
    var(pm, nw, :pg_bus)[id] = var(pm, nw, :pg, id)
end


"Only support wye-connected, constant-power loads."
function constraint_mc_load_power(pm::_PM.AbstractActivePowerModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
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
