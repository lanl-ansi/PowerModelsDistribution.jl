# # BUS

# # Constraints

""
function constraint_mc_power_balance(pm::AbstractMultiConductorACRModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)
    p    = get(var(pm, nw), :p_bus,   Dict()); _check_var_keys(p,   bus_arcs,       "active power",   "branch")
    q    = get(var(pm, nw), :q_bus,   Dict()); _check_var_keys(q,   bus_arcs,       "reactive power", "branch")
    pg   = get(var(pm, nw), :pg_bus,  Dict()); _check_var_keys(pg,  bus_gens,       "active power",   "generator")
    qg   = get(var(pm, nw), :qg_bus,  Dict()); _check_var_keys(qg,  bus_gens,       "reactive power", "generator")
    # ps   = get(var(pm, nw), :ps_bus,  Dict()); _check_var_keys(ps,  bus_storage,    "active power",   "storage")
    # qs   = get(var(pm, nw), :qs_bus,  Dict()); _check_var_keys(qs,  bus_storage,    "reactive power", "storage")
    # psw  = get(var(pm, nw), :psw_bus, Dict()); _check_var_keys(psw, bus_arcs_sw,    "active power",   "switch")
    # qsw  = get(var(pm, nw), :qsw_bus, Dict()); _check_var_keys(qsw, bus_arcs_sw,    "reactive power", "switch")
    pt   = get(var(pm, nw), :pt_bus,  Dict()); _check_var_keys(pt,  bus_arcs_trans, "active power",   "transformer")
    qt   = get(var(pm, nw), :qt_bus,  Dict()); _check_var_keys(qt,  bus_arcs_trans, "reactive power", "transformer")
    pd   = get(var(pm, nw), :pd_bus,  Dict()); _check_var_keys(pd,  bus_loads,      "active power",   "load")
    qd   = get(var(pm, nw), :qd_bus,  Dict()); _check_var_keys(pd,  bus_loads,      "reactive power", "load")

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    # pd/qd can be NLexpressions, so cannot be vectorized
    for (idx, t) in ungrounded_terminals
        # cp = @smart_constraint(pm.model, [p, q, pg, qg, ps, qs, psw, qsw, pt, qt, pd, qd, vr, vi],
        cp = @smart_constraint(pm.model, [p, q, pg, qg, pt, qt, pd, qd, vr, vi],
              sum(  p[arc][t] for (arc, conns) in bus_arcs if t in conns)
            # + sum(psw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
            + sum( pt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
            ==
              sum(pg[gen][t] for (gen, conns) in bus_gens if t in conns)
            # - sum(ps[strg][t] for (strg, conns) in bus_storage if t in conns)
            - sum(pd[load][t] for (load, conns) in bus_loads if t in conns)
            + ( -vr[t] * sum(Gt[idx,jdx]*vr[u]-Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals)
                -vi[t] * sum(Gt[idx,jdx]*vi[u]+Bt[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals)
            )
        )
        push!(cstr_p, cp)

        # cq = @smart_constraint(pm.model, [p, q, pg, qg, ps, qs, psw, qsw, pt, qt, pd, qd, vr, vi],
        cq = @smart_constraint(pm.model, [p, q, pg, qg, pt, qt, pd, qd, vr, vi],
              sum(  q[arc][t] for (arc, conns) in bus_arcs if t in conns)
            # + sum(qsw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
            + sum( qt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
            ==
              sum(qg[gen][t] for (gen, conns) in bus_gens if t in conns)
            - sum(qd[load][t] for (load, conns) in bus_loads if t in conns)
            # - sum(qs[strg][t] for (strg, conns) in bus_storage if t in conns)
            + ( vr[t] * sum(Gt[idx,jdx]*vi[u]+Bt[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals)
               -vi[t] * sum(Gt[idx,jdx]*vr[u]-Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals)
            )
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


# # GENERATOR

# # GENERATOR - Variables

# ""
# function variable_mc_generator_current_real(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(i => _infer_int_dim_unit(gen, false) for (i,gen) in ref(pm, nw, :gen))
#     crg = var(pm, nw)[:crg] = Dict(i => JuMP.@variable(pm.model,
#             [c in 1:int_dim[i]], base_name="$(nw)_crg_$(i)",
#             start = comp_start_value(ref(pm, nw, :gen, i), "crg_start", c, 0.0)
#         ) for i in ids(pm, nw, :gen)
#     )
#     # if bounded
#     #     for (i, g) in ref(pm, nw, :gen) if haskey(g, "c_rating")
#     #         cmax = c["c_rating"]
#     #         set_lower_bound.(crg[i], -cmax)
#     #         set_upper_bound.(crg[i],  cmax)
#     #     end
#     # end

#     report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :crg, ids(pm, nw, :gen), crg)
# end


# ""
# function variable_mc_generator_current_imaginary(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(i => _infer_int_dim_unit(gen, false) for (i,gen) in ref(pm, nw, :gen))
#     cig = var(pm, nw)[:cig] = Dict(i => JuMP.@variable(pm.model,
#             [c in 1:int_dim[i]], base_name="$(nw)_cig_$(i)",
#             start = comp_start_value(ref(pm, nw, :gen, i), "cig_start", c, 0.0)
#         ) for i in ids(pm, nw, :gen)
#     )
#     # if bounded
#     #     for (i, g) in ref(pm, nw, :gen) if haskey(g, "c_rating")
#     #         cmax = c["c_rating"]
#     #         set_lower_bound.(cig[i], -cmax)
#     #         set_upper_bound.(cig[i],  cmax)
#     #     end
#     # end

#     report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :cig, ids(pm, nw, :gen), cig)
# end


""
function variable_mc_generator_power(pm::AbstractMultiConductorACRModel; nw::Int=nw_id_default, kwargs...)
    variable_mc_generator_power_real(pm; nw=nw, kwargs...)
    variable_mc_generator_power_imaginary(pm; nw=nw, kwargs...)
    var(pm, nw)[:pg_bus] = Dict{Int, Any}()
    var(pm, nw)[:qg_bus] = Dict{Int, Any}()
end


# ""
# function variable_mc_generator_power_real(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(i => _infer_int_dim_unit(gen, false) for (i,gen) in ref(pm, :gen))
#     pg = var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
#             [c in 1:int_dim[i]], base_name="$(nw)_pg_$(i)",
#             start = comp_start_value(ref(pm, nw, :gen, i), "pg_start", c, 0.0)
#         ) for i in ids(pm, nw, :gen)
#     )

#     if bounded
#         for (i,gen) in ref(pm, nw, :gen)
#             set_lower_bound.(pg[i], gen["pmin"])
#             set_upper_bound.(pg[i], gen["pmax"])
#         end
#     end

#     report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :pg, ids(pm, nw, :gen), pg)
# end


# ""
# function variable_mc_generator_power_imaginary(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(i => _infer_int_dim_unit(gen, false) for (i,gen) in ref(pm, :gen))
#     qg = var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
#             [c in 1:int_dim[i]], base_name="$(nw)_qg_$(i)",
#             start = comp_start_value(ref(pm, nw, :gen, i), "qg_start", c, 0.0)
#         ) for i in ids(pm, nw, :gen)
#     )

#     if bounded
#         for (i,gen) in ref(pm, nw, :gen)
#             set_lower_bound.(qg[i], gen["qmin"])
#             set_upper_bound.(qg[i], gen["qmax"])
#         end
#     end

#     report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :qg, ids(pm, nw, :gen), qg)
# end


# # GENERATOR - Variables - Non-linear

# ""
# function variable_mc_generator_power(pm::AbstractNLMultiConductorIVRModel; kwargs...)
#     #TODO docs
# end

# # GENERATOR - Variables - Quadratic

# ""
# function variable_mc_generator_power(pm::AbstractQuadraticMultiConductorIVRModel; kwargs...)
#     variable_mc_generator_power_real(pm; kwargs...)
#     variable_mc_generator_power_imaginary(pm; kwargs...)
# end


# # GENERATOR - Constraints


# function constraint_mc_generator_current(pm::AbstractMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)
#     generator = ref(pm, nw, :gen, id)

#     nphases = _infer_int_dim_unit(generator, false)
#     if get(generator, "configuration", WYE) == WYE || nphases==1
#         constraint_mc_generator_current_wye(pm, nw, id, generator["connections"]; report=report, bounded=bounded)
#     else
#         constraint_mc_generator_current_delta(pm, nw, id, generator["connections"]; report=report, bounded=bounded)
#     end
# end


# ""
# function constraint_mc_generator_current_wye(pm::AbstractMultiConductorIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
#     crg = var(pm, nw, :crg, id)
#     cig = var(pm, nw, :cig, id)
#     var(pm, nw, :crg_bus)[id] = _merge_bus_flows(pm, [crg..., -sum(crg)], connections)
#     var(pm, nw, :cig_bus)[id] = _merge_bus_flows(pm, [cig..., -sum(cig)], connections)
# end


# ""
# function constraint_mc_generator_current_delta(pm::AbstractMultiConductorIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
#     crg = var(pm, nw, :crg, id)
#     cig = var(pm, nw, :cig, id)
#     Md = _get_delta_transformation_matrix(length(connections))
#     var(pm, nw, :crg_bus)[id] = _merge_bus_flows(pm, Md'*crg, connections)
#     var(pm, nw, :cig_bus)[id] = _merge_bus_flows(pm, Md'*cig, connections)
# end


# # GENERATOR - Constraints - Non-linear

""
function constraint_mc_generator_power_wye(pm::AbstractMultiConductorACRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    pg = var(pm, nw, :pg, id)
    qg = var(pm, nw, :qg, id)

    phases = connections[1:end-1]
    n      = connections[end]
    P      = length(phases)
    
    if iszero(vr[n]) && iszero(vi[n])
        pg_bus_unmerged = [pg..., 0.0]
        qg_bus_unmerged = [qg..., 0.0]
    else
        crg = [JuMP.@NLexpression(pm.model, 
            (pg[idx]*(vr[p]-vr[n])+qg[idx]*(vi[p]-vi[n]))/((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)
        ) for (idx,p) in enumerate(phases)]
        cig = [JuMP.@NLexpression(pm.model, 
            (-pg[idx]*(vi[p]-vi[n])+qg[idx]*(vr[p]-vr[n]))/((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)
        ) for (idx,p) in enumerate(phases)]
        pg_bus_unmerged = [
            [JuMP.@NLexpression(pm.model, vr[p]*crg[idx]+vi[p]*cig[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, vr[n]*sum(-crg[idx] for idx in 1:P)+vi[n]*sum(-cig[idx] for idx in 1:P))
        ]
        qg_bus_unmerged = [
            [JuMP.@NLexpression(pm.model, -vr[p]*cig[idx]+vi[p]*crg[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, -vr[n]*sum(-cig[idx] for idx in 1:P)+vi[n]*sum(-crg[idx] for idx in 1:P))
        ]
    end
    var(pm, nw, :pg_bus)[id] = pg_bus = _merge_bus_flows(pm, pg_bus_unmerged, connections)
    var(pm, nw, :qg_bus)[id] = qg_bus = _merge_bus_flows(pm, qg_bus_unmerged, connections)

    if report
        sol(pm, nw, :gen, id)[:pg_bus] = pg_bus
        sol(pm, nw, :gen, id)[:qg_bus] = qg_bus
    end
end


# ""
# function constraint_mc_generator_power_delta(pm::AbstractNLMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
#     vr = var(pm, nw, :vr, bus_id)
#     vi = var(pm, nw, :vi, bus_id)
#     crg = var(pm, nw, :crg, id)
#     cig = var(pm, nw, :cig, id)

#     nph = length(pmin)

#     prev = Dict(c=>connections[(idx+nph-2)%nph+1] for (idx,c) in enumerate(connections))
#     next = Dict(c=>connections[idx%nph+1] for (idx,c) in enumerate(connections))

#     vrg = Dict()
#     vig = Dict()
#     for c in connections
#         vrg[c] = JuMP.@NLexpression(pm.model, vr[c]-vr[next[c]])
#         vig[c] = JuMP.@NLexpression(pm.model, vi[c]-vi[next[c]])
#     end

#     pg = Vector{JuMP.NonlinearExpression}([])
#     qg = Vector{JuMP.NonlinearExpression}([])
#     for (idx,c) in enumerate(connections)
#         push!(pg, JuMP.@NLexpression(pm.model,  vrg[c]*crg[idx]+vig[c]*cig[idx]))
#         push!(qg, JuMP.@NLexpression(pm.model, -vrg[c]*cig[idx]+vig[c]*crg[idx]))
#     end

#     if bounded
#         JuMP.@NLconstraint(pm.model, [i in 1:nph], pmin[i] <= pg[i])
#         JuMP.@NLconstraint(pm.model, [i in 1:nph], pmax[i] >= pg[i])
#         JuMP.@NLconstraint(pm.model, [i in 1:nph], qmin[i] <= qg[i])
#         JuMP.@NLconstraint(pm.model, [i in 1:nph], qmax[i] >= qg[i])
#     end

#     var(pm, nw, :pg)[id] = JuMP.Containers.DenseAxisArray(pg, connections)
#     var(pm, nw, :qg)[id] = JuMP.Containers.DenseAxisArray(qg, connections)

#     if report
#         sol(pm, nw, :gen, id)[:pg] = pg
#         sol(pm, nw, :gen, id)[:qg] = qg
#     end
# end


# # GENERATOR - Constraints - Quadratic

# ""
# function constraint_mc_generator_power_wye(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
#     vr = var(pm, nw, :vr, bus_id)
#     vi = var(pm, nw, :vi, bus_id)
#     crg = var(pm, nw, :crg, id)
#     cig = var(pm, nw, :cig, id)

#     phases = connections[1:end-1]
#     n      = connections[end]

#     vr_pn = [vr[p]-vr[n] for p in phases]
#     vi_pn = [vi[p]-vi[n] for p in phases]

#     pg = var(pm, nw, :pg, id)
#     qg = var(pm, nw, :qg, id)

#     JuMP.@constraint(pm.model, pg .==  vr_pn.*crg .+ vi_pn.*cig)
#     JuMP.@constraint(pm.model, qg .== -vr_pn.*cig .+ vi_pn.*crg)
# end


# ""
# function constraint_mc_generator_power_delta(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
#     vr = var(pm, nw, :vr, bus_id)
#     vi = var(pm, nw, :vi, bus_id)
#     crg = var(pm, nw, :crg, id)
#     cig = var(pm, nw, :cig, id)

#     nph = length(pmin)

#     prev = Dict(c=>connections[(idx+nph-2)%nph+1] for (idx,c) in enumerate(connections))
#     next = Dict(c=>connections[idx%nph+1] for (idx,c) in enumerate(connections))

#     vrg = [vr[c]-vr[next[c]] for p in connections]
#     vig = [vi[c]-vi[next[c]] for p in connections]

#     pg = var(pm, nw, :pg, id)
#     qg = var(pm, nw, :qg, id)

#     JuMP.@constraint(pm.model, pg .==  vrg.*crg .+ vig.*cig)
#     JuMP.@constraint(pm.model, qg .== -vrg.*cig .+ vig.*crg)
# end


# # LOAD

# # LOAD - Variables

# ""
# function variable_mc_load_current(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
#     var(pm, nw)[:crd] = Dict{Int, Any}()
#     var(pm, nw)[:cid] = Dict{Int, Any}()
#     var(pm, nw)[:crd_bus] = Dict{Int, Any}()
#     var(pm, nw)[:cid_bus] = Dict{Int, Any}()
# end


""
function variable_mc_load_power(pm::AbstractMultiConductorACRModel; nw=nw_id_default, bounded::Bool=true, report::Bool=true)
    var(pm, nw)[:pd] = Dict{Int, Any}()
    var(pm, nw)[:qd] = Dict{Int, Any}()
    var(pm, nw)[:pd_bus] = Dict{Int, Any}()
    var(pm, nw)[:qd_bus] = Dict{Int, Any}()
end


# ""
# function variable_mc_load_current(pm::AbstractQuadraticMultiConductorIVRModel; nw::Int=nw_id_default, kwargs...)
#     load_ids_exponential = [id for (id,load) in ref(pm, nw, :load) if load["model"]==EXPONENTIAL]
#     @assert isempty(load_ids_exponential) "Exponential loads cannot be represented quadratically."

#     variable_mc_load_current_real(pm; nw=nw, kwargs...)
#     variable_mc_load_current_imaginary(pm; nw=nw, kwargs...)

#     var(pm, nw)[:crd_bus] = Dict{Int,Any}()
#     var(pm, nw)[:cid_bus] = Dict{Int,Any}()
# end


# ""
# function variable_mc_load_current_real(pm::AbstractQuadraticMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(i => _infer_int_dim_unit(load, false) for (i,load) in ref(pm, :load))

#     crd = var(pm, nw)[:crd] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
#             [c in 1:int_dim[i]], base_name="$(nw)_crd_$(i)",
#             start = comp_start_value(ref(pm, nw, :load, i), "crd_start", c, 0.0)
#         ) for i in ids(pm, nw, :load)
#     )

#     # if bounded
#     #     for (i,load) in ref(pm, nw, :gen)
#     #         if haskey(gen, "pmin")
#     #             set_lower_bound.(pg[i], gen["pmin"])
#     #         end
#     #         if haskey(gen, "pmax")
#     #             set_upper_bound.(pg[i], gen["pmax"])
#     #         end
#     #     end
#     # end

#     report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :crd, ids(pm, nw, :load), crd)
# end


# ""
# function variable_mc_load_current_imaginary(pm::AbstractQuadraticMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(i => _infer_int_dim_unit(load, false) for (i,load) in ref(pm, :load))
#     load_ids_current = [id for (id,load) in ref(pm, nw, :load) if load["model"]==CURRENT]
    
#     cid = var(pm, nw)[:cid] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
#             [c in 1:int_dim[i]], base_name="$(nw)_cid_$(i)",
#             start = comp_start_value(ref(pm, nw, :load, i), "cid_start", c, 0.0)
#         ) for i in ids(pm, nw, :load)
#     )

#     # if bounded
#     #     for (i,gen) in ref(pm, nw, :gen)
#     #         if haskey(gen, "qmin")
#     #             set_lower_bound.(qg[i], gen["qmin"])
#     #         end
#     #         if haskey(gen, "qmax")
#     #             set_upper_bound.(qg[i], gen["qmax"])
#     #         end
#     #     end
#     # end

#     report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :cid, ids(pm, nw, :load), cid)
# end


# # LOAD - Variables - Quadratic

# ""
# function variable_mc_load_power(pm::AbstractQuadraticMultiConductorIVRModel; nw::Int=nw_id_default, kwargs...)
#     load_ids_exponential = [id for (id,load) in ref(pm, nw, :load) if load["model"]==EXPONENTIAL]
#     @assert isempty(load_ids_exponential) "Exponential loads cannot be represented quadratically."

#     variable_mc_load_power_real(pm; nw=nw, kwargs...)
#     variable_mc_load_power_imaginary(pm; nw=nw, kwargs...)
# end


# ""
# function variable_mc_load_power_real(pm::AbstractQuadraticMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(i => _infer_int_dim_unit(load, false) for (i,load) in ref(pm, :load))
#     load_ids_current = [id for (id,load) in ref(pm, nw, :load) if load["model"]==CURRENT]
    
#     pd = var(pm, nw)[:pd] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
#             [c in 1:int_dim[i]], base_name="$(nw)_pd_$(i)",
#             start = comp_start_value(ref(pm, nw, :load, i), "pd_start", c, 0.0)
#         ) for i in load_ids_current
#     )

#     report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :pd, load_ids_current, pd)
# end


# ""
# function variable_mc_load_power_imaginary(pm::AbstractQuadraticMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(i => _infer_int_dim_unit(load, false) for (i,load) in ref(pm, :load))
#     load_ids_current = [id for (id,load) in ref(pm, nw, :load) if load["model"]==CURRENT]
    
#     qd = var(pm, nw)[:qd] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
#             [c in 1:int_dim[i]], base_name="$(nw)_qd_$(i)",
#             start = comp_start_value(ref(pm, nw, :load, i), "qd_start", c, 0.0)
#         ) for i in load_ids_current
#     )

#     report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :qd, load_ids_current, qd)
# end


# # LOAD - Constraints

function constraint_mc_load_power(pm::MultiConductorModels, id::Int; nw::Int=nw_id_default, report::Bool=true)
    load = ref(pm, nw, :load, id)
    bus = ref(pm, nw,:bus, load["load_bus"])

    configuration = load["configuration"]

    a, alpha, b, beta = _load_expmodel_params(load, bus)

    if configuration==WYE || length(a)==1
        constraint_mc_load_power_wye(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    else
        constraint_mc_load_power_delta(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    end
end


function constraint_mc_load_power_wye(pm::AbstractMultiConductorACRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    phases = connections[1:end-1]
    n      = connections[end]
    P      = length(phases)

    vr_pn = [vr[p]-vr[n] for p in phases]
    vi_pn = [vi[p]-vi[n] for p in phases]

    crd = JuMP.@NLexpression(pm.model, [idx in 1:P], a[idx]*vr_pn[idx]*(vr_pn[idx]^2+vi_pn[idx]^2)^(alpha[idx]/2-1)+b[idx]*vi_pn[idx]*(vr_pn[idx]^2+vi_pn[idx]^2)^(beta[idx]/2 -1))
    cid = JuMP.@NLexpression(pm.model, [idx in 1:P], a[idx]*vi_pn[idx]*(vr_pn[idx]^2+vi_pn[idx]^2)^(alpha[idx]/2-1)-b[idx]*vr_pn[idx]*(vr_pn[idx]^2+vi_pn[idx]^2)^(beta[idx]/2 -1))


    # if constant power load
    if all(alpha.==0) && all(beta.==0)
        pd = a
        qd = b
    else
        pd = JuMP.@NLexpression(pm.model, [idx in 1:P],  vr_pn[idx]*crd[idx]+vi_pn[idx]*cid[idx])
        qd = JuMP.@NLexpression(pm.model, [idx in 1:P], -vr_pn[idx]*cid[idx]+vi_pn[idx]*crd[idx])
    end

    if iszero(vr[n]) && iszero(vi[n])
        pd_bus_unmerged = [pd..., 0.0]
        qd_bus_unmerged = [qd..., 0.0]
    else
        pd_bus_unmerged = [
            [JuMP.@NLexpression(pm.model, vr[p]*crd[idx]+vi[p]*cid[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, vr[n]*sum(-crd[idx] for idx in 1:P)+vi[n]*sum(-cid[idx] for idx in 1:P))
        ]
        qd_bus_unmerged = [
            [JuMP.@NLexpression(pm.model, -vr[p]*cid[idx]+vi[p]*crd[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, -vr[n]*sum(-cid[idx] for idx in 1:P)+vi[n]*sum(-crd[idx] for idx in 1:P))
        ]
    end
    var(pm, nw, :pd_bus)[id] = pd_bus = _merge_bus_flows(pm, pd_bus_unmerged, connections)
    var(pm, nw, :qd_bus)[id] = qd_bus = _merge_bus_flows(pm, qd_bus_unmerged, connections)

    if report
        sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        sol(pm, nw, :load, id)[:pd] = pd
        sol(pm, nw, :load, id)[:qd] = qd
    end
end


""
function constraint_mc_load_power_delta(pm::AbstractUnbalancedACRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    nph = length(a)
    @assert nph == 3 "only phases == 3 delta loads are currently supported"
    next(i)  = i%nph+1 # 1=>2, 2=>3, 3=>1
    prev(i)  = (i-2+nph)%nph+1  # 1=>3, 2=>1, 3=>2

    vrd = [vr[t]-vr[connections[next(idx)]] for (idx,t) in enumerate(connections)]
    vid = [vi[t]-vi[connections[next(idx)]] for (idx,t) in enumerate(connections)]

    crd = JuMP.@NLexpression(pm.model, [idx in 1:nph], a[idx]*vrd[idx]*(vrd[idx]^2+vid[idx]^2)^(alpha[idx]/2-1)+b[idx]*vid[idx]*(vrd[idx]^2+vid[idx]^2)^(beta[idx]/2 -1))
    cid = JuMP.@NLexpression(pm.model, [idx in 1:nph], a[idx]*vid[idx]*(vrd[idx]^2+vid[idx]^2)^(alpha[idx]/2-1)-b[idx]*vrd[idx]*(vrd[idx]^2+vid[idx]^2)^(beta[idx]/2 -1))

    crd_bus_unmerged = JuMP.@NLexpression(pm.model, [idx in 1:nph], crd[idx]-crd[prev(idx)])
    cid_bus_unmerged = JuMP.@NLexpression(pm.model, [idx in 1:nph], cid[idx]-cid[prev(idx)])

    pd_bus_unmerged = [JuMP.@NLexpression(pm.model,  vr[p]*crd_bus_unmerged[idx]+vi[p]*cid_bus_unmerged[idx]) for (idx,p) in enumerate(connections)]
    qd_bus_unmerged = [JuMP.@NLexpression(pm.model, -vr[p]*cid_bus_unmerged[idx]+vi[p]*crd_bus_unmerged[idx]) for (idx,p) in enumerate(connections)]

    var(pm, nw, :pd_bus)[id] = pd_bus = _merge_bus_flows(pm, pd_bus_unmerged, connections)
    var(pm, nw, :qd_bus)[id] = qd_bus = _merge_bus_flows(pm, qd_bus_unmerged, connections)

    if report
        sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        pd = Vector{JuMP.NonlinearExpression}([])
        qd = Vector{JuMP.NonlinearExpression}([])
        for idx in 1:nph
            push!(pd, JuMP.@NLexpression(pm.model, a[idx]*(vrd[idx]^2+vid[idx]^2)^(alpha[idx]/2) ))
            push!(qd, JuMP.@NLexpression(pm.model, b[idx]*(vrd[idx]^2+vid[idx]^2)^(beta[idx]/2)  ))
        end
        sol(pm, nw, :load, id)[:pd] = JuMP.Containers.DenseAxisArray(pd, connections)
        sol(pm, nw, :load, id)[:qd] = JuMP.Containers.DenseAxisArray(qd, connections)
    end
end


# ""
# function constraint_mc_load_current(pm::AbstractMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
#     load = ref(pm, nw, :load, id)
#     bus = ref(pm, nw,:bus, load["load_bus"])

#     configuration = load["configuration"]

#     a, alpha, b, beta = _load_expmodel_params(load, bus)

#     int_dim = _infer_int_dim_unit(load, false)
#     if configuration==WYE || int_dim==1
#         constraint_mc_load_current_wye(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
#     else
#         constraint_mc_load_current_delta(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
#     end
# end


# ""
# function constraint_mc_load_current_wye(pm::AbstractMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
#     vr = var(pm, nw, :vr, bus_id)
#     vi = var(pm, nw, :vi, bus_id)

#     crd = Vector{JuMP.NonlinearExpression}([])
#     cid = Vector{JuMP.NonlinearExpression}([])

#     phases = connections[1:end-1]
#     n      = connections[end]

#     for (idx, p) in enumerate(phases)
#         push!(crd, JuMP.@NLexpression(pm.model,
#              a[idx]*(vr[p]-vr[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(alpha[idx]/2-1)
#             +b[idx]*(vi[p]-vi[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(beta[idx]/2 -1)
#         ))
#         push!(cid, JuMP.@NLexpression(pm.model,
#              a[idx]*(vi[p]-vi[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(alpha[idx]/2-1)
#             -b[idx]*(vr[p]-vr[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(beta[idx]/2 -1)
#         ))
#     end

#     var(pm, nw, :crd)[id] = crd
#     var(pm, nw, :cid)[id] = cid

#     crd_bus_n = JuMP.@NLexpression(pm.model, -sum(crd[i] for i in 1:length(phases)))
#     cid_bus_n = JuMP.@NLexpression(pm.model, -sum(cid[i] for i in 1:length(phases)))

#     var(pm, nw, :crd_bus)[id] = crd_bus = _merge_bus_flows(pm, [crd..., crd_bus_n], connections)
#     var(pm, nw, :cid_bus)[id] = cid_bus = _merge_bus_flows(pm, [cid..., cid_bus_n], connections)

#     if report
#         pd_bus = Vector{JuMP.NonlinearExpression}([])
#         qd_bus = Vector{JuMP.NonlinearExpression}([])
#         for (idx,c) in enumerate(connections)
#             push!(pd_bus, JuMP.@NLexpression(pm.model,  vr[c]*crd_bus[c]+vi[c]*cid_bus[c]))
#             push!(qd_bus, JuMP.@NLexpression(pm.model, -vr[c]*cid_bus[c]+vi[c]*crd_bus[c]))
#         end

#         sol(pm, nw, :load, id)[:pd_bus] = JuMP.Containers.DenseAxisArray(pd_bus, connections)
#         sol(pm, nw, :load, id)[:qd_bus] = JuMP.Containers.DenseAxisArray(qd_bus, connections)

#         sol(pm, nw, :load, id)[:crd] = JuMP.Containers.DenseAxisArray(crd, connections)
#         sol(pm, nw, :load, id)[:cid] = JuMP.Containers.DenseAxisArray(cid, connections)

#         sol(pm, nw, :load, id)[:crd_bus] = crd_bus
#         sol(pm, nw, :load, id)[:cid_bus] = cid_bus

#         pd = Vector{JuMP.NonlinearExpression}([])
#         qd = Vector{JuMP.NonlinearExpression}([])
#         for (idx, p) in enumerate(phases)
#             push!(pd, JuMP.@NLexpression(pm.model, a[idx]*(vr[p]^2+vi[p]^2)^(alpha[idx]/2) ))
#             push!(qd, JuMP.@NLexpression(pm.model, b[idx]*(vr[p]^2+vi[p]^2)^(beta[idx]/2)  ))
#         end
#         sol(pm, nw, :load, id)[:pd] = JuMP.Containers.DenseAxisArray(pd, connections)
#         sol(pm, nw, :load, id)[:qd] = JuMP.Containers.DenseAxisArray(qd, connections)
#     end
# end


# ""
# function constraint_mc_load_current_delta(pm::AbstractMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
#     vr = var(pm, nw, :vr, bus_id)
#     vi = var(pm, nw, :vi, bus_id)

#     nph = length(connections)
#     prev = Dict(i=>(i+nph-2)%nph+1 for i in 1:nph)
#     next = Dict(i=>i%nph+1 for i in 1:nph)

#     vrd = JuMP.@NLexpression(pm.model, [i in 1:nph], vr[i]-vr[next[i]])
#     vid = JuMP.@NLexpression(pm.model, [i in 1:nph], vi[i]-vi[next[i]])

#     crd = JuMP.@NLexpression(pm.model, [i in 1:nph],
#         a[i]*vrd[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2-1)
#        +b[i]*vid[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2 -1)
#     )
#     cid = JuMP.@NLexpression(pm.model, [i in 1:nph],
#         a[i]*vid[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2-1)
#        -b[i]*vrd[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2 -1)
#     )

#     crd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], crd[i]-crd[prev[i]])
#     cid_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], cid[i]-cid[prev[i]])

#     var(pm, nw, :crd_bus)[id] = _merge_bus_flows(pm, crd_bus, connections)
#     var(pm, nw, :cid_bus)[id] = _merge_bus_flows(pm, cid_bus, connections)

#     if report
#         pd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph],  vr[i]*crd_bus[i]+vi[i]*cid_bus[i])
#         qd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], -vr[i]*cid_bus[i]+vi[i]*crd_bus[i])

#         sol(pm, nw, :load, id)[:pd_bus] = pd_bus
#         sol(pm, nw, :load, id)[:qd_bus] = qd_bus

#         pd = JuMP.@NLexpression(pm.model, [i in 1:nph], a[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2) )
#         qd = JuMP.@NLexpression(pm.model, [i in 1:nph], b[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2)  )
#         sol(pm, nw, :load, id)[:pd] = pd
#         sol(pm, nw, :load, id)[:qd] = qd
#     end
# end


# # LOAD - Constraints - Quadratic

# ""
# function constraint_mc_load_power(pm::AbstractQuadraticMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
#     load = ref(pm, nw, :load, id)
#     bus = ref(pm, nw,:bus, load["load_bus"])

#     configuration = load["configuration"]
#     int_dim = _infer_int_dim_unit(load, false)
#     a, alpha, b, beta = _load_expmodel_params(load, bus)

#     if configuration==WYE || int_dim==1
#         constraint_mc_load_power_wye(pm, nw, id, load["load_bus"], load["connections"], load["model"], a, b; report=report)
#     else
#         constraint_mc_load_power_delta(pm, nw, id, load["load_bus"], load["connections"], load["model"], a, b; report=report)
#     end
# end


# ""
# function constraint_mc_load_power_wye(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, model::LoadModel, a::Vector{<:Real}, b::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
#     vr = var(pm, nw, :vr, bus_id)
#     vi = var(pm, nw, :vi, bus_id)

#     crd = var(pm, nw, :crd, id)
#     cid = var(pm, nw, :cid, id)

#     phases = connections[1:end-1]
#     n      = connections[end]

#     vr_pn = [vr[p]-vr[n] for p in phases]
#     vi_pn = [vi[p]-vi[n] for p in phases]

#     if model==POWER
#         pd = a
#         qd = b
#     elseif model==IMPEDANCE
#         pd = a .* (vr_pn.^2 .+ vi_pn.^2)
#         qd = b .* (vr_pn.^2 .+ vi_pn.^2)
#     elseif model==CURRENT
#         pd = var(pm, nw, :pd, id)
#         qd = var(pm, nw, :qd, id)
#         JuMP.@constraint(pm.model, pd.^2 .== a.^2 .* (vr_pn.^2 .+ vi_pn.^2))
#         JuMP.@constraint(pm.model, sign.(a).*pd .>= 0) #TODO handle in bounds instead
#         JuMP.@constraint(pm.model, qd.^2 .== b.^2 .* (vr_pn.^2 .+ vi_pn.^2))
#         JuMP.@constraint(pm.model, sign.(b).*qd .>= 0) #TODO handle in bounds instead
#     else
#         error("Load model $model for load $id is not supported by this formulation.")
#     end

#     JuMP.@constraint(pm.model, pd .==  vr_pn.*crd .+ vi_pn.*cid)
#     JuMP.@constraint(pm.model, qd .== -vr_pn.*cid .+ vi_pn.*crd)

#     # constant current loads are already reported through variable function
#     if report && model!=CURRENT
#         sol(pm, nw, :load, id)[:pd] = pd
#         sol(pm, nw, :load, id)[:qd] = qd
#     end
# end


# ""
# function constraint_mc_load_power_delta(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, model::LoadModel, a::Vector{<:Real}, b::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
#     vr = var(pm, nw, :vr, bus_id)
#     vi = var(pm, nw, :vi, bus_id)

#     crd = var(pm, nw, :crd, id)
#     cid = var(pm, nw, :cid, id)

#     phases = connections[1:end-1]
#     n      = connections[end]

#     Md = _get_delta_transformation_matrix(length(connections))
#     vrd = Md*[vr[p] for p in phases]
#     vid = Md*[vi[p] for p in phases]

#     if model==POWER
#         pd = a
#         qd = b
#     elseif model==IMPEDANCE
#         pd = a .* (vrd.^2 .+ vid.^2)
#         qd = b .* (vrd.^2 .+ vid.^2)
#     elseif model==CURRENT
#         pd = var(pm, nw, :pd, id)
#         qd = var(pm, nw, :qd, id)
#         JuMP.@constraint(pm.model, pd.^2 .== a.^2 .* (vrd.^2 .+ vid.^2))
#         JuMP.@constraint(pm.model, sign.(a).*pd .>= 0) #TODO handle in bounds instead
#         JuMP.@constraint(pm.model, qd.^2 .== b.^2 .* (vrd.^2 .+ vid.^2))
#         JuMP.@constraint(pm.model, sign.(b).*qd .>= 0) #TODO handle in bounds instead
#     else
#         error("Load model $model for load $id is not supported by this formulation.")
#     end

#     JuMP.@constraint(pm.model, pd .==  vrd.*crd .+ vid.*cid)
#     JuMP.@constraint(pm.model, qd .== -vrd.*cid .+ vid.*crd)

#     # constant current loads are already reported through variable function
#     if report && model!=CURRENT
#         sol(pm, nw, :load, id)[:pd] = pd
#         sol(pm, nw, :load, id)[:qd] = qd
#     end
# end


# function constraint_mc_load_current(pm::AbstractQuadraticMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)
#     load = ref(pm, nw, :load, id)

#     int_dim = _infer_int_dim_unit(load, false)
#     if get(load, "configuration", WYE) == WYE || int_dim==1
#         constraint_mc_load_current_wye(pm, nw, id, load["connections"]; report=report, bounded=bounded)
#     else
#         constraint_mc_load_current_delta(pm, nw, id, load["connections"]; report=report, bounded=bounded)
#     end
# end


# ""
# function constraint_mc_load_current_wye(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
#     crd = var(pm, nw, :crd, id)
#     cid = var(pm, nw, :cid, id)
#     var(pm, nw, :crd_bus)[id] = _merge_bus_flows(pm, [crd..., -sum(crd)], connections)
#     var(pm, nw, :cid_bus)[id] = _merge_bus_flows(pm, [cid..., -sum(cid)], connections)
# end


# ""
# function constraint_mc_load_current_delta(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
#     crd = var(pm, nw, :crd, id)
#     cid = var(pm, nw, :cid, id)
#     Md = _get_delta_transformation_matrix(length(connections))
#     var(pm, nw, :crd_bus)[id] = _merge_bus_flows(pm, Md'*crd, connections)
#     var(pm, nw, :cid_bus)[id] = _merge_bus_flows(pm, Md'*cid, connections)
# end


# # TRANSFORMER

# # Variables

# ""
# function variable_mc_transformer_current(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
#     variable_mc_transformer_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
#     variable_mc_transformer_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

#     var(pm)[:crt_bus] = Dict{Tuple{Int,Int,Int}, Any}()
#     var(pm)[:cit_bus] = Dict{Tuple{Int,Int,Int}, Any}()
# end


# "variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
# function variable_mc_transformer_current_real(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in ref(pm, :transformer))
#     crt = var(pm, nw)[:crt] = Dict((l,i,j) => JuMP.@variable(pm.model,
#             [c in 1:int_dim[l]], base_name="$(nw)_crt_$((l,i,j))",
#             start = comp_start_value(ref(pm, nw, :transformer, l), "crt_start", c, 0.0)
#         ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
#     )

#     report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :crt_fr, :crt_to, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), crt)
# end


# "variable: `ci[l,i,j] ` for `(l,i,j)` in `arcs`"
# function variable_mc_transformer_current_imaginary(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in ref(pm, :transformer))
#     cit = var(pm, nw)[:cit] = Dict((l,i,j) => JuMP.@variable(pm.model,
#             [c in 1:int_dim[l]], base_name="$(nw)_cit_$((l,i,j))",
#             start = comp_start_value(ref(pm, nw, :transformer, l), "cit_start", c, 0.0)
#         ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
#     )

#     report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :cit_fr, :cit_to, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), cit)
# end


# ""
# function variable_mc_transformer_power_active(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in ref(pm, :transformer))
#     pt = var(pm, nw)[:pt] = Dict((l,i,j) => JuMP.@variable(pm.model,
#             [c in 1:int_dim[l]], base_name="$(nw)_pt_$((l,i,j))",
#             start = comp_start_value(ref(pm, nw, :transformer, l), "pt_start", c, 0.0)
#         ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
#     )
    
#     if bounded
#         for (l,i,j) in ref(pm, nw, :arcs_transformer_from)
#             trans = ref(pm, nw, :transformer, l)
#             f_bus = ref(pm, nw, :bus, i)
#             t_bus = ref(pm, nw, :bus, j)
#             sm_ub = trans["sm_ub"]
#             set_lower_bound(pt[(l,i,j)], -sm_ub)
#             set_upper_bound(pt[(l,i,j)],  sm_ub)
#             set_lower_bound(pt[(l,j,i)], -sm_ub)
#             set_upper_bound(pt[(l,j,i)],  sm_ub)
#         end
#     end

#     report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :pt_fr, :pt_to, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), pt)
# end


# ""
# function variable_mc_transformer_power_reactive(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in ref(pm, :transformer))
#     qt = var(pm, nw)[:qt] = Dict((l,i,j) => JuMP.@variable(pm.model,
#             [c in 1:int_dim[l]], base_name="$(nw)_qt_$((l,i,j))",
#             start = comp_start_value(ref(pm, nw, :transformer, l), "qt_start", c, 0.0)
#         ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
#     )
    
#     if bounded
#         for (l,i,j) in ref(pm, nw, :arcs_transformer_from)
#             trans = ref(pm, nw, :transformer, l)
#             f_bus = ref(pm, nw, :bus, i)
#             t_bus = ref(pm, nw, :bus, j)
#             sm_ub = trans["sm_ub"]
#             set_lower_bound(qt[(l,i,j)], -sm_ub)
#             set_upper_bound(qt[(l,i,j)],  sm_ub)
#             set_lower_bound(qt[(l,j,i)], -sm_ub)
#             set_upper_bound(qt[(l,j,i)],  sm_ub)
#         end
#     end

#     report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :qt_fr, :qt_to, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), qt)
# end


# # TRANSFORMER -  Variable - Non-linear

# ""
# function variable_mc_transformer_power(pm::AbstractNLMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
#     # do nothing    
# end


# # TRANSFORMER - Variable - Quadratic

""
function variable_mc_transformer_power(pm::AbstractMultiConductorACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_transformer_power_active(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_transformer_power_reactive(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    var(pm, nw)[:pt_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:qt_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end

# # Constraints

"This function adds all constraints required to model a two-winding, wye-wye connected transformer."
function constraint_mc_transformer_power_yy(pm::AbstractMultiConductorACRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    pt_fr = var(pm, nw, :pt, f_idx)
    pt_to = var(pm, nw, :pt, t_idx)
    qt_fr = var(pm, nw, :qt, f_idx)
    qt_to = var(pm, nw, :qt, t_idx)
    vr_fr = var(pm, nw, :vr, f_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    JuMP.@constraint(pm.model, pt_fr + pt_to .== 0)
    JuMP.@constraint(pm.model, qt_fr + qt_to .== 0)

    f_phases = f_connections[1:end-1]
    f_n      = f_connections[end]
    t_phases = t_connections[1:end-1]
    t_n      = t_connections[end]
    P        = length(f_phases)

    # from-side
    if iszero(vr_fr[f_n]) && iszero(vi_fr[f_n])
        pt_bus_fr_unmerged = [pt_fr..., 0.0]
        qt_bus_fr_unmerged = [qt_fr..., 0.0]
    else
        vr_fr_pn = [vr_fr[p]-vr_fr[f_n] for p in f_phases]
        vi_fr_pn = [vi_fr[p]-vi_fr[f_n] for p in f_phases]
        
        crt_fr = [JuMP.@NLexpression(pm.model, ( pt_fr[idx]*vr_fr_pn[idx]+qt_fr[idx]*vi_fr_pn[idx])/(vr_fr_pn[idx]^2+vi_fr_pn[idx]^2)) for idx in 1:P]
        cit_fr = [JuMP.@NLexpression(pm.model, (-pt_fr[idx]*vi_fr_pn[idx]+qt_fr[idx]*vr_fr_pn[idx])/(vr_fr_pn[idx]^2+vi_fr_pn[idx]^2)) for idx in 1:P]

        crt_fr_bus_unmerged = 

        pt_bus_fr_unmerged = [
            [JuMP.@NLexpression(pm.model, vr_fr[p]*crt_fr[idx]+vi_fr[p]*cit_fr[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, vr_fr[n]*sum(-crt_fr[idx] for idx in 1:P)+vi_fr[n]*sum(-cit_fr[idx] for idx in 1:P))
        ]
        qt_bus_fr_unmerged = [
            [JuMP.@NLexpression(pm.model, -vr_fr[p]*cit_fr[idx]+vi_fr[p]*crt_fr[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, -vr_fr[n]*sum(-cit_fr[idx] for idx in 1:P)+vi_fr[n]*sum(-crt_fr[idx] for idx in 1:P))
        ]
    end
    var(pm, nw, :pt_bus)[f_idx] = pt_bus_fr = _merge_bus_flows(pm, pt_bus_fr_unmerged, f_connections)
    var(pm, nw, :qt_bus)[f_idx] = qt_bus_fr = _merge_bus_flows(pm, qt_bus_fr_unmerged, f_connections)

    # to-side
    if iszero(vr_to[t_n]) && iszero(vi_to[t_n])
        pt_bus_to_unmerged = [pt_to..., 0.0]
        qt_bus_to_unmerged = [qt_to..., 0.0]
    else
        vr_to_pn = [vr_to[p]-vr_to[t_n] for p in t_phases]
        vi_to_pn = [vi_to[p]-vi_to[t_n] for p in t_phases]
        
        crt_to = [JuMP.@NLexpression(pm.model, ( pt_to[idx]*vr_to_pn[idx]+qt_to[idx]*vi_to_pn[idx])/(vr_to_pn[idx]^2+vi_to_pn[idx]^2)) for idx in 1:P]
        cit_to = [JuMP.@NLexpression(pm.model, (-pt_to[idx]*vi_to_pn[idx]+qt_to[idx]*vr_to_pn[idx])/(vr_to_pn[idx]^2+vi_to_pn[idx]^2)) for idx in 1:P]

        pt_bus_to_unmerged = [
            [JuMP.@NLexpression(pm.model, vr_to[p]*crt_to[idx]+vi_to[p]*cit_to[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, vr_to[n]*sum(-crt_to[idx] for idx in 1:P)+vi_to[n]*sum(-cit_to[idx] for idx in 1:P))
        ]
        qt_bus_to_unmerged = [
            [JuMP.@NLexpression(pm.model, -vr_to[p]*cit_to[idx]+vi_to[p]*crt_to[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, -vr_to[n]*sum(-cit_to[idx] for idx in 1:P)+vi_to[n]*sum(-crt_to[idx] for idx in 1:P))
        ]
    end
    var(pm, nw, :pt_bus)[t_idx] = pt_bus_to = _merge_bus_flows(pm, pt_bus_to_unmerged, t_connections)
    var(pm, nw, :qt_bus)[t_idx] = qt_bus_to = _merge_bus_flows(pm, qt_bus_to_unmerged, t_connections)
end

"This function adds all constraints required to model a two-winding, wye-wye connected transformer."
function constraint_mc_transformer_power_dy(pm::AbstractMultiConductorACRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    pt_fr = var(pm, nw, :pt, f_idx)
    pt_to = var(pm, nw, :pt, t_idx)
    qt_fr = var(pm, nw, :qt, f_idx)
    qt_to = var(pm, nw, :qt, t_idx)
    vr_fr = var(pm, nw, :vr, f_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    JuMP.@constraint(pm.model, pt_fr + pt_to .== 0)
    JuMP.@constraint(pm.model, qt_fr + qt_to .== 0)

    t_phases = t_connections[1:end-1]
    t_n      = t_connections[end]
    P        = length(t_phases)
    next(i)  = i%P+1 # 1=>2, 2=>3, 3=>1
    prev(i)  = (i-2+P)%P+1  # 1=>3, 2=>1, 3=>2

    # from-side
    vrg_fr = [vr_fr[p]-vr_fr[f_connections[next(idx)]] for (idx,p) in enumerate(f_connections)]
    vig_fr = [vi_fr[p]-vi_fr[f_connections[next(idx)]] for (idx,p) in enumerate(f_connections)]
    
    crt_fr = JuMP.@NLexpression(pm.model, [idx in 1:P], (pt_fr[idx]*vrg_fr[idx]+qt_fr[idx]*vig_fr[idx])/(vrg_fr[idx]^2+vig_fr[idx]^2))
    cit_fr = JuMP.@NLexpression(pm.model, [idx in 1:P], (pt_fr[idx]*vig_fr[idx]-qt_fr[idx]*vrg_fr[idx])/(vrg_fr[idx]^2+vig_fr[idx]^2))

    crt_bus_fr_unmerged = JuMP.@NLexpression(pm.model, [idx in 1:P], crt_fr[idx]-crt_fr[prev(idx)])
    cit_bus_fr_unmerged = JuMP.@NLexpression(pm.model, [idx in 1:P], cit_fr[idx]-cit_fr[prev(idx)])

    pt_bus_fr_unmerged = [JuMP.@NLexpression(pm.model,  vr_fr[p]*crt_bus_fr_unmerged[idx]+vi_fr[p]*cit_bus_fr_unmerged[idx]) for (idx,p) in enumerate(f_connections)]
    qt_bus_fr_unmerged = [JuMP.@NLexpression(pm.model, -vr_fr[p]*cit_bus_fr_unmerged[idx]+vi_fr[p]*crt_bus_fr_unmerged[idx]) for (idx,p) in enumerate(f_connections)]

    var(pm, nw, :pt_bus)[f_idx] = pt_bus_fr = _merge_bus_flows(pm, pt_bus_fr_unmerged, f_connections)
    var(pm, nw, :qt_bus)[f_idx] = qt_bus_fr = _merge_bus_flows(pm, qt_bus_fr_unmerged, f_connections)

    # to-side
    if iszero(vr_to[t_n]) && iszero(vi_to[t_n])
        pt_bus_to_unmerged = [pt_to..., 0.0]
        qt_bus_to_unmerged = [qt_to..., 0.0]
    else
        vr_to_pn = [vr_to[p]-vr_to[t_n] for p in t_phases]
        vi_to_pn = [vi_to[p]-vi_to[t_n] for p in t_phases]
        
        crt_to = [JuMP.@NLexpression(pm.model, ( pt_to[idx]*vr_to_pn[idx]+qt_to[idx]*vi_to_pn[idx])/(vr_to_pn[idx]^2+vi_to_pn[idx]^2)) for idx in 1:P]
        cit_to = [JuMP.@NLexpression(pm.model, (-pt_to[idx]*vi_to_pn[idx]+qt_to[idx]*vr_to_pn[idx])/(vr_to_pn[idx]^2+vi_to_pn[idx]^2)) for idx in 1:P]

        pt_bus_to_unmerged = [
            [JuMP.@NLexpression(pm.model, vr_to[p]*crt_to[idx]+vi_to[p]*cit_to[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, vr_to[n]*sum(-crt_to[idx] for idx in 1:P)+vi_to[n]*sum(-cit_to[idx] for idx in 1:P))
        ]
        qt_bus_to_unmerged = [
            [JuMP.@NLexpression(pm.model, -vr_to[p]*cit_to[idx]+vi_to[p]*crt_to[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, -vr_to[n]*sum(-cit_to[idx] for idx in 1:P)+vi_to[n]*sum(-crt_to[idx] for idx in 1:P))
        ]
    end
    var(pm, nw, :pt_bus)[t_idx] = pt_bus_to = _merge_bus_flows(pm, pt_bus_to_unmerged, t_connections)
    var(pm, nw, :qt_bus)[t_idx] = qt_bus_to = _merge_bus_flows(pm, qt_bus_to_unmerged, t_connections)
end


""
function constraint_mc_transformer_power_rating(pm::AbstractMultiConductorACRModel, nw::Int, id::Int, f_idx::Tuple, t_idx::Tuple, f_bus::Int, t_bus::Int, f_connections::Vector, t_connections::Vector, config::ConnConfig, sm_ub::Real; report::Bool=true)
    pt_fr = var(pm, nw, :pt, f_idx)
    qt_fr = var(pm, nw, :qt, f_idx)
    pt_to = var(pm, nw, :pt, t_idx)
    qt_to = var(pm, nw, :qt, t_idx)

    if sm_ub<Inf
        JuMP.@NLconstraint(pm.model, sum(pt_fr[i] for i in idxs)^2 + sum(qt_fr[i] for i in idxs)^2 <= sm_ub^2)
        JuMP.@NLconstraint(pm.model, sum(pt_to[i] for i in idxs)^2 + sum(qt_to[i] for i in idxs)^2 <= sm_ub^2)
    end
end


# ""
# function constraint_mc_transformer_power_rating(pm::AbstractMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
#     trans = ref(pm, nw, :transformer, id)
#     f_bus = trans["f_bus"]
#     t_bus = trans["t_bus"]
#     f_idx = (id,f_bus,t_bus)
#     t_idx = (id,t_bus,f_bus)
#     f_conns = trans["f_connections"]
#     t_conns = trans["t_connections"]
#     config = trans["configuration"]
#     sm_ub = trans["sm_ub"]

#     constraint_mc_transformer_power_rating(pm, nw, id, f_idx, t_idx, f_bus, t_bus, f_conns, t_conns, config, sm_ub)
# end


# ""
# function constraint_mc_transformer_power_rating(pm::AbstractMultiConductorIVRModel, nw::Int, id::Int, f_idx::Tuple, t_idx::Tuple, f_bus::Int, t_bus::Int, f_connections::Vector, t_connections::Vector, config::ConnConfig, sm_ub::Real; report::Bool=true)
#     vr_fr = var(pm, nw, :vr, f_bus)
#     vi_fr = var(pm, nw, :vi, f_bus)
#     vr_to = var(pm, nw, :vr, t_bus)
#     vi_to = var(pm, nw, :vi, t_bus)

#     crt_fr = var(pm, nw, :crt, f_idx)
#     cit_fr = var(pm, nw, :cit, f_idx)
#     crt_to = var(pm, nw, :crt, t_idx)
#     cit_to = var(pm, nw, :cit, t_idx)

#     if config==WYE || length(crt_fr)==1
#         P_fr = f_connections[1:end-1]
#         n_fr = f_connections[end]
#         vrt_fr = [vr_fr[p]-vr_fr[n_fr] for p in P_fr]
#         vit_fr = [vi_fr[p]-vi_fr[n_fr] for p in P_fr]
#     elseif config==DELTA && length(crt_fr)==3
#         M = _get_delta_transformation_matrix(3)
#         vrt_fr = M*[vr_to[p] for p in f_connections]
#         vit_fr = M*[vi_to[p] for p in f_connections]
#     else
#         error("The configuration $config of dimension $(length(crt)) is not supported.")
#     end

#     P_to = t_connections[1:end-1]
#     n_to = t_connections[end]
#     vrt_to = [vr_to[p]-vr_to[n_to] for p in P_to]
#     vit_to = [vi_to[p]-vi_to[n_to] for p in P_to]

#     idxs = [1:length(vrt_fr)...]
#     pt_fr = JuMP.@NLexpression(pm.model, [i in idxs],  vrt_fr[i]*crt_fr[i] + vit_fr[i]*cit_fr[i])
#     qt_fr = JuMP.@NLexpression(pm.model, [i in idxs], -vrt_fr[i]*cit_fr[i] + vit_fr[i]*crt_fr[i])
#     pt_to = JuMP.@NLexpression(pm.model, [i in idxs],  vrt_to[i]*crt_to[i] + vit_to[i]*cit_to[i])
#     qt_to = JuMP.@NLexpression(pm.model, [i in idxs], -vrt_to[i]*cit_to[i] + vit_to[i]*crt_to[i])

#     if sm_ub<Inf
#         JuMP.@NLconstraint(pm.model, sum(pt_fr[i] for i in idxs)^2 + sum(qt_fr[i] for i in idxs)^2 <= sm_ub^2)
#         JuMP.@NLconstraint(pm.model, sum(pt_to[i] for i in idxs)^2 + sum(qt_to[i] for i in idxs)^2 <= sm_ub^2)
#     end

#     if report
#         sol(pm, nw, :transformer, id)[:pt_fr] = pt_fr
#         sol(pm, nw, :transformer, id)[:qt_fr] = qt_fr
#         sol(pm, nw, :transformer, id)[:pt_to] = pt_to
#         sol(pm, nw, :transformer, id)[:qt_to] = qt_to
#         sol(pm, nw, :transformer, id)[:smtot_fr] = JuMP.@NLexpression(pm.model, sqrt(sum(pt_fr[i] for i in idxs)^2 + sum(qt_fr[i] for i in idxs)^2))
#         sol(pm, nw, :transformer, id)[:smtot_to] = JuMP.@NLexpression(pm.model, sqrt(sum(pt_to[i] for i in idxs)^2 + sum(qt_to[i] for i in idxs)^2))
#     end
# end


# # TRANSFORMER - Constraint - Quadratic

# ""
# function constraint_mc_transformer_power_rating(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, f_idx::Tuple, t_idx::Tuple, f_bus::Int, t_bus::Int, f_connections::Vector, t_connections::Vector, config::ConnConfig, sm_ub::Real; report::Bool=true)
#     vr_fr = var(pm, nw, :vr, f_bus)
#     vi_fr = var(pm, nw, :vi, f_bus)
#     vr_to = var(pm, nw, :vr, t_bus)
#     vi_to = var(pm, nw, :vi, t_bus)

#     crt_fr = var(pm, nw, :crt, f_idx)
#     cit_fr = var(pm, nw, :cit, f_idx)
#     crt_to = var(pm, nw, :crt, t_idx)
#     cit_to = var(pm, nw, :cit, t_idx)

#     if config==WYE || length(crt_fr)==1
#         P_fr = f_connections[1:end-1]
#         n_fr = f_connections[end]
#         vrt_fr = [vr_fr[p]-vr_fr[n_fr] for p in P_fr]
#         vit_fr = [vi_fr[p]-vi_fr[n_fr] for p in P_fr]
#     elseif config==DELTA && length(crt_fr)==3
#         M = _get_delta_transformation_matrix(3)
#         vrt_fr = M*[vr_to[p] for p in f_connections]
#         vit_fr = M*[vi_to[p] for p in f_connections]
#     else
#         error("The configuration $config of dimension $(length(crt)) is not supported.")
#     end

#     P_to = t_connections[1:end-1]
#     n_to = t_connections[end]
#     vrt_to = [vr_to[p]-vr_to[n_to] for p in P_to]
#     vit_to = [vi_to[p]-vi_to[n_to] for p in P_to]

#     pt_fr = var(pm, nw, :pt, f_idx)
#     qt_fr = var(pm, nw, :qt, f_idx)
#     pt_to = var(pm, nw, :pt, t_idx)
#     qt_to = var(pm, nw, :qt, t_idx)

#     JuMP.@constraint(pm.model, pt_fr .==  vrt_fr.*crt_fr .+ vit_fr.*cit_fr)
#     JuMP.@constraint(pm.model, qt_fr .== -vrt_fr.*cit_fr .+ vit_fr.*crt_fr)
#     JuMP.@constraint(pm.model, pt_to .==  vrt_to.*crt_to .+ vit_to.*cit_to)
#     JuMP.@constraint(pm.model, qt_to .== -vrt_to.*cit_to .+ vit_to.*crt_to)

#     if sm_ub<Inf
#         JuMP.@constraint(pm.model, sum(pt_fr)^2+sum(qt_fr)^2 <= sm_ub^2)
#         JuMP.@constraint(pm.model, sum(pt_to)^2+sum(qt_to)^2 <= sm_ub^2)
#     end
    
#     # if report
#     #     sol(pm, nw, :transformer, id)[:sm_fr] = JuMP.@NLexpression(pm.model, sqrt(pt_fr^2+qt_fr^2))
#     #     sol(pm, nw, :transformer, id)[:sm_to] = JuMP.@NLexpression(pm.model, sqrt(pt_to^2+qt_to^2))
#     # end
# end


# # BRANCH

# # BRANCH - Variables

""
function variable_mc_branch_power(pm::AbstractMultiConductorACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_branch_power_active(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_branch_power_reactive(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    var(pm, nw)[:p_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:q_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end


""
function variable_mc_branch_power_active(pm::MultiConductorModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref(pm, nw, :branch))
    p = var(pm, nw)[:p] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_p_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "p_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch, :p_fr, :p_to, ref(pm, nw, :arcs_branch_from), ref(pm, nw, :arcs_branch_to), p)
end


""
function variable_mc_branch_power_reactive(pm::MultiConductorModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref(pm, nw, :branch))
    q = var(pm, nw)[:q] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_q_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "q_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch, :q_fr, :q_to, ref(pm, nw, :arcs_branch_from), ref(pm, nw, :arcs_branch_to), q)
end


# "variable: `csr[l]` for `l` in `branch`"
# function variable_mc_branch_current_series_real(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref(pm, nw, :branch))
#     csr = var(pm, nw)[:csr] = Dict(l => JuMP.@variable(pm.model,
#             [c in 1:nconds[l]], base_name="$(nw)_csr_$(l)",
#             start = comp_start_value(ref(pm, nw, :branch, l), "csr_start", c, 0.0)
#         ) for (l,i,j) in ref(pm, nw, :arcs_branch)
#     )

#     # if bounded
#     #     for (l,i,j) in ref(pm, nw, :arcs_branch_from)
#     #         cmax = _calc_branch_series_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
#     #         for (idx,c) in enumerate(connections[(l,i,j)])
#     #             set_upper_bound(csr[l][c],  cmax[idx])
#     #             set_lower_bound(csr[l][c], -cmax[idx])
#     #         end
#     #     end
#     # end

#     report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch, :csr_fr, ids(pm, nw, :branch), csr)
# end


# "variable: `csi[l]` for `l` in `branch`"
# function variable_mc_branch_current_series_imaginary(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref(pm, nw, :branch))
#     csi = var(pm, nw)[:csi] = Dict(l => JuMP.@variable(pm.model,
#             [c in 1:nconds[l]], base_name="$(nw)_csi_$(l)",
#             start = comp_start_value(ref(pm, nw, :branch, l), "csi_start", c, 0.0)
#         ) for (l,i,j) in ref(pm, nw, :arcs_branch)
#     )

#     # if bounded
#     #     for (l,i,j) in ref(pm, nw, :arcs_branch_from)
#     #         cmax = _calc_branch_series_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
#     #         for (idx,c) in enumerate(connections[(l,i,j)])
#     #             set_upper_bound(csi[l][c],  cmax[idx])
#     #             set_lower_bound(csi[l][c], -cmax[idx])
#     #         end
#     #     end
#     # end

#     report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch, :csi_fr, ids(pm, nw, :branch), csi)
# end


# # BRANCH - Variables - Reduced

# ""
# function variable_mc_branch_current(pm::ReducedMultiConductorIVRModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
#     variable_mc_branch_current_series_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
#     variable_mc_branch_current_series_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
#     var(pm, nw)[:cr] = Dict{Tuple{Int,Int,Int}, Any}()
#     var(pm, nw)[:ci] = Dict{Tuple{Int,Int,Int}, Any}()
#     var(pm, nw)[:cr_bus] = Dict{Tuple{Int,Int,Int}, Any}()
#     var(pm, nw)[:ci_bus] = Dict{Tuple{Int,Int,Int}, Any}()
# end


# # BRANCH - Constraints

"""
Creates Ohms constraints

s_fr = v_fr.*conj(Y*(v_fr-v_to))
s_fr = (vr_fr+im*vi_fr).*(G-im*B)*([vr_fr-vr_to]-im*[vi_fr-vi_to])
s_fr = (vr_fr+im*vi_fr).*([G*vr_fr-G*vr_to-B*vi_fr+B*vi_to]-im*[G*vi_fr-G*vi_to+B*vr_fr-B*vr_to])
"""
function constraint_mc_ohms_yt_from(pm::AbstractMultiConductorACRModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})
    p_fr  = var(pm, nw, :p, f_idx)
    q_fr  = var(pm, nw, :q, f_idx)
    vr_fr = [var(pm, nw, :vr, f_bus)[t] for t in f_connections]
    vr_to = [var(pm, nw, :vr, t_bus)[t] for t in t_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[t] for t in f_connections]
    vi_to = [var(pm, nw, :vi, t_bus)[t] for t in t_connections]

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

    var(pm, nw, :p_bus)[f_idx] = _merge_bus_flows(pm, p_fr, f_connections)
    var(pm, nw, :q_bus)[f_idx] = _merge_bus_flows(pm, q_fr, f_connections)
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_mc_ohms_yt_to(pm::AbstractMultiConductorACRModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix, B::Matrix, G_to::Matrix, B_to::Matrix)
    constraint_mc_ohms_yt_from(pm, nw, t_bus, f_bus, t_idx, f_idx, t_connections, f_connections, G, B, G_to, B_to)
end


""
function constraint_mc_branch_current_rating(pm::AbstractMultiConductorACRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int},  f_connections::Vector, t_connections::Vector, c_rating::Vector{<:Real}; report::Bool=true)
    p_fr = var(pm, nw, :p, f_idx)
    q_fr = var(pm, nw, :q, f_idx)
    p_to = var(pm, nw, :p, t_idx)
    q_to = var(pm, nw, :q, t_idx)
    vr_fr   = var(pm, nw, :vr, f_idx[2])
    vi_fr   = var(pm, nw, :vi, f_idx[2])
    vr_to   = var(pm, nw, :vr, t_idx[2])
    vi_to   = var(pm, nw, :vi, t_idx[2])

    for (idx, f_terminal, t_terminal, r) in zip(1:length(f_connections), f_connections, t_connections, c_rating)
        if r<Inf
            JuMP.@constraint(pm.model, p_fr[idx]^2+q_fr[idx]^2 <= r^2*(vr_fr[f_terminal]^2+vi_fr[f_terminal]^2))
            JuMP.@constraint(pm.model, p_to[idx]^2+q_to[idx]^2 <= r^2*(vr_to[t_terminal]^2+vi_to[t_terminal]^2))
        end
    end
end


# # BRANCH - Constraints - Reduced

# "Defines how current distributes over series and shunt impedances of a pi-model branch"
# function constraint_mc_current_from(pm::ReducedMultiConductorIVRModels, nw::Int, f_bus::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real}; report::Bool=true)
#     vr_fr = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
#     vi_fr = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]

#     csr_fr =  var(pm, nw, :csr, f_idx[1])
#     csi_fr =  var(pm, nw, :csi, f_idx[1])

#     cr_fr = csr_fr + g_sh_fr*vr_fr - b_sh_fr*vi_fr
#     ci_fr = csi_fr + g_sh_fr*vi_fr + b_sh_fr*vr_fr

#     var(pm, nw, :cr)[f_idx] = cr_fr
#     var(pm, nw, :ci)[f_idx] = ci_fr

#     var(pm, nw, :cr_bus)[f_idx] = cr_bus_fr = _merge_bus_flows(pm, cr_fr, f_connections)
#     var(pm, nw, :ci_bus)[f_idx] = ci_bus_fr = _merge_bus_flows(pm, ci_fr, f_connections)

#     if report
#         sol(pm, nw, :branch, f_idx[1])[:cr_fr] = cr_fr
#         sol(pm, nw, :branch, f_idx[1])[:ci_fr] = ci_fr
#     end
# end


# "Defines how current distributes over series and shunt impedances of a pi-model branch"
# function constraint_mc_current_to(pm::ReducedMultiConductorIVRModels, nw::Int, t_bus, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, g_sh_to::Matrix{<:Real}, b_sh_to::Matrix{<:Real}; report::Bool=true)
#     vr_to = [var(pm, nw, :vr, t_bus)[c] for c in t_connections]
#     vi_to = [var(pm, nw, :vi, t_bus)[c] for c in t_connections]

#     csr_to = -var(pm, nw, :csr, f_idx[1])
#     csi_to = -var(pm, nw, :csi, f_idx[1])

#     cr_to = csr_to + g_sh_to*vr_to - b_sh_to*vi_to
#     ci_to = csi_to + g_sh_to*vi_to + b_sh_to*vr_to

#     var(pm, nw, :cr)[t_idx] = cr_to
#     var(pm, nw, :ci)[t_idx] = ci_to

#     var(pm, nw, :cr_bus)[t_idx] = cr_bus_to = _merge_bus_flows(pm, cr_to, t_connections)
#     var(pm, nw, :ci_bus)[t_idx] = ci_bus_to = _merge_bus_flows(pm, ci_to, t_connections)

#     if report
#         sol(pm, nw, :branch, f_idx[1])[:cr_to] = cr_to
#         sol(pm, nw, :branch, f_idx[1])[:ci_to] = ci_to
#     end
# end


# ## SHARED 

# ""
# function _merge_bus_flows(pm::AbstractMultiConductorIVRModel, flows::Vector, connections::Vector)::JuMP.Containers.DenseAxisArray
#     flows_merged = []
#     conns_unique = unique(connections)
#     for t in conns_unique
#         idxs = findall(connections.==t)
#         flows_t = flows[idxs]
#         if length(flows_t)==1
#             flows_merged_t = flows_t[1]
#         elseif any(isa(a, JuMP.NonlinearExpression) for a in flows_t)
#             flows_merged_t = JuMP.@NLexpression(pm.model, sum(flows_t[i] for i in 1:length(flows_t)))
#         else
#             flows_merged_t = sum(flows_t)
#         end
#         push!(flows_merged, flows_merged_t)
#     end
#     JuMP.Containers.DenseAxisArray(flows_merged, conns_unique)
# end


# # STORAGE

# # STORAGE - Variable

# "variables for modeling storage units, includes grid injection and internal variables"
# function variable_mc_storage(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, kwargs...)
#     variable_mc_storage_power_real(pm; bounded=bounded, kwargs...)
#     variable_mc_storage_power_imaginary(pm; bounded=bounded, kwargs...)
#     variable_mc_storage_power_control_imaginary(pm; bounded=false, kwargs...)
#     var(pm, nw)[:crs_bus] = Dict{Int,Any}()
#     var(pm, nw)[:cis_bus] = Dict{Int,Any}()
#     variable_mc_storage_current_real(pm; bounded=bounded, kwargs...)
#     variable_mc_storage_current_imaginary(pm; bounded=bounded, kwargs...)
#     variable_storage_energy(pm; bounded=bounded, kwargs...)
#     variable_storage_charge(pm; bounded=bounded, kwargs...)
#     variable_storage_discharge(pm; bounded=bounded, kwargs...)
# end


# ""
# function variable_mc_storage_power_real(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(i =>_infer_int_dim_unit(strg, false) for (i,strg) in ref(pm, nw, :storage))
#     ps = var(pm, nw)[:ps] = Dict(i => JuMP.@variable(pm.model,
#             [c in 1:int_dim[i]], base_name="$(nw)_ps_$(i)",
#             start = comp_start_value(ref(pm, nw, :storage, i), "ps_start", c, 0.0)
#         ) for i in ids(pm, nw, :storage)
#     )

#     report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :ps, ids(pm, nw, :storage), ps)
# end


# ""
# function variable_mc_storage_power_imaginary(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(i =>_infer_int_dim_unit(strg, false) for (i,strg) in ref(pm, nw, :storage))
#     qs = var(pm, nw)[:qs] = Dict(i => JuMP.@variable(pm.model,
#             [c in int_dim[i]], base_name="$(nw)_qs_$(i)",
#             start = comp_start_value(ref(pm, nw, :storage, i), "qs_start", c, 0.0)
#         ) for i in ids(pm, nw, :storage)
#     )

#     report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :qs, ids(pm, nw, :storage), qs)
# end


# ""
# function variable_mc_storage_current_real(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(i =>_infer_int_dim_unit(strg, false) for (i,strg) in ref(pm, nw, :storage))
#     crs = var(pm, nw)[:crs] = Dict(i => JuMP.@variable(pm.model,
#             [c in 1:int_dim[i]], base_name="$(nw)_crs_$(i)",
#             start = comp_start_value(ref(pm, nw, :storage, i), "crs_start", c, 0.0)
#         ) for i in ids(pm, nw, :storage)
#     )

#     report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :crs, ids(pm, nw, :storage), crs)
# end


# ""
# function variable_mc_storage_current_imaginary(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     int_dim = Dict(i =>_infer_int_dim_unit(strg, false) for (i,strg) in ref(pm, nw, :storage))
#     cis = var(pm, nw)[:cis] = Dict(i => JuMP.@variable(pm.model,
#             [c in 1:int_dim[i]], base_name="$(nw)_cis_$(i)",
#             start = comp_start_value(ref(pm, nw, :storage, i), "cis_start", c, 0.0)
#         ) for i in ids(pm, nw, :storage)
#     )

#     report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :cis, ids(pm, nw, :storage), cis)
# end


# # STORAGE - Constraints

# function constraint_mc_storage_current(pm::AbstractMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)
#     strg = ref(pm, nw, :storage, id)

#     nphases = _infer_int_dim_unit(strg, false)
#     if get(strg, "configuration", WYE) == WYE || nphases==1
#         constraint_mc_storage_current_wye(pm, nw, id, strg["connections"]; report=report, bounded=bounded)
#     else
#         constraint_mc_storage_current_delta(pm, nw, id, strg["connections"]; report=report, bounded=bounded)
#     end
# end


# ""
# function constraint_mc_storage_current_wye(pm::AbstractMultiConductorIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
#     crs = var(pm, nw, :crs, id)
#     cis = var(pm, nw, :cis, id)

#     var(pm, nw, :crs_bus)[id] = _merge_bus_flows(pm, [crs..., -sum(crs)], connections)
#     var(pm, nw, :cis_bus)[id] = _merge_bus_flows(pm, [cis..., -sum(cis)], connections)
# end


# ""
# function constraint_mc_storage_current_delta(pm::AbstractMultiConductorIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
#     crs = var(pm, nw, :crs, id)
#     cis = var(pm, nw, :cis, id)
#     Md = _get_delta_transformation_matrix(length(connections))
#     var(pm, nw, :crs_bus)[id] = _merge_bus_flows(pm, Md'*crs, connections)
#     var(pm, nw, :cis_bus)[id] = _merge_bus_flows(pm, Md'*cis, connections)
# end


# ""
# function constraint_mc_storage_losses(pm::AbstractMultiConductorIVRModel, i::Int; nw::Int=nw_id_default, kwargs...)
#     strg = ref(pm, nw, :storage, i)

#     vr  = var(pm, nw,  :vr, strg["storage_bus"])
#     vi  = var(pm, nw,  :vi, strg["storage_bus"])
#     ps  = var(pm, nw,  :ps, i)
#     qs  = var(pm, nw,  :qs, i)
#     crs = var(pm, nw,  :crs, i)
#     cis = var(pm, nw,  :cis, i)
#     sc  = var(pm, nw,  :sc, i)
#     sd  = var(pm, nw,  :sd, i)
#     qsc = var(pm, nw, :qsc, i)

#     p_loss = strg["p_loss"]
#     q_loss = strg["q_loss"]
#     r = strg["r"]
#     x = strg["x"]

#     JuMP.@constraint(pm.model,
#         sum(ps) + (sd - sc)
#         ==
#         p_loss + sum(r .* (crs.^2 .+ cis.^2))
#     )

#     JuMP.@constraint(pm.model,
#         sum(qs)
#         ==
#         qsc + q_loss + sum(x .* (crs.^2 .+ cis.^2))
#     )
# end


# ""
# function constraint_mc_storage_thermal_limit(pm::AbstractMultiConductorIVRModel, nw::Int, i::Int, connections::Vector{Int}, rating::Vector{<:Real})
#     crs = var(pm, nw, :crs, i)
#     cis = var(pm, nw, :cis, i)

    
#     #TODO is this a current or power bound?
#     JuMP.@constraint(pm.model, crs.^2 + cis.^2 .<= rating.^2)
# end


# ""
# function constraint_mc_storage_power(pm::AbstractUnbalancedPowerModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)
#     strg = ref(pm, nw, :gen, id)
#     bus = ref(pm, nw, :bus, strg["storage_bus"])

#     if get(strg, "configuration", WYE) == WYE
#         constraint_mc_storage_power_wye(pm, nw, id, bus["index"], strg["connections"]; report=report, bounded=bounded)
#     else
#         constraint_mc_storage_power_delta(pm, nw, id, bus["index"], strg["connections"]; report=report, bounded=bounded)
#     end
# end


# ""
# function constraint_mc_storage_power_wye(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
#     vr = var(pm, nw, :vr, bus_id)
#     vi = var(pm, nw, :vi, bus_id)
#     crs = var(pm, nw, :crs, id)
#     cis = var(pm, nw, :cis, id)

#     phases = connections[1:end-1]
#     n      = connections[end]

#     vr_pn = [vr[p]-vr[n] for p in phases]
#     vi_pn = [vi[p]-vi[n] for p in phases]

#     ps = var(pm, nw, :ps, id)
#     qs = var(pm, nw, :qs, id)

#     JuMP.@constraint(pm.model, ps .==  vr_pn.*crs .+ vi_pn.*cis)
#     JuMP.@constraint(pm.model, qs .== -vr_pn.*cis .+ vi_pn.*crs)
# end


# ""
# function constraint_mc_storage_power_delta(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
#     vr = var(pm, nw, :vr, bus_id)
#     vi = var(pm, nw, :vi, bus_id)
#     crs = var(pm, nw, :crs, id)
#     cis = var(pm, nw, :cis, id)

#     nph = length(pmin)

#     prev = Dict(c=>connections[(idx+nph-2)%nph+1] for (idx,c) in enumerate(connections))
#     next = Dict(c=>connections[idx%nph+1] for (idx,c) in enumerate(connections))

#     vrs = [vr[c]-vr[next[c]] for p in connections]
#     vis = [vi[c]-vi[next[c]] for p in connections]

#     ps = var(pm, nw, :ps, id)
#     qs = var(pm, nw, :qs, id)

#     JuMP.@constraint(pm.model, ps .==  vrs.*crs .+ vis.*cis)
#     JuMP.@constraint(pm.model, qs .== -vrs.*cis .+ vis.*crs)
# end

