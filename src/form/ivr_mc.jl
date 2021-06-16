# TODO 
# - reorganize load constraints across power and current
# - add switch current constraints
# - rating to limit rename

# BUS

# Variables

""
function variable_mc_bus_voltage(pm::AbstractMultiConductorIVRModel; nw=nw_id_default, bounded::Bool=true, kwargs...)
    variable_mc_bus_voltage_real(pm; nw=nw, bounded=bounded, kwargs...)
    variable_mc_bus_voltage_imaginary(pm; nw=nw, bounded=bounded, kwargs...)
end


""
function variable_mc_bus_voltage_real(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref(pm, nw, :bus))

    vr = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_vr_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "vr_start", t, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in ref(pm, nw, :bus)
            for (idx,t) in enumerate(terminals[i])
                set_lower_bound(vr[i][t], -bus["vmax"][idx])
                set_upper_bound(vr[i][t],  bus["vmax"][idx])
            end
        end
    end

    vr = var(pm, nw)[:vr] = Dict(i=>JuMP.Containers.DenseAxisArray(
        Vector{JuMP.AffExpr}([t in vr[i].axes[1] ? vr[i][t] : 0.0 for t in bus["terminals"]]), bus["terminals"]
        ) for (i, bus) in ref(pm, nw, :bus))

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :vr, ids(pm, nw, :bus), vr)
end


""
function variable_mc_bus_voltage_imaginary(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i,bus) in ref(pm, nw, :bus))
    vi = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_vi_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "vi_start", t, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in ref(pm, nw, :bus)
            for (idx,t) in enumerate(terminals[i])
                set_lower_bound(vi[i][t], -bus["vmax"][idx])
                set_upper_bound(vi[i][t],  bus["vmax"][idx])
            end
        end
    end

    
    vi = var(pm, nw)[:vi] = Dict(i=>JuMP.Containers.DenseAxisArray(
        Vector{JuMP.AffExpr}([t in vi[i].axes[1] ? vi[i][t] : 0.0 for t in bus["terminals"]]), bus["terminals"]
        ) for (i, bus) in ref(pm, nw, :bus))

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :vi, ids(pm, nw, :bus), vi)
end


# Constraints

"""
Kirchhoff's current law applied to buses
`sum(cr + im*ci) = 0`
"""
function constraint_mc_current_balance(pm::AbstractMultiConductorIVRModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    cr    = get(var(pm, nw),    :cr_bus, Dict()); _check_var_keys(cr, bus_arcs, "real current", "branch")
    ci    = get(var(pm, nw),    :ci_bus, Dict()); _check_var_keys(ci, bus_arcs, "imaginary current", "branch")
    crd   = get(var(pm, nw),   :crd_bus, Dict()); _check_var_keys(crd, bus_loads, "real current", "load")
    cid   = get(var(pm, nw),   :cid_bus, Dict()); _check_var_keys(cid, bus_loads, "imaginary current", "load")
    crg   = get(var(pm, nw),   :crg_bus, Dict()); _check_var_keys(crg, bus_gens, "real current", "generator")
    cig   = get(var(pm, nw),   :cig_bus, Dict()); _check_var_keys(cig, bus_gens, "imaginary current", "generator")
    crs   = get(var(pm, nw),   :crs_bus, Dict()); _check_var_keys(crs, bus_storage, "real currentr", "storage")
    cis   = get(var(pm, nw),   :cis_bus, Dict()); _check_var_keys(cis, bus_storage, "imaginary current", "storage")
    crsw  = get(var(pm, nw),  :crsw, Dict()); _check_var_keys(crsw, bus_arcs_sw, "real current", "switch")
    cisw  = get(var(pm, nw),  :cisw, Dict()); _check_var_keys(cisw, bus_arcs_sw, "imaginary current", "switch")
    crt   = get(var(pm, nw),   :crt_bus, Dict()); _check_var_keys(crt, bus_arcs_trans, "real current", "transformer")
    cit   = get(var(pm, nw),   :cit_bus, Dict()); _check_var_keys(cit, bus_arcs_trans, "imaginary current", "transformer")

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx, t) in ungrounded_terminals
        # @smart_constraint(pm.model,  [cr, crd, crg, crs, crsw, crt, vr],
        @smart_constraint(pm.model,  [cr, crd, crg, crt, vr],
                                      sum(cr[a][t] for (a, conns) in bus_arcs if t in conns)
                                    + sum(crsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                    + sum(crt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                    ==
                                      sum(crg[g][t]         for (g, conns) in bus_gens if t in conns)
                                    - sum(crs[s][t]         for (s, conns) in bus_storage if t in conns)
                                    - sum(crd[d][t]         for (d, conns) in bus_loads if t in conns)
                                    - sum( Gt[idx,jdx]*vr[u] -Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals) # shunts
                                    )
        # @smart_constraint(pm.model, [ci, cid, cig, cis, cisw, cit, vi],
        @smart_constraint(pm.model, [ci, cid, cig, cit, vi],
                                      sum(ci[a][t] for (a, conns) in bus_arcs if t in conns)
                                    + sum(cisw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                    + sum(cit[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                    ==
                                      sum(cig[g][t]         for (g, conns) in bus_gens if t in conns)
                                    - sum(cis[s][t]         for (s, conns) in bus_storage if t in conns)
                                    - sum(cid[d][t]         for (d, conns) in bus_loads if t in conns)
                                    - sum( Gt[idx,jdx]*vi[u] +Bt[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals) # shunts
                                    )
    end
end


""
function constraint_mc_voltage_absolute(pm::AbstractMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    bus = ref(pm, nw, :bus, id)

    constraint_mc_voltage_absolute(pm, nw, id, bus["terminals"], bus["grounded"], bus["vmin"], bus["vmax"])
end


""
function constraint_mc_voltage_absolute(pm::AbstractMultiConductorIVRModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, vmin::Vector{<:Real}, vmax::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    ungrounded_terminals = terminals[(!).(grounded)]
    JuMP.@constraint(pm.model, [t in ungrounded_terminals], vmin[t]^2 <= vr[t]^2+vi[t]^2 <= vmax[t]^2)
end


""
function constraint_mc_voltage_pairwise(pm::AbstractMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    bus = ref(pm, nw, :bus, id)

    vm_pair_lb = bus["vm_pair_lb"]
    vm_pair_ub = bus["vm_pair_ub"]

    constraint_mc_voltage_pairwise(pm, nw, id, bus["terminals"], bus["grounded"], vm_pair_lb, vm_pair_ub)
end


""
function constraint_mc_voltage_pairwise(pm::AbstractMultiConductorIVRModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, vm_pair_lb::Vector, vm_pair_ub::Vector; report::Bool=true)
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    for (a,b,lb) in vm_pair_lb
        JuMP.@constraint(pm.model, (vr[a]-vr[b])^2 + (vi[a]-vi[b])^2 <= lb^2)
    end

    for (a,b,ub) in vm_pair_ub
        JuMP.@constraint(pm.model, (vr[a]-vr[b])^2 + (vi[a]-vi[b])^2 <= ub^2)
    end
end

# GENERATOR

# GENERATOR - Variables

""
function variable_mc_generator_current_real(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(gen, false) for (i,gen) in ref(pm, nw, :gen))
    crg = var(pm, nw)[:crg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_crg_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "crg_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )
    # if bounded
    #     for (i, g) in ref(pm, nw, :gen) if haskey(g, "c_rating")
    #         cmax = c["c_rating"]
    #         set_lower_bound.(crg[i], -cmax)
    #         set_upper_bound.(crg[i],  cmax)
    #     end
    # end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :crg, ids(pm, nw, :gen), crg)
end


""
function variable_mc_generator_current_imaginary(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(gen, false) for (i,gen) in ref(pm, nw, :gen))
    cig = var(pm, nw)[:cig] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_cig_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "cig_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )
    # if bounded
    #     for (i, g) in ref(pm, nw, :gen) if haskey(g, "c_rating")
    #         cmax = c["c_rating"]
    #         set_lower_bound.(cig[i], -cmax)
    #         set_upper_bound.(cig[i],  cmax)
    #     end
    # end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :cig, ids(pm, nw, :gen), cig)
end


""
function variable_mc_generator_power_real(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(gen, false) for (i,gen) in ref(pm, :gen))
    pg = var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_pg_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "pg_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )

    if bounded
        for (i,gen) in ref(pm, nw, :gen)
            set_lower_bound.(pg[i], gen["pmin"])
            set_upper_bound.(pg[i], gen["pmax"])
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :pg, ids(pm, nw, :gen), pg)
end


""
function variable_mc_generator_power_imaginary(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(gen, false) for (i,gen) in ref(pm, :gen))
    qg = var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_qg_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "qg_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )

    if bounded
        for (i,gen) in ref(pm, nw, :gen)
            set_lower_bound.(qg[i], gen["qmin"])
            set_upper_bound.(qg[i], gen["qmax"])
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :qg, ids(pm, nw, :gen), qg)
end


# GENERATOR - Variables - Non-linear

""
function variable_mc_generator_power(pm::AbstractNLMultiConductorIVRModel; kwargs...)
    #TODO docs
end

# GENERATOR - Variables - Quadratic

""
function variable_mc_generator_power(pm::AbstractQuadraticMultiConductorIVRModel; kwargs...)
    variable_mc_generator_power_real(pm; kwargs...)
    variable_mc_generator_power_imaginary(pm; kwargs...)
end


# GENERATOR - Constraints


function constraint_mc_generator_current(pm::AbstractMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)
    generator = ref(pm, nw, :gen, id)

    nphases = _infer_int_dim_unit(generator, false)
    if get(generator, "configuration", WYE) == WYE || nphases==1
        constraint_mc_generator_current_wye(pm, nw, id, generator["connections"]; report=report, bounded=bounded)
    else
        constraint_mc_generator_current_delta(pm, nw, id, generator["connections"]; report=report, bounded=bounded)
    end
end


""
function constraint_mc_generator_current_wye(pm::AbstractMultiConductorIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)
    var(pm, nw, :crg_bus)[id] = _merge_bus_flows(pm, [crg..., -sum(crg)], connections)
    var(pm, nw, :cig_bus)[id] = _merge_bus_flows(pm, [cig..., -sum(cig)], connections)
end


""
function constraint_mc_generator_current_delta(pm::AbstractMultiConductorIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)
    Md = _get_delta_transformation_matrix(length(connections))
    var(pm, nw, :crg_bus)[id] = _merge_bus_flows(pm, Md'*crg, connections)
    var(pm, nw, :cig_bus)[id] = _merge_bus_flows(pm, Md'*cig, connections)
end


# GENERATOR - Constraints - Non-linear

""
function constraint_mc_generator_power_wye(pm::AbstractNLMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)

    phases = connections[1:end-1]
    n      = connections[end]

    pg = Vector{JuMP.NonlinearExpression}([])
    qg = Vector{JuMP.NonlinearExpression}([])

    for (idx, p) in enumerate(phases)
        push!(pg, JuMP.@NLexpression(pm.model,  (vr[p]-vr[n])*crg[idx]+(vi[p]-vi[n])*cig[idx]))
        push!(qg, JuMP.@NLexpression(pm.model, -(vr[p]-vr[n])*cig[idx]+(vi[p]-vi[n])*crg[idx]))
    end

    if bounded
        for (idx, p) in enumerate(phases)
            if pmin[idx]>-Inf
                JuMP.@constraint(pm.model, pmin[idx] .<= (vr[p]-vr[n])*crg[idx]  + (vi[p]-vi[n])*cig[idx])
            end
            if pmax[idx]< Inf
                JuMP.@constraint(pm.model, pmax[idx] .>= (vr[p]-vr[n])*crg[idx]  + (vi[p]-vi[n])*cig[idx])
            end
            if qmin[idx]>-Inf
                JuMP.@constraint(pm.model, qmin[idx] .<= (vi[p]-vi[n])*crg[idx]  - (vr[p]-vr[n])*cig[idx])
            end
            if qmax[idx]< Inf
                JuMP.@constraint(pm.model, qmax[idx] .>= (vi[p]-vi[n])*crg[idx]  - (vr[p]-vr[n])*cig[idx])
            end
        end
    end

    var(pm, nw, :pg)[id] = pg
    var(pm, nw, :qg)[id] = qg

    if report
        sol(pm, nw, :gen, id)[:pg] = pg
        sol(pm, nw, :gen, id)[:qg] = qg
    end
end


""
function constraint_mc_generator_power_delta(pm::AbstractNLMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)

    nph = length(pmin)

    prev = Dict(c=>connections[(idx+nph-2)%nph+1] for (idx,c) in enumerate(connections))
    next = Dict(c=>connections[idx%nph+1] for (idx,c) in enumerate(connections))

    vrg = Dict()
    vig = Dict()
    for c in connections
        vrg[c] = JuMP.@NLexpression(pm.model, vr[c]-vr[next[c]])
        vig[c] = JuMP.@NLexpression(pm.model, vi[c]-vi[next[c]])
    end

    pg = Vector{JuMP.NonlinearExpression}([])
    qg = Vector{JuMP.NonlinearExpression}([])
    for (idx,c) in enumerate(connections)
        push!(pg, JuMP.@NLexpression(pm.model,  vrg[c]*crg[idx]+vig[c]*cig[idx]))
        push!(qg, JuMP.@NLexpression(pm.model, -vrg[c]*cig[idx]+vig[c]*crg[idx]))
    end

    if bounded
        JuMP.@NLconstraint(pm.model, [i in 1:nph], pmin[i] <= pg[i])
        JuMP.@NLconstraint(pm.model, [i in 1:nph], pmax[i] >= pg[i])
        JuMP.@NLconstraint(pm.model, [i in 1:nph], qmin[i] <= qg[i])
        JuMP.@NLconstraint(pm.model, [i in 1:nph], qmax[i] >= qg[i])
    end

    var(pm, nw, :pg)[id] = JuMP.Containers.DenseAxisArray(pg, connections)
    var(pm, nw, :qg)[id] = JuMP.Containers.DenseAxisArray(qg, connections)

    if report
        sol(pm, nw, :gen, id)[:pg] = pg
        sol(pm, nw, :gen, id)[:qg] = qg
    end
end


# GENERATOR - Constraints - Quadratic

""
function constraint_mc_generator_power_wye(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)

    phases = connections[1:end-1]
    n      = connections[end]

    vr_pn = [vr[p]-vr[n] for p in phases]
    vi_pn = [vi[p]-vi[n] for p in phases]

    pg = var(pm, nw, :pg, id)
    qg = var(pm, nw, :qg, id)

    JuMP.@constraint(pm.model, pg .==  vr_pn.*crg .+ vi_pn.*cig)
    JuMP.@constraint(pm.model, qg .== -vr_pn.*cig .+ vi_pn.*crg)
end


""
function constraint_mc_generator_power_delta(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)

    nph = length(pmin)

    prev = Dict(c=>connections[(idx+nph-2)%nph+1] for (idx,c) in enumerate(connections))
    next = Dict(c=>connections[idx%nph+1] for (idx,c) in enumerate(connections))

    vrg = [vr[c]-vr[next[c]] for p in connections]
    vig = [vi[c]-vi[next[c]] for p in connections]

    pg = var(pm, nw, :pg, id)
    qg = var(pm, nw, :qg, id)

    JuMP.@constraint(pm.model, pg .==  vrg.*crg .+ vig.*cig)
    JuMP.@constraint(pm.model, qg .== -vrg.*cig .+ vig.*crg)
end


# LOAD

# LOAD - Variables

""
function variable_mc_load_current(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    var(pm, nw)[:crd] = Dict{Int, Any}()
    var(pm, nw)[:cid] = Dict{Int, Any}()
    var(pm, nw)[:crd_bus] = Dict{Int, Any}()
    var(pm, nw)[:cid_bus] = Dict{Int, Any}()
end


""
function variable_mc_load_power(pm::AbstractNLMultiConductorIVRModel; nw=nw_id_default, bounded::Bool=true, report::Bool=true)
    var(pm, nw)[:pd] = Dict{Int, Any}()
    var(pm, nw)[:qd] = Dict{Int, Any}()
    var(pm, nw)[:pd_bus] = Dict{Int, Any}()
    var(pm, nw)[:qd_bus] = Dict{Int, Any}()
end


""
function variable_mc_load_current(pm::AbstractQuadraticMultiConductorIVRModel; nw::Int=nw_id_default, kwargs...)
    load_ids_exponential = [id for (id,load) in ref(pm, nw, :load) if load["model"]==EXPONENTIAL]
    @assert isempty(load_ids_exponential) "Exponential loads cannot be represented quadratically."

    variable_mc_load_current_real(pm; nw=nw, kwargs...)
    variable_mc_load_current_imaginary(pm; nw=nw, kwargs...)

    var(pm, nw)[:crd_bus] = Dict{Int,Any}()
    var(pm, nw)[:cid_bus] = Dict{Int,Any}()
end


""
function variable_mc_load_current_real(pm::AbstractQuadraticMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(load, false) for (i,load) in ref(pm, :load))

    crd = var(pm, nw)[:crd] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_crd_$(i)",
            start = comp_start_value(ref(pm, nw, :load, i), "crd_start", c, 0.0)
        ) for i in ids(pm, nw, :load)
    )

    # if bounded
    #     for (i,load) in ref(pm, nw, :gen)
    #         if haskey(gen, "pmin")
    #             set_lower_bound.(pg[i], gen["pmin"])
    #         end
    #         if haskey(gen, "pmax")
    #             set_upper_bound.(pg[i], gen["pmax"])
    #         end
    #     end
    # end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :crd, ids(pm, nw, :load), crd)
end


""
function variable_mc_load_current_imaginary(pm::AbstractQuadraticMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(load, false) for (i,load) in ref(pm, :load))
    load_ids_current = [id for (id,load) in ref(pm, nw, :load) if load["model"]==CURRENT]
    
    cid = var(pm, nw)[:cid] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_cid_$(i)",
            start = comp_start_value(ref(pm, nw, :load, i), "cid_start", c, 0.0)
        ) for i in ids(pm, nw, :load)
    )

    # if bounded
    #     for (i,gen) in ref(pm, nw, :gen)
    #         if haskey(gen, "qmin")
    #             set_lower_bound.(qg[i], gen["qmin"])
    #         end
    #         if haskey(gen, "qmax")
    #             set_upper_bound.(qg[i], gen["qmax"])
    #         end
    #     end
    # end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :cid, ids(pm, nw, :load), cid)
end


# LOAD - Variables - Quadratic

""
function variable_mc_load_power(pm::AbstractQuadraticMultiConductorIVRModel; nw::Int=nw_id_default, kwargs...)
    load_ids_exponential = [id for (id,load) in ref(pm, nw, :load) if load["model"]==EXPONENTIAL]
    @assert isempty(load_ids_exponential) "Exponential loads cannot be represented quadratically."

    variable_mc_load_power_real(pm; nw=nw, kwargs...)
    variable_mc_load_power_imaginary(pm; nw=nw, kwargs...)
end


""
function variable_mc_load_power_real(pm::AbstractQuadraticMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(load, false) for (i,load) in ref(pm, :load))
    load_ids_current = [id for (id,load) in ref(pm, nw, :load) if load["model"]==CURRENT]
    
    pd = var(pm, nw)[:pd] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_pd_$(i)",
            start = comp_start_value(ref(pm, nw, :load, i), "pd_start", c, 0.0)
        ) for i in load_ids_current
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :pd, load_ids_current, pd)
end


""
function variable_mc_load_power_imaginary(pm::AbstractQuadraticMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(load, false) for (i,load) in ref(pm, :load))
    load_ids_current = [id for (id,load) in ref(pm, nw, :load) if load["model"]==CURRENT]
    
    qd = var(pm, nw)[:qd] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_qd_$(i)",
            start = comp_start_value(ref(pm, nw, :load, i), "qd_start", c, 0.0)
        ) for i in load_ids_current
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :qd, load_ids_current, qd)
end


# LOAD - Constraints

""
function constraint_mc_load_power(pm::AbstractMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
    # nothing to do
end

""
function constraint_mc_load_current(pm::AbstractMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
    load = ref(pm, nw, :load, id)
    bus = ref(pm, nw,:bus, load["load_bus"])

    configuration = load["configuration"]

    a, alpha, b, beta = _load_expmodel_params(load, bus)

    int_dim = _infer_int_dim_unit(load, false)
    if configuration==WYE || int_dim==1
        constraint_mc_load_current_wye(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    else
        constraint_mc_load_current_delta(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    end
end


""
function constraint_mc_load_current_wye(pm::AbstractMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    crd = Vector{JuMP.NonlinearExpression}([])
    cid = Vector{JuMP.NonlinearExpression}([])

    phases = connections[1:end-1]
    n      = connections[end]

    for (idx, p) in enumerate(phases)
        push!(crd, JuMP.@NLexpression(pm.model,
             a[idx]*(vr[p]-vr[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(alpha[idx]/2-1)
            +b[idx]*(vi[p]-vi[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(beta[idx]/2 -1)
        ))
        push!(cid, JuMP.@NLexpression(pm.model,
             a[idx]*(vi[p]-vi[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(alpha[idx]/2-1)
            -b[idx]*(vr[p]-vr[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(beta[idx]/2 -1)
        ))
    end

    var(pm, nw, :crd)[id] = crd
    var(pm, nw, :cid)[id] = cid

    crd_bus_n = JuMP.@NLexpression(pm.model, -sum(crd[i] for i in 1:length(phases)))
    cid_bus_n = JuMP.@NLexpression(pm.model, -sum(cid[i] for i in 1:length(phases)))

    var(pm, nw, :crd_bus)[id] = crd_bus = _merge_bus_flows(pm, [crd..., crd_bus_n], connections)
    var(pm, nw, :cid_bus)[id] = cid_bus = _merge_bus_flows(pm, [cid..., cid_bus_n], connections)

    if report
        pd_bus = Vector{JuMP.NonlinearExpression}([])
        qd_bus = Vector{JuMP.NonlinearExpression}([])
        for (idx,c) in enumerate(connections)
            push!(pd_bus, JuMP.@NLexpression(pm.model,  vr[c]*crd_bus[c]+vi[c]*cid_bus[c]))
            push!(qd_bus, JuMP.@NLexpression(pm.model, -vr[c]*cid_bus[c]+vi[c]*crd_bus[c]))
        end

        sol(pm, nw, :load, id)[:pd_bus] = JuMP.Containers.DenseAxisArray(pd_bus, connections)
        sol(pm, nw, :load, id)[:qd_bus] = JuMP.Containers.DenseAxisArray(qd_bus, connections)

        sol(pm, nw, :load, id)[:crd] = JuMP.Containers.DenseAxisArray(crd, connections)
        sol(pm, nw, :load, id)[:cid] = JuMP.Containers.DenseAxisArray(cid, connections)

        sol(pm, nw, :load, id)[:crd_bus] = crd_bus
        sol(pm, nw, :load, id)[:cid_bus] = cid_bus

        pd = Vector{JuMP.NonlinearExpression}([])
        qd = Vector{JuMP.NonlinearExpression}([])
        for (idx, p) in enumerate(phases)
            push!(pd, JuMP.@NLexpression(pm.model, a[idx]*(vr[p]^2+vi[p]^2)^(alpha[idx]/2) ))
            push!(qd, JuMP.@NLexpression(pm.model, b[idx]*(vr[p]^2+vi[p]^2)^(beta[idx]/2)  ))
        end
        sol(pm, nw, :load, id)[:pd] = JuMP.Containers.DenseAxisArray(pd, connections)
        sol(pm, nw, :load, id)[:qd] = JuMP.Containers.DenseAxisArray(qd, connections)
    end
end


""
function constraint_mc_load_current_delta(pm::AbstractMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    nph = length(connections)
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

    var(pm, nw, :crd_bus)[id] = _merge_bus_flows(pm, crd_bus, connections)
    var(pm, nw, :cid_bus)[id] = _merge_bus_flows(pm, cid_bus, connections)

    if report
        pd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph],  vr[i]*crd_bus[i]+vi[i]*cid_bus[i])
        qd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], -vr[i]*cid_bus[i]+vi[i]*crd_bus[i])

        sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        pd = JuMP.@NLexpression(pm.model, [i in 1:nph], a[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2) )
        qd = JuMP.@NLexpression(pm.model, [i in 1:nph], b[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2)  )
        sol(pm, nw, :load, id)[:pd] = pd
        sol(pm, nw, :load, id)[:qd] = qd
    end
end


# LOAD - Constraints - Quadratic

""
function constraint_mc_load_power(pm::AbstractQuadraticMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
    load = ref(pm, nw, :load, id)
    bus = ref(pm, nw,:bus, load["load_bus"])

    configuration = load["configuration"]
    int_dim = _infer_int_dim_unit(load, false)
    a, alpha, b, beta = _load_expmodel_params(load, bus)

    if configuration==WYE || int_dim==1
        constraint_mc_load_power_wye(pm, nw, id, load["load_bus"], load["connections"], load["model"], a, b; report=report)
    else
        constraint_mc_load_power_delta(pm, nw, id, load["load_bus"], load["connections"], load["model"], a, b; report=report)
    end
end


""
function constraint_mc_load_power_wye(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, model::LoadModel, a::Vector{<:Real}, b::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    crd = var(pm, nw, :crd, id)
    cid = var(pm, nw, :cid, id)

    phases = connections[1:end-1]
    n      = connections[end]

    vr_pn = [vr[p]-vr[n] for p in phases]
    vi_pn = [vi[p]-vi[n] for p in phases]

    if model==POWER
        pd = a
        qd = b
    elseif model==IMPEDANCE
        pd = a .* (vr_pn.^2 .+ vi_pn.^2)
        qd = b .* (vr_pn.^2 .+ vi_pn.^2)
    elseif model==CURRENT
        pd = var(pm, nw, :pd, id)
        qd = var(pm, nw, :qd, id)
        JuMP.@constraint(pm.model, pd.^2 .== a.^2 .* (vr_pn.^2 .+ vi_pn.^2))
        JuMP.@constraint(pm.model, sign.(a).*pd .>= 0) #TODO handle in bounds instead
        JuMP.@constraint(pm.model, qd.^2 .== b.^2 .* (vr_pn.^2 .+ vi_pn.^2))
        JuMP.@constraint(pm.model, sign.(b).*qd .>= 0) #TODO handle in bounds instead
    else
        error("Load model $model for load $id is not supported by this formulation.")
    end

    JuMP.@constraint(pm.model, pd .==  vr_pn.*crd .+ vi_pn.*cid)
    JuMP.@constraint(pm.model, qd .== -vr_pn.*cid .+ vi_pn.*crd)

    # constant current loads are already reported through variable function
    if report && model!=CURRENT
        sol(pm, nw, :load, id)[:pd] = pd
        sol(pm, nw, :load, id)[:qd] = qd
    end
end


""
function constraint_mc_load_power_delta(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, model::LoadModel, a::Vector{<:Real}, b::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    crd = var(pm, nw, :crd, id)
    cid = var(pm, nw, :cid, id)

    phases = connections[1:end-1]
    n      = connections[end]

    Md = _get_delta_transformation_matrix(length(connections))
    vrd = Md*[vr[p] for p in phases]
    vid = Md*[vi[p] for p in phases]

    if model==POWER
        pd = a
        qd = b
    elseif model==IMPEDANCE
        pd = a .* (vrd.^2 .+ vid.^2)
        qd = b .* (vrd.^2 .+ vid.^2)
    elseif model==CURRENT
        pd = var(pm, nw, :pd, id)
        qd = var(pm, nw, :qd, id)
        JuMP.@constraint(pm.model, pd.^2 .== a.^2 .* (vrd.^2 .+ vid.^2))
        JuMP.@constraint(pm.model, sign.(a).*pd .>= 0) #TODO handle in bounds instead
        JuMP.@constraint(pm.model, qd.^2 .== b.^2 .* (vrd.^2 .+ vid.^2))
        JuMP.@constraint(pm.model, sign.(b).*qd .>= 0) #TODO handle in bounds instead
    else
        error("Load model $model for load $id is not supported by this formulation.")
    end

    JuMP.@constraint(pm.model, pd .==  vrd.*crd .+ vid.*cid)
    JuMP.@constraint(pm.model, qd .== -vrd.*cid .+ vid.*crd)

    # constant current loads are already reported through variable function
    if report && model!=CURRENT
        sol(pm, nw, :load, id)[:pd] = pd
        sol(pm, nw, :load, id)[:qd] = qd
    end
end


function constraint_mc_load_current(pm::AbstractQuadraticMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)
    load = ref(pm, nw, :load, id)

    int_dim = _infer_int_dim_unit(load, false)
    if get(load, "configuration", WYE) == WYE || int_dim==1
        constraint_mc_load_current_wye(pm, nw, id, load["connections"]; report=report, bounded=bounded)
    else
        constraint_mc_load_current_delta(pm, nw, id, load["connections"]; report=report, bounded=bounded)
    end
end


""
function constraint_mc_load_current_wye(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crd = var(pm, nw, :crd, id)
    cid = var(pm, nw, :cid, id)
    var(pm, nw, :crd_bus)[id] = _merge_bus_flows(pm, [crd..., -sum(crd)], connections)
    var(pm, nw, :cid_bus)[id] = _merge_bus_flows(pm, [cid..., -sum(cid)], connections)
end


""
function constraint_mc_load_current_delta(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crd = var(pm, nw, :crd, id)
    cid = var(pm, nw, :cid, id)
    Md = _get_delta_transformation_matrix(length(connections))
    var(pm, nw, :crd_bus)[id] = _merge_bus_flows(pm, Md'*crd, connections)
    var(pm, nw, :cid_bus)[id] = _merge_bus_flows(pm, Md'*cid, connections)
end


# TRANSFORMER

# Variables

""
function variable_mc_transformer_current(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_transformer_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_transformer_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    var(pm)[:crt_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm)[:cit_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end


"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_transformer_current_real(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in ref(pm, :transformer))
    crt = var(pm, nw)[:crt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:int_dim[l]], base_name="$(nw)_crt_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :transformer, l), "crt_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
    )

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :crt_fr, :crt_to, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), crt)
end


"variable: `ci[l,i,j] ` for `(l,i,j)` in `arcs`"
function variable_mc_transformer_current_imaginary(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in ref(pm, :transformer))
    cit = var(pm, nw)[:cit] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:int_dim[l]], base_name="$(nw)_cit_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :transformer, l), "cit_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
    )

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :cit_fr, :cit_to, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), cit)
end


""
function variable_mc_transformer_power_active(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in ref(pm, :transformer))
    pt = var(pm, nw)[:pt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:int_dim[l]], base_name="$(nw)_pt_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :transformer, l), "pt_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
    )
    
    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_transformer_from)
            trans = ref(pm, nw, :transformer, l)
            f_bus = ref(pm, nw, :bus, i)
            t_bus = ref(pm, nw, :bus, j)
            sm_ub = trans["sm_ub"]
            set_lower_bound(pt[(l,i,j)], -sm_ub)
            set_upper_bound(pt[(l,i,j)],  sm_ub)
            set_lower_bound(pt[(l,j,i)], -sm_ub)
            set_upper_bound(pt[(l,j,i)],  sm_ub)
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :pt_fr, :pt_to, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), pt)
end


""
function variable_mc_transformer_power_reactive(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in ref(pm, :transformer))
    qt = var(pm, nw)[:qt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:int_dim[l]], base_name="$(nw)_qt_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :transformer, l), "qt_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
    )
    
    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_transformer_from)
            trans = ref(pm, nw, :transformer, l)
            f_bus = ref(pm, nw, :bus, i)
            t_bus = ref(pm, nw, :bus, j)
            sm_ub = trans["sm_ub"]
            set_lower_bound(qt[(l,i,j)], -sm_ub)
            set_upper_bound(qt[(l,i,j)],  sm_ub)
            set_lower_bound(qt[(l,j,i)], -sm_ub)
            set_upper_bound(qt[(l,j,i)],  sm_ub)
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :qt_fr, :qt_to, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), qt)
end


# TRANSFORMER -  Variable - Non-linear

""
function variable_mc_transformer_power(pm::AbstractNLMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    # do nothing    
end


# TRANSFORMER - Variable - Quadratic

""
function variable_mc_transformer_power(pm::AbstractQuadraticMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_transformer_power_active(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_transformer_power_reactive(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# Constraints

""
function constraint_mc_transformer_voltage(pm::AbstractMultiConductorIVRModel, i::Int; nw::Int=nw_id_default, fix_taps::Bool=true)
    # if ref(pm, nw_id_default, :conductors)!=3
    #     error("Transformers only work with networks with three conductors.")
    # end

    transformer = ref(pm, nw, :transformer, i)
    f_bus = transformer["f_bus"]
    t_bus = transformer["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    configuration = transformer["configuration"]
    f_connections = transformer["f_connections"]
    t_connections = transformer["t_connections"]
    tm_set = transformer["tm_set"]
    tm_fixed = fix_taps ? ones(Bool, length(tm_set)) : transformer["tm_fix"]
    tm_scale = calculate_tm_scale(transformer, ref(pm, nw, :bus, f_bus), ref(pm, nw, :bus, t_bus))

    #TODO change data model
    # there is redundancy in specifying polarity seperately on from and to side
    #TODO change this once migrated to new data model
    pol = transformer["polarity"]

    if configuration == WYE
        constraint_mc_transformer_voltage_yy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == DELTA
        constraint_mc_transformer_voltage_dy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == "zig-zag"
        error("Zig-zag not yet supported.")
    end
end


""
function constraint_mc_transformer_voltage_yy(pm::AbstractMultiConductorIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    vr_fr_P = [var(pm, nw, :vr, f_bus)[c] for c in f_connections[1:end-1]]
    vi_fr_P = [var(pm, nw, :vi, f_bus)[c] for c in f_connections[1:end-1]]
    vr_fr_n = var(pm, nw, :vr, f_bus)[f_connections[end]]
    vi_fr_n = var(pm, nw, :vi, f_bus)[f_connections[end]]
    vr_to_P = [var(pm, nw, :vr, t_bus)[c] for c in t_connections[1:end-1]]
    vi_to_P = [var(pm, nw, :vi, t_bus)[c] for c in t_connections[1:end-1]]
    vr_to_n = var(pm, nw, :vr, t_bus)[t_connections[end]]
    vi_to_n = var(pm, nw, :vi, t_bus)[t_connections[end]]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(f_connections)-1]
    scale = (tm_scale*pol).*tm_set

    JuMP.@constraint(pm.model, (vr_fr_P.-vr_fr_n) .== scale.*(vr_to_P.-vr_to_n))
    JuMP.@constraint(pm.model, (vi_fr_P.-vi_fr_n) .== scale.*(vi_to_P.-vi_to_n))
end


""
function constraint_mc_transformer_voltage_dy(pm::AbstractMultiConductorIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    vr_fr_P = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr_P = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]
    vr_to_P = [var(pm, nw, :vr, t_bus)[c] for c in t_connections[1:end-1]]
    vi_to_P = [var(pm, nw, :vi, t_bus)[c] for c in t_connections[1:end-1]]
    vr_to_n = var(pm, nw, :vr, t_bus)[t_connections[end]]
    vi_to_n = var(pm, nw, :vi, t_bus)[t_connections[end]]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))]
    scale = (tm_scale*pol).*tm_set

    n_phases = length(tm)
    Md = _get_delta_transformation_matrix(n_phases)

    JuMP.@constraint(pm.model, Md*vr_fr_P .== scale.*(vr_to_P .- vr_to_n))
    JuMP.@constraint(pm.model, Md*vi_fr_P .== scale.*(vi_to_P .- vi_to_n))
end


""
function constraint_mc_transformer_current(pm::AbstractMultiConductorIVRModel, i::Int; nw::Int=nw_id_default, fix_taps::Bool=true)
    # if ref(pm, nw_id_default, :conductors)!=3
    #     error("Transformers only work with networks with three conductors.")
    # end

    transformer = ref(pm, nw, :transformer, i)
    f_bus = transformer["f_bus"]
    t_bus = transformer["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    configuration = transformer["configuration"]
    f_connections = transformer["f_connections"]
    t_connections = transformer["t_connections"]
    tm_set = transformer["tm_set"]
    tm_fixed = fix_taps ? ones(Bool, length(tm_set)) : transformer["tm_fix"]
    tm_scale = calculate_tm_scale(transformer, ref(pm, nw, :bus, f_bus), ref(pm, nw, :bus, t_bus))

    #TODO change data model
    # there is redundancy in specifying polarity seperately on from and to side
    #TODO change this once migrated to new data model
    pol = transformer["polarity"]

    if configuration == WYE
        constraint_mc_transformer_current_yy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == DELTA
        constraint_mc_transformer_current_dy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == "zig-zag"
        error("Zig-zag not yet supported.")
    end
end


""
function constraint_mc_transformer_current_yy(pm::AbstractMultiConductorIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    cr_fr_P = var(pm, nw, :crt, f_idx)
    ci_fr_P = var(pm, nw, :cit, f_idx)
    cr_to_P = var(pm, nw, :crt, t_idx)
    ci_to_P = var(pm, nw, :cit, t_idx)

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(f_connections)-1]
    scale = (tm_scale*pol).*tm_set

    JuMP.@constraint(pm.model, scale.*cr_fr_P .+ cr_to_P .== 0)
    JuMP.@constraint(pm.model, scale.*ci_fr_P .+ ci_to_P .== 0)

    var(pm, nw, :crt_bus)[f_idx] = _merge_bus_flows(pm, [cr_fr_P..., -sum(cr_fr_P)], f_connections)
    var(pm, nw, :cit_bus)[f_idx] = _merge_bus_flows(pm, [ci_fr_P..., -sum(ci_fr_P)], f_connections)
    var(pm, nw, :crt_bus)[t_idx] = _merge_bus_flows(pm, [cr_to_P..., -sum(cr_to_P)], t_connections)
    var(pm, nw, :cit_bus)[t_idx] = _merge_bus_flows(pm, [ci_to_P..., -sum(ci_to_P)], t_connections)
end


""
function constraint_mc_transformer_current_dy(pm::AbstractMultiConductorIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    cr_fr_P = var(pm, nw, :crt, f_idx)
    ci_fr_P = var(pm, nw, :cit, f_idx)
    cr_to_P = var(pm, nw, :crt, t_idx)
    ci_to_P = var(pm, nw, :cit, t_idx)

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))]
    scale = (tm_scale*pol).*tm_set

    n_phases = length(tm)
    Md = _get_delta_transformation_matrix(n_phases)

    JuMP.@constraint(pm.model, scale.*cr_fr_P .+ cr_to_P .== 0)
    JuMP.@constraint(pm.model, scale.*ci_fr_P .+ ci_to_P .== 0)

    var(pm, nw, :crt_bus)[f_idx] = _merge_bus_flows(pm, Md'*cr_fr_P, f_connections)
    var(pm, nw, :cit_bus)[f_idx] = _merge_bus_flows(pm, Md'*ci_fr_P, f_connections)
    var(pm, nw, :crt_bus)[t_idx] = _merge_bus_flows(pm, [cr_to_P..., -sum(cr_to_P)], t_connections)
    var(pm, nw, :cit_bus)[t_idx] = _merge_bus_flows(pm, [ci_to_P..., -sum(ci_to_P)], t_connections)
end


""
function constraint_mc_transformer_power_rating(pm::AbstractMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    trans = ref(pm, nw, :transformer, id)
    f_bus = trans["f_bus"]
    t_bus = trans["t_bus"]
    f_idx = (id,f_bus,t_bus)
    t_idx = (id,t_bus,f_bus)
    f_conns = trans["f_connections"]
    t_conns = trans["t_connections"]
    config = trans["configuration"]
    sm_ub = trans["sm_ub"]

    constraint_mc_transformer_power_rating(pm, nw, id, f_idx, t_idx, f_bus, t_bus, f_conns, t_conns, config, sm_ub)
end


""
function constraint_mc_transformer_power_rating(pm::AbstractMultiConductorIVRModel, nw::Int, id::Int, f_idx::Tuple, t_idx::Tuple, f_bus::Int, t_bus::Int, f_connections::Vector, t_connections::Vector, config::ConnConfig, sm_ub::Real; report::Bool=true)
    vr_fr = var(pm, nw, :vr, f_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    crt_fr = var(pm, nw, :crt, f_idx)
    cit_fr = var(pm, nw, :cit, f_idx)
    crt_to = var(pm, nw, :crt, t_idx)
    cit_to = var(pm, nw, :cit, t_idx)

    if config==WYE || length(crt_fr)==1
        P_fr = f_connections[1:end-1]
        n_fr = f_connections[end]
        vrt_fr = [vr_fr[p]-vr_fr[n_fr] for p in P_fr]
        vit_fr = [vi_fr[p]-vi_fr[n_fr] for p in P_fr]
    elseif config==DELTA && length(crt_fr)==3
        M = _get_delta_transformation_matrix(3)
        vrt_fr = M*[vr_to[p] for p in f_connections]
        vit_fr = M*[vi_to[p] for p in f_connections]
    else
        error("The configuration $config of dimension $(length(crt)) is not supported.")
    end

    P_to = t_connections[1:end-1]
    n_to = t_connections[end]
    vrt_to = [vr_to[p]-vr_to[n_to] for p in P_to]
    vit_to = [vi_to[p]-vi_to[n_to] for p in P_to]

    idxs = [1:length(vrt_fr)...]
    pt_fr = JuMP.@NLexpression(pm.model, [i in idxs],  vrt_fr[i]*crt_fr[i] + vit_fr[i]*cit_fr[i])
    qt_fr = JuMP.@NLexpression(pm.model, [i in idxs], -vrt_fr[i]*cit_fr[i] + vit_fr[i]*crt_fr[i])
    pt_to = JuMP.@NLexpression(pm.model, [i in idxs],  vrt_to[i]*crt_to[i] + vit_to[i]*cit_to[i])
    qt_to = JuMP.@NLexpression(pm.model, [i in idxs], -vrt_to[i]*cit_to[i] + vit_to[i]*crt_to[i])

    if sm_ub<Inf
        JuMP.@NLconstraint(pm.model, sum(pt_fr[i] for i in idxs)^2 + sum(qt_fr[i] for i in idxs)^2 <= sm_ub^2)
        JuMP.@NLconstraint(pm.model, sum(pt_to[i] for i in idxs)^2 + sum(qt_to[i] for i in idxs)^2 <= sm_ub^2)
    end

    if report
        sol(pm, nw, :transformer, id)[:pt_fr] = pt_fr
        sol(pm, nw, :transformer, id)[:qt_fr] = qt_fr
        sol(pm, nw, :transformer, id)[:pt_to] = pt_to
        sol(pm, nw, :transformer, id)[:qt_to] = qt_to
        sol(pm, nw, :transformer, id)[:smtot_fr] = JuMP.@NLexpression(pm.model, sqrt(sum(pt_fr[i] for i in idxs)^2 + sum(qt_fr[i] for i in idxs)^2))
        sol(pm, nw, :transformer, id)[:smtot_to] = JuMP.@NLexpression(pm.model, sqrt(sum(pt_to[i] for i in idxs)^2 + sum(qt_to[i] for i in idxs)^2))
    end
end


# TRANSFORMER - Constraint - Quadratic

""
function constraint_mc_transformer_power_rating(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, f_idx::Tuple, t_idx::Tuple, f_bus::Int, t_bus::Int, f_connections::Vector, t_connections::Vector, config::ConnConfig, sm_ub::Real; report::Bool=true)
    vr_fr = var(pm, nw, :vr, f_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    crt_fr = var(pm, nw, :crt, f_idx)
    cit_fr = var(pm, nw, :cit, f_idx)
    crt_to = var(pm, nw, :crt, t_idx)
    cit_to = var(pm, nw, :cit, t_idx)

    if config==WYE || length(crt_fr)==1
        P_fr = f_connections[1:end-1]
        n_fr = f_connections[end]
        vrt_fr = [vr_fr[p]-vr_fr[n_fr] for p in P_fr]
        vit_fr = [vi_fr[p]-vi_fr[n_fr] for p in P_fr]
    elseif config==DELTA && length(crt_fr)==3
        M = _get_delta_transformation_matrix(3)
        vrt_fr = M*[vr_to[p] for p in f_connections]
        vit_fr = M*[vi_to[p] for p in f_connections]
    else
        error("The configuration $config of dimension $(length(crt)) is not supported.")
    end

    P_to = t_connections[1:end-1]
    n_to = t_connections[end]
    vrt_to = [vr_to[p]-vr_to[n_to] for p in P_to]
    vit_to = [vi_to[p]-vi_to[n_to] for p in P_to]

    pt_fr = var(pm, nw, :pt, f_idx)
    qt_fr = var(pm, nw, :qt, f_idx)
    pt_to = var(pm, nw, :pt, t_idx)
    qt_to = var(pm, nw, :qt, t_idx)

    JuMP.@constraint(pm.model, pt_fr .==  vrt_fr.*crt_fr .+ vit_fr.*cit_fr)
    JuMP.@constraint(pm.model, qt_fr .== -vrt_fr.*cit_fr .+ vit_fr.*crt_fr)
    JuMP.@constraint(pm.model, pt_to .==  vrt_to.*crt_to .+ vit_to.*cit_to)
    JuMP.@constraint(pm.model, qt_to .== -vrt_to.*cit_to .+ vit_to.*crt_to)

    if sm_ub<Inf
        JuMP.@constraint(pm.model, sum(pt_fr)^2+sum(qt_fr)^2 <= sm_ub^2)
        JuMP.@constraint(pm.model, sum(pt_to)^2+sum(qt_to)^2 <= sm_ub^2)
    end
    
    # if report
    #     sol(pm, nw, :transformer, id)[:sm_fr] = JuMP.@NLexpression(pm.model, sqrt(pt_fr^2+qt_fr^2))
    #     sol(pm, nw, :transformer, id)[:sm_to] = JuMP.@NLexpression(pm.model, sqrt(pt_to^2+qt_to^2))
    # end
end


# BRANCH

# BRANCH - Variables

""
function variable_mc_branch_current(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    variable_mc_branch_current_series_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_branch_current_series_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    var(pm, nw)[:cr_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:ci_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end


"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_current_real(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref(pm, nw, :branch))
    cr = var(pm, nw)[:cr] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_cr_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "cr_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    # if bounded
    #     for (l,i,j) in ref(pm, nw, :arcs_branch)
    #         cmax = _calc_branch_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
    #         for (idx,c) in enumerate(connections[(l,i,j)])
    #             set_upper_bound(cr[(l,i,j)][c],  cmax[idx])
    #             set_lower_bound(cr[(l,i,j)][c], -cmax[idx])
    #         end
    #     end
    # end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch, :cr_fr, :cr_to, ref(pm, nw, :arcs_branch_from), ref(pm, nw, :arcs_branch_to), cr)
end


"variable: `ci[l,i,j] ` for `(l,i,j)` in `arcs`"
function variable_mc_branch_current_imaginary(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref(pm, nw, :branch))
    ci = var(pm, nw)[:ci] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_ci_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "ci_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    # if bounded
    #     for (l,i,j) in ref(pm, nw, :arcs_branch)
    #         cmax = _calc_branch_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
    #         for (idx,c) in enumerate(connections[(l,i,j)])
    #             set_upper_bound(ci[(l,i,j)][c],  cmax[idx])
    #             set_lower_bound(ci[(l,i,j)][c], -cmax[idx])
    #         end
    #     end
    # end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch, :ci_fr, :ci_to, ref(pm, nw, :arcs_branch_from), ref(pm, nw, :arcs_branch_to), ci)
end


"variable: `csr[l]` for `l` in `branch`"
function variable_mc_branch_current_series_real(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref(pm, nw, :branch))
    csr = var(pm, nw)[:csr] = Dict(l => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_csr_$(l)",
            start = comp_start_value(ref(pm, nw, :branch, l), "csr_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    # if bounded
    #     for (l,i,j) in ref(pm, nw, :arcs_branch_from)
    #         cmax = _calc_branch_series_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
    #         for (idx,c) in enumerate(connections[(l,i,j)])
    #             set_upper_bound(csr[l][c],  cmax[idx])
    #             set_lower_bound(csr[l][c], -cmax[idx])
    #         end
    #     end
    # end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch, :csr_fr, ids(pm, nw, :branch), csr)
end


"variable: `csi[l]` for `l` in `branch`"
function variable_mc_branch_current_series_imaginary(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref(pm, nw, :branch))
    csi = var(pm, nw)[:csi] = Dict(l => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_csi_$(l)",
            start = comp_start_value(ref(pm, nw, :branch, l), "csi_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    # if bounded
    #     for (l,i,j) in ref(pm, nw, :arcs_branch_from)
    #         cmax = _calc_branch_series_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
    #         for (idx,c) in enumerate(connections[(l,i,j)])
    #             set_upper_bound(csi[l][c],  cmax[idx])
    #             set_lower_bound(csi[l][c], -cmax[idx])
    #         end
    #     end
    # end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch, :csi_fr, ids(pm, nw, :branch), csi)
end


# BRANCH - Variables - Reduced

""
function variable_mc_branch_current(pm::ReducedMultiConductorIVRModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_branch_current_series_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_branch_current_series_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    var(pm, nw)[:cr] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:ci] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:cr_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:ci_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end


# BRANCH - Constraints

"Defines how current distributes over series and shunt impedances of a pi-model branch"
function constraint_mc_current_from(pm::AbstractMultiConductorIVRModel, nw::Int, f_bus::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real}; report::Bool=true)
    vr_fr = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]

    cr_fr =  var(pm, nw, :cr, f_idx)
    ci_fr =  var(pm, nw, :ci, f_idx)

    csr_fr =  var(pm, nw, :csr, f_idx[1])
    csi_fr =  var(pm, nw, :csi, f_idx[1])

    JuMP.@constraint(pm.model, cr_fr .== csr_fr + g_sh_fr*vr_fr - b_sh_fr*vi_fr)
    JuMP.@constraint(pm.model, ci_fr .== csi_fr + g_sh_fr*vi_fr + b_sh_fr*vr_fr)

    var(pm, nw, :cr_bus)[f_idx] = cr_bus_fr = _merge_bus_flows(pm, cr_fr, f_connections)
    var(pm, nw, :ci_bus)[f_idx] = ci_bus_fr = _merge_bus_flows(pm, ci_fr, f_connections)

    if report
        sol(pm, nw, :branch, f_idx[1])[:cr_fr] = cr_fr
        sol(pm, nw, :branch, f_idx[1])[:ci_fr] = ci_fr
    end
    
end


"Defines how current distributes over series and shunt impedances of a pi-model branch"
function constraint_mc_current_to(pm::AbstractMultiConductorIVRModel, nw::Int, t_bus, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, g_sh_to::Matrix{<:Real}, b_sh_to::Matrix{<:Real}; report::Bool=true)
    vr_to = [var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_to = [var(pm, nw, :vi, t_bus)[c] for c in t_connections]

    cr_to = var(pm, nw, :cr, t_idx)
    ci_to = var(pm, nw, :ci, t_idx)

    csr_to = -var(pm, nw, :csr, f_idx[1])
    csi_to = -var(pm, nw, :csi, f_idx[1])

    JuMP.@constraint(pm.model, cr_to .== csr_to + g_sh_to*vr_to - b_sh_to*vi_to)
    JuMP.@constraint(pm.model, ci_to .== csi_to + g_sh_to*vi_to + b_sh_to*vr_to)

    var(pm, nw, :cr_bus)[t_idx] = cr_bus_to = _merge_bus_flows(pm, cr_to, t_connections)
    var(pm, nw, :ci_bus)[t_idx] = ci_bus_to = _merge_bus_flows(pm, ci_to, t_connections)

    if report
        sol(pm, nw, :branch, f_idx[1])[:cr_to] = cr_to
        sol(pm, nw, :branch, f_idx[1])[:ci_to] = ci_to
        sol(pm, nw, :branch, f_idx[1])[:cr_bus_to] = cr_bus_to
        sol(pm, nw, :branch, f_idx[1])[:ci_bus_to] = ci_bus_to
    end
end


"Defines voltage drop over a branch, linking from and to side complex voltage"
function constraint_mc_bus_voltage_drop(pm::AbstractMultiConductorIVRModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, r::Matrix{<:Real}, x::Matrix{<:Real})
    vr_fr = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]

    vr_to = [var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_to = [var(pm, nw, :vi, t_bus)[c] for c in t_connections]

    csr_fr = var(pm, nw, :csr, f_idx[1])
    csi_fr = var(pm, nw, :csi, f_idx[1])

    JuMP.@constraint(pm.model, vr_to .== vr_fr - r*csr_fr + x*csi_fr)
    JuMP.@constraint(pm.model, vi_to .== vi_fr - r*csi_fr - x*csr_fr)
end


function constraint_mc_branch_current_rating(pm::AbstractMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    branch = ref(pm, nw, :branch, id)
    f_idx = (id,branch["f_bus"],branch["t_bus"])
    t_idx = (id,branch["t_bus"],branch["f_bus"])

    constraint_mc_branch_current_rating(pm, nw, f_idx, t_idx, branch["c_rating_a"])
end


""
function constraint_mc_branch_current_rating(pm::AbstractMultiConductorIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, c_rating::Vector{<:Real}; report::Bool=true)
    cr_fr = var(pm, nw, :cr, f_idx)
    ci_fr = var(pm, nw, :ci, f_idx)
    cr_to = var(pm, nw, :cr, t_idx)
    ci_to = var(pm, nw, :ci, t_idx)

    cnds_finite_rating = [c for (c,r) in enumerate(c_rating) if r<Inf]
    JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_fr[c]^2+ci_fr[c]^2 <= c_rating[c]^2)
    JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_to[c]^2+ci_to[c]^2 <= c_rating[c]^2)
end


# BRANCH - Constraints - Reduced

"Defines how current distributes over series and shunt impedances of a pi-model branch"
function constraint_mc_current_from(pm::ReducedMultiConductorIVRModels, nw::Int, f_bus::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real}; report::Bool=true)
    vr_fr = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]

    csr_fr =  var(pm, nw, :csr, f_idx[1])
    csi_fr =  var(pm, nw, :csi, f_idx[1])

    cr_fr = csr_fr + g_sh_fr*vr_fr - b_sh_fr*vi_fr
    ci_fr = csi_fr + g_sh_fr*vi_fr + b_sh_fr*vr_fr

    var(pm, nw, :cr)[f_idx] = cr_fr
    var(pm, nw, :ci)[f_idx] = ci_fr

    var(pm, nw, :cr_bus)[f_idx] = cr_bus_fr = _merge_bus_flows(pm, cr_fr, f_connections)
    var(pm, nw, :ci_bus)[f_idx] = ci_bus_fr = _merge_bus_flows(pm, ci_fr, f_connections)

    if report
        sol(pm, nw, :branch, f_idx[1])[:cr_fr] = cr_fr
        sol(pm, nw, :branch, f_idx[1])[:ci_fr] = ci_fr
    end
end


"Defines how current distributes over series and shunt impedances of a pi-model branch"
function constraint_mc_current_to(pm::ReducedMultiConductorIVRModels, nw::Int, t_bus, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, g_sh_to::Matrix{<:Real}, b_sh_to::Matrix{<:Real}; report::Bool=true)
    vr_to = [var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_to = [var(pm, nw, :vi, t_bus)[c] for c in t_connections]

    csr_to = -var(pm, nw, :csr, f_idx[1])
    csi_to = -var(pm, nw, :csi, f_idx[1])

    cr_to = csr_to + g_sh_to*vr_to - b_sh_to*vi_to
    ci_to = csi_to + g_sh_to*vi_to + b_sh_to*vr_to

    var(pm, nw, :cr)[t_idx] = cr_to
    var(pm, nw, :ci)[t_idx] = ci_to

    var(pm, nw, :cr_bus)[t_idx] = cr_bus_to = _merge_bus_flows(pm, cr_to, t_connections)
    var(pm, nw, :ci_bus)[t_idx] = ci_bus_to = _merge_bus_flows(pm, ci_to, t_connections)

    if report
        sol(pm, nw, :branch, f_idx[1])[:cr_to] = cr_to
        sol(pm, nw, :branch, f_idx[1])[:ci_to] = ci_to
    end
end


## SHARED 

""
function _merge_bus_flows(pm::AbstractMultiConductorIVRModel, flows::Vector, connections::Vector)::JuMP.Containers.DenseAxisArray
    flows_merged = []
    conns_unique = unique(connections)
    for t in conns_unique
        idxs = findall(connections.==t)
        flows_t = flows[idxs]
        if length(flows_t)==1
            flows_merged_t = flows_t[1]
        elseif any(isa(a, JuMP.NonlinearExpression) for a in flows_t)
            flows_merged_t = JuMP.@NLexpression(pm.model, sum(flows_t[i] for i in 1:length(flows_t)))
        else
            flows_merged_t = sum(flows_t)
        end
        push!(flows_merged, flows_merged_t)
    end
    JuMP.Containers.DenseAxisArray(flows_merged, conns_unique)
end


# STORAGE

# STORAGE - Variable

"variables for modeling storage units, includes grid injection and internal variables"
function variable_mc_storage(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, kwargs...)
    variable_mc_storage_power_real(pm; bounded=bounded, kwargs...)
    variable_mc_storage_power_imaginary(pm; bounded=bounded, kwargs...)
    variable_mc_storage_power_control_imaginary(pm; bounded=false, kwargs...)
    var(pm, nw)[:crs_bus] = Dict{Int,Any}()
    var(pm, nw)[:cis_bus] = Dict{Int,Any}()
    variable_mc_storage_current_real(pm; bounded=bounded, kwargs...)
    variable_mc_storage_current_imaginary(pm; bounded=bounded, kwargs...)
    variable_storage_energy(pm; bounded=bounded, kwargs...)
    variable_storage_charge(pm; bounded=bounded, kwargs...)
    variable_storage_discharge(pm; bounded=bounded, kwargs...)
end


""
function variable_mc_storage_power_real(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i =>_infer_int_dim_unit(strg, false) for (i,strg) in ref(pm, nw, :storage))
    ps = var(pm, nw)[:ps] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_ps_$(i)",
            start = comp_start_value(ref(pm, nw, :storage, i), "ps_start", c, 0.0)
        ) for i in ids(pm, nw, :storage)
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :ps, ids(pm, nw, :storage), ps)
end


""
function variable_mc_storage_power_imaginary(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i =>_infer_int_dim_unit(strg, false) for (i,strg) in ref(pm, nw, :storage))
    qs = var(pm, nw)[:qs] = Dict(i => JuMP.@variable(pm.model,
            [c in int_dim[i]], base_name="$(nw)_qs_$(i)",
            start = comp_start_value(ref(pm, nw, :storage, i), "qs_start", c, 0.0)
        ) for i in ids(pm, nw, :storage)
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :qs, ids(pm, nw, :storage), qs)
end


""
function variable_mc_storage_current_real(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i =>_infer_int_dim_unit(strg, false) for (i,strg) in ref(pm, nw, :storage))
    crs = var(pm, nw)[:crs] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_crs_$(i)",
            start = comp_start_value(ref(pm, nw, :storage, i), "crs_start", c, 0.0)
        ) for i in ids(pm, nw, :storage)
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :crs, ids(pm, nw, :storage), crs)
end


""
function variable_mc_storage_current_imaginary(pm::AbstractMultiConductorIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i =>_infer_int_dim_unit(strg, false) for (i,strg) in ref(pm, nw, :storage))
    cis = var(pm, nw)[:cis] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_cis_$(i)",
            start = comp_start_value(ref(pm, nw, :storage, i), "cis_start", c, 0.0)
        ) for i in ids(pm, nw, :storage)
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :cis, ids(pm, nw, :storage), cis)
end


# STORAGE - Constraints

function constraint_mc_storage_current(pm::AbstractMultiConductorIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)
    strg = ref(pm, nw, :storage, id)

    nphases = _infer_int_dim_unit(strg, false)
    if get(strg, "configuration", WYE) == WYE || nphases==1
        constraint_mc_storage_current_wye(pm, nw, id, strg["connections"]; report=report, bounded=bounded)
    else
        constraint_mc_storage_current_delta(pm, nw, id, strg["connections"]; report=report, bounded=bounded)
    end
end


""
function constraint_mc_storage_current_wye(pm::AbstractMultiConductorIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crs = var(pm, nw, :crs, id)
    cis = var(pm, nw, :cis, id)

    var(pm, nw, :crs_bus)[id] = _merge_bus_flows(pm, [crs..., -sum(crs)], connections)
    var(pm, nw, :cis_bus)[id] = _merge_bus_flows(pm, [cis..., -sum(cis)], connections)
end


""
function constraint_mc_storage_current_delta(pm::AbstractMultiConductorIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crs = var(pm, nw, :crs, id)
    cis = var(pm, nw, :cis, id)
    Md = _get_delta_transformation_matrix(length(connections))
    var(pm, nw, :crs_bus)[id] = _merge_bus_flows(pm, Md'*crs, connections)
    var(pm, nw, :cis_bus)[id] = _merge_bus_flows(pm, Md'*cis, connections)
end


""
function constraint_mc_storage_losses(pm::AbstractMultiConductorIVRModel, i::Int; nw::Int=nw_id_default, kwargs...)
    strg = ref(pm, nw, :storage, i)

    vr  = var(pm, nw,  :vr, strg["storage_bus"])
    vi  = var(pm, nw,  :vi, strg["storage_bus"])
    ps  = var(pm, nw,  :ps, i)
    qs  = var(pm, nw,  :qs, i)
    crs = var(pm, nw,  :crs, i)
    cis = var(pm, nw,  :cis, i)
    sc  = var(pm, nw,  :sc, i)
    sd  = var(pm, nw,  :sd, i)
    qsc = var(pm, nw, :qsc, i)

    p_loss = strg["p_loss"]
    q_loss = strg["q_loss"]
    r = strg["r"]
    x = strg["x"]

    JuMP.@constraint(pm.model,
        sum(ps) + (sd - sc)
        ==
        p_loss + sum(r .* (crs.^2 .+ cis.^2))
    )

    JuMP.@constraint(pm.model,
        sum(qs)
        ==
        qsc + q_loss + sum(x .* (crs.^2 .+ cis.^2))
    )
end


""
function constraint_mc_storage_thermal_limit(pm::AbstractMultiConductorIVRModel, nw::Int, i::Int, connections::Vector{Int}, rating::Vector{<:Real})
    crs = var(pm, nw, :crs, i)
    cis = var(pm, nw, :cis, i)

    
    #TODO is this a current or power bound?
    JuMP.@constraint(pm.model, crs.^2 + cis.^2 .<= rating.^2)
end


""
function constraint_mc_storage_power(pm::AbstractUnbalancedPowerModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)
    strg = ref(pm, nw, :gen, id)
    bus = ref(pm, nw, :bus, strg["storage_bus"])

    if get(strg, "configuration", WYE) == WYE
        constraint_mc_storage_power_wye(pm, nw, id, bus["index"], strg["connections"]; report=report, bounded=bounded)
    else
        constraint_mc_storage_power_delta(pm, nw, id, bus["index"], strg["connections"]; report=report, bounded=bounded)
    end
end


""
function constraint_mc_storage_power_wye(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crs = var(pm, nw, :crs, id)
    cis = var(pm, nw, :cis, id)

    phases = connections[1:end-1]
    n      = connections[end]

    vr_pn = [vr[p]-vr[n] for p in phases]
    vi_pn = [vi[p]-vi[n] for p in phases]

    ps = var(pm, nw, :ps, id)
    qs = var(pm, nw, :qs, id)

    JuMP.@constraint(pm.model, ps .==  vr_pn.*crs .+ vi_pn.*cis)
    JuMP.@constraint(pm.model, qs .== -vr_pn.*cis .+ vi_pn.*crs)
end


""
function constraint_mc_storage_power_delta(pm::AbstractQuadraticMultiConductorIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crs = var(pm, nw, :crs, id)
    cis = var(pm, nw, :cis, id)

    nph = length(pmin)

    prev = Dict(c=>connections[(idx+nph-2)%nph+1] for (idx,c) in enumerate(connections))
    next = Dict(c=>connections[idx%nph+1] for (idx,c) in enumerate(connections))

    vrs = [vr[c]-vr[next[c]] for p in connections]
    vis = [vi[c]-vi[next[c]] for p in connections]

    ps = var(pm, nw, :ps, id)
    qs = var(pm, nw, :qs, id)

    JuMP.@constraint(pm.model, ps .==  vrs.*crs .+ vis.*cis)
    JuMP.@constraint(pm.model, qs .== -vrs.*cis .+ vis.*crs)
end

