""
function variable_mc_bus_voltage(pm::AbstractUnbalancedDCPModel; nw=nw_id_default, kwargs...)
    variable_mc_bus_voltage_angle(pm; nw=nw, kwargs...)
end


######## AbstractDCPForm Models (has va but assumes vm is 1.0) ########

"nothing to do, these models do not have complex voltage constraints"
function constraint_mc_model_voltage(pm::AbstractUnbalancedDCPModel, n::Int, c::Int)
end


""
function variable_mc_bus_voltage_on_off(pm::AbstractUnbalancedDCPModel; kwargs...)
    variable_mc_bus_voltage_angle(pm; kwargs...)
end


### DC Power Flow Approximation ###

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] == -b*(t[f_bus] - t[t_bus])
```
"""
function constraint_mc_ohms_yt_from(pm::AbstractUnbalancedDCPModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})
    p_fr  = var(pm, nw,  :p, f_idx)
    va_fr = var(pm, nw, :va, f_bus)
    va_to = var(pm, nw, :va, t_bus)

    for (idx, (fc,tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, p_fr[fc] == -sum(B[idx,jdx]*(va_fr[fc] - va_to[td]) for (jdx, (fd,td)) in enumerate(zip(f_connections, t_connections))))
    end
end


"power balance constraint with line shunts and transformers for load shed problem, DCP formulation"
function constraint_mc_power_balance_shed(pm::AbstractUnbalancedDCPModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    p    = get(var(pm, nw), :p,      Dict()); _check_var_keys(p, bus_arcs, "active power", "branch")
    pg   = get(var(pm, nw), :pg_bus, Dict()); _check_var_keys(pg, bus_gens, "active power", "generator")
    ps   = get(var(pm, nw), :ps,     Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    psw  = get(var(pm, nw), :psw,    Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    pt   = get(var(pm, nw), :pt,     Dict()); _check_var_keys(pt, bus_arcs_trans, "active power", "transformer")

    z_demand = var(pm, nw, :z_demand)
    z_gen = haskey(var(pm, nw), :z_gen) ? var(pm, nw, :z_gen) : Dict(i => 1.0 for i in ids(pm, nw, :gen))
    z_storage = haskey(var(pm, nw), :z_storage) ? var(pm, nw, :z_storage) : Dict(i => 1.0 for i in ids(pm, nw, :storage))
    z_shunt  = haskey(var(pm, nw), :z_shunt) ? var(pm, nw, :z_shunt) : Dict(i => 1.0 for i in ids(pm, nw, :shunt))

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
              sum(p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(pt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
            + sum(ref(pm, nw, :load, d, "pd")[findfirst(isequal(t), conns)]*z_demand[d] for (d, conns) in bus_loads if t in conns)
            - sum(pg[g][t]*z_gen[g] for (g, conns) in bus_gens if t in conns)
            - sum(ps[s][t]*z_storage[s] for (s, conns) in bus_storage if t in conns)
            + sum(LinearAlgebra.diag(Gt)[idx]*z_shunt[sh] for (sh, conns) in bus_shunts if t in conns)
            == 0
        )
        push!(cstr_p, cp)
    end

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = []

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
    end
end


"power balance constraint with line shunts and transformers for load shed problem, DCP formulation"
function constraint_mc_power_balance_shed_simple(pm::AbstractUnbalancedDCPModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    p    = get(var(pm, nw), :p,      Dict()); _check_var_keys(p, bus_arcs, "active power", "branch")
    pg   = get(var(pm, nw), :pg_bus, Dict()); _check_var_keys(pg, bus_gens, "active power", "generator")
    ps   = get(var(pm, nw), :ps,     Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    psw  = get(var(pm, nw), :psw,    Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    pt   = get(var(pm, nw), :pt,     Dict()); _check_var_keys(pt, bus_arcs_trans, "active power", "transformer")

    z_demand = var(pm, nw, :z_demand)

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
              sum(p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(pt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
            + sum(ref(pm, nw, :load, d, "pd")[findfirst(isequal(t), conns)]*z_demand[d] for (d, conns) in bus_loads if t in conns)
            - sum(pg[g][t]*z_gen[g] for (g, conns) in bus_gens if t in conns)
            - sum(ps[s][t]*z_storage[s] for (s, conns) in bus_storage if t in conns)
            + sum(LinearAlgebra.diag(Gt)[idx]*z_shunt[sh] for (sh, conns) in bus_shunts if t in conns)
            == 0
        )
        push!(cstr_p, cp)
    end

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = []

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
    end
end



"on/off bus voltage constraint for DCP formulation, nothing to do"
function constraint_mc_bus_voltage_on_off(pm::AbstractUnbalancedDCPModel; nw::Int=nw_id_default, kwargs...)
end


""
function variable_mc_branch_power_real(pm::AbstractUnbalancedAPLossLessModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)
    p = Dict((l,i,j) => JuMP.@variable(pm.model,
        [c in connections[(l,i,j)]], base_name="$(nw)_($l,$i,$j)_p",
        start = comp_start_value(ref(pm, nw, :branch, l), "p_start", c, 0.0)
    ) for (l,i,j) in ref(pm, nw, :arcs_branch_from))

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_branch_from)
            smax = _calc_branch_power_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
            for (k,t) in enumerate(connections[(l,i,j)])
                if !ismissing(smax)
                    set_upper_bound.(p[(l,i,j)][t],  smax[k])
                    set_lower_bound.(p[(l,i,j)][t], -smax[k])
                end
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
    p_expr = Dict{Any,Any}( ((l,i,j), p[(l,i,j)]) for (l,i,j) in ref(pm, nw, :arcs_branch_from) )
    p_expr = merge(p_expr, Dict( ((l,j,i), -1.0*p[(l,i,j)]) for (l,i,j) in ref(pm, nw, :arcs_branch_from)))
    var(pm, nw)[:p] = p_expr
end


""
function constraint_mc_switch_state_closed(pm::AbstractUnbalancedDCPModel, nw::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int})
    va_fr = var(pm, nw, :va, f_bus)
    va_to = var(pm, nw, :va, t_bus)

    for (idx,(fc,tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, va_fr[fc] == va_to[fc])
    end
end


""
function constraint_mc_storage_losses(pm::AbstractUnbalancedAPLossLessModels, n::Int, i, bus, r, x, p_loss, q_loss; conductors=[1])
    ps = var(pm, n, :ps, i)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)

    JuMP.@constraint(pm.model, sum(ps[c] for c in conductors) + (sd - sc) == p_loss)
end
