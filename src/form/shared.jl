import LinearAlgebra: diag


"`vm[i] == vmref`"
function constraint_mc_voltage_magnitude_only(pm::AbstractUnbalancedWModels, nw::Int, i::Int, vm_ref::Vector{<:Real})
    w = [var(pm, nw, :w, i)[t] for t in ref(pm, nw, :bus, i)["terminals"]]
    JuMP.@constraint(pm.model, w .== vm_ref.^2)
end


""
function constraint_mc_switch_state_closed(pm::AbstractUnbalancedWModels, nw::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int})
    w_fr = var(pm, nw, :w, f_bus)
    w_to = var(pm, nw, :w, t_bus)

    for (idx, (fc, tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, w_fr[fc] == w_to[tc])
    end
end


""
function constraint_mc_switch_state_on_off(pm::AbstractUnbalancedWModels, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int}; relax::Bool=false)
    w_fr = var(pm, nw, :w, f_bus)
    w_to = var(pm, nw, :w, t_bus)

    z = var(pm, nw, :switch_state, i)

    for (fc, tc) in zip(f_connections, t_connections)
        if relax
            M = 1e20
            JuMP.@constraint(pm.model, w_fr[fc] - w_to[tc] <=  M * (1-z))
            JuMP.@constraint(pm.model, w_fr[fc] - w_to[tc] >= -M * (1-z))
        else
            JuMP.@constraint(pm.model, z => {w_fr[fc] == w_to[tc]})
        end
    end
end


""
function constraint_mc_power_balance_slack(pm::AbstractUnbalancedWModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{<:Tuple{Tuple{Int,Int,Int},Vector{Union{String,Int}}}}, bus_arcs_sw::Vector{<:Tuple{Tuple{Int,Int,Int},Vector{Union{String,Int}}}}, bus_arcs_trans::Vector{<:Tuple{Tuple{Int,Int,Int},Vector{Union{String,Int}}}}, bus_gens::Vector{<:Tuple{Int,Vector{Union{String,Int}}}}, bus_storage::Vector{<:Tuple{Int,Vector{Union{String,Int}}}}, bus_loads::Vector{<:Tuple{Int,Vector{Union{String,Int}}}}, bus_shunts::Vector{<:Tuple{Int,Vector{Union{String,Int}}}})
    w    = var(pm, nw, :w, i)
    p    = get(var(pm, nw),    :p, Dict()); _check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(var(pm, nw),    :q, Dict()); _check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, nw),   :pg, Dict()); _check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, nw),   :qg, Dict()); _check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, nw),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, nw),   :qs, Dict()); _check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, nw),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, nw),  :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(var(pm, nw),   :pt, Dict()); _check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(var(pm, nw),   :qt, Dict()); _check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    p_slack = var(pm, nw, :p_slack, i)
    q_slack = var(pm, nw, :q_slack, i)

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx, t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
              sum(p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(pt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
            ==
            sum(pg[g][t] for (g, conns) in bus_gens if t in conns)
            - sum(ps[s][t] for (s, conns) in bus_storage if t in conns)
            - sum(ref(pm, nw, :load, l, "pd")[findfirst(isequal(t), conns)] for (l, conns) in bus_loads if t in conns)
            - sum(w[t] * diag(Gt')[idx] for (sh, conns) in bus_shunts if t in conns)
            + p_slack[t]
        )
        push!(cstr_p, cp)
        cq = JuMP.@constraint(pm.model,
              sum(q[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(qsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(qt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
            ==
            sum(qg[g][t] for (g, conns) in bus_gens if t in conns)
            - sum(qs[s][t] for (s, conns) in bus_storage if t in conns)
            - sum(ref(pm, nw, :load, l, "qd")[findfirst(isequal(t), conns)] for (l, conns) in bus_loads if t in conns)
            - sum(-w[t] * diag(Bt')[idx] for (sh, conns) in bus_shunts if t in conns)
            + q_slack[t]
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


"do nothing, no way to represent this in these variables"
function constraint_mc_theta_ref(pm::AbstractUnbalancedWModels, n::Int, d::Int, va_ref)
end


"Creates phase angle constraints at reference buses"
function constraint_mc_theta_ref(pm::AbstractUnbalancedPolarModels, nw::Int, i::Int, va_ref::Vector{<:Real})
    terminals = ref(pm, nw, :bus, i)["terminals"]

    va = [var(pm, nw, :va, i)[t] for t in terminals]

    JuMP.@constraint(pm.model, va .== va_ref)
end


"""
For a variable tap transformer, fix the tap variables which are fixed. For
example, an OLTC where the third phase is fixed, will have tap variables for
all phases, but the third tap variable should be fixed.
"""
function constraint_mc_oltc_tap_fix(pm::AbstractUnbalancedPowerModel, i::Int, fixed::Vector, tm::Vector; nw=nw_id_default)
    for (c,fixed) in enumerate(fixed)
        if fixed
            JuMP.@constraint(pm.model, var(pm, nw, c, :tap)[i]==tm[c])
        end
    end
end


"KCL for load shed problem with transformers (AbstractWForms)"
function constraint_mc_power_balance_shed(pm::AbstractUnbalancedWModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    w        = var(pm, nw, :w, i)
    p        = get(var(pm, nw),    :p, Dict()); _check_var_keys(p, bus_arcs, "active power", "branch")
    q        = get(var(pm, nw),    :q, Dict()); _check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg       = get(var(pm, nw),   :pg, Dict()); _check_var_keys(pg, bus_gens, "active power", "generator")
    qg       = get(var(pm, nw),   :qg, Dict()); _check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps       = get(var(pm, nw),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    qs       = get(var(pm, nw),   :qs, Dict()); _check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw      = get(var(pm, nw),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw      = get(var(pm, nw),  :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt       = get(var(pm, nw),   :pt, Dict()); _check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt       = get(var(pm, nw),   :qt, Dict()); _check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    z_demand = var(pm, nw, :z_demand)
    z_shunt  = var(pm, nw, :z_shunt)

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx, t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
              sum(p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(pt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
            ==
            sum(pg[g][t] for (g, conns) in bus_gens if t in conns)
            - sum(ps[s][t] for (s, conns) in bus_storage if t in conns)
            - sum(ref(pm, nw, :load, l, "pd")[findfirst(isequal(t), conns)] * z_demand[l] for (l, conns) in bus_loads if t in conns)
            - sum(z_shunt[sh] *(w[t] * diag(Gt')[idx]) for (sh, conns) in bus_shunts if t in conns)
        )
        push!(cstr_p, cp)
        cq = JuMP.@constraint(pm.model,
              sum(q[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(qsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(qt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
            ==
            sum(qg[g][t] for (g, conns) in bus_gens if t in conns)
            - sum(qs[s][t] for (s, conns) in bus_storage if t in conns)
            - sum(ref(pm, nw, :load, l, "qd")[findfirst(isequal(t), conns)]*z_demand[l] for (l, conns) in bus_loads if t in conns)
            - sum(z_shunt[sh] * (-w[t] * diag(Bt')[idx]) for (sh, conns) in bus_shunts if t in conns)
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


""
function constraint_mc_power_balance(pm::AbstractUnbalancedWModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    Wr = var(pm, nw, :Wr, i)
    Wi = var(pm, nw, :Wi, i)
    P = get(var(pm, nw), :P, Dict()); _check_var_keys(P, bus_arcs, "active power", "branch")
    Q = get(var(pm, nw), :Q, Dict()); _check_var_keys(Q, bus_arcs, "reactive power", "branch")
    Psw  = get(var(pm, nw),  :Psw, Dict()); _check_var_keys(Psw, bus_arcs_sw, "active power", "switch")
    Qsw  = get(var(pm, nw),  :Qsw, Dict()); _check_var_keys(Qsw, bus_arcs_sw, "reactive power", "switch")
    Pt   = get(var(pm, nw),   :Pt, Dict()); _check_var_keys(Pt, bus_arcs_trans, "active power", "transformer")
    Qt   = get(var(pm, nw),   :Qt, Dict()); _check_var_keys(Qt, bus_arcs_trans, "reactive power", "transformer")

    pd = get(var(pm, nw), :pd_bus, Dict()); _check_var_keys(pd, bus_loads, "active power", "load")
    qd = get(var(pm, nw), :qd_bus, Dict()); _check_var_keys(qd, bus_loads, "reactive power", "load")
    pg = get(var(pm, nw), :pg_bus, Dict()); _check_var_keys(pg, bus_gens, "active power", "generator")
    qg = get(var(pm, nw), :qg_bus, Dict()); _check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, nw),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, nw),   :qs, Dict()); _check_var_keys(qs, bus_storage, "reactive power", "storage")

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
            sum(diag(P[a])[findfirst(isequal(t), conns)] for (a, conns) in bus_arcs if t in conns)
            + sum(diag(Psw[a_sw])[findfirst(isequal(t), conns)] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(diag(Pt[a_trans])[findfirst(isequal(t), conns)] for (a_trans, conns) in bus_arcs_trans if t in conns)
            ==
            sum(pg[g][t] for (g, conns) in bus_gens if t in conns)
            - sum(ps[s][t] for (s, conns) in bus_storage if t in conns)
            - sum(pd[d][t] for (d, conns) in bus_loads if t in conns)
            - diag(Wr*Gt'+Wi*Bt')[idx]
        )
        push!(cstr_p, cp)

        cq = JuMP.@constraint(pm.model,
            sum(diag(Q[a])[findfirst(isequal(t), conns)] for (a, conns) in bus_arcs if t in conns)
            + sum(diag(Qsw[a_sw])[findfirst(isequal(t), conns)] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(diag(Qt[a_trans])[findfirst(isequal(t), conns)] for (a_trans, conns) in bus_arcs_trans if t in conns)
            ==
            sum(qg[g][t] for (g, conns) in bus_gens if t in conns)
            - sum(qs[s][t] for (s, conns) in bus_storage if t in conns)
            - sum(qd[d][t] for (d, conns) in bus_loads if t in conns)
            - diag(-Wr*Bt'+Wi*Gt')[idx]
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


"Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)"
function constraint_mc_ohms_yt_from(pm::AbstractUnbalancedWModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    w_fr = var(pm, n, :w, f_bus)
    wr   = var(pm, n, :wr, (f_bus, t_bus))
    wi   = var(pm, n, :wi, (f_bus, t_bus))

    JuMP.@constraint(pm.model, p_fr ==  (g+g_fr)/tm^2*w_fr + (-g*tr+b*ti)/tm^2*wr + (-b*tr-g*ti)/tm^2*wi )
    JuMP.@constraint(pm.model, q_fr == -(b+b_fr)/tm^2*w_fr - (-b*tr-g*ti)/tm^2*wr + (-g*tr+b*ti)/tm^2*wi )
end


"Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)"
function constraint_mc_ohms_yt_to(pm::AbstractUnbalancedWModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    q_to = var(pm, n, :q, t_idx)
    p_to = var(pm, n, :p, t_idx)
    w_to = var(pm, n, :w, t_bus)
    wr   = var(pm, n, :wr, (f_bus, t_bus))
    wi   = var(pm, n, :wi, (f_bus, t_bus))

    JuMP.@constraint(pm.model, p_to ==  (g+g_to)*w_to + (-g*tr-b*ti)/tm^2*wr + (-b*tr+g*ti)/tm^2*-wi )
    JuMP.@constraint(pm.model, q_to == -(b+b_to)*w_to - (-b*tr+g*ti)/tm^2*wr + (-g*tr-b*ti)/tm^2*-wi )
end


"on/off bus voltage constraint for relaxed forms"
function constraint_mc_bus_voltage_on_off(pm::AbstractUnbalancedWModels, n::Int; kwargs...)
    for (i, bus) in ref(pm, n, :bus)
        constraint_mc_bus_voltage_magnitude_sqr_on_off(pm, i, nw=n)
    end
end


""
function constraint_mc_voltage_angle_difference(pm::AbstractUnbalancedPolarModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})
    i, f_bus, t_bus = f_idx

    va_fr = [var(pm, nw, :va, f_bus)[fc] for fc in f_connections]
    va_to = [var(pm, nw, :va, t_bus)[tc] for tc in t_connections]

    JuMP.@constraint(pm.model, va_fr .- va_to .<= angmax)
    JuMP.@constraint(pm.model, va_fr .- va_to .>= angmin)
end


""
function constraint_mc_storage_on_off(pm::AbstractUnbalancedPowerModel, nw::Int, i::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}, charge_ub, discharge_ub)
    z_storage =var(pm, nw, :z_storage, i)
    ps = [var(pm, nw, :ps, i)[c] for c in connections]
    qs = [var(pm, nw, :qs, i)[c] for c in connections]

    JuMP.@constraint(pm.model, ps .<= z_storage.*pmax)
    JuMP.@constraint(pm.model, ps .>= z_storage.*pmin)

    JuMP.@constraint(pm.model, qs .<= z_storage.*qmax)
    JuMP.@constraint(pm.model, qs .>= z_storage.*qmin)
end


""
function constraint_mc_generator_power_wye(pm::AbstractUnbalancedPowerModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    var(pm, nw, :pg_bus)[id] = var(pm, nw, :pg, id)
    var(pm, nw, :qg_bus)[id] = var(pm, nw, :qg, id)

    if report
        sol(pm, nw, :gen, id)[:pg_bus] = var(pm, nw, :pg_bus, id)
        sol(pm, nw, :gen, id)[:qg_bus] = var(pm, nw, :qg_bus, id)
    end
end


"do nothing by default but some formulations require this"
function variable_mc_storage_current(pm::AbstractUnbalancedWConvexModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i, strg) in ref(pm, nw, :storage))
    ccms = var(pm, nw)[:ccms] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_ccms_$(i)",
            start = comp_start_value(ref(pm, nw, :storage, i), "ccms_start", c, 0.0)
        ) for i in ids(pm, nw, :storage)
    )

    if bounded
        bus = ref(pm, nw, :bus)
        for (i, storage) in ref(pm, nw, :storage)
            for (idx, c) in enumerate(connections[i])
                ub = Inf
                if haskey(storage, "thermal_rating")
                    sb = bus[storage["storage_bus"]]
                    ub = (storage["thermal_rating"][idx]/sb["vmin"][findfirst(isequal(c), bus[storage["storage_bus"]]["terminals"])])^2
                end

                set_lower_bound(ccms[i][c], 0.0)
                if !isinf(ub)
                    set_upper_bound(ccms[i][c], ub)
                end
            end
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :storage, :ccms, ids(pm, nw, :storage), ccms)
end


""
function constraint_mc_storage_losses(pm::AbstractUnbalancedWConvexModels, n::Int, i::Int, bus::Int, connections::Vector{Int}, r::Real, x::Real, p_loss::Real, q_loss::Real)
    w = var(pm, n, :w, bus)
    ccms = var(pm, n, :ccms, i)
    ps = var(pm, n, :ps, i)
    qs = var(pm, n, :qs, i)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    qsc = var(pm, n, :qsc, i)

    for c in connections
        JuMP.@constraint(pm.model, ps[c]^2 + qs[c]^2 <= w[c]*ccms[c])
    end

    JuMP.@constraint(pm.model,
        sum(ps[c] for c in connections) + (sd - sc)
        ==
        p_loss + r * sum(ccms[c] for c in connections)
    )

    JuMP.@constraint(pm.model,
        sum(qs[c] for c in connections)
        ==
        qsc + q_loss + x * sum(ccms[c] for c in connections)
    )
end


@doc raw"""
    constraint_mc_ampacity_from(pm::AbstractUnbalancedWModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches from-side

math```
p_{fr}^2 + q_{fr}^2 \leq w_{fr} i_{max}^2
```
"""
function constraint_mc_ampacity_from(pm::AbstractUnbalancedWModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_fr = [var(pm, nw, :p, f_idx)[c] for c in f_connections]
    q_fr = [var(pm, nw, :q, f_idx)[c] for c in f_connections]
    w_fr = [var(pm, nw, :w, f_idx[2])[c] for c in f_connections]

    con(pm, nw, :mu_cm_branch)[f_idx] = mu_cm_fr = [JuMP.@constraint(pm.model, p_fr[idx]^2 + q_fr[idx]^2 .<= w_fr[idx] * c_rating[idx]^2) for idx in findall(c_rating .< Inf)]

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, f_idx[1])[:mu_cm_fr] = mu_cm_fr
    end

    nothing
end


@doc raw"""
    constraint_mc_ampacity_to(pm::AbstractUnbalancedWModels, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches to-side

math```
p_{to}^2 + q_{to}^2 \leq w_{to} i_{max}^2
```
"""
function constraint_mc_ampacity_to(pm::AbstractUnbalancedWModels, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_to = [var(pm, nw, :p, t_idx)[c] for c in t_connections]
    q_to = [var(pm, nw, :q, t_idx)[c] for c in t_connections]
    w_to = [var(pm, nw, :w, t_idx[2])[c] for c in t_connections]

    con(pm, nw, :mu_cm_branch)[t_idx] = mu_cm_to = [JuMP.@constraint(pm.model, p_to[idx]^2 + q_to[idx]^2 .<= w_to[idx] * c_rating[idx]^2) for idx in findall(c_rating .< Inf)]

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, t_idx[1])[:mu_cm_to] = mu_cm_to
    end

    nothing
end


@doc raw"""
    constraint_mc_switch_ampacity(pm::AbstractUnbalancedWModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on switches from-side

math```
p_{fr}^2 + q_{fr}^2 \leq w_{fr} i_{max}^2
```
"""
function constraint_mc_switch_ampacity(pm::AbstractUnbalancedWModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    psw_fr = [var(pm, nw, :psw, f_idx)[c] for c in f_connections]
    qsw_fr = [var(pm, nw, :qsw, f_idx)[c] for c in f_connections]
    w_fr = [var(pm, nw, :w, f_idx[2])[c] for c in f_connections]

    con(pm, nw, :mu_cm_switch)[f_idx] = mu_cm_fr = [JuMP.@constraint(pm.model, psw_fr[idx]^2 + qsw_fr[idx]^2 .<= w_fr[idx] * c_rating[idx]^2) for idx in findall(c_rating .< Inf)]

    if _IM.report_duals(pm)
        sol(pm, nw, :switch, f_idx[1])[:mu_cm_fr] = mu_cm_fr
    end

    nothing
end
