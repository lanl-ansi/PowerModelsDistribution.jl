import LinearAlgebra: diag


"`vm[i] == vmref`"
function constraint_mc_voltage_magnitude_only(pm::_PM.AbstractWModels, nw::Int, i::Int, vm_ref::Vector{<:Real})
    w = [var(pm, nw, :w, i)[t] for t in ref(pm, nw, :bus, i)["terminals"]]
    JuMP.@constraint(pm.model, w .== vm_ref.^2)
end


""
function constraint_mc_switch_state_closed(pm::_PM.AbstractWModels, nw::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int})
    w_fr = var(pm, nw, :w, f_bus)
    w_to = var(pm, nw, :w, t_bus)

    for (idx, (fc, tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, w_fr[fc] == w_to[tc])
    end
end


""
function constraint_mc_switch_state_on_off(pm::_PM.AbstractWModels, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int}; relax::Bool=false)
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
function constraint_mc_power_balance_slack(pm::_PM.AbstractWModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{<:Tuple{Tuple{Int,Int,Int},Vector{Union{String,Int}}}}, bus_arcs_sw::Vector{<:Tuple{Tuple{Int,Int,Int},Vector{Union{String,Int}}}}, bus_arcs_trans::Vector{<:Tuple{Tuple{Int,Int,Int},Vector{Union{String,Int}}}}, bus_gens::Vector{<:Tuple{Int,Vector{Union{String,Int}}}}, bus_storage::Vector{<:Tuple{Int,Vector{Union{String,Int}}}}, bus_loads::Vector{<:Tuple{Int,Vector{Union{String,Int}}}}, bus_shunts::Vector{<:Tuple{Int,Vector{Union{String,Int}}}})
    w    = var(pm, nw, :w, i)
    p    = get(var(pm, nw),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(var(pm, nw),    :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, nw),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, nw),   :qg, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, nw),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, nw),   :qs, Dict()); _PM._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, nw),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, nw),  :qsw, Dict()); _PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(var(pm, nw),   :pt, Dict()); _PM._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(var(pm, nw),   :qt, Dict()); _PM._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
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
function constraint_mc_theta_ref(pm::_PM.AbstractWModels, n::Int, d::Int, va_ref)
end


"Creates phase angle constraints at reference buses"
function constraint_mc_theta_ref(pm::_PM.AbstractPolarModels, nw::Int, i::Int, va_ref::Vector{<:Real})
    terminals = ref(pm, nw, :bus, i)["terminals"]

    va = [var(pm, nw, :va, i)[t] for t in terminals]

    JuMP.@constraint(pm.model, va .== va_ref)
end


"""
For a variable tap transformer, fix the tap variables which are fixed. For
example, an OLTC where the third phase is fixed, will have tap variables for
all phases, but the third tap variable should be fixed.
"""
function constraint_mc_oltc_tap_fix(pm::_PM.AbstractPowerModel, i::Int, fixed::Vector, tm::Vector; nw=pm.cnw)
    for (c,fixed) in enumerate(fixed)
        if fixed
            JuMP.@constraint(pm.model, var(pm, nw, c, :tap)[i]==tm[c])
        end
    end
end


"KCL for load shed problem with transformers (AbstractWForms)"
function constraint_mc_power_balance_shed(pm::_PM.AbstractWModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    w        = var(pm, nw, :w, i)
    p        = get(var(pm, nw),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q        = get(var(pm, nw),    :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg       = get(var(pm, nw),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg       = get(var(pm, nw),   :qg, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps       = get(var(pm, nw),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    qs       = get(var(pm, nw),   :qs, Dict()); _PM._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw      = get(var(pm, nw),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw      = get(var(pm, nw),  :qsw, Dict()); _PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt       = get(var(pm, nw),   :pt, Dict()); _PM._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt       = get(var(pm, nw),   :qt, Dict()); _PM._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
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
function constraint_mc_power_balance(pm::_PM.AbstractWModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    Wr = var(pm, nw, :Wr, i)
    Wi = var(pm, nw, :Wi, i)
    P = get(var(pm, nw), :P, Dict()); _PM._check_var_keys(P, bus_arcs, "active power", "branch")
    Q = get(var(pm, nw), :Q, Dict()); _PM._check_var_keys(Q, bus_arcs, "reactive power", "branch")
    Psw  = get(var(pm, nw),  :Psw, Dict()); _PM._check_var_keys(Psw, bus_arcs_sw, "active power", "switch")
    Qsw  = get(var(pm, nw),  :Qsw, Dict()); _PM._check_var_keys(Qsw, bus_arcs_sw, "reactive power", "switch")
    Pt   = get(var(pm, nw),   :Pt, Dict()); _PM._check_var_keys(Pt, bus_arcs_trans, "active power", "transformer")
    Qt   = get(var(pm, nw),   :Qt, Dict()); _PM._check_var_keys(Qt, bus_arcs_trans, "reactive power", "transformer")

    pd = get(var(pm, nw), :pd_bus, Dict()); _PM._check_var_keys(pd, bus_loads, "active power", "load")
    qd = get(var(pm, nw), :qd_bus, Dict()); _PM._check_var_keys(qd, bus_loads, "reactive power", "load")
    pg = get(var(pm, nw), :pg_bus, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg = get(var(pm, nw), :qg_bus, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, nw),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, nw),   :qs, Dict()); _PM._check_var_keys(qs, bus_storage, "reactive power", "storage")

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


"delegate back to PowerModels"
function constraint_mc_ohms_yt_from(pm::_PM.AbstractWModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    _PM.constraint_ohms_yt_from(pm, n, c, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


"delegate back to PowerModels"
function constraint_mc_ohms_yt_to(pm::_PM.AbstractWModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    _PM.constraint_ohms_yt_to(pm, n, c, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end


"on/off bus voltage constraint for relaxed forms"
function constraint_mc_bus_voltage_on_off(pm::_PM.AbstractWModels, n::Int; kwargs...)
    for (i, bus) in ref(pm, n, :bus)
        constraint_mc_bus_voltage_magnitude_sqr_on_off(pm, i, nw=n)
    end
end


""
function constraint_mc_voltage_angle_difference(pm::_PM.AbstractPolarModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})
    i, f_bus, t_bus = f_idx

    va_fr = [var(pm, nw, :va, f_bus)[fc] for fc in f_connections]
    va_to = [var(pm, nw, :va, t_bus)[tc] for tc in t_connections]

    JuMP.@constraint(pm.model, va_fr .- va_to .<= angmax)
    JuMP.@constraint(pm.model, va_fr .- va_to .>= angmin)
end


""
function constraint_mc_voltage_angle_difference(pm::_PM.AbstractWModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})
    i, f_bus, t_bus = f_idx

    w_fr = var(pm, nw, :w, f_bus)
    w_to = var(pm, nw, :w, t_bus)
    wr   = [var(pm, nw, :wr)[(f_bus, t_bus, fc, tc)] for (fc,tc) in zip(f_connections,t_connections)]
    wi   = [var(pm, nw, :wi)[(f_bus, t_bus, fc, tc)] for (fc,tc) in zip(f_connections,t_connections)]

    JuMP.@constraint(pm.model, wi .<= tan.(angmax).*wr)
    JuMP.@constraint(pm.model, wi .>= tan.(angmin).*wr)

    for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))
        _PM.cut_complex_product_and_angle_difference(pm.model, w_fr[fc], w_to[tc], wr[idx], wi[idx], angmin[idx], angmax[idx])
    end
end


""
function constraint_mc_storage_on_off(pm::_PM.AbstractPowerModel, nw::Int, i::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}, charge_ub, discharge_ub)
    z_storage =var(pm, nw, :z_storage, i)
    ps = [var(pm, nw, :ps, i)[c] for c in connections]
    qs = [var(pm, nw, :qs, i)[c] for c in connections]

    JuMP.@constraint(pm.model, ps .<= z_storage.*pmax)
    JuMP.@constraint(pm.model, ps .>= z_storage.*pmin)

    JuMP.@constraint(pm.model, qs .<= z_storage.*qmax)
    JuMP.@constraint(pm.model, qs .>= z_storage.*qmin)
end


""
function constraint_mc_generator_power_wye(pm::_PM.AbstractPowerModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    var(pm, nw, :pg_bus)[id] = var(pm, nw, :pg, id)
    var(pm, nw, :qg_bus)[id] = var(pm, nw, :qg, id)

    if report
        sol(pm, nw, :gen, id)[:pg_bus] = var(pm, nw, :pg_bus, id)
        sol(pm, nw, :gen, id)[:qg_bus] = var(pm, nw, :qg_bus, id)
    end
end


"do nothing by default but some formulations require this"
function variable_mc_storage_current(pm::_PM.AbstractWConvexModels; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
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

    report && _IM.sol_component_value(pm, nw, :storage, :ccms, ids(pm, nw, :storage), ccms)
end
