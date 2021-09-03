""
function variable_mc_bus_voltage(pm::AbstractUnbalancedACRModel; nw::Int=nw_id_default, bounded::Bool=true, kwargs...)
    variable_mc_bus_voltage_real(pm; nw=nw, bounded=bounded, kwargs...)
    variable_mc_bus_voltage_imaginary(pm; nw=nw, bounded=bounded, kwargs...)

    # local infeasbility issues without proper initialization;
    # convergence issues start when the equivalent angles of the starting point
    # are further away than 90 degrees from the solution (as given by ACP)
    # this is the default behaviour of _PM, initialize all phases as (1,0)
    # the magnitude seems to have little effect on the convergence (>0.05)
    # updating the starting point to a balanced phasor does the job
    for id in ids(pm, nw, :bus)
        busref = ref(pm, nw, :bus, id)
        terminals = busref["terminals"]
        grounded = busref["grounded"]

        ncnd = length(terminals)

        if haskey(busref, "vr_start") && haskey(busref, "vm_start")
            vr = busref["vr_start"]
            vi = busref["vi_start"]
        else
            vm = haskey(busref, "vm_start") ? busref["vm_start"] : fill(0.0, ncnd)
            vm[.!grounded] .= 1.0

            # TODO how to do this more generally?
            nph = 3
            va = haskey(busref, "va_start") ? busref["va_start"] : [c <= nph ? _wrap_to_pi(2 * pi / nph * (1-c)) : 0.0 for c in terminals]

            vr = vm .* cos.(va)
            vi = vm .* sin.(va)
        end

        for (idx,t) in enumerate(terminals)
            JuMP.set_start_value(var(pm, nw, :vr, id)[t], vr[idx])
            JuMP.set_start_value(var(pm, nw, :vi, id)[t], vi[idx])
        end
    end

    # apply bounds if bounded
    if bounded
        for i in ids(pm, nw, :bus)
            constraint_mc_voltage_magnitude_bounds(pm, i, nw=nw)
        end
    end
end


""
function variable_mc_bus_voltage_on_off(pm::AbstractUnbalancedACRModel; nw::Int=nw_id_default, bounded::Bool=true, kwargs...)
    variable_mc_bus_voltage_real_on_off(pm; nw=nw, bounded=bounded, kwargs...)
    variable_mc_bus_voltage_imaginary_on_off(pm; nw=nw, bounded=bounded, kwargs...)

    for id in ids(pm, nw, :bus)
        busref = ref(pm, nw, :bus, id)
        terminals = busref["terminals"]
        grounded = busref["grounded"]

        ncnd = length(terminals)

        if haskey(busref, "vr_start") && haskey(busref, "vm_start")
            vr = busref["vr_start"]
            vi = busref["vi_start"]
        else
            vm = haskey(busref, "vm_start") ? busref["vm_start"] : fill(0.0, ncnd)
            vm[.!grounded] .= 1.0

            # TODO how to do this more generally?
            nph = 3
            va = haskey(busref, "va_start") ? busref["va_start"] : [c <= nph ? _wrap_to_pi(2 * pi / nph * (1-c)) : 0.0 for c in terminals]

            vr = vm .* cos.(va)
            vi = vm .* sin.(va)
        end

        for (idx,t) in enumerate(terminals)
            JuMP.set_start_value(var(pm, nw, :vr, id)[t], vr[idx])
            JuMP.set_start_value(var(pm, nw, :vi, id)[t], vi[idx])
        end
        # apply bounds if bounded
        if bounded
            for i in ids(pm, nw, :bus)
                constraint_mc_voltage_magnitude_bounds(pm, i, nw=nw)
            end
        end
    end
end


""
function constraint_mc_switch_state_closed(pm::AbstractUnbalancedACRModel, nw::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int})
    vr_fr = var(pm, nw, :vr, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)

    vi_fr = var(pm, nw, :vi, f_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    for (idx,(fc,tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, vr_fr[fc] == vr_to[tc])
        JuMP.@constraint(pm.model, vi_fr[fc] == vi_to[tc])
    end
end


""
function constraint_mc_switch_state_on_off(pm::AbstractUnbalancedACRModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int}; relax::Bool=false)
    vr_fr = var(pm, nw, :vr, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)

    vi_fr = var(pm, nw, :vi, f_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    z = var(pm, nw, :switch_state, i)

    for (idx,(fc,tc)) in enumerate(zip(f_connections, t_connections))
        if relax
            M = 1e20
            JuMP.@constraint(pm.model, vr_fr[fc] - vr_to[tc] <=  M * (1-z))
            JuMP.@constraint(pm.model, vr_fr[fc] - vr_to[tc] >= -M * (1-z))

            JuMP.@constraint(pm.model, vi_fr[fc] - vi_to[tc] <=  M * (1-z))
            JuMP.@constraint(pm.model, vi_fr[fc] - vi_to[tc] >= -M * (1-z))
        else
            JuMP.@constraint(pm.model, z => {vr_fr[fc] == vr_to[tc]})
            JuMP.@constraint(pm.model, z => {vi_fr[fc] == vi_to[tc]})
        end
    end
end


""
function variable_mc_bus_voltage_real_on_off(pm::AbstractUnbalancedACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    vr = var(pm, nw)[:vr] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_vr_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "vr_start", t, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in ref(pm, nw, :bus)
            if haskey(bus, "vmax")
                for (idx, t) in enumerate(terminals[i])
                    set_lower_bound(vr[i][t], -bus["vmax"][idx])
                    set_upper_bound(vr[i][t],  bus["vmax"][idx])
                end
            end
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :vr, ids(pm, nw, :bus), vr)
end


""
function variable_mc_bus_voltage_imaginary_on_off(pm::AbstractUnbalancedACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    vi = var(pm, nw)[:vi] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_vi_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "vi_start", t, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in ref(pm, nw, :bus)
            if haskey(bus, "vmax")
                for (idx, t) in enumerate(terminals[i])
                    set_lower_bound(vi[i][t], -bus["vmax"][idx])
                    set_upper_bound(vi[i][t],  bus["vmax"][idx])
                end
            end
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :vi, ids(pm, nw, :bus), vi)
end


"`vmin <= vm[i] <= vmax`"
function constraint_mc_voltage_magnitude_bounds(pm::AbstractUnbalancedACRModel, nw::Int, i::Int, vmin::Vector{<:Real}, vmax::Vector{<:Real})
    @assert all(vmin .<= vmax)
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    for (idx,t) in enumerate(ref(pm, nw, :bus, i)["terminals"])
        JuMP.@constraint(pm.model, vmin[idx]^2 <= vr[t]^2 + vi[t]^2)
        if vmax[idx] < Inf
            JuMP.@constraint(pm.model, vmax[idx]^2 >= vr[t]^2 + vi[t]^2)
        end
    end
end


"bus voltage on/off constraint for load shed problem"
function constraint_mc_bus_voltage_on_off(pm::AbstractUnbalancedACRModel; nw::Int=nw_id_default, kwargs...)
    for (i,bus) in ref(pm, nw, :bus)
        constraint_mc_bus_voltage_magnitude_on_off(pm, i; nw=nw)
    end
end


"on/off bus voltage magnitude constraint"
function constraint_mc_bus_voltage_magnitude_on_off(pm::AbstractUnbalancedACRModel, nw::Int, i::Int, vmin::Vector{<:Real}, vmax::Vector{<:Real})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)
    z_voltage = var(pm, nw, :z_voltage, i)

    # TODO: non-convex constraints, look into ways to avoid in the future
    # z_voltage*vr_lb[c] <= vr[c] <= z_voltage*vr_ub[c]
    # z_voltage*vi_lb[c] <= vi[c] <= z_voltage*vi_ub[c]
    for (idx, t) in enumerate(ref(pm, nw, :bus, i)["terminals"])
        if isfinite(vmax[idx])
            JuMP.@constraint(pm.model, vr[t]^2 + vi[t]^2 <= vmax[idx]^2*z_voltage)
        end

        if isfinite(vmin[idx])
            JuMP.@constraint(pm.model, vr[t]^2 + vi[t]^2 >= vmin[idx]^2*z_voltage)
        end
    end
end


"Creates phase angle constraints at reference buses"
function constraint_mc_theta_ref(pm::AbstractUnbalancedACRModel, nw::Int, i::Int, va_ref::Vector{<:Real})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    # deal with cases first where tan(theta)==Inf or tan(theta)==0
    for (idx, t) in enumerate(ref(pm, nw, :bus, i)["terminals"])
        if va_ref[t] == pi/2
            JuMP.@constraint(pm.model, vr[t] == 0)
            JuMP.@constraint(pm.model, vi[t] >= 0)
        elseif va_ref[t] == -pi/2
            JuMP.@constraint(pm.model, vr[t] == 0)
            JuMP.@constraint(pm.model, vi[t] <= 0)
        elseif va_ref[t] == 0
            JuMP.@constraint(pm.model, vr[t] >= 0)
            JuMP.@constraint(pm.model, vi[t] == 0)
        elseif va_ref[t] == pi
            JuMP.@constraint(pm.model, vr[t] >= 0)
            JuMP.@constraint(pm.model, vi[t] == 0)
        else
            JuMP.@constraint(pm.model, vi[t] == tan(va_ref[t])*vr[t])
            # va_ref also implies a sign for vr, vi
            if 0<=va_ref[t] && va_ref[t] <= pi
                JuMP.@constraint(pm.model, vi[t] >= 0)
            else
                JuMP.@constraint(pm.model, vi[t] <= 0)
            end
        end
    end
end


""
function constraint_mc_voltage_angle_difference(pm::AbstractUnbalancedACRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})
    i, f_bus, t_bus = f_idx

    vr_fr = var(pm, nw, :vr, f_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    for (idx, (fc,tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, (vi_fr[fc] * vr_to[tc] .- vr_fr[fc] * vi_to[tc]) <= tan(angmax[idx]) * (vr_fr[fc] * vr_to[tc] .+ vi_fr[fc] * vi_to[tc]))
        JuMP.@constraint(pm.model, (vi_fr[fc] * vr_to[tc] .- vr_fr[fc] * vi_to[tc]) >= tan(angmin[idx]) * (vr_fr[fc] * vr_to[tc] .+ vi_fr[fc] * vi_to[tc]))
    end
end


""
function constraint_mc_power_balance_slack(pm::AbstractUnbalancedACRModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)
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

    ncnds = length(terminals)
    Pd = fill(0.0, ncnds)
    Qd = fill(0.0, ncnds)
    for (ld_i, connections) in bus_loads
        load = ref(pm, nw, :load, ld_i)
        for (idx, c) in enumerate(connections)
            Pd[findfirst(isequal(c), terminals)] += load["pd"][idx]
            Qd[findfirst(isequal(c), terminals)] += load["qd"][idx]
        end
    end

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx, t) for (idx, t) in enumerate(terminals) if !grounded[idx]]

    for (idx, t) in ungrounded_terminals
        crsh = sum(Gt[idx,jdx]*vr[s]-Bt[idx,jdx]*vi[s] for (jdx,s) in enumerate(terminals) if !grounded[jdx])
        cish = sum(Gt[idx,jdx]*vi[s]+Bt[idx,jdx]*vr[s] for (jdx,s) in enumerate(terminals) if !grounded[jdx])

        cp = JuMP.@constraint(pm.model,
            sum(p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(pt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
            ==
            sum(pg[g][t] for (g, conns) in bus_gens if t in conns)
            - sum(ps[s][t] for (s, conns) in bus_storage if t in conns)
            - Pd[idx]
            - (vr.*(Gt*vr-Bt*vi) + vi.*(Gt*vi+Bt*vr))
            + p_slack[t]
        )
        push!(cstr_p, cp)

        cq = JuMP.@constraint(pm.model,
            sum(q[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(qsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(qt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
            ==
              sum(qg[g][t] for (g, conns) in bus_gens if t in conns)
            + sum(qs[s][t] for (s, conns) in bus_storage if t in conns)
            - Qd[idx]
            + ( vr[t] * cish - vi[t] * crsh)
            + q_slack[t]
        )
        push!(cstr_q, cq)
    end

    con(pm, nw, :lam_kcl_r)[i] = isa(cstr_p, Array) ? cstr_p : [cstr_p]
    con(pm, nw, :lam_kcl_i)[i] = isa(cstr_q, Array) ? cstr_q : [cstr_q]

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_mc_power_balance_simple(pm::AbstractUnbalancedACRModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)
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


    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    ncnds = length(terminals)
    Pd = fill(0.0, ncnds)
    Qd = fill(0.0, ncnds)
    for (ld_i, connections) in bus_loads
        load = ref(pm, nw, :load, ld_i)
        for (idx, c) in enumerate(connections)
            Pd[findfirst(isequal(c), terminals)] += load["pd"][idx]
            Qd[findfirst(isequal(c), terminals)] += load["qd"][idx]
        end
    end

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
            - Pd[idx]
            + sum(-vr[t] * sum(Gt[idx,jdx]*vr[s]-Bt[idx,jdx]*vi[s] for (jdx,s) in ungrounded_terminals)
                  -vi[t] * sum(Gt[idx,jdx]*vi[s]+Bt[idx,jdx]*vr[s] for (jdx,s) in ungrounded_terminals))
        )
        push!(cstr_p, cp)

        cq = JuMP.@constraint(pm.model,
              sum(q[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(qsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(qt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
            ==
              sum(qg[g][t] for (g, conns) in bus_gens if t in conns)
            - sum(qs[s][t] for (s, conns) in bus_storage if t in conns)
            - Qd[idx]
            + ( vr[t] * sum(Gt[idx,jdx]*vi[s]+Bt[idx,jdx]*vr[s] for (jdx,s) in ungrounded_terminals)
               -vi[t] * sum(Gt[idx,jdx]*vr[s]-Bt[idx,jdx]*vi[s] for (jdx,s) in ungrounded_terminals) )
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
function constraint_mc_power_balance(pm::AbstractUnbalancedACRModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)
    p    = get(var(pm, nw), :p,      Dict()); _check_var_keys(p,   bus_arcs,       "active power",   "branch")
    q    = get(var(pm, nw), :q,      Dict()); _check_var_keys(q,   bus_arcs,       "reactive power", "branch")
    pg   = get(var(pm, nw), :pg_bus, Dict()); _check_var_keys(pg,  bus_gens,       "active power",   "generator")
    qg   = get(var(pm, nw), :qg_bus, Dict()); _check_var_keys(qg,  bus_gens,       "reactive power", "generator")
    ps   = get(var(pm, nw), :ps,     Dict()); _check_var_keys(ps,  bus_storage,    "active power",   "storage")
    qs   = get(var(pm, nw), :qs,     Dict()); _check_var_keys(qs,  bus_storage,    "reactive power", "storage")
    psw  = get(var(pm, nw), :psw,    Dict()); _check_var_keys(psw, bus_arcs_sw,    "active power",   "switch")
    qsw  = get(var(pm, nw), :qsw,    Dict()); _check_var_keys(qsw, bus_arcs_sw,    "reactive power", "switch")
    pt   = get(var(pm, nw), :pt,     Dict()); _check_var_keys(pt,  bus_arcs_trans, "active power",   "transformer")
    qt   = get(var(pm, nw), :qt,     Dict()); _check_var_keys(qt,  bus_arcs_trans, "reactive power", "transformer")
    pd   = get(var(pm, nw), :pd_bus, Dict()); _check_var_keys(pd,  bus_loads,      "active power",   "load")
    qd   = get(var(pm, nw), :qd_bus, Dict()); _check_var_keys(pd,  bus_loads,      "reactive power", "load")

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    # pd/qd can be NLexpressions, so cannot be vectorized
    for (idx, t) in ungrounded_terminals
        cp = @smart_constraint(pm.model, [p, q, pg, qg, ps, qs, psw, qsw, pt, qt, pd, qd, vr, vi],
              sum(  p[arc][t] for (arc, conns) in bus_arcs if t in conns)
            + sum(psw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
            + sum( pt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
            ==
              sum(pg[gen][t] for (gen, conns) in bus_gens if t in conns)
            - sum(ps[strg][t] for (strg, conns) in bus_storage if t in conns)
            - sum(pd[load][t] for (load, conns) in bus_loads if t in conns)
            + ( -vr[t] * sum(Gt[idx,jdx]*vr[u]-Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals)
                -vi[t] * sum(Gt[idx,jdx]*vi[u]+Bt[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals)
            )
        )
        push!(cstr_p, cp)

        cq = @smart_constraint(pm.model, [p, q, pg, qg, ps, qs, psw, qsw, pt, qt, pd, qd, vr, vi],
              sum(  q[arc][t] for (arc, conns) in bus_arcs if t in conns)
            + sum(qsw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
            + sum( qt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
            ==
              sum(qg[gen][t] for (gen, conns) in bus_gens if t in conns)
            - sum(qd[load][t] for (load, conns) in bus_loads if t in conns)
            - sum(qs[strg][t] for (strg, conns) in bus_storage if t in conns)
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


""
function constraint_mc_power_balance_shed(pm::AbstractUnbalancedACRModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)
    p    = get(var(pm, nw), :p,      Dict()); _check_var_keys(p,   bus_arcs,       "active power",   "branch")
    q    = get(var(pm, nw), :q,      Dict()); _check_var_keys(q,   bus_arcs,       "reactive power", "branch")
    pg   = get(var(pm, nw), :pg, Dict()); _check_var_keys(pg,  bus_gens,       "active power",   "generator")
    qg   = get(var(pm, nw), :qg, Dict()); _check_var_keys(qg,  bus_gens,       "reactive power", "generator")
    ps   = get(var(pm, nw), :ps,     Dict()); _check_var_keys(ps,  bus_storage,    "active power",   "storage")
    qs   = get(var(pm, nw), :qs,     Dict()); _check_var_keys(qs,  bus_storage,    "reactive power", "storage")
    psw  = get(var(pm, nw), :psw,    Dict()); _check_var_keys(psw, bus_arcs_sw,    "active power",   "switch")
    qsw  = get(var(pm, nw), :qsw,    Dict()); _check_var_keys(qsw, bus_arcs_sw,    "reactive power", "switch")
    pt   = get(var(pm, nw), :pt,     Dict()); _check_var_keys(pt,  bus_arcs_trans, "active power",   "transformer")
    qt   = get(var(pm, nw), :qt,     Dict()); _check_var_keys(qt,  bus_arcs_trans, "reactive power", "transformer")

    zd = var(pm, nw, :z_demand)
    z_shunt  = var(pm, nw, :z_shunt)  # TODO add support for z_shunt in power balance shed
    zg = haskey(var(pm, nw), :z_gen) ? var(pm, nw, :z_gen) : Dict(i => 1.0 for i in ids(pm, nw, :gen))
    zs = haskey(var(pm, nw), :z_storage) ? var(pm, nw, :z_storage) : Dict(i => 1.0 for i in ids(pm, nw, :storage))


    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    # pd/qd can be NLexpressions, so cannot be vectorized
    for (idx, t) in ungrounded_terminals
        cp = @smart_constraint(pm.model, [p, pg, ps, psw, pt],
              sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
            + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
            - sum(pg[g][t]*zg[g] for (g, conns) in bus_gens if t in conns)
            + sum(ps[s][t]*zs[s] for (s, conns) in bus_storage if t in conns)
            + sum(ref(pm, nw, :load, d, "pd")[findfirst(isequal(t), conns)]*zd[d] for (d, conns) in bus_loads if t in conns)
            + (+vr[t] * sum(Gt[idx,jdx]*vr[u]-Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals)
               +vi[t] * sum(Gt[idx,jdx]*vi[u]+Bt[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals))
            ==
            0.0
        )
        push!(cstr_p, cp)

        cq = @smart_constraint(pm.model, [q, qg, qs, qsw, qt],
              sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
            + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
            - sum(qg[g][t]*zg[g] for (g, conns) in bus_gens if t in conns)
            + sum(qs[s][t]*zs[s] for (s, conns) in bus_storage if t in conns)
            + sum(ref(pm, nw, :load, d, "qd")[findfirst(isequal(t), conns)]*zd[d] for (d, conns) in bus_loads if t in conns)
            + (-vr[t] * sum(Gt[idx,jdx]*vi[u]+Bt[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals)
               +vi[t] * sum(Gt[idx,jdx]*vr[u]-Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals))
            ==
            0.0
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


"""
Creates Ohms constraints

s_fr = v_fr.*conj(Y*(v_fr-v_to))
s_fr = (vr_fr+im*vi_fr).*(G-im*B)*([vr_fr-vr_to]-im*[vi_fr-vi_to])
s_fr = (vr_fr+im*vi_fr).*([G*vr_fr-G*vr_to-B*vi_fr+B*vi_to]-im*[G*vi_fr-G*vi_to+B*vr_fr-B*vr_to])
"""
function constraint_mc_ohms_yt_from(pm::AbstractUnbalancedACRModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})
    p_fr  = [var(pm, nw, :p, f_idx)[t] for t in f_connections]
    q_fr  = [var(pm, nw, :q, f_idx)[t] for t in f_connections]
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
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_mc_ohms_yt_to(pm::AbstractUnbalancedACRModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix, B::Matrix, G_to::Matrix, B_to::Matrix)
    constraint_mc_ohms_yt_from(pm, nw, t_bus, f_bus, t_idx, f_idx, t_connections, f_connections, G, B, G_to, B_to)
end


""
function constraint_mc_load_power_wye(pm::AbstractUnbalancedACRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    # if constant power load
    if all(alpha.==0) && all(beta.==0)
        pd_bus = a
        qd_bus = b
    else
        pd_bus = Vector{JuMP.NonlinearExpression}([])
        qd_bus = Vector{JuMP.NonlinearExpression}([])

        for (idx,c) in enumerate(connections)
            crd = JuMP.@NLexpression(pm.model, a[idx]*vr[c]*(vr[c]^2+vi[c]^2)^(alpha[idx]/2-1)+b[idx]*vi[c]*(vr[c]^2+vi[c]^2)^(beta[idx]/2 -1))
            cid = JuMP.@NLexpression(pm.model, a[idx]*vi[c]*(vr[c]^2+vi[c]^2)^(alpha[idx]/2-1)-b[idx]*vr[c]*(vr[c]^2+vi[c]^2)^(beta[idx]/2 -1))

            push!(pd_bus, JuMP.@NLexpression(pm.model,  vr[c]*crd+vi[c]*cid))
            push!(qd_bus, JuMP.@NLexpression(pm.model, -vr[c]*cid+vi[c]*crd))
        end
    end

    pd_bus = JuMP.Containers.DenseAxisArray(pd_bus, connections)
    qd_bus = JuMP.Containers.DenseAxisArray(qd_bus, connections)

    var(pm, nw, :pd_bus)[id] = pd_bus
    var(pm, nw, :qd_bus)[id] = qd_bus

    if report
        sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        pd = Vector{JuMP.NonlinearExpression}([])
        qd = Vector{JuMP.NonlinearExpression}([])

        for (idx,c) in enumerate(connections)
            push!(pd, JuMP.@NLexpression(pm.model, a[idx]*(vr[c]^2+vi[c]^2)^(alpha[idx]/2) ))
            push!(qd, JuMP.@NLexpression(pm.model, b[idx]*(vr[c]^2+vi[c]^2)^(beta[idx]/2)  ))
        end
        sol(pm, nw, :load, id)[:pd] = JuMP.Containers.DenseAxisArray(pd, connections)
        sol(pm, nw, :load, id)[:qd] = JuMP.Containers.DenseAxisArray(qd, connections)
    end
end


""
function constraint_mc_load_power_delta(pm::AbstractUnbalancedACRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    nph = length(a)
    @assert nph == 3 "only phases == 3 delta loads are currently supported"

    prev = Dict(c=>connections[(idx+nph-2)%nph+1] for (idx,c) in enumerate(connections))
    next = Dict(c=>connections[idx%nph+1] for (idx,c) in enumerate(connections))

    vrd = Dict()
    vid = Dict()
    for (idx, c) in enumerate(connections)
        vrd[c] = JuMP.@NLexpression(pm.model, vr[c]-vr[next[c]])
        vid[c] = JuMP.@NLexpression(pm.model, vi[c]-vi[next[c]])
    end

    crd = Dict()
    cid = Dict()
    for (idx, c) in enumerate(connections)
        crd[c] = JuMP.@NLexpression(pm.model, a[idx]*vrd[c]*(vrd[c]^2+vid[c]^2)^(alpha[idx]/2-1)+b[idx]*vid[c]*(vrd[c]^2+vid[c]^2)^(beta[idx]/2 -1))
        cid[c] = JuMP.@NLexpression(pm.model, a[idx]*vid[c]*(vrd[c]^2+vid[c]^2)^(alpha[idx]/2-1)-b[idx]*vrd[c]*(vrd[c]^2+vid[c]^2)^(beta[idx]/2 -1))
    end

    crd_bus = Dict()
    cid_bus = Dict()
    for (idx, c) in enumerate(connections)
        crd_bus[c] = JuMP.@NLexpression(pm.model, crd[c]-crd[prev[c]])
        cid_bus[c] = JuMP.@NLexpression(pm.model, cid[c]-cid[prev[c]])
    end

    pd_bus = Vector{JuMP.NonlinearExpression}([])
    qd_bus = Vector{JuMP.NonlinearExpression}([])
    for (idx,c) in enumerate(connections)
        push!(pd_bus, JuMP.@NLexpression(pm.model,  vr[c]*crd_bus[c]+vi[c]*cid_bus[c]))
        push!(qd_bus, JuMP.@NLexpression(pm.model, -vr[c]*cid_bus[c]+vi[c]*crd_bus[c]))
    end

    pd_bus = JuMP.Containers.DenseAxisArray(pd_bus, connections)
    qd_bus = JuMP.Containers.DenseAxisArray(qd_bus, connections)

    var(pm, nw, :pd_bus)[id] = pd_bus
    var(pm, nw, :qd_bus)[id] = qd_bus

    if report
        sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        pd = Vector{JuMP.NonlinearExpression}([])
        qd = Vector{JuMP.NonlinearExpression}([])
        for (idx,c) in enumerate(connections)
            push!(pd, JuMP.@NLexpression(pm.model, a[idx]*(vrd[c]^2+vid[c]^2)^(alpha[idx]/2) ))
            push!(qd, JuMP.@NLexpression(pm.model, b[idx]*(vrd[c]^2+vid[c]^2)^(beta[idx]/2)  ))
        end
        sol(pm, nw, :load, id)[:pd] = JuMP.Containers.DenseAxisArray(pd, connections)
        sol(pm, nw, :load, id)[:qd] = JuMP.Containers.DenseAxisArray(qd, connections)
    end
end


"`vm[i] == vmref`"
function constraint_mc_voltage_magnitude_only(pm::AbstractUnbalancedACRModel, nw::Int, i::Int, vm_ref::Vector{<:Real})
    vr = [var(pm, nw, :vr, i)[t] for t in ref(pm, nw, :bus, i)["terminals"]]
    vi = [var(pm, nw, :vi, i)[t] for t in ref(pm, nw, :bus, i)["terminals"]]

    JuMP.@constraint(pm.model, vr.^2 + vi.^2  .== vm_ref.^2)
end


""
function constraint_mc_generator_power_delta(pm::AbstractUnbalancedACRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    pg = var(pm, nw, :pg, id)
    qg = var(pm, nw, :qg, id)

    nph = length(pmin)
    @assert nph == 3 "only phases == 3 delta generators are currently supported"

    prev = Dict(c=>connections[(idx+nph-2)%nph+1] for (idx,c) in enumerate(connections))
    next = Dict(c=>connections[idx%nph+1] for (idx,c) in enumerate(connections))

    vrg = Dict()
    vig = Dict()
    for c in connections
        vrg[c] = JuMP.@NLexpression(pm.model, vr[c]-vr[next[c]])
        vig[c] = JuMP.@NLexpression(pm.model, vi[c]-vi[next[c]])
    end

    crg = Dict()
    cig = Dict()
    for c in connections
        crg[c] = JuMP.@NLexpression(pm.model, (pg[c]*vrg[c]+qg[c]*vig[c])/(vrg[c]^2+vig[c]^2) )
        cig[c] = JuMP.@NLexpression(pm.model, (pg[c]*vig[c]-qg[c]*vrg[c])/(vrg[c]^2+vig[c]^2) )
    end

    crg_bus = Dict()
    cig_bus = Dict()
    for c in connections
        crg_bus[c] = JuMP.@NLexpression(pm.model, crg[c]-crg[prev[c]])
        cig_bus[c] = JuMP.@NLexpression(pm.model, cig[c]-cig[prev[c]])
    end

    pg_bus = Vector{JuMP.NonlinearExpression}([])
    qg_bus = Vector{JuMP.NonlinearExpression}([])
    for (idx,c) in enumerate(connections)
        push!(pg_bus, JuMP.@NLexpression(pm.model,  vr[c]*crg_bus[c]+vi[c]*cig_bus[c]))
        push!(qg_bus, JuMP.@NLexpression(pm.model, -vr[c]*cig_bus[c]+vi[c]*crg_bus[c]))
    end
    pd_bus = JuMP.Containers.DenseAxisArray(pg_bus, connections)
    qd_bus = JuMP.Containers.DenseAxisArray(qg_bus, connections)

    var(pm, nw, :pg_bus)[id] = pg_bus
    var(pm, nw, :qg_bus)[id] = qg_bus

    if report
        sol(pm, nw, :gen, id)[:pg_bus] = pg_bus
        sol(pm, nw, :gen, id)[:qg_bus] = qg_bus
    end
end


"This function adds all constraints required to model a two-winding, wye-wye connected transformer."
function constraint_mc_transformer_power_yy(pm::AbstractUnbalancedACRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    transformer = ref(pm, nw, :transformer, trans_id)

    vr_fr = var(pm, nw, :vr, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    p_fr = [var(pm, nw, :pt, f_idx)[c] for c in f_connections]
    p_to = [var(pm, nw, :pt, t_idx)[c] for c in t_connections]
    q_fr = [var(pm, nw, :qt, f_idx)[c] for c in f_connections]
    q_to = [var(pm, nw, :qt, t_idx)[c] for c in t_connections]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))]

    for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))
        if tm_fixed[idx]
            JuMP.@constraint(pm.model, vr_fr[fc] == pol*tm_scale*tm[idx]*vr_to[tc])
            JuMP.@constraint(pm.model, vi_fr[fc] == pol*tm_scale*tm[idx]*vi_to[tc])
        else
            # transformer taps without regcontrol, tap variable not required in regcontrol formulation
            JuMP.@constraint(pm.model, vr_fr[fc] == pol*tm_scale*tm[idx]*vr_to[tc])
            JuMP.@constraint(pm.model, vi_fr[fc] == pol*tm_scale*tm[idx]*vi_to[tc])

            # with regcontrol
            if haskey(transformer,"controls")
                v_ref = transformer["controls"]["vreg"][idx]
                δ = transformer["controls"]["band"][idx]
                r = transformer["controls"]["r"][idx]
                x = transformer["controls"]["x"][idx]

                # (cr+jci) = (p-jq)/(vr-j⋅vi)
                cr = JuMP.@NLexpression(pm.model, ( p_to[idx]*vr_to[tc] + q_to[idx]*vi_to[tc])/(vr_to[tc]^2+vi_to[tc]^2))
                ci = JuMP.@NLexpression(pm.model, (-q_to[idx]*vr_to[tc] + p_to[idx]*vi_to[tc])/(vr_to[tc]^2+vi_to[tc]^2))
                # v_drop = (cr+jci)⋅(r+jx)
                vr_drop = JuMP.@NLexpression(pm.model, r*cr-x*ci)
                vi_drop = JuMP.@NLexpression(pm.model, r*ci+x*cr)

                # (v_ref-δ)^2 ≤ (vr_fr-vr_drop)^2 + (vi_fr-vi_drop)^2 ≤ (v_ref+δ)^2
                # (vr_fr^2 + vi_fr^2)/1.1^2 ≤ (vr_to^2 + vi_to^2) ≤ (vr_fr^2 + vi_fr^2)/0.9^2
                JuMP.@NLconstraint(pm.model, (vr_fr[fc]-vr_drop)^2 + (vi_fr[fc]-vi_drop)^2 ≥ (v_ref - δ)^2)
                JuMP.@NLconstraint(pm.model, (vr_fr[fc]-vr_drop)^2 + (vi_fr[fc]-vi_drop)^2 ≤ (v_ref + δ)^2)
                JuMP.@constraint(pm.model, (vr_fr[fc]^2 + vi_fr[fc]^2)/1.1^2 ≤ vr_to[tc]^2 + vi_to[tc]^2)
                JuMP.@constraint(pm.model, (vr_fr[fc]^2 + vi_fr[fc]^2)/0.9^2 ≥ vr_to[tc]^2 + vi_to[tc]^2)
            end
        end
    end

    JuMP.@constraint(pm.model, p_fr + p_to .== 0)
    JuMP.@constraint(pm.model, q_fr + q_to .== 0)
end


"This function adds all constraints required to model a two-winding, delta-wye connected transformer."
function constraint_mc_transformer_power_dy(pm::AbstractUnbalancedACRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    vr_p_fr = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vr_p_to = [var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_p_fr = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]
    vi_p_to = [var(pm, nw, :vi, t_bus)[c] for c in t_connections]

    nph = length(tm_set)
    @assert length(f_connections) == length(t_connections) && nph == 3 "only phases == 3 dy transformers are currently supported"

    M = _get_delta_transformation_matrix(nph)

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[fc] for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))]

    # introduce auxialiary variable vd = Md*v_fr
    vrd = M*vr_p_fr
    vid = M*vi_p_fr

    JuMP.@constraint(pm.model, vrd .== (pol*tm_scale)*tm.*vr_p_to)
    JuMP.@constraint(pm.model, vid .== (pol*tm_scale)*tm.*vi_p_to)

    p_fr = var(pm, nw, :pt, f_idx)
    p_to = var(pm, nw, :pt, t_idx)
    q_fr = var(pm, nw, :qt, f_idx)
    q_to = var(pm, nw, :qt, t_idx)

    id_re = Array{Any,1}(undef, nph)
    id_im = Array{Any,1}(undef, nph)
    # s/v      = (p+jq)/|v|^2*conj(v)
    #          = (p+jq)/|v|*(cos(va)-j*sin(va))
    # Re(s/v)  = (p*cos(va)+q*sin(va))/|v|
    # -Im(s/v) = -(q*cos(va)-p*sin(va))/|v|
    for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))
        # id = conj(s_to/v_to)./tm
        id_re[idx] = JuMP.@NLexpression(pm.model, (p_to[tc]*vr_p_to[idx]+q_to[tc]*vi_p_to[idx])/(tm_scale*tm[idx]*pol)/(vr_p_to[idx]^2+vi_p_to[idx]^2))
        id_im[idx] = JuMP.@NLexpression(pm.model, (p_to[tc]*vi_p_to[idx]-q_to[tc]*vr_p_to[idx])/(tm_scale*tm[idx]*pol)/(vr_p_to[idx]^2+vi_p_to[idx]^2))
    end
    for (idx,(fc,tc)) in enumerate(zip(f_connections, t_connections))
        jdx = (idx-1+nph-1)%nph+1
        # s_fr  = v_fr*conj(i_fr)
        #       = v_fr*conj(id[q]-id[p])
        #       = v_fr*(id_re[q]-j*id_im[q]-id_re[p]+j*id_im[p])
        JuMP.@NLconstraint(pm.model, p_fr[fc] ==
             vr_p_fr[idx]*( id_re[jdx]-id_re[idx])
            -vi_p_fr[idx]*(-id_im[jdx]+id_im[idx])
        )
        JuMP.@NLconstraint(pm.model, q_fr[fc] ==
             vr_p_fr[idx]*(-id_im[jdx]+id_im[idx])
            +vi_p_fr[idx]*( id_re[jdx]-id_re[idx])
        )
    end
end


""
function constraint_mc_storage_losses(pm::AbstractUnbalancedACRModel, i::Int; nw::Int=nw_id_default, kwargs...)
    storage = ref(pm, nw, :storage, i)

    vr  = var(pm, nw,  :vr, storage["storage_bus"])
    vi  = var(pm, nw,  :vi, storage["storage_bus"])
    ps  = var(pm, nw,  :ps, i)
    qs  = var(pm, nw,  :qs, i)
    sc  = var(pm, nw,  :sc, i)
    sd  = var(pm, nw,  :sd, i)
    qsc = var(pm, nw, :qsc, i)

    p_loss = storage["p_loss"]
    q_loss = storage["q_loss"]
    r = storage["r"]
    x = storage["x"]

    JuMP.@NLconstraint(pm.model,
        sum(ps[c] for c in storage["connections"]) + (sd - sc)
        ==
        p_loss + sum(r[idx]*(ps[c]^2 + qs[c]^2)/(vr[c]^2 + vi[c]^2) for (idx,c) in enumerate(storage["connections"]))
    )

    JuMP.@NLconstraint(pm.model,
        sum(qs[c] for c in storage["connections"])
        ==
        qsc + q_loss + sum(x[idx]*(ps[c]^2 + qs[c]^2)/(vr[c]^2 + vi[c]^2) for (idx,c) in enumerate(storage["connections"]))
    )
end


""
function constraint_storage_losses(pm::AbstractUnbalancedACRModel, n::Int, i, bus, r, x, p_loss, q_loss; conductors=[1])
    vr = var(pm, n, :vr, bus)
    vi = var(pm, n, :vi, bus)
    ps = var(pm, n, :ps, i)
    qs = var(pm, n, :qs, i)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    qsc = var(pm, n, :qsc, i)

    JuMP.@NLconstraint(pm.model,
        sum(ps[c] for c in conductors) + (sd - sc)
        ==
        p_loss + sum(r[c]*(ps[c]^2 + qs[c]^2)/(vr[c]^2 + vi[c]^2) for c in conductors)
    )

    JuMP.@NLconstraint(pm.model,
        sum(qs[c] for c in conductors)
        ==
        qsc + q_loss + sum(x[c]*(ps[c]^2 + qs[c]^2)/(vr[c]^2 + vi[c]^2) for c in conductors)
    )
end
