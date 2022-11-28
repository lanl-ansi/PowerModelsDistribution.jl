# The model in rectangular coordinates is linearized around an initial operating point using a
# first order Taylor approximation (FOT) method


"""
    variable_mc_bus_voltage(pm::FOTRUPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

Voltage variables are defined in rectangular coordinates similar to ACRUPowerModel.
An initial operating point is specified for linearization similar to FBSUBFPowerModel.
"""
function variable_mc_bus_voltage(pm::FOTRUPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_bus_voltage_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_bus_voltage_imaginary(pm; nw=nw, bounded=bounded, report=report)

    # initial operating point for linearization (using flat-start)
    var(pm, nw)[:vr0] = Dict{Int,Vector{Float64}}()
    var(pm, nw)[:vi0] = Dict{Int,Vector{Float64}}()

    for id in ids(pm, nw, :bus)
        busref = ref(pm, nw, :bus, id)
        terminals = busref["terminals"]
        grounded = busref["grounded"]

        ncnd = length(terminals)

        vm_start = fill(1.0, 3)
        for t in 1:3
            if t in terminals
                vmax = busref["vmax"][findfirst(isequal(t), terminals)]
                vm_start[t] = min(vm_start[t], vmax)

                vmin = busref["vmin"][findfirst(isequal(t), terminals)]
                vm_start[t] = max(vm_start[t], vmin)
            end
        end
        default_va = [[_wrap_to_pi(2 * pi / 3 * (1-t)) for t in 1:3]..., zeros(length(terminals))...][terminals]
        vm = haskey(busref, "vm_start") ? busref["vm_start"] : haskey(busref, "vm") ? busref["vm"] : [vm_start..., fill(0.0, ncnd)...][terminals]
        va = haskey(busref, "va_start") ? busref["va_start"] : haskey(busref, "va") ? busref["va"] : default_va

        vr = vm.*cos.(va)
        vi = vm.*sin.(va)

        JuMP.set_start_value.(var(pm, nw, :vr, id), vr)
        JuMP.set_start_value.(var(pm, nw, :vi, id), vi)

        # TODO: update initial operating point with warm-start (causes infeasbility if not flat start)
        var(pm, nw, :vr0)[id] = fill(1.0, ncnd) .* cos.(default_va) # vr
        var(pm, nw, :vi0)[id] = fill(1.0, ncnd) .* sin.(default_va) # vi
    end

    # apply bounds if bounded
    if bounded
        for i in ids(pm, nw, :bus)
            constraint_mc_voltage_magnitude_bounds(pm, i, nw=nw)
        end
    end
end


"""
    constraint_mc_voltage_magnitude_bounds(pm::FOTRUPowerModel, nw::Int, i::Int, vmin::Vector{<:Real}, vmax::Vector{<:Real})

Linearized voltage magnitude limits similar to FBSUBFPowerModel.
"""
function constraint_mc_voltage_magnitude_bounds(pm::FOTRUPowerModel, nw::Int, i::Int, vmin::Vector{<:Real}, vmax::Vector{<:Real})
    @assert all(vmin .<= vmax)
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)
    vr0 = var(pm, nw, :vr0, i)
    vi0 = var(pm, nw, :vi0, i)

    for (idx,t) in enumerate(ref(pm, nw, :bus, i)["terminals"])
        # linearized lower voltage magnitude limits at reference bus
        if ref(pm, nw, :bus, i)["bus_type"] == 3.0 && vmax[idx] == vmin[idx]
            vr_ref = ref(pm, nw, :bus, i)["vm"][idx]*cos(ref(pm, nw, :bus, i)["va"][idx])
            vi_ref = ref(pm, nw, :bus, i)["vm"][idx]*sin(ref(pm, nw, :bus, i)["va"][idx])
            JuMP.@constraint(pm.model, 2*vr[t]*vr_ref + 2*vi[t]*vi_ref - vr_ref^2 - vi_ref^2 >= vmin[idx]^2)
        # linearized lower voltage magnitude limits at at all other buses
        else
            JuMP.@constraint(pm.model, 2*vr[t]*vr0[idx] + 2*vi[t]*vi0[idx] - vr0[idx]^2 - vi0[idx]^2 >= vmin[idx]^2)
        end
        # Outer approximation of upper voltage magnitude limits
        if vmax[idx] < Inf
            JuMP.@constraint(pm.model, -vmax[idx] <= vr[t])
            JuMP.@constraint(pm.model,  vmax[idx] >= vr[t])
            JuMP.@constraint(pm.model, -vmax[idx] <= vi[t])
            JuMP.@constraint(pm.model,  vmax[idx] >= vi[t])
            JuMP.@constraint(pm.model, -sqrt(2)*vmax[idx] <= vr[t] + vi[t])
            JuMP.@constraint(pm.model,  sqrt(2)*vmax[idx] >= vr[t] + vi[t])
            JuMP.@constraint(pm.model, -sqrt(2)*vmax[idx] <= vr[t] - vi[t])
            JuMP.@constraint(pm.model,  sqrt(2)*vmax[idx] >= vr[t] - vi[t])
        end
    end
end


"""
    constraint_mc_voltage_angle_difference(pm::FOTRUPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})

Nothing to do, this model ignores angle difference constraints"
"""
function constraint_mc_voltage_angle_difference(pm::FOTRUPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})
end


@doc raw"""
    constraint_mc_power_balance(pm::FOTRUPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})

Power balance constraints similar to ACRUPowerModel with shunt current linearized around initial operating point.

```math
\begin{align}
&\text{Initial operating point: }  (v_{r0} + j ⋅ v_{i0})\\
& v_{r} ⋅ v_{i} = v_{r0} ⋅ v_{i0} + v_{r} ⋅ v_{i0} + v_{r0} ⋅ v_{i} - 2 ⋅ v_{r0} ⋅ v_{i0}
\end{align}
```
"""
function constraint_mc_power_balance(pm::FOTRUPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)
    vr0 = var(pm, nw, :vr0, i)
    vi0 = var(pm, nw, :vi0, i)

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
        cp = JuMP.@constraint(pm.model,
              sum(  p[arc][t] for (arc, conns) in bus_arcs if t in conns)
            + sum(psw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
            + sum( pt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
            ==
              sum(pg[gen][t] for (gen, conns) in bus_gens if t in conns)
            - sum(ps[strg][t] for (strg, conns) in bus_storage if t in conns)
            - sum(pd[load][t] for (load, conns) in bus_loads if t in conns)
            + ( -sum(Gt[idx,jdx]*vr0[idx]*vr0[jdx]-Bt[idx,jdx]*vr0[idx]*vi0[jdx] for (jdx,u) in ungrounded_terminals)
                -sum(Gt[idx,jdx]*vi0[idx]*vi0[jdx]+Bt[idx,jdx]*vi0[idx]*vr0[jdx] for (jdx,u) in ungrounded_terminals)
            )
            + ( -sum(Gt[idx,jdx]*(vr[t]*vr0[jdx]+vr0[idx]*vr[u]-2*vr0[idx]*vr0[jdx])-Bt[idx,jdx]*(vr[t]*vi0[jdx]+vr0[idx]*vi[u]-2*vr0[idx]*vi0[jdx]) for (jdx,u) in ungrounded_terminals)
                -sum(Gt[idx,jdx]*(vi[t]*vi0[jdx]+vi0[idx]*vi[u]-2*vi0[idx]*vi0[jdx])+Bt[idx,jdx]*(vi[t]*vr0[jdx]+vi0[idx]*vr[u]-2*vi0[idx]*vr0[jdx]) for (jdx,u) in ungrounded_terminals)
            )
        )
        push!(cstr_p, cp)

        cq = JuMP.@constraint(pm.model,
              sum(  q[arc][t] for (arc, conns) in bus_arcs if t in conns)
            + sum(qsw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
            + sum( qt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
            ==
              sum(qg[gen][t] for (gen, conns) in bus_gens if t in conns)
            - sum(qd[load][t] for (load, conns) in bus_loads if t in conns)
            - sum(qs[strg][t] for (strg, conns) in bus_storage if t in conns)
            + ( sum(Gt[idx,jdx]*vr0[idx]*vi0[jdx]+Bt[idx,jdx]*vr0[idx]*vr0[jdx] for (jdx,u) in ungrounded_terminals)
               -sum(Gt[idx,jdx]*vi0[idx]*vr0[jdx]-Bt[idx,jdx]*vi0[idx]*vi0[jdx] for (jdx,u) in ungrounded_terminals)
            )
            + ( sum(Gt[idx,jdx]*(vr[t]*vi0[jdx]+vr0[idx]*vi[u]-2*vr0[idx]*vi0[jdx])+Bt[idx,jdx]*(vr[t]*vr0[jdx]+vr0[idx]*vr[u]-2*vr0[idx]*vr0[jdx]) for (jdx,u) in ungrounded_terminals)
               -sum(Gt[idx,jdx]*(vi[t]*vr0[jdx]+vi0[idx]*vr[u]-2*vi0[idx]*vr0[jdx])-Bt[idx,jdx]*(vi[t]*vi0[jdx]+vi0[idx]*vi[u]-2*vi0[idx]*vi0[jdx]) for (jdx,u) in ungrounded_terminals)
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


@doc raw"""
    constraint_mc_power_balance_capc(pm::FOTRUPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})

Power balance constraints with capacitor control with shunt current calculated using initial operating point.

```math
\begin{align}
    & B_t = b_s ⋅ z,~~ cq_{sh} = B_t ⋅ v, \\
    &\text{FOT approximation: }  B_t ⋅ v_r ⋅ v_i = B_{t0} ⋅ v_{r0} ⋅ v_{i0} + B_{t} ⋅ v_{r0} ⋅ v_{i0} + B_{t0} ⋅ v_{r} ⋅ v_{i0} + B_{t0} ⋅ v_{r0} ⋅ v_{i} - 3 ⋅ B_{t0} ⋅ v_{r0} ⋅ v_{i0}
\end{align}
```
"""
function constraint_mc_power_balance_capc(pm::FOTRUPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)
    vr0 = var(pm, nw, :vr0, i)
    vi0 = var(pm, nw, :vi0, i)

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

    # calculate Gs, Bs
    ncnds = length(terminals)
    Gt = fill(0.0, ncnds, ncnds)
    Bt0 = fill(0.0, ncnds, ncnds)
    Bt = convert(Matrix{JuMP.AffExpr}, JuMP.@expression(pm.model, [idx=1:ncnds, jdx=1:ncnds], 0.0))
    for (val, connections) in bus_shunts
        shunt = ref(pm,nw,:shunt,val)
        for (idx,c) in enumerate(connections)
            cap_state = haskey(shunt,"controls") ? var(pm, nw, :capacitor_state, val)[c] : 1.0
            for (jdx,d) in enumerate(connections)
                Gt[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] += shunt["gs"][idx,jdx]
                Bt0[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] += shunt["bs"][idx,jdx]
                Bt[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] = JuMP.@expression(pm.model, Bt[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] + shunt["bs"][idx,jdx]*cap_state)
            end
        end
    end

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    # add constraints to model capacitor switching
    if !isempty(bus_shunts) && haskey(ref(pm, nw, :shunt, bus_shunts[1][1]), "controls")
        constraint_capacitor_on_off(pm, nw, i, bus_shunts)

        for (idx, t) in ungrounded_terminals
            cp = JuMP.@constraint(pm.model,
                sum(  p[arc][t] for (arc, conns) in bus_arcs if t in conns)
                + sum(psw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
                + sum( pt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
                ==
                sum(pg[gen][t] for (gen, conns) in bus_gens if t in conns)
                - sum(ps[strg][t] for (strg, conns) in bus_storage if t in conns)
                - sum(pd[load][t] for (load, conns) in bus_loads if t in conns)
                + ( -sum(Gt[idx,jdx]*vr0[idx]*vr0[jdx]-Bt0[idx,jdx]*vr0[idx]*vi0[jdx] for (jdx,u) in ungrounded_terminals)
                    -sum(Gt[idx,jdx]*vi0[idx]*vi0[jdx]+Bt0[idx,jdx]*vi0[idx]*vr0[jdx] for (jdx,u) in ungrounded_terminals)
                )
                + ( -sum(Gt[idx,jdx]*(vr[t]*vr0[jdx]+vr0[idx]*vr[u]-2*vr0[idx]*vr0[jdx])-(Bt[idx,jdx]*vr0[idx]*vi0[jdx]+Bt0[idx,jdx]*vr[t]*vi0[jdx]+Bt0[idx,jdx]*vr0[idx]*vi[u]-3*Bt0[idx,jdx]*vr0[idx]*vi0[jdx]) for (jdx,u) in ungrounded_terminals)
                    -sum(Gt[idx,jdx]*(vi[t]*vi0[jdx]+vi0[idx]*vi[u]-2*vi0[idx]*vi0[jdx])+(Bt[idx,jdx]*vi0[idx]*vr0[jdx]+Bt0[idx,jdx]*vi[t]*vr0[jdx]+Bt0[idx,jdx]*vi0[idx]*vr[u]-3*Bt0[idx,jdx]*vi0[idx]*vr0[jdx]) for (jdx,u) in ungrounded_terminals)
                )
            )
            push!(cstr_p, cp)

            cq = JuMP.@constraint(pm.model,
                sum(  q[arc][t] for (arc, conns) in bus_arcs if t in conns)
                + sum(qsw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
                + sum( qt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
                ==
                sum(qg[gen][t] for (gen, conns) in bus_gens if t in conns)
                - sum(qd[load][t] for (load, conns) in bus_loads if t in conns)
                - sum(qs[strg][t] for (strg, conns) in bus_storage if t in conns)
                + ( sum(Gt[idx,jdx]*vr0[idx]*vi0[jdx]+Bt0[idx,jdx]*vr0[idx]*vr0[jdx] for (jdx,u) in ungrounded_terminals)
                   -sum(Gt[idx,jdx]*vi0[idx]*vr0[jdx]-Bt0[idx,jdx]*vi0[idx]*vi0[jdx] for (jdx,u) in ungrounded_terminals)
                )
                + ( sum(Gt[idx,jdx]*(vr[t]*vi0[jdx]+vr0[idx]*vi[u]-2*vr0[idx]*vi0[jdx])+(Bt[idx,jdx]*vr0[idx]*vr0[jdx]+Bt0[idx,jdx]*vr[t]*vr0[jdx]+Bt0[idx,jdx]*vr0[idx]*vr[u]-3*Bt0[idx,jdx]*vr0[idx]*vr0[jdx]) for (jdx,u) in ungrounded_terminals)
                   -sum(Gt[idx,jdx]*(vi[t]*vr0[jdx]+vi0[idx]*vr[u]-2*vi0[idx]*vr0[jdx])-(Bt[idx,jdx]*vi0[idx]*vi0[jdx]+Bt0[idx,jdx]*vi[t]*vi0[jdx]+Bt0[idx,jdx]*vi0[idx]*vi[u]-3*Bt0[idx,jdx]*vi0[idx]*vi0[jdx]) for (jdx,u) in ungrounded_terminals)
                )
            )
            push!(cstr_q, cq)
        end
    else
        for (idx, t) in ungrounded_terminals
            cp = JuMP.@constraint(pm.model,
                sum(  p[arc][t] for (arc, conns) in bus_arcs if t in conns)
                + sum(psw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
                + sum( pt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
                ==
                sum(pg[gen][t] for (gen, conns) in bus_gens if t in conns)
                - sum(ps[strg][t] for (strg, conns) in bus_storage if t in conns)
                - sum(pd[load][t] for (load, conns) in bus_loads if t in conns)
                + ( -sum(Gt[idx,jdx]*vr0[idx]*vr0[jdx]-Bt[idx,jdx]*vr0[idx]*vi0[jdx] for (jdx,u) in ungrounded_terminals)
                    -sum(Gt[idx,jdx]*vi0[idx]*vi0[jdx]+Bt[idx,jdx]*vi0[idx]*vr0[jdx] for (jdx,u) in ungrounded_terminals)
                )
                + ( -sum(Gt[idx,jdx]*(vr[t]*vr0[jdx]+vr0[idx]*vr[u]-2*vr0[idx]*vr0[jdx])-Bt[idx,jdx]*(vr[t]*vi0[jdx]+vr0[idx]*vi[u]-2*vr0[idx]*vi0[jdx]) for (jdx,u) in ungrounded_terminals)
                    -sum(Gt[idx,jdx]*(vi[t]*vi0[jdx]+vi0[idx]*vi[u]-2*vi0[idx]*vi0[jdx])+Bt[idx,jdx]*(vi[t]*vr0[jdx]+vi0[idx]*vr[u]-2*vi0[idx]*vr0[jdx]) for (jdx,u) in ungrounded_terminals)
                )
            )
            push!(cstr_p, cp)

            cq = JuMP.@constraint(pm.model,
                sum(  q[arc][t] for (arc, conns) in bus_arcs if t in conns)
                + sum(qsw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
                + sum( qt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
                ==
                sum(qg[gen][t] for (gen, conns) in bus_gens if t in conns)
                - sum(qd[load][t] for (load, conns) in bus_loads if t in conns)
                - sum(qs[strg][t] for (strg, conns) in bus_storage if t in conns)
                + ( sum(Gt[idx,jdx]*vr0[idx]*vi0[jdx]+Bt[idx,jdx]*vr0[idx]*vr0[jdx] for (jdx,u) in ungrounded_terminals)
                -sum(Gt[idx,jdx]*vi0[idx]*vr0[jdx]-Bt[idx,jdx]*vi0[idx]*vi0[jdx] for (jdx,u) in ungrounded_terminals)
                )
                + ( sum(Gt[idx,jdx]*(vr[t]*vi0[jdx]+vr0[idx]*vi[u]-2*vr0[idx]*vi0[jdx])+Bt[idx,jdx]*(vr[t]*vr0[jdx]+vr0[idx]*vr[u]-2*vr0[idx]*vr0[jdx]) for (jdx,u) in ungrounded_terminals)
                -sum(Gt[idx,jdx]*(vi[t]*vr0[jdx]+vi0[idx]*vr[u]-2*vi0[idx]*vr0[jdx])-Bt[idx,jdx]*(vi[t]*vi0[jdx]+vi0[idx]*vi[u]-2*vi0[idx]*vi0[jdx]) for (jdx,u) in ungrounded_terminals)
                )
            )
            push!(cstr_q, cq)
        end
    end

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = cstr_q

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


"""
    constraint_capacitor_on_off(pm::FOTRUPowerModel, i::Int; nw::Int=nw_id_default)

Add constraints to model capacitor switching similar to FBSUBFPowerModel
"""
function constraint_capacitor_on_off(pm::FOTRUPowerModel, nw::Int, i::Int, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    cap_state = var(pm, nw, :capacitor_state, bus_shunts[1][1])
    shunt = ref(pm, nw, :shunt, bus_shunts[1][1])
    ϵ = 1e-5
    M_q = 1e5
    M_v = 2
    elem_type = shunt["controls"]["element"]["type"]
    if shunt["controls"]["type"] == CAP_REACTIVE_POWER
        bus_idx = shunt["controls"]["terminal"] == 1 ? (shunt["controls"]["element"]["index"], shunt["controls"]["element"]["f_bus"], shunt["controls"]["element"]["t_bus"]) : (shunt["controls"]["element"]["index"], shunt["controls"]["element"]["t_bus"], shunt["controls"]["element"]["f_bus"])
        q_fr = elem_type == "branch" ? var(pm, nw, :q)[bus_idx] : elem_type == "switch" ? var(pm, nw, :qsw) : var(pm, nw, :qt, bus_idx)
        JuMP.@constraint(pm.model, sum(q_fr) - shunt["controls"]["onsetting"] ≤ M_q*cap_state[shunt["connections"][1]] - ϵ*(1-cap_state[shunt["connections"][1]]))
        JuMP.@constraint(pm.model, sum(q_fr) - shunt["controls"]["offsetting"] ≥ -M_q*(1-cap_state[shunt["connections"][1]]) - ϵ*cap_state[shunt["connections"][1]])
        JuMP.@constraint(pm.model, cap_state .== cap_state[shunt["connections"][1]])
        if shunt["controls"]["voltoverride"]
            for (idx,val) in enumerate(shunt["connections"])
                vr_cap = var(pm, nw, :vr, i)[val]
                vi_cap = var(pm, nw, :vi, i)[val]
                vr0_cap = var(pm, nw, :vr0, i)[findfirst(isequal(val), ref(pm, nw, :bus, i, "terminals"))]
                vi0_cap = var(pm, nw, :vi0, i)[findfirst(isequal(val), ref(pm, nw, :bus, i, "terminals"))]
                JuMP.@constraint(pm.model, 2*vr_cap*vr0_cap + 2*vi_cap*vi0_cap - vr0_cap^2 - vi0_cap^2 - shunt["controls"]["vmin"]^2 ≥ -M_v*cap_state[val] + ϵ*(1-cap_state[val]))
                JuMP.@constraint(pm.model, 2*vr_cap*vr0_cap + 2*vi_cap*vi0_cap - vr0_cap^2 - vi0_cap^2 - shunt["controls"]["vmax"]^2 ≤ M_v*(1-cap_state[val]) - ϵ*cap_state[val])
            end
        end
    else
        for (idx,val) in enumerate(shunt["connections"])
            if shunt["controls"]["type"][idx] == CAP_VOLTAGE
                bus_idx = shunt["controls"]["terminal"][idx] == 1 ? shunt["controls"]["element"]["f_bus"] : shunt["controls"]["element"]["t_bus"]
                vr_cap = var(pm, nw, :vr, i)[val]
                vi_cap = var(pm, nw, :vi, i)[val]
                vr0_cap = var(pm, nw, :vr0, i)[findfirst(isequal(val), ref(pm, nw, :bus, i, "terminals"))]
                vi0_cap = var(pm, nw, :vi0, i)[findfirst(isequal(val), ref(pm, nw, :bus, i, "terminals"))]
                JuMP.@constraint(pm.model, 2*vr_cap*vr0_cap + 2*vi_cap*vi0_cap - vr0_cap^2 - vi0_cap^2 - shunt["controls"]["onsetting"][idx]^2 ≤ M_v*cap_state[val] - ϵ*(1-cap_state[val]))
                JuMP.@constraint(pm.model, 2*vr_cap*vr0_cap + 2*vi_cap*vi0_cap - vr0_cap^2 - vi0_cap^2 - shunt["controls"]["offsetting"][idx]^2 ≥ -M_v*(1-cap_state[val]) - ϵ*cap_state[val])
            end
            if shunt["controls"]["voltoverride"][idx]
                vr_cap = var(pm, nw, :vr, i)[val]
                vi_cap = var(pm, nw, :vi, i)[val]
                vr0_cap = var(pm, nw, :vr0, i)[findfirst(isequal(val), ref(pm, nw, :bus, i, "terminals"))]
                vi0_cap = var(pm, nw, :vi0, i)[findfirst(isequal(val), ref(pm, nw, :bus, i, "terminals"))]
                JuMP.@constraint(pm.model, 2*vr_cap*vr0_cap + 2*vi_cap*vi0_cap - vr0_cap^2 - vi0_cap^2 - shunt["controls"]["vmin"][idx]^2 ≥ -M_v*cap_state[val] + ϵ*(1-cap_state[val]))
                JuMP.@constraint(pm.model, 2*vr_cap*vr0_cap + 2*vi_cap*vi0_cap - vr0_cap^2 - vi0_cap^2 - shunt["controls"]["vmax"][idx]^2 ≤ M_v*(1-cap_state[val]) - ϵ*cap_state[val])
            end
            if shunt["controls"]["type"][idx] == CAP_DISABLED
                JuMP.@constraint(pm.model, cap_state[val] == 1 )
            end
        end
    end
end


"""
    constraint_mc_ohms_yt_from(pm::FOTRUPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})

Creates Ohms constraints by linearizing (similar to power balance constraints) around initial operating point.
"""
function constraint_mc_ohms_yt_from(pm::FOTRUPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})
    p_fr  = [var(pm, nw, :p, f_idx)[t] for t in f_connections]
    q_fr  = [var(pm, nw, :q, f_idx)[t] for t in f_connections]
    vr_fr = [var(pm, nw, :vr, f_bus)[t] for t in f_connections]
    vr_to = [var(pm, nw, :vr, t_bus)[t] for t in t_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[t] for t in f_connections]
    vi_to = [var(pm, nw, :vi, t_bus)[t] for t in t_connections]
    vr0_fr = [var(pm, nw, :vr0, f_bus)[findfirst(isequal(t), ref(pm, nw, :bus, f_bus, "terminals"))] for t in f_connections]
    vr0_to = [var(pm, nw, :vr0, t_bus)[findfirst(isequal(t), ref(pm, nw, :bus, t_bus, "terminals"))] for t in t_connections]
    vi0_fr = [var(pm, nw, :vi0, f_bus)[findfirst(isequal(t), ref(pm, nw, :bus, f_bus, "terminals"))] for t in f_connections]
    vi0_to = [var(pm, nw, :vi0, t_bus)[findfirst(isequal(t), ref(pm, nw, :bus, t_bus, "terminals"))] for t in t_connections]

    p0_fr =  vr0_fr.*(G*vr0_fr-G*vr0_to-B*vi0_fr+B*vi0_to)
        +vi0_fr.*(G*vi0_fr-G*vi0_to+B*vr0_fr-B*vr0_to)
        +vr0_fr.*(G_fr*vr0_fr-B_fr*vi0_fr)
        +vi0_fr.*(G_fr*vi0_fr+B_fr*vr0_fr)
    q0_fr = -vr0_fr.*(G*vi0_fr-G*vi0_to+B*vr0_fr-B*vr0_to)
        +vi0_fr.*(G*vr0_fr-G*vr0_to-B*vi0_fr+B*vi0_to)
        -vr0_fr.*(G_fr*vi0_fr+B_fr*vr0_fr)
        +vi0_fr.*(G_fr*vr0_fr-B_fr*vi0_fr)

    for (idx,fc) in enumerate(f_connections)
        JuMP.@constraint(pm.model,
                p_fr[idx] ==  p0_fr[idx] + sum(G[idx,jdx]*(vr_fr[jdx]*vr0_fr[idx]+vr0_fr[jdx]*vr_fr[idx]-2*vr0_fr[jdx]*vr0_fr[idx]) - G[idx,jdx]*(vr_to[jdx]*vr0_fr[idx]+vr0_to[jdx]*vr_fr[idx]-2*vr0_to[jdx]*vr0_fr[idx])
                            -B[idx,jdx]*(vi_fr[jdx]*vr0_fr[idx]+vi0_fr[jdx]*vr_fr[idx]-2*vi0_fr[jdx]*vr0_fr[idx]) + B[idx,jdx]*(vi_to[jdx]*vr0_fr[idx]+vi0_to[jdx]*vr_fr[idx]-2*vi0_to[jdx]*vr0_fr[idx]) for (jdx,tc) in enumerate(t_connections))
                        +sum(G[idx,jdx]*(vi_fr[jdx]*vi0_fr[idx]+vi0_fr[jdx]*vi_fr[idx]-2*vi0_fr[jdx]*vi0_fr[idx]) - G[idx,jdx]*(vi_to[jdx]*vi0_fr[idx]+vi0_to[jdx]*vi_fr[idx]-2*vi0_to[jdx]*vi0_fr[idx])
                            +B[idx,jdx]*(vr_fr[jdx]*vi0_fr[idx]+vr0_fr[jdx]*vi_fr[idx]-2*vr0_fr[jdx]*vi0_fr[idx]) - B[idx,jdx]*(vr_to[jdx]*vi0_fr[idx]+vr0_to[jdx]*vi_fr[idx]-2*vr0_to[jdx]*vi0_fr[idx]) for (jdx,tc) in enumerate(t_connections))
                        +sum(G_fr[idx,jdx]*(vr_fr[jdx]*vr0_fr[idx]+vr0_fr[jdx]*vr_fr[idx]-2*vr0_fr[jdx]*vr0_fr[idx]) - B_fr[idx,jdx]*(vi_fr[jdx]*vr0_fr[idx]+vi0_fr[jdx]*vr_fr[idx]-2*vi0_fr[jdx]*vr0_fr[idx]) for (jdx,tc) in enumerate(t_connections))
                        +sum(G_fr[idx,jdx]*(vi_fr[jdx]*vi0_fr[idx]+vi0_fr[jdx]*vi_fr[idx]-2*vi0_fr[jdx]*vi0_fr[idx]) + B_fr[idx,jdx]*(vr_fr[jdx]*vi0_fr[idx]+vr0_fr[jdx]*vi_fr[idx]-2*vr0_fr[jdx]*vi0_fr[idx]) for (jdx,tc) in enumerate(t_connections))
        )

        JuMP.@constraint(pm.model,
                q_fr[idx] == q0_fr[idx] - sum(G[idx,jdx]*(vi_fr[jdx]*vr0_fr[idx]+vi0_fr[jdx]*vr_fr[idx]-2*vi0_fr[jdx]*vr0_fr[idx]) - G[idx,jdx]*(vi_to[jdx]*vr0_fr[idx]+vi0_to[jdx]*vr_fr[idx]-2*vi0_to[jdx]*vr0_fr[idx])
                            +B[idx,jdx]*(vr_fr[jdx]*vr0_fr[idx]+vr0_fr[jdx]*vr_fr[idx]-2*vr0_fr[jdx]*vr0_fr[idx]) - B[idx,jdx]*(vr_to[jdx]*vr0_fr[idx]+vr0_to[jdx]*vr_fr[idx]-2*vr0_to[jdx]*vr0_fr[idx]) for (jdx,tc) in enumerate(t_connections))
                        +sum(G[idx,jdx]*(vr_fr[jdx]*vi0_fr[idx]+vr0_fr[jdx]*vi_fr[idx]-2*vr0_fr[jdx]*vi0_fr[idx]) - G[idx,jdx]*(vr_to[jdx]*vi0_fr[idx]+vr0_to[jdx]*vi_fr[idx]-2*vr0_to[jdx]*vi0_fr[idx])
                            -B[idx,jdx]*(vi_fr[jdx]*vi0_fr[idx]+vi0_fr[jdx]*vi_fr[idx]-2*vi0_fr[jdx]*vi0_fr[idx]) + B[idx,jdx]*(vi_to[jdx]*vi0_fr[idx]+vi0_to[jdx]*vi_fr[idx]-2*vi0_to[jdx]*vi0_fr[idx]) for (jdx,tc) in enumerate(t_connections))
                        -sum(G_fr[idx,jdx]*(vi_fr[jdx]*vr0_fr[idx]+vi0_fr[jdx]*vr_fr[idx]-2*vi0_fr[jdx]*vr0_fr[idx]) + B_fr[idx,jdx]*(vr_fr[jdx]*vr0_fr[idx]+vr0_fr[jdx]*vr_fr[idx]-2*vr0_fr[jdx]*vr0_fr[idx]) for (jdx,tc) in enumerate(t_connections))
                        +sum(G_fr[idx,jdx]*(vr_fr[jdx]*vi0_fr[idx]+vr0_fr[jdx]*vi_fr[idx]-2*vr0_fr[jdx]*vi0_fr[idx]) - B_fr[idx,jdx]*(vi_fr[jdx]*vi0_fr[idx]+vi0_fr[jdx]*vi_fr[idx]-2*vi0_fr[jdx]*vi0_fr[idx]) for (jdx,tc) in enumerate(t_connections))
        )
    end
end


"""
    constraint_mc_ohms_yt_to(pm::FOTRUPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix, B::Matrix, G_to::Matrix, B_to::Matrix)

Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_mc_ohms_yt_to(pm::FOTRUPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix, B::Matrix, G_to::Matrix, B_to::Matrix)
    constraint_mc_ohms_yt_from(pm, nw, t_bus, f_bus, t_idx, f_idx, t_connections, f_connections, G, B, G_to, B_to)
end


"""
    constraint_mc_load_power(pm::FOTRUPowerModel, load_id::Int; nw::Int=nw_id_default, report::Bool=true)

Load model is linearized around initial operating point similar to FBSUBFPowerModel.
"""
function constraint_mc_load_power(pm::FOTRUPowerModel, load_id::Int; nw::Int=nw_id_default, report::Bool=true)
    # shared variables and parameters
    load = ref(pm, nw, :load, load_id)
    connections = load["connections"]
    pd0 = load["pd"]
    qd0 = load["qd"]
    bus_id = load["load_bus"]
    bus = ref(pm, nw, :bus, bus_id)
    terminals = bus["terminals"]

    # calculate load params
    a, alpha, b, beta = _load_expmodel_params(load, bus)

    # first-order approximation
    if load["configuration"]==WYE
        if load["model"]==POWER
            pd_bus = a
            qd_bus = b
        elseif load["model"]==IMPEDANCE
            vr = var(pm, nw, :vr, bus_id)
            vi = var(pm, nw, :vi, bus_id)
            vr0 = [var(pm, nw, :vr0, bus_id)[findfirst(isequal(c), terminals)] for c in connections]
            vi0 = [var(pm, nw, :vi0, bus_id)[findfirst(isequal(c), terminals)] for c in connections]

            pd_bus = []
            qd_bus = []

            for (idx,c) in enumerate(connections)
                vm0 = vr0[idx]^2 + vi0[idx]^2
                push!(pd_bus, a[idx]*(-vm0 + 2*vr[c]*vr0[idx] + 2*vi[c]*vi0[idx]))
                push!(qd_bus, b[idx]*(-vm0 + 2*vr[c]*vr0[idx] + 2*vi[c]*vi0[idx]))
            end
        else
            vr = var(pm, nw, :vr, bus_id)
            vi = var(pm, nw, :vi, bus_id)
            vr0 = [var(pm, nw, :vr0, bus_id)[findfirst(isequal(c), terminals)] for c in connections]
            vi0 = [var(pm, nw, :vi0, bus_id)[findfirst(isequal(c), terminals)] for c in connections]

            pd_bus = []
            qd_bus = []

            for (idx,c) in enumerate(connections)
                vm0 = sqrt(vr0[idx]^2 + vi0[idx]^2)
                push!(pd_bus, a[idx]*(vm0 + (vr[c]*vr0[idx] + vi[c]*vi0[idx]-vm0^2)/vm0))
                push!(qd_bus, b[idx]*(vm0 + (vr[c]*vr0[idx] + vi[c]*vi0[idx]-vm0^2)/vm0))
            end
        end

        if report
            sol(pm, nw, :load, load_id)[:pd_bus] = pd_bus
            sol(pm, nw, :load, load_id)[:qd_bus] = qd_bus
            sol(pm, nw, :load, load_id)[:pd] = pd_bus
            sol(pm, nw, :load, load_id)[:qd] = qd_bus
        end
        pd_bus = JuMP.Containers.DenseAxisArray(pd_bus, connections)
        qd_bus = JuMP.Containers.DenseAxisArray(qd_bus, connections)

        var(pm, nw, :pd_bus)[load_id] = pd_bus
        var(pm, nw, :qd_bus)[load_id] = qd_bus

    # zero-order approximation
    elseif load["configuration"]==DELTA
        vr0 = var(pm, nw, :vr0, bus_id)
        vi0 = var(pm, nw, :vi0, bus_id)

        nph = length(a)

        prev = Dict(c=>connections[(idx+nph-2)%nph+1] for (idx,c) in enumerate(connections))
        next = Dict(c=>connections[idx%nph+1] for (idx,c) in enumerate(connections))

        vrd0 = [vr0[idx]-vr0[next[idx]] for (idx, c) in enumerate(connections)]
        vid0 = [vi0[idx]-vi0[next[idx]] for (idx, c) in enumerate(connections)]

        crd0 = Array{Any,1}(undef, nph)
        cid0 = Array{Any,1}(undef, nph)
        for (idx, c) in enumerate(connections)
            crd0[c] = a[idx]*vrd0[c]*(vrd0[c]^2+vid0[c]^2)^(alpha[idx]/2-1)+b[idx]*vid0[c]*(vrd0[c]^2+vid0[c]^2)^(beta[idx]/2 -1)
            cid0[c] = a[idx]*vid0[c]*(vrd0[c]^2+vid0[c]^2)^(alpha[idx]/2-1)-b[idx]*vrd0[c]*(vrd0[c]^2+vid0[c]^2)^(beta[idx]/2 -1)
        end
        
        crd0_bus = [crd0[idx]-crd0[prev[idx]] for (idx, c) in enumerate(connections)]
        cid0_bus = [cid0[idx]-cid0[prev[idx]] for (idx, c) in enumerate(connections)]

        pd_bus = [ vr0[c]*crd0_bus[c]+vi0[c]*cid0_bus[c] for (idx,c) in enumerate(connections)]
        qd_bus = [-vr0[c]*cid0_bus[c]+vi0[c]*crd0_bus[c] for (idx,c) in enumerate(connections)]

        pd_bus = JuMP.Containers.DenseAxisArray(pd_bus, connections)
        qd_bus = JuMP.Containers.DenseAxisArray(qd_bus, connections)

        var(pm, nw, :pd_bus)[load_id] = pd_bus
        var(pm, nw, :qd_bus)[load_id] = qd_bus

        if report
            sol(pm, nw, :load, load_id)[:pd_bus] = pd_bus
            sol(pm, nw, :load, load_id)[:qd_bus] = qd_bus

            pd = []
            qd = []
            for (idx,c) in enumerate(connections)
                push!(pd, JuMP.@expression(pm.model, a[idx]*(vrd0[c]^2+vid0[c]^2)^(alpha[idx]/2) ))
                push!(qd, JuMP.@expression(pm.model, b[idx]*(vrd0[c]^2+vid0[c]^2)^(beta[idx]/2)  ))
            end
            sol(pm, nw, :load, load_id)[:pd] = pd
            sol(pm, nw, :load, load_id)[:qd] = qd
        end
    end
end


"""
    constraint_mc_transformer_power_yy(pm::FOTRUPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)

Add all constraints required to model a two-winding, wye-wye connected transformer similar to ACRUPowerModel.
"""
function constraint_mc_transformer_power_yy(pm::FOTRUPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    transformer = ref(pm, nw, :transformer, trans_id)

    vr_fr = var(pm, nw, :vr, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    vr0_fr = [var(pm, nw, :vr0, f_bus)[findfirst(isequal(c), ref(pm, nw, :bus, f_bus, "terminals"))] for c in f_connections]
    vr0_to = [var(pm, nw, :vr0, t_bus)[findfirst(isequal(c), ref(pm, nw, :bus, t_bus, "terminals"))] for c in t_connections]
    vi0_fr = [var(pm, nw, :vi0, f_bus)[findfirst(isequal(c), ref(pm, nw, :bus, f_bus, "terminals"))] for c in f_connections]
    vi0_to = [var(pm, nw, :vi0, t_bus)[findfirst(isequal(c), ref(pm, nw, :bus, t_bus, "terminals"))] for c in t_connections]

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

                # (cr+jci) = (p-jq)/(vr0-j⋅vi0)
                cr = JuMP.@expression(pm.model, ( p_to[idx]*vr0_to[idx] + q_to[idx]*vi0_to[idx])/(vr0_to[idx]^2+vi0_to[idx]^2))
                ci = JuMP.@expression(pm.model, (-q_to[idx]*vr0_to[idx] + p_to[idx]*vi0_to[idx])/(vr0_to[idx]^2+vi0_to[idx]^2))
                # linearized v_drop = (cr+jci)⋅(r+jx)
                vr_drop = JuMP.@expression(pm.model, r*cr-x*ci)
                vi_drop = JuMP.@expression(pm.model, r*ci+x*cr)

                # linearized voltage magnitude squared v_lin_sq = 2⋅vr⋅vr0 + 2⋅vi⋅vi0 - (vr0^2+vi0^2)
                # outer approximation of upper limits: -(v_ref+δ) ≤ (vr_fr-vr_drop) ≤ (v_ref+δ)
                #                                      -(v_ref+δ) ≤ (vi_fr-vi_drop) ≤ (v_ref+δ)
                #                              -\sqrt(2)(v_ref+δ) ≤ (vr_fr-vr_drop) + (vi_fr-vi_drop) ≤ \sqrt(2)(v_ref+δ)
                #                              -\sqrt(2)(v_ref+δ) ≤ (vr_fr-vr_drop) - (vi_fr-vi_drop) ≤ \sqrt(2)(v_ref+δ)
                JuMP.@constraint(pm.model, (vr_fr[fc]-vr_drop) ≤  (v_ref + δ))
                JuMP.@constraint(pm.model, (vr_fr[fc]-vr_drop) ≥ -(v_ref + δ))
                JuMP.@constraint(pm.model, (vi_fr[fc]-vi_drop) ≤  (v_ref + δ))
                JuMP.@constraint(pm.model, (vi_fr[fc]-vi_drop) ≥ -(v_ref + δ))
                JuMP.@constraint(pm.model, (vr_fr[fc]-vr_drop) + (vi_fr[fc]-vi_drop) ≤  sqrt(2)*(v_ref + δ)^2)
                JuMP.@constraint(pm.model, (vr_fr[fc]-vr_drop) + (vi_fr[fc]-vi_drop) ≥ -sqrt(2)*(v_ref + δ)^2)
                JuMP.@constraint(pm.model, (vr_fr[fc]-vr_drop) - (vi_fr[fc]-vi_drop) ≥ -sqrt(2)*(v_ref + δ)^2)
                JuMP.@constraint(pm.model, (vr_fr[fc]-vr_drop) - (vi_fr[fc]-vi_drop) ≤  sqrt(2)*(v_ref + δ)^2)
                JuMP.@constraint(pm.model, (vr_fr[fc]-vr_drop)^2 + (vi_fr[fc]-vi_drop)^2 ≥ (v_ref - δ)^2)
                # TODO: linearized lower limits: (v_ref-δ)^2 ≤ v_lin_sq
                # JuMP.@constraint(pm.model, 2*vr0_fr[fc]*(vr_fr[fc]-vr_drop) + 2*vi0_fr[fc]*(vi_fr[fc]-vi_drop) - vr_fr[fc]^2 - vi_fr[fc]^2 ≥ (v_ref - δ)^2)
                JuMP.@constraint(pm.model, (2*vr_fr[fc]*vr0_fr[idx] + 2*vi_fr[fc]*vi0_fr[idx] - vr_fr[fc]^2 - vi_fr[fc]^2)/1.1^2 ≤ 2*vr_to[tc]*vr0_to[idx] + 2*vi_to[tc]*vi0_to[idx] - vr_to[tc]^2 - vi_to[tc]^2)
                JuMP.@constraint(pm.model, (2*vr_fr[fc]*vr0_fr[idx] + 2*vi_fr[fc]*vi0_fr[idx] - vr_fr[fc]^2 - vi_fr[fc]^2)/0.9^2 ≥ 2*vr_to[tc]*vr0_to[idx] + 2*vi_to[tc]*vi0_to[idx] - vr_to[tc]^2 - vi_to[tc]^2)
            end
        end
    end

    JuMP.@constraint(pm.model, p_fr + p_to .== 0)
    JuMP.@constraint(pm.model, q_fr + q_to .== 0)
end


"""
    constraint_mc_transformer_power_dy(pm::FOTRUPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)

Add all constraints required to model a two-winding, delta-wye connected transformer similar to ACRUPowerModel
with power constraints using initial operating point voltage instead of actual voltage variables.
"""
function constraint_mc_transformer_power_dy(pm::FOTRUPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
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
    vr0_p_fr = [var(pm, nw, :vr0, f_bus)[findfirst(isequal(c), ref(pm, nw, :bus, f_bus, "terminals"))] for c in f_connections]
    vr0_p_to = [var(pm, nw, :vr0, t_bus)[findfirst(isequal(c), ref(pm, nw, :bus, t_bus, "terminals"))] for c in t_connections]
    vi0_p_fr = [var(pm, nw, :vi0, f_bus)[findfirst(isequal(c), ref(pm, nw, :bus, f_bus, "terminals"))] for c in f_connections]
    vi0_p_to = [var(pm, nw, :vi0, t_bus)[findfirst(isequal(c), ref(pm, nw, :bus, t_bus, "terminals"))] for c in t_connections]

    id_re = Array{Any,1}(undef, nph)
    id_im = Array{Any,1}(undef, nph)

    for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))
        id_re[idx] = JuMP.@expression(pm.model, (p_to[tc]*vr0_p_to[idx]+q_to[tc]*vi0_p_to[idx])/(tm_scale*tm[idx]*pol)/(vr0_p_to[idx]^2+vi0_p_to[idx]^2))
        id_im[idx] = JuMP.@expression(pm.model, (p_to[tc]*vi0_p_to[idx]-q_to[tc]*vr0_p_to[idx])/(tm_scale*tm[idx]*pol)/(vr0_p_to[idx]^2+vi0_p_to[idx]^2))
    end
    for (idx,(fc,tc)) in enumerate(zip(f_connections, t_connections))
        jdx = (idx-1+nph-1)%nph+1

        JuMP.@constraint(pm.model, p_fr[fc] ==
             vr0_p_fr[idx]*( id_re[jdx]-id_re[idx])
            -vi0_p_fr[idx]*(-id_im[jdx]+id_im[idx])
        )
        JuMP.@constraint(pm.model, q_fr[fc] ==
             vr0_p_fr[idx]*(-id_im[jdx]+id_im[idx])
            +vi0_p_fr[idx]*( id_re[jdx]-id_re[idx])
        )
    end
end
