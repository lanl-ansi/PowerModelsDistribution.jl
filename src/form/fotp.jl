# The model in polar coordinates is linearized around an initial operating point using a
# first order Taylor approximation (FOT) method


"""
    constraint_mc_load_current_delta(pm::FOTPUPowerModel, nw::Int, load_id::Int, load_bus_id::Int, cp::Vector, cq::Vector)

No loads require a current variable. Delta loads are zero-order approximations and
wye loads are first-order approximations around the initial operating point.
"""
function constraint_mc_load_current_delta(pm::FOTPUPowerModel, nw::Int, load_id::Int, load_bus_id::Int, cp::Vector, cq::Vector)
end


"""
    variable_mc_bus_voltage(pm::FOTPUPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

Voltage variables are defined in polar coordinates similar to ACPUPowerModel.
An initial operating point is specified for linearization.
"""
function variable_mc_bus_voltage(pm::FOTPUPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_bus_voltage_angle(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_bus_voltage_magnitude_only(pm; nw=nw, bounded=bounded, report=report)

    # initial operating point for linearization (using flat-start)
    var(pm, nw)[:vm0] = Dict{Int,Vector{Float64}}()  #(i => [1 1 1] for i in ids(pm, nw, :bus))
    var(pm, nw)[:va0] = Dict{Int,Vector{Float64}}()  #(i => [0 -2*pi/3 2*pi/3] for i in ids(pm, nw, :bus))

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
        va = haskey(busref, "va_start") ? busref["va_start"] : haskey(busref, "va") ? busref["va"] : [deg2rad.([0, -120, 120])..., zeros(length(terminals))...][terminals]

        JuMP.set_start_value.(var(pm, nw, :vm, id), vm)
        JuMP.set_start_value.(var(pm, nw, :va, id), va)

        # TODO: update initial operating point with warm-start (causes infeasbility if not flat start)
        var(pm, nw, :vm0)[id] = fill(1.0, ncnd)  # vm
        var(pm, nw, :va0)[id] = get(busref, "va_start", default_va)  # va
    end
end


@doc raw"""
    constraint_mc_power_balance(pm::FOTPUPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})

Power balance equations similar to ACPUPowerModel.
The nonlinear functions are approximated around initial operating point.

```math
\begin{align}
&\text{Initial operating points: }  v_{m0}^{t} \angle v_{a0}^t,~v_{m0}^u \angle v_{a0}^u\\
& {v_m^t}^2 \Rightarrow {v_{m0}^t}^2+2 \cdot v_{m0}^t \cdot (v_m^t-v_{m0}^t)\\
& v_m^t \cdot v_m^u \cdot \cos(v_a^t-v_a^u) \Rightarrow v_{m0}^t \cdot v_{m0}^u \cdot \cos(v_{a0}^t-v_{a0}^u) +
\begin{bmatrix}
v_{m0}^u \cdot \cos(v_{a0}^t-v_{a0}^u) \\
v_{m0}^t \cdot \cos(v_{a0}^t-v_{a0}^u) \\
-v_{m0}^t \cdot v_{m0}^u \cdot \sin(v_{a0}^t-v_{a0}^u) \\
v_{m0}^t \cdot v_{m0}^u \cdot \sin(v_{a0}^t-v_{a0}^u)
\end{bmatrix}^\top
\begin{bmatrix}
v_m^t-v_{m0}^t \\
v_m^u-v_{m0}^u \\
v_a^t-v_{a0}^t \\
v_a^u-v_{a0}^u
\end{bmatrix} \\
& v_m^t \cdot v_m^u \cdot \sin(v_a^t-v_a^u) \Rightarrow v_{m0}^t \cdot v_{m0}^u \cdot \sin(v_{a0}^t-v_{a0}^u) +
\begin{bmatrix}
v_{m0}^u \cdot \sin(v_{a0}^t-v_{a0}^u) \\
v_{m0}^t \cdot \sin(v_{a0}^t-v_{a0}^u) \\
v_{m0}^t \cdot v_{m0}^u \cdot \cos(v_{a0}^t-v_{a0}^u) \\
-v_{m0}^t \cdot v_{m0}^u \cdot \cos(v_{a0}^t-v_{a0}^u)
\end{bmatrix}^\top
\begin{bmatrix}
v_m^t-v_{m0}^t \\
v_m^u-v_{m0}^u \\
v_a^t-v_{a0}^t \\
v_a^u-v_{a0}^u
\end{bmatrix}
\end{align}
```
"""
function constraint_mc_power_balance(pm::FOTPUPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vm   = var(pm, nw, :vm, i)
    va   = var(pm, nw, :va, i)
    vm0   = var(pm, nw, :vm0, i)
    va0   = var(pm, nw, :va0, i)
    p    = get(var(pm, nw),      :p, Dict()); _check_var_keys(  p, bus_arcs, "active power", "branch")
    q    = get(var(pm, nw),      :q, Dict()); _check_var_keys(  q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, nw), :pg_bus, Dict()); _check_var_keys( pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, nw), :qg_bus, Dict()); _check_var_keys( qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, nw),     :ps, Dict()); _check_var_keys( ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, nw),     :qs, Dict()); _check_var_keys( qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, nw),    :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, nw),    :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(var(pm, nw),     :pt, Dict()); _check_var_keys( pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(var(pm, nw),     :qt, Dict()); _check_var_keys( qt, bus_arcs_trans, "reactive power", "transformer")
    pd   = get(var(pm, nw), :pd_bus, Dict()); _check_var_keys( pd, bus_loads, "active power", "load")
    qd   = get(var(pm, nw), :qd_bus, Dict()); _check_var_keys( pd, bus_loads, "reactive power", "load")

    Gs, Bs = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []
    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        if any(Bs[idx,jdx] != 0 for (jdx, u) in ungrounded_terminals if idx != jdx) || any(Gs[idx,jdx] != 0 for (jdx, u) in ungrounded_terminals if idx != jdx)
            cp = JuMP.@constraint(pm.model,
                  sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( pd[l][t] for (l, conns) in bus_loads if t in conns)
                + ( Gs[idx,idx]*(vm0[idx]^2+2*vm0[idx]*(vm[t]-vm0[idx]))
                    +sum( Gs[idx,jdx] * vm0[idx]*vm0[jdx] * cos(va0[idx]-va0[jdx])
                         +Bs[idx,jdx] * vm0[idx]*vm0[jdx] * sin(va0[idx]-va0[jdx])
                         +[Gs[idx,jdx]*vm0[jdx]*cos(va0[idx]-va0[jdx]) Gs[idx,jdx]*vm0[idx]*cos(va0[idx]-va0[jdx]) -Gs[idx,jdx]*vm0[idx]*vm0[jdx]*sin(va0[idx]-va0[jdx])  Gs[idx,jdx]*vm0[idx]*vm0[jdx]*sin(va0[idx]-va0[jdx])]*[vm[t]-vm0[idx];vm[u]-vm0[jdx];va[t]-va0[idx];va[u]-va0[jdx]]
                         +[Bs[idx,jdx]*vm0[jdx]*sin(va0[idx]-va0[jdx]) Bs[idx,jdx]*vm0[idx]*sin(va0[idx]-va0[jdx])  Bs[idx,jdx]*vm0[idx]*vm0[jdx]*cos(va0[idx]-va0[jdx]) -Bs[idx,jdx]*vm0[idx]*vm0[jdx]*cos(va0[idx]-va0[jdx])]*[vm[t]-vm0[idx];vm[u]-vm0[jdx];va[t]-va0[idx];va[u]-va0[jdx]]
                        for (jdx,u) in ungrounded_terminals if idx != jdx)
                )
                ==
                0.0
            )
            push!(cstr_p, cp)

            cq = JuMP.@constraint(pm.model,
                  sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( qd[l][t] for (l, conns) in bus_loads if t in conns)
                + ( -Bs[idx,idx]*(vm0[idx]^2+2*vm0[idx]*(vm[t]-vm0[idx]))
                    -sum( Bs[idx,jdx] * vm0[idx]*vm0[jdx] * cos(va0[idx]-va0[jdx])
                         -Gs[idx,jdx] * vm0[idx]*vm0[jdx] * sin(va0[idx]-va0[jdx])
                         +[Bs[idx,jdx]*vm0[jdx]*cos(va0[idx]-va0[jdx])   Bs[idx,jdx]*vm0[idx]*cos(va0[idx]-va0[jdx]) -Bs[idx,jdx]*vm0[idx]*vm0[jdx]*sin(va0[idx]-va0[jdx]) Bs[idx,jdx]*vm0[idx]*vm0[jdx]*sin(va0[idx]-va0[jdx])]*[vm[t]-vm0[idx];vm[u]-vm0[jdx];va[t]-va0[idx];va[u]-va0[jdx]]
                         +[-Gs[idx,jdx]*vm0[jdx]*sin(va0[idx]-va0[jdx]) -Gs[idx,jdx]*vm0[idx]*sin(va0[idx]-va0[jdx]) -Gs[idx,jdx]*vm0[idx]*vm0[jdx]*cos(va0[idx]-va0[jdx]) Gs[idx,jdx]*vm0[idx]*vm0[jdx]*cos(va0[idx]-va0[jdx])]*[vm[t]-vm0[idx];vm[u]-vm0[jdx];va[t]-va0[idx];va[u]-va0[jdx]]
                         for (jdx,u) in ungrounded_terminals if idx != jdx)
                )
                ==
                0.0
            )
            push!(cstr_q, cq)
        else
            cp = JuMP.@constraint(pm.model,
                  sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( pd[l][t] for (l, conns) in bus_loads if t in conns)
                + Gs[idx,idx]*(vm0[idx]^2+2*vm0[idx]*(vm[t]-vm0[idx]))
                ==
                0.0
            )
            push!(cstr_p, cp)

            cq = JuMP.@constraint(pm.model,
                  sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( qd[l][t] for (l, conns) in bus_loads if t in conns)
                - Bs[idx,idx]*(vm0[idx]^2+2*vm0[idx]*(vm[t]-vm0[idx]))
                ==
                0.0
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


@doc raw"""
    constraint_mc_power_balance_capc(pm::FOTPUPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})

Power balance constraints with capacitor control with shunt current calculated using initial operating point.

```math
\begin{align}
    & B_s = b_s ⋅ z,~~ cq_{sh} = B_s ⋅ v, \\
    & B_s \cdot v_m^t \cdot v_m^u \cdot \cos(v_a^t-v_a^u) \Rightarrow B_{s0} \cdot v_{m0}^t \cdot v_{m0}^u \cdot \cos(v_{a0}^t-v_{a0}^u) +
\begin{bmatrix}
B_{s0} \cdot v_{m0}^u \cdot \cos(v_{a0}^t-v_{a0}^u) \\
B_{s0} \cdot v_{m0}^t \cdot \cos(v_{a0}^t-v_{a0}^u) \\
-B_{s0} \cdot v_{m0}^t \cdot v_{m0}^u \cdot \sin(v_{a0}^t-v_{a0}^u) \\
B_{s0} \cdot v_{m0}^t \cdot v_{m0}^u \cdot \sin(v_{a0}^t-v_{a0}^u) \\
v_{m0}^t \cdot v_{m0}^u \cdot \cos(v_{a0}^t-v_{a0}^u)
\end{bmatrix}^\top
\begin{bmatrix}
v_m^t-v_{m0}^t \\
v_m^u-v_{m0}^u \\
v_a^t-v_{a0}^t \\
v_a^u-v_{a0}^u \\
B_{s} -B_{s0}
\end{bmatrix} \\
& B_s \cdot v_m^t \cdot v_m^u \cdot \sin(v_a^t-v_a^u) \Rightarrow B_{s0} \cdot v_{m0}^t \cdot v_{m0}^u \cdot \sin(v_{a0}^t-v_{a0}^u) +
\begin{bmatrix}
 B_{s0} \cdot v_{m0}^u \cdot \sin(v_{a0}^t-v_{a0}^u) \\
 B_{s0} \cdot v_{m0}^t \cdot \sin(v_{a0}^t-v_{a0}^u) \\
 B_{s0} \cdot v_{m0}^t \cdot v_{m0}^u \cdot \cos(v_{a0}^t-v_{a0}^u) \\
 -B_{s0} \cdot v_{m0}^t \cdot v_{m0}^u \cdot \cos(v_{a0}^t-v_{a0}^u) \\
 v_{m0}^t \cdot v_{m0}^u \cdot \sin(v_{a0}^t-v_{a0}^u)
\end{bmatrix}^\top
\begin{bmatrix}
v_m^t-v_{m0}^t \\
v_m^u-v_{m0}^u \\
v_a^t-v_{a0}^t \\
v_a^u-v_{a0}^u \\
B_{s} -B_{s0}
\end{bmatrix}

\end{align}
```
"""
function constraint_mc_power_balance_capc(pm::FOTPUPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vm   = var(pm, nw, :vm, i)
    va   = var(pm, nw, :va, i)
    vm0   = var(pm, nw, :vm0, i)
    va0   = var(pm, nw, :va0, i)
    p    = get(var(pm, nw),      :p, Dict()); _check_var_keys(  p, bus_arcs, "active power", "branch")
    q    = get(var(pm, nw),      :q, Dict()); _check_var_keys(  q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, nw), :pg_bus, Dict()); _check_var_keys( pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, nw), :qg_bus, Dict()); _check_var_keys( qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, nw),     :ps, Dict()); _check_var_keys( ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, nw),     :qs, Dict()); _check_var_keys( qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, nw),    :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, nw),    :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(var(pm, nw),     :pt, Dict()); _check_var_keys( pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(var(pm, nw),     :qt, Dict()); _check_var_keys( qt, bus_arcs_trans, "reactive power", "transformer")
    pd   = get(var(pm, nw), :pd_bus, Dict()); _check_var_keys( pd, bus_loads, "active power", "load")
    qd   = get(var(pm, nw), :qd_bus, Dict()); _check_var_keys( pd, bus_loads, "reactive power", "load")

    # calculate Gs, Bs
    ncnds = length(terminals)
    Gs = fill(0.0, ncnds, ncnds)
    Bs0 = fill(0.0, ncnds, ncnds)
    Bs = convert(Matrix{JuMP.AffExpr}, JuMP.@expression(pm.model, [idx=1:ncnds, jdx=1:ncnds], 0.0))
    for (val, connections) in bus_shunts
        shunt = ref(pm,nw,:shunt,val)
        for (idx,c) in enumerate(connections)
            cap_state = haskey(shunt,"controls") ? var(pm, nw, :capacitor_state, val)[c] : 1.0
            for (jdx,d) in enumerate(connections)
                Gs[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] += shunt["gs"][idx,jdx]
                Bs0[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] += shunt["bs"][idx,jdx]
                Bs[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] = JuMP.@expression(pm.model, Bs[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] + shunt["bs"][idx,jdx]*cap_state)
            end
        end
    end

    cstr_p = []
    cstr_q = []
    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    # add constraints to model capacitor switching
    if !isempty(bus_shunts) && haskey(ref(pm, nw, :shunt, bus_shunts[1][1]), "controls")
        constraint_capacitor_on_off(pm, nw, i, bus_shunts)

        for (idx,t) in ungrounded_terminals
            if any(Bs[idx,jdx] != 0 for (jdx, u) in ungrounded_terminals if idx != jdx) || any(Gs[idx,jdx] != 0 for (jdx, u) in ungrounded_terminals if idx != jdx)
                cp = JuMP.@constraint(pm.model,
                    sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
                    + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                    + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                    - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
                    + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
                    + sum( pd[l][t] for (l, conns) in bus_loads if t in conns)
                    + ( Gs[idx,idx]*(vm0[idx]^2+2*vm0[idx]*(vm[t]-vm0[idx]))
                        +sum( Gs[idx,jdx] * vm0[idx]*vm0[jdx] * cos(va0[idx]-va0[jdx])
                            +Bs0[idx,jdx] * vm0[idx]*vm0[jdx] * sin(va0[idx]-va0[jdx])
                            +[Gs[idx,jdx]*vm0[jdx]*cos(va0[idx]-va0[jdx]) Gs[idx,jdx]*vm0[idx]*cos(va0[idx]-va0[jdx]) -Gs[idx,jdx]*vm0[idx]*vm0[jdx]*sin(va0[idx]-va0[jdx])  Gs[idx,jdx]*vm0[idx]*vm0[jdx]*sin(va0[idx]-va0[jdx])]*[vm[t]-vm0[idx];vm[u]-vm0[jdx];va[t]-va0[idx];va[u]-va0[jdx]]
                            +[Bs0[idx,jdx]*vm0[jdx]*sin(va0[idx]-va0[jdx]) Bs0[idx,jdx]*vm0[idx]*sin(va0[idx]-va0[jdx])  Bs0[idx,jdx]*vm0[idx]*vm0[jdx]*cos(va0[idx]-va0[jdx]) -Bs0[idx,jdx]*vm0[idx]*vm0[jdx]*cos(va0[idx]-va0[jdx]) vm0[t]*vm0[u]*sin(va0[t]-va0[u])]*[vm[t]-vm0[idx];vm[u]-vm0[jdx];va[t]-va0[idx];va[u]-va0[jdx];Bs[idx,jdx]-Bs0[idx,jdx]]
                            for (jdx,u) in ungrounded_terminals if idx != jdx)
                    )
                    ==
                    0.0
                )
                push!(cstr_p, cp)

                cq = JuMP.@constraint(pm.model,
                    sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
                    + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                    + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                    - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
                    + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
                    + sum( qd[l][t] for (l, conns) in bus_loads if t in conns)
                    + ( -Bs0[idx,idx]*vm0[idx]^2 - Bs0[idx,idx]*2*vm0[idx]*(vm[t]-vm0[idx]) - vm0[idx]^2*(Bs[idx,idx]-Bs0[idx,idx])
                        -sum( Bs0[idx,jdx] * vm0[idx]*vm0[jdx] * cos(va0[idx]-va0[jdx])
                            -Gs[idx,jdx] * vm0[idx]*vm0[jdx] * sin(va0[idx]-va0[jdx])
                            +[Bs0[idx,jdx]*vm0[jdx]*cos(va0[idx]-va0[jdx])   Bs0[idx,jdx]*vm0[idx]*cos(va0[idx]-va0[jdx]) -Bs0[idx,jdx]*vm0[idx]*vm0[jdx]*sin(va0[idx]-va0[jdx]) Bs0[idx,jdx]*vm0[idx]*vm0[jdx]*sin(va0[idx]-va0[jdx]) vm0[t]*vm0[u]*cos(va0[t]-va0[u])]*[vm[t]-vm0[idx];vm[u]-vm0[jdx];va[t]-va0[idx];va[u]-va0[jdx];Bs[idx,jdx]-Bs0[idx,jdx]]
                            +[-Gs[idx,jdx]*vm0[jdx]*sin(va0[idx]-va0[jdx]) -Gs[idx,jdx]*vm0[idx]*sin(va0[idx]-va0[jdx]) -Gs[idx,jdx]*vm0[idx]*vm0[jdx]*cos(va0[idx]-va0[jdx]) Gs[idx,jdx]*vm0[idx]*vm0[jdx]*cos(va0[idx]-va0[jdx])]*[vm[t]-vm0[idx];vm[u]-vm0[jdx];va[t]-va0[idx];va[u]-va0[jdx]]
                            for (jdx,u) in ungrounded_terminals if idx != jdx)
                    )
                    ==
                    0.0
                )
                push!(cstr_q, cq)
            else
                cp = JuMP.@constraint(pm.model,
                    sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
                    + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                    + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                    - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
                    + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
                    + sum( pd[l][t] for (l, conns) in bus_loads if t in conns)
                    + Gs[idx,idx]*(vm0[idx]^2+2*vm0[idx]*(vm[t]-vm0[idx]))
                    ==
                    0.0
                )
                push!(cstr_p, cp)

                cq = JuMP.@constraint(pm.model,
                    sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
                    + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                    + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                    - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
                    + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
                    + sum( qd[l][t] for (l, conns) in bus_loads if t in conns)
                    - Bs[idx,idx]*(vm0[idx]^2+2*vm0[idx]*(vm[t]-vm0[idx]))
                    ==
                    0.0
                )
                push!(cstr_q, cq)
            end
        end
    else
        for (idx,t) in ungrounded_terminals
            if any(Bs[idx,jdx] != 0 for (jdx, u) in ungrounded_terminals if idx != jdx) || any(Gs[idx,jdx] != 0 for (jdx, u) in ungrounded_terminals if idx != jdx)
                cp = JuMP.@constraint(pm.model,
                    sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
                    + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                    + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                    - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
                    + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
                    + sum( pd[l][t] for (l, conns) in bus_loads if t in conns)
                    + ( Gs[idx,idx]*(vm0[idx]^2+2*vm0[idx]*(vm[t]-vm0[idx]))
                        +sum( Gs[idx,jdx] * vm0[idx]*vm0[jdx] * cos(va0[idx]-va0[jdx])
                            +Bs[idx,jdx] * vm0[idx]*vm0[jdx] * sin(va0[idx]-va0[jdx])
                            +[Gs[idx,jdx]*vm0[jdx]*cos(va0[idx]-va0[jdx]) Gs[idx,jdx]*vm0[idx]*cos(va0[idx]-va0[jdx]) -Gs[idx,jdx]*vm0[idx]*vm0[jdx]*sin(va0[idx]-va0[jdx])  Gs[idx,jdx]*vm0[idx]*vm0[jdx]*sin(va0[idx]-va0[jdx])]*[vm[t]-vm0[idx];vm[u]-vm0[jdx];va[t]-va0[idx];va[u]-va0[jdx]]
                            +[Bs[idx,jdx]*vm0[jdx]*sin(va0[idx]-va0[jdx]) Bs[idx,jdx]*vm0[idx]*sin(va0[idx]-va0[jdx])  Bs[idx,jdx]*vm0[idx]*vm0[jdx]*cos(va0[idx]-va0[jdx]) -Bs[idx,jdx]*vm0[idx]*vm0[jdx]*cos(va0[idx]-va0[jdx])]*[vm[t]-vm0[idx];vm[u]-vm0[jdx];va[t]-va0[idx];va[u]-va0[jdx]]
                            for (jdx,u) in ungrounded_terminals if idx != jdx)
                    )
                    ==
                    0.0
                )
                push!(cstr_p, cp)

                cq = JuMP.@constraint(pm.model,
                    sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
                    + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                    + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                    - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
                    + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
                    + sum( qd[l][t] for (l, conns) in bus_loads if t in conns)
                    + ( -Bs[idx,idx]*(vm0[idx]^2+2*vm0[idx]*(vm[t]-vm0[idx]))
                        -sum( Bs[idx,jdx] * vm0[idx]*vm0[jdx] * cos(va0[idx]-va0[jdx])
                            -Gs[idx,jdx] * vm0[idx]*vm0[jdx] * sin(va0[idx]-va0[jdx])
                            +[Bs[idx,jdx]*vm0[jdx]*cos(va0[idx]-va0[jdx])   Bs[idx,jdx]*vm0[idx]*cos(va0[idx]-va0[jdx]) -Bs[idx,jdx]*vm0[idx]*vm0[jdx]*sin(va0[idx]-va0[jdx]) Bs[idx,jdx]*vm0[idx]*vm0[jdx]*sin(va0[idx]-va0[jdx])]*[vm[t]-vm0[idx];vm[u]-vm0[jdx];va[t]-va0[idx];va[u]-va0[jdx]]
                            +[-Gs[idx,jdx]*vm0[jdx]*sin(va0[idx]-va0[jdx]) -Gs[idx,jdx]*vm0[idx]*sin(va0[idx]-va0[jdx]) -Gs[idx,jdx]*vm0[idx]*vm0[jdx]*cos(va0[idx]-va0[jdx]) Gs[idx,jdx]*vm0[idx]*vm0[jdx]*cos(va0[idx]-va0[jdx])]*[vm[t]-vm0[idx];vm[u]-vm0[jdx];va[t]-va0[idx];va[u]-va0[jdx]]
                            for (jdx,u) in ungrounded_terminals if idx != jdx)
                    )
                    ==
                    0.0
                )
                push!(cstr_q, cq)
            else
                cp = JuMP.@constraint(pm.model,
                    sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
                    + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                    + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                    - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
                    + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
                    + sum( pd[l][t] for (l, conns) in bus_loads if t in conns)
                    + Gs[idx,idx]*(vm0[idx]^2+2*vm0[idx]*(vm[t]-vm0[idx]))
                    ==
                    0.0
                )
                push!(cstr_p, cp)

                cq = JuMP.@constraint(pm.model,
                    sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
                    + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                    + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                    - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
                    + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
                    + sum( qd[l][t] for (l, conns) in bus_loads if t in conns)
                    - Bs[idx,idx]*(vm0[idx]^2+2*vm0[idx]*(vm[t]-vm0[idx]))
                    ==
                    0.0
                )
                push!(cstr_q, cq)
            end
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
    constraint_mc_ohms_yt_from(pm::FOTPUPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})

Ohm constraints similar to ACPUPowerModel.
The nonlinear functions are approximated around initial operating points.
"""
function constraint_mc_ohms_yt_from(pm::FOTPUPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})
    p_fr  = var(pm, nw,  :p, f_idx)
    q_fr  = var(pm, nw,  :q, f_idx)
    vm_fr = var(pm, nw, :vm, f_bus)
    vm_to = var(pm, nw, :vm, t_bus)
    va_fr = var(pm, nw, :va, f_bus)
    va_to = var(pm, nw, :va, t_bus)
    vm0_fr = [var(pm, nw, :vm0, f_bus)[findfirst(isequal(c), ref(pm, nw, :bus, f_bus, "terminals"))] for c in f_connections]
    vm0_to = [var(pm, nw, :vm0, t_bus)[findfirst(isequal(c), ref(pm, nw, :bus, t_bus, "terminals"))] for c in t_connections]
    va0_fr = [var(pm, nw, :va0, f_bus)[findfirst(isequal(c), ref(pm, nw, :bus, f_bus, "terminals"))] for c in f_connections]
    va0_to = [var(pm, nw, :va0, t_bus)[findfirst(isequal(c), ref(pm, nw, :bus, t_bus, "terminals"))] for c in t_connections]

    for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))
        p_s_fr = []
        q_s_fr = []
        p_s_fr_dg = []
        q_s_fr_dg = []
        for (jdx, (fd,td)) in enumerate(zip(f_connections,t_connections))
            if idx != jdx
                p0_s_fr = (G[idx,jdx]+G_fr[idx,jdx])*vm0_fr[idx]*vm0_fr[jdx]*cos(va0_fr[idx]-va0_fr[jdx]) + (B[idx,jdx]+B_fr[idx,jdx])*vm0_fr[idx]*vm0_fr[jdx]*sin(va0_fr[idx]-va0_fr[jdx])
                q0_s_fr = (B[idx,jdx]+B_fr[idx,jdx])*vm0_fr[idx]*vm0_fr[jdx]*cos(va0_fr[idx]-va0_fr[jdx]) - (G[idx,jdx]+G_fr[idx,jdx])*vm0_fr[idx]*vm0_fr[jdx]*sin(va0_fr[idx]-va0_fr[jdx])
                push!(p_s_fr, JuMP.@expression(pm.model, p0_s_fr + (G[idx,jdx]+G_fr[idx,jdx])*(
                 (vm0_fr[jdx]*cos(va0_fr[idx]-va0_fr[jdx]))*(vm_fr[fc]-vm0_fr[idx])
                +(vm0_fr[idx]*cos(va0_fr[idx]-va0_fr[jdx]))*(vm_fr[fd]-vm0_fr[jdx])
                +(-vm0_fr[idx]*vm0_fr[jdx]*sin(va0_fr[idx]-va0_fr[jdx]))*(va_fr[fc]-va0_fr[idx])
                +( vm0_fr[idx]*vm0_fr[jdx]*sin(va0_fr[idx]-va0_fr[jdx]))*(va_fr[fd]-va0_fr[jdx])) + (B[idx,jdx]+B_fr[idx,jdx])*(
                 (vm0_fr[jdx]*sin(va0_fr[idx]-va0_fr[jdx]))*(vm_fr[fc]-vm0_fr[idx])
                +(vm0_fr[idx]*sin(va0_fr[idx]-va0_fr[jdx]))*(vm_fr[fd]-vm0_fr[jdx])
                +( vm0_fr[idx]*vm0_fr[jdx]*cos(va0_fr[idx]-va0_fr[jdx]))*(va_fr[fc]-va0_fr[idx])
                +(-vm0_fr[idx]*vm0_fr[jdx]*cos(va0_fr[idx]-va0_fr[jdx]))*(va_fr[fd]-va0_fr[jdx])))
                )
                push!(q_s_fr, JuMP.@expression(pm.model, q0_s_fr + (B[idx,jdx]+B_fr[idx,jdx])*(
                    (vm0_fr[jdx]*cos(va0_fr[idx]-va0_fr[jdx]))*(vm_fr[fc]-vm0_fr[idx])
                   +(vm0_fr[idx]*cos(va0_fr[idx]-va0_fr[jdx]))*(vm_fr[fd]-vm0_fr[jdx])
                   +(-vm0_fr[idx]*vm0_fr[jdx]*sin(va0_fr[idx]-va0_fr[jdx]))*(va_fr[fc]-va0_fr[idx])
                   +( vm0_fr[idx]*vm0_fr[jdx]*sin(va0_fr[idx]-va0_fr[jdx]))*(va_fr[fd]-va0_fr[jdx])) - (G[idx,jdx]+G_fr[idx,jdx])*(
                    (vm0_fr[jdx]*sin(va0_fr[idx]-va0_fr[jdx]))*(vm_fr[fc]-vm0_fr[idx])
                   +(vm0_fr[idx]*sin(va0_fr[idx]-va0_fr[jdx]))*(vm_fr[fd]-vm0_fr[jdx])
                   +( vm0_fr[idx]*vm0_fr[jdx]*cos(va0_fr[idx]-va0_fr[jdx]))*(va_fr[fc]-va0_fr[idx])
                   +(-vm0_fr[idx]*vm0_fr[jdx]*cos(va0_fr[idx]-va0_fr[jdx]))*(va_fr[fd]-va0_fr[jdx])))
                   )
            end
            p0_s_fr_dg = -G[idx,jdx]*vm0_fr[idx]*vm0_to[jdx]*cos(va0_fr[idx]-va0_to[jdx])-B[idx,jdx]*vm0_fr[idx]*vm0_to[jdx]*sin(va0_fr[idx]-va0_to[jdx])
            q0_s_fr_dg = -B[idx,jdx]*vm0_fr[idx]*vm0_to[jdx]*cos(va0_fr[idx]-va0_to[jdx])+G[idx,jdx]*vm0_fr[idx]*vm0_to[jdx]*sin(va0_fr[idx]-va0_to[jdx])
            push!(p_s_fr_dg, JuMP.@expression(pm.model, p0_s_fr_dg - G[idx,jdx]*(
                 (vm0_to[jdx]*cos(va0_fr[idx]-va0_to[jdx]))*(vm_fr[fc]-vm0_fr[idx])
                +(vm0_fr[idx]*cos(va0_fr[idx]-va0_to[jdx]))*(vm_to[td]-vm0_to[jdx])
                +(-vm0_fr[idx]*vm0_to[jdx]*sin(va0_fr[idx]-va0_to[jdx]))*(va_fr[fc]-va0_fr[idx])
                +( vm0_fr[idx]*vm0_to[jdx]*sin(va0_fr[idx]-va0_to[jdx]))*(va_to[td]-va0_to[jdx])) - B[idx,jdx]*(
                 (vm0_to[jdx]*sin(va0_fr[idx]-va0_to[jdx]))*(vm_fr[fc]-vm0_fr[idx])
                +(vm0_fr[idx]*sin(va0_fr[idx]-va0_to[jdx]))*(vm_to[td]-vm0_to[jdx])
                +( vm0_fr[idx]*vm0_to[jdx]*cos(va0_fr[idx]-va0_to[jdx]))*(va_fr[fc]-va0_fr[idx])
                +(-vm0_fr[idx]*vm0_to[jdx]*cos(va0_fr[idx]-va0_to[jdx]))*(va_to[td]-va0_to[jdx])))
                )
            push!(q_s_fr_dg, JuMP.@expression(pm.model, q0_s_fr_dg - B[idx,jdx]*(
                    (vm0_to[jdx]*cos(va0_fr[idx]-va0_to[jdx]))*(vm_fr[fc]-vm0_fr[idx])
                   +(vm0_fr[idx]*cos(va0_fr[idx]-va0_to[jdx]))*(vm_to[td]-vm0_to[jdx])
                   +(-vm0_fr[idx]*vm0_to[jdx]*sin(va0_fr[idx]-va0_to[jdx]))*(va_fr[fc]-va0_fr[idx])
                   +(vm0_fr[idx]*vm0_to[jdx]*sin(va0_fr[idx]-va0_to[jdx]))*(va_to[td]-va0_to[jdx])) + G[idx,jdx]*(
                    (vm0_to[jdx]*sin(va0_fr[idx]-va0_to[jdx]))*(vm_fr[fc]-vm0_fr[idx])
                   +(vm0_fr[idx]*sin(va0_fr[idx]-va0_to[jdx]))*(vm_to[td]-vm0_to[jdx])
                   +(vm0_fr[idx]*vm0_to[jdx]*cos(va0_fr[idx]-va0_to[jdx]))*(va_fr[fc]-va0_fr[idx])
                   +(-vm0_fr[idx]*vm0_to[jdx]*cos(va0_fr[idx]-va0_to[jdx]))*(va_to[td]-va0_to[jdx])))
                   )
        end
        JuMP.@constraint(pm.model, p_fr[fc] == (G[idx,idx]+G_fr[idx,idx])*(2*vm0_fr[idx]*vm_fr[fc]-vm0_fr[idx]^2)
            +sum(p_s_fr[j] for j=1:length(p_s_fr)) + sum(p_s_fr_dg[j] for j=1:length(p_s_fr_dg))
        )

        JuMP.@constraint(pm.model, q_fr[fc] == -(B[idx,idx]+B_fr[idx,idx])*(2*vm0_fr[idx]*vm_fr[fc]-vm0_fr[idx]^2)
            -sum(q_s_fr[j] for j=1:length(q_s_fr)) - sum(q_s_fr_dg[j] for j=1:length(q_s_fr_dg))
        )
    end
end


"""
    constraint_mc_ohms_yt_to(pm::FOTPUPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_to::Matrix{<:Real}, B_to::Matrix{<:Real})

Ohm constraints similar to ACPUPowerModel.
The nonlinear functions are approximated around initial operating points.
"""
function constraint_mc_ohms_yt_to(pm::FOTPUPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_to::Matrix{<:Real}, B_to::Matrix{<:Real})
    constraint_mc_ohms_yt_from(pm, nw, t_bus, f_bus, t_idx, f_idx, t_connections, f_connections, G, B, G_to, B_to)
end


"""
    constraint_mc_transformer_power_yy(pm::FOTPUPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)

Add all constraints required to model a two-winding, wye-wye connected transformer similar to ACPUPowerModel.
"""
function constraint_mc_transformer_power_yy(pm::FOTPUPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    transformer = ref(pm, nw, :transformer, trans_id)

    vm_fr = var(pm, nw, :vm, f_bus)
    vm_to = var(pm, nw, :vm, t_bus)
    va_fr = var(pm, nw, :va, f_bus)
    va_to = var(pm, nw, :va, t_bus)
    vm0_fr = [var(pm, nw, :vm0, f_bus)[findfirst(isequal(c), ref(pm, nw, :bus, f_bus, "terminals"))] for c in f_connections]
    vm0_to = [var(pm, nw, :vm0, t_bus)[findfirst(isequal(c), ref(pm, nw, :bus, t_bus, "terminals"))] for c in t_connections]
    va0_fr = [var(pm, nw, :va0, f_bus)[findfirst(isequal(c), ref(pm, nw, :bus, f_bus, "terminals"))] for c in f_connections]
    va0_to = [var(pm, nw, :va0, t_bus)[findfirst(isequal(c), ref(pm, nw, :bus, t_bus, "terminals"))] for c in t_connections]

    p_fr = [var(pm, nw, :pt, f_idx)[c] for c in f_connections]
    p_to = [var(pm, nw, :pt, t_idx)[c] for c in t_connections]
    q_fr = [var(pm, nw, :qt, f_idx)[c] for c in f_connections]
    q_to = [var(pm, nw, :qt, t_idx)[c] for c in t_connections]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))]

    for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))
        if tm_fixed[idx]
            JuMP.@constraint(pm.model, vm_fr[fc] == tm_scale*tm[idx]*vm_to[tc])
        else
            # transformer taps without regcontrol, tap variable not required in regcontrol formulation
            JuMP.@constraint(pm.model, vm_fr[fc] == tm_scale*tm[idx]*vm_to[tc])

            # with regcontrol
            if haskey(transformer,"controls")
                v_ref = transformer["controls"]["vreg"][idx]
                δ = transformer["controls"]["band"][idx]
                r = transformer["controls"]["r"][idx]
                x = transformer["controls"]["x"][idx]

                # linearized voltage: vm_drop = (r⋅p+x⋅q)/vm0
                vm_drop = JuMP.@expression(pm.model, (r*p_to[idx] + x*q_to[idx])/vm0_fr[idx])

                # v_ref-δ ≤ vm_fr-vm_drop ≤ v_ref+δ
                # vm_fr/1.1 ≤ vm_to ≤ vm_fr/0.9
                JuMP.@constraint(pm.model, vm_fr[fc] - vm_drop ≥ v_ref - δ)
                JuMP.@constraint(pm.model, vm_fr[fc] - vm_drop ≤ v_ref + δ)
                JuMP.@constraint(pm.model, vm_fr[fc]/1.1 ≤ vm_to[tc])
                JuMP.@constraint(pm.model, vm_fr[fc]/0.9 ≥ vm_to[tc])
            end
        end
        pol_angle = pol == 1 ? 0 : pi
        JuMP.@constraint(pm.model, va_fr[fc] == va_to[tc] + pol_angle)
    end

    JuMP.@constraint(pm.model, p_fr + p_to .== 0)
    JuMP.@constraint(pm.model, q_fr + q_to .== 0)
end


"""
    constraint_mc_transformer_power_dy(pm::FOTPUPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)

Add all constraints required to model a two-winding, delta-wye connected transformer similar to ACPUPowerModel
with voltage constraints linearized using first-order Taylor approximation and
power constraints simplified using initial operating point voltage instead of actual voltage variables.
"""
function constraint_mc_transformer_power_dy(pm::FOTPUPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    vm_fr = var(pm, nw, :vm, f_bus)
    vm_to = var(pm, nw, :vm, t_bus)
    va_fr = var(pm, nw, :va, f_bus)
    va_to = var(pm, nw, :va, t_bus)
    vm0_fr = [var(pm, nw, :vm0, f_bus)[findfirst(isequal(c), ref(pm, nw, :bus, f_bus, "terminals"))] for c in f_connections]
    vm0_to = [var(pm, nw, :vm0, t_bus)[findfirst(isequal(c), ref(pm, nw, :bus, t_bus, "terminals"))] for c in t_connections]
    va0_fr = [var(pm, nw, :va0, f_bus)[findfirst(isequal(c), ref(pm, nw, :bus, f_bus, "terminals"))] for c in f_connections]
    va0_to = [var(pm, nw, :va0, t_bus)[findfirst(isequal(c), ref(pm, nw, :bus, t_bus, "terminals"))] for c in t_connections]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[fc] for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))]

    nph = length(tm_set)

    # introduce auxialiary variable vd = Md*v_fr
    vd_re = Array{Any,1}(undef, nph)
    vd_im = Array{Any,1}(undef, nph)
    for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))
        # rotate by 1 to get 'previous' phase
        # e.g., for nph=3: 1->3, 2->1, 3->2
        jdx = (idx-1+1)%nph+1
        fd = f_connections[jdx]
        vd_re[idx] = JuMP.@expression(pm.model, vm0_fr[idx]*cos(va0_fr[idx])-vm0_fr[jdx]*cos(va0_fr[jdx])
                                             + (vm_fr[fc]-vm0_fr[idx])*cos(va0_fr[idx]) - (va_fr[fc]-va0_fr[idx])*sin(va0_fr[idx])
                                             - (vm_fr[fd]-vm0_fr[jdx])*cos(va0_fr[jdx]) + (va_fr[fd]-va0_fr[jdx])*sin(va0_fr[jdx])
                                    )
        vd_im[idx] = JuMP.@expression(pm.model, vm0_fr[idx]*sin(va0_fr[idx])-vm0_fr[jdx]*sin(va0_fr[jdx])
                                             + (vm_fr[fc]-vm0_fr[idx])*sin(va0_fr[idx]) + (va_fr[fc]-va0_fr[idx])*cos(va0_fr[idx])
                                             - (vm_fr[fd]-vm0_fr[jdx])*sin(va0_fr[jdx]) - (va_fr[fd]-va0_fr[jdx])*cos(va0_fr[jdx])
                                    )
        JuMP.@constraint(pm.model, vd_re[idx] == pol*tm_scale*tm[idx]*vm0_to[idx]*cos(va0_to[idx])
                                            + (vm_to[tc]-vm0_to[idx])*cos(va0_to[idx]) - (va_to[tc]-va0_to[idx])*sin(va0_to[idx]))
        JuMP.@constraint(pm.model, vd_im[idx] == pol*tm_scale*tm[idx]*vm0_to[idx]*sin(va0_to[idx])
                                            + (vm_to[tc]-vm0_to[idx])*sin(va0_to[idx]) + (va_to[tc]-va0_to[idx])*cos(va0_to[idx]))
    end

    p_fr = var(pm, nw, :pt, f_idx)
    p_to = var(pm, nw, :pt, t_idx)
    q_fr = var(pm, nw, :qt, f_idx)
    q_to = var(pm, nw, :qt, t_idx)

    id_re = Array{Any,1}(undef, nph)
    id_im = Array{Any,1}(undef, nph)

    for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))
        id_re[idx] = JuMP.@expression(pm.model,  (p_to[tc]*cos(va0_to[idx])+q_to[tc]*sin(va0_to[idx]))/vm0_to[idx]/(tm_scale*tm[idx])/pol)
        id_im[idx] = JuMP.@expression(pm.model, -(q_to[tc]*cos(va0_to[idx])-p_to[tc]*sin(va0_to[idx]))/vm0_to[idx]/(tm_scale*tm[idx])/pol)
    end
    for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))
        jdx = (idx-1+nph-1)%nph+1
        JuMP.@constraint(pm.model, p_fr[fc] ==
             vm0_fr[idx]*cos(va0_fr[idx])*( id_re[jdx]-id_re[idx])
            -vm0_fr[idx]*sin(va0_fr[idx])*(-id_im[jdx]+id_im[idx])
        )
        JuMP.@constraint(pm.model, q_fr[fc] ==
             vm0_fr[idx]*cos(va0_fr[idx])*(-id_im[jdx]+id_im[idx])
            +vm0_fr[idx]*sin(va0_fr[idx])*( id_re[jdx]-id_re[idx])
        )
    end
end


@doc raw"""
    constraint_mc_load_power(pm::FOTPUPowerModel, load_id::Int; nw::Int=nw_id_default, report::Bool=true)

Load model is linearized around initial operating point.
Wye loads are first-order and delta loads are zero-order approximations.

```math
\begin{align}
&\text{Initial operating point: }   v_{m0} \angle v_{a0}\\
&\text{Constant power: }  P^d = P^{d0},~Q^d = Q^{d0} \\
&\text{Constant impedance: }  P^d = a \cdot \left({v_{m0}}^2+2 \cdot v_{m0} \cdot (v_m-v_{m0})\right),\\
&  Q^d = b \cdot \left({v_{m0}}^2+2 \cdot v_{m0} \cdot (v_m-v_{m0})\right),  \\
&\text{Constant current: }  P^d = a \cdot v_m,\\
& Q^d = b \cdot v_m.
\end{align}
```
"""
function constraint_mc_load_power(pm::FOTPUPowerModel, load_id::Int; nw::Int=nw_id_default, report::Bool=true)
    # shared variables and parameters
    load = ref(pm, nw, :load, load_id)
    connections = load["connections"]
    bus_id = load["load_bus"]
    bus = ref(pm, nw, :bus, bus_id)

    # calculate load params
    a, alpha, b, beta = _load_expmodel_params(load, bus)

    # first-order approximation
    if load["configuration"]==WYE
        if load["model"]==POWER
            pd_bus = a
            qd_bus = b
        elseif load["model"]==IMPEDANCE
            vm = var(pm, nw, :vm, bus_id)
            vm0 = [var(pm, nw, :vm0, bus_id)[findfirst(isequal(c), bus["terminals"])] for c in connections]

            pd_bus = []
            qd_bus = []

            for (idx,c) in enumerate(connections)
                push!(pd_bus, a[idx]*(vm0[idx]^2 + 2*vm0[idx]*(vm[c]-vm0[idx])))
                push!(qd_bus, b[idx]*(vm0[idx]^2 + 2*vm0[idx]*(vm[c]-vm0[idx])))
            end
        else
            vm = var(pm, nw, :vm, bus_id)
            pd_bus = []
            qd_bus = []

            for (idx,c) in enumerate(connections)
                push!(pd_bus, a[idx]*vm[c])
                push!(qd_bus, b[idx]*vm[c])
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
        vm0 = var(pm, nw, :vm0, bus_id)
        va0 = var(pm, nw, :va0, bus_id)
        vr0 = vm0.*cos.(va0)
        vi0 = vm0.*sin.(va0)

        nph = length(a)
        is_triplex = nph < 3
        conn_bus = is_triplex ? ref(pm, nw, :bus, bus_id)["terminals"] : connections

        prev = Dict(c=>conn_bus[(idx+nph-2)%nph+1] for (idx,c) in enumerate(conn_bus))
        next = is_triplex ? conn_bus[2] : Dict(c=>conn_bus[idx%nph+1] for (idx,c) in enumerate(conn_bus))

        vrd0 = [vr0[idx]-vr0[next[idx]] for (idx, c) in enumerate(connections)]
        vid0 = [vi0[idx]-vi0[next[idx]] for (idx, c) in enumerate(connections)]

        crd0 = Array{Any,1}(undef, nph)
        cid0 = Array{Any,1}(undef, nph)
        for (idx, c) in enumerate(connections)
            crd0[c] = a[idx]*vrd0[c]*(vrd0[c]^2+vid0[c]^2)^(alpha[idx]/2-1)+b[idx]*vid0[c]*(vrd0[c]^2+vid0[c]^2)^(beta[idx]/2 -1)
            cid0[c] = a[idx]*vid0[c]*(vrd0[c]^2+vid0[c]^2)^(alpha[idx]/2-1)-b[idx]*vrd0[c]*(vrd0[c]^2+vid0[c]^2)^(beta[idx]/2 -1)
        end

        crd0_bus = is_triplex ? [(-1.0)^(c-1)*crd0[1] for (idx, c) in enumerate(conn_bus)] : [crd0[idx]-crd0[prev[idx]] for (idx, c) in enumerate(conn_bus)]
        cid0_bus = is_triplex ? [(-1.0)^(c-1)*cid0[1] for (idx, c) in enumerate(conn_bus)] : [cid0[idx]-cid0[prev[idx]] for (idx, c) in enumerate(conn_bus)]

        pd_bus = [ vr0[c]*crd0_bus[c]+vi0[c]*cid0_bus[c] for (idx,c) in enumerate(conn_bus)]
        qd_bus = [-vr0[c]*cid0_bus[c]+vi0[c]*crd0_bus[c] for (idx,c) in enumerate(conn_bus)]

        pd_bus = JuMP.Containers.DenseAxisArray(pd_bus, conn_bus)
        qd_bus = JuMP.Containers.DenseAxisArray(qd_bus, conn_bus)

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

@doc raw"""
    constraint_mc_generator_power_delta(pm::FOTPUPowerModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)

Adds constraints for delta-connected generators similar to delta-connected loads (zero-order approximation).

```math
\begin{align}
&\text{Initial line-neutral voltage: }   V_0 = V_{m0} \angle V_{a0}\\
&\text{Three-phase delta transformation matrix: }  M^\Delta = \begin{bmatrix}\;\;\;1 & -1 & \;\;0\\ \;\;\;0 & \;\;\;1 & -1\\ -1 & \;\;\;0 & \;\;\;1\end{bmatrix} \\
&\text{Single-phase delta transformation matrix (triple nodes): }  M^\Delta = \begin{bmatrix}\;1 & -1 \end{bmatrix} \\
&\text{Initial line-line voltage: }  V_0^\Delta = M^\Delta V_0 \\
&\text{Line-line current: }  (I^\Delta)^* = S^\Delta \oslash V_0^\Delta \\
&\text{Line-neutral current: }  I_{bus} = (M^\Delta)^T I^\Delta \\
&\text{Line-neutral generation power: }  S_{bus} = V_0 \oslash I_{bus}^*
\end{align}
```
"""
function constraint_mc_generator_power_delta(pm::FOTPUPowerModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vm0 = var(pm, nw, :vm0, bus_id)
    va0 = var(pm, nw, :va0, bus_id)
    vr0 = vm0.*cos.(va0)
    vi0 = vm0.*sin.(va0)
    pg = var(pm, nw, :pg, id)
    qg = var(pm, nw, :qg, id)

    nph = length(connections)
    is_triplex = nph < 3
    conn_bus = is_triplex ? ref(pm, nw, :bus, bus_id)["terminals"] : connections

    prev = Dict(c=>conn_bus[(idx+nph-2)%nph+1] for (idx,c) in enumerate(conn_bus))
    next = is_triplex ? conn_bus[2] : Dict(c=>conn_bus[idx%nph+1] for (idx,c) in enumerate(conn_bus))

    vrg0 = [vr0[idx]-vr0[next[idx]] for (idx, c) in enumerate(connections)]
    vig0 = [vi0[idx]-vi0[next[idx]] for (idx, c) in enumerate(connections)]

    crg = Dict() # Re(s/v)  = (p*vr+q*vi)/|v|^2
    cig = Dict() # Im(s/v) = -(q*vr-p*vi)/|v|^2
    for c in connections
        crg[c] = JuMP.@expression(pm.model, (pg[c]*vrg0[c]+qg[c]*vig0[c])/(vrg0[c]^2+vig0[c]^2) )
        cig[c] = JuMP.@expression(pm.model, (pg[c]*vig0[c]-qg[c]*vrg0[c])/(vrg0[c]^2+vig0[c]^2) )
    end

    crg_bus = Dict()
    cig_bus = Dict()
    for c in conn_bus
        if is_triplex
            crg_bus[c] = JuMP.@expression(pm.model, (-1.0)^(c-1)*crg[1])
            cig_bus[c] = JuMP.@expression(pm.model, (-1.0)^(c-1)*cig[1])
        else
            crg_bus[c] = JuMP.@expression(pm.model, crg[c]-crg[prev[c]])
            cig_bus[c] = JuMP.@expression(pm.model, cig[c]-cig[prev[c]])
        end
    end

    pg_bus = []
    qg_bus = []
    for c in conn_bus
        push!(pg_bus, JuMP.@expression(pm.model,  vr0[c]*crg_bus[c]+vi0[c]*cig_bus[c]))
        push!(qg_bus, JuMP.@expression(pm.model, -vr0[c]*cig_bus[c]+vi0[c]*crg_bus[c]))
    end

    pg_bus = JuMP.Containers.DenseAxisArray(pg_bus, conn_bus)
    qg_bus = JuMP.Containers.DenseAxisArray(qg_bus, conn_bus)

    var(pm, nw, :pg_bus)[id] = pg_bus
    var(pm, nw, :qg_bus)[id] = qg_bus

    if report
        sol(pm, nw, :gen, id)[:pg_bus] = pg_bus
        sol(pm, nw, :gen, id)[:qg_bus] = qg_bus
    end
end
