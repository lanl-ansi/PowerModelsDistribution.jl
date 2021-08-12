# The model in polar coordinates is linearized around an initial operating point using a 
# first order Taylor approximation (FOT) method


"""
    constraint_mc_load_current_delta(pm::AbstractLPUModel, nw::Int, load_id::Int, load_bus_id::Int, cp::Vector, cq::Vector)

No loads require a current variable. Delta loads are zero-order approximations and 
wye loads are first-order approximations around the initial operating point.
"""
function constraint_mc_load_current_delta(pm::AbstractLPUModel, nw::Int, load_id::Int, load_bus_id::Int, cp::Vector, cq::Vector)
end


"""
    variable_mc_bus_voltage(pm::AbstractLPUModel; nw=nw_id_default, kwargs...)

Voltage variables are defined in polar coordinates similar to ACPUPowerModel.
An initial operating point is specified for linearization.
"""
function variable_mc_bus_voltage(pm::AbstractLPUModel; nw=nw_id_default, kwargs...)
    variable_mc_bus_voltage_angle(pm; nw=nw, kwargs...)
    variable_mc_bus_voltage_magnitude_only(pm; nw=nw, kwargs...)

    # initial operating point for linearization (using flat-start)
    vm0 = var(pm, nw)[:vm0] = Dict(i => [1 1 1] for i in ids(pm, nw, :bus))
    va0 = var(pm, nw)[:va0] = Dict(i => [0 -2*pi/3 2*pi/3] for i in ids(pm, nw, :bus))

    for id in ids(pm, nw, :bus)
        busref = ref(pm, nw, :bus, id)
        terminals = busref["terminals"]
        grounded = busref["grounded"]

        ncnd = length(terminals)

        vm = haskey(busref, "vm_start") ? busref["vm_start"] : fill(0.0, ncnd)
        vm[.!grounded] .= 1.0

        # TODO how to do this more generally
        nph = 3
        va = haskey(busref, "va_start") ? busref["va_start"] : [c <= nph ? _wrap_to_pi(2 * pi / nph * (1-c)) : 0.0 for c in terminals]

        for (idx,t) in enumerate(terminals)
            JuMP.set_start_value(var(pm, nw, :vm, id)[t], vm[idx])
            JuMP.set_start_value(var(pm, nw, :va, id)[t], va[idx])
            # update initial operating point with warm-start
            var(pm, nw, :vm0, id)[t] = vm[idx]
            var(pm, nw, :va0, id)[t] = va[idx]
        end
    end
end


@doc raw"""
    constraint_mc_power_balance(pm::FOTUPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})

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
function constraint_mc_power_balance(pm::FOTUPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
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
                + ( Gs[idx,idx]*(vm0[t]^2+2*vm0[t]*(vm[t]-vm0[t]))
                    +sum( Gs[idx,jdx] * vm0[t]*vm0[u] * cos(va0[t]-va0[u])
                         +Bs[idx,jdx] * vm0[t]*vm0[u] * sin(va0[t]-va0[u])
                         +[Gs[idx,jdx]*vm0[u]*cos(va0[t]-va0[u]) Gs[idx,jdx]*vm0[t]*cos(va0[t]-va0[u]) -Gs[idx,jdx]*vm0[t]*vm0[u]*sin(va0[t]-va0[u])  Gs[idx,jdx]*vm0[t]*vm0[u]*sin(va0[t]-va0[u])]*[vm[t]-vm0[t];vm[u]-vm0[u];va[t]-va0[t];va[u]-va0[u]]
                         +[Bs[idx,jdx]*vm0[u]*sin(va0[t]-va0[u]) Bs[idx,jdx]*vm0[t]*sin(va0[t]-va0[u])  Bs[idx,jdx]*vm0[t]*vm0[u]*cos(va0[t]-va0[u]) -Bs[idx,jdx]*vm0[t]*vm0[u]*cos(va0[t]-va0[u])]*[vm[t]-vm0[t];vm[u]-vm0[u];va[t]-va0[t];va[u]-va0[u]]
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
                + ( -Bs[idx,idx]*(vm0[t]^2+2*vm0[t]*(vm[t]-vm0[t]))
                    -sum( Bs[idx,jdx] * vm0[t]*vm0[u] * cos(va0[t]-va0[u])
                         -Gs[idx,jdx] * vm0[t]*vm0[u] * sin(va0[t]-va0[u])
                         +[Bs[idx,jdx]*vm0[u]*cos(va0[t]-va0[u])   Bs[idx,jdx]*vm0[t]*cos(va0[t]-va0[u]) -Bs[idx,jdx]*vm0[t]*vm0[u]*sin(va0[t]-va0[u]) Bs[idx,jdx]*vm0[t]*vm0[u]*sin(va0[t]-va0[u])]*[vm[t]-vm0[t];vm[u]-vm0[u];va[t]-va0[t];va[u]-va0[u]]
                         +[-Gs[idx,jdx]*vm0[u]*sin(va0[t]-va0[u]) -Gs[idx,jdx]*vm0[t]*sin(va0[t]-va0[u]) -Gs[idx,jdx]*vm0[t]*vm0[u]*cos(va0[t]-va0[u]) Gs[idx,jdx]*vm0[t]*vm0[u]*cos(va0[t]-va0[u])]*[vm[t]-vm0[t];vm[u]-vm0[u];va[t]-va0[t];va[u]-va0[u]]
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
                + Gs[idx,idx]*(vm0[t]^2+2*vm0[t]*(vm[t]-vm0[t]))
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
                - Bs[idx,idx]*(vm0[t]^2+2*vm0[t]*(vm[t]-vm0[t]))
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


"""
    constraint_mc_ohms_yt_from(pm::FOTUPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})

Ohm constraints similar to ACPUPowerModel. 
The nonlinear functions are approximated around initial operating points.
"""
function constraint_mc_ohms_yt_from(pm::FOTUPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})
    p_fr  = var(pm, nw,  :p, f_idx)
    q_fr  = var(pm, nw,  :q, f_idx)
    vm_fr = var(pm, nw, :vm, f_bus)
    vm_to = var(pm, nw, :vm, t_bus)
    va_fr = var(pm, nw, :va, f_bus)
    va_to = var(pm, nw, :va, t_bus)
    vm0_fr = var(pm, nw, :vm0, f_bus)
    vm0_to = var(pm, nw, :vm0, t_bus)
    va0_fr = var(pm, nw, :va0, f_bus)
    va0_to = var(pm, nw, :va0, t_bus)

    for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))
        p_s_fr = []
        q_s_fr = []
        p_s_fr_dg = []
        q_s_fr_dg = []
        for (jdx, (fd,td)) in enumerate(zip(f_connections,t_connections)) 
            if idx != jdx
                p0_s_fr = (G[idx,jdx]+G_fr[idx,jdx])*vm0_fr[fc]*vm0_fr[fd]*cos(va0_fr[fc]-va0_fr[fd]) + (B[idx,jdx]+B_fr[idx,jdx])*vm0_fr[fc]*vm0_fr[fd]*sin(va0_fr[fc]-va0_fr[fd])
                q0_s_fr = (B[idx,jdx]+B_fr[idx,jdx])*vm0_fr[fc]*vm0_fr[fd]*cos(va0_fr[fc]-va0_fr[fd]) - (G[idx,jdx]+G_fr[idx,jdx])*vm0_fr[fc]*vm0_fr[fd]*sin(va0_fr[fc]-va0_fr[fd])
                push!(p_s_fr, JuMP.@expression(pm.model, p0_s_fr + (G[idx,jdx]+G_fr[idx,jdx])*(  
                 (vm0_fr[fd]*cos(va0_fr[fc]-va0_fr[fd]))*(vm_fr[fc]-vm0_fr[fc]) 
                +(vm0_fr[fc]*cos(va0_fr[fc]-va0_fr[fd]))*(vm_fr[fd]-vm0_fr[fd]) 
                +(-vm0_fr[fc]*vm0_fr[fd]*sin(va0_fr[fc]-va0_fr[fd]))*(va_fr[fc]-va0_fr[fc]) 
                +( vm0_fr[fc]*vm0_fr[fd]*sin(va0_fr[fc]-va0_fr[fd]))*(va_fr[fd]-va0_fr[fd])) + (B[idx,jdx]+B_fr[idx,jdx])*(
                 (vm0_fr[fd]*sin(va0_fr[fc]-va0_fr[fd]))*(vm_fr[fc]-vm0_fr[fc]) 
                +(vm0_fr[fc]*sin(va0_fr[fc]-va0_fr[fd]))*(vm_fr[fd]-vm0_fr[fd]) 
                +( vm0_fr[fc]*vm0_fr[fd]*cos(va0_fr[fc]-va0_fr[fd]))*(va_fr[fc]-va0_fr[fc]) 
                +(-vm0_fr[fc]*vm0_fr[fd]*cos(va0_fr[fc]-va0_fr[fd]))*(va_fr[fd]-va0_fr[fd])))
                )
                push!(q_s_fr, JuMP.@expression(pm.model, q0_s_fr + (B[idx,jdx]+B_fr[idx,jdx])*(  
                    (vm0_fr[fd]*cos(va0_fr[fc]-va0_fr[fd]))*(vm_fr[fc]-vm0_fr[fc]) 
                   +(vm0_fr[fc]*cos(va0_fr[fc]-va0_fr[fd]))*(vm_fr[fd]-vm0_fr[fd]) 
                   +(-vm0_fr[fc]*vm0_fr[fd]*sin(va0_fr[fc]-va0_fr[fd]))*(va_fr[fc]-va0_fr[fc]) 
                   +( vm0_fr[fc]*vm0_fr[fd]*sin(va0_fr[fc]-va0_fr[fd]))*(va_fr[fd]-va0_fr[fd])) - (G[idx,jdx]+G_fr[idx,jdx])*(
                    (vm0_fr[fd]*sin(va0_fr[fc]-va0_fr[fd]))*(vm_fr[fc]-vm0_fr[fc]) 
                   +(vm0_fr[fc]*sin(va0_fr[fc]-va0_fr[fd]))*(vm_fr[fd]-vm0_fr[fd]) 
                   +( vm0_fr[fc]*vm0_fr[fd]*cos(va0_fr[fc]-va0_fr[fd]))*(va_fr[fc]-va0_fr[fc]) 
                   +(-vm0_fr[fc]*vm0_fr[fd]*cos(va0_fr[fc]-va0_fr[fd]))*(va_fr[fd]-va0_fr[fd])))
                   )
            end
            p0_s_fr_dg = -G[idx,jdx]*vm0_fr[fc]*vm0_to[td]*cos(va0_fr[fc]-va0_to[td])-B[idx,jdx]*vm0_fr[fc]*vm0_to[td]*sin(va0_fr[fc]-va0_to[td])
            q0_s_fr_dg = -B[idx,jdx]*vm0_fr[fc]*vm0_to[td]*cos(va0_fr[fc]-va0_to[td])+G[idx,jdx]*vm0_fr[fc]*vm0_to[td]*sin(va0_fr[fc]-va0_to[td])
            push!(p_s_fr_dg, JuMP.@expression(pm.model, p0_s_fr_dg - G[idx,jdx]*(  
                 (vm0_to[td]*cos(va0_fr[fc]-va0_to[td]))*(vm_fr[fc]-vm0_fr[fc]) 
                +(vm0_fr[fc]*cos(va0_fr[fc]-va0_to[td]))*(vm_to[td]-vm0_to[td]) 
                +(-vm0_fr[fc]*vm0_to[td]*sin(va0_fr[fc]-va0_to[td]))*(va_fr[fc]-va0_fr[fc]) 
                +( vm0_fr[fc]*vm0_to[td]*sin(va0_fr[fc]-va0_to[td]))*(va_to[td]-va0_to[td])) - B[idx,jdx]*(
                 (vm0_to[td]*sin(va0_fr[fc]-va0_to[td]))*(vm_fr[fc]-vm0_fr[fc]) 
                +(vm0_fr[fc]*sin(va0_fr[fc]-va0_to[td]))*(vm_to[td]-vm0_to[td]) 
                +( vm0_fr[fc]*vm0_to[td]*cos(va0_fr[fc]-va0_to[td]))*(va_fr[fc]-va0_fr[fc]) 
                +(-vm0_fr[fc]*vm0_to[td]*cos(va0_fr[fc]-va0_to[td]))*(va_to[td]-va0_to[td])))
                )
            push!(q_s_fr_dg, JuMP.@expression(pm.model, q0_s_fr_dg - B[idx,jdx]*(  
                    (vm0_to[td]*cos(va0_fr[fc]-va0_to[td]))*(vm_fr[fc]-vm0_fr[fc]) 
                   +(vm0_fr[fc]*cos(va0_fr[fc]-va0_to[td]))*(vm_to[td]-vm0_to[td]) 
                   +(-vm0_fr[fc]*vm0_to[td]*sin(va0_fr[fc]-va0_to[td]))*(va_fr[fc]-va0_fr[fc]) 
                   +(vm0_fr[fc]*vm0_to[td]*sin(va0_fr[fc]-va0_to[td]))*(va_to[td]-va0_to[td])) + G[idx,jdx]*(
                    (vm0_to[td]*sin(va0_fr[fc]-va0_to[td]))*(vm_fr[fc]-vm0_fr[fc]) 
                   +(vm0_fr[fc]*sin(va0_fr[fc]-va0_to[td]))*(vm_to[td]-vm0_to[td]) 
                   +(vm0_fr[fc]*vm0_to[td]*cos(va0_fr[fc]-va0_to[td]))*(va_fr[fc]-va0_fr[fc]) 
                   +(-vm0_fr[fc]*vm0_to[td]*cos(va0_fr[fc]-va0_to[td]))*(va_to[td]-va0_to[td])))
                   )
        end
        JuMP.@constraint(pm.model, p_fr[fc] == (G[idx,idx]+G_fr[idx,idx])*(2*vm0_fr[fc]*vm_fr[fc]-vm0_fr[fc]^2)
            +sum(p_s_fr[j] for j=1:length(p_s_fr)) + sum(p_s_fr_dg[j] for j=1:length(p_s_fr_dg))
        )

        JuMP.@constraint(pm.model, q_fr[fc] == -(B[idx,idx]+B_fr[idx,idx])*(2*vm0_fr[fc]*vm_fr[fc]-vm0_fr[fc]^2)
            -sum(q_s_fr[j] for j=1:length(q_s_fr)) - sum(q_s_fr_dg[j] for j=1:length(q_s_fr_dg))
        )
    end
end


"""
    constraint_mc_ohms_yt_to(pm::FOTUPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_to::Matrix{<:Real}, B_to::Matrix{<:Real})

Ohm constraints similar to ACPUPowerModel. 
The nonlinear functions are approximated around initial operating points.
"""
function constraint_mc_ohms_yt_to(pm::FOTUPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_to::Matrix{<:Real}, B_to::Matrix{<:Real})
    constraint_mc_ohms_yt_from(pm, nw, t_bus, f_bus, t_idx, f_idx, t_connections, f_connections, G, B, G_to, B_to)
end


"""
    constraint_mc_transformer_power_yy(pm::FOTUPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)

Add all constraints required to model a two-winding, wye-wye connected transformer similar to ACPUPowerModel.
"""
function constraint_mc_transformer_power_yy(pm::FOTUPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    transformer = ref(pm, nw, :transformer, trans_id)

    vm_fr = var(pm, nw, :vm, f_bus)
    vm_to = var(pm, nw, :vm, t_bus)
    va_fr = var(pm, nw, :va, f_bus)
    va_to = var(pm, nw, :va, t_bus)
    vm0_fr = var(pm, nw, :vm0, f_bus)
    vm0_to = var(pm, nw, :vm0, t_bus)
    va0_fr = var(pm, nw, :va0, f_bus)
    va0_to = var(pm, nw, :va0, t_bus)

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
                vm_drop = JuMP.@expression(pm.model, (r*p_to[idx] + x*q_to[idx])/vm0_fr[fc])

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


@doc raw"""
    constraint_mc_transformer_power_dy(pm::FOTUPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)

Add all constraints required to model a two-winding, delta-wye connected transformer similar to ACPUPowerModel
with voltage constraints linearized using first-order Taylor approximation and
power constraints simplified using initial operating point voltage instead of actual voltage variables.

```math
\begin{align}
&\text{Initial operating point: }  v_{m0} \angle v_{a0}\\
& v_m \cdot \cos(v_a) \Rightarrow v_{m0} \cdot \cos(v_{a0}) + \cos(v_{a0}) \cdot(v_{m}-v_{m0}) -\sin(v_{a0}) \cdot(v_{a}-v_{a0})\\
& v_m \cdot \sin(v_a) \Rightarrow v_{m0} \cdot \sin(v_{a0}) + \sin(v_{a0}) \cdot(v_{m}-v_{m0}) +\cos(v_{a0}) \cdot(v_{a}-v_{a0})
\end{align}
The load power is
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
function constraint_mc_transformer_power_dy(pm::FOTUPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    vm_fr = var(pm, nw, :vm, f_bus)
    vm_to = var(pm, nw, :vm, t_bus)
    va_fr = var(pm, nw, :va, f_bus)
    va_to = var(pm, nw, :va, t_bus)
    vm0_fr = var(pm, nw, :vm0, f_bus)
    vm0_to = var(pm, nw, :vm0, t_bus)
    va0_fr = var(pm, nw, :va0, f_bus)
    va0_to = var(pm, nw, :va0, t_bus)

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
        vd_re[idx] = JuMP.@expression(pm.model, vm0_fr[fc]*cos(va0_fr[fc])-vm0_fr[fd]*cos(va0_fr[fd]) 
                                             + (vm_fr[fc]-vm0_fr[fc])*cos(va0_fr[fc]) - (va_fr[fc]-va0_fr[fc])*sin(va0_fr[fc])
                                             - (vm_fr[fd]-vm0_fr[fd])*cos(va0_fr[fd]) + (va_fr[fd]-va0_fr[fd])*sin(va0_fr[fd])
                                    )
        vd_im[idx] = JuMP.@expression(pm.model, vm0_fr[fc]*sin(va0_fr[fc])-vm0_fr[fd]*sin(va0_fr[fd]) 
                                             + (vm_fr[fc]-vm0_fr[fc])*sin(va0_fr[fc]) + (va_fr[fc]-va0_fr[fc])*cos(va0_fr[fc])
                                             - (vm_fr[fd]-vm0_fr[fd])*sin(va0_fr[fd]) - (va_fr[fd]-va0_fr[fd])*cos(va0_fr[fd])
                                    )
        JuMP.@constraint(pm.model, vd_re[idx] == pol*tm_scale*tm[idx]*vm0_to[tc]*cos(va0_to[tc]) 
                                            + (vm_to[tc]-vm0_to[tc])*cos(va0_to[tc]) - (va_to[tc]-va0_to[tc])*sin(va0_to[tc]))
        JuMP.@constraint(pm.model, vd_im[idx] == pol*tm_scale*tm[idx]*vm0_to[tc]*sin(va0_to[tc]) 
                                            + (vm_to[tc]-vm0_to[tc])*sin(va0_to[tc]) + (va_to[tc]-va0_to[tc])*cos(va0_to[tc]))
    end

    p_fr = var(pm, nw, :pt, f_idx)
    p_to = var(pm, nw, :pt, t_idx)
    q_fr = var(pm, nw, :qt, f_idx)
    q_to = var(pm, nw, :qt, t_idx)

    id_re = Array{Any,1}(undef, nph)
    id_im = Array{Any,1}(undef, nph)

    for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))
        id_re[idx] = JuMP.@expression(pm.model,  (p_to[tc]*cos(va0_to[tc])+q_to[tc]*sin(va0_to[tc]))/vm0_to[tc]/(tm_scale*tm[idx])/pol)
        id_im[idx] = JuMP.@expression(pm.model, -(q_to[tc]*cos(va0_to[tc])-p_to[tc]*sin(va0_to[tc]))/vm0_to[tc]/(tm_scale*tm[idx])/pol)
    end
    for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))
        jdx = (idx-1+nph-1)%nph+1
        JuMP.@constraint(pm.model, p_fr[fc] ==
             vm0_fr[fc]*cos(va0_fr[fc])*( id_re[jdx]-id_re[idx])
            -vm0_fr[fc]*sin(va0_fr[fc])*(-id_im[jdx]+id_im[idx])
        )
        JuMP.@constraint(pm.model, q_fr[fc] ==
             vm0_fr[fc]*cos(va0_fr[fc])*(-id_im[jdx]+id_im[idx])
            +vm0_fr[fc]*sin(va0_fr[fc])*( id_re[jdx]-id_re[idx])
        )
    end
end


@doc raw"""
    constraint_mc_load_power(pm::FOTUPowerModel, load_id::Int; nw::Int=nw_id_default, report::Bool=true)

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
function constraint_mc_load_power(pm::FOTUPowerModel, load_id::Int; nw::Int=nw_id_default, report::Bool=true)
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
            vm0 = var(pm, nw, :vm0, bus_id)

            pd_bus = []
            qd_bus = []

            for (idx,c) in enumerate(connections)
                push!(pd_bus, a[idx]*(vm0[c]^2 + 2*vm0[c]*(vm[c]-vm0[c])))
                push!(qd_bus, b[idx]*(vm0[c]^2 + 2*vm0[c]*(vm[c]-vm0[c])))
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

        nph = length(connections)
        prev = Dict(i=>(i+nph-2)%nph+1 for i in 1:nph)
        next = Dict(i=>i%nph+1 for i in 1:nph)

        vrd0 = [vr0[i]-vr0[next[i]] for i in 1:nph]
        vid0 = [vi0[i]-vi0[next[i]] for i in 1:nph]

        crd0 = [a[i]*vrd0[i]*(vrd0[i]^2+vid0[i]^2)^(alpha[i]/2-1)+b[i]*vid0[i]*(vrd0[i]^2+vid0[i]^2)^(beta[i]/2 -1) for i in 1:nph]
        cid0 = [a[i]*vid0[i]*(vrd0[i]^2+vid0[i]^2)^(alpha[i]/2-1)-b[i]*vrd0[i]*(vrd0[i]^2+vid0[i]^2)^(beta[i]/2 -1) for i in 1:nph]
        crd0_bus = [crd0[i]-crd0[prev[i]] for i in 1:nph]
        cid0_bus = [cid0[i]-cid0[prev[i]] for i in 1:nph]

        pd_bus = [ vr0[i]*crd0_bus[i]+vi0[i]*cid0_bus[i] for i in 1:nph]
        qd_bus = [-vr0[i]*cid0_bus[i]+vi0[i]*crd0_bus[i] for i in 1:nph] 
        var(pm, nw, :pd_bus)[load_id] = pd_bus
        var(pm, nw, :qd_bus)[load_id] = qd_bus
        
        if report
            sol(pm, nw, :load, load_id)[:pd_bus] = pd_bus
            sol(pm, nw, :load, load_id)[:qd_bus] = qd_bus

            pd = JuMP.@expression(pm.model, [i in 1:nph], a[i]*(vrd0[i]^2+vid0[i]^2)^(alpha[i]/2) )
            qd = JuMP.@expression(pm.model, [i in 1:nph], b[i]*(vrd0[i]^2+vid0[i]^2)^(beta[i]/2) )
          
            sol(pm, nw, :load, load_id)[:pd] = pd
            sol(pm, nw, :load, load_id)[:qd] = qd
        end
    end
end

