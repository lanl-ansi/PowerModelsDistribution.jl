# The branch flow model is linearized around an initial operating point using a single iteration
# of  the forward-backward sweep (FBS) method


"""
    variable_mc_load_power_delta_aux(pm::FBSUBFPowerModel, load_ids::Vector{Int}; nw::Int=nw_id_default, eps::Real=0.1, bounded::Bool=true, report::Bool=true)

Auxiliary variables are not required since delta loads are zero-order approximations
calculated using the initial operating point.
"""
function variable_mc_load_power_delta_aux(pm::FBSUBFPowerModel, load_ids::Vector{Int}; nw::Int=nw_id_default, eps::Real=0.1, bounded::Bool=true, report::Bool=true)
end


"""
    variable_mc_load_current(pm::FBSUBFPowerModel, load_ids::Vector{Int}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

No loads require a current variable. Delta loads are zero-order approximations and
wye loads are first-order approximations around the initial operating point.
"""
function variable_mc_load_current(pm::FBSUBFPowerModel, load_ids::Vector{Int}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
end


"""
    variable_mc_branch_power(pm::FBSUBFPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

Branch flow variables similar to LPUBFDiagModel
"""
function variable_mc_branch_power(pm::FBSUBFPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_branch_power_real(pm, nw=nw, bounded=bounded)
    variable_mc_branch_power_imaginary(pm, nw=nw, bounded=bounded)
end


"""
    variable_mc_bus_voltage(pm::FBSUBFPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

Voltage variables are defined in rectangular coordinates similar to ACRUPowerModel.
An initial operating point is specified for linearization.
"""
function variable_mc_bus_voltage(pm::FBSUBFPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_bus_voltage_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_bus_voltage_imaginary(pm; nw=nw, bounded=bounded, report=report)

    # initial operating point for linearization (using flat-start)
    var(pm, nw)[:vr0] = Dict{Int,Vector{Float64}}()  #(i => [cosd(0), cosd(-120), cosd(120)] for i in ids(pm, nw, :bus))
    var(pm, nw)[:vi0] = Dict{Int,Vector{Float64}}()  #(i => [sind(0), sind(-120), sind(120)] for i in ids(pm, nw, :bus))

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

        vm = haskey(busref, "vm_start") ? busref["vm_start"] : haskey(busref, "vm") ? busref["vm"] : [vm_start..., fill(0.0, ncnd)...][terminals]
        if busref["vbase"] == 0.12 # different angle initializations for center-tap transformer secondary nodes
            txfr_bus_id = [data["t_bus"] for (idx,data) in ref(pm, nw, :transformer) if data["f_bus"] == busref["index"]][1]
            txfr_bus = ref(pm, nw, :bus, txfr_bus_id)
            terminals = txfr_bus["terminals"]
            default_va = [[_wrap_to_pi(2 * pi / 3 * (1-t)) for t in 1:3]..., zeros(length(terminals))...][terminals][1]
            va = [default_va, _wrap_to_pi(default_va+pi)]
        else
            default_va = [[_wrap_to_pi(2 * pi / 3 * (1-t)) for t in 1:3]..., zeros(length(terminals))...][terminals]
            va = haskey(busref, "va_start") ? busref["va_start"] : haskey(busref, "va") ? busref["va"] : default_va
        end
        
        vr = vm.*cos.(va)
        vi = vm.*sin.(va)

        JuMP.set_start_value.(var(pm, nw, :vr, id), vr)
        JuMP.set_start_value.(var(pm, nw, :vi, id), vi)

        # update initial operating point with warm-start (causes infeasbility if not flat start)
        var(pm, nw, :vr0)[id] = busref["vbase"] == 0.12 ? vr : fill(1.0, ncnd) .* cos.(default_va) # vr
        var(pm, nw, :vi0)[id] = busref["vbase"] == 0.12 ? vi : fill(1.0, ncnd) .* sin.(default_va) # vi
    end
    # apply bounds if bounded
    if bounded
        for i in ids(pm, nw, :bus)
            constraint_mc_voltage_magnitude_bounds(pm, i, nw=nw)
        end
    end
end


"""
    variable_mc_capcontrol(pm::FBSUBFPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)

Capacitor switching variables.
"""
function variable_mc_capcontrol(pm::FBSUBFPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    variable_mc_capacitor_switch_state(pm; nw=nw, relax=relax, report=report)
end


"Creates variables for both `active` and `reactive` power flow at each transformer."
function variable_mc_transformer_power(pm::FBSUBFPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_transformer_power_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_transformer_power_imaginary(pm; nw=nw, bounded=bounded, report=report)
end


@doc raw"""
    constraint_mc_voltage_magnitude_bounds(pm::FBSUBFPowerModel, nw::Int, i::Int, vmin::Vector{<:Real}, vmax::Vector{<:Real})

Upper voltage magnitude limits are linearized using outer approximation.
Lower voltage magnitude limits are linearized around initial operating point.

```math
\begin{align}
&\text{Initial operating point: } ⇒ v_{r}^0 + j ⋅ v_{i}^0~\text{where}~{(v_m^0)}^2 = {(v_{r}^0)}^2 + {(v_{i}^0)}^2\\
&\text{Lower limits: }  2 ⋅ v_{r} ⋅ v_{r}^0 + 2 ⋅ v_{i} ⋅ v_{i}^0 - {(v_{m}^0)}^2 ≥ v_{min}^2,\\
&\text{Upper limits: } -v_{max} ≤  v_{r} ≤ v_{max},\\
& -v_{max} ≤  v_{i} ≤ v_{max},\\
&-\sqrt{2} ⋅ v_{max} ≤  v_{r} + v_{i} ≤ \sqrt{2} ⋅ v_{max},\\
& -\sqrt{2} ⋅ v_{max} ≤  v_{r} - v_{i} ≤ \sqrt{2} ⋅ v_{max}.
\end{align}
```
"""
function constraint_mc_voltage_magnitude_bounds(pm::FBSUBFPowerModel, nw::Int, i::Int, vmin::Vector{<:Real}, vmax::Vector{<:Real})
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
            JuMP.@constraint(pm.model, 2*(vr[t]*vr0[idx] + vi[t]*vi0[idx]) - (vr0[idx]^2 + vi0[idx]^2) >=  vmin[idx]^2)
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
    constraint_mc_voltage_angle_difference(pm::FBSUBFPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})

Nothing to do, this model ignores angle difference constraints"
"""
function constraint_mc_voltage_angle_difference(pm::FBSUBFPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})
end


@doc raw"""
    constraint_mc_power_losses(pm::FBSUBFPowerModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, r::Matrix{<:Real}, x::Matrix{<:Real}, g_sh_fr::Matrix{<:Real}, g_sh_to::Matrix{<:Real}, b_sh_fr::Matrix{<:Real}, b_sh_to::Matrix{<:Real})

Branch flow model power flow equation linearized around initial operating point (backward sweep)

```math
\begin{align}
&\text{Initial operating points: }  (v_{r0}^{fr} + j ⋅ v_{i0}^{fr}),~ (v_{r0}^{to} + j ⋅ v_{i0}^{to})\\
&\text{Voltage drop: }  v_{drop} = (v_{r0}^{fr} + j ⋅ v_{i0}^{fr}) - (v_{r0}^{to} + j ⋅ v_{i0}^{to}),\\
&\text{Line series admittance: } y = (r+j ⋅ x)^{-1},\\
&\text{Power loss: }  s_{loss} = v_{drop} ⋅ (y ⋅ v_{drop})^*,\\
&\text{Active power flow: }  p^{fr} + p^{to} = g_{sh}^{fr} ⋅ {(v_{m0}^{fr})}^2 +  g_{sh}^{to} ⋅ {(v_{m0}^{to})}^2 + \Re(s_{loss}),\\
&\text{Reactive power flow: }  q^{fr} + q^{to} = -b_{sh}^{fr} ⋅ {(v_{m0}^{fr})}^2 -  b_{sh}^{to} ⋅ {(v_{m0}^{to})}^2 + \Im(s_{loss}).
\end{align}
```
"""
function constraint_mc_power_losses(pm::FBSUBFPowerModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, r::Matrix{<:Real}, x::Matrix{<:Real}, g_sh_fr::Matrix{<:Real}, g_sh_to::Matrix{<:Real}, b_sh_fr::Matrix{<:Real}, b_sh_to::Matrix{<:Real})
    branch = ref(pm, nw, :branch, i)
    f_connections = branch["f_connections"]
    t_connections = branch["t_connections"]

    p_fr = var(pm, nw, :p)[f_idx]
    q_fr = var(pm, nw, :q)[f_idx]

    p_to = var(pm, nw, :p)[t_idx]
    q_to = var(pm, nw, :q)[t_idx]

    vr0_fr = [var(pm, nw, :vr0)[f_bus][findfirst(isequal(c), ref(pm, nw, :bus, f_bus, "terminals"))] for c in f_connections]
    vi0_fr = [var(pm, nw, :vi0)[f_bus][findfirst(isequal(c), ref(pm, nw, :bus, f_bus, "terminals"))] for c in f_connections]

    vr0_to = [var(pm, nw, :vr0)[t_bus][findfirst(isequal(c), ref(pm, nw, :bus, t_bus, "terminals"))] for c in t_connections]
    vi0_to = [var(pm, nw, :vi0)[t_bus][findfirst(isequal(c), ref(pm, nw, :bus, t_bus, "terminals"))] for c in t_connections]

    # linearized line voltage drop phasor: v_drop = v0_fr  - v0_to
    v_drop = vr0_fr .+ im.*vi0_fr .- vr0_to .- im.*vi0_to
    y = LinearAlgebra.pinv(r+im*x)
    N = length(f_connections)

    for (idx, (fc,tc)) in enumerate(zip(f_connections, t_connections))
        # linearized line power loss phasor: s_loss = v_drop*conj(y*v_drop)
        s_loss = v_drop[idx]*conj(sum(y[idx,j]*v_drop[j] for j=1:N))
        JuMP.@constraint(pm.model, p_fr[fc] + p_to[tc] ==  g_sh_fr[idx,idx]*(vr0_fr[idx]^2+vi0_fr[idx]^2) +  g_sh_to[idx,idx]*(vr0_to[idx]^2+vi0_to[idx]^2) + real(s_loss))
        JuMP.@constraint(pm.model, q_fr[fc] + q_to[tc] == -b_sh_fr[idx,idx]*(vr0_fr[idx]^2+vi0_fr[idx]^2) + -b_sh_to[idx,idx]*(vr0_to[idx]^2+vi0_to[idx]^2) + imag(s_loss))
    end
end


@doc raw"""
    constraint_mc_model_voltage_magnitude_difference(pm::FBSUBFPowerModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, r::Matrix{<:Real}, x::Matrix{<:Real}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real})

Voltage drop over a branch linearized around initial operating point (forward sweep)

```math
\begin{align}
&\text{Initial operating points: }  (v_{r0}^{fr} + j ⋅ v_{i0}^{fr}),~ (v_{r0}^{to} + j ⋅ v_{i0}^{to})\\
&\text{Series active power flow: }  p_s^{fr} =  p^{fr} -  g_{sh}^{fr} ⋅ {(v_{m0}^{fr})}^2,\\
&\text{Series reactive power flow: }  q_s^{fr} =  q^{fr} +  b_{sh}^{fr} ⋅ {(v_{m0}^{fr})}^2,\\
&\text{Series real current flow: }  cr_s^{fr} =  \frac{(p_s^{fr} ⋅ v_{r0}^{fr} + q_s^{fr} ⋅ v_{i0}^{fr})}{{(v_{m0}^{fr})}^2},\\
&\text{Series imaginary current flow: }  ci_s^{fr} =  \frac{(-q_s^{fr} ⋅ v_{r0}^{fr} + p_s^{fr} ⋅ v_{i0}^{fr})}{{(v_{m0}^{fr})}^2},\\
&\text{Series real voltage drop: } v_{r}^{to} = v_{r}^{fr} - r ⋅ cr_s^{fr} + x ⋅ ci_s^{fr} ,\\
&\text{Series imaginary voltage drop: } v_{i}^{to} = v_{i}^{fr} - x ⋅ cr_s^{fr} - r ⋅ ci_s^{fr}.
\end{align}
```
"""
function constraint_mc_model_voltage_magnitude_difference(pm::FBSUBFPowerModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, r::Matrix{<:Real}, x::Matrix{<:Real}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real})
    f_connections = ref(pm, nw, :branch, i)["f_connections"]
    t_connections = ref(pm, nw, :branch, i)["t_connections"]

    vr_fr = var(pm, nw, :vr)[f_bus]
    vi_fr = var(pm, nw, :vi)[f_bus]
    vr_to = var(pm, nw, :vr)[t_bus]
    vi_to = var(pm, nw, :vi)[t_bus]

    vr0_fr = [var(pm, nw, :vr0)[f_bus][findfirst(isequal(c), ref(pm, nw, :bus, f_bus, "terminals"))] for c in f_connections]
    vi0_fr = [var(pm, nw, :vi0)[f_bus][findfirst(isequal(c), ref(pm, nw, :bus, f_bus, "terminals"))] for c in f_connections]

    p_fr = var(pm, nw, :p)[f_idx]
    q_fr = var(pm, nw, :q)[f_idx]

    p_s_fr = [p_fr[fc]- LinearAlgebra.diag(g_sh_fr)[idx].*(vr0_fr[idx]^2+vi0_fr[idx]^2) for (idx,fc) in enumerate(f_connections)]
    q_s_fr = [q_fr[fc]+ LinearAlgebra.diag(b_sh_fr)[idx].*(vr0_fr[idx]^2+vi0_fr[idx]^2) for (idx,fc) in enumerate(f_connections)]

    # linearized current injection at from node: c0_s_fr = (p_s_fr - j.q_s_fr)/conj(v0_fr)
    cr_s_fr = [( p_s_fr[idx]*vr0_fr[idx] + q_s_fr[idx]*vi0_fr[idx])/(vr0_fr[idx]^2 + vi0_fr[idx]^2) for (idx,fc) in enumerate(f_connections)]
    ci_s_fr = [(-q_s_fr[idx]*vr0_fr[idx] + p_s_fr[idx]*vi0_fr[idx])/(vr0_fr[idx]^2 + vi0_fr[idx]^2) for (idx,fc) in enumerate(f_connections)]

    N = length(f_connections)
    for (idx, (fc, tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, vr_to[tc] == vr_fr[fc] - sum(r[idx,j]*cr_s_fr[j] for j in 1:N) + sum(x[idx,j]*ci_s_fr[j] for j in 1:N))
        JuMP.@constraint(pm.model, vi_to[tc] == vi_fr[fc] - sum(x[idx,j]*cr_s_fr[j] for j in 1:N) - sum(r[idx,j]*ci_s_fr[j] for j in 1:N))
    end
end


"""
    constraint_mc_theta_ref(pm::FBSUBFPowerModel, nw::Int, i::Int, va_ref::Vector{<:Real})

Creates phase angle constraints at reference buses similar to ACRUPowerModel.
"""
function constraint_mc_theta_ref(pm::FBSUBFPowerModel, nw::Int, i::Int, va_ref::Vector{<:Real})
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


"""
    constraint_mc_power_balance(pm::FBSUBFPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})

Power balance constraints similar to ACRUPowerModel with shunt current calculated using initial operating point.
"""
function constraint_mc_power_balance(pm::FBSUBFPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
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

    for (idx, t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
              sum(  p[arc][t] for (arc, conns) in bus_arcs if t in conns)
            + sum(psw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
            + sum( pt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
            ==
              sum(pg[gen][t] for (gen, conns) in bus_gens if t in conns)
            - sum(ps[strg][t] for (strg, conns) in bus_storage if t in conns)
            - sum(pd[load][t] for (load, conns) in bus_loads if t in conns)
            + ( -vr0[idx] * sum(Gt[idx,jdx]*vr0[jdx]-Bt[idx,jdx]*vi0[jdx] for (jdx,u) in ungrounded_terminals)
                -vi0[idx] * sum(Gt[idx,jdx]*vi0[jdx]+Bt[idx,jdx]*vr0[jdx] for (jdx,u) in ungrounded_terminals)
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
            + ( vr0[idx] * sum(Gt[idx,jdx]*vi0[jdx]+Bt[idx,jdx]*vr0[jdx] for (jdx,u) in ungrounded_terminals)
               -vi0[idx] * sum(Gt[idx,jdx]*vr0[jdx]-Bt[idx,jdx]*vi0[jdx] for (jdx,u) in ungrounded_terminals)
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


"""
    constraint_mc_power_balance_capc(pm::FBSUBFPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})

Power balance constraints with capacitor control similar to ACRUPowerModel with shunt current calculated using initial operating point.
"""
function constraint_mc_power_balance_capc(pm::FBSUBFPowerModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
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

    # add constraints to model capacitor switching
    if !isempty(bus_shunts) && haskey(ref(pm, nw, :shunt, bus_shunts[1][1]), "controls")
        constraint_capacitor_on_off(pm, nw, i, bus_shunts)
    end

    # calculate Gs, Bs
    ncnds = length(terminals)
    Gt = fill(0.0, ncnds, ncnds)
    Bt = convert(Matrix{JuMP.AffExpr}, JuMP.@expression(pm.model, [idx=1:ncnds, jdx=1:ncnds], 0.0))
    for (val, connections) in bus_shunts
        shunt = ref(pm,nw,:shunt,val)
        for (idx,c) in enumerate(connections)
            cap_state = haskey(shunt,"controls") ? var(pm, nw, :capacitor_state, val)[c] : 1.0
            for (jdx,d) in enumerate(connections)
                Gt[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] += shunt["gs"][idx,jdx]
                Bt[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] = JuMP.@expression(pm.model, Bt[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] + shunt["bs"][idx,jdx]*cap_state)
            end
        end
    end

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx, t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
              sum(  p[arc][t] for (arc, conns) in bus_arcs if t in conns)
            + sum(psw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
            + sum( pt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
            ==
              sum(pg[gen][t] for (gen, conns) in bus_gens if t in conns)
            - sum(ps[strg][t] for (strg, conns) in bus_storage if t in conns)
            - sum(pd[load][t] for (load, conns) in bus_loads if t in conns)
            + ( -vr0[idx] * sum(Gt[idx,jdx]*vr0[jdx]-Bt[idx,jdx]*vi0[jdx] for (jdx,u) in ungrounded_terminals)
                -vi0[idx] * sum(Gt[idx,jdx]*vi0[jdx]+Bt[idx,jdx]*vr0[jdx] for (jdx,u) in ungrounded_terminals)
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
            + ( vr0[idx] * sum(Gt[idx,jdx]*vi0[jdx]+Bt[idx,jdx]*vr0[jdx] for (jdx,u) in ungrounded_terminals)
               -vi0[idx] * sum(Gt[idx,jdx]*vr0[jdx]-Bt[idx,jdx]*vi0[jdx] for (jdx,u) in ungrounded_terminals)
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
    constraint_capacitor_on_off(pm::FBSUBFPowerModel, i::Int; nw::Int=nw_id_default)

Add constraints to model capacitor switching

```math
\begin{align}
&\text{kvar control (ON): }  q-q_\text{on} ≤ M_q ⋅ z - ϵ ⋅ (1-z), \\
&\text{kvar control (OFF): } q-q_\text{off} ≥ -M_q ⋅ (1-z) - ϵ ⋅ z, \\
&\text{voltage control (ON): }  2 ⋅ v_{r0} ⋅ v_r + 2 ⋅ v_{i0} ⋅ v_i - v_{r0}^2 - v_{i0}^2 - v_\text{min}^2 ≥ -M_v ⋅ z + ϵ ⋅ (1-z), \\
&\text{voltage control (OFF): } 2 ⋅ v_{r0} ⋅ v_r + 2 ⋅ v_{i0} ⋅ v_i - v_{r0}^2 - v_{i0}^2 - v_\text{max}^2 ≤ M_v ⋅ (1-z) - ϵ ⋅ z.
\end{align}
```
"""
function constraint_capacitor_on_off(pm::FBSUBFPowerModel, nw::Int, i::Int, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
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
                JuMP.@constraint(pm.model, 2*vr_cap*vr0_cap + 2*vi_cap*vi0_cap - vr0_cap^2 - vi0_cap^2 - shunt["controls"]["vmin"][idx]^2 ≥ -M_v*cap_state[val] + ϵ*(1-cap_state[val]))
                JuMP.@constraint(pm.model, 2*vr_cap*vr0_cap + 2*vi_cap*vi0_cap - vr0_cap^2 - vi0_cap^2 - shunt["controls"]["vmax"][idx]^2 ≤ M_v*(1-cap_state[val]) - ϵ*cap_state[val])
            end
            if shunt["controls"]["type"][idx] == CAP_DISABLED
                JuMP.@constraint(pm.model, cap_state[val] == 1 )
            end
        end
    end
end


@doc raw"""
    constraint_mc_load_power(pm::FBSUBFPowerModel, load_id::Int; nw::Int=nw_id_default, report::Bool=true)

Load model is linearized around initial operating point.
Wye loads are first-order and delta loads are zero-order approximations.

```math
\begin{align}
&\text{Initial operating point: }  v_{rd}^0 + j ⋅ v_{id}^0~\text{where}~{(v_m^0)}^2 = {(v_{rd}^0)}^2 + {(v_{id}^0)}^2\\
&\text{Constant power: }  P^d = P^{d0},~Q^d = Q^{d0} \\
&\text{Constant impedance: }  P^d = a ⋅ \left(2\cdot v_{rd} ⋅ v_{rd}^0+2 ⋅ v_{id}*v_{id}^0-{(v_{m}^0)}^2\right),\\
&  Q^d = b ⋅ \left(2\cdot v_{rd} ⋅ v_{rd}^0+2 ⋅ v_{id}*v_{id}^0-{(v_{m}^0)}^2\right),  \\
&\text{Constant current: }  P^d = a ⋅ \left(v_{m}^0 + \frac{v_{rd} ⋅ v_{rd}^0+ v_{id}*v_{id}^0-{(v_{m}^0)}^2}{v_{m}^0} \right),\\
& Q^d = b ⋅ \left(v_{m}^0 + \frac{v_{rd} ⋅ v_{rd}^0+ v_{id}*v_{id}^0-{(v_{m}^0)}^2}{v_{m}^0} \right).
\end{align}
```
"""
function constraint_mc_load_power(pm::FBSUBFPowerModel, load_id::Int; nw::Int=nw_id_default, report::Bool=true)
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
    constraint_mc_switch_state_closed(pm::FBSUBFPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int})

Voltage constraints for closed switches similar to ACRUPowerModel.
"""
function constraint_mc_switch_state_closed(pm::FBSUBFPowerModel, nw::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int})
    vr_fr = var(pm, nw, :vr, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)

    vi_fr = var(pm, nw, :vi, f_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    for (idx,(fc,tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, vr_fr[fc] == vr_to[tc])
        JuMP.@constraint(pm.model, vi_fr[fc] == vi_to[tc])
    end
end


"""
    constraint_mc_transformer_power_yy(pm::FBSUBFPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)

Add all constraints required to model a two-winding, wye-wye connected transformer similar to ACRUPowerModel.
"""
function constraint_mc_transformer_power_yy(pm::FBSUBFPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
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
                # JuMP.@constraint(pm.model, 2*vr0_fr[idx]*(vr_fr[fc]-vr_drop) + 2*vi0_fr[idx]*(vi_fr[fc]-vi_drop) - vr_fr[fc]^2 - vi_fr[fc]^2 ≥ (v_ref - δ)^2)
                JuMP.@constraint(pm.model, (2*vr_fr[fc]*vr0_fr[idx] + 2*vi_fr[fc]*vi0_fr[idx] - vr_fr[fc]^2 - vi_fr[fc]^2)/1.1^2 ≤ 2*vr_to[tc]*vr0_to[idx] + 2*vi_to[tc]*vi0_to[idx] - vr_to[tc]^2 - vi_to[tc]^2)
                JuMP.@constraint(pm.model, (2*vr_fr[fc]*vr0_fr[idx] + 2*vi_fr[fc]*vi0_fr[idx] - vr_fr[fc]^2 - vi_fr[fc]^2)/0.9^2 ≥ 2*vr_to[tc]*vr0_to[idx] + 2*vi_to[tc]*vi0_to[idx] - vr_to[tc]^2 - vi_to[tc]^2)
            end
        end
    end

    JuMP.@constraint(pm.model, p_fr + p_to .== 0)
    JuMP.@constraint(pm.model, q_fr + q_to .== 0)
end


"""
    constraint_mc_transformer_power_dy(pm::FBSUBFPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)

Add all constraints required to model a two-winding, delta-wye connected transformer similar to ACRUPowerModel
with power constraints using initial operating point voltage instead of actual voltage variables.
"""
function constraint_mc_transformer_power_dy(pm::FBSUBFPowerModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
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
