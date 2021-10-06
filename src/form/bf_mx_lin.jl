import LinearAlgebra: diagm


"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::AbstractLPUBFModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
end


""
function variable_mc_branch_current(pm::AbstractLPUBFModel; kwargs...)
end


""
function variable_mc_bus_voltage(pm::LPUBFDiagModel; nw::Int=nw_id_default, bounded::Bool=true)
    variable_mc_bus_voltage_magnitude_sqr(pm, nw=nw)
end


""
function variable_mc_branch_power(pm::LPUBFDiagModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_branch_power_real(pm, nw=nw, bounded=bounded)
    variable_mc_branch_power_imaginary(pm, nw=nw, bounded=bounded)
end


"""
    variable_mc_capcontrol(pm::AbstractLPUBFModel; nw::Int=nw_id_default, relax::Bool=false)

Capacitor switching and relaxed power variables.
"""
function variable_mc_capcontrol(pm::AbstractLPUBFModel; nw::Int=nw_id_default, relax::Bool=false)
    variable_mc_capacitor_switch_state(pm; nw=nw, relax=relax)
    variable_mc_capacitor_reactive_power(pm; nw=nw)
end


"Defines branch flow model power flow equations"
function constraint_mc_power_losses(pm::LPUBFDiagModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, r::Matrix{<:Real}, x::Matrix{<:Real}, g_sh_fr::Matrix{<:Real}, g_sh_to::Matrix{<:Real}, b_sh_fr::Matrix{<:Real}, b_sh_to::Matrix{<:Real})
    p_fr = var(pm, nw, :p)[f_idx]
    q_fr = var(pm, nw, :q)[f_idx]

    p_to = var(pm, nw, :p)[t_idx]
    q_to = var(pm, nw, :q)[t_idx]

    w_fr = var(pm, nw, :w)[f_bus]
    w_to = var(pm, nw, :w)[t_bus]


    fb = ref(pm, nw, :bus, f_bus)
    tb = ref(pm, nw, :bus, t_bus)

    branch = ref(pm, nw, :branch, i)
    f_connections = branch["f_connections"]
    t_connections = branch["t_connections"]

    for (idx, (fc,tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, p_fr[fc] + p_to[tc] == g_sh_fr[idx,idx]*w_fr[fc] +  g_sh_to[idx,idx]*w_to[tc])
        JuMP.@constraint(pm.model, q_fr[fc] + q_to[tc] == -b_sh_fr[idx,idx]*w_fr[fc] + -b_sh_to[idx,idx]*w_to[tc])
    end
end


"Defines voltage drop over a branch, linking from and to side voltage"
function constraint_mc_model_voltage_magnitude_difference(pm::LPUBFDiagModel, n::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, r::Matrix{<:Real}, x::Matrix{<:Real}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real})
    f_connections = ref(pm, n, :branch, i)["f_connections"]
    t_connections = ref(pm, n, :branch, i)["t_connections"]

    w_fr = var(pm, n, :w)[f_bus]
    w_to = var(pm, n, :w)[t_bus]

    p_fr = var(pm, n, :p)[f_idx]
    q_fr = var(pm, n, :q)[f_idx]

    p_s_fr = [p_fr[fc]- diag(g_sh_fr)[idx].*w_fr[fc] for (idx,fc) in enumerate(f_connections)]
    q_s_fr = [q_fr[fc]+ diag(b_sh_fr)[idx].*w_fr[fc] for (idx,fc) in enumerate(f_connections)]

    alpha = exp(-im*2*pi/3)
    Gamma = [1 alpha^2 alpha; alpha 1 alpha^2; alpha^2 alpha 1][f_connections,t_connections]

    MP = 2*(real(Gamma).*r + imag(Gamma).*x)
    MQ = 2*(real(Gamma).*x - imag(Gamma).*r)

    N = length(f_connections)

    for (idx, (fc, tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, w_to[tc]== w_fr[fc] - sum(MP[idx,j]*p_s_fr[j] for j in 1:N) - sum(MQ[idx,j]*q_s_fr[j] for j in 1:N))
    end
end


"balanced three-phase phasor"
function constraint_mc_theta_ref(pm::LPUBFDiagModel, nw::Int, i::Int, va_ref::Vector{<:Real})
    w = [var(pm, nw, :w, i)[t] for t in ref(pm, nw, :bus, i)["terminals"]]

    JuMP.@constraint(pm.model, w[2:end] .== w[1])
end


"Creates variables for both `active` and `reactive` power flow at each transformer."
function variable_mc_transformer_power(pm::LPUBFDiagModel; kwargs...)
    variable_mc_transformer_power_real(pm; kwargs...)
    variable_mc_transformer_power_imaginary(pm; kwargs...)
end


""
function constraint_mc_power_balance(pm::LPUBFDiagModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    w = var(pm, nw, :w, i)
    p   = get(var(pm, nw),      :p,   Dict()); _check_var_keys(p,   bus_arcs, "active power", "branch")
    q   = get(var(pm, nw),      :q,   Dict()); _check_var_keys(q,   bus_arcs, "reactive power", "branch")
    psw = get(var(pm, nw),    :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw = get(var(pm, nw),    :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt  = get(var(pm, nw),     :pt,  Dict()); _check_var_keys(pt,  bus_arcs_trans, "active power", "transformer")
    qt  = get(var(pm, nw),     :qt,  Dict()); _check_var_keys(qt,  bus_arcs_trans, "reactive power", "transformer")
    pg  = get(var(pm, nw),     :pg,  Dict()); _check_var_keys(pg,  bus_gens, "active power", "generator")
    qg  = get(var(pm, nw),     :qg,  Dict()); _check_var_keys(qg,  bus_gens, "reactive power", "generator")
    ps  = get(var(pm, nw),     :ps,  Dict()); _check_var_keys(ps,  bus_storage, "active power", "storage")
    qs  = get(var(pm, nw),     :qs,  Dict()); _check_var_keys(qs,  bus_storage, "reactive power", "storage")
    pd  = get(var(pm, nw), :pd_bus,  Dict()); _check_var_keys(pd,  bus_loads, "active power", "load")
    qd  = get(var(pm, nw), :qd_bus,  Dict()); _check_var_keys(qd,  bus_loads, "reactive power", "load")

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
              sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
            + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
            - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
            + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
            + sum( pd[d][t] for (d, conns) in bus_loads if t in conns)
            + sum(diag(ref(pm, nw, :shunt, sh, "gs"))[findfirst(isequal(t), conns)]*w[t] for (sh, conns) in bus_shunts if t in conns)
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
            + sum( qd[d][t] for (d, conns) in bus_loads if t in conns)
            - sum(diag(ref(pm, nw, :shunt, sh, "bs"))[findfirst(isequal(t), conns)]*w[t] for (sh, conns) in bus_shunts if t in conns)
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


@doc raw"""
    constraint_mc_power_balance_capc(pm::LPUBFDiagModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})

Power balance constraints with capacitor control linearized using McCormick envelopes

```math
\begin{align}
    & B_s = b_s ⋅ z,~~ cq_{sh} = B_s ⋅ w, \\
    &\text{Underestimator: }  cq_{sh} ≥ B_s ⋅ w_\text{min},~~ cq_{sh} ≥ b_s ⋅ w  + B_s ⋅ w_\text{max} - b_s ⋅ w_\text{max}\\
    &\text{Overestimator: }   cq_{sh} ≤ B_s ⋅ w_\text{max},~~ cq_{sh} ≤ b_s ⋅ w  + B_s ⋅ w_\text{min} - b_s ⋅ w_\text{min}\\
\end{align}
```
"""
function constraint_mc_power_balance_capc(pm::LPUBFDiagModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    w = var(pm, nw, :w, i)
    p   = get(var(pm, nw),      :p,   Dict()); _check_var_keys(p,   bus_arcs, "active power", "branch")
    q   = get(var(pm, nw),      :q,   Dict()); _check_var_keys(q,   bus_arcs, "reactive power", "branch")
    psw = get(var(pm, nw),    :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw = get(var(pm, nw),    :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt  = get(var(pm, nw),     :pt,  Dict()); _check_var_keys(pt,  bus_arcs_trans, "active power", "transformer")
    qt  = get(var(pm, nw),     :qt,  Dict()); _check_var_keys(qt,  bus_arcs_trans, "reactive power", "transformer")
    pg  = get(var(pm, nw),     :pg,  Dict()); _check_var_keys(pg,  bus_gens, "active power", "generator")
    qg  = get(var(pm, nw),     :qg,  Dict()); _check_var_keys(qg,  bus_gens, "reactive power", "generator")
    ps  = get(var(pm, nw),     :ps,  Dict()); _check_var_keys(ps,  bus_storage, "active power", "storage")
    qs  = get(var(pm, nw),     :qs,  Dict()); _check_var_keys(qs,  bus_storage, "reactive power", "storage")
    pd  = get(var(pm, nw), :pd_bus,  Dict()); _check_var_keys(pd,  bus_loads, "active power", "load")
    qd  = get(var(pm, nw), :qd_bus,  Dict()); _check_var_keys(qd,  bus_loads, "reactive power", "load")

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
              sum(  p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
            + sum( pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
            - sum( pg[g][t] for (g, conns) in bus_gens if t in conns)
            + sum( ps[s][t] for (s, conns) in bus_storage if t in conns)
            + sum( pd[d][t] for (d, conns) in bus_loads if t in conns)
            + sum(diag(ref(pm, nw, :shunt, sh, "gs"))[findfirst(isequal(t), conns)]*w[t] for (sh, conns) in bus_shunts if t in conns)
            ==
            0.0
        )
        push!(cstr_p, cp)
    end

    # add constraints to model capacitor switching
    if !isempty(bus_shunts) && haskey(ref(pm, nw, :shunt, bus_shunts[1][1]), "controls")
        constraint_capacitor_on_off(pm, i, bus_shunts)

        ncnds = length(bus_shunts[1][2])
        cq_sh = convert(Vector{JuMP.AffExpr}, JuMP.@expression(pm.model, [idx=1:ncnds], 0.0))
        for (idx,t) in ungrounded_terminals
            for (sh, conns) in bus_shunts
                if t in conns
                    cq_cap = var(pm, nw, :capacitor_reactive_power, sh)[t]
                    cap_state = var(pm, nw, :capacitor_state, sh)[t]
                    bs = diag(ref(pm, nw, :shunt, sh, "bs"))[findfirst(isequal(t), conns)]
                    w_min = 0.9^2
                    w_max = 1.1^2
                    # McCormick envelope constraints
                    JuMP.@constraint(pm.model, cq_cap ≥ bs*cap_state*w_min)
                    JuMP.@constraint(pm.model, cq_cap ≥ bs*w[t] + bs*cap_state*w_max - bs*w_max)
                    JuMP.@constraint(pm.model, cq_cap ≤ bs*cap_state*w_max)
                    JuMP.@constraint(pm.model, cq_cap ≤ bs*w[t] + bs*cap_state*w_min - bs*w_min)
                end
            end
            cq = JuMP.@constraint(pm.model,
                sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( qd[d][t] for (d, conns) in bus_loads if t in conns)
                - sum(var(pm, nw, :capacitor_reactive_power, bus_shunts[1][1]))
                ==
                0.0
            )
            push!(cstr_q, cq)
        end
    else
        for (idx,t) in ungrounded_terminals
            cq = JuMP.@constraint(pm.model,
                sum(  q[a][t] for (a, conns) in bus_arcs if t in conns)
                + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
                + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
                - sum( qg[g][t] for (g, conns) in bus_gens if t in conns)
                + sum( qs[s][t] for (s, conns) in bus_storage if t in conns)
                + sum( qd[d][t] for (d, conns) in bus_loads if t in conns)
                - sum(diag(ref(pm, nw, :shunt, sh, "bs"))[findfirst(isequal(t), conns)]*w[t] for (sh, conns) in bus_shunts if t in conns)
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
    constraint_capacitor_on_off(pm::LPUBFDiagModel, i::Int; nw::Int=nw_id_default)

Add constraints to model capacitor switching

```math
\begin{align}
&\text{kvar control (ON): }  q-q_\text{on} ≤ M_q ⋅ z - ϵ ⋅ (1-z), \\
&\text{kvar control (OFF): } q-q_\text{off} ≥ -M_q ⋅ (1-z) - ϵ ⋅ z, \\
&\text{voltage control (ON): }  w - v_\text{min}^2 ≥ -M_v ⋅ z + ϵ ⋅ (1-z), \\
&\text{voltage control (OFF): } w - v_\text{max}^2 ≤ M_v ⋅ (1-z) - ϵ ⋅ z.
\end{align}
```
"""
function constraint_capacitor_on_off(pm::LPUBFDiagModel, i::Int, bus_shunts::Vector{Tuple{Int,Vector{Int}}}; nw::Int=nw_id_default)
    cap_state = var(pm, nw, :capacitor_state, bus_shunts[1][1])
    shunt = ref(pm, nw, :shunt, bus_shunts[1][1])
    ϵ = 1e-5
    M_q = 1e5
    M_v = 2
    if shunt["controls"]["type"] == CAP_REACTIVE_POWER
        bus_idx = shunt["controls"]["terminal"] == 1 ? (shunt["controls"]["element"]["index"], shunt["controls"]["element"]["f_bus"], shunt["controls"]["element"]["t_bus"]) : (shunt["controls"]["element"]["index"], shunt["controls"]["element"]["t_bus"], shunt["controls"]["element"]["f_bus"])
        q_fr =  shunt["controls"]["element"]["type"] == "branch" ? var(pm, nw, :q)[bus_idx] : var(pm, nw, :qt, bus_idx)
        JuMP.@constraint(pm.model, sum(q_fr) - shunt["controls"]["onsetting"] ≤ M_q*cap_state[shunt["connections"][1]] - ϵ*(1-cap_state[shunt["connections"][1]]))
        JuMP.@constraint(pm.model, sum(q_fr) - shunt["controls"]["offsetting"] ≥ -M_q*(1-cap_state[shunt["connections"][1]]) - ϵ*cap_state[shunt["connections"][1]])
        JuMP.@constraint(pm.model, cap_state .== cap_state[shunt["connections"][1]])
        if shunt["controls"]["voltoverride"]
            for (idx,val) in enumerate(shunt["connections"])
                w = var(pm, nw, :w, i)[val]
                JuMP.@constraint(pm.model, w - shunt["controls"]["vmin"]^2 ≥ -M_v*cap_state[val] + ϵ*(1-cap_state[val]))
                JuMP.@constraint(pm.model, w - shunt["controls"]["vmax"]^2 ≤ M_v*(1-cap_state[val]) - ϵ*cap_state[val])
            end
        end
    else
        for (idx,val) in enumerate(shunt["connections"])
            if shunt["controls"]["type"][idx] == CAP_VOLTAGE
                bus_idx = shunt["controls"]["terminal"][idx] == 1 ? shunt["controls"]["element"]["f_bus"] : shunt["controls"]["element"]["t_bus"]
                w = var(pm, nw, :w, i)[val]
                JuMP.@constraint(pm.model, w - shunt["controls"]["onsetting"][idx]^2 ≤ M_v*cap_state[val] - ϵ*(1-cap_state[val]))
                JuMP.@constraint(pm.model, w - shunt["controls"]["offsetting"][idx]^2 ≥ -M_v*(1-cap_state[val]) - ϵ*cap_state[val])
            end
            if shunt["controls"]["voltoverride"][idx]
                w = var(pm, nw, :w, i)[val]
                JuMP.@constraint(pm.model, w - shunt["controls"]["vmin"][idx]^2 ≥ -M_v*cap_state[val] + ϵ*(1-cap_state[val]))
                JuMP.@constraint(pm.model, w - shunt["controls"]["vmax"][idx]^2 ≤ M_v*(1-cap_state[val]) - ϵ*cap_state[val])
            end
            if shunt["controls"]["type"][idx] == CAP_DISABLED
                JuMP.@constraint(pm.model, cap_state[val] == 1 )
            end
        end
    end
end


"Neglects the active and reactive loss terms associated with the squared current magnitude."
function constraint_mc_storage_losses(pm::AbstractUBFModels, i::Int; nw::Int=nw_id_default, kwargs...)
    storage = ref(pm, nw, :storage, i)

    p_loss, q_loss = storage["p_loss"], storage["q_loss"]
    conductors = storage["connections"]

    ps = var(pm, nw, :ps, i)
    qs = var(pm, nw, :qs, i)
    sc = var(pm, nw, :sc, i)
    sd = var(pm, nw, :sd, i)
    qsc = var(pm, nw, :qsc, i)


    JuMP.@constraint(pm.model,
        sum(ps[c] for c in conductors) + (sd - sc)
        ==
        p_loss
    )

    JuMP.@constraint(pm.model,
        sum(qs[c] for c in conductors)
        ==
        qsc + q_loss
    )
end

@doc raw"""
    constraint_mc_load_power(pm::LPUBFDiagModel, load_id::Int; nw::Int=nw_id_default, report::Bool=true)

Delta/voltage-dependent load models for LPUBFDiagModel. Delta loads use the auxilary power variable (X). The constant current load model is derived by linearizing around the flat-start voltage solution.

```math
\begin{align}
&\text{Constant power:} \Rightarrow P_i^d = P_i^{d0},~Q_i^d = Q_i^{d0} ~\forall i \in L \\
&\text{Constant impedance (Wye):} \Rightarrow P_i^d = a_i \cdot w_i,~Q_i^d = b_i \cdot w_i ~\forall i \in L \\
&\text{Constant impedance (Delta):} \Rightarrow P_i^d = 3\cdot a_i \cdot w_i,~Q_i^d = 3\cdot b_i \cdot w_i ~\forall i \in L \\
&\text{Constant current (Wye):} \Rightarrow P_i^d = \frac{a_i}{2}\cdot \left( 1+w_i \right),~Q_i^d = \frac{b_i}{2}\cdot \left( 1+w_i \right) \forall i \in L \\
&\text{Constant current (Delta):} \Rightarrow P_i^d = \frac{\sqrt{3} \cdot a_i}{2}\cdot \left( 1+w_i \right),~Q_i^d = \frac{\sqrt{3} \cdot b_i}{2}\cdot \left( 1+w_i \right) \forall i \in L
\end{align}
```
"""
function constraint_mc_load_power(pm::LPUBFDiagModel, load_id::Int; nw::Int=nw_id_default, report::Bool=true)
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
    vmin, vmax = _calc_load_vbounds(load, bus)
    wmin = vmin.^2
    wmax = vmax.^2
    pmin, pmax, qmin, qmax = _calc_load_pq_bounds(load, bus)

    # take care of connections
    if load["configuration"]==WYE
        if load["model"]==POWER
            var(pm, nw, :pd)[load_id] = JuMP.Containers.DenseAxisArray(pd0, connections)
            var(pm, nw, :qd)[load_id] = JuMP.Containers.DenseAxisArray(qd0, connections)
        elseif load["model"]==IMPEDANCE
            w = var(pm, nw, :w)[bus_id][[c for c in connections]]
            var(pm, nw, :pd)[load_id] = a.*w
            var(pm, nw, :qd)[load_id] = b.*w
        # in this case, :pd has a JuMP variable
        else
            w = var(pm, nw, :w)[bus_id][[c for c in connections]]
            pd = var(pm, nw, :pd)[load_id]
            qd = var(pm, nw, :qd)[load_id]
            for (idx,c) in enumerate(connections)
                JuMP.@constraint(pm.model, pd[c]==1/2*a[idx]*(w[c]+1))
                JuMP.@constraint(pm.model, qd[c]==1/2*b[idx]*(w[c]+1))
            end
        end
        # :pd_bus is identical to :pd now
        var(pm, nw, :pd_bus)[load_id] = var(pm, nw, :pd)[load_id]
        var(pm, nw, :qd_bus)[load_id] = var(pm, nw, :qd)[load_id]

        ## reporting
        if report
            sol(pm, nw, :load, load_id)[:pd] = var(pm, nw, :pd)[load_id]
            sol(pm, nw, :load, load_id)[:qd] = var(pm, nw, :qd)[load_id]
            sol(pm, nw, :load, load_id)[:pd_bus] = var(pm, nw, :pd_bus)[load_id]
            sol(pm, nw, :load, load_id)[:qd_bus] = var(pm, nw, :qd_bus)[load_id]
        end
    elseif load["configuration"]==DELTA
        Xdr = var(pm, nw, :Xdr, load_id)
        Xdi = var(pm, nw, :Xdi, load_id)
        Td = [1 -1 0; 0 1 -1; -1 0 1]  # TODO
        # define pd/qd and pd_bus/qd_bus as affine transformations of X
        pd_bus = LinearAlgebra.diag(Xdr*Td)
        qd_bus = LinearAlgebra.diag(Xdi*Td)
        pd = LinearAlgebra.diag(Td*Xdr)
        qd = LinearAlgebra.diag(Td*Xdi)
        # Equate missing edge parameters to zero
        for (idx, c) in enumerate(connections)
            if abs(pd0[idx]+im*qd0[idx]) == 0.0
                JuMP.@constraint(pm.model, Xdr[:,idx] .== 0)
                JuMP.@constraint(pm.model, Xdi[:,idx] .== 0)
            end
        end

        var(pm, nw, :pd_bus)[load_id] = pd_bus
        var(pm, nw, :qd_bus)[load_id] = qd_bus
        var(pm, nw, :pd)[load_id] = pd
        var(pm, nw, :qd)[load_id] = qd
        if load["model"]==POWER
            for (idx, c) in enumerate(connections)
                JuMP.@constraint(pm.model, pd[idx]==pd0[idx])
                JuMP.@constraint(pm.model, qd[idx]==qd0[idx])
            end
        elseif load["model"]==IMPEDANCE
            w = var(pm, nw, :w)[bus_id][[c for c in connections]]
            for (idx,c) in enumerate(connections)
                JuMP.@constraint(pm.model, pd[idx]==3*a[idx]*w[idx])
                JuMP.@constraint(pm.model, qd[idx]==3*b[idx]*w[idx])
            end
        else
            w = var(pm, nw, :w)[bus_id][[c for c in connections]]
            for (idx,c) in enumerate(connections)
                JuMP.@constraint(pm.model, pd[idx]==sqrt(3)/2*a[idx]*(w[idx]+1))
                JuMP.@constraint(pm.model, qd[idx]==sqrt(3)/2*b[idx]*(w[idx]+1))
            end
        end

        ## reporting; for delta these are not available as saved variables!
        if report
            sol(pm, nw, :load, load_id)[:pd] = pd
            sol(pm, nw, :load, load_id)[:qd] = qd
            sol(pm, nw, :load, load_id)[:pd_bus] = pd_bus
            sol(pm, nw, :load, load_id)[:qd_bus] = qd_bus
        end
    end
end


"""
	function constraint_mc_switch_thermal_limit(
		pm::AbstractNLExplicitNeutralIVRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		rating::Vector{<:Real}
	)

This method is not yet implemented for AbstractLPUBFModel.
If the limit is finite, a warning is thrown.
"""
function constraint_mc_switch_thermal_limit(pm::AbstractLPUBFModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{<:Real})
    if any(rating.<=Inf)
        @warn "Encountered a finite switch thermal limit; these are not yet implemented for AbstractLPUBFModel."
    end
end
