import LinearAlgebra: diagm


"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::AbstractLPUBFModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
end


""
function variable_mc_branch_current(pm::AbstractLPUBFModel; kwargs...)
end


""
function variable_mc_bus_voltage(pm::LPUBFDiagModel; nw::Int=pm.cnw, bounded::Bool=true)
    variable_mc_bus_voltage_magnitude_sqr(pm, nw=nw)
end


""
function variable_mc_branch_power(pm::LPUBFDiagModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    variable_mc_branch_power_real(pm, nw=nw, bounded=bounded)
    variable_mc_branch_power_imaginary(pm, nw=nw, bounded=bounded)
end


"Defines branch flow model power flow equations"
function constraint_mc_power_losses(pm::LPUBFDiagModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, r::Array{<:Real}, x::Array{<:Real}, g_sh_fr::Array{<:Real}, g_sh_to::Array{<:Real}, b_sh_fr::Array{<:Real}, b_sh_to::Array{<:Real})
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
function constraint_mc_model_voltage_magnitude_difference(pm::LPUBFDiagModel, n::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, r::Array{<:Real}, x::Array{<:Real}, g_sh_fr::Array{<:Real}, b_sh_fr::Array{<:Real})
    f_connections = ref(pm, n, :branch, i)["f_connections"]
    t_connections = ref(pm, n, :branch, i)["t_connections"]

    w_fr = var(pm, n, :w)[f_bus]
    w_to = var(pm, n, :w)[t_bus]

    p_fr = var(pm, n, :p)[f_idx]
    q_fr = var(pm, n, :q)[f_idx]

    dg = size(f_connections)[1] == 1 ? :diagm : :diag 
    p_s_fr = [p_fr[fc]- eval.(dg)(g_sh_fr)[idx].*w_fr[fc] for (idx,fc) in enumerate(f_connections)]
    q_s_fr = [q_fr[fc]+ eval.(dg)(b_sh_fr)[idx].*w_fr[fc] for (idx,fc) in enumerate(f_connections)]
    
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


""
function constraint_mc_power_balance(pm::LPUBFDiagModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    w = var(pm, nw, :w, i)
    p   = get(var(pm, nw),      :p,   Dict()); _PM._check_var_keys(p,   bus_arcs, "active power", "branch")
    q   = get(var(pm, nw),      :q,   Dict()); _PM._check_var_keys(q,   bus_arcs, "reactive power", "branch")
    psw = get(var(pm, nw),    :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw = get(var(pm, nw),    :qsw, Dict()); _PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt  = get(var(pm, nw),     :pt,  Dict()); _PM._check_var_keys(pt,  bus_arcs_trans, "active power", "transformer")
    qt  = get(var(pm, nw),     :qt,  Dict()); _PM._check_var_keys(qt,  bus_arcs_trans, "reactive power", "transformer")
    pg  = get(var(pm, nw),     :pg,  Dict()); _PM._check_var_keys(pg,  bus_gens, "active power", "generator")
    qg  = get(var(pm, nw),     :qg,  Dict()); _PM._check_var_keys(qg,  bus_gens, "reactive power", "generator")
    ps  = get(var(pm, nw),     :ps,  Dict()); _PM._check_var_keys(ps,  bus_storage, "active power", "storage")
    qs  = get(var(pm, nw),     :qs,  Dict()); _PM._check_var_keys(qs,  bus_storage, "reactive power", "storage")
    pd  = get(var(pm, nw), :pd_bus,  Dict()); _PM._check_var_keys(pd,  bus_loads, "active power", "load")
    qd  = get(var(pm, nw), :qd_bus,  Dict()); _PM._check_var_keys(qd,  bus_loads, "reactive power", "load")


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


"Neglects the active and reactive loss terms associated with the squared current magnitude."
function constraint_mc_storage_losses(pm::AbstractUBFModels, i::Int; nw::Int=pm.cnw, kwargs...)
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
