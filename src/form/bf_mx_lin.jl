import LinearAlgebra: diagm


"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::AbstractLPUBFModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
end


""
function variable_mc_branch_current(pm::AbstractLPUBFModel; kwargs...)
end


""
function variable_mc_voltage_prod_hermitian(pm::LPdiagUBFModel; n_cond::Int=3, nw::Int=pm.cnw, bounded = true)
    variable_mc_voltage_magnitude_sqr(pm, nw=nw)
end


""
function variable_mc_branch_flow(pm::LPfullUBFModel; n_cond::Int=3, nw::Int=pm.cnw, bounded::Bool=true,  report::Bool=true)
    @assert(n_cond == 3)
    variable_mc_branch_flow_active(pm, nw=nw, bounded=bounded)
    variable_mc_branch_flow_reactive(pm, nw=nw, bounded=bounded)

    #Store dictionary with matrix variables by arc
    p_mat_dict = Dict{Tuple{Int64,Int64,Int64}, Any}()
    q_mat_dict = Dict{Tuple{Int64,Int64,Int64}, Any}()

    for i in _PMs.ref(pm, nw, :arcs)
        l = i[1]
        f_bus = i[2]
        branch = _PMs.ref(pm, nw, :branch, l)
        p_d =  _PMs.var(pm, nw, :p)[i]
        q_d =  _PMs.var(pm, nw, :q)[i]

        alpha = exp(-im*2*pi/3)
        Gamma = [1 alpha^2 alpha; alpha 1 alpha^2; alpha^2 alpha 1]
        p_d_mx = [p_d[1] 0 0; 0 p_d[2] 0; 0 0 p_d[3]]
        q_d_mx = [q_d[1] 0 0; 0 q_d[2] 0; 0 0 q_d[3]]


        ps_mat = real(Gamma)*p_d_mx - imag(Gamma)*q_d_mx
        qs_mat = imag(Gamma)*p_d_mx + real(Gamma)*q_d_mx

        g_sh_fr = branch["g_fr"].values
        b_sh_fr = branch["b_fr"].values

        w_fr_re = _PMs.var(pm, nw, :Wr)[f_bus]
        w_fr_im = _PMs.var(pm, nw, :Wi)[f_bus]

        p_mat = ps_mat + (w_fr_re*(g_sh_fr)' + w_fr_im*(b_sh_fr)')
        q_mat = qs_mat + (w_fr_im*(g_sh_fr)' - w_fr_re*(b_sh_fr)')

        p_mat_dict[i] = p_mat
        q_mat_dict[i] = q_mat
    end
    _PMs.var(pm, nw)[:P] = p_mat_dict
    _PMs.var(pm, nw)[:Q] = q_mat_dict
end


""
function variable_mc_branch_flow(pm::LPdiagUBFModel; n_cond::Int=3, nw::Int=pm.cnw, bounded::Bool=true,  report::Bool=true)
    @assert(n_cond == 3)
    variable_mc_branch_flow_active(pm, nw=nw, bounded=bounded)
    variable_mc_branch_flow_reactive(pm, nw=nw, bounded=bounded)
end


"Defines branch flow model power flow equations"
function constraint_mc_flow_losses(pm::LPfullUBFModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to, tm)
    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    p_to = _PMs.var(pm, n, :P)[t_idx]
    q_to = _PMs.var(pm, n, :Q)[t_idx]

    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_to_re = _PMs.var(pm, n, :Wr)[t_bus]

    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]
    w_to_im = _PMs.var(pm, n, :Wi)[t_bus]

    JuMP.@constraint(pm.model, diag(p_fr) + diag(p_to) .==  diag(w_fr_re*(g_sh_fr)' + w_fr_im*(b_sh_fr)' +  w_to_re*(g_sh_to)'  + w_to_im*(b_sh_to)'))
    JuMP.@constraint(pm.model, diag(q_fr) + diag(q_to) .==  diag(w_fr_im*(g_sh_fr)' - w_fr_re*(b_sh_fr)' +  w_to_im*(g_sh_to)'  - w_to_re*(b_sh_to)'))
end


"Defines branch flow model power flow equations"
function constraint_mc_flow_losses(pm::LPdiagUBFModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to, tm)
    p_fr = _PMs.var(pm, n, :p)[f_idx]
    q_fr = _PMs.var(pm, n, :q)[f_idx]

    p_to = _PMs.var(pm, n, :p)[t_idx]
    q_to = _PMs.var(pm, n, :q)[t_idx]

    w_fr = _PMs.var(pm, n, :w)[f_bus]
    w_to = _PMs.var(pm, n, :w)[t_bus]

    JuMP.@constraint(pm.model, p_fr + p_to .== diag( g_sh_fr).*w_fr + diag( g_sh_to).*w_to)
    JuMP.@constraint(pm.model, q_fr + q_to .== diag(-b_sh_fr).*w_fr + diag(-b_sh_to).*w_to)
end


"Defines voltage drop over a branch, linking from and to side voltage"
function constraint_mc_model_voltage_magnitude_difference(pm::LPfullUBFModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]

    w_to_re = _PMs.var(pm, n, :Wr)[t_bus]
    w_to_im = _PMs.var(pm, n, :Wi)[t_bus]

    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    p_s_fr = p_fr - (w_fr_re*(g_sh_fr)' + w_fr_im*(b_sh_fr)')
    q_s_fr = q_fr - (w_fr_im*(g_sh_fr)' - w_fr_re*(b_sh_fr)')

        #KVL over the line:
    JuMP.@constraint(pm.model, diag(w_to_re)        .==        diag(w_fr_re   - p_s_fr  *r' - q_s_fr*x'        - r*p_s_fr'    - x*q_s_fr'))
end


"Defines voltage drop over a branch, linking from and to side voltage"
function constraint_mc_model_voltage_magnitude_difference(pm::LPdiagUBFModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
    w_fr = _PMs.var(pm, n, :w)[f_bus]
    w_to = _PMs.var(pm, n, :w)[t_bus]

    p_fr = _PMs.var(pm, n, :p)[f_idx]
    q_fr = _PMs.var(pm, n, :q)[f_idx]

    p_s_fr = p_fr - diag(g_sh_fr).*w_fr
    q_s_fr = q_fr + diag(b_sh_fr).*w_to

    alpha = exp(-im*2*pi/3)
    Gamma = [1 alpha^2 alpha; alpha 1 alpha^2; alpha^2 alpha 1]

    MP = 2*(real(Gamma).*r + imag(Gamma).*x)
    MQ = 2*(real(Gamma).*x - imag(Gamma).*r)

    JuMP.@constraint(pm.model, w_to .== w_fr - MP*p_s_fr - MQ*q_s_fr)

end


"balanced three-phase phasor"
function constraint_mc_theta_ref(pm::LPdiagUBFModel, n::Int, i)
    ncnds = length(_PMs.conductor_ids(pm))
    @assert(ncnds >= 2)

    w = _PMs.var(pm, n, :w)[i]
    JuMP.@constraint(pm.model, w[2:ncnds]   .== w[1])
end


""
function constraint_mc_power_balance(pm::LPdiagUBFModel, nw::Int, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    w = _PMs.var(pm, nw, :w, i)

    p = get(_PMs.var(pm, nw), :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q = get(_PMs.var(pm, nw), :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")

    psw  = get(_PMs.var(pm, nw),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(_PMs.var(pm, nw),  :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(_PMs.var(pm, nw),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(_PMs.var(pm, nw),   :qt, Dict()); _PMs._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")

    pg = get(_PMs.var(pm, nw), :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg = get(_PMs.var(pm, nw), :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PMs.var(pm, nw),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PMs.var(pm, nw),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")

    cstr_p = []
    cstr_q = []

    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        .==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        - sum(gs for gs in values(bus_gs))*w
    )

    cstr_q = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(qt[a_trans] for a_trans in bus_arcs_trans)
        .==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd for qd in values(bus_qd))
        + sum(bs for bs in values(bus_bs))*w
    )

    if _PMs.report_duals(pm)
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMs.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end
