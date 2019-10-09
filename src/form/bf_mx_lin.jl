import LinearAlgebra: diagm


"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::AbstractLPUBFModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
end


""
function variable_mc_branch_current(pm::AbstractLPUBFModel; kwargs...)
end


""
function variable_mc_voltage_prod_hermitian(pm::LPdiagUBFModel; n_cond::Int=3, nw::Int=pm.cnw, bounded = true)
    n_diag_el = n_cond
    n_lower_triangle_el = Int((n_cond^2 - n_cond)/2)
    for c in 1:n_diag_el
        PowerModels.variable_voltage_magnitude_sqr(pm, nw=nw, cnd=c)
    end
    #Store dictionary with matrix variables by bus
    w_re_dict = Dict{Int64, Any}()
    for i in _PMs.ids(pm, nw, :bus)
        w =  [_PMs.var(pm, nw, h, :w,  i) for h in 1:n_diag_el]
        w_re_dict[i] = w
    end
    _PMs.var(pm, nw)[:Wr] = w_re_dict
end


""
function variable_mc_branch_flow(pm::LPfullUBFModel; n_cond::Int=3, nw::Int=pm.cnw, bounded = true)
    n_diag_el = n_cond
    n_lower_triangle_el = Int((n_cond^2 - n_cond)/2)

    for i in 1:n_diag_el
        _PMs.variable_active_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
        _PMs.variable_reactive_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
    end

    #Store dictionary with matrix variables by arc
    p_mat_dict = Dict{Tuple{Int64,Int64,Int64}, Any}()
    q_mat_dict = Dict{Tuple{Int64,Int64,Int64}, Any}()

    for i in _PMs.ref(pm, nw, :arcs)
        l = i[1]
        f_bus = i[2]
        branch = _PMs.ref(pm, nw, :branch, l)
        p_d =  [_PMs.var(pm, nw, c, :p,    i) for c in 1:n_diag_el]
        q_d =  [_PMs.var(pm, nw, c, :q,    i) for c in 1:n_diag_el]

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
function variable_mc_branch_flow(pm::LPdiagUBFModel; n_cond::Int=3, nw::Int=pm.cnw, bounded = true)
    n_diag_el = n_cond

    for i in 1:n_diag_el
        _PMs.variable_active_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
        _PMs.variable_reactive_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
    end

    #Store dictionary with matrix variables by arc
    p_mat_dict = Dict{Tuple{Int64,Int64,Int64}, Any}()
    q_mat_dict = Dict{Tuple{Int64,Int64,Int64}, Any}()

    for i in _PMs.ref(pm, nw, :arcs)
        l = i[1]
        f_bus = i[2]
        branch = _PMs.ref(pm, nw, :branch, l)
        p_d =  [_PMs.var(pm, nw, c, :p,    i) for c in 1:n_diag_el]
        q_d =  [_PMs.var(pm, nw, c, :q,    i) for c in 1:n_diag_el]

        ps_mat = p_d
        qs_mat = q_d

        g_sh_fr = branch["g_fr"].values
        b_sh_fr = branch["b_fr"].values

        w_fr_re = _PMs.var(pm, nw, :Wr)[f_bus]
        w_fr_re_mx = [w_fr_re[1] 0 0; 0 w_fr_re[2] 0; 0 0 w_fr_re[3]]

        p_mat = ps_mat + diag(( w_fr_re_mx*(g_sh_fr)'))
        q_mat = qs_mat + diag((-w_fr_re_mx*(b_sh_fr)'))


        p_mat_dict[i] = p_mat
        q_mat_dict[i] = q_mat
    end
    _PMs.var(pm, nw)[:P] = p_mat_dict
    _PMs.var(pm, nw)[:Q] = q_mat_dict
end


"Defines branch flow model power flow equations"
function constraint_mc_flow_losses(pm::LPfullUBFModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to)
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
function constraint_mc_flow_losses(pm::LPdiagUBFModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to)
    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    p_to = _PMs.var(pm, n, :P)[t_idx]
    q_to = _PMs.var(pm, n, :Q)[t_idx]

    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_to_re = _PMs.var(pm, n, :Wr)[t_bus]

    JuMP.@constraint(pm.model, p_fr + p_to .== diag( g_sh_fr).*w_fr_re + diag( g_sh_to).*w_to_re)
    JuMP.@constraint(pm.model, q_fr + q_to .== diag(-b_sh_fr).*w_fr_re + diag(-b_sh_to).*w_to_re)
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
    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_to_re = _PMs.var(pm, n, :Wr)[t_bus]

    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]



    p_s_fr = p_fr - diag(g_sh_fr).*w_fr_re
    q_s_fr = q_fr + diag(b_sh_fr).*w_fr_re

    alpha = exp(-im*2*pi/3)
    Gamma = [1 alpha^2 alpha; alpha 1 alpha^2; alpha^2 alpha 1]

    MP = 2*(real(Gamma).*r + imag(Gamma).*x)
    MQ = 2*(real(Gamma).*x - imag(Gamma).*r)

    JuMP.@constraint(pm.model, w_to_re .== w_fr_re - MP*p_s_fr - MQ*q_s_fr)

end


"balanced three-phase phasor"
function constraint_mc_theta_ref(pm::LPdiagUBFModel, n::Int, i)
    nconductors = length(_PMs.conductor_ids(pm))

    w_re = _PMs.var(pm, n, :Wr)[i]
    JuMP.@constraint(pm.model, w_re[2:3]   .== w_re[1])
end
