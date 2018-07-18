""
function variable_branch_flow(pm::GenericPowerModel{T}; n_cond::Int=3, nw::Int=pm.cnw, bounded = true) where T<:LPUBFForm
    n_diag_el = n_cond
    n_lower_triangle_el = Int((n_cond^2 - n_cond)/2)
    assert(n_cond==3)#TODO generalize this variable defintion

    for i in 1:n_diag_el
        PMs.variable_active_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
        PMs.variable_reactive_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
    end

    #Store dictionary with matrix variables by arc
    p_mat_dict = Dict{Tuple{Int64,Int64,Int64}, Any}()
    q_mat_dict = Dict{Tuple{Int64,Int64,Int64}, Any}()

    for i in ref(pm, nw, :arcs)
        l = i[1]
        f_bus = i[2]
        branch = ref(pm, nw, :branch, l)
        p_d =  [var(pm, nw, c, :p,    i) for c in 1:n_diag_el]
        q_d =  [var(pm, nw, c, :q,    i) for c in 1:n_diag_el]

        alpha = exp(-im*2*pi/3)
        Gamma = [1 alpha^2 alpha; alpha 1 alpha^2; alpha^2 alpha 1]

        ps_mat = real(Gamma)*diagm(p_d) - imag(Gamma)*diagm(q_d)
        qs_mat = imag(Gamma)*diagm(p_d) + real(Gamma)*diagm(q_d)

        g_fr = diagm(branch["g_fr"].values) #TODO remove .values
        b_fr = diagm(branch["b_fr"].values) #TODO remove .values

        w_fr_re = var(pm, nw, :w_re)[f_bus]

        p_mat = ps_mat + g_fr*w_fr_re
        q_mat = qs_mat - b_fr*w_fr_re

        p_mat_dict[i] = p_mat
        q_mat_dict[i] = q_mat
    end
    var(pm, nw)[:p_mat] = p_mat_dict
    var(pm, nw)[:q_mat] = q_mat_dict
end

"""
Defines branch flow model power flow equations
"""
function constraint_tp_flow_losses_mat(pm::GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to) where T <: LPUBFForm
    p_fr = var(pm, n, :p_mat)[f_idx]
    q_fr = var(pm, n, :q_mat)[f_idx]

    p_to = var(pm, n, :p_mat)[t_idx]
    q_to = var(pm, n, :q_mat)[t_idx]

    w_fr_re = var(pm, n, :w_re)[f_bus]
    w_to_re = var(pm, n, :w_re)[t_bus]

    @constraint(pm.model, diag(p_fr) + diag(p_to) .==  diag( g_sh_fr*(w_fr_re) +  g_sh_to*w_to_re))
    @constraint(pm.model, diag(q_fr) + diag(q_to) .==  diag(-b_sh_fr*(w_fr_re) + -b_sh_to*w_to_re))
end

"""
Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude
"""
function constraint_tp_branch_current_mat(pm::GenericPowerModel{T}, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr) where T <: LPUBFForm
end

function constraint_tp_valid_inquality_branch(pm::GenericPowerModel{T}, n::Int, i, f_idx, t_idx) where T <: LPUBFForm
    #not needed, line is lossless already
end


function constraint_tp_valid_inquality_bus(pm::GenericPowerModel{T}, n::Int, i, vmax) where T <: LPUBFForm
end

"""
Defines voltage drop over a branch, linking from and to side voltage
"""
function constraint_tp_voltage_magnitude_difference_mat(pm::GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm) where T <: LPfullUBFForm
    w_fr_re = var(pm, n, :w_re)[f_bus]
    w_fr_im = var(pm, n, :w_im)[f_bus]

    w_to_re = var(pm, n, :w_re)[t_bus]
    w_to_im = var(pm, n, :w_im)[t_bus]

    p_fr = var(pm, n, :p_mat)[f_idx]
    q_fr = var(pm, n, :q_mat)[f_idx]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

        #KVL over the line:
    @constraint(pm.model, diag(w_to_re) .== diag(w_fr_re   - p_s_fr  *r' - q_s_fr*x'        - r*p_s_fr'    - x*q_s_fr'))
    @constraint(pm.model, mat2utrivec(w_to_re) .== mat2utrivec(w_fr_re   - p_s_fr  *r' - q_s_fr*x'        - r*p_s_fr'    - x*q_s_fr'))
    @constraint(pm.model, mat2utrivec(w_to_im) .== mat2utrivec(w_fr_im   - q_s_fr  *r' + p_s_fr*x'        - x*p_s_fr'    + r*q_s_fr'))
end

"""
Defines voltage drop over a branch, linking from and to side voltage
"""
function constraint_tp_voltage_magnitude_difference_mat(pm::GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm) where T <: LPdiagUBFForm
    w_fr_re = var(pm, n, :w_re)[f_bus]

    w_to_re = var(pm, n, :w_re)[t_bus]

    p_fr = var(pm, n, :p_mat)[f_idx]
    q_fr = var(pm, n, :q_mat)[f_idx]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

        #KVL over the line:
    @constraint(pm.model, diag(w_to_re) .== diag(w_fr_re   - p_s_fr  *r' - q_s_fr*x'        - r*p_s_fr'    - x*q_s_fr'))
end
