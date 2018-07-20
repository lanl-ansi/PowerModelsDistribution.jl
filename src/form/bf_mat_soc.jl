
"""
Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude
"""
function constraint_tp_branch_current(pm::GenericPowerModel{T}, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr) where T <: SOCUBFForm
    p_fr = var(pm, n, :p_mat)[f_idx]
    q_fr = var(pm, n, :q_mat)[f_idx]

    w_fr_re = var(pm, n, :w_re)[f_bus]
    w_fr_im = var(pm, n, :w_im)[f_bus]

    ccm_re =  var(pm, n, :ccm_re)[i]
    ccm_im =  var(pm, n, :ccm_im)[i]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

    mat = [
    w_fr_re    p_s_fr   -w_fr_im   -q_s_fr;
    p_s_fr'    ccm_re    q_s_fr'   -ccm_im;
    w_fr_im    q_s_fr    w_fr_re    p_s_fr;
    -q_s_fr'   ccm_im    p_s_fr'    ccm_re
    ]

    psd_to_soc(pm, mat)
    # psd_to_soc_diag(pm, p_s_fr, q_s_fr, w_fr_re, ccm_re)

    #TODO Required to make the SOC-NLP formulation more accurate
    (l,i,j) = f_idx
    t_idx = (l,j,i)
    p_to = var(pm, n, :p_mat)[t_idx]
    # total losses are positive when g_fr, g_to and r are positive
    # not guaranteed for individual phases though when matrix obtained through Kron's reduction
    # @constraint(pm.model, trace(p_fr) + trace(p_to) >= 0)
end
