"""
Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude
"""
function constraint_tp_branch_current(pm::PMs.GenericPowerModel{T}, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr) where T <: SDPUBFForm
    p_fr = PMs.var(pm, n, :P_mx)[f_idx]
    q_fr = PMs.var(pm, n, :Q_mx)[f_idx]

    w_fr_re = PMs.var(pm, n, :W_re)[f_bus]
    w_fr_im = PMs.var(pm, n, :W_im)[f_bus]

    ccm_re =  PMs.var(pm, n, :CC_re)[i]
    ccm_im =  PMs.var(pm, n, :CC_im)[i]

    p_s_fr = p_fr - (w_fr_re*(g_sh_fr)' + w_fr_im*(b_sh_fr)')
    q_s_fr = q_fr - (w_fr_im*(g_sh_fr)' - w_fr_re*(b_sh_fr)')


    mat_real = [
    w_fr_re     p_s_fr  ;
    p_s_fr'    ccm_re  ;
    ]

    mat_imag = [
    w_fr_im     q_s_fr  ;
    -q_s_fr'    ccm_im  ;
    ]

    JuMP.@SDconstraint(pm.model,
    [
    mat_real  -mat_imag;
    mat_imag   mat_real
    ] >=0)

end
