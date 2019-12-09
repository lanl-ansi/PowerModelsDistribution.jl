"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SDPUBFModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]

    ccm_re =  _PMs.var(pm, n, :CCr)[i]
    ccm_im =  _PMs.var(pm, n, :CCi)[i]

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

    JuMP.@constraint(pm.model,
    [
    mat_real  -mat_imag;
    mat_imag   mat_real
    ] in JuMP.PSDCone())

end
