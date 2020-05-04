"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SDPUBFModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
    p_fr = var(pm, n, :P)[f_idx]
    q_fr = var(pm, n, :Q)[f_idx]

    w_fr_re = var(pm, n, :Wr)[f_bus]
    w_fr_im = var(pm, n, :Wi)[f_bus]

    ccm_re =  var(pm, n, :CCr)[i]
    ccm_im =  var(pm, n, :CCi)[i]

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


"Add explicit PSD-ness of W for nodes where it is not implied"
function constraint_mc_voltage_psd(pm::SDPUBFModel, n::Int, i)
    Wr = var(pm, n, :Wr)[i]
    Wi = var(pm, n, :Wi)[i]

    JuMP.@constraint(pm.model,
    [
    Wr  -Wi;
    Wi   Wr
    ] in JuMP.PSDCone())
end
