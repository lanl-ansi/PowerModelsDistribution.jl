"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SOCUBFModels, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]

    ccm_re =  _PMs.var(pm, n, :CCr)[i]
    ccm_im =  _PMs.var(pm, n, :CCi)[i]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

    mat_real = [
    w_fr_re     p_s_fr  ;
    p_s_fr'    ccm_re  ;
    ]

    mat_imag = [
    w_fr_im     q_s_fr  ;
    -q_s_fr'    ccm_im  ;
    ]

    relaxation_psd_to_soc(pm.model, mat_real, mat_imag)

    # code below useful for debugging: valid inequality equired to make the SOC-NLP formulation more accurate
    (l,i,j) = f_idx
    t_idx = (l,j,i)
    p_to = _PMs.var(pm, n, :P)[t_idx]
    # total losses are positive when g_fr, g_to and r are positive
    # not guaranteed for individual phases though when matrix obtained through Kron's reduction
    # JuMP.@constraint(pm.model, tr(p_fr) + tr(p_to) >= 0)
    # JuMP.@constraint(pm.model, sum(diag(p_fr)) + sum(diag(p_to)) .>= 0)
end


"Add explicit PSD-ness of W for nodes where it is not implied"
function constraint_mc_voltage_psd(pm::SOCUBFModels, n::Int, i)
    Wr = _PMs.var(pm, n, :Wr)[i]
    Wi = _PMs.var(pm, n, :Wi)[i]

    relaxation_psd_to_soc(pm.model, Wr, Wi)
end



"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SOCConicUBFModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]

    ccm_re =  _PMs.var(pm, n, :CCr)[i]
    ccm_im =  _PMs.var(pm, n, :CCi)[i]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

    mat_real = [
    w_fr_re     p_s_fr  ;
    p_s_fr'    ccm_re  ;
    ]

    mat_imag = [
    w_fr_im     q_s_fr  ;
    -q_s_fr'    ccm_im  ;
    ]

    relaxation_psd_to_soc_conic(pm.model, mat_real, mat_imag)

    (l,i,j) = f_idx
    t_idx = (l,j,i)
    p_to = _PMs.var(pm, n, :P)[t_idx]
    # JuMP.@constraint(pm.model, sum(diag(p_fr)) + sum(diag(p_to)) .>= 0)

end

"Add explicit PSD-ness of W for nodes where it is not implied"
function constraint_mc_voltage_psd(pm::SOCConicUBFModel, n::Int, i)
    Wr = _PMs.var(pm, n, :Wr)[i]
    Wi = _PMs.var(pm, n, :Wi)[i]

    relaxation_psd_to_soc_conic(pm.model, Wr, Wi)
end


"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::QPUBFPowerModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]

    ccm_re =  _PMs.var(pm, n, :CCr)[i]
    ccm_im =  _PMs.var(pm, n, :CCi)[i]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

    mat_real = [
    w_fr_re     p_s_fr  ;
    p_s_fr'    ccm_re  ;
    ]

    mat_imag = [
    w_fr_im     q_s_fr  ;
    -q_s_fr'    ccm_im  ;
    ]

    reformulation_psd_rank1_to_quadratic(pm.model, mat_real, mat_imag)
end


"do nothing, not needed in this variable space"
function constraint_mc_voltage_psd(pm::QPUBFPowerModel, n::Int, i)
end


"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SOCNLPUBFKimKojimaPowerModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]

    ccm_re =  _PMs.var(pm, n, :CCr)[i]
    ccm_im =  _PMs.var(pm, n, :CCi)[i]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

    mat_real = [
    w_fr_re     p_s_fr  ;
    p_s_fr'    ccm_re  ;
    ]

    mat_imag = [
    w_fr_im     q_s_fr  ;
    -q_s_fr'    ccm_im  ;
    ]

    relaxation_psd_to_soc(pm.model, mat_real, mat_imag)
    relaxation_psd_to_soc_complex_kim_kojima_3x3(pm.model, w_fr_re, w_fr_im)
    relaxation_psd_to_soc_complex_kim_kojima_3x3(pm.model, ccm_re, ccm_im)
end


"Add explicit PSD-ness of W for nodes where it is not implied"
function constraint_mc_voltage_psd(pm::SOCNLPUBFKimKojimaPowerModel, n::Int, i)
    Wr = _PMs.var(pm, n, :Wr)[i]
    Wi = _PMs.var(pm, n, :Wi)[i]

    relaxation_psd_to_soc(pm.model, Wr, Wi)
    relaxation_psd_to_soc_complex_kim_kojima_3x3(pm.model, Wr, Wi)
end

"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SOCConicUBFKimKojimaPowerModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]

    ccm_re =  _PMs.var(pm, n, :CCr)[i]
    ccm_im =  _PMs.var(pm, n, :CCi)[i]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

    mat_real = [
    w_fr_re     p_s_fr  ;
    p_s_fr'    ccm_re  ;
    ]

    mat_imag = [
    w_fr_im     q_s_fr  ;
    -q_s_fr'    ccm_im  ;
    ]

    relaxation_psd_to_soc_conic(pm.model, mat_real, mat_imag, complex=true)
    relaxation_psd_to_soc_complex_kim_kojima_3x3_conic(pm.model, w_fr_re, w_fr_im)
    relaxation_psd_to_soc_complex_kim_kojima_3x3_conic(pm.model, ccm_re, ccm_im)
end


"Add explicit PSD-ness of W for nodes where it is not implied"
function constraint_mc_voltage_psd(pm::SOCConicUBFKimKojimaPowerModel, n::Int, i)
    Wr = _PMs.var(pm, n, :Wr)[i]
    Wi = _PMs.var(pm, n, :Wi)[i]

    relaxation_psd_to_soc_conic(pm.model, Wr, Wi)
    relaxation_psd_to_soc_complex_kim_kojima_3x3_conic(pm.model, Wr, Wi)
end

"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SOCNLPUBFKCLMXKimKojimaPowerModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]

    ccm_re =  _PMs.var(pm, n, :CCr)[i]
    ccm_im =  _PMs.var(pm, n, :CCi)[i]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

    mat_real = [
    w_fr_re     p_s_fr  ;
    p_s_fr'    ccm_re  ;
    ]

    mat_imag = [
    w_fr_im     q_s_fr  ;
    -q_s_fr'    ccm_im  ;
    ]

    relaxation_psd_to_soc(pm.model, mat_real, mat_imag, complex=true)
    relaxation_psd_to_soc_complex_kim_kojima_3x3(pm.model, w_fr_re, w_fr_im)
    relaxation_psd_to_soc_complex_kim_kojima_3x3(pm.model, ccm_re, ccm_im)
end


"Add explicit PSD-ness of W for nodes where it is not implied"
function constraint_mc_voltage_psd(pm::SOCNLPUBFKCLMXKimKojimaPowerModel, n::Int, i)
    Wr = _PMs.var(pm, n, :Wr)[i]
    Wi = _PMs.var(pm, n, :Wi)[i]

    relaxation_psd_to_soc(pm.model, Wr, Wi)
    relaxation_psd_to_soc_complex_kim_kojima_3x3(pm.model, Wr, Wi)
end


"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SOCConicUBFKCLMXKimKojimaPowerModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]

    ccm_re =  _PMs.var(pm, n, :CCr)[i]
    ccm_im =  _PMs.var(pm, n, :CCi)[i]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

    mat_real = [
    w_fr_re     p_s_fr  ;
    p_s_fr'    ccm_re  ;
    ]

    mat_imag = [
    w_fr_im     q_s_fr  ;
    -q_s_fr'    ccm_im  ;
    ]

    relaxation_psd_to_soc_conic(pm.model, mat_real, mat_imag, complex=true)
    relaxation_psd_to_soc_complex_kim_kojima_3x3_conic(pm.model, w_fr_re, w_fr_im)
    relaxation_psd_to_soc_complex_kim_kojima_3x3_conic(pm.model, ccm_re, ccm_im)
end


"Add explicit PSD-ness of W for nodes where it is not implied"
function constraint_mc_voltage_psd(pm::SOCConicUBFKCLMXKimKojimaPowerModel, n::Int, i)
    Wr = _PMs.var(pm, n, :Wr)[i]
    Wi = _PMs.var(pm, n, :Wi)[i]

    relaxation_psd_to_soc_conic(pm.model, Wr, Wi)
    relaxation_psd_to_soc_complex_kim_kojima_3x3_conic(pm.model, Wr, Wi)
end

"""
Link the current and power withdrawn by a generator at the bus through a PSD
constraint. The rank-1 constraint is dropped in this formulation.
"""
function constraint_mc_generation(pm::SOCNLPUBFKCLMXKimKojimaPowerModel, gen_id::Int; nw::Int=pm.cnw)
    Pg = _PMs.var(pm, nw, :Pg, gen_id)
    Qg = _PMs.var(pm, nw, :Qg, gen_id)
    bus_id = _PMs.ref(pm, nw, :gen, gen_id)["gen_bus"]
    Wr = _PMs.var(pm, nw, :Wr, bus_id)
    Wi = _PMs.var(pm, nw, :Wi, bus_id)
    CCgr = _PMs.var(pm, nw, :CCgr, gen_id)
    CCgi = _PMs.var(pm, nw, :CCgi, gen_id)
    # constraint_SWL_psd(pm.model, Pg, Qg, Wr, Wi, CCgr, CCgi)

    mat_real = [
    Wr     Pg  ;
    Pg'    CCgr
    ]

    mat_imag = [
    Wi     Qg  ;
    -Qg'    CCgi
    ]
    relaxation_psd_to_soc(pm.model, mat_real, mat_imag, complex=true)
    relaxation_psd_to_soc_complex_kim_kojima_3x3(pm.model, CCgr, CCgi)
end

"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::QPUBFKCLMXPowerModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]

    ccm_re =  _PMs.var(pm, n, :CCr)[i]
    ccm_im =  _PMs.var(pm, n, :CCi)[i]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

    mat_real = [
    w_fr_re     p_s_fr  ;
    p_s_fr'    ccm_re  ;
    ]

    mat_imag = [
    w_fr_im     q_s_fr  ;
    -q_s_fr'    ccm_im  ;
    ]

    reformulation_psd_rank1_to_quadratic(pm.model, mat_real, mat_imag)
end

"Add explicit PSD-ness of W for nodes where it is not implied"
function constraint_mc_voltage_psd(pm::QPUBFKCLMXPowerModel, n::Int, i)

end

"""
Link the current and power withdrawn by a generator at the bus through a PSD
constraint. The rank-1 constraint is dropped in this formulation.
"""
function constraint_mc_generation(pm::QPUBFKCLMXPowerModel, gen_id::Int; nw::Int=pm.cnw)
    Pg = _PMs.var(pm, nw, :Pg, gen_id)
    Qg = _PMs.var(pm, nw, :Qg, gen_id)
    bus_id = _PMs.ref(pm, nw, :gen, gen_id)["gen_bus"]
    Wr = _PMs.var(pm, nw, :Wr, bus_id)
    Wi = _PMs.var(pm, nw, :Wi, bus_id)
    CCgr = _PMs.var(pm, nw, :CCgr, gen_id)
    CCgi = _PMs.var(pm, nw, :CCgi, gen_id)
    # constraint_SWL_psd(pm.model, Pg, Qg, Wr, Wi, CCgr, CCgi)

    mat_real = [
    Wr     Pg  ;
    Pg'    CCgr
    ]

    mat_imag = [
    Wi     Qg  ;
    -Qg'    CCgi
    ]
    reformulation_psd_rank1_to_quadratic(pm.model, mat_real, mat_imag)
end
