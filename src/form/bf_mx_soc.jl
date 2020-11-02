"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SOCUBFModels, nw::Int, i::Int, f_bus::Int, f_idx::Tuple{Int,Int,Int}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real})
    branch = ref(pm, nw, :branch, f_idx[1])
    f_connections = branch["f_connections"]
    t_connections = branch["t_connections"]
    bus = ref(pm, nw, :bus, f_bus)
    terminals = bus["terminals"]

    p_fr = var(pm, nw, :P)[f_idx]
    q_fr = var(pm, nw, :Q)[f_idx]

    w_fr_re = var(pm, nw, :Wr, f_bus)[[findfirst(isequal(fc), terminals) for fc in f_connections],[findfirst(isequal(tc), terminals) for tc in t_connections]]
    w_fr_im = var(pm, nw, :Wi, f_bus)[[findfirst(isequal(fc), terminals) for fc in f_connections],[findfirst(isequal(tc), terminals) for tc in t_connections]]

    ccm_re =  var(pm, nw, :CCr)[i]
    ccm_im =  var(pm, nw, :CCi)[i]

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

    # code below useful for debugging: valid inequality equired to make the SOC-NLP formulation more accurate
    # (l,i,j) = f_idx
    # t_idx = (l,j,i)
    # p_to = var(pm, n, :P)[t_idx]
    # total losses are positive when g_fr, g_to and r are positive
    # not guaranteed for individual phases though when matrix obtained through Kron's reduction
    # JuMP.@constraint(pm.model, tr(p_fr) + tr(p_to) >= 0)
end


"Add explicit PSD-ness of W for nodes where it is not implied"
function constraint_mc_voltage_psd(pm::SOCUBFModels, nw::Int, i::Int)
    Wr = var(pm, nw, :Wr)[i]
    Wi = var(pm, nw, :Wi)[i]

    relaxation_psd_to_soc(pm.model, Wr, Wi)
end


"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SOCConicUBFModel, nw::Int, i::Int, f_bus::Int, f_idx::Tuple{Int,Int,Int}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real})
    branch = ref(pm, nw, :branch, f_idx[1])
    f_connections = branch["f_connections"]
    t_connections = branch["t_connections"]
    bus = ref(pm, nw, :bus, f_bus)
    terminals = bus["terminals"]

    p_fr = var(pm, nw, :P)[f_idx]
    q_fr = var(pm, nw, :Q)[f_idx]

    w_fr_re = var(pm, nw, :Wr, f_bus)[[findfirst(isequal(fc), terminals) for fc in f_connections],[findfirst(isequal(tc), terminals) for tc in t_connections]]
    w_fr_im = var(pm, nw, :Wi, f_bus)[[findfirst(isequal(fc), terminals) for fc in f_connections],[findfirst(isequal(tc), terminals) for tc in t_connections]]

    ccm_re =  var(pm, nw, :CCr)[i]
    ccm_im =  var(pm, nw, :CCi)[i]

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
end


"Add explicit PSD-ness of W for nodes where it is not implied"
function constraint_mc_voltage_psd(pm::SOCConicUBFModel, nw::Int, i::Int)
    Wr = var(pm, nw, :Wr)[i]
    Wi = var(pm, nw, :Wi)[i]

    relaxation_psd_to_soc_conic(pm.model, Wr, Wi)
end
