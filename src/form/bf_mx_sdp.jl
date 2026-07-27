"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SDPUBFModel, nw::Int, i::Int, f_bus::Int, f_idx::Tuple{Int,Int,Int}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real})
    branch = ref(pm, nw, :branch, f_idx[1])
    f_connections = branch["f_connections"]
    t_connections = branch["t_connections"]
    bus = ref(pm, nw, :bus, f_bus)
    terminals = bus["terminals"]

    p_fr = var(pm, nw, :P, f_idx)
    q_fr = var(pm, nw, :Q, f_idx)

    w_fr_re = var(pm, nw, :Wr, f_bus)[[findfirst(isequal(fc), terminals) for fc in f_connections],[findfirst(isequal(tc), terminals) for tc in t_connections]]
    w_fr_im = var(pm, nw, :Wi, f_bus)[[findfirst(isequal(fc), terminals) for fc in f_connections],[findfirst(isequal(tc), terminals) for tc in t_connections]]

    ccm_re =  var(pm, nw, :CCr, i)
    ccm_im =  var(pm, nw, :CCi, i)

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
function constraint_mc_voltage_psd(pm::SDPUBFModel, nw::Int, i::Int)
    Wr = var(pm, nw, :Wr, i)
    Wi = var(pm, nw, :Wi, i)

    JuMP.@constraint(pm.model,
    [
    Wr  -Wi;
    Wi   Wr
    ] in JuMP.PSDCone())
end
