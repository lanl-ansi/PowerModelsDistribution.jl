import LinearAlgebra: diagm


"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::AbstractLPUBFModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
end


""
function variable_mc_branch_current(pm::AbstractLPUBFModel; kwargs...)
end


""
function variable_mc_bus_voltage(pm::LPUBFDiagModel; n_cond::Int=3, nw::Int=pm.cnw, bounded::Bool=true)
    variable_mc_bus_voltage_magnitude_sqr(pm, nw=nw)
end


""
function variable_mc_branch_power(pm::LPUBFDiagModel; n_cond::Int=3, nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    @assert(n_cond == 3)
    variable_mc_branch_power_real(pm, nw=nw, bounded=bounded)
    variable_mc_branch_power_imaginary(pm, nw=nw, bounded=bounded)
end


"Defines branch flow model power flow equations"
function constraint_mc_power_losses(pm::LPUBFDiagModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to, tm)
  p_fr = var(pm, n, :p)[f_idx]
  q_fr = var(pm, n, :q)[f_idx]

  p_to = var(pm, n, :p)[t_idx]
  q_to = var(pm, n, :q)[t_idx]

  w_fr = var(pm, n, :w)[f_bus]
  w_to = var(pm, n, :w)[t_bus]


  fb = ref(pm, n, :bus, f_bus)
  tb = ref(pm, n, :bus, t_bus)
  vmax1=fb["vmax"]
  vmax2=tb["vmax"]

  cnds = conductor_ids(pm; nw=n)
  ncnds = length(cnds)

  for c in 1:ncnds
    if((vmax1[c]!=0) && (vmax2[c]!=0))
      JuMP.@constraint(pm.model, p_fr[c] + p_to[c] == g_sh_fr[c,c]*w_fr[c] +  g_sh_to[c,c]*w_to[c])
      JuMP.@constraint(pm.model, q_fr[c] + q_to[c] == -b_sh_fr[c,c]*w_fr[c] + -b_sh_to[c,c]*w_to[c])
    else
      JuMP.@constraint(pm.model, p_fr[c] + p_to[c] == 0)
      JuMP.@constraint(pm.model, q_fr[c] + q_to[c] == 0)
    end
  end
end


"Defines voltage drop over a branch, linking from and to side voltage"
function constraint_mc_model_voltage_magnitude_difference(pm::LPUBFDiagModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
    w_fr = var(pm, n, :w)[f_bus]
    w_to = var(pm, n, :w)[t_bus]

    p_fr = var(pm, n, :p)[f_idx]
    q_fr = var(pm, n, :q)[f_idx]

    p_s_fr = p_fr- diag(g_sh_fr).*w_fr
    q_s_fr = q_fr+ diag(b_sh_fr).*w_fr

    fb = ref(pm, n, :bus, f_bus)
    tb = ref(pm, n, :bus, t_bus)

    vmax1=fb["vmax"]
    vmax2=tb["vmax"]

    cnds = conductor_ids(pm; nw=n)
    ncnds = length(cnds)

    alpha = exp(-im*2*pi/3)
    Gamma = [1 alpha^2 alpha; alpha 1 alpha^2; alpha^2 alpha 1]

    MP = 2*(real(Gamma).*r + imag(Gamma).*x)
    MQ = 2*(real(Gamma).*x - imag(Gamma).*r)

    for c in 1:ncnds
      if((vmax1[c]!=0) && (vmax2[c]!=0))
        JuMP.@constraint(pm.model, w_to[c]== w_fr[c] - (MP[c,1]*p_s_fr[1]+MP[c,2]*p_s_fr[2]+MP[c,3]*p_s_fr[3]) - (MQ[c,1]*q_s_fr[1]+MQ[c,2]*q_s_fr[2]+MQ[c,3]*q_s_fr[3]))
      else
        JuMP.@constraint(pm.model, p_s_fr[c]==0)
        JuMP.@constraint(pm.model, q_s_fr[c]==0)
      end
   end
end


"balanced three-phase phasor"
function constraint_mc_theta_ref(pm::LPUBFDiagModel, n::Int, i::Int, va_ref)
    ncnds = length(conductor_ids(pm))
    @assert(ncnds >= 2)
    w = var(pm, n, :w)[i]
    JuMP.@constraint(pm.model, w[2:ncnds]   .== w[1])
end


""
function constraint_mc_load_power_balance(pm::LPUBFDiagModel, nw::Int, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_converter, bus_loads, bus_gs, bus_bs)
    w = var(pm, nw, :w, i)
    bus = ref(pm, nw, :bus, i)
    vmax=bus["vmax"]

    p = get(var(pm, nw), :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q = get(var(pm, nw), :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")

    psw  = get(var(pm, nw),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, nw),  :qsw, Dict()); _PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(var(pm, nw),   :pt, Dict()); _PM._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(var(pm, nw),   :qt, Dict()); _PM._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")

    pd = get(var(pm, nw), :pd, Dict()); _PM._check_var_keys(pd, bus_loads, "active power", "load")
    qd = get(var(pm, nw), :qd, Dict()); _PM._check_var_keys(qd, bus_loads, "reactive power", "load")
    pg = get(var(pm, nw), :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg = get(var(pm, nw), :qg, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, nw),   :ps, Dict()); _PM._check_var_keys(ps, bus_converter, "active power", "converter")
    qs   = get(var(pm, nw),   :qs, Dict()); _PM._check_var_keys(qs, bus_converter, "reactive power", "converter")

    cstr_p = []
    cstr_q = []

    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    for c in 1:ncnds
        if(vmax[c]!=0)
            cp = JuMP.@constraint(pm.model,
            sum((p[a])[c] for a in bus_arcs)
            + sum((psw[a_sw])[c] for a_sw in bus_arcs_sw)
            + sum((pt[a_trans])[c] for a_trans in bus_arcs_trans)
            ==
            sum((pg[g])[c] for g in bus_gens)
            - sum((ps[s])[c] for s in bus_converter)
            - sum((pd[d])[c] for d in bus_loads)
            - sum(gs[c]*w[c] for gs in values(bus_gs)))
            push!(cstr_p, cp)
        end
    end

    for c in 1:ncnds
        if(vmax[c]!=0)
            cq = JuMP.@constraint(pm.model,
            sum((q[a])[c] for a in bus_arcs)
            + sum((qsw[a_sw])[c] for a_sw in bus_arcs_sw)
            + sum((qt[a_trans])[c] for a_trans in bus_arcs_trans)
            ==
            sum((qg[g])[c] for g in bus_gens)
            - sum((qs[s])[c] for s in bus_converter)
            - sum((qd[d])[c] for d in bus_loads)
            + sum(bs[c]*w[c] for bs in values(bus_bs)))
            push!(cstr_q, cq)
        end
   end

    con(pm, nw, :lam_kcl_r)[i] = isa(cstr_p, Array) ? cstr_p : [cstr_p]
    con(pm, nw, :lam_kcl_i)[i] = isa(cstr_q, Array) ? cstr_q : [cstr_q]

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end
