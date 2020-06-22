# Kirchhoff's circuit laws as defined the current-voltage variable space.
# Even though the branch model is linear, the feasible set is non-convex
# in the context of constant-power loads or generators

""
function variable_mc_branch_current(pm::_PM.AbstractIVRModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    # store expressions in rectangular power variable space
    p = Dict()
    q = Dict()

    for (l,i,j) in ref(pm, nw, :arcs_from)
        vr_fr = var(pm, nw, :vr, i)
        vi_fr = var(pm, nw, :vi, i)
        cr_fr = var(pm, nw, :cr, (l,i,j))
        ci_fr = var(pm, nw, :ci, (l,i,j))

        vr_to = var(pm, nw, :vr, j)
        vi_to = var(pm, nw, :vi, j)
        cr_to = var(pm, nw, :cr, (l,j,i))
        ci_to = var(pm, nw, :ci, (l,j,i))
        p[(l,i,j)] = vr_fr.*cr_fr  + vi_fr.*ci_fr
        q[(l,i,j)] = vi_fr.*cr_fr  - vr_fr.*ci_fr
        p[(l,j,i)] = vr_to.*cr_to  + vi_to.*ci_to
        q[(l,j,i)] = vi_to.*cr_to  - vr_to.*ci_to
    end

    var(pm, nw)[:p] = p
    var(pm, nw)[:q] = q
    report && _IM.sol_component_value_edge(pm, nw, :branch, :pf, :pt, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), p)
    report && _IM.sol_component_value_edge(pm, nw, :branch, :qf, :qt, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), q)

    variable_mc_branch_current_series_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_branch_current_series_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end


""
function variable_mc_transformer_current(pm::_PM.AbstractIVRModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_transformer_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_transformer_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    # store expressions in rectangular power variable space
    p = Dict()
    q = Dict()

    f_cnd = [1, 2, 3] #TODO extend to this being passed in constraint template
    t_cnd = [1, 2, 3] #TODO extend to this being passed in constraint template

    for (l,i,j) in ref(pm, nw, :arcs_from_trans)
        vr_fr = var(pm, nw, :vr, i)[f_cnd]
        vi_fr = var(pm, nw, :vi, i)[f_cnd]
        cr_fr = var(pm, nw, :crt, (l,i,j))
        ci_fr = var(pm, nw, :cit, (l,i,j))

        vr_to = var(pm, nw, :vr, j)[t_cnd]
        vi_to = var(pm, nw, :vi, j)[t_cnd]
        cr_to = var(pm, nw, :crt, (l,j,i))
        ci_to = var(pm, nw, :cit, (l,j,i))
        p[(l,i,j)] = vr_fr.*cr_fr  + vi_fr.*ci_fr
        q[(l,i,j)] = vi_fr.*cr_fr  - vr_fr.*ci_fr
        p[(l,j,i)] = vr_to.*cr_to  + vi_to.*ci_to
        q[(l,j,i)] = vi_to.*cr_to  - vr_to.*ci_to
    end

    var(pm, nw)[:p] = p
    var(pm, nw)[:q] = q
    report && _IM.sol_component_value_edge(pm, nw, :transformer, :pf, :pt, ref(pm, nw, :arcs_from_trans), ref(pm, nw, :arcs_to_trans), p)
    report && _IM.sol_component_value_edge(pm, nw, :transformer, :qf, :qt, ref(pm, nw, :arcs_from_trans), ref(pm, nw, :arcs_to_trans), q)
end


""
function variable_mc_load_setpoint(pm::_PM.AbstractIVRModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true, kwargs...)
    var(pm, nw)[:crd] = Dict{Int, Any}()
    var(pm, nw)[:cid] = Dict{Int, Any}()
    var(pm, nw)[:crd_bus] = Dict{Int, Any}()
    var(pm, nw)[:cid_bus] = Dict{Int, Any}()
end


""
function variable_mc_gen_power_setpoint(pm::_PM.AbstractIVRModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_gen_current_setpoint_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_gen_current_setpoint_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    var(pm, nw)[:crg_bus] = Dict{Int, Any}()
    var(pm, nw)[:cig_bus] = Dict{Int, Any}()

    # store active and reactive power expressions for use in objective + post processing
    var(pm, nw)[:pg] = Dict{Int, Any}()
    var(pm, nw)[:qg] = Dict{Int, Any}()
end


""
function variable_mc_bus_voltage(pm::_PM.AbstractIVRModel; nw=pm.cnw, bounded::Bool=true, kwargs...)
    variable_mc_bus_voltage_real(pm; nw=nw, bounded=bounded, kwargs...)
    variable_mc_bus_voltage_imaginary(pm; nw=nw, bounded=bounded, kwargs...)

    # local infeasbility issues without proper initialization;
    # convergence issues start when the equivalent angles of the starting point
    # are further away than 90 degrees from the solution (as given by ACP)
    # this is the default behaviour of _PM, initialize all phases as (1,0)
    # the magnitude seems to have little effect on the convergence (>0.05)
    # updating the starting point to a balanced phasor does the job

    ncnd = length(conductor_ids(pm))
    theta = [_wrap_to_pi(2 * pi / ncnd * (1-c)) for c in 1:ncnd]

    for id in ids(pm, nw, :bus)
        busref = ref(pm, nw, :bus, id)
        vm = haskey(busref, "vm_start") ? busref["vm_start"] : fill(1.0, ncnd)
        if !haskey(busref, "va_start")
            for c in 1:ncnd
                vr = vm[c]*cos(theta[c])
                vi = vm[c]*sin(theta[c])
                JuMP.set_start_value(var(pm, nw, :vr, id)[c], vr)
                JuMP.set_start_value(var(pm, nw, :vi, id)[c], vi)
            end
        end
    end

    # apply bounds if bounded
    if bounded
        for i in ids(pm, nw, :bus)
            constraint_mc_voltage_magnitude_bounds(pm, i, nw=nw)
        end
    end
end


"Defines how current distributes over series and shunt impedances of a pi-model branch"
function constraint_mc_current_from(pm::_PM.AbstractIVRModel, n::Int, f_bus, f_idx, g_sh_fr, b_sh_fr, tr, ti, tm)
    vr_fr = var(pm, n, :vr, f_bus)
    vi_fr = var(pm, n, :vi, f_bus)

    csr_fr =  var(pm, n, :csr, f_idx[1])
    csi_fr =  var(pm, n, :csi, f_idx[1])

    cr_fr =  var(pm, n, :cr, f_idx)
    ci_fr =  var(pm, n, :ci, f_idx)

    tr = tr
    ti = ti
    tm = tm
    g_sh_fr = g_sh_fr
    b_sh_fr = b_sh_fr

    JuMP.@constraint(pm.model, cr_fr .== (tr.*csr_fr - ti.*csi_fr + g_sh_fr*vr_fr - b_sh_fr*vi_fr)./tm.^2)
    JuMP.@constraint(pm.model, ci_fr .== (tr.*csi_fr + ti.*csr_fr + g_sh_fr*vi_fr + b_sh_fr*vr_fr)./tm.^2)
end


"Defines how current distributes over series and shunt impedances of a pi-model branch"
function constraint_mc_current_to(pm::_PM.AbstractIVRModel, n::Int, t_bus, f_idx, t_idx, g_sh_to, b_sh_to)
    vr_to = var(pm, n, :vr, t_bus)
    vi_to = var(pm, n, :vi, t_bus)

    csr_to =  -var(pm, n, :csr, f_idx[1])
    csi_to =  -var(pm, n, :csi, f_idx[1])

    cr_to =  var(pm, n, :cr, t_idx)
    ci_to =  var(pm, n, :ci, t_idx)

    g_sh_to = g_sh_to
    b_sh_to = b_sh_to

    JuMP.@constraint(pm.model, cr_to .== csr_to + g_sh_to*vr_to - b_sh_to*vi_to)
    JuMP.@constraint(pm.model, ci_to .== csi_to + g_sh_to*vi_to + b_sh_to*vr_to)
end


"Defines voltage drop over a branch, linking from and to side complex voltage"
function constraint_mc_bus_voltage_drop(pm::_PM.AbstractIVRModel, n::Int, i, f_bus, t_bus, f_idx, r, x, tr, ti, tm)
    vr_fr = var(pm, n, :vr, f_bus)
    vi_fr = var(pm, n, :vi, f_bus)

    vr_to = var(pm, n, :vr, t_bus)
    vi_to = var(pm, n, :vi, t_bus)

    csr_fr =  var(pm, n, :csr, f_idx[1])
    csi_fr =  var(pm, n, :csi, f_idx[1])

    r = r
    x = x

    JuMP.@constraint(pm.model, vr_to .== (vr_fr.*tr + vi_fr.*ti)./tm.^2 - r*csr_fr + x*csi_fr)
    JuMP.@constraint(pm.model, vi_to .== (vi_fr.*tr - vr_fr.*ti)./tm.^2 - r*csi_fr - x*csr_fr)
end


"Bounds the voltage angle difference between bus pairs"
function constraint_mc_voltage_angle_difference(pm::_PM.AbstractIVRModel, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx

    vr_fr = var(pm, n, :vr, f_bus)
    vi_fr = var(pm, n, :vi, f_bus)
    vr_to = var(pm, n, :vr, t_bus)
    vi_to = var(pm, n, :vi, t_bus)
    vvr = vr_fr.*vr_to + vi_fr.*vi_to
    vvi = vi_fr.*vr_to - vr_fr.*vi_to

    JuMP.@constraint(pm.model, tan.(angmin).*vvr .<= vvi)
    JuMP.@constraint(pm.model, tan.(angmax).*vvr .>= vvi)
end


"""
Kirchhoff's current law applied to buses
`sum(cr + im*ci) = 0`
"""
function constraint_mc_load_current_balance(pm::_PM.AbstractIVRModel, n::Int, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_gs, bus_bs)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    cr    = get(var(pm, n),    :cr, Dict()); _PM._check_var_keys(cr, bus_arcs, "real current", "branch")
    ci    = get(var(pm, n),    :ci, Dict()); _PM._check_var_keys(ci, bus_arcs, "imaginary current", "branch")
    crd   = get(var(pm, n),   :crd_bus, Dict()); _PM._check_var_keys(crd, bus_loads, "real current", "load")
    cid   = get(var(pm, n),   :cid_bus, Dict()); _PM._check_var_keys(cid, bus_loads, "imaginary current", "load")
    crg   = get(var(pm, n),   :crg_bus, Dict()); _PM._check_var_keys(crg, bus_gens, "real current", "generator")
    cig   = get(var(pm, n),   :cig_bus, Dict()); _PM._check_var_keys(cig, bus_gens, "imaginary current", "generator")
    crs   = get(var(pm, n),   :crs, Dict()); _PM._check_var_keys(crs, bus_storage, "real currentr", "storage")
    cis   = get(var(pm, n),   :cis, Dict()); _PM._check_var_keys(cis, bus_storage, "imaginary current", "storage")
    crsw  = get(var(pm, n),  :crsw, Dict()); _PM._check_var_keys(crsw, bus_arcs_sw, "real current", "switch")
    cisw  = get(var(pm, n),  :cisw, Dict()); _PM._check_var_keys(cisw, bus_arcs_sw, "imaginary current", "switch")
    crt   = get(var(pm, n),   :crt, Dict()); _PM._check_var_keys(crt, bus_arcs_trans, "real current", "transformer")
    cit   = get(var(pm, n),   :cit, Dict()); _PM._check_var_keys(cit, bus_arcs_trans, "imaginary current", "transformer")

    cnds = conductor_ids(pm; nw=n)
    ncnds = length(cnds)

    Gt = isempty(bus_gs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_gs))
    Bt = isempty(bus_bs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_bs))

    for c in cnds
        JuMP.@NLconstraint(pm.model,  sum(cr[a][c] for a in bus_arcs)
                                    + sum(crsw[a_sw][c] for a_sw in bus_arcs_sw)
                                    + sum(crt[a_trans][c] for a_trans in bus_arcs_trans)
                                    ==
                                      sum(crg[g][c]         for g in bus_gens)
                                    - sum(crs[s][c]         for s in bus_storage)
                                    - sum(crd[d][c]         for d in bus_loads)
                                    - sum( Gt[c,d]*vr[d] -Bt[c,d]*vi[d] for d in cnds) # shunts
                                    )
        JuMP.@NLconstraint(pm.model, sum(ci[a][c] for a in bus_arcs)
                                    + sum(cisw[a_sw][c] for a_sw in bus_arcs_sw)
                                    + sum(cit[a_trans][c] for a_trans in bus_arcs_trans)
                                    ==
                                      sum(cig[g][c]         for g in bus_gens)
                                    - sum(cis[s][c]         for s in bus_storage)
                                    - sum(cid[d][c]         for d in bus_loads)
                                    - sum( Gt[c,d]*vi[d] +Bt[c,d]*vr[d] for d in cnds) # shunts
                                    )
    end
end


"`p[f_idx]^2 + q[f_idx]^2 <= rate_a^2`"
function constraint_mc_thermal_limit_from(pm::_PM.AbstractIVRModel, n::Int, f_idx, rate_a)
    (l, f_bus, t_bus) = f_idx

    vr = var(pm, n, :vr, f_bus)
    vi = var(pm, n, :vi, f_bus)
    crf = var(pm, n, :cr, f_idx)
    cif = var(pm, n, :ci, f_idx)

    cnds = conductor_ids(pm; nw=n)
    ncnds = length(cnds)

    for c in cnds
        JuMP.@NLconstraint(pm.model, (vr[c]^2 + vi[c]^2)*(crf[c]^2 + cif[c]^2) <= rate_a[c]^2)
    end
end


"`p[t_idx]^2 + q[t_idx]^2 <= rate_a^2`"
function constraint_mc_thermal_limit_to(pm::_PM.AbstractIVRModel, n::Int, t_idx, rate_a)
    (l, t_bus, f_bus) = t_idx

    vr = var(pm, n, :vr, t_bus)
    vi = var(pm, n, :vi, t_bus)
    crt = var(pm, n, :cr, t_idx)
    cit = var(pm, n, :ci, t_idx)

    cnds = conductor_ids(pm; nw=n)
    ncnds = length(cnds)

    for c in cnds
        JuMP.@NLconstraint(pm.model, (vr[c]^2 + vi[c]^2)*(crt[c]^2 + cit[c]^2) <= rate_a[c]^2)
    end
end


"""
Bounds the current magnitude at both from and to side of a branch
`cr[f_idx]^2 + ci[f_idx]^2 <= c_rating_a^2`
`cr[t_idx]^2 + ci[t_idx]^2 <= c_rating_a^2`
"""
function constraint_mc_current_limit(pm::_PM.AbstractIVRModel, n::Int, f_idx, c_rating_a)
    (l, f_bus, t_bus) = f_idx
    t_idx = (l, t_bus, f_bus)

    crf =  var(pm, n, :cr, f_idx)
    cif =  var(pm, n, :ci, f_idx)

    crt =  var(pm, n, :cr, t_idx)
    cit =  var(pm, n, :ci, t_idx)

    JuMP.@constraint(pm.model, crf.^2 + cif.^2 .<= c_rating_a.^2)
    JuMP.@constraint(pm.model, crt.^2 + cit.^2 .<= c_rating_a.^2)
end


"""
`pmin <= Re(v*cg') <= pmax`
"""
function constraint_mc_gen_active_bounds(pm::_PM.AbstractIVRModel, n::Int, i, bus, pmax, pmin)
    @assert pmin <= pmax

    vr = var(pm, n, :vr, bus)
    vi = var(pm, n, :vi, bus)
    cr = var(pm, n, :crg, i)
    ci = var(pm, n, :cig, i)

    JuMP.@constraint(pm.model, pmin .<= vr.*cr  + vi.*ci)
    JuMP.@constraint(pm.model, pmax .>= vr.*cr  + vi.*ci)
end


"""
`qmin <= Im(v*cg') <= qmax`
"""
function constraint_mc_gen_reactive_bounds(pm::_PM.AbstractIVRModel, n::Int, i, bus, qmax, qmin)
    @assert qmin <= qmax

    vr = var(pm, n, :vr, bus)
    vi = var(pm, n, :vi, bus)
    cr = var(pm, n, :crg, i)
    ci = var(pm, n, :cig, i)

    JuMP.@constraint(pm.model, qmin .<= vi.*cr  - vr.*ci)
    JuMP.@constraint(pm.model, qmax .>= vi.*cr  - vr.*ci)
end


"`pg[i] == pg`"
function constraint_mc_gen_power_setpoint_real(pm::_PM.AbstractIVRModel, n::Int, i, pgref)
    gen = ref(pm, n, :gen, i)
    bus = gen["gen_bus"]
    vr = var(pm, n, :vr, bus)
    vi = var(pm, n, :vi, bus)
    cr = var(pm, n, :crg, i)
    ci = var(pm, n, :cig, i)

    JuMP.@constraint(pm.model, pgref .== vr.*cr  + vi.*ci)
end


"`qq[i] == qq`"
function constraint_mc_regen_setpoint_active(pm::_PM.AbstractIVRModel, n::Int, i, qgref)
    gen = ref(pm, n, :gen, i)
    bus = gen["gen_bus"]
    vr = var(pm, n, :vr, bus)
    vi = var(pm, n, :vi, bus)
    cr = var(pm, n, :crg, i)
    ci = var(pm, n, :cig, i)

    JuMP.@constraint(pm.model, qgref .== vi.*cr  - vr.*ci)
end


"wye-wye transformer power constraint for IVR formulation"
function constraint_mc_transformer_power_yy(pm::_PM.AbstractIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
    vr_fr_P = var(pm, nw, :vr, f_bus)[f_cnd]
    vi_fr_P = var(pm, nw, :vi, f_bus)[f_cnd]
    vr_fr_n = 0
    vi_fr_n = 0
    vr_to_P = var(pm, nw, :vr, t_bus)[t_cnd]
    vi_to_P = var(pm, nw, :vi, t_bus)[t_cnd]
    vr_to_n = 0
    vi_to_n = 0

    cr_fr_P = var(pm, nw, :crt, f_idx)[f_cnd]
    ci_fr_P = var(pm, nw, :cit, f_idx)[f_cnd]
    cr_to_P = var(pm, nw, :crt, t_idx)[t_cnd]
    ci_to_P = var(pm, nw, :cit, t_idx)[t_cnd]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[p] ? tm_set[p] : var(pm, nw, :tap, trans_id)[p] for p in conductor_ids(pm)]
    scale = (tm_scale*pol).*tm_set

    JuMP.@constraint(pm.model, (vr_fr_P.-vr_fr_n) .== scale.*(vr_to_P.-vr_to_n))
    JuMP.@constraint(pm.model, (vi_fr_P.-vi_fr_n) .== scale.*(vi_to_P.-vi_to_n))
    JuMP.@constraint(pm.model, scale.*cr_fr_P .+ cr_to_P .== 0)
    JuMP.@constraint(pm.model, scale.*ci_fr_P .+ ci_to_P .== 0)
end


"delta-wye transformer power constraint for IVR formulation"
function constraint_mc_transformer_power_dy(pm::_PM.AbstractIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
    vr_fr_P = var(pm, nw, :vr, f_bus)[f_cnd]
    vi_fr_P = var(pm, nw, :vi, f_bus)[f_cnd]
    vr_to_P = var(pm, nw, :vr, t_bus)[t_cnd]
    vi_to_P = var(pm, nw, :vi, t_bus)[t_cnd]
    vr_to_n = 0
    vi_to_n = 0

    cr_fr_P = var(pm, nw, :crt, f_idx)[f_cnd]
    ci_fr_P = var(pm, nw, :cit, f_idx)[f_cnd]
    cr_to_P = var(pm, nw, :crt, t_idx)[t_cnd]
    ci_to_P = var(pm, nw, :cit, t_idx)[t_cnd]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[p] ? tm_set[p] : var(pm, nw, :tap, trans_id)[p] for p in conductor_ids(pm)]
    scale = (tm_scale*pol).*tm_set

    n_phases = length(tm)
    Md = _get_delta_transformation_matrix(n_phases)

    JuMP.@constraint(pm.model, Md*vr_fr_P .== scale.*(vr_to_P .- vr_to_n))
    JuMP.@constraint(pm.model, Md*vi_fr_P .== scale.*(vi_to_P .- vi_to_n))

    JuMP.@constraint(pm.model, scale.*cr_fr_P .+ Md'*cr_to_P .== 0)
    JuMP.@constraint(pm.model, scale.*ci_fr_P .+ Md'*ci_to_P .== 0)
end


"wye connected load setpoint constraint for IVR formulation"
function constraint_mc_load_setpoint_wye(pm::_PM.IVRPowerModel, nw::Int, id::Int, bus_id::Int, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    nph = 3

    crd = JuMP.@NLexpression(pm.model, [i in 1:nph],
        a[i]*vr[i]*(vr[i]^2+vi[i]^2)^(alpha[i]/2-1)
       +b[i]*vi[i]*(vr[i]^2+vi[i]^2)^(beta[i]/2 -1)
    )
    cid = JuMP.@NLexpression(pm.model, [i in 1:nph],
        a[i]*vi[i]*(vr[i]^2+vi[i]^2)^(alpha[i]/2-1)
       -b[i]*vr[i]*(vr[i]^2+vi[i]^2)^(beta[i]/2 -1)
    )

    var(pm, nw, :crd_bus)[id] = crd
    var(pm, nw, :cid_bus)[id] = cid

    if report
        pd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph],  vr[i]*crd[i]+vi[i]*cid[i])
        qd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], -vr[i]*cid[i]+vi[i]*crd[i])

        sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        sol(pm, nw, :load, id)[:crd_bus] = crd
        sol(pm, nw, :load, id)[:cid_bus] = cid

        pd = JuMP.@NLexpression(pm.model, [i in 1:nph], a[i]*(vr[i]^2+vi[i]^2)^(alpha[i]/2) )
        qd = JuMP.@NLexpression(pm.model, [i in 1:nph], b[i]*(vr[i]^2+vi[i]^2)^(beta[i]/2)  )
        sol(pm, nw, :load, id)[:pd] = pd
        sol(pm, nw, :load, id)[:qd] = qd
    end
end


"delta connected load setpoint constraint for IVR formulation"
function constraint_mc_load_setpoint_delta(pm::_PM.IVRPowerModel, nw::Int, id::Int, bus_id::Int, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    nph = 3
    prev = Dict(i=>(i+nph-2)%nph+1 for i in 1:nph)
    next = Dict(i=>i%nph+1 for i in 1:nph)

    vrd = JuMP.@NLexpression(pm.model, [i in 1:nph], vr[i]-vr[next[i]])
    vid = JuMP.@NLexpression(pm.model, [i in 1:nph], vi[i]-vi[next[i]])

    crd = JuMP.@NLexpression(pm.model, [i in 1:nph],
        a[i]*vrd[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2-1)
       +b[i]*vid[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2 -1)
    )
    cid = JuMP.@NLexpression(pm.model, [i in 1:nph],
        a[i]*vid[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2-1)
       -b[i]*vrd[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2 -1)
    )

    crd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], crd[i]-crd[prev[i]])
    cid_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], cid[i]-cid[prev[i]])

    var(pm, nw, :crd_bus)[id] = crd_bus
    var(pm, nw, :cid_bus)[id] = cid_bus

    if report
        pd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph],  vr[i]*crd_bus[i]+vi[i]*cid_bus[i])
        qd_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], -vr[i]*cid_bus[i]+vi[i]*crd_bus[i])

        sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        pd = JuMP.@NLexpression(pm.model, [i in 1:nph], a[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2) )
        qd = JuMP.@NLexpression(pm.model, [i in 1:nph], b[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2)  )
        sol(pm, nw, :load, id)[:pd] = pd
        sol(pm, nw, :load, id)[:qd] = qd
    end
end


"wye connected generator setpoint constraint for IVR formulation"
function constraint_mc_gen_setpoint_wye(pm::_PM.IVRPowerModel, nw::Int, id::Int, bus_id::Int, pmin::Vector, pmax::Vector, qmin::Vector, qmax::Vector; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)

    nph = 3

    pg = JuMP.@NLexpression(pm.model, [i in 1:nph],  vr[i]*crg[i]+vi[i]*cig[i])
    qg = JuMP.@NLexpression(pm.model, [i in 1:nph], -vr[i]*cig[i]+vi[i]*crg[i])

    if bounded
        for c in 1:nph
            if pmin[c]>-Inf
                JuMP.@constraint(pm.model, pmin[c] .<= vr[c]*crg[c]  + vi[c]*cig[c])
            end
            if pmax[c]< Inf
                JuMP.@constraint(pm.model, pmax[c] .>= vr[c]*crg[c]  + vi[c]*cig[c])
            end
            if qmin[c]>-Inf
                JuMP.@constraint(pm.model, qmin[c] .<= vi[c]*crg[c]  - vr[c]*cig[c])
            end
            if qmax[c]< Inf
                JuMP.@constraint(pm.model, qmax[c] .>= vi[c]*crg[c]  - vr[c]*cig[c])
            end
        end
    end

    var(pm, nw, :crg_bus)[id] = crg
    var(pm, nw, :cig_bus)[id] = cig
    var(pm, nw, :pg)[id] = pg
    var(pm, nw, :qg)[id] = qg

    if report
        sol(pm, nw, :gen, id)[:crg_bus] = var(pm, nw, :crg_bus, id)
        sol(pm, nw, :gen, id)[:cig_bus] = var(pm, nw, :cig_bus, id)

        sol(pm, nw, :gen, id)[:pg] = pg
        sol(pm, nw, :gen, id)[:qg] = qg
    end
end


"delta connected generator setpoint constraint for IVR formulation"
function constraint_mc_gen_setpoint_delta(pm::_PM.IVRPowerModel, nw::Int, id::Int, bus_id::Int, pmin::Vector, pmax::Vector, qmin::Vector, qmax::Vector; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)

    nph = 3
    prev = Dict(i=>(i+nph-2)%nph+1 for i in 1:nph)
    next = Dict(i=>i%nph+1 for i in 1:nph)

    vrg = JuMP.@NLexpression(pm.model, [i in 1:nph], vr[i]-vr[next[i]])
    vig = JuMP.@NLexpression(pm.model, [i in 1:nph], vi[i]-vi[next[i]])

    pg = JuMP.@NLexpression(pm.model, [i in 1:nph],  vrg[i]*crg[i]+vig[i]*cig[i])
    qg = JuMP.@NLexpression(pm.model, [i in 1:nph], -vrg[i]*cig[i]+vig[i]*crg[i])

    if bounded
        JuMP.@NLconstraint(pm.model, [i in 1:nph], pmin[i] <= pg[i])
        JuMP.@NLconstraint(pm.model, [i in 1:nph], pmax[i] >= pg[i])
        JuMP.@NLconstraint(pm.model, [i in 1:nph], qmin[i] <= qg[i])
        JuMP.@NLconstraint(pm.model, [i in 1:nph], qmax[i] >= qg[i])
    end

    crg_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], crg[i]-crg[prev[i]])
    cig_bus = JuMP.@NLexpression(pm.model, [i in 1:nph], cig[i]-cig[prev[i]])

    var(pm, nw, :crg_bus)[id] = crg_bus
    var(pm, nw, :cig_bus)[id] = cig_bus
    var(pm, nw, :pg)[id] = pg
    var(pm, nw, :qg)[id] = qg

    if report
        sol(pm, nw, :gen, id)[:crg_bus] = crg_bus
        sol(pm, nw, :gen, id)[:cig_bus] = cig_bus
        sol(pm, nw, :gen, id)[:pg] = pg
        sol(pm, nw, :gen, id)[:qg] = qg
    end
end
