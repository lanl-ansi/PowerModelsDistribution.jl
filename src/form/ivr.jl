# Kirchhoff's circuit laws as defined the current-voltage variable space.
# Even though the branch model is linear, the feasible set is non-convex
# in the context of constant-power loads or generators

""
function variable_mc_branch_current(pm::_PMs.AbstractIVRModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    # store expressions in rectangular power variable space
    p = Dict()
    q = Dict()

    for (l,i,j) in _PMs.ref(pm, nw, :arcs_from)
        vr_fr = _PMs.var(pm, nw, :vr, i)
        vi_fr = _PMs.var(pm, nw, :vi, i)
        cr_fr = _PMs.var(pm, nw, :cr, (l,i,j))
        ci_fr = _PMs.var(pm, nw, :ci, (l,i,j))

        vr_to = _PMs.var(pm, nw, :vr, j)
        vi_to = _PMs.var(pm, nw, :vi, j)
        cr_to = _PMs.var(pm, nw, :cr, (l,j,i))
        ci_to = _PMs.var(pm, nw, :ci, (l,j,i))
        p[(l,i,j)] = vr_fr.*cr_fr  + vi_fr.*ci_fr
        q[(l,i,j)] = vi_fr.*cr_fr  - vr_fr.*ci_fr
        p[(l,j,i)] = vr_to.*cr_to  + vi_to.*ci_to
        q[(l,j,i)] = vi_to.*cr_to  - vr_to.*ci_to
    end

    _PMs.var(pm, nw)[:p] = p
    _PMs.var(pm, nw)[:q] = q
    report && _PMs.sol_component_value_edge(pm, nw, :branch, :pf, :pt, _PMs.ref(pm, nw, :arcs_from), _PMs.ref(pm, nw, :arcs_to), p)
    report && _PMs.sol_component_value_edge(pm, nw, :branch, :qf, :qt, _PMs.ref(pm, nw, :arcs_from), _PMs.ref(pm, nw, :arcs_to), q)

    variable_mc_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end


""
function variable_mc_transformer_current(pm::_PMs.AbstractIVRModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_transformer_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_transformer_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    # store expressions in rectangular power variable space
    p = Dict()
    q = Dict()

    f_cnd = [1, 2, 3] #TODO extend to this being passed in constraint template
    t_cnd = [1, 2, 3] #TODO extend to this being passed in constraint template

    for (l,i,j) in _PMs.ref(pm, nw, :arcs_from_trans)
        vr_fr = _PMs.var(pm, nw, :vr, i)[f_cnd]
        vi_fr = _PMs.var(pm, nw, :vi, i)[f_cnd]
        cr_fr = _PMs.var(pm, nw, :crt, (l,i,j))
        ci_fr = _PMs.var(pm, nw, :cit, (l,i,j))

        vr_to = _PMs.var(pm, nw, :vr, j)[t_cnd]
        vi_to = _PMs.var(pm, nw, :vi, j)[t_cnd]
        cr_to = _PMs.var(pm, nw, :crt, (l,j,i))
        ci_to = _PMs.var(pm, nw, :cit, (l,j,i))
        p[(l,i,j)] = vr_fr.*cr_fr  + vi_fr.*ci_fr
        q[(l,i,j)] = vi_fr.*cr_fr  - vr_fr.*ci_fr
        p[(l,j,i)] = vr_to.*cr_to  + vi_to.*ci_to
        q[(l,j,i)] = vi_to.*cr_to  - vr_to.*ci_to
    end

    _PMs.var(pm, nw)[:p] = p
    _PMs.var(pm, nw)[:q] = q
    report && _PMs.sol_component_value_edge(pm, nw, :transformer, :pf, :pt, _PMs.ref(pm, nw, :arcs_from_trans), _PMs.ref(pm, nw, :arcs_to_trans), p)
    report && _PMs.sol_component_value_edge(pm, nw, :transformer, :qf, :qt, _PMs.ref(pm, nw, :arcs_from_trans), _PMs.ref(pm, nw, :arcs_to_trans), q)
end


""
function variable_mc_generation(pm::_PMs.AbstractIVRModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_generation_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_generation_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    # store active and reactive power expressions for use in objective + post processing
    pg = Dict()
    qg = Dict()
    for (i,gen) in _PMs.ref(pm, nw, :gen)
        busid = gen["gen_bus"]
        vr = _PMs.var(pm, nw, :vr, busid)
        vi = _PMs.var(pm, nw, :vi, busid)
        crg = _PMs.var(pm, nw, :crg, i)
        cig = _PMs.var(pm, nw, :cig, i)

        pg[i] = []
        qg[i] = []
        for c in 1:ncnds
            push!(pg[i], JuMP.@NLexpression(pm.model, vr[c]*crg[c]  + vi[c]*cig[c]))
            push!(qg[i], JuMP.@NLexpression(pm.model, vi[c]*crg[c]  - vr[c]*cig[c]))
        end
    end
    _PMs.var(pm, nw)[:pg] = pg
    _PMs.var(pm, nw)[:qg] = qg
    report && _PMs.sol_component_value(pm, nw, :gen, :pg, _PMs.ids(pm, nw, :gen), pg)
    report && _PMs.sol_component_value(pm, nw, :gen, :qg, _PMs.ids(pm, nw, :gen), qg)

    if bounded
        for (i,gen) in _PMs.ref(pm, nw, :gen)
            constraint_mc_generation_active_power_limits(pm, i, nw=nw)
            constraint_mc_generation_reactive_power_limits(pm, i, nw=nw)
        end
    end
end


""
function variable_mc_voltage(pm::_PMs.AbstractIVRModel; nw=pm.cnw, bounded::Bool=true, kwargs...)
    variable_mc_voltage_real(pm; nw=nw, bounded=bounded, kwargs...)
    variable_mc_voltage_imaginary(pm; nw=nw, bounded=bounded, kwargs...)

    # local infeasbility issues without proper initialization;
    # convergence issues start when the equivalent angles of the starting point
    # are further away than 90 degrees from the solution (as given by ACP)
    # this is the default behaviour of _PMs, initialize all phases as (1,0)
    # the magnitude seems to have little effect on the convergence (>0.05)
    # updating the starting point to a balanced phasor does the job

    ncnd = length(_PMs.conductor_ids(pm))
    theta = [_wrap_to_pi(2 * pi / ncnd * (1-c)) for c in 1:ncnd]
    vm = 1
    for id in _PMs.ids(pm, nw, :bus)
        busref = _PMs.ref(pm, nw, :bus, id)
        if !haskey(busref, "va_start")
            for c in 1:ncnd
                vr = vm*cos(theta[c])
                vi = vm*sin(theta[c])
                JuMP.set_start_value(_PMs.var(pm, nw, :vr, id)[c], vr)
                JuMP.set_start_value(_PMs.var(pm, nw, :vi, id)[c], vi)
            end
        end
    end

    # apply bounds if bounded
    if bounded
        for i in _PMs.ids(pm, nw, :bus)
            constraint_mc_voltage_magnitude_bounds(pm, i, nw=nw)
        end
    end
end

"""
Defines how current distributes over series and shunt impedances of a pi-model branch
"""
function constraint_mc_current_from(pm::_PMs.AbstractIVRModel, n::Int, f_bus, f_idx, g_sh_fr, b_sh_fr, tr, ti, tm)
    vr_fr = _PMs.var(pm, n, :vr, f_bus)
    vi_fr = _PMs.var(pm, n, :vi, f_bus)

    csr_fr =  _PMs.var(pm, n, :csr, f_idx[1])
    csi_fr =  _PMs.var(pm, n, :csi, f_idx[1])

    cr_fr =  _PMs.var(pm, n, :cr, f_idx)
    ci_fr =  _PMs.var(pm, n, :ci, f_idx)

    tr = tr
    ti = ti
    tm = tm
    g_sh_fr = g_sh_fr
    b_sh_fr = b_sh_fr

    JuMP.@constraint(pm.model, cr_fr .== (tr.*csr_fr - ti.*csi_fr + g_sh_fr*vr_fr - b_sh_fr*vi_fr)./tm.^2)
    JuMP.@constraint(pm.model, ci_fr .== (tr.*csi_fr + ti.*csr_fr + g_sh_fr*vi_fr + b_sh_fr*vr_fr)./tm.^2)
end

"""
Defines how current distributes over series and shunt impedances of a pi-model branch
"""
function constraint_mc_current_to(pm::_PMs.AbstractIVRModel, n::Int, t_bus, f_idx, t_idx, g_sh_to, b_sh_to)
    vr_to = _PMs.var(pm, n, :vr, t_bus)
    vi_to = _PMs.var(pm, n, :vi, t_bus)

    csr_to =  -_PMs.var(pm, n, :csr, f_idx[1])
    csi_to =  -_PMs.var(pm, n, :csi, f_idx[1])

    cr_to =  _PMs.var(pm, n, :cr, t_idx)
    ci_to =  _PMs.var(pm, n, :ci, t_idx)

    g_sh_to = g_sh_to
    b_sh_to = b_sh_to

    JuMP.@constraint(pm.model, cr_to .== csr_to + g_sh_to*vr_to - b_sh_to*vi_to)
    JuMP.@constraint(pm.model, ci_to .== csi_to + g_sh_to*vi_to + b_sh_to*vr_to)
end


"""
Defines voltage drop over a branch, linking from and to side complex voltage
"""
function constraint_mc_voltage_drop(pm::_PMs.AbstractIVRModel, n::Int, i, f_bus, t_bus, f_idx, r, x, tr, ti, tm)
    vr_fr = _PMs.var(pm, n, :vr, f_bus)
    vi_fr = _PMs.var(pm, n, :vi, f_bus)

    vr_to = _PMs.var(pm, n, :vr, t_bus)
    vi_to = _PMs.var(pm, n, :vi, t_bus)

    csr_fr =  _PMs.var(pm, n, :csr, f_idx[1])
    csi_fr =  _PMs.var(pm, n, :csi, f_idx[1])

    r = r
    x = x

    JuMP.@constraint(pm.model, vr_to .== (vr_fr.*tr + vi_fr.*ti)./tm.^2 - r*csr_fr + x*csi_fr)
    JuMP.@constraint(pm.model, vi_to .== (vi_fr.*tr - vr_fr.*ti)./tm.^2 - r*csi_fr - x*csr_fr)
end

"""
Bounds the voltage angle difference between bus pairs
"""
function constraint_mc_voltage_angle_difference(pm::_PMs.AbstractIVRModel, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx

    vr_fr = _PMs.var(pm, n, :vr, f_bus)
    vi_fr = _PMs.var(pm, n, :vi, f_bus)
    vr_to = _PMs.var(pm, n, :vr, t_bus)
    vi_to = _PMs.var(pm, n, :vi, t_bus)
    vvr = vr_fr.*vr_to + vi_fr.*vi_to
    vvi = vi_fr.*vr_to - vr_fr.*vi_to

    JuMP.@constraint(pm.model, tan.(angmin).*vvr .<= vvi)
    JuMP.@constraint(pm.model, tan.(angmax).*vvr .>= vvi)
end

"""
Kirchhoff's current law applied to buses
`sum(cr + im*ci) = 0`
"""
function constraint_mc_current_balance(pm::_PMs.AbstractIVRModel, n::Int, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)

    cr    = get(_PMs.var(pm, n),    :cr, Dict()); _PMs._check_var_keys(cr, bus_arcs, "real current", "branch")
    ci    = get(_PMs.var(pm, n),    :ci, Dict()); _PMs._check_var_keys(ci, bus_arcs, "imaginary current", "branch")
    crg   = get(_PMs.var(pm, n),   :crg, Dict()); _PMs._check_var_keys(crg, bus_gens, "real current", "generator")
    cig   = get(_PMs.var(pm, n),   :cig, Dict()); _PMs._check_var_keys(cig, bus_gens, "imaginary current", "generator")
    crs   = get(_PMs.var(pm, n),   :crs, Dict()); _PMs._check_var_keys(crs, bus_storage, "real currentr", "storage")
    cis   = get(_PMs.var(pm, n),   :cis, Dict()); _PMs._check_var_keys(cis, bus_storage, "imaginary current", "storage")
    crsw  = get(_PMs.var(pm, n),  :crsw, Dict()); _PMs._check_var_keys(crsw, bus_arcs_sw, "real current", "switch")
    cisw  = get(_PMs.var(pm, n),  :cisw, Dict()); _PMs._check_var_keys(cisw, bus_arcs_sw, "imaginary current", "switch")
    crt   = get(_PMs.var(pm, n),   :crt, Dict()); _PMs._check_var_keys(crt, bus_arcs_trans, "real current", "transformer")
    cit   = get(_PMs.var(pm, n),   :cit, Dict()); _PMs._check_var_keys(cit, bus_arcs_trans, "imaginary current", "transformer")

    cnds = _PMs.conductor_ids(pm; nw=n)
    ncnds = length(cnds)

    for c in cnds
        JuMP.@NLconstraint(pm.model, sum(cr[a][c] for a in bus_arcs)
                                    + sum(crsw[a_sw][c] for a_sw in bus_arcs_sw)
                                    + sum(crt[a_trans][c] for a_trans in bus_arcs_trans)
                                    ==
                                    sum(crg[g][c] for g in bus_gens)
                                    - sum(crs[s][c] for s in bus_storage)
                                    - (sum(pd[c] for pd in values(bus_pd))*vr[c] + sum(qd[c] for qd in values(bus_qd))*vi[c])/(vr[c]^2 + vi[c]^2)
                                    - sum(gs[c] for gs in values(bus_gs))*vr[c] + sum(bs[c] for bs in values(bus_bs))*vi[c]
                                    )
        JuMP.@NLconstraint(pm.model, sum(ci[a][c] for a in bus_arcs)
                                    + sum(cisw[a_sw][c] for a_sw in bus_arcs_sw)
                                    + sum(cit[a_trans][c] for a_trans in bus_arcs_trans)
                                    ==
                                    sum(cig[g][c] for g in bus_gens)
                                    - sum(cis[s][c] for s in bus_storage)
                                    - (sum(pd[c] for pd in values(bus_pd))*vi[c] - sum(qd[c] for qd in values(bus_qd))*vr[c])/(vr[c]^2 + vi[c]^2)
                                    - sum(gs[c] for gs in values(bus_gs))*vi[c] - sum(bs[c] for bs in values(bus_bs))*vr[c]
                                    )
    end
end

"`p[f_idx]^2 + q[f_idx]^2 <= rate_a^2`"
function constraint_mc_thermal_limit_from(pm::_PMs.AbstractIVRModel, n::Int, f_idx, rate_a)
    (l, f_bus, t_bus) = f_idx

    vr = _PMs.var(pm, n, :vr, f_bus)
    vi = _PMs.var(pm, n, :vi, f_bus)
    crf = _PMs.var(pm, n, :cr, f_idx)
    cif = _PMs.var(pm, n, :ci, f_idx)

    cnds = _PMs.conductor_ids(pm; nw=n)
    ncnds = length(cnds)

    for c in cnds
        JuMP.@NLconstraint(pm.model, (vr[c]^2 + vi[c]^2)*(crf[c]^2 + cif[c]^2) <= rate_a[c]^2)
    end
end

"`p[t_idx]^2 + q[t_idx]^2 <= rate_a^2`"
function constraint_mc_thermal_limit_to(pm::_PMs.AbstractIVRModel, n::Int, t_idx, rate_a)
    (l, t_bus, f_bus) = t_idx

    vr = _PMs.var(pm, n, :vr, t_bus)
    vi = _PMs.var(pm, n, :vi, t_bus)
    crt = _PMs.var(pm, n, :cr, t_idx)
    cit = _PMs.var(pm, n, :ci, t_idx)

    cnds = _PMs.conductor_ids(pm; nw=n)
    ncnds = length(cnds)

    for c in cnds
        JuMP.@NLconstraint(pm.model, (vr[c]^2 + vi[c]^2)*(crt[c]^2 + cit[c]^2) <= rate_a[c]^2)
    end
end

"""
Bounds the current magnitude at both from and to side of a branch
`cr[f_idx]^2 + ci[f_idx]^2 <= c_rating^2`
`cr[t_idx]^2 + ci[t_idx]^2 <= c_rating^2`
"""
function constraint_mc_current_limit(pm::_PMs.AbstractIVRModel, n::Int, f_idx, c_rating)
    (l, f_bus, t_bus) = f_idx
    t_idx = (l, t_bus, f_bus)

    crf =  _PMs.var(pm, n, :cr, f_idx)
    cif =  _PMs.var(pm, n, :ci, f_idx)

    crt =  _PMs.var(pm, n, :cr, t_idx)
    cit =  _PMs.var(pm, n, :ci, t_idx)

    JuMP.@constraint(pm.model, crf.^2 + cif.^2 .<= c_rating.^2)
    JuMP.@constraint(pm.model, crt.^2 + cit.^2 .<= c_rating.^2)
end

"""
`pmin <= Re(v*cg') <= pmax`
"""
function constraint_mc_generation_active_power_limits(pm::_PMs.AbstractIVRModel, n::Int, i, bus, pmax, pmin)
    @assert pmin <= pmax

    vr = _PMs.var(pm, n, :vr, bus)
    vi = _PMs.var(pm, n, :vi, bus)
    cr = _PMs.var(pm, n, :crg, i)
    ci = _PMs.var(pm, n, :cig, i)

    JuMP.@constraint(pm.model, pmin .<= vr.*cr  + vi.*ci)
    JuMP.@constraint(pm.model, pmax .>= vr.*cr  + vi.*ci)
end

"""
`qmin <= Im(v*cg') <= qmax`
"""
function constraint_mc_generation_reactive_power_limits(pm::_PMs.AbstractIVRModel, n::Int, i, bus, qmax, qmin)
    @assert qmin <= qmax

    vr = _PMs.var(pm, n, :vr, bus)
    vi = _PMs.var(pm, n, :vi, bus)
    cr = _PMs.var(pm, n, :crg, i)
    ci = _PMs.var(pm, n, :cig, i)

    JuMP.@constraint(pm.model, qmin .<= vi.*cr  - vr.*ci)
    JuMP.@constraint(pm.model, qmax .>= vi.*cr  - vr.*ci)
end

"`pg[i] == pg`"
function constraint_mc_active_gen_setpoint(pm::_PMs.AbstractIVRModel, n::Int, i, pgref)
    gen = _PMs.ref(pm, n, :gen, i)
    bus = gen["gen_bus"]
    vr = _PMs.var(pm, n, :vr, bus)
    vi = _PMs.var(pm, n, :vi, bus)
    cr = _PMs.var(pm, n, :crg, i)
    ci = _PMs.var(pm, n, :cig, i)

    JuMP.@constraint(pm.model, pgref .== vr.*cr  + vi.*ci)
end

"`qq[i] == qq`"
function constraint_mc_reactive_gen_setpoint(pm::_PMs.AbstractIVRModel, n::Int, i, qgref)
    gen = _PMs.ref(pm, n, :gen, i)
    bus = gen["gen_bus"]
    vr = _PMs.var(pm, n, :vr, bus)
    vi = _PMs.var(pm, n, :vi, bus)
    cr = _PMs.var(pm, n, :crg, i)
    ci = _PMs.var(pm, n, :cig, i)

    JuMP.@constraint(pm.model, qgref .== vi.*cr  - vr.*ci)
end


function _PMs._objective_min_fuel_cost_polynomial_linquad(pm::_PMs.AbstractIVRModel; report::Bool=true)
    gen_cost = Dict()
    dcline_cost = Dict()

    for (n, nw_ref) in _PMs.nws(pm)
        for (i,gen) in nw_ref[:gen]
            bus = gen["gen_bus"]

            #to avoid function calls inside of @NLconstraint:
            pg = _PMs.var(pm, n, :pg, i)
            nc = length(_PMs.conductor_ids(pm, n))
            if length(gen["cost"]) == 1
                gen_cost[(n,i)] = gen["cost"][1]
            elseif length(gen["cost"]) == 2
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*sum(pg[c] for c in 1:nc) + gen["cost"][2])
            elseif length(gen["cost"]) == 3
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*sum(pg[c] for c in 1:nc)^2 + gen["cost"][2]*sum(pg[c] for c in 1:nc) + gen["cost"][3])
            else
                gen_cost[(n,i)] = 0.0
            end
        end
    end

    return JuMP.@NLobjective(pm.model, Min,
        sum(
            sum(    gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] )
            + sum( dcline_cost[(n,i)] for (i,dcline) in nw_ref[:dcline] )
        for (n, nw_ref) in _PMs.nws(pm))
    )
end


"adds pg_cost variables and constraints"
function objective_variable_pg_cost(pm::_PMs.AbstractIVRModel; report::Bool=true)
    for (n, nw_ref) in _PMs.nws(pm)
        gen_lines = calc_cost_pwl_lines(nw_ref[:gen])

        #to avoid function calls inside of @NLconstraint
        pg_cost = _PMs.var(pm, n)[:pg_cost] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, n, :gen)], base_name="$(n)_pg_cost",
        )
        report && _PMs.sol_component_value(pm, n, :gen, :pg_cost, _PMs.ids(pm, n, :gen), pg_cost)

        nc = length(conductor_ids(pm, n))

        # gen pwl cost
        for (i, gen) in nw_ref[:gen]
            pg = var(pm, n, :pg, i)
            for line in gen_lines[i]
                JuMP.@NLconstraint(pm.model, pg_cost[i] >= line.slope*sum(pg[c] for c in 1:nc) + line.intercept)
            end
        end
    end
end


function constraint_mc_trans_yy(pm::_PMs.AbstractIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
    vr_fr_P = _PMs.var(pm, nw, :vr, f_bus)[f_cnd]
    vi_fr_P = _PMs.var(pm, nw, :vi, f_bus)[f_cnd]
    vr_fr_n = 0
    vi_fr_n = 0
    vr_to_P = _PMs.var(pm, nw, :vr, t_bus)[t_cnd]
    vi_to_P = _PMs.var(pm, nw, :vi, t_bus)[t_cnd]
    vr_to_n = 0
    vi_to_n = 0

    cr_fr_P = _PMs.var(pm, nw, :crt, f_idx)[f_cnd]
    ci_fr_P = _PMs.var(pm, nw, :cit, f_idx)[f_cnd]
    cr_to_P = _PMs.var(pm, nw, :crt, t_idx)[t_cnd]
    ci_to_P = _PMs.var(pm, nw, :cit, t_idx)[t_cnd]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[p] ? tm_set[p] : _PMs.var(pm, nw, :tap, trans_id)[p] for p in _PMs.conductor_ids(pm)]
    scale = (tm_scale*pol).*tm_set

    JuMP.@constraint(pm.model, (vr_fr_P.-vr_fr_n) .== scale.*(vr_to_P.-vr_to_n))
    JuMP.@constraint(pm.model, (vi_fr_P.-vi_fr_n) .== scale.*(vi_to_P.-vi_to_n))
    JuMP.@constraint(pm.model, scale.*cr_fr_P .+ cr_to_P .== 0)
    JuMP.@constraint(pm.model, scale.*ci_fr_P .+ ci_to_P .== 0)
end


function constraint_mc_trans_dy(pm::_PMs.AbstractIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
    vr_fr_P = _PMs.var(pm, nw, :vr, f_bus)[f_cnd]
    vi_fr_P = _PMs.var(pm, nw, :vi, f_bus)[f_cnd]
    vr_to_P = _PMs.var(pm, nw, :vr, t_bus)[t_cnd]
    vi_to_P = _PMs.var(pm, nw, :vi, t_bus)[t_cnd]
    vr_to_n = 0
    vi_to_n = 0

    cr_fr_P = _PMs.var(pm, nw, :crt, f_idx)[f_cnd]
    ci_fr_P = _PMs.var(pm, nw, :cit, f_idx)[f_cnd]
    cr_to_P = _PMs.var(pm, nw, :crt, t_idx)[t_cnd]
    ci_to_P = _PMs.var(pm, nw, :cit, t_idx)[t_cnd]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[p] ? tm_set[p] : _PMs.var(pm, nw, :tap, trans_id)[p] for p in _PMs.conductor_ids(pm)]
    scale = (tm_scale*pol).*tm_set

    n_phases = length(tm)
    Md = _get_delta_transformation_matrix(n_phases)

    JuMP.@constraint(pm.model, Md*vr_fr_P .== scale.*(vr_to_P .- vr_to_n))
    JuMP.@constraint(pm.model, Md*vi_fr_P .== scale.*(vi_to_P .- vi_to_n))

    JuMP.@constraint(pm.model, scale.*cr_fr_P .+ Md'*cr_to_P .== 0)
    JuMP.@constraint(pm.model, scale.*ci_fr_P .+ Md'*ci_to_P .== 0)
end
