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
        f_connections = ref(pm, nw, :branch, l, "f_connections")
        t_connections = ref(pm, nw, :branch, l, "t_connections")

        vr_fr = [var(pm, nw, :vr, i)[c] for c in f_connections]
        vi_fr = [var(pm, nw, :vi, i)[c] for c in f_connections]
        cr_fr = [var(pm, nw, :cr, (l,i,j))[c] for c in f_connections]
        ci_fr = [var(pm, nw, :ci, (l,i,j))[c] for c in f_connections]

        vr_to = [var(pm, nw, :vr, j)[c] for c in t_connections]
        vi_to = [var(pm, nw, :vi, j)[c] for c in t_connections]
        cr_to = [var(pm, nw, :cr, (l,j,i))[c] for c in t_connections]
        ci_to = [var(pm, nw, :ci, (l,j,i))[c] for c in t_connections]
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

    for (l,i,j) in ref(pm, nw, :arcs_from_trans)
        f_connections = ref(pm, nw, :transformer, l, "f_connections")
        t_connections = ref(pm, nw, :transformer, l, "t_connections")

        vr_fr = [var(pm, nw, :vr, i)[c] for c in f_connections]
        vi_fr = [var(pm, nw, :vi, i)[c] for c in f_connections]
        cr_fr = [var(pm, nw, :crt, (l,i,j))[c] for c in f_connections]
        ci_fr = [var(pm, nw, :cit, (l,i,j))[c] for c in f_connections]

        vr_to = [var(pm, nw, :vr, j)[c] for c in t_connections]
        vi_to = [var(pm, nw, :vi, j)[c] for c in t_connections]
        cr_to = [var(pm, nw, :crt, (l,j,i))[c] for c in t_connections]
        ci_to = [var(pm, nw, :cit, (l,j,i))[c] for c in t_connections]

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
function variable_mc_load_current(pm::_PM.AbstractIVRModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true, kwargs...)
    var(pm, nw)[:crd] = Dict{Int, Any}()
    var(pm, nw)[:cid] = Dict{Int, Any}()
    var(pm, nw)[:crd_bus] = Dict{Int, Any}()
    var(pm, nw)[:cid_bus] = Dict{Int, Any}()
end


""
function variable_mc_generator_current(pm::_PM.AbstractIVRModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_generator_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_generator_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

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

    for id in ids(pm, nw, :bus)
        busref = ref(pm, nw, :bus, id)
        terminals = busref["terminals"]
        grounded = busref["grounded"]

        ncnd = length(terminals)
        vm = haskey(busref, "vm_start") ? busref["vm_start"] : fill(0.0, ncnd)
        vm[.!grounded] .= 1.0

        # TODO how to do this more generally
        nph = 3
        va = haskey(busref, "va_start") ? busref["va_start"] : [c <= nph ? _wrap_to_pi(2 * pi / nph * (1-c)) : 0.0 for c in terminals]

        for (idx,t) in enumerate(terminals)
            vr = vm[idx]*cos(va[idx])
            vi = vm[idx]*sin(va[idx])
            JuMP.set_start_value(var(pm, nw, :vr, id)[t], vr)
            JuMP.set_start_value(var(pm, nw, :vi, id)[t], vi)
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
function constraint_mc_current_from(pm::_PM.AbstractIVRModel, nw::Int, f_bus::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real})
    vr_fr = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]

    csr_fr =  [var(pm, nw, :csr, f_idx[1])[c] for c in f_connections]
    csi_fr =  [var(pm, nw, :csi, f_idx[1])[c] for c in f_connections]

    cr_fr =  [var(pm, nw, :cr, f_idx)[c] for c in f_connections]
    ci_fr =  [var(pm, nw, :ci, f_idx)[c] for c in f_connections]

    JuMP.@constraint(pm.model, cr_fr .== (csr_fr + g_sh_fr*vr_fr - b_sh_fr*vi_fr))
    JuMP.@constraint(pm.model, ci_fr .== (csi_fr + g_sh_fr*vi_fr + b_sh_fr*vr_fr))
end


"Defines how current distributes over series and shunt impedances of a pi-model branch"
function constraint_mc_current_to(pm::_PM.AbstractIVRModel, n::Int, t_bus, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, g_sh_to::Matrix{<:Real}, b_sh_to::Matrix{<:Real})
    vr_to = [var(pm, n, :vr, t_bus)[c] for c in t_connections]
    vi_to = [var(pm, n, :vi, t_bus)[c] for c in t_connections]

    csr_to = [-var(pm, n, :csr, f_idx[1])[c] for c in f_connections]
    csi_to = [-var(pm, n, :csi, f_idx[1])[c] for c in f_connections]

    cr_to = [var(pm, n, :cr, t_idx)[c] for c in t_connections]
    ci_to = [var(pm, n, :ci, t_idx)[c] for c in t_connections]

    JuMP.@constraint(pm.model, cr_to .== csr_to + g_sh_to*vr_to - b_sh_to*vi_to)
    JuMP.@constraint(pm.model, ci_to .== csi_to + g_sh_to*vi_to + b_sh_to*vr_to)
end


"Defines voltage drop over a branch, linking from and to side complex voltage"
function constraint_mc_bus_voltage_drop(pm::_PM.AbstractIVRModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, r::Matrix{<:Real}, x::Matrix{<:Real})
    vr_fr = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]

    vr_to = [var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_to = [var(pm, nw, :vi, t_bus)[c] for c in t_connections]

    csr_fr = [var(pm, nw, :csr, f_idx[1])[c] for c in f_connections]
    csi_fr = [var(pm, nw, :csi, f_idx[1])[c] for c in f_connections]

    r = r
    x = x

    JuMP.@constraint(pm.model, vr_to .== vr_fr - r*csr_fr + x*csi_fr)
    JuMP.@constraint(pm.model, vi_to .== vi_fr - r*csi_fr - x*csr_fr)
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
function constraint_mc_current_balance(pm::_PM.AbstractIVRModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    cr    = get(var(pm, nw),    :cr, Dict()); _PM._check_var_keys(cr, bus_arcs, "real current", "branch")
    ci    = get(var(pm, nw),    :ci, Dict()); _PM._check_var_keys(ci, bus_arcs, "imaginary current", "branch")
    crd   = get(var(pm, nw),   :crd_bus, Dict()); _PM._check_var_keys(crd, bus_loads, "real current", "load")
    cid   = get(var(pm, nw),   :cid_bus, Dict()); _PM._check_var_keys(cid, bus_loads, "imaginary current", "load")
    crg   = get(var(pm, nw),   :crg_bus, Dict()); _PM._check_var_keys(crg, bus_gens, "real current", "generator")
    cig   = get(var(pm, nw),   :cig_bus, Dict()); _PM._check_var_keys(cig, bus_gens, "imaginary current", "generator")
    crs   = get(var(pm, nw),   :crs, Dict()); _PM._check_var_keys(crs, bus_storage, "real currentr", "storage")
    cis   = get(var(pm, nw),   :cis, Dict()); _PM._check_var_keys(cis, bus_storage, "imaginary current", "storage")
    crsw  = get(var(pm, nw),  :crsw, Dict()); _PM._check_var_keys(crsw, bus_arcs_sw, "real current", "switch")
    cisw  = get(var(pm, nw),  :cisw, Dict()); _PM._check_var_keys(cisw, bus_arcs_sw, "imaginary current", "switch")
    crt   = get(var(pm, nw),   :crt, Dict()); _PM._check_var_keys(crt, bus_arcs_trans, "real current", "transformer")
    cit   = get(var(pm, nw),   :cit, Dict()); _PM._check_var_keys(cit, bus_arcs_trans, "imaginary current", "transformer")

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx, t) in ungrounded_terminals
        JuMP.@NLconstraint(pm.model,  sum(cr[a][t] for (a, conns) in bus_arcs if t in conns)
                                    + sum(crsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                    + sum(crt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                    ==
                                      sum(crg[g][t]         for (g, conns) in bus_gens if t in conns)
                                    - sum(crs[s][t]         for (s, conns) in bus_storage if t in conns)
                                    - sum(crd[d][t]         for (d, conns) in bus_loads if t in conns)
                                    - sum( Gt[idx,jdx]*vr[u] -Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals) # shunts
                                    )
        JuMP.@NLconstraint(pm.model, sum(ci[a][t] for (a, conns) in bus_arcs if t in conns)
                                    + sum(cisw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                    + sum(cit[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                    ==
                                      sum(cig[g][t]         for (g, conns) in bus_gens if t in conns)
                                    - sum(cis[s][t]         for (s, conns) in bus_storage if t in conns)
                                    - sum(cid[d][t]         for (d, conns) in bus_loads if t in conns)
                                    - sum( Gt[idx,jdx]*vi[u] +Bt[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals) # shunts
                                    )
    end
end


"`p[f_idx]^2 + q[f_idx]^2 <= rate_a^2`"
function constraint_mc_thermal_limit_from(pm::_PM.AbstractIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})
    (l, f_bus, t_bus) = f_idx

    vr = var(pm, nw, :vr, f_bus)
    vi = var(pm, nw, :vi, f_bus)
    crf = var(pm, nw, :cr, f_idx)
    cif = var(pm, nw, :ci, f_idx)

    for (idx, fc) in enumerate(f_connections)
        JuMP.@NLconstraint(pm.model, (vr[fc]^2 + vi[fc]^2)*(crf[fc]^2 + cif[fc]^2) <= rate_a[idx]^2)
    end
end


"`p[t_idx]^2 + q[t_idx]^2 <= rate_a^2`"
function constraint_mc_thermal_limit_to(pm::_PM.AbstractIVRModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})
    (l, t_bus, f_bus) = t_idx

    vr = var(pm, nw, :vr, t_bus)
    vi = var(pm, nw, :vi, t_bus)
    crt = var(pm, nw, :cr, t_idx)
    cit = var(pm, nw, :ci, t_idx)

    for (idx, tc) in enumerate(t_connections)
        JuMP.@NLconstraint(pm.model, (vr[tc]^2 + vi[tc]^2)*(crt[tc]^2 + cit[tc]^2) <= rate_a[idx]^2)
    end
end


"""
Bounds the current magnitude at both from and to side of a branch
`cr[f_idx]^2 + ci[f_idx]^2 <= c_rating_a^2`
`cr[t_idx]^2 + ci[t_idx]^2 <= c_rating_a^2`
"""
function constraint_mc_current_limit(pm::_PM.AbstractIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, c_rating_a::Vector{<:Real})
    (l, f_bus, t_bus) = f_idx
    t_idx = (l, t_bus, f_bus)

    crf =  [var(pm, nw, :cr, f_idx)[c] for c in f_connections]
    cif =  [var(pm, nw, :ci, f_idx)[c] for c in f_connections]

    crt =  [var(pm, nw, :cr, t_idx)[c] for c in t_connections]
    cit =  [var(pm, nw, :ci, t_idx)[c] for c in t_connections]

    JuMP.@constraint(pm.model, crf.^2 + cif.^2 .<= c_rating_a.^2)
    JuMP.@constraint(pm.model, crt.^2 + cit.^2 .<= c_rating_a.^2)
end


"""
`pmin <= Re(v*cg') <= pmax`
"""
function constraint_mc_gen_active_bounds(pm::_PM.AbstractIVRModel, nw::Int, i::Int, bus::Int, connections::Vector{Int}, pmax::Vector{<:Real}, pmin::Vector{<:Real})
    @assert all(pmin .<= pmax)

    vr = [var(pm, nw, :vr, bus)[c] for c in connections]
    vi = [var(pm, nw, :vi, bus)[c] for c in connections]
    cr = [var(pm, nw, :crg, i)[c] for c in connections]
    ci = [var(pm, nw, :cig, i)[c] for c in connections]

    JuMP.@constraint(pm.model, pmin .<= vr.*cr  + vi.*ci)
    JuMP.@constraint(pm.model, pmax .>= vr.*cr  + vi.*ci)
end


"""
`qmin <= Im(v*cg') <= qmax`
"""
function constraint_mc_gen_reactive_bounds(pm::_PM.AbstractIVRModel, nw::Int, i::Int, bus::Int, connections::Vector{Int}, qmax::Vector{<:Real}, qmin::Vector{<:Real})
    @assert all(qmin .<= qmax)

    vr = [var(pm, nw, :vr, bus)[c] for c in connections]
    vi = [var(pm, nw, :vi, bus)[c] for c in connections]
    cr = [var(pm, nw, :crg, i)[c] for c in connections]
    ci = [var(pm, nw, :cig, i)[c] for c in connections]

    JuMP.@constraint(pm.model, qmin .<= vi.*cr  - vr.*ci)
    JuMP.@constraint(pm.model, qmax .>= vi.*cr  - vr.*ci)
end


"`pg[i] == pg`"
function constraint_mc_gen_power_setpoint_real(pm::_PM.AbstractIVRModel, nw::Int, i::Int, pg_ref::Vector{<:Real})
    gen = ref(pm, nw, :gen, i)
    bus = gen["gen_bus"]
    connections = gen["connections"]
    vr = [var(pm, nw, :vr, bus)[c] for c in connections]
    vi = [var(pm, nw, :vi, bus)[c] for c in connections]
    cr = [var(pm, nw, :crg, i)[c] for c in connections]
    ci = [var(pm, nw, :cig, i)[c] for c in connections]

    JuMP.@constraint(pm.model, pg_ref .== vr.*cr  + vi.*ci)
end


"`qq[i] == qq`"
function constraint_mc_regen_setpoint_active(pm::_PM.AbstractIVRModel, n::Int, i, qg_ref)
    gen = ref(pm, n, :gen, i)
    bus = gen["gen_bus"]
    connections = gen["connections"]
    vr = [var(pm, n, :vr, bus)[c] for c in connections]
    vi = [var(pm, n, :vi, bus)[c] for c in connections]
    cr = [var(pm, n, :crg, i)[c] for c in connections]
    ci = [var(pm, n, :cig, i)[c] for c in connections]

    JuMP.@constraint(pm.model, qg_ref .== vi.*cr  - vr.*ci)
end


"wye-wye transformer power constraint for IVR formulation"
function constraint_mc_transformer_power_yy(pm::_PM.AbstractIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    vr_fr_P = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr_P = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]
    vr_fr_n = 0
    vi_fr_n = 0
    vr_to_P = [var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_to_P = [var(pm, nw, :vi, t_bus)[c] for c in t_connections]
    vr_to_n = 0
    vi_to_n = 0

    cr_fr_P = [var(pm, nw, :crt, f_idx)[c] for c in f_connections]
    ci_fr_P = [var(pm, nw, :cit, f_idx)[c] for c in f_connections]
    cr_to_P = [var(pm, nw, :crt, t_idx)[c] for c in t_connections]
    ci_to_P = [var(pm, nw, :cit, t_idx)[c] for c in t_connections]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))]
    scale = (tm_scale*pol).*tm_set

    JuMP.@constraint(pm.model, (vr_fr_P.-vr_fr_n) .== scale.*(vr_to_P.-vr_to_n))
    JuMP.@constraint(pm.model, (vi_fr_P.-vi_fr_n) .== scale.*(vi_to_P.-vi_to_n))
    JuMP.@constraint(pm.model, scale.*cr_fr_P .+ cr_to_P .== 0)
    JuMP.@constraint(pm.model, scale.*ci_fr_P .+ ci_to_P .== 0)
end


"delta-wye transformer power constraint for IVR formulation"
function constraint_mc_transformer_power_dy(pm::_PM.AbstractIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    vr_fr_P = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr_P = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]
    vr_to_P = [var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_to_P = [var(pm, nw, :vi, t_bus)[c] for c in t_connections]
    vr_to_n = 0
    vi_to_n = 0

    cr_fr_P = [var(pm, nw, :crt, f_idx)[c] for c in f_connections]
    ci_fr_P = [var(pm, nw, :cit, f_idx)[c] for c in f_connections]
    cr_to_P = [var(pm, nw, :crt, t_idx)[c] for c in t_connections]
    ci_to_P = [var(pm, nw, :cit, t_idx)[c] for c in t_connections]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))]
    scale = (tm_scale*pol).*tm_set

    n_phases = length(tm)
    Md = _get_delta_transformation_matrix(n_phases)

    JuMP.@constraint(pm.model, Md*vr_fr_P .== scale.*(vr_to_P .- vr_to_n))
    JuMP.@constraint(pm.model, Md*vi_fr_P .== scale.*(vi_to_P .- vi_to_n))

    JuMP.@constraint(pm.model, scale.*cr_fr_P .+ Md'*cr_to_P .== 0)
    JuMP.@constraint(pm.model, scale.*ci_fr_P .+ Md'*ci_to_P .== 0)
end


"wye connected load setpoint constraint for IVR formulation"
function constraint_mc_load_power_wye(pm::_PM.IVRPowerModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    crd = Vector{JuMP.NonlinearExpression}([])
    cid = Vector{JuMP.NonlinearExpression}([])

    for (idx, c) in enumerate(connections)
        push!(crd, JuMP.@NLexpression(pm.model,
             a[idx]*vr[c]*(vr[c]^2+vi[c]^2)^(alpha[idx]/2-1)
            +b[idx]*vi[c]*(vr[c]^2+vi[c]^2)^(beta[idx]/2 -1)
        ))
        push!(cid, JuMP.@NLexpression(pm.model,
             a[idx]*vi[c]*(vr[c]^2+vi[c]^2)^(alpha[idx]/2-1)
            -b[idx]*vr[c]*(vr[c]^2+vi[c]^2)^(beta[idx]/2 -1)
        ))
    end

    var(pm, nw, :crd_bus)[id] = JuMP.Containers.DenseAxisArray(crd, connections)
    var(pm, nw, :cid_bus)[id] = JuMP.Containers.DenseAxisArray(cid, connections)

    if report
        pd_bus = Vector{JuMP.NonlinearExpression}([])
        qd_bus = Vector{JuMP.NonlinearExpression}([])
        for (idx,c) in enumerate(connections)
            push!(pd_bus, JuMP.@NLexpression(pm.model,  vr[c]*crd[idx]+vi[c]*cid[idx]))
            push!(qd_bus, JuMP.@NLexpression(pm.model, -vr[c]*cid[idx]+vi[c]*crd[idx]))
        end

        sol(pm, nw, :load, id)[:pd_bus] = JuMP.Containers.DenseAxisArray(pd_bus, connections)
        sol(pm, nw, :load, id)[:qd_bus] = JuMP.Containers.DenseAxisArray(qd_bus, connections)

        sol(pm, nw, :load, id)[:crd_bus] = JuMP.Containers.DenseAxisArray(crd, connections)
        sol(pm, nw, :load, id)[:cid_bus] = JuMP.Containers.DenseAxisArray(cid, connections)

        pd = Vector{JuMP.NonlinearExpression}([])
        qd = Vector{JuMP.NonlinearExpression}([])
        for (idx, c) in enumerate(connections)
            push!(pd, JuMP.@NLexpression(pm.model, a[idx]*(vr[c]^2+vi[c]^2)^(alpha[idx]/2) ))
            push!(qd, JuMP.@NLexpression(pm.model, b[idx]*(vr[c]^2+vi[c]^2)^(beta[idx]/2)  ))
        end
        sol(pm, nw, :load, id)[:pd] = JuMP.Containers.DenseAxisArray(pd, connections)
        sol(pm, nw, :load, id)[:qd] = JuMP.Containers.DenseAxisArray(qd, connections)
    end
end


"delta connected load setpoint constraint for IVR formulation"
function constraint_mc_load_power_delta(pm::_PM.IVRPowerModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    nph = length(connections)
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
function constraint_mc_generator_power_wye(pm::_PM.IVRPowerModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)

    pg = Vector{JuMP.NonlinearExpression}([])
    qg = Vector{JuMP.NonlinearExpression}([])

    for (idx, c) in enumerate(connections)
        push!(pg, JuMP.@NLexpression(pm.model,  vr[c]*crg[c]+vi[c]*cig[c]))
        push!(qg, JuMP.@NLexpression(pm.model, -vr[c]*cig[c]+vi[c]*crg[c]))
    end

    if bounded
        for (idx,c) in enumerate(connections)
            if pmin[idx]>-Inf
                JuMP.@constraint(pm.model, pmin[idx] .<= vr[c]*crg[c]  + vi[c]*cig[c])
            end
            if pmax[idx]< Inf
                JuMP.@constraint(pm.model, pmax[idx] .>= vr[c]*crg[c]  + vi[c]*cig[c])
            end
            if qmin[idx]>-Inf
                JuMP.@constraint(pm.model, qmin[idx] .<= vi[c]*crg[c]  - vr[c]*cig[c])
            end
            if qmax[idx]< Inf
                JuMP.@constraint(pm.model, qmax[idx] .>= vi[c]*crg[c]  - vr[c]*cig[c])
            end
        end
    end

    var(pm, nw, :crg_bus)[id] = crg
    var(pm, nw, :cig_bus)[id] = cig
    var(pm, nw, :pg)[id] = JuMP.Containers.DenseAxisArray(pg, connections)
    var(pm, nw, :qg)[id] = JuMP.Containers.DenseAxisArray(qg, connections)

    if report
        sol(pm, nw, :gen, id)[:crg_bus] = var(pm, nw, :crg_bus, id)
        sol(pm, nw, :gen, id)[:cig_bus] = var(pm, nw, :cig_bus, id)

        sol(pm, nw, :gen, id)[:pg] = JuMP.Containers.DenseAxisArray(pg, connections)
        sol(pm, nw, :gen, id)[:qg] = JuMP.Containers.DenseAxisArray(qg, connections)
    end
end


"delta connected generator setpoint constraint for IVR formulation"
function constraint_mc_generator_power_delta(pm::_PM.IVRPowerModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)

    nph = length(pmin)

    prev = Dict(c=>connections[(idx+nph-2)%nph+1] for (idx,c) in enumerate(connections))
    next = Dict(c=>connections[idx%nph+1] for (idx,c) in enumerate(connections))

    vrg = Dict()
    vig = Dict()
    for c in connections
        vrg[c] = JuMP.@NLexpression(pm.model, vr[c]-vr[next[c]])
        vig[c] = JuMP.@NLexpression(pm.model, vi[c]-vi[next[c]])
    end

    pg = Vector{JuMP.NonlinearExpression}([])
    qg = Vector{JuMP.NonlinearExpression}([])
    for c in connections
        push!(pg, JuMP.@NLexpression(pm.model,  vrg[c]*crg[c]+vig[c]*cig[c]))
        push!(qg, JuMP.@NLexpression(pm.model, -vrg[c]*cig[c]+vig[c]*crg[c]))
    end

    if bounded
        JuMP.@NLconstraint(pm.model, [i in 1:nph], pmin[i] <= pg[i])
        JuMP.@NLconstraint(pm.model, [i in 1:nph], pmax[i] >= pg[i])
        JuMP.@NLconstraint(pm.model, [i in 1:nph], qmin[i] <= qg[i])
        JuMP.@NLconstraint(pm.model, [i in 1:nph], qmax[i] >= qg[i])
    end

    crg_bus = Vector{JuMP.NonlinearExpression}([])
    cig_bus = Vector{JuMP.NonlinearExpression}([])
    for c in connections
        push!(crg_bus, JuMP.@NLexpression(pm.model, crg[c]-crg[prev[c]]))
        push!(cig_bus, JuMP.@NLexpression(pm.model, cig[c]-cig[prev[c]]))
    end

    var(pm, nw, :crg_bus)[id] = JuMP.Containers.DenseAxisArray(crg_bus, connections)
    var(pm, nw, :cig_bus)[id] = JuMP.Containers.DenseAxisArray(cig_bus, connections)
    var(pm, nw, :pg)[id] = JuMP.Containers.DenseAxisArray(pg, connections)
    var(pm, nw, :qg)[id] = JuMP.Containers.DenseAxisArray(qg, connections)

    if report
        sol(pm, nw, :gen, id)[:crg_bus] = JuMP.Containers.DenseAxisArray(crg_bus, connections)
        sol(pm, nw, :gen, id)[:cig_bus] = JuMP.Containers.DenseAxisArray(cig_bus, connections)
        sol(pm, nw, :gen, id)[:pg] = JuMP.Containers.DenseAxisArray(pg, connections)
        sol(pm, nw, :gen, id)[:qg] = JuMP.Containers.DenseAxisArray(qg, connections)
    end
end
