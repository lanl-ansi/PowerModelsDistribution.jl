# BUS

# BUS - Variables

"""
	function variable_mc_bus_voltage(
		pm::RectangularVoltageExplicitNeutralModels;
		nw=nw_id_default,
		bounded::Bool=true,
	)

Creates rectangular voltage variables `:vr` and `:vi` for models with explicit neutrals
"""
function variable_mc_bus_voltage(pm::RectangularVoltageExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_bus_voltage_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_bus_voltage_imaginary(pm; nw=nw, bounded=bounded, report=report)
end


"""
	function variable_mc_bus_voltage_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates real voltage variables `:vr` for models with explicit neutrals
"""
function variable_mc_bus_voltage_real(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref(pm, nw, :bus))

    vr = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_vr_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "vr_start", t, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in ref(pm, nw, :bus)
            for (idx,t) in enumerate(terminals[i])
                set_lower_bound.(vr[i][t], -bus["vmax"][idx])
                set_upper_bound.(vr[i][t],  bus["vmax"][idx])
            end
        end
    end

    # perfectly grounded terminals get a constant 0 instead of a variable
    vr = var(pm, nw)[:vr] = Dict(i=>JuMP.Containers.DenseAxisArray(
        Vector{JuMP.AffExpr}([t in vr[i].axes[1] ? vr[i][t] : 0.0 for t in bus["terminals"]]), bus["terminals"]
        ) for (i, bus) in ref(pm, nw, :bus))

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :vr, ids(pm, nw, :bus), vr)
end


"""
	function variable_mc_bus_voltage_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates imaginary voltage variables `:vr` for models with explicit neutrals
"""
function variable_mc_bus_voltage_imaginary(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i,bus) in ref(pm, nw, :bus))
    vi = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_vi_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "vi_start", t, 0.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in ref(pm, nw, :bus)
            for (idx,t) in enumerate(terminals[i])
                set_lower_bound.(vi[i][t], -bus["vmax"][idx])
                set_upper_bound.(vi[i][t],  bus["vmax"][idx])
            end
        end
    end

    # perfectly grounded terminals get a constant 0 instead of a variable
    vi = var(pm, nw)[:vi] = Dict(i=>JuMP.Containers.DenseAxisArray(
        Vector{JuMP.AffExpr}([t in vi[i].axes[1] ? vi[i][t] : 0.0 for t in bus["terminals"]]), bus["terminals"]
        ) for (i, bus) in ref(pm, nw, :bus))

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :vi, ids(pm, nw, :bus), vi)
end


# Constraints

"""
	function constraint_mc_current_balance(
		pm::RectangularVoltageExplicitNeutralModels,
		nw::Int,
		i::Int,
		terminals::Vector{Int},
		grounded::Vector{Bool},
		bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
		bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
		bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
		bus_gens::Vector{Tuple{Int,Vector{Int}}},
		bus_storage::Vector{Tuple{Int,Vector{Int}}},
		bus_loads::Vector{Tuple{Int,Vector{Int}}},
		bus_shunts::Vector{Tuple{Int,Vector{Int}}}
	)

Kirchhoff's current law applied to buses
`sum(cr + im*ci) = 0`
"""
function constraint_mc_current_balance(pm::RectangularVoltageExplicitNeutralModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    cr    = get(var(pm, nw), :cr_bus,   Dict()); _check_var_keys(cr, bus_arcs, "real current", "branch")
    ci    = get(var(pm, nw), :ci_bus,   Dict()); _check_var_keys(ci, bus_arcs, "imaginary current", "branch")
    crd   = get(var(pm, nw), :crd_bus,  Dict()); _check_var_keys(crd, bus_loads, "real current", "load")
    cid   = get(var(pm, nw), :cid_bus,  Dict()); _check_var_keys(cid, bus_loads, "imaginary current", "load")
    crg   = get(var(pm, nw), :crg_bus,  Dict()); _check_var_keys(crg, bus_gens, "real current", "generator")
    cig   = get(var(pm, nw), :cig_bus,  Dict()); _check_var_keys(cig, bus_gens, "imaginary current", "generator")
    crs   = get(var(pm, nw), :crs_bus,  Dict()); _check_var_keys(crs, bus_storage, "real currentr", "storage")
    cis   = get(var(pm, nw), :cis_bus,  Dict()); _check_var_keys(cis, bus_storage, "imaginary current", "storage")
    crsw  = get(var(pm, nw), :crsw_bus, Dict()); _check_var_keys(crsw, bus_arcs_sw, "real current", "switch")
    cisw  = get(var(pm, nw), :cisw_bus, Dict()); _check_var_keys(cisw, bus_arcs_sw, "imaginary current", "switch")
    crt   = get(var(pm, nw), :crt_bus,  Dict()); _check_var_keys(crt, bus_arcs_trans, "real current", "transformer")
    cit   = get(var(pm, nw), :cit_bus,  Dict()); _check_var_keys(cit, bus_arcs_trans, "imaginary current", "transformer")

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx, t) in ungrounded_terminals
        JuMP.@constraint(pm.model,
                                      sum(cr[a][t] for (a, conns) in bus_arcs if t in conns)
                                    + sum(crsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                    + sum(crt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                    ==
                                      sum(crg[g][t]         for (g, conns) in bus_gens if t in conns)
                                    - sum(crs[s][t]         for (s, conns) in bus_storage if t in conns)
                                    - sum(crd[d][t]         for (d, conns) in bus_loads if t in conns)
                                    - sum( Gt[idx,jdx]*vr[u] -Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals) # shunts
                                    )
        JuMP.@constraint(pm.model,
                                      sum(ci[a][t] for (a, conns) in bus_arcs if t in conns)
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


"""
	function constraint_mc_power_balance(
		pm::RectangularVoltageExplicitNeutralModels,
		nw::Int,
		i::Int,
		terminals::Vector{Int},
		grounded::Vector{Bool},
		bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
		bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
		bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
		bus_gens::Vector{Tuple{Int,Vector{Int}}},
		bus_storage::Vector{Tuple{Int,Vector{Int}}},
		bus_loads::Vector{Tuple{Int,Vector{Int}}},
		bus_shunts::Vector{Tuple{Int,Vector{Int}}}
	)

Imposes power balance constraints at each ungrounded terminal of bus `i` for rectangular voltage models with explicit neutrals.
`sum(p + im*q) = 0`
"""
function constraint_mc_power_balance(pm::RectangularVoltageExplicitNeutralModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)
    p    = get(var(pm, nw), :p_bus,   Dict()); _check_var_keys(p,   bus_arcs,       "active power",   "branch")
    q    = get(var(pm, nw), :q_bus,   Dict()); _check_var_keys(q,   bus_arcs,       "reactive power", "branch")
    pg   = get(var(pm, nw), :pg_bus,  Dict()); _check_var_keys(pg,  bus_gens,       "active power",   "generator")
    qg   = get(var(pm, nw), :qg_bus,  Dict()); _check_var_keys(qg,  bus_gens,       "reactive power", "generator")
    ps   = get(var(pm, nw), :ps_bus,  Dict()); _check_var_keys(ps,  bus_storage,    "active power",   "storage")
    qs   = get(var(pm, nw), :qs_bus,  Dict()); _check_var_keys(qs,  bus_storage,    "reactive power", "storage")
    psw  = get(var(pm, nw), :psw_bus, Dict()); _check_var_keys(psw, bus_arcs_sw,    "active power",   "switch")
    qsw  = get(var(pm, nw), :qsw_bus, Dict()); _check_var_keys(qsw, bus_arcs_sw,    "reactive power", "switch")
    pt   = get(var(pm, nw), :pt_bus,  Dict()); _check_var_keys(pt,  bus_arcs_trans, "active power",   "transformer")
    qt   = get(var(pm, nw), :qt_bus,  Dict()); _check_var_keys(qt,  bus_arcs_trans, "reactive power", "transformer")
    pd   = get(var(pm, nw), :pd_bus,  Dict()); _check_var_keys(pd,  bus_loads,      "active power",   "load")
    qd   = get(var(pm, nw), :qd_bus,  Dict()); _check_var_keys(pd,  bus_loads,      "reactive power", "load")

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    # pd/qd can be NLexpressions, so cannot be vectorized
    for (idx, t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
              sum(  p[arc][t] for (arc, conns) in bus_arcs if t in conns)
            + sum(psw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
            + sum( pt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
            ==
              sum(pg[gen][t] for (gen, conns) in bus_gens if t in conns)
            - sum(ps[strg][t] for (strg, conns) in bus_storage if t in conns)
            - sum(pd[load][t] for (load, conns) in bus_loads if t in conns)
            + ( -vr[t] * sum(Gt[idx,jdx]*vr[u]-Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals)
                -vi[t] * sum(Gt[idx,jdx]*vi[u]+Bt[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals)
            )
        )
        push!(cstr_p, cp)

        cq = JuMP.@constraint(pm.model,
              sum(  q[arc][t] for (arc, conns) in bus_arcs if t in conns)
            + sum(qsw[arc][t] for (arc, conns) in bus_arcs_sw if t in conns)
            + sum( qt[arc][t] for (arc, conns) in bus_arcs_trans if t in conns)
            ==
              sum(qg[gen][t] for (gen, conns) in bus_gens if t in conns)
            - sum(qd[load][t] for (load, conns) in bus_loads if t in conns)
            - sum(qs[strg][t] for (strg, conns) in bus_storage if t in conns)
            + ( vr[t] * sum(Gt[idx,jdx]*vi[u]+Bt[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals)
               -vi[t] * sum(Gt[idx,jdx]*vr[u]-Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals)
            )
        )
        push!(cstr_q, cq)
    end

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = cstr_q

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


"""
	function constraint_mc_voltage_absolute(
		pm::RectangularVoltageExplicitNeutralModels,
		nw::Int,
		i::Int,
		terminals::Vector{Int},
		grounded::Vector{Bool},
		vmin::Vector{<:Real},
		vmax::Vector{<:Real};
		report::Bool=true
	)

Imposes absolute voltage magnitude bounds for models with explicit neutrals
"""
function constraint_mc_voltage_absolute(pm::RectangularVoltageExplicitNeutralModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, vmin::Vector{<:Real}, vmax::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    ungrounded_terminals = terminals[(!).(grounded)]
    for (idx,t) in enumerate(terminals)
        if !grounded[idx]
            if vmax[idx] < Inf
                JuMP.@constraint(pm.model, vr[t]^2+vi[t]^2 <= vmax[idx]^2)
            end
            if vmin[idx] > 0.0
                JuMP.@constraint(pm.model, vr[t]^2+vi[t]^2 >= vmin[idx]^2)
            end
        end
    end
end


"""
	function constraint_mc_voltage_pairwise(
		pm::RectangularVoltageExplicitNeutralModels,
		nw::Int,
		i::Int,
		terminals::Vector{Int},
		grounded::Vector{Bool},
		vm_pair_lb::Vector,
		vm_pair_ub::Vector;
		report::Bool=true
	)

Imposes pairwise voltage magnitude bounds, i.e. magnitude bounds on the voltage between to terminals, for models with explicit neutrals
"""
function constraint_mc_voltage_pairwise(pm::RectangularVoltageExplicitNeutralModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, vm_pair_lb::Vector{<:Tuple{Any,Any,Real}}, vm_pair_ub::Vector{<:Tuple{Any,Any,Real}}; report::Bool=true)
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    for (a,b,lb) in vm_pair_lb
        if lb > 0.0
            JuMP.@constraint(pm.model, (vr[a]-vr[b])^2 + (vi[a]-vi[b])^2 >= lb^2)
        end
    end

    for (a,b,ub) in vm_pair_ub
        if ub < Inf
            JuMP.@constraint(pm.model, (vr[a]-vr[b])^2 + (vi[a]-vi[b])^2 <= ub^2)
        end
    end
end


"""
	function constraint_mc_voltage_fixed(
		pm::RectangularVoltageExplicitNeutralModels,
		nw::Int,
		i::Int,
		vm::Vector{<:Real},
		va::Vector{<:Real},
		terminals::Vector{Int},
		grounded::Vector{Bool}
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
"""
function constraint_mc_voltage_fixed(pm::RectangularVoltageExplicitNeutralModels, nw::Int, i::Int, vm::Vector{<:Real}, va::Vector{<:Real}, terminals::Vector{Int}, grounded::Vector{Bool})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    idxs = findall((!).(grounded))

    JuMP.@constraint(pm.model, [i in idxs], vr[terminals[i]]==vm[i]*cos(va[i]))
    JuMP.@constraint(pm.model, [i in idxs], vi[terminals[i]]==vm[i]*sin(va[i]))
end


"""
	function constraint_mc_voltage_magnitude_fixed(
		pm::RectangularVoltageExplicitNeutralModels,
		nw::Int,
		i::Int,
		vm::Vector{<:Real},
		va::Vector{<:Real},
		terminals::Vector{Int},
		grounded::Vector{Bool}
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
"""
function constraint_mc_voltage_magnitude_fixed(pm::RectangularVoltageExplicitNeutralModels, nw::Int, i::Int, vm::Vector{<:Real}, terminals::Vector{Int}, grounded::Vector{Bool})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    @assert iszero(vm[grounded]) "Infeasible model; the voltage magnitude of a grounded terminal is fixed to a non-zero value."
    idxs = findall((!).(grounded))

    JuMP.@constraint(pm.model, [i in idxs], vr[terminals[i]].^2 + vi[terminals[i]].^2 == vm[i].^2)
end


"""
	function constraint_mc_theta_ref(
		pm::RectangularVoltageExplicitNeutralModels,
		nw::Int,
		i::Int,
		va::Vector{<:Real},
		terminals::Vector{Int},
		grounded::Vector{Bool}
	)

Creates phase angle constraints at bus `i`
"""
function constraint_mc_theta_ref(pm::RectangularVoltageExplicitNeutralModels, nw::Int, i::Int, va::Vector{<:Real}, terminals::Vector{Int}, grounded::Vector{Bool})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    # deal with cases first where tan(theta)==Inf or tan(theta)==0
    for idx in findall[(!).(grounded)]
        t = terminals[idx]
        if va[t] == pi/2
            JuMP.@constraint(pm.model, vr[t] == 0)
            JuMP.@constraint(pm.model, vi[t] >= 0)
        elseif va[t] == -pi/2
            JuMP.@constraint(pm.model, vr[t] == 0)
            JuMP.@constraint(pm.model, vi[t] <= 0)
        elseif va[t] == 0
            JuMP.@constraint(pm.model, vr[t] >= 0)
            JuMP.@constraint(pm.model, vi[t] == 0)
        elseif va[t] == pi
            JuMP.@constraint(pm.model, vr[t] >= 0)
            JuMP.@constraint(pm.model, vi[t] == 0)
        else
            JuMP.@constraint(pm.model, vi[t] == tan(va[t])*vr[t])
            # va_ref also implies a sign for vr, vi
            if 0<=va[t] && va[t] <= pi
                JuMP.@constraint(pm.model, vi[t] >= 0)
            else
                JuMP.@constraint(pm.model, vi[t] <= 0)
            end
        end
    end
end


# GENERATOR

# GENERATOR - Variables

"""
	function variable_mc_generator_current_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates generator real current variables `:crg` for models with explicit neutrals
"""
function variable_mc_generator_current_real(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(gen, false) for (i,gen) in ref(pm, nw, :gen))
    crg = var(pm, nw)[:crg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_crg_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "crg_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :crg, ids(pm, nw, :gen), crg)
end


"""
	function variable_mc_generator_current_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates generator imaginary current variables `:cig` for models with explicit neutrals
"""
function variable_mc_generator_current_imaginary(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(gen, false) for (i,gen) in ref(pm, nw, :gen))
    cig = var(pm, nw)[:cig] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_cig_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "cig_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :cig, ids(pm, nw, :gen), cig)
end


"""
	function variable_mc_generator_power_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates generator active power variables `:pg` for models with explicit neutrals
"""
function variable_mc_generator_power_real(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(gen, false) for (i,gen) in ref(pm, nw, :gen))
    pg = var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_pg_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "pg_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )

    if bounded
        for (i,gen) in ref(pm, nw, :gen)
            set_lower_bound.(pg[i], gen["pmin"])
            set_upper_bound.(pg[i], gen["pmax"])
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :pg, ids(pm, nw, :gen), pg)
end


"""
	function variable_mc_generator_power_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates generator reactive power variables `:qg` for models with explicit neutrals
"""
function variable_mc_generator_power_imaginary(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(gen, false) for (i,gen) in ref(pm, nw, :gen))
    qg = var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_qg_$(i)",
            start = comp_start_value(ref(pm, nw, :gen, i), "qg_start", c, 0.0)
        ) for i in ids(pm, nw, :gen)
    )

    if bounded
        for (i,gen) in ref(pm, nw, :gen)
            set_lower_bound.(qg[i], gen["qmin"])
            set_upper_bound.(qg[i], gen["qmax"])
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :qg, ids(pm, nw, :gen), qg)
end


# LOAD

# LOAD - Variables

"""
	function variable_mc_load_current_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates load real current variables `:crd` for models with explicit neutrals
"""
function variable_mc_load_current_real(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(load, false) for (i,load) in ref(pm, nw, :load))

    crd = var(pm, nw)[:crd] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_crd_$(i)",
            start = comp_start_value(ref(pm, nw, :load, i), "crd_start", c, 0.0)
        ) for i in ids(pm, nw, :load)
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :crd, ids(pm, nw, :load), crd)
end


"""
	function variable_mc_load_current_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates load imaginary current variables `:cid` for models with explicit neutrals
"""
function variable_mc_load_current_imaginary(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(load, false) for (i,load) in ref(pm, nw, :load))
    load_ids_current = [id for (id,load) in ref(pm, nw, :load) if load["model"]==CURRENT]

    cid = var(pm, nw)[:cid] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_cid_$(i)",
            start = comp_start_value(ref(pm, nw, :load, i), "cid_start", c, 0.0)
        ) for i in ids(pm, nw, :load)
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :cid, ids(pm, nw, :load), cid)
end


# TRANSFORMER

# TRANSFORMER - Variables

"""
	function variable_mc_transformer_current_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates transformer real current variables `:crt` for models with explicit neutrals
"""
function variable_mc_transformer_current_real(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in ref(pm, nw, :transformer))
    crt = var(pm, nw)[:crt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:int_dim[l]], base_name="$(nw)_crt_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :transformer, l), "crt_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
    )

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :cr_fr, :cr_to, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), crt)
end


"""
	function variable_mc_transformer_current_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates transformer imaginary current variables `:cit` for models with explicit neutrals
"""
function variable_mc_transformer_current_imaginary(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in ref(pm, nw, :transformer))
    cit = var(pm, nw)[:cit] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:int_dim[l]], base_name="$(nw)_cit_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :transformer, l), "cit_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
    )

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :ci_fr, :ci_to, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), cit)
end


"""
	function variable_mc_transformer_power_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates transformer active power variables `:pt` for models with explicit neutrals
"""
function variable_mc_transformer_power_real(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in ref(pm, nw, :transformer))
    pt = var(pm, nw)[:pt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:int_dim[l]], base_name="$(nw)_pt_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :transformer, l), "pt_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_transformer_from)
            trans = ref(pm, nw, :transformer, l)
            f_bus = ref(pm, nw, :bus, i)
            t_bus = ref(pm, nw, :bus, j)
            sm_ub = trans["sm_ub"]
            set_lower_bound(pt[(l,i,j)], -sm_ub)
            set_upper_bound(pt[(l,i,j)],  sm_ub)
            set_lower_bound(pt[(l,j,i)], -sm_ub)
            set_upper_bound(pt[(l,j,i)],  sm_ub)
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :pf, :pt, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), pt)
end


"""
	function variable_mc_transformer_power_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates transformer reactive power variables `:qt` for models with explicit neutrals
"""
function variable_mc_transformer_power_imaginary(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in ref(pm, nw, :transformer))
    qt = var(pm, nw)[:qt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:int_dim[l]], base_name="$(nw)_qt_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :transformer, l), "qt_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_transformer)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_transformer_from)
            trans = ref(pm, nw, :transformer, l)
            f_bus = ref(pm, nw, :bus, i)
            t_bus = ref(pm, nw, :bus, j)
            sm_ub = trans["sm_ub"]
            set_lower_bound(qt[(l,i,j)], -sm_ub)
            set_upper_bound(qt[(l,i,j)],  sm_ub)
            set_lower_bound(qt[(l,j,i)], -sm_ub)
            set_upper_bound(qt[(l,j,i)],  sm_ub)
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer, :qf, :qt, ref(pm, nw, :arcs_transformer_from), ref(pm, nw, :arcs_transformer_to), qt)
end


# TRANSFORMER - Constraints


"""
	function constraint_mc_transformer_voltage_yy(
		pm::RectangularVoltageExplicitNeutralModels,
		nw::Int,
		trans_id::Int,
		f_bus::Int,
		t_bus::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int},
		pol::Int,
		tm_set::Vector{<:Real},
		tm_fixed::Vector{Bool},
		tm_scale::Real
	)

For rectangular voltage models with explicit neutrals,
links the voltage of the from-side and to-side transformer windings
for wye-wye connected transformers

```
(vr_fr_P-vr_fr_n) == scale * (vr_to_P.-vr_to_n)
(vi_fr_P-vi_fr_n) == scale * (vi_to_P.-vi_to_n)
```
"""
function constraint_mc_transformer_voltage_yy(pm::RectangularVoltageExplicitNeutralModels, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    vr_fr_P = [var(pm, nw, :vr, f_bus)[c] for c in f_connections[1:end-1]]
    vi_fr_P = [var(pm, nw, :vi, f_bus)[c] for c in f_connections[1:end-1]]
    vr_fr_n = var(pm, nw, :vr, f_bus)[f_connections[end]]
    vi_fr_n = var(pm, nw, :vi, f_bus)[f_connections[end]]
    vr_to_P = [var(pm, nw, :vr, t_bus)[c] for c in t_connections[1:end-1]]
    vi_to_P = [var(pm, nw, :vi, t_bus)[c] for c in t_connections[1:end-1]]
    vr_to_n = var(pm, nw, :vr, t_bus)[t_connections[end]]
    vi_to_n = var(pm, nw, :vi, t_bus)[t_connections[end]]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
    scale = (tm_scale*pol).*tm_set

    JuMP.@constraint(pm.model, (vr_fr_P.-vr_fr_n) .== scale.*(vr_to_P.-vr_to_n))
    JuMP.@constraint(pm.model, (vi_fr_P.-vi_fr_n) .== scale.*(vi_to_P.-vi_to_n))
end


"""
	function constraint_mc_transformer_voltage_dy(
		pm::RectangularVoltageExplicitNeutralModels,
		nw::Int,
		trans_id::Int,
		f_bus::Int,
		t_bus::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int},
		pol::Int,
		tm_set::Vector{<:Real},
		tm_fixed::Vector{Bool},
		tm_scale::Real
	)

For rectangular voltage models with explicit neutrals,
links the voltage of the from-side and to-side transformer windings
for delta-wye connected transformers

```
Md*vr_fr_P == scale * (vr_to_P - vr_to_n)
Md*vi_fr_P == scale * (vi_to_P - vi_to_n)
```
"""
function constraint_mc_transformer_voltage_dy(pm::RectangularVoltageExplicitNeutralModels, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    vr_fr_P = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr_P = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]
    vr_to_P = [var(pm, nw, :vr, t_bus)[c] for c in t_connections[1:end-1]]
    vi_to_P = [var(pm, nw, :vi, t_bus)[c] for c in t_connections[1:end-1]]
    vr_to_n = var(pm, nw, :vr, t_bus)[t_connections[end]]
    vi_to_n = var(pm, nw, :vi, t_bus)[t_connections[end]]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
    scale = (tm_scale*pol).*tm_set

    n_phases = length(tm)
    Md = _get_delta_transformation_matrix(n_phases)

    JuMP.@constraint(pm.model, Md*vr_fr_P .== scale.*(vr_to_P .- vr_to_n))
    JuMP.@constraint(pm.model, Md*vi_fr_P .== scale.*(vi_to_P .- vi_to_n))
end


# BRANCH

# BRANCH - Variables

"""
	function variable_mc_branch_current_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates branch real current variables `:cr` for models with explicit neutrals
"""
function variable_mc_branch_current_real(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref(pm, nw, :branch))
    cr = var(pm, nw)[:cr] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_cr_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "cr_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_branch)
            cmax = ref(pm, nw, :branch, l)["c_rating_a"]
            set_upper_bound.(cr[(l,i,j)],  cmax)
            set_lower_bound.(cr[(l,i,j)], -cmax)
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch, :cr_fr, :cr_to, ref(pm, nw, :arcs_branch_from), ref(pm, nw, :arcs_branch_to), cr)
end


"""
	function variable_mc_branch_current_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates branch imaginary current variables `:ci` for models with explicit neutrals
"""
function variable_mc_branch_current_imaginary(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref(pm, nw, :branch))
    ci = var(pm, nw)[:ci] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_ci_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "ci_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_branch)
            cmax = ref(pm, nw, :branch, l)["c_rating_a"]
            set_upper_bound.(ci[(l,i,j)],  cmax)
            set_lower_bound.(ci[(l,i,j)], -cmax)
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch, :ci_fr, :ci_to, ref(pm, nw, :arcs_branch_from), ref(pm, nw, :arcs_branch_to), ci)
end


"""
	function variable_mc_branch_current_series_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates branch real series current variables `:csr` for models with explicit neutrals
"""
function variable_mc_branch_current_series_real(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref(pm, nw, :branch))
    csr = var(pm, nw)[:csr] = Dict(l => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_csr_$(l)",
            start = comp_start_value(ref(pm, nw, :branch, l), "csr_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch, :csr_fr, ids(pm, nw, :branch), csr)
end


"""
	function variable_mc_branch_current_series_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates branch imaginary series current variables `:csi` for models with explicit neutrals
"""
function variable_mc_branch_current_series_imaginary(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref(pm, nw, :branch))
    csi = var(pm, nw)[:csi] = Dict(l => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_csi_$(l)",
            start = comp_start_value(ref(pm, nw, :branch, l), "csi_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch, :csi_fr, ids(pm, nw, :branch), csi)
end


"""
	function variable_mc_branch_power_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates branch active power variables `:p` for models with explicit neutrals
"""
function variable_mc_branch_power_real(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref(pm, nw, :branch))
    p = var(pm, nw)[:p] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_p_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "p_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch, :pf, :pt, ref(pm, nw, :arcs_branch_from), ref(pm, nw, :arcs_branch_to), p)
end


"""
	function variable_mc_branch_power_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates branch reactive power variables `:q` for models with explicit neutrals
"""
function variable_mc_branch_power_imaginary(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref(pm, nw, :branch))
    q = var(pm, nw)[:q] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_q_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "q_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_branch)
    )

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch, :qf, :qt, ref(pm, nw, :arcs_branch_from), ref(pm, nw, :arcs_branch_to), q)
end


# BRANCH - Constraints

"""
	function constraint_mc_branch_current_limit(
		pm::ExplicitNeutralModels,
		id::Int;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
	)

For models with explicit neutrals,
imposes a bound on the current magnitude per conductor
at both ends of the branch (total current, i.e. including shunt contributions)
"""
function constraint_mc_branch_current_limit(pm::ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
branch = ref(pm, nw, :branch, id)
f_idx = (id,branch["f_bus"],branch["t_bus"])
t_idx = (id,branch["t_bus"],branch["f_bus"])

constraint_mc_branch_current_limit(pm, nw, f_idx, t_idx, branch["f_connections"], branch["t_connections"], branch["c_rating_a"])
end


# SWITCH

# SWITCH - Variables

"""
	function variable_mc_switch_current_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For models with explicit neutrals,
creates switch real current variables `:crsw` for models with explicit neutrals.
"""
function variable_mc_switch_current_real(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(switch["f_connections"]) for (l,switch) in ref(pm, nw, :switch))
    crsw = var(pm, nw)[:crsw] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_crsw_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :switch, l), "crsw_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_switch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_switch)
            cmax = ref(pm, nw, :switch, l)["current_rating"]
            set_upper_bound.(crsw[(l,i,j)],  cmax)
            set_lower_bound.(crsw[(l,i,j)], -cmax)
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :switch, :cr_fr, :cr_to, ref(pm, nw, :arcs_switch_from), ref(pm, nw, :arcs_switch_to), crsw)
end


"""
	function variable_mc_switch_current_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For models with explicit neutrals,
creates switch imaginary current variables `:cisw` for models with explicit neutrals.
"""
function variable_mc_switch_current_imaginary(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(switch["f_connections"]) for (l,switch) in ref(pm, nw, :switch))
    cisw = var(pm, nw)[:cisw] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_cisw_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :switch, l), "cisw_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_switch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_switch)
            cmax = ref(pm, nw, :switch, l)["current_rating"]
            set_upper_bound.(cisw[(l,i,j)],  cmax)
            set_lower_bound.(cisw[(l,i,j)], -cmax)
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :switch, :ci_fr, :ci_to, ref(pm, nw, :arcs_switch_from), ref(pm, nw, :arcs_switch_to), cisw)
end


"""
	function variable_mc_switch_power_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For models with explicit neutrals,
creates switch active power variables `:psw` for models with explicit neutrals.
This is defined per arc, i.e. with a variable for the from-side and to-side power.
"""
function variable_mc_switch_power_real(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(switch["f_connections"]) for (l,switch) in ref(pm, nw, :switch))
    psw = var(pm, nw)[:psw] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_psw_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :switch, l), "psw_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_switch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_switch)
            smax = ref(pm, nw, :switch, l)["thermal_rating"]
            set_upper_bound.(psw[(l,i,j)],  smax)
            set_lower_bound.(psw[(l,i,j)], -smax)
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :switch, :pf, :pt, ref(pm, nw, :arcs_switch_from), ref(pm, nw, :arcs_switch_to), psw)
end


"""
	function variable_mc_switch_power_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For models with explicit neutrals,
creates switch reactive power variables `:qsw` for models with explicit neutrals.
This is defined per arc, i.e. with a variable for the from-side and to-side power.
"""
function variable_mc_switch_power_imaginary(pm::ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(switch["f_connections"]) for (l,switch) in ref(pm, nw, :switch))
    qsw = var(pm, nw)[:qsw] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_qsw_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :switch, l), "qsw_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs_switch)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs_switch)
            smax = ref(pm, nw, :switch, l)["thermal_rating"]
            set_upper_bound.(qsw[(l,i,j)],  smax)
            set_lower_bound.(qsw[(l,i,j)], -smax)
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :switch, :qf, :qt, ref(pm, nw, :arcs_switch_from), ref(pm, nw, :arcs_switch_to), qsw)
end


