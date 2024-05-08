# GENERATOR

# GENERATOR - Variables - Non-linear

"""
	function variable_mc_generator_power(
		pm::AbstractNLExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true
	)

For IVR models with explicit neutrals,
no power variables are required
"""
function variable_mc_generator_power(pm::AbstractNLExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # do nothing
end


# GENERATOR - Variables - Quadratic

"""
	function variable_mc_generator_power(
		pm::AbstractQuadraticExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true
	)

For quadratic IVR models with explicit neutrals,
creates generator power variables `:pg` and `:qg`
"""
function variable_mc_generator_power(pm::AbstractQuadraticExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_generator_power_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_generator_power_imaginary(pm; nw=nw, bounded=bounded, report=report)
end


# GENERATOR - Constraints

"""
	function constraint_mc_generator_current(
		pm::AbstractExplicitNeutralIVRModel,
		id::Int;
		nw::Int=nw_id_default,
		report::Bool=true,
		bounded::Bool=true
	)

For IVR models with explicit neutrals,
creates expressions for the terminal current flows `:crg_bus` and `:cig_bus`.
"""
function constraint_mc_generator_current(pm::AbstractExplicitNeutralIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)
    generator = ref(pm, nw, :gen, id)

    nphases = _infer_int_dim_unit(generator, false)
    # Note that one-dimensional delta generators are handled as wye-connected generators.
    # The distinction between one-dimensional wye and delta generators is purely semantic
    # when neutrals are modeled explicitly.
    if get(generator, "configuration", WYE) == WYE || nphases==1
        constraint_mc_generator_current_wye(pm, nw, id, generator["connections"]; report=report, bounded=bounded)
    else
        constraint_mc_generator_current_delta(pm, nw, id, generator["connections"]; report=report, bounded=bounded)
    end
end


"""
	function constraint_mc_generator_current_wye(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		connections::Vector{Int};
		report::Bool=true,
		bounded::Bool=true
	)

For IVR models with explicit neutrals,
creates expressions for the terminal current flows `:crg_bus` and `:cig_bus` of wye-connected generators
"""
function constraint_mc_generator_current_wye(pm::AbstractExplicitNeutralIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)
    var(pm, nw, :crg_bus)[id] = _merge_bus_flows(pm, [crg..., -sum(crg)], connections)
    var(pm, nw, :cig_bus)[id] = _merge_bus_flows(pm, [cig..., -sum(cig)], connections)
end


"""
	function constraint_mc_generator_current_delta(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		connections::Vector{Int};
		report::Bool=true,
		bounded::Bool=true
	)

For IVR models with explicit neutrals,
creates expressions for the terminal current flows `:crg_bus` and `:cig_bus` of delta-connected generators
"""
function constraint_mc_generator_current_delta(pm::AbstractExplicitNeutralIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)
    Md = _get_delta_transformation_matrix(length(connections))
    var(pm, nw, :crg_bus)[id] = _merge_bus_flows(pm, Md'*crg, connections)
    var(pm, nw, :cig_bus)[id] = _merge_bus_flows(pm, Md'*cig, connections)
end


# GENERATOR - Constraints - Non-linear

"""
	function constraint_mc_generator_power_wye(
		pm::AbstractNLExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
		pmin::Vector{<:Real},
		pmax::Vector{<:Real},
		qmin::Vector{<:Real},
		qmax::Vector{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
creates non-linear expressions for the generator power `:pd` and `:qd`
of wye-connected generators as a function of voltage and current
"""
function constraint_mc_generator_power_wye(pm::AbstractNLExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)

    phases = connections[1:end-1]
    n      = connections[end]

    pg = Vector{JuMP.NonlinearExpr}([])
    qg = Vector{JuMP.NonlinearExpr}([])

    for (idx, p) in enumerate(phases)
        push!(pg, JuMP.@expression(pm.model,  ((vr[p]-vr[n])*crg[idx])+((vi[p]-vi[n])*cig[idx])))
        push!(qg, JuMP.@expression(pm.model, -((vr[p]-vr[n])*cig[idx])+((vi[p]-vi[n])*crg[idx])))
    end

    for (idx, p) in enumerate(phases)
        if pmin[idx]>-Inf
            JuMP.@constraint(pm.model, pmin[idx] .<= (vr[p]-vr[n])*crg[idx]  + (vi[p]-vi[n])*cig[idx])
        end
        if pmax[idx]< Inf
            JuMP.@constraint(pm.model, pmax[idx] .>= (vr[p]-vr[n])*crg[idx]  + (vi[p]-vi[n])*cig[idx])
        end
        if qmin[idx]>-Inf
            JuMP.@constraint(pm.model, qmin[idx] .<= (vi[p]-vi[n])*crg[idx]  - (vr[p]-vr[n])*cig[idx])
        end
        if qmax[idx]< Inf
            JuMP.@constraint(pm.model, qmax[idx] .>= (vi[p]-vi[n])*crg[idx]  - (vr[p]-vr[n])*cig[idx])
        end
    end

    var(pm, nw, :pg)[id] = pg
    var(pm, nw, :qg)[id] = qg

    if report
        sol(pm, nw, :gen, id)[:pg] = pg
        sol(pm, nw, :gen, id)[:qg] = qg
    end
end


"""
	function constraint_mc_generator_power_delta(
		pm::AbstractNLExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
		pmin::Vector{<:Real},
		pmax::Vector{<:Real},
		qmin::Vector{<:Real},
		qmax::Vector{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
creates non-linear expressions for the generator power `:pd` and `:qd`
of delta-connected generators as a function of voltage and current
"""
function constraint_mc_generator_power_delta(pm::AbstractNLExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)

    nph = length(pmin)

    vrg = Dict()
    vig = Dict()
    for (idx,c,d) in zip(1:nph, connections, [connections[2:end]..., connections[1]])
        vrg[idx] = JuMP.@expression(pm.model, vr[c]-vr[d])
        vig[idx] = JuMP.@expression(pm.model, vi[c]-vi[d])
    end

    pg = Vector{JuMP.NonlinearExpr}([])
    qg = Vector{JuMP.NonlinearExpr}([])
    for idx in 1:nph
        push!(pg, JuMP.@expression(pm.model,  (vrg[idx]*crg[idx])+(vig[idx]*cig[idx])))
        push!(qg, JuMP.@expression(pm.model, -(vrg[idx]*cig[idx])+(vig[idx]*crg[idx])))
    end

    JuMP.@constraint(pm.model, [i in 1:nph], pmin[i] <= pg[i])
    JuMP.@constraint(pm.model, [i in 1:nph], pmax[i] >= pg[i])
    JuMP.@constraint(pm.model, [i in 1:nph], qmin[i] <= qg[i])
    JuMP.@constraint(pm.model, [i in 1:nph], qmax[i] >= qg[i])

    var(pm, nw, :pg)[id] = JuMP.Containers.DenseAxisArray(pg, connections)
    var(pm, nw, :qg)[id] = JuMP.Containers.DenseAxisArray(qg, connections)

    if report
        sol(pm, nw, :gen, id)[:pg] = pg
        sol(pm, nw, :gen, id)[:qg] = qg
    end
end


# GENERATOR - Constraints - Quadratic

"""
	function constraint_mc_generator_power_wye(
		pm::AbstractQuadraticExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
		pmin::Vector{<:Real},
		pmax::Vector{<:Real},
		qmin::Vector{<:Real},
		qmax::Vector{<:Real};
		report::Bool=true,
		bounded::Bool=true
	)

For quadratic IVR models with explicit neutrals,
links the generator power variables `:pd` and `:qd`
of wye-connected generators to the voltage and current
"""
function constraint_mc_generator_power_wye(pm::AbstractQuadraticExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)

    phases = connections[1:end-1]
    n      = connections[end]

    vr_pn = [vr[p]-vr[n] for p in phases]
    vi_pn = [vi[p]-vi[n] for p in phases]

    pg = var(pm, nw, :pg, id)
    qg = var(pm, nw, :qg, id)

    JuMP.@constraint(pm.model, pg .==  vr_pn.*crg .+ vi_pn.*cig)
    JuMP.@constraint(pm.model, qg .== -vr_pn.*cig .+ vi_pn.*crg)
end


"""
	function constraint_mc_generator_power_delta(
		pm::AbstractQuadraticExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
		pmin::Vector{<:Real},
		pmax::Vector{<:Real},
		qmin::Vector{<:Real},
		qmax::Vector{<:Real};
		report::Bool=true,
		bounded::Bool=true
	)

For quadratic IVR models with explicit neutrals,
links the generator power variables `:pd` and `:qd`
of delta-connected generators to the voltage and current
"""
function constraint_mc_generator_power_delta(pm::AbstractQuadraticExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)

    ph = connections
    ph_next = [connections[2:end]..., connections[1]]

    vrg = [vr[c]-vr[d] for (c,d) in zip(ph,ph_next)]
    vig = [vi[c]-vi[d] for (c,d) in zip(ph,ph_next)]

    pg = var(pm, nw, :pg, id)
    qg = var(pm, nw, :qg, id)

    JuMP.@constraint(pm.model, pg .==  vrg.*crg .+ vig.*cig)
    JuMP.@constraint(pm.model, qg .== -vrg.*cig .+ vig.*crg)
end


# LOAD

# LOAD - Variables

"""
	function variable_mc_load_current(
		pm::AbstractExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For IVR models with explicit neutrals,
creates placeholder dictionaries for the load current `:crd` and `:cid`,
and for the terminal current flows `:crd_bus` and `:cid_bus`
"""
function variable_mc_load_current(pm::AbstractExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    var(pm, nw)[:crd] = Dict{Int, Any}()
    var(pm, nw)[:cid] = Dict{Int, Any}()
    var(pm, nw)[:crd_bus] = Dict{Int, Any}()
    var(pm, nw)[:cid_bus] = Dict{Int, Any}()
end


"""
	function variable_mc_load_power(
		pm::AbstractNLExplicitNeutralIVRModel;
		nw=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For non-linear IVR models with explicit neutrals,
creates placeholder dictionaries for the load power `:pd` and `:qd`,
and for the terminal power flows `:pd_bus` and `:qd_bus`
"""
function variable_mc_load_power(pm::AbstractNLExplicitNeutralIVRModel; nw=nw_id_default, bounded::Bool=true, report::Bool=true)
    var(pm, nw)[:pd] = Dict{Int, Any}()
    var(pm, nw)[:qd] = Dict{Int, Any}()
    var(pm, nw)[:pd_bus] = Dict{Int, Any}()
    var(pm, nw)[:qd_bus] = Dict{Int, Any}()
end


# LOAD - Variables - Quadratic


"""
	function variable_mc_load_current(
		pm::AbstractQuadraticExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true
	)

For quadratic IVR models with explicit neutrals,
creates load current variables `:crd` and `:cid`,
and placeholder dictionaries for the terminal current flows `:crd_bus` and `:cid_bus`
"""
function variable_mc_load_current(pm::AbstractQuadraticExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    load_ids_exponential = [id for (id,load) in ref(pm, nw, :load) if load["model"]==EXPONENTIAL]
    @assert isempty(load_ids_exponential) "Exponential loads cannot be represented quadratically."

    variable_mc_load_current_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_load_current_imaginary(pm; nw=nw, bounded=bounded, report=report)

    var(pm, nw)[:crd_bus] = Dict{Int,Any}()
    var(pm, nw)[:cid_bus] = Dict{Int,Any}()
end

"""
	function variable_mc_load_power(
		pm::AbstractQuadraticExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true
	)

For quadratic IVR models with explicit neutrals,
creates load power variables `:pd` and `:qd`
"""
function variable_mc_load_power(pm::AbstractQuadraticExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    load_ids_exponential = [id for (id,load) in ref(pm, nw, :load) if load["model"]==EXPONENTIAL]
    @assert isempty(load_ids_exponential) "Exponential loads cannot be represented quadratically."

    variable_mc_load_power_active(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_load_power_reactive(pm; nw=nw, bounded=bounded, report=report)
end


"""
	function variable_mc_load_power_active(
		pm::AbstractQuadraticExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates load active power variables `:pd` for models with explicit neutrals
"""
function variable_mc_load_power_active(pm::AbstractQuadraticExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(load, false) for (i,load) in ref(pm, nw, :load))
    load_ids_current = [id for (id,load) in ref(pm, nw, :load) if load["model"]==CURRENT]

    pd = var(pm, nw)[:pd] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_pd_$(i)",
            start = comp_start_value(ref(pm, nw, :load, i), "pd_start", c, 0.0)
        ) for i in load_ids_current
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :pd, load_ids_current, pd)
end


"""
	function variable_mc_load_power_reactive(
		pm::AbstractQuadraticExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates load reactive power variables `:qd` for models with explicit neutrals
"""
function variable_mc_load_power_reactive(pm::AbstractQuadraticExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _infer_int_dim_unit(load, false) for (i,load) in ref(pm, nw, :load))
    load_ids_current = [id for (id,load) in ref(pm, nw, :load) if load["model"]==CURRENT]

    qd = var(pm, nw)[:qd] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_qd_$(i)",
            start = comp_start_value(ref(pm, nw, :load, i), "qd_start", c, 0.0)
        ) for i in load_ids_current
    )

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :qd, load_ids_current, qd)
end


# LOAD - Constraints

"""
	function constraint_mc_load_power(
		pm::AbstractExplicitNeutralIVRModel,
		id::Int;
		nw::Int=nw_id_default,
		report::Bool=true
	)

For IVR models with explicit neutrals,
the load power does not require any constraints.
"""
function constraint_mc_load_power(pm::AbstractExplicitNeutralIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
    # nothing to do
end

"""
	function constraint_mc_load_current(
		pm::AbstractExplicitNeutralIVRModel,
		id::Int;
		nw::Int=nw_id_default,
		report::Bool=true
	)

For IVR models with explicit neutrals,
create non-linear expressions for the terminal current flows `:crd_bus` and `:cid_bus`
"""
function constraint_mc_load_current(pm::AbstractExplicitNeutralIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
    load = ref(pm, nw, :load, id)
    bus = ref(pm, nw,:bus, load["load_bus"])

    configuration = load["configuration"]

    a, alpha, b, beta = _load_expmodel_params(load, bus)

    int_dim = _infer_int_dim_unit(load, false)
    if configuration==WYE || int_dim==1
        constraint_mc_load_current_wye(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    else
        constraint_mc_load_current_delta(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    end
end


"""
	function constraint_mc_load_current_wye(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
		a::Vector{<:Real},
		alpha::Vector{<:Real},
		b::Vector{<:Real},
		beta::Vector{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
create non-linear expressions for the terminal current flows `:crd_bus` and `:cid_bus`
of wye-connected loads
"""
function constraint_mc_load_current_wye(pm::AbstractExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    crd = Vector{JuMP.NonlinearExpr}([])
    cid = Vector{JuMP.NonlinearExpr}([])

    phases = connections[1:end-1]
    n      = connections[end]

    for (idx, p) in enumerate(phases)
        push!(crd, JuMP.@expression(pm.model,
             (a[idx]*(vr[p]-vr[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(alpha[idx]/2-1))
            +(b[idx]*(vi[p]-vi[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(beta[idx]/2 -1))
        ))
        push!(cid, JuMP.@expression(pm.model,
             (a[idx]*(vi[p]-vi[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(alpha[idx]/2-1))
            -(b[idx]*(vr[p]-vr[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(beta[idx]/2 -1))
        ))
    end

    var(pm, nw, :crd)[id] = crd
    var(pm, nw, :cid)[id] = cid

    crd_bus_n = JuMP.@expression(pm.model, -sum(crd[i] for i in 1:length(phases)))
    cid_bus_n = JuMP.@expression(pm.model, -sum(cid[i] for i in 1:length(phases)))

    var(pm, nw, :crd_bus)[id] = crd_bus = _merge_bus_flows(pm, [crd..., crd_bus_n], connections)
    var(pm, nw, :cid_bus)[id] = cid_bus = _merge_bus_flows(pm, [cid..., cid_bus_n], connections)

    if report
        pd_bus = Vector{JuMP.NonlinearExpr}([])
        qd_bus = Vector{JuMP.NonlinearExpr}([])
        for (idx,c) in enumerate(connections)
            push!(pd_bus, JuMP.@expression(pm.model,  (vr[c]*crd_bus[c])+(vi[c]*cid_bus[c])))
            push!(qd_bus, JuMP.@expression(pm.model, -(vr[c]*cid_bus[c])+(vi[c]*crd_bus[c])))
        end

        sol(pm, nw, :load, id)[:pd_bus] = JuMP.Containers.DenseAxisArray(pd_bus, connections)
        sol(pm, nw, :load, id)[:qd_bus] = JuMP.Containers.DenseAxisArray(qd_bus, connections)

        sol(pm, nw, :load, id)[:crd] = JuMP.Containers.DenseAxisArray(crd, connections)
        sol(pm, nw, :load, id)[:cid] = JuMP.Containers.DenseAxisArray(cid, connections)

        sol(pm, nw, :load, id)[:crd_bus] = crd_bus
        sol(pm, nw, :load, id)[:cid_bus] = cid_bus

        pd = Vector{JuMP.NonlinearExpr}([])
        qd = Vector{JuMP.NonlinearExpr}([])
        for (idx, p) in enumerate(phases)
            push!(pd, JuMP.@expression(pm.model, a[idx]*(vr[p]^2+vi[p]^2)^(alpha[idx]/2) ))
            push!(qd, JuMP.@expression(pm.model, b[idx]*(vr[p]^2+vi[p]^2)^(beta[idx]/2)  ))
        end
        sol(pm, nw, :load, id)[:pd] = JuMP.Containers.DenseAxisArray(pd, connections)
        sol(pm, nw, :load, id)[:qd] = JuMP.Containers.DenseAxisArray(qd, connections)
    end
end


"""
	function constraint_mc_load_current_delta(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
		a::Vector{<:Real},
		alpha::Vector{<:Real},
		b::Vector{<:Real},
		beta::Vector{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
create non-linear expressions for the terminal current flows `:crd_bus` and `:cid_bus`
of delta-connected loads
"""
function constraint_mc_load_current_delta(pm::AbstractExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)


    ph = connections
    ph_next = [connections[2:end]..., connections[1]]
    P = length(ph)
    idxs = 1:P
    idxs_prev = [idxs[end], idxs[1:end-1]...]

    vrd = [vr[c]-vr[d] for (c,d) in zip(ph,ph_next)]
    vid = [vi[c]-vi[d] for (c,d) in zip(ph,ph_next)]

    crd = JuMP.@expression(pm.model, [i in 1:P],
        (a[i]*vrd[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2-1))
       +(b[i]*vid[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2 -1))
    )
    cid = JuMP.@expression(pm.model, [i in 1:P],
        (a[i]*vid[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2-1))
       -(b[i]*vrd[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2 -1))
    )

    crd_bus = JuMP.@expression(pm.model, [i in 1:P], crd[i]-crd[idxs_prev[i]])
    cid_bus = JuMP.@expression(pm.model, [i in 1:P], cid[i]-cid[idxs_prev[i]])

    var(pm, nw, :crd_bus)[id] = _merge_bus_flows(pm, crd_bus, connections)
    var(pm, nw, :cid_bus)[id] = _merge_bus_flows(pm, cid_bus, connections)

    if report
        pd_bus = JuMP.@expression(pm.model, [i in 1:P],  (vr[i]*crd_bus[i])+(vi[i]*cid_bus[i]))
        qd_bus = JuMP.@expression(pm.model, [i in 1:P], -(vr[i]*cid_bus[i])+(vi[i]*crd_bus[i]))

        sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        pd = JuMP.@expression(pm.model, [i in 1:P], a[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2) )
        qd = JuMP.@expression(pm.model, [i in 1:P], b[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2)  )
        sol(pm, nw, :load, id)[:pd] = pd
        sol(pm, nw, :load, id)[:qd] = qd
    end
end


# LOAD - Constraints - Quadratic

"""
	function constraint_mc_load_power(
		pm::AbstractQuadraticExplicitNeutralIVRModel,
		id::Int;
		nw::Int=nw_id_default,
		report::Bool=true
	)

For quadratic IVR models with explicit neutrals,
link the load power variables `:pd` and `:qd` to the voltage,
and link together the power, voltage and current variables
"""
function constraint_mc_load_power(pm::AbstractQuadraticExplicitNeutralIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
    load = ref(pm, nw, :load, id)
    bus = ref(pm, nw,:bus, load["load_bus"])

    configuration = load["configuration"]
    int_dim = _infer_int_dim_unit(load, false)
    a, alpha, b, beta = _load_expmodel_params(load, bus)

    # Note that one-dimensional delta loads are handled as wye-connected loads.
    # The distinction between one-dimensional wye and delta loads is purely semantic
    # when neutrals are modeled explicitly.
    if configuration==WYE || int_dim==1
        constraint_mc_load_power_wye(pm, nw, id, load["load_bus"], load["connections"], load["model"], a, b; report=report)
    else
        constraint_mc_load_power_delta(pm, nw, id, load["load_bus"], load["connections"], load["model"], a, b; report=report)
    end
end


"""
	function constraint_mc_load_power_wye(
		pm::AbstractQuadraticExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
		model::LoadModel,
		a::Vector{<:Real},
		b::Vector{<:Real};
		report::Bool=true,
		bounded::Bool=true
	)

For quadratic IVR models with explicit neutrals,
link the load power variables `:pd` and `:qd` to the voltage,
and link together the power, voltage and current variables
for wye-connected loads
"""
function constraint_mc_load_power_wye(pm::AbstractQuadraticExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, model::LoadModel, a::Vector{<:Real}, b::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    crd = var(pm, nw, :crd, id)
    cid = var(pm, nw, :cid, id)

    phases = connections[1:end-1]
    n      = connections[end]

    vr_pn = [vr[p]-vr[n] for p in phases]
    vi_pn = [vi[p]-vi[n] for p in phases]

    if model==POWER
        pd = a
        qd = b
    elseif model==IMPEDANCE
        pd = a .* (vr_pn.^2 .+ vi_pn.^2)
        qd = b .* (vr_pn.^2 .+ vi_pn.^2)
    elseif model==CURRENT
        pd = var(pm, nw, :pd, id)
        qd = var(pm, nw, :qd, id)
        JuMP.@constraint(pm.model, pd.^2 .== a.^2 .* (vr_pn.^2 .+ vi_pn.^2))
        JuMP.@constraint(pm.model, sign.(a).*pd .>= 0)
        JuMP.@constraint(pm.model, qd.^2 .== b.^2 .* (vr_pn.^2 .+ vi_pn.^2))
        JuMP.@constraint(pm.model, sign.(b).*qd .>= 0)
    else
        error("Load model $model for load $id is not supported by this formulation.")
    end

    JuMP.@constraint(pm.model, pd .==  vr_pn.*crd .+ vi_pn.*cid)
    JuMP.@constraint(pm.model, qd .== -vr_pn.*cid .+ vi_pn.*crd)

    # constant current loads are already reported through variable function
    if report && model!=CURRENT
        sol(pm, nw, :load, id)[:pd] = pd
        sol(pm, nw, :load, id)[:qd] = qd
    end
end


"""
	function constraint_mc_load_power_delta(
		pm::AbstractQuadraticExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
		model::LoadModel,
		a::Vector{<:Real},
		b::Vector{<:Real};
		report::Bool=true,
		bounded::Bool=true
	)

For quadratic IVR models with explicit neutrals,
link the load power variables `:pd` and `:qd` to the voltage,
and link together the power, voltage and current variables
for delta-connected loads
"""
function constraint_mc_load_power_delta(pm::AbstractQuadraticExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, model::LoadModel, a::Vector{<:Real}, b::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    crd = var(pm, nw, :crd, id)
    cid = var(pm, nw, :cid, id)

    phases = connections

    Md = _get_delta_transformation_matrix(length(connections))
    vrd = Md*[vr[p] for p in phases]
    vid = Md*[vi[p] for p in phases]

    if model==POWER
        pd = a
        qd = b
    elseif model==IMPEDANCE
        pd = a .* (vrd.^2 .+ vid.^2)
        qd = b .* (vrd.^2 .+ vid.^2)
    elseif model==CURRENT
        pd = var(pm, nw, :pd, id)
        qd = var(pm, nw, :qd, id)
        JuMP.@constraint(pm.model, pd.^2 .== a.^2 .* (vrd.^2 .+ vid.^2))
        JuMP.@constraint(pm.model, sign.(a).*pd .>= 0)
        JuMP.@constraint(pm.model, qd.^2 .== b.^2 .* (vrd.^2 .+ vid.^2))
        JuMP.@constraint(pm.model, sign.(b).*qd .>= 0)
    else
        error("Load model $model for load $id is not supported by this formulation.")
    end

    JuMP.@constraint(pm.model, pd .==  vrd.*crd .+ vid.*cid)
    JuMP.@constraint(pm.model, qd .== -vrd.*cid .+ vid.*crd)

    # constant current loads are already reported through variable function
    if report && model!=CURRENT
        sol(pm, nw, :load, id)[:pd] = pd
        sol(pm, nw, :load, id)[:qd] = qd
    end
end


"""
	function constraint_mc_load_current(
		pm::AbstractQuadraticExplicitNeutralIVRModel,
		id::Int;
		nw::Int=nw_id_default,
		report::Bool=true,
		bounded::Bool=true
	)

For quadratic IVR models with explicit neutrals,
create expressions for the terminal current flows `:crd_bus` and `:cid_bus`
"""
function constraint_mc_load_current(pm::AbstractQuadraticExplicitNeutralIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)
    load = ref(pm, nw, :load, id)

    int_dim = _infer_int_dim_unit(load, false)
    # Note that one-dimensional delta loads are handled as wye-connected loads.
    # The distinction between one-dimensional wye and delta loads is purely semantic
    # when neutrals are modeled explicitly.
    if get(load, "configuration", WYE) == WYE || int_dim==1
        constraint_mc_load_current_wye(pm, nw, id, load["connections"]; report=report, bounded=bounded)
    else
        constraint_mc_load_current_delta(pm, nw, id, load["connections"]; report=report, bounded=bounded)
    end
end


"""
	function constraint_mc_load_current_wye(
		pm::AbstractQuadraticExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		connections::Vector{Int};
		report::Bool=true,
		bounded::Bool=true
	)

For quadratic IVR models with explicit neutrals,
create expressions for the terminal current flows `:crd_bus` and `:cid_bus`
for wye-connected loads
"""
function constraint_mc_load_current_wye(pm::AbstractQuadraticExplicitNeutralIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crd = var(pm, nw, :crd, id)
    cid = var(pm, nw, :cid, id)
    var(pm, nw, :crd_bus)[id] = _merge_bus_flows(pm, [crd..., -sum(crd)], connections)
    var(pm, nw, :cid_bus)[id] = _merge_bus_flows(pm, [cid..., -sum(cid)], connections)
end


"""
	function constraint_mc_load_current_delta(
		pm::AbstractQuadraticExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		connections::Vector{Int};
		report::Bool=true,
		bounded::Bool=true
	)

For quadratic IVR models with explicit neutrals,
create expressions for the terminal current flows `:crd_bus` and `:cid_bus`
for delta-connected loads
"""
function constraint_mc_load_current_delta(pm::AbstractQuadraticExplicitNeutralIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crd = var(pm, nw, :crd, id)
    cid = var(pm, nw, :cid, id)
    Md = _get_delta_transformation_matrix(length(connections))
    var(pm, nw, :crd_bus)[id] = _merge_bus_flows(pm, Md'*crd, connections)
    var(pm, nw, :cid_bus)[id] = _merge_bus_flows(pm, Md'*cid, connections)
end


# TRANSFORMER

# TRANSFORMER - Variables

"""
	function variable_mc_transformer_current(
		pm::AbstractExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
	)

For IVR models with explicit neutrals,
create transformer current variables `:crt` and `:cit`,
and placeholder dictionaries for the terminal current flows `:crt_bus` and `:cit_bus`
"""
function variable_mc_transformer_current(pm::AbstractExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_transformer_current_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_transformer_current_imaginary(pm; nw=nw, bounded=bounded, report=report)

    var(pm, nw)[:crt_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:cit_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end


# TRANSFORMER -  Variable - Non-linear

"""
	function variable_mc_transformer_power(
		pm::AbstractNLExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
	)

For non-linear IVR models with explicit neutrals,
no power variables are required.
"""
function variable_mc_transformer_power(pm::AbstractNLExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # do nothing
end


# TRANSFORMER - Variable - Quadratic

"""
	function variable_mc_transformer_power(
		pm::AbstractQuadraticExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
	)

For quadratic IVR models with explicit neutrals,
creates transformer power variables `:pt` and `:qt`
"""
function variable_mc_transformer_power(pm::AbstractQuadraticExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_transformer_power_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_transformer_power_imaginary(pm; nw=nw, bounded=bounded, report=report)
end


# TRANSFORMER - Constraints

"""
	function constraint_mc_transformer_current(
		pm::AbstractExplicitNeutralIVRModel,
		i::Int;
		nw::Int=nw_id_default,
		fix_taps::Bool=true
	)

For IVR models with explicit neutrals,
links the current variables of the from-side and to-side transformer windings,
and creates expressions for the terminal current flows
"""
function constraint_mc_transformer_current(pm::AbstractExplicitNeutralIVRModel, i::Int; nw::Int=nw_id_default, fix_taps::Bool=true)
    # if ref(pm, nw_id_default, :conductors)!=3
    #     error("Transformers only work with networks with three conductors.")
    # end

    transformer = ref(pm, nw, :transformer, i)
    f_bus = transformer["f_bus"]
    t_bus = transformer["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    configuration = transformer["configuration"]
    f_connections = transformer["f_connections"]
    t_connections = transformer["t_connections"]
    tm_set = transformer["tm_set"]
    tm_fixed = fix_taps ? ones(Bool, length(tm_set)) : transformer["tm_fix"]
    tm_scale = calculate_tm_scale(transformer, ref(pm, nw, :bus, f_bus), ref(pm, nw, :bus, t_bus))

    #TODO change data model
    # there is redundancy in specifying polarity seperately on from and to side
    #TODO change this once migrated to new data model
    pol = transformer["polarity"]

    if configuration == WYE
        constraint_mc_transformer_current_yy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == DELTA
        constraint_mc_transformer_current_dy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == "zig-zag"
        error("Zig-zag not yet supported.")
    end
end


"""
	function constraint_mc_transformer_current_yy(
		pm::AbstractExplicitNeutralIVRModel,
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

For IVR models with explicit neutrals,
links the current variables of the from-side and to-side transformer windings,
and creates expressions for the terminal current flows
for wye-wye connected transformers

```
scale*cr_fr_P + cr_to_P == 0
scale*ci_fr_P + ci_to_P == 0
```
"""
function constraint_mc_transformer_current_yy(pm::AbstractExplicitNeutralIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    cr_fr_P = var(pm, nw, :crt, f_idx)
    ci_fr_P = var(pm, nw, :cit, f_idx)
    cr_to_P = var(pm, nw, :crt, t_idx)
    ci_to_P = var(pm, nw, :cit, t_idx)

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
    scale = (tm_scale*pol).*tm_set

    JuMP.@constraint(pm.model, scale.*cr_fr_P .+ cr_to_P .== 0)
    JuMP.@constraint(pm.model, scale.*ci_fr_P .+ ci_to_P .== 0)

    var(pm, nw, :crt_bus)[f_idx] = _merge_bus_flows(pm, [cr_fr_P..., -sum(cr_fr_P)], f_connections)
    var(pm, nw, :cit_bus)[f_idx] = _merge_bus_flows(pm, [ci_fr_P..., -sum(ci_fr_P)], f_connections)
    var(pm, nw, :crt_bus)[t_idx] = _merge_bus_flows(pm, [cr_to_P..., -sum(cr_to_P)], t_connections)
    var(pm, nw, :cit_bus)[t_idx] = _merge_bus_flows(pm, [ci_to_P..., -sum(ci_to_P)], t_connections)
end


"""
	function constraint_mc_transformer_current_dy(
		pm::AbstractExplicitNeutralIVRModel,
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

For IVR models with explicit neutrals,
links the current variables of the from-side and to-side transformer windings,
and creates expressions for the terminal current flows
for delta-wye connected transformers

```
scale*cr_fr_P + cr_to_P == 0
scale*ci_fr_P + ci_to_P == 0
```
"""
function constraint_mc_transformer_current_dy(pm::AbstractExplicitNeutralIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    cr_fr_P = var(pm, nw, :crt, f_idx)
    ci_fr_P = var(pm, nw, :cit, f_idx)
    cr_to_P = var(pm, nw, :crt, t_idx)
    ci_to_P = var(pm, nw, :cit, t_idx)

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
    scale = (tm_scale*pol).*tm_set

    n_phases = length(tm)
    Md = _get_delta_transformation_matrix(n_phases)

    JuMP.@constraint(pm.model, scale.*cr_fr_P .+ cr_to_P .== 0)
    JuMP.@constraint(pm.model, scale.*ci_fr_P .+ ci_to_P .== 0)

    var(pm, nw, :crt_bus)[f_idx] = _merge_bus_flows(pm, Md'*cr_fr_P, f_connections)
    var(pm, nw, :cit_bus)[f_idx] = _merge_bus_flows(pm, Md'*ci_fr_P, f_connections)
    var(pm, nw, :crt_bus)[t_idx] = _merge_bus_flows(pm, [cr_to_P..., -sum(cr_to_P)], t_connections)
    var(pm, nw, :cit_bus)[t_idx] = _merge_bus_flows(pm, [ci_to_P..., -sum(ci_to_P)], t_connections)
end


"""
	function constraint_mc_transformer_thermal_limit(
		pm::AbstractNLExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		f_idx::Tuple,
		t_idx::Tuple,
		f_bus::Int,
		t_bus::Int,
		f_connections::Vector,
		t_connections::Vector,
		config::ConnConfig,
		sm_ub::Real;
		report::Bool=true
	)

For non-linear IVR models with explicit neutrals,
imposes a bound on the magnitude of the total apparent power at both windings.
Expressions are created for the transformer power variables.

```
sum(pt_fr)^2 + sum(qt_fr)^2 <= sm_ub^2
sum(pt_to)^2 + sum(qt_to)^2 <= sm_ub^2
```
"""
function constraint_mc_transformer_thermal_limit(pm::AbstractNLExplicitNeutralIVRModel, nw::Int, id::Int, f_idx::Tuple, t_idx::Tuple, f_bus::Int, t_bus::Int, f_connections::Vector, t_connections::Vector, config::ConnConfig, sm_ub::Real; report::Bool=true)
    vr_fr = var(pm, nw, :vr, f_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    crt_fr = var(pm, nw, :crt, f_idx)
    cit_fr = var(pm, nw, :cit, f_idx)
    crt_to = var(pm, nw, :crt, t_idx)
    cit_to = var(pm, nw, :cit, t_idx)

    if config==WYE || length(crt_fr)==1
        P_fr = f_connections[1:end-1]
        n_fr = f_connections[end]
        vrt_fr = [vr_fr[p]-vr_fr[n_fr] for p in P_fr]
        vit_fr = [vi_fr[p]-vi_fr[n_fr] for p in P_fr]
    elseif config==DELTA && length(crt_fr)==3
        M = _get_delta_transformation_matrix(3)
        vrt_fr = M*[vr_to[p] for p in f_connections]
        vit_fr = M*[vi_to[p] for p in f_connections]
    else
        error("The configuration $config of dimension $(length(crt)) is not supported.")
    end

    P_to = t_connections[1:end-1]
    n_to = t_connections[end]
    vrt_to = [vr_to[p]-vr_to[n_to] for p in P_to]
    vit_to = [vi_to[p]-vi_to[n_to] for p in P_to]

    idxs = [1:length(vrt_fr)...]
    pt_fr = JuMP.@expression(pm.model, [i in idxs],  (vrt_fr[i]*crt_fr[i]) + (vit_fr[i]*cit_fr[i]))
    qt_fr = JuMP.@expression(pm.model, [i in idxs], -(vrt_fr[i]*cit_fr[i]) + (vit_fr[i]*crt_fr[i]))
    pt_to = JuMP.@expression(pm.model, [i in idxs],  (vrt_to[i]*crt_to[i]) + (vit_to[i]*cit_to[i]))
    qt_to = JuMP.@expression(pm.model, [i in idxs], -(vrt_to[i]*cit_to[i]) + (vit_to[i]*crt_to[i]))

    if sm_ub<Inf
        JuMP.@constraint(pm.model, sum(pt_fr[i] for i in idxs)^2 + sum(qt_fr[i] for i in idxs)^2 <= sm_ub^2)
        JuMP.@constraint(pm.model, sum(pt_to[i] for i in idxs)^2 + sum(qt_to[i] for i in idxs)^2 <= sm_ub^2)
    end

    if report
        sol(pm, nw, :transformer, id)[:pf] = pt_fr
        sol(pm, nw, :transformer, id)[:qf] = qt_fr
        sol(pm, nw, :transformer, id)[:pt] = pt_to
        sol(pm, nw, :transformer, id)[:qt] = qt_to
        sol(pm, nw, :transformer, id)[:smtot_fr] = JuMP.@expression(pm.model, sqrt(sum(pt_fr[i] for i in idxs)^2 + sum(qt_fr[i] for i in idxs)^2))
        sol(pm, nw, :transformer, id)[:smtot_to] = JuMP.@expression(pm.model, sqrt(sum(pt_to[i] for i in idxs)^2 + sum(qt_to[i] for i in idxs)^2))
    end
end


# TRANSFORMER - Constraint - Quadratic

"""
	function constraint_mc_transformer_thermal_limit(
		pm::AbstractQuadraticExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		f_idx::Tuple,
		t_idx::Tuple,
		f_bus::Int,
		t_bus::Int,
		f_connections::Vector,
		t_connections::Vector,
		config::ConnConfig,
		sm_ub::Real;
		report::Bool=true
	)

For quadratic IVR models with explicit neutrals,
imposes a bound on the magnitude of the total apparent power at both windings.

```
sum(pt_fr)^2 + sum(qt_fr)^2 <= sm_ub^2
sum(pt_to)^2 + sum(qt_to)^2 <= sm_ub^2
```
"""
function constraint_mc_transformer_thermal_limit(pm::AbstractQuadraticExplicitNeutralIVRModel, nw::Int, id::Int, f_idx::Tuple, t_idx::Tuple, f_bus::Int, t_bus::Int, f_connections::Vector, t_connections::Vector, config::ConnConfig, sm_ub::Real; report::Bool=true)
    vr_fr = var(pm, nw, :vr, f_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    crt_fr = var(pm, nw, :crt, f_idx)
    cit_fr = var(pm, nw, :cit, f_idx)
    crt_to = var(pm, nw, :crt, t_idx)
    cit_to = var(pm, nw, :cit, t_idx)

    if config==WYE || length(crt_fr)==1
        P_fr = f_connections[1:end-1]
        n_fr = f_connections[end]
        vrt_fr = [vr_fr[p]-vr_fr[n_fr] for p in P_fr]
        vit_fr = [vi_fr[p]-vi_fr[n_fr] for p in P_fr]
    elseif config==DELTA && length(crt_fr)==3
        M = _get_delta_transformation_matrix(3)
        vrt_fr = M*[vr_to[p] for p in f_connections]
        vit_fr = M*[vi_to[p] for p in f_connections]
    else
        error("The configuration $config of dimension $(length(crt)) is not supported.")
    end

    P_to = t_connections[1:end-1]
    n_to = t_connections[end]
    vrt_to = [vr_to[p]-vr_to[n_to] for p in P_to]
    vit_to = [vi_to[p]-vi_to[n_to] for p in P_to]

    pt_fr = var(pm, nw, :pt, f_idx)
    qt_fr = var(pm, nw, :qt, f_idx)
    pt_to = var(pm, nw, :pt, t_idx)
    qt_to = var(pm, nw, :qt, t_idx)

    JuMP.@constraint(pm.model, pt_fr .==  vrt_fr.*crt_fr .+ vit_fr.*cit_fr)
    JuMP.@constraint(pm.model, qt_fr .== -vrt_fr.*cit_fr .+ vit_fr.*crt_fr)
    JuMP.@constraint(pm.model, pt_to .==  vrt_to.*crt_to .+ vit_to.*cit_to)
    JuMP.@constraint(pm.model, qt_to .== -vrt_to.*cit_to .+ vit_to.*crt_to)

    if sm_ub < Inf
        JuMP.@constraint(pm.model, sum(pt_fr)^2+sum(qt_fr)^2 <= sm_ub^2)
        JuMP.@constraint(pm.model, sum(pt_to)^2+sum(qt_to)^2 <= sm_ub^2)
    end
end


# BRANCH

# BRANCH - Variables

"""
	function variable_mc_branch_current(
		pm::AbstractExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For IVR models with explicit neutrals,
creates total current variables `:cr` and `:ci`,
series current variables `:csr` and `:csi`,
and placeholder dictionaries for the terminal current flows `:cr_bus` and `:ci_bus`
"""
function variable_mc_branch_current(pm::AbstractExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_branch_current_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_branch_current_imaginary(pm; nw=nw, bounded=bounded, report=report)

    variable_mc_branch_current_series_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_branch_current_series_imaginary(pm; nw=nw, bounded=bounded, report=report)

    var(pm, nw)[:cr_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:ci_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end


# BRANCH - Variables - Reduced

"""
	function variable_mc_branch_current(
		pm::ReducedExplicitNeutralIVRModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For branch-reduced IVR models with explicit neutrals,
creates series current variables `:csr` and `:csi`,
placeholder dictionaries for the total current `:cr` and `:ci`,
and placeholder dictionaries for the terminal current flows `:cr_bus` and `:ci_bus`
"""
function variable_mc_branch_current(pm::ReducedExplicitNeutralIVRModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_branch_current_series_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_branch_current_series_imaginary(pm; nw=nw, bounded=bounded, report=report)

    var(pm, nw)[:cr] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:ci] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:cr_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:ci_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end


# BRANCH - Constraints

"""
	function constraint_mc_current_from(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		f_bus::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		g_sh_fr::Matrix{<:Real},
		b_sh_fr::Matrix{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
defines how current distributes over series and shunt impedances of a pi-model branch.

```
cr_fr == csr_fr + g_sh_fr*vr_fr - b_sh_fr*vi_fr
ci_fr == csi_fr + g_sh_fr*vi_fr + b_sh_fr*vr_fr
```
"""
function constraint_mc_current_from(pm::AbstractExplicitNeutralIVRModel, nw::Int, f_bus::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real}; report::Bool=true)
    vr_fr = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]

    cr_fr =  var(pm, nw, :cr, f_idx)
    ci_fr =  var(pm, nw, :ci, f_idx)

    csr_fr =  var(pm, nw, :csr, f_idx[1])
    csi_fr =  var(pm, nw, :csi, f_idx[1])

    JuMP.@constraint(pm.model, cr_fr .== csr_fr + g_sh_fr*vr_fr - b_sh_fr*vi_fr)
    JuMP.@constraint(pm.model, ci_fr .== csi_fr + g_sh_fr*vi_fr + b_sh_fr*vr_fr)

    var(pm, nw, :cr_bus)[f_idx] = cr_bus_fr = _merge_bus_flows(pm, cr_fr, f_connections)
    var(pm, nw, :ci_bus)[f_idx] = ci_bus_fr = _merge_bus_flows(pm, ci_fr, f_connections)

    if report
        sol(pm, nw, :branch, f_idx[1])[:pf] =  cr_fr.*vr_fr .+ ci_fr.*vi_fr
        sol(pm, nw, :branch, f_idx[1])[:qf] = -cr_fr.*vi_fr .+ ci_fr.*vr_fr
    end

end


"""
	function constraint_mc_current_to(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		t_bus,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int},
		g_sh_to::Matrix{<:Real},
		b_sh_to::Matrix{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
defines how current distributes over series and shunt impedances of a pi-model branch.

```
cr_to == csr_to + g_sh_to*vr_to - b_sh_to*vi_to
ci_to == csi_to + g_sh_to*vi_to + b_sh_to*vr_to
```
"""
function constraint_mc_current_to(pm::AbstractExplicitNeutralIVRModel, nw::Int, t_bus, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, g_sh_to::Matrix{<:Real}, b_sh_to::Matrix{<:Real}; report::Bool=true)
    vr_to = [var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_to = [var(pm, nw, :vi, t_bus)[c] for c in t_connections]

    cr_to = var(pm, nw, :cr, t_idx)
    ci_to = var(pm, nw, :ci, t_idx)

    csr_to = -var(pm, nw, :csr, f_idx[1])
    csi_to = -var(pm, nw, :csi, f_idx[1])

    JuMP.@constraint(pm.model, cr_to .== csr_to + g_sh_to*vr_to - b_sh_to*vi_to)
    JuMP.@constraint(pm.model, ci_to .== csi_to + g_sh_to*vi_to + b_sh_to*vr_to)

    var(pm, nw, :cr_bus)[t_idx] = cr_bus_to = _merge_bus_flows(pm, cr_to, t_connections)
    var(pm, nw, :ci_bus)[t_idx] = ci_bus_to = _merge_bus_flows(pm, ci_to, t_connections)

    if report
        sol(pm, nw, :branch, t_idx[1])[:pt] =  cr_to.*vr_to .+ ci_to.*vi_to
        sol(pm, nw, :branch, t_idx[1])[:qt] = -cr_to.*vi_to .+ ci_to.*vr_to
    end
end


"""
	function constraint_mc_bus_voltage_drop(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		i::Int,
		f_bus::Int,
		t_bus::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int},
		r::Matrix{<:Real},
		x::Matrix{<:Real}
	)

For IVR models with explicit neutrals,
defines voltage drop over a branch, linking from and to side complex voltage.

```
vr_to == vr_fr - r*csr_fr + x*csi_fr
vi_to == vi_fr - r*csi_fr - x*csr_fr
```
"""
function constraint_mc_bus_voltage_drop(pm::AbstractExplicitNeutralIVRModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, r::Matrix{<:Real}, x::Matrix{<:Real})
    vr_fr = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]

    vr_to = [var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_to = [var(pm, nw, :vi, t_bus)[c] for c in t_connections]

    csr_fr = var(pm, nw, :csr, f_idx[1])
    csi_fr = var(pm, nw, :csi, f_idx[1])

    JuMP.@constraint(pm.model, vr_to .== vr_fr - r*csr_fr + x*csi_fr)
    JuMP.@constraint(pm.model, vi_to .== vi_fr - r*csi_fr - x*csr_fr)
end


"""
	function constraint_mc_branch_current_limit(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector,
		t_connections::Vector,
		c_rating::Vector{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
imposes a bound on the current magnitude per conductor
at both ends of the branch (total current, i.e. including shunt contributions).

```
cr_fr^2 + ci_fr^2 <= c_rating^2
cr_to^2 + ci_to^2 <= c_rating^2
```
"""
function constraint_mc_branch_current_limit(pm::AbstractExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector, t_connections::Vector, c_rating::Vector{<:Real}; report::Bool=true)
    cr_fr = var(pm, nw, :cr, f_idx)
    ci_fr = var(pm, nw, :ci, f_idx)
    cr_to = var(pm, nw, :cr, t_idx)
    ci_to = var(pm, nw, :ci, t_idx)

    cnds_finite_rating = [c for (c,r) in enumerate(c_rating) if r<Inf]
    JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_fr[c]^2+ci_fr[c]^2 <= c_rating[c]^2)
    JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_to[c]^2+ci_to[c]^2 <= c_rating[c]^2)
end


"""
	function constraint_mc_thermal_limit_from(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		rate_a::Vector{<:Real}
	)

For IVR models with explicit neutrals,
imposes a bound on the from-side line power magnitude.
"""
function constraint_mc_thermal_limit_from(pm::AbstractExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})
    vr_fr = [var(pm, nw, :vr, f_idx[1])[t] for t in f_connections]
    vi_fr = [var(pm, nw, :vi, f_idx[1])[t] for t in f_connections]
    cr_fr = var(pm, nw, :cr, f_idx)
    ci_fr = var(pm, nw, :ci, f_idx)

    for idx in 1:length(rate_a)
        if rate_a[idx]<Inf
            pf_idx = JuMP.@expression(pm.model,  (vr_fr[idx]*cr_fr[idx]) + (vi_fr[idx]*ci_fr[idx]))
            qf_idx = JuMP.@expression(pm.model, -(vr_fr[idx]*ci_fr[idx]) + (vi_fr[idx]*cr_fr[idx]))

            JuMP.@constraint(pm.model, pf_idx^2 + qf_idx^2 <= rate_a[idx]^2)
        end
    end
end


"""
	function constraint_mc_thermal_limit_to(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		t_idx::Tuple{Int,Int,Int},
		t_connections::Vector{Int},
		rate_a::Vector{<:Real}
	)

For IVR models with explicit neutrals,
imposes a bound on the to-side line power magnitude.
"""
function constraint_mc_thermal_limit_to(pm::AbstractExplicitNeutralIVRModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})
    vr_to = [var(pm, nw, :vr, t_idx[1])[t] for t in t_connections]
    vi_to = [var(pm, nw, :vi, t_idx[1])[t] for t in t_connections]
    cr_to = var(pm, nw, :cr, t_idx)
    ci_to = var(pm, nw, :ci, t_idx)

    for idx in 1:length(rate_a)
        if rate_a[idx]<Inf
            pt_idx = JuMP.@expression(pm.model,  (vr_to[idx]*cr_to[idx]) + (vi_to[idx]*ci_to[idx]))
            qt_idx = JuMP.@expression(pm.model, -(vr_to[idx]*ci_to[idx]) + (vi_to[idx]*cr_to[idx]))

            JuMP.@constraint(pm.model, pt_idx^2 + qt_idx^2 <= rate_a[idx]^2)
        end
    end
end


# BRANCH - Constraints - Quadratic

"""
	function constraint_mc_thermal_limit_from(
		pm::AbstractQuadraticExplicitNeutralIVRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		rate_a::Vector{<:Real}
	)

For quadratic IVR models with explicit neutrals,
throw an error because this cannot be represented quadratically
without introducing explicit power variables.
"""
function constraint_mc_thermal_limit_from(pm::AbstractQuadraticExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})
    if any(rate_a.<Inf)
        @warn("""
            A branch power bound cannot be represented quadratically in the default AbstractQuadraticExplicitNeutralIVRModel.
            Either extend this quadratic formulation by including explicit branch power variables, or use AbstractNLExplicitNeutralIVRModel instead.""")
    end
end


"""
	function constraint_mc_thermal_limit_to(
		pm::AbstractQuadraticExplicitNeutralIVRModel,
		nw::Int,
		t_idx::Tuple{Int,Int,Int},
		t_connections::Vector{Int},
		rate_a::Vector{<:Real}
	)

For quadratic IVR models with explicit neutrals,
throw an error because this cannot be represented quadratically
without introducing explicit power variables.
"""
function constraint_mc_thermal_limit_to(pm::AbstractQuadraticExplicitNeutralIVRModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})
    @warn("""
        A branch power bound cannot be represented quadratically in the default AbstractQuadraticExplicitNeutralIVRModel.
        Either extend this quadratic formulation by including explicit branch power variables, or use AbstractNLExplicitNeutralIVRModel instead.
        """)
end


# BRANCH - Constraints - Reduced

"""
	function constraint_mc_current_from(
		pm::ReducedExplicitNeutralIVRModels,
		nw::Int,
		f_bus::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		g_sh_fr::Matrix{<:Real},
		b_sh_fr::Matrix{<:Real};
		report::Bool=true
	)

For branch-reduced IVR models with explicit neutrals,
defines how current distributes over series and shunt impedances of a pi-model branch.

```
cr_fr = csr_fr + g_sh_fr*vr_fr - b_sh_fr*vi_fr
ci_fr = csi_fr + g_sh_fr*vi_fr + b_sh_fr*vr_fr
```
"""
function constraint_mc_current_from(pm::ReducedExplicitNeutralIVRModels, nw::Int, f_bus::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real}; report::Bool=true)
    vr_fr = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]

    csr_fr =  var(pm, nw, :csr, f_idx[1])
    csi_fr =  var(pm, nw, :csi, f_idx[1])

    cr_fr = csr_fr + g_sh_fr*vr_fr - b_sh_fr*vi_fr
    ci_fr = csi_fr + g_sh_fr*vi_fr + b_sh_fr*vr_fr

    var(pm, nw, :cr)[f_idx] = cr_fr
    var(pm, nw, :ci)[f_idx] = ci_fr

    var(pm, nw, :cr_bus)[f_idx] = cr_bus_fr = _merge_bus_flows(pm, cr_fr, f_connections)
    var(pm, nw, :ci_bus)[f_idx] = ci_bus_fr = _merge_bus_flows(pm, ci_fr, f_connections)

    if report
        sol(pm, nw, :branch, f_idx[1])[:cr_fr] = cr_fr
        sol(pm, nw, :branch, f_idx[1])[:ci_fr] = ci_fr
        sol(pm, nw, :branch, f_idx[1])[:pf] =  cr_fr.*vr_fr .+ ci_fr.*vi_fr
        sol(pm, nw, :branch, f_idx[1])[:qf] = -cr_fr.*vi_fr .+ ci_fr.*vr_fr
    end
end


"""
	function constraint_mc_current_to(
		pm::ReducedExplicitNeutralIVRModels,
		nw::Int,
		t_bus,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int},
		g_sh_to::Matrix{<:Real},
		b_sh_to::Matrix{<:Real};
		report::Bool=true
	)

For branch-reduced IVR models with explicit neutrals,
defines how current distributes over series and shunt impedances of a pi-model branch.

```
cr_to = csr_to + g_sh_to*vr_to - b_sh_to*vi_to
ci_to = csi_to + g_sh_to*vi_to + b_sh_to*vr_to
```
"""
function constraint_mc_current_to(pm::ReducedExplicitNeutralIVRModels, nw::Int, t_bus, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, g_sh_to::Matrix{<:Real}, b_sh_to::Matrix{<:Real}; report::Bool=true)
    vr_to = [var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_to = [var(pm, nw, :vi, t_bus)[c] for c in t_connections]

    csr_to = -var(pm, nw, :csr, f_idx[1])
    csi_to = -var(pm, nw, :csi, f_idx[1])

    cr_to = csr_to + g_sh_to*vr_to - b_sh_to*vi_to
    ci_to = csi_to + g_sh_to*vi_to + b_sh_to*vr_to

    var(pm, nw, :cr)[t_idx] = cr_to
    var(pm, nw, :ci)[t_idx] = ci_to

    var(pm, nw, :cr_bus)[t_idx] = cr_bus_to = _merge_bus_flows(pm, cr_to, t_connections)
    var(pm, nw, :ci_bus)[t_idx] = ci_bus_to = _merge_bus_flows(pm, ci_to, t_connections)

    if report
        sol(pm, nw, :branch, f_idx[1])[:cr_to] = cr_to
        sol(pm, nw, :branch, f_idx[1])[:ci_to] = ci_to
        sol(pm, nw, :branch, t_idx[1])[:pt] =  cr_to.*vr_to .+ ci_to.*vi_to
        sol(pm, nw, :branch, t_idx[1])[:qt] = -cr_to.*vi_to .+ ci_to.*vr_to
    end
end


# SWITCH

# SWITCH - Variables

"""
	function variable_mc_switch_current(
		pm::AbstractExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
	)

For IVR models with explicit neutrals,
creates switch current variables `:crs` and `:cis`,
and placeholder dictionaries for the terminal current flows `:crsw_bus` and `:cisw_bus`
"""
function variable_mc_switch_current(pm::AbstractExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_switch_current_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_switch_current_imaginary(pm; nw=nw, bounded=bounded, report=report)

    var(pm, nw)[:crsw_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:cisw_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end


"""
	function constraint_mc_switch_current(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int};
		report::Bool=true
	)

For IVR models with explicit neutrals,
create expressions for the terminal current flows `:crsw_bus` and `cisw_bus`,
and link the from-side to the to-side switch current
"""
function constraint_mc_switch_current(pm::AbstractExplicitNeutralIVRModel, nw::Int, id::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}; report::Bool=true)
    crsw_fr = var(pm, nw, :crsw, f_idx)
    cisw_fr = var(pm, nw, :cisw, f_idx)
    crsw_to = var(pm, nw, :crsw, t_idx)
    cisw_to = var(pm, nw, :cisw, t_idx)

    JuMP.@constraint(pm.model, crsw_fr .+ crsw_to .== 0)
    JuMP.@constraint(pm.model, cisw_fr .+ cisw_to .== 0)

    var(pm, nw, :crsw_bus)[f_idx] = _merge_bus_flows(pm, crsw_fr, f_connections)
    var(pm, nw, :cisw_bus)[f_idx] = _merge_bus_flows(pm, cisw_fr, f_connections)
    var(pm, nw, :crsw_bus)[t_idx] = _merge_bus_flows(pm, crsw_to, t_connections)
    var(pm, nw, :cisw_bus)[t_idx] = _merge_bus_flows(pm, cisw_to, t_connections)
end


"""
	function constraint_mc_switch_current_limit(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		connections::Vector{Int},
		rating::Vector{<:Real}
	)

For IVR models with explicit neutrals,
imposes a bound on the switch current magnitude per conductor.
Note that a bound on the from-side implies the same bound on the to-side current,
so it suffices to apply this only explicitly at the from-side.
"""
function constraint_mc_switch_current_limit(pm::AbstractExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, connections::Vector{Int}, rating::Vector{<:Real})::Nothing
    crsw = var(pm, nw, :crsw, f_idx)
    cisw = var(pm, nw, :cisw, f_idx)

    mu_cm_fr = JuMP.ConstraintRef[]
    for idx in 1:length(rating)
        if rating[idx] < Inf
            push!(mu_cm_fr, JuMP.@constraint(pm.model, crsw[idx]^2 + cisw[idx]^2 <= rating[idx]^2))
        end
    end

    con(pm, nw, :mu_cm_switch)[f_idx] = mu_cm_fr

    nothing
end


"""
	function constraint_mc_switch_thermal_limit(
		pm::AbstractNLExplicitNeutralIVRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		rating::Vector{<:Real}
	)

For IVR models with explicit neutrals,
imposes a bound on the switch power magnitude per conductor.
Note that a bound on the from-side implies the same bound on the to-side power
when the switch is closed (equal voltages), and also when it is open since the
power then equals zero on both ends.
"""
function constraint_mc_switch_thermal_limit(pm::AbstractNLExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{<:Real})::Nothing
    vr_fr = [var(pm, nw, :vr, f_idx[2])[t] for t in f_connections]
    vi_fr = [var(pm, nw, :vi, f_idx[2])[t] for t in f_connections]
    crsw_fr = var(pm, nw, :crsw, f_idx)
    cisw_fr = var(pm, nw, :cisw, f_idx)

    mu_sm_fr = JuMP.ConstraintRef[]
    for idx in 1:length(rating)
        if rating[idx] < Inf
            psw_fr_idx = JuMP.@expression(pm.model,  (vr_fr[idx]*crsw_fr[idx]) + (vi_fr[idx]*cisw_fr[idx]))
            qsw_fr_idx = JuMP.@expression(pm.model, -(vr_fr[idx]*cisw_fr[idx]) + (vi_fr[idx]*crsw_fr[idx]))

            push!(mu_sm_fr, JuMP.@constraint(pm.model, psw_fr_idx^2 + qsw_fr_idx^2 <= rating[idx]^2))
        end
    end

    con(pm, nw, :mu_sm_switch)[f_idx] = mu_sm_fr

    nothing
end


"""
	function constraint_mc_switch_thermal_limit(
		pm::AbstractQuadraticExplicitNeutralIVRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		rating::Vector{<:Real}
	)

For quadratic IVR models with explicit neutrals,
throw an error because this cannot be represented quadratically
without introducing explicit power variables.
"""
function constraint_mc_switch_thermal_limit(pm::AbstractQuadraticExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{<:Real})
    @warn("""
        A switch power bound cannot be represented quadratically in the default AbstractQuadraticExplicitNeutralIVRModel.
        Either extend this quadratic formulation by including explicit switch power variables, or use AbstractNLExplicitNeutralIVRModel instead.
        """)
end
