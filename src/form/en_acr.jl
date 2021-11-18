# GENERATOR

# GENERATOR - Variables

"""
	function variable_mc_generator_power(
		pm::AbstractExplicitNeutralACRModel;
		nw::Int=nw_id_default,
		kwargs...
	)

For ACR models with explicit neutrals,
creates generator power variables `:pg` and `:qg`,
and placeholder dictionaries for terminal power flows `:pg_bus` and `:qg_bus`,
"""
function variable_mc_generator_power(pm::AbstractExplicitNeutralACRModel; nw::Int=nw_id_default, kwargs...)
    variable_mc_generator_power_real(pm; nw=nw, kwargs...)
    variable_mc_generator_power_imaginary(pm; nw=nw, kwargs...)
    var(pm, nw)[:pg_bus] = Dict{Int, Any}()
    var(pm, nw)[:qg_bus] = Dict{Int, Any}()
end


# GENERATOR - Constraints

"""
	function constraint_mc_generator_power_wye(
		pm::AbstractExplicitNeutralACRModel,
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

For ACR models with explicit neutrals,
links the terminal power flows `:pg_bus` and `:qg_bus` to the power variables `:pg` and `:qg` for wye-connected generators
"""
function constraint_mc_generator_power_wye(pm::AbstractExplicitNeutralACRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    pg = var(pm, nw, :pg, id)
    qg = var(pm, nw, :qg, id)

    phases = connections[1:end-1]
    n      = connections[end]
    P      = length(phases)

    if iszero(vr[n]) && iszero(vi[n])
        pg_bus_unmerged = [pg..., 0.0]
        qg_bus_unmerged = [qg..., 0.0]
    else
        vr_pn = [vr[p]-vr[n] for p in phases]
        vi_pn = [vi[p]-vi[n] for p in phases]

        crg = JuMP.@NLexpression(pm.model, [idx in 1:P],
            ( pg[idx]*vr_pn[idx] + qg[idx]*vi_pn[idx] )/( vr_pn[idx]^2 + vi_pn[idx]^2 )
        )
        cig = JuMP.@NLexpression(pm.model, [idx in 1:P],
            ( pg[idx]*vi_pn[idx] - qg[idx]*vr_pn[idx] )/( vr_pn[idx]^2 + vi_pn[idx]^2 )
        )
        pg_bus_unmerged = [
            [JuMP.@NLexpression(pm.model, vr[p]*crg[idx]+vi[p]*cig[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, vr[n]*sum(-crg[idx] for idx in 1:P)+vi[n]*sum(-cig[idx] for idx in 1:P))
        ]
        qg_bus_unmerged = [
            [JuMP.@NLexpression(pm.model, -vr[p]*cig[idx]+vi[p]*crg[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, -vr[n]*sum(-cig[idx] for idx in 1:P)+vi[n]*sum(-crg[idx] for idx in 1:P))
        ]
    end
    var(pm, nw, :pg_bus)[id] = pg_bus = _merge_bus_flows(pm, pg_bus_unmerged, connections)
    var(pm, nw, :qg_bus)[id] = qg_bus = _merge_bus_flows(pm, qg_bus_unmerged, connections)

    if report
        sol(pm, nw, :gen, id)[:pg_bus] = pg_bus
        sol(pm, nw, :gen, id)[:qg_bus] = qg_bus
    end
end


"""
	function constraint_mc_generator_power_delta(
		pm::AbstractExplicitNeutralACRModel,
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

For ACR models with explicit neutrals,
links the terminal power flows `:pg_bus` and `:qg_bus` to the power variables `:pg` and `:qg` for delta-connected generators
"""
function constraint_mc_generator_power_delta(pm::AbstractExplicitNeutralACRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    pg = var(pm, nw, :pg, id)
    qg = var(pm, nw, :qg, id)

    ph = connections
    ph_next = [connections[2:end]..., connections[1]]
    P = length(ph)
    idxs = 1:length(ph)
    idxs_prev = [idxs[end], idxs[1:end-1]...]

    vr_pp = [vr[c]-vr[d] for (c,d) in zip(ph,ph_next)]
    vi_pp = [vi[c]-vi[d] for (c,d) in zip(ph,ph_next)]

    crg = JuMP.@NLexpression(pm.model, [idx in 1:P],
        ( pg[idx]*vr_pp[idx] + qg[idx]*vi_pp[idx] )/( vr_pp[idx]^2 + vi_pp[idx]^2 )
    )
    cig = JuMP.@NLexpression(pm.model, [idx in 1:P],
        ( pg[idx]*vi_pp[idx] - qg[idx]*vr_pp[idx] )/( vr_pp[idx]^2 + vi_pp[idx]^2 )
    )
    crg_bus = JuMP.@NLexpression(pm.model, [idx in 1:P], crg[idx] - crg[idxs_prev[idx]])
    cig_bus = JuMP.@NLexpression(pm.model, [idx in 1:P], cig[idx] - cig[idxs_prev[idx]])

    pg_bus_unmerged = [JuMP.@NLexpression(pm.model,  vr[p]*crg_bus[idx]+vi[p]*cig_bus[idx]) for (idx,p) in enumerate(ph)]
    qg_bus_unmerged = [JuMP.@NLexpression(pm.model, -vr[p]*cig_bus[idx]+vi[p]*crg_bus[idx]) for (idx,p) in enumerate(ph)]

    var(pm, nw, :pg_bus)[id] = pg_bus = _merge_bus_flows(pm, pg_bus_unmerged, connections)
    var(pm, nw, :qg_bus)[id] = qg_bus = _merge_bus_flows(pm, qg_bus_unmerged, connections)

    if report
        sol(pm, nw, :gen, id)[:pg_bus] = pg_bus
        sol(pm, nw, :gen, id)[:qg_bus] = qg_bus
    end
end


# LOAD

# LOAD - Variables

"""
	function variable_mc_load_power(
		pm::AbstractExplicitNeutralACRModel;
		nw=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For ACR models with explicit neutrals,
creates placeholder dictionaries for load power expressions `:pd` and `:qd`
"""
function variable_mc_load_power(pm::AbstractExplicitNeutralACRModel; nw=nw_id_default, bounded::Bool=true, report::Bool=true)
    var(pm, nw)[:pd] = Dict{Int, Any}()
    var(pm, nw)[:qd] = Dict{Int, Any}()
    var(pm, nw)[:pd_bus] = Dict{Int, Any}()
    var(pm, nw)[:qd_bus] = Dict{Int, Any}()
end


# LOAD - Constraints

"""
	function constraint_mc_load_power_wye(
		pm::AbstractExplicitNeutralACRModel,
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

For ACR models with explicit neutrals,
creates non-linear expressions for terminal power flows ':pd_bus' and ':qd_bus' of wye-connected loads
"""
function constraint_mc_load_power_wye(pm::AbstractExplicitNeutralACRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    phases = connections[1:end-1]
    n      = connections[end]
    P      = length(phases)

    vr_pn = [vr[p]-vr[n] for p in phases]
    vi_pn = [vi[p]-vi[n] for p in phases]

    crd = JuMP.@NLexpression(pm.model, [idx in 1:P], a[idx]*vr_pn[idx]*(vr_pn[idx]^2+vi_pn[idx]^2)^(alpha[idx]/2-1)+b[idx]*vi_pn[idx]*(vr_pn[idx]^2+vi_pn[idx]^2)^(beta[idx]/2 -1))
    cid = JuMP.@NLexpression(pm.model, [idx in 1:P], a[idx]*vi_pn[idx]*(vr_pn[idx]^2+vi_pn[idx]^2)^(alpha[idx]/2-1)-b[idx]*vr_pn[idx]*(vr_pn[idx]^2+vi_pn[idx]^2)^(beta[idx]/2 -1))


    # if constant power load
    if all(alpha.==0) && all(beta.==0)
        pd = a
        qd = b
    else
        pd = JuMP.@NLexpression(pm.model, [idx in 1:P],  vr_pn[idx]*crd[idx]+vi_pn[idx]*cid[idx])
        qd = JuMP.@NLexpression(pm.model, [idx in 1:P], -vr_pn[idx]*cid[idx]+vi_pn[idx]*crd[idx])
    end

    if iszero(vr[n]) && iszero(vi[n])
        pd_bus_unmerged = [pd..., 0.0]
        qd_bus_unmerged = [qd..., 0.0]
    else
        pd_bus_unmerged = [
            [JuMP.@NLexpression(pm.model, vr[p]*crd[idx]+vi[p]*cid[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, vr[n]*sum(-crd[idx] for idx in 1:P)+vi[n]*sum(-cid[idx] for idx in 1:P))
        ]
        qd_bus_unmerged = [
            [JuMP.@NLexpression(pm.model, -vr[p]*cid[idx]+vi[p]*crd[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, -vr[n]*sum(-cid[idx] for idx in 1:P)+vi[n]*sum(-crd[idx] for idx in 1:P))
        ]
    end
    var(pm, nw, :pd_bus)[id] = pd_bus = _merge_bus_flows(pm, pd_bus_unmerged, connections)
    var(pm, nw, :qd_bus)[id] = qd_bus = _merge_bus_flows(pm, qd_bus_unmerged, connections)

    if report
        sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        sol(pm, nw, :load, id)[:pd] = pd
        sol(pm, nw, :load, id)[:qd] = qd
    end
end


"""
	function constraint_mc_load_power_delta(
		pm::AbstractExplicitNeutralACRModel,
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

For ACR models with explicit neutrals,
creates non-linear expressions for terminal power flows ':pd_bus' and ':qd_bus' of delta-connected loads
"""
function constraint_mc_load_power_delta(pm::AbstractExplicitNeutralACRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    nph = length(a)
    @assert nph == 3 "only phases == 3 delta loads are currently supported"
    next(i)  = i%nph+1 # 1=>2, 2=>3, 3=>1
    prev(i)  = (i-2+nph)%nph+1  # 1=>3, 2=>1, 3=>2

    vrd = [vr[t]-vr[connections[next(idx)]] for (idx,t) in enumerate(connections)]
    vid = [vi[t]-vi[connections[next(idx)]] for (idx,t) in enumerate(connections)]

    crd = JuMP.@NLexpression(pm.model, [idx in 1:nph], a[idx]*vrd[idx]*(vrd[idx]^2+vid[idx]^2)^(alpha[idx]/2-1)+b[idx]*vid[idx]*(vrd[idx]^2+vid[idx]^2)^(beta[idx]/2 -1))
    cid = JuMP.@NLexpression(pm.model, [idx in 1:nph], a[idx]*vid[idx]*(vrd[idx]^2+vid[idx]^2)^(alpha[idx]/2-1)-b[idx]*vrd[idx]*(vrd[idx]^2+vid[idx]^2)^(beta[idx]/2 -1))

    crd_bus_unmerged = JuMP.@NLexpression(pm.model, [idx in 1:nph], crd[idx]-crd[prev(idx)])
    cid_bus_unmerged = JuMP.@NLexpression(pm.model, [idx in 1:nph], cid[idx]-cid[prev(idx)])

    pd_bus_unmerged = [JuMP.@NLexpression(pm.model,  vr[p]*crd_bus_unmerged[idx]+vi[p]*cid_bus_unmerged[idx]) for (idx,p) in enumerate(connections)]
    qd_bus_unmerged = [JuMP.@NLexpression(pm.model, -vr[p]*cid_bus_unmerged[idx]+vi[p]*crd_bus_unmerged[idx]) for (idx,p) in enumerate(connections)]

    var(pm, nw, :pd_bus)[id] = pd_bus = _merge_bus_flows(pm, pd_bus_unmerged, connections)
    var(pm, nw, :qd_bus)[id] = qd_bus = _merge_bus_flows(pm, qd_bus_unmerged, connections)

    if report
        sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        pd = Vector{JuMP.NonlinearExpression}([])
        qd = Vector{JuMP.NonlinearExpression}([])
        for idx in 1:nph
            push!(pd, JuMP.@NLexpression(pm.model, a[idx]*(vrd[idx]^2+vid[idx]^2)^(alpha[idx]/2) ))
            push!(qd, JuMP.@NLexpression(pm.model, b[idx]*(vrd[idx]^2+vid[idx]^2)^(beta[idx]/2)  ))
        end
        sol(pm, nw, :load, id)[:pd] = JuMP.Containers.DenseAxisArray(pd, connections)
        sol(pm, nw, :load, id)[:qd] = JuMP.Containers.DenseAxisArray(qd, connections)
    end
end


# TRANSFORMER - Variable

"""
	function variable_mc_transformer_power(
		pm::AbstractExplicitNeutralACRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
		kwargs...
	)

For ACR models with explicit neutrals,
creates transfomer power variables `:pt` and `:qt`, and placeholder dictionaries for transformer terminal power flows `:pt_bus` and `:qt_bus`
"""
function variable_mc_transformer_power(pm::AbstractExplicitNeutralACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_transformer_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_transformer_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    var(pm, nw)[:pt_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:qt_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end

# TRANSFORMER - Constraints

"""
	function constraint_mc_transformer_power_yy(
		pm::AbstractExplicitNeutralACRModel,
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

For ACR models with explicit neutrals,
links the from-side and to-side power variables of wye-wye connected transformers.
Expressions for the terminal power flow variables are also added.
"""
function constraint_mc_transformer_power_yy(pm::AbstractExplicitNeutralACRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    pt_fr = var(pm, nw, :pt, f_idx)
    pt_to = var(pm, nw, :pt, t_idx)
    qt_fr = var(pm, nw, :qt, f_idx)
    qt_to = var(pm, nw, :qt, t_idx)
    vr_fr = var(pm, nw, :vr, f_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    JuMP.@constraint(pm.model, pt_fr + pt_to .== 0)
    JuMP.@constraint(pm.model, qt_fr + qt_to .== 0)

    f_phases = f_connections[1:end-1]
    f_n      = f_connections[end]
    t_phases = t_connections[1:end-1]
    t_n      = t_connections[end]
    P        = length(f_phases)

    # from-side
    if iszero(vr_fr[f_n]) && iszero(vi_fr[f_n])
        pt_bus_fr_unmerged = [pt_fr..., 0.0]
        qt_bus_fr_unmerged = [qt_fr..., 0.0]
    else
        vr_fr_pn = [vr_fr[p]-vr_fr[f_n] for p in f_phases]
        vi_fr_pn = [vi_fr[p]-vi_fr[f_n] for p in f_phases]

        crt_fr = [JuMP.@NLexpression(pm.model, ( pt_fr[idx]*vr_fr_pn[idx]+qt_fr[idx]*vi_fr_pn[idx])/(vr_fr_pn[idx]^2+vi_fr_pn[idx]^2)) for idx in 1:P]
        cit_fr = [JuMP.@NLexpression(pm.model, (-pt_fr[idx]*vi_fr_pn[idx]+qt_fr[idx]*vr_fr_pn[idx])/(vr_fr_pn[idx]^2+vi_fr_pn[idx]^2)) for idx in 1:P]

        pt_bus_fr_unmerged = [
            [JuMP.@NLexpression(pm.model, vr_fr[p]*crt_fr[idx]+vi_fr[p]*cit_fr[idx]) for (idx,p) in enumerate(f_phases)]...,
            JuMP.@NLexpression(pm.model, vr_fr[f_n]*sum(-crt_fr[idx] for idx in 1:P)+vi_fr[f_n]*sum(-cit_fr[idx] for idx in 1:P))
        ]
        qt_bus_fr_unmerged = [
            [JuMP.@NLexpression(pm.model, -vr_fr[p]*cit_fr[idx]+vi_fr[p]*crt_fr[idx]) for (idx,p) in enumerate(f_phases)]...,
            JuMP.@NLexpression(pm.model, -vr_fr[f_n]*sum(-cit_fr[idx] for idx in 1:P)+vi_fr[f_n]*sum(-crt_fr[idx] for idx in 1:P))
        ]
    end
    var(pm, nw, :pt_bus)[f_idx] = pt_bus_fr = _merge_bus_flows(pm, pt_bus_fr_unmerged, f_connections)
    var(pm, nw, :qt_bus)[f_idx] = qt_bus_fr = _merge_bus_flows(pm, qt_bus_fr_unmerged, f_connections)

    # to-side
    if iszero(vr_to[t_n]) && iszero(vi_to[t_n])
        pt_bus_to_unmerged = [pt_to..., 0.0]
        qt_bus_to_unmerged = [qt_to..., 0.0]
    else
        vr_to_pn = [vr_to[p]-vr_to[t_n] for p in t_phases]
        vi_to_pn = [vi_to[p]-vi_to[t_n] for p in t_phases]

        crt_to = [JuMP.@NLexpression(pm.model, ( pt_to[idx]*vr_to_pn[idx]+qt_to[idx]*vi_to_pn[idx])/(vr_to_pn[idx]^2+vi_to_pn[idx]^2)) for idx in 1:P]
        cit_to = [JuMP.@NLexpression(pm.model, (-pt_to[idx]*vi_to_pn[idx]+qt_to[idx]*vr_to_pn[idx])/(vr_to_pn[idx]^2+vi_to_pn[idx]^2)) for idx in 1:P]

        pt_bus_to_unmerged = [
            [JuMP.@NLexpression(pm.model, vr_to[p]*crt_to[idx]+vi_to[p]*cit_to[idx]) for (idx,p) in enumerate(t_phases)]...,
            JuMP.@NLexpression(pm.model, vr_to[t_n]*sum(-crt_to[idx] for idx in 1:P)+vi_to[t_n]*sum(-cit_to[idx] for idx in 1:P))
        ]
        qt_bus_to_unmerged = [
            [JuMP.@NLexpression(pm.model, -vr_to[p]*cit_to[idx]+vi_to[p]*crt_to[idx]) for (idx,p) in enumerate(t_phases)]...,
            JuMP.@NLexpression(pm.model, -vr_to[t_n]*sum(-cit_to[idx] for idx in 1:P)+vi_to[t_n]*sum(-crt_to[idx] for idx in 1:P))
        ]
    end
    var(pm, nw, :pt_bus)[t_idx] = pt_bus_to = _merge_bus_flows(pm, pt_bus_to_unmerged, t_connections)
    var(pm, nw, :qt_bus)[t_idx] = qt_bus_to = _merge_bus_flows(pm, qt_bus_to_unmerged, t_connections)
end

"""
	function constraint_mc_transformer_power_dy(
		pm::AbstractExplicitNeutralACRModel,
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

For ACR models with explicit neutrals,
links the from-side and to-side power variables of delta-wye connected transformers.
Expressions for the terminal power flow variables are also added.
"""
function constraint_mc_transformer_power_dy(pm::AbstractExplicitNeutralACRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    pt_fr = var(pm, nw, :pt, f_idx)
    pt_to = var(pm, nw, :pt, t_idx)
    qt_fr = var(pm, nw, :qt, f_idx)
    qt_to = var(pm, nw, :qt, t_idx)
    vr_fr = var(pm, nw, :vr, f_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    JuMP.@constraint(pm.model, pt_fr + pt_to .== 0)
    JuMP.@constraint(pm.model, qt_fr + qt_to .== 0)

    t_phases = t_connections[1:end-1]
    t_n      = t_connections[end]
    P        = length(t_phases)
    next(i)  = i%P+1 # 1=>2, 2=>3, 3=>1
    prev(i)  = (i-2+P)%P+1  # 1=>3, 2=>1, 3=>2

    # from-side
    vrg_fr = [vr_fr[p]-vr_fr[f_connections[next(idx)]] for (idx,p) in enumerate(f_connections)]
    vig_fr = [vi_fr[p]-vi_fr[f_connections[next(idx)]] for (idx,p) in enumerate(f_connections)]

    crt_fr = JuMP.@NLexpression(pm.model, [idx in 1:P], (pt_fr[idx]*vrg_fr[idx]+qt_fr[idx]*vig_fr[idx])/(vrg_fr[idx]^2+vig_fr[idx]^2))
    cit_fr = JuMP.@NLexpression(pm.model, [idx in 1:P], (pt_fr[idx]*vig_fr[idx]-qt_fr[idx]*vrg_fr[idx])/(vrg_fr[idx]^2+vig_fr[idx]^2))

    crt_bus_fr_unmerged = JuMP.@NLexpression(pm.model, [idx in 1:P], crt_fr[idx]-crt_fr[prev(idx)])
    cit_bus_fr_unmerged = JuMP.@NLexpression(pm.model, [idx in 1:P], cit_fr[idx]-cit_fr[prev(idx)])

    pt_bus_fr_unmerged = [JuMP.@NLexpression(pm.model,  vr_fr[p]*crt_bus_fr_unmerged[idx]+vi_fr[p]*cit_bus_fr_unmerged[idx]) for (idx,p) in enumerate(f_connections)]
    qt_bus_fr_unmerged = [JuMP.@NLexpression(pm.model, -vr_fr[p]*cit_bus_fr_unmerged[idx]+vi_fr[p]*crt_bus_fr_unmerged[idx]) for (idx,p) in enumerate(f_connections)]

    var(pm, nw, :pt_bus)[f_idx] = pt_bus_fr = _merge_bus_flows(pm, pt_bus_fr_unmerged, f_connections)
    var(pm, nw, :qt_bus)[f_idx] = qt_bus_fr = _merge_bus_flows(pm, qt_bus_fr_unmerged, f_connections)

    # to-side
    if iszero(vr_to[t_n]) && iszero(vi_to[t_n])
        pt_bus_to_unmerged = [pt_to..., 0.0]
        qt_bus_to_unmerged = [qt_to..., 0.0]
    else
        vr_to_pn = [vr_to[p]-vr_to[t_n] for p in t_phases]
        vi_to_pn = [vi_to[p]-vi_to[t_n] for p in t_phases]

        crt_to = [JuMP.@NLexpression(pm.model, ( pt_to[idx]*vr_to_pn[idx]+qt_to[idx]*vi_to_pn[idx])/(vr_to_pn[idx]^2+vi_to_pn[idx]^2)) for idx in 1:P]
        cit_to = [JuMP.@NLexpression(pm.model, (-pt_to[idx]*vi_to_pn[idx]+qt_to[idx]*vr_to_pn[idx])/(vr_to_pn[idx]^2+vi_to_pn[idx]^2)) for idx in 1:P]

        pt_bus_to_unmerged = [
            [JuMP.@NLexpression(pm.model, vr_to[p]*crt_to[idx]+vi_to[p]*cit_to[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, vr_to[n]*sum(-crt_to[idx] for idx in 1:P)+vi_to[n]*sum(-cit_to[idx] for idx in 1:P))
        ]
        qt_bus_to_unmerged = [
            [JuMP.@NLexpression(pm.model, -vr_to[p]*cit_to[idx]+vi_to[p]*crt_to[idx]) for (idx,p) in enumerate(phases)]...,
            JuMP.@NLexpression(pm.model, -vr_to[n]*sum(-cit_to[idx] for idx in 1:P)+vi_to[n]*sum(-crt_to[idx] for idx in 1:P))
        ]
    end
    var(pm, nw, :pt_bus)[t_idx] = pt_bus_to = _merge_bus_flows(pm, pt_bus_to_unmerged, t_connections)
    var(pm, nw, :qt_bus)[t_idx] = qt_bus_to = _merge_bus_flows(pm, qt_bus_to_unmerged, t_connections)
end


"""
	function constraint_mc_transformer_thermal_limit(
		pm::AbstractExplicitNeutralACRModel,
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

For ACR models with explicit neutrals,
imposes a bound on the magnitude of the total apparent power at each winding of the transformer.

```
sum(pt_fr)^2 + sum(qt_fr)^2 <= sm_ub^2
sum(pt_to)^2 + sum(qt_to)^2 <= sm_ub^2
```
"""
function constraint_mc_transformer_thermal_limit(pm::AbstractExplicitNeutralACRModel, nw::Int, id::Int, f_idx::Tuple, t_idx::Tuple, f_bus::Int, t_bus::Int, f_connections::Vector, t_connections::Vector, config::ConnConfig, sm_ub::Real; report::Bool=true)
    pt_fr = var(pm, nw, :pt, f_idx)
    qt_fr = var(pm, nw, :qt, f_idx)
    pt_to = var(pm, nw, :pt, t_idx)
    qt_to = var(pm, nw, :qt, t_idx)

    if sm_ub < Inf
        JuMP.@constraint(pm.model, sum(pt_fr)^2 + sum(qt_fr)^2 <= sm_ub^2)
        JuMP.@constraint(pm.model, sum(pt_to)^2 + sum(qt_to)^2 <= sm_ub^2)
    end
end


# BRANCH

# BRANCH - Variables

"""
	function variable_mc_branch_power(
		pm::AbstractExplicitNeutralACRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
		kwargs...
	)

For ACR models with explicit neutrals,
creates branch power variables `:p` and `:q` and placeholder dictionaries for the terminal power flows `:p_bus` and `:q_bus`.
"""
function variable_mc_branch_power(pm::AbstractExplicitNeutralACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_branch_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_branch_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    var(pm, nw)[:p_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:q_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end


# BRANCH - Constraints

"""
	function constraint_mc_ohms_yt_from(
		pm::AbstractExplicitNeutralACRModel,
		nw::Int,
		f_bus::Int,
		t_bus::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int},
		G::Matrix{<:Real},
		B::Matrix{<:Real},
		G_fr::Matrix{<:Real},
		B_fr::Matrix{<:Real}
	)

For ACR models with explicit neutrals,
creates Ohms constraints for ACR models with explicit neutrals.

```
s_fr = v_fr.*conj(Y*(v_fr-v_to))
s_fr = (vr_fr+im*vi_fr).*(G-im*B)*([vr_fr-vr_to]-im*[vi_fr-vi_to])
s_fr = (vr_fr+im*vi_fr).*([G*vr_fr-G*vr_to-B*vi_fr+B*vi_to]-im*[G*vi_fr-G*vi_to+B*vr_fr-B*vr_to])
```
"""
function constraint_mc_ohms_yt_from(pm::AbstractExplicitNeutralACRModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real})
    p_fr  = var(pm, nw, :p, f_idx)
    q_fr  = var(pm, nw, :q, f_idx)
    vr_fr = [var(pm, nw, :vr, f_bus)[t] for t in f_connections]
    vr_to = [var(pm, nw, :vr, t_bus)[t] for t in t_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[t] for t in f_connections]
    vi_to = [var(pm, nw, :vi, t_bus)[t] for t in t_connections]

    JuMP.@constraint(pm.model,
            p_fr .==  vr_fr.*(G*vr_fr-G*vr_to-B*vi_fr+B*vi_to)
                     +vi_fr.*(G*vi_fr-G*vi_to+B*vr_fr-B*vr_to)
                     # shunt
                     +vr_fr.*(G_fr*vr_fr-B_fr*vi_fr)
                     +vi_fr.*(G_fr*vi_fr+B_fr*vr_fr)
    )
    JuMP.@constraint(pm.model,
            q_fr .== -vr_fr.*(G*vi_fr-G*vi_to+B*vr_fr-B*vr_to)
                     +vi_fr.*(G*vr_fr-G*vr_to-B*vi_fr+B*vi_to)
                     # shunt
                     -vr_fr.*(G_fr*vi_fr+B_fr*vr_fr)
                     +vi_fr.*(G_fr*vr_fr-B_fr*vi_fr)
    )

    var(pm, nw, :p_bus)[f_idx] = _merge_bus_flows(pm, p_fr, f_connections)
    var(pm, nw, :q_bus)[f_idx] = _merge_bus_flows(pm, q_fr, f_connections)
end


"""
	function constraint_mc_ohms_yt_to(
		pm::AbstractExplicitNeutralACRModel,
		nw::Int,
		f_bus::Int,
		t_bus::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int},
		G::Matrix,
		B::Matrix,
		G_to::Matrix,
		B_to::Matrix
	)

For ACR models with explicit neutrals,
creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form).
```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_mc_ohms_yt_to(pm::AbstractExplicitNeutralACRModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix, B::Matrix, G_to::Matrix, B_to::Matrix)
    constraint_mc_ohms_yt_from(pm, nw, t_bus, f_bus, t_idx, f_idx, t_connections, f_connections, G, B, G_to, B_to)
end


"""
	function constraint_mc_branch_current_limit(
		pm::AbstractExplicitNeutralACRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector,
		t_connections::Vector,
		c_rating::Vector{<:Real};
		report::Bool=true
	)

For ACR models with explicit neutrals,
imposes a bound on the total current magnitude per conductor.

```
p_fr^2 + q_fr^2 <= r^2 * (vr_fr^2 + vi_fr^2)
p_to^2 + q_to^2 <= r^2 * (vr_to^2 + vi_to^2)
```
"""
function constraint_mc_branch_current_limit(pm::AbstractExplicitNeutralACRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int},  f_connections::Vector, t_connections::Vector, c_rating::Vector{<:Real}; report::Bool=true)
    p_fr = var(pm, nw, :p, f_idx)
    q_fr = var(pm, nw, :q, f_idx)
    p_to = var(pm, nw, :p, t_idx)
    q_to = var(pm, nw, :q, t_idx)
    vr_fr   = var(pm, nw, :vr, f_idx[2])
    vi_fr   = var(pm, nw, :vi, f_idx[2])
    vr_to   = var(pm, nw, :vr, t_idx[2])
    vi_to   = var(pm, nw, :vi, t_idx[2])

    for (idx, f_terminal, t_terminal, rating) in zip(1:length(f_connections), f_connections, t_connections, c_rating)
        if rating<Inf
            JuMP.@constraint(pm.model, p_fr[idx]^2+q_fr[idx]^2 <= rating^2*(vr_fr[f_terminal]^2+vi_fr[f_terminal]^2))
            JuMP.@constraint(pm.model, p_to[idx]^2+q_to[idx]^2 <= rating^2*(vr_to[t_terminal]^2+vi_to[t_terminal]^2))
        end
    end
end


"""
	function constraint_mc_thermal_limit_from(
		pm::AbstractExplicitNeutralACRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		rate_a::Vector{<:Real}
	)

For ACR models with explicit neutrals,
imposes a bound on the from-side line power magnitude.
"""
function constraint_mc_thermal_limit_from(pm::AbstractExplicitNeutralACRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})
    p_fr = var(pm, nw, :p, f_idx)
    q_fr = var(pm, nw, :q, f_idx)

    for idx in 1:length(rate_a)
        if rate_a[idx]<Inf
            # for branch-reduced models, p_fr and q_fr can be a sum of several terms
            # therefore, use a NLconstraint to handle these cases as well
            JuMP.@NLconstraint(pm.model, p_fr[idx]^2 + q_fr[idx]^2 <= rate_a[idx]^2)
        end
    end
end


"""
	function constraint_mc_thermal_limit_to(
		pm::AbstractExplicitNeutralACRModel,
		nw::Int,
		t_idx::Tuple{Int,Int,Int},
		t_connections::Vector{Int},
		rate_a::Vector{<:Real}
	)

For ACR models with explicit neutrals,
imposes a bound on the from-side line power magnitude.
"""
function constraint_mc_thermal_limit_to(pm::AbstractExplicitNeutralACRModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})
    p_to = var(pm, nw, :p, t_idx)
    q_to = var(pm, nw, :q, t_idx)

    for idx in 1:length(rate_a)
        if rate_a[idx]<Inf
            # for branch-reduced models, p_fr and q_fr can be a sum of several terms
            # therefore, use a NLconstraint to handle these cases as well
            JuMP.@NLconstraint(pm.model, p_to[idx]^2 + q_to[idx]^2 <= rate_a[idx]^2)
        end
    end
end


# SWITCH

# SWITCH - Variables

"""
	function variable_mc_switch_power(
		pm::AbstractExplicitNeutralACRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
		kwargs...
	)

For ACR models with explicit neutrals,
creates switch power variables `:p` and `:q` and placeholder dictionaries for the terminal power flows `:ps_bus` and `:qs_bus`.
"""
function variable_mc_switch_power(pm::AbstractExplicitNeutralACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_switch_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_switch_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    var(pm, nw)[:psw_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:qsw_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end


# SWITCH - Constraints

"""
	function constraint_mc_switch_power(
		pm::AbstractExplicitNeutralACRModel,
		nw::Int,
		id::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int};
		report::Bool=true
	)

constraint_mc_switch_power(
    pm::ReducedExplicitNeutralIVRModels,
    nw::Int,
    id::Int,
    f_idx::Tuple{Int,Int,Int},
    t_idx::Tuple{Int,Int,Int},
    f_connections::Vector{Int},
    t_connections::Vector{Int};
    report::Bool=true
)

For IVR models with explicit neutrals,
create expressions for the terminal power flows `:psw_bus` and `qsw_bus`,
and link the from-side to the to-side switch power
"""
function constraint_mc_switch_power(pm::AbstractExplicitNeutralACRModel, nw::Int, id::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}; report::Bool=true)
    psw_fr = var(pm, nw, :psw, f_idx)
    qsw_fr = var(pm, nw, :qsw, f_idx)
    psw_to = var(pm, nw, :psw, t_idx)
    qsw_to = var(pm, nw, :qsw, t_idx)

    JuMP.@constraint(pm.model, psw_fr .+ psw_to .== 0)
    JuMP.@constraint(pm.model, qsw_fr .+ qsw_to .== 0)

    var(pm, nw, :psw_bus)[f_idx] = _merge_bus_flows(pm, psw_fr, f_connections)
    var(pm, nw, :qsw_bus)[f_idx] = _merge_bus_flows(pm, qsw_fr, f_connections)
    var(pm, nw, :psw_bus)[t_idx] = _merge_bus_flows(pm, psw_to, t_connections)
    var(pm, nw, :qsw_bus)[t_idx] = _merge_bus_flows(pm, qsw_to, t_connections)
end


"""
	function constraint_mc_switch_current_limit(
		pm::AbstractExplicitNeutralACRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		rating::Vector{<:Real}
	)

For ACR models with explicit neutrals,
imposes a bound on the switch current magnitude per conductor.
Note that a bound on the from-side implies the same bound on the to-side current,
so it suffices to apply this only explicitly at the from-side.
"""
function constraint_mc_switch_current_limit(pm::AbstractExplicitNeutralACRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{<:Real})::Nothing
    vr_fr = [var(pm, nw, :vr, f_idx[2])[t] for t in f_connections]
    vi_fr = [var(pm, nw, :vi, f_idx[2])[t] for t in f_connections]
    psw = var(pm, nw, :psw, f_idx)
    qsw = var(pm, nw, :qsw, f_idx)

    mu_cm_fr = JuMP.ConstraintRef[]
    for idx in 1:length(rating)
        if rating[idx] < Inf
            push!(mu_cm_fr, JuMP.@constraint(pm.model, psw[idx]^2 + qsw[idx]^2 <= rating[idx]^2 * (vr_fr[idx]^2 + vi_fr[idx]^2)))
        end
    end

    con(pm, nw, :mu_cm_switch)[f_idx] = mu_cm_fr

    nothing
end


"""
	function constraint_mc_switch_thermal_limit(
		pm::AbstractExplicitNeutralACRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		rating::Vector{<:Real}
	)

For ACR models with explicit neutrals,
imposes a bound on the switch power magnitude per conductor.
Note that a bound on the from-side implies the same bound on the to-side power
when the switch is closed (equal voltages), and also when it is open since the
power then equals zero on both ends.
"""
function constraint_mc_switch_thermal_limit(pm::AbstractExplicitNeutralACRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{<:Real})::Nothing
    psw_fr = var(pm, nw, :psw, f_idx)
    qsw_fr = var(pm, nw, :qsw, f_idx)

    mu_sm_fr = JuMP.ConstraintRef[]

    for idx in 1:length(rating)
        if rating[idx] < Inf
            push!(mu_sm_fr, JuMP.@constraint(pm.model, psw_fr[idx]^2 + qsw_fr[idx]^2 <= rating[idx]^2))
        end
    end

    con(pm, nw, :mu_sm_switch)[f_idx] = mu_sm_fr

    nothing
end
