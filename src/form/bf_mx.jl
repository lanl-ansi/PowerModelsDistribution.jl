import LinearAlgebra: diag, diagm


""
function variable_tp_branch_current(pm::_PMs.GenericPowerModel{T}; kwargs...) where T <: AbstractUBFForm
    variable_tp_branch_series_current_prod_hermitian(pm; kwargs...)
end


""
function variable_tp_voltage(pm::_PMs.GenericPowerModel{T}; kwargs...) where T <: AbstractUBFForm
    variable_tp_voltage_prod_hermitian(pm; kwargs...)
end


""
function variable_tp_voltage_prod_hermitian(pm::_PMs.GenericPowerModel{T}; n_cond::Int=3, nw::Int=pm.cnw, bounded = true) where T <: AbstractUBFForm
    bus_ids = collect(_PMs.ids(pm, nw, :bus))

    if bounded
        # get bounds
        vmax = Dict([(id, _PMs.ref(pm, nw, :bus, id, "vmax").values) for id in bus_ids])
        vmin = Dict([(id, _PMs.ref(pm, nw, :bus, id, "vmin").values) for id in bus_ids])
        # create bounded Hermitian matrix variables
        (Wre,Wim) = variable_mx_hermitian_sqrt_bounds(pm.model, bus_ids, n_cond,
            vmax, vmin; name="W", prefix="$nw")
    else
        # create unbounded Hermitian matrix variables
        (Wre,Wim) = variable_mx_hermitian(pm.model, bus_ids, n_cond; name="W", prefix="$nw", lb_diag_zero=0)
    end

    # save references in dict
    _PMs.var(pm, nw)[:W_re] = Wre
    _PMs.var(pm, nw)[:W_im] = Wim
    for c in 1:n_cond
        _PMs.var(pm, nw, c)[:w] = Dict{Int, Any}([(id, Wre[id][c,c]) for id in bus_ids])
    end
end


""
function variable_tp_branch_series_current_prod_hermitian(pm::_PMs.GenericPowerModel{T}; n_cond::Int=3, nw::Int=pm.cnw, bounded = true) where T <: AbstractUBFForm
    branches = _PMs.ref(pm, nw, :branch)
    buses = _PMs.ref(pm, nw, :bus)

    branch_ids = collect(keys(branches))

    if bounded
        # calculate max series current for each branch
        cmax = Dict{eltype(branch_ids), Array{Real,1}}()
        for (key, branch) in branches
            bus_fr = buses[branch["f_bus"]]
            bus_to = buses[branch["t_bus"]]

            vmin_fr = bus_fr["vmin"].values
            vmin_to = bus_to["vmin"].values

            vmax_fr = bus_fr["vmax"].values
            vmax_to = bus_to["vmax"].values

            # assumed to be matrices already
            # temportary fix by shunts_diag2mat!
            y_fr = branch["g_fr"].values + im* branch["b_fr"].values
            y_to = branch["g_to"].values + im* branch["b_to"].values

            y_fr = diagm(0=>y_fr)
            y_to = diagm(0=>y_to)

            smax = branch["rate_a"].values
            cmaxfr = smax./vmin_fr + abs.(y_fr)*vmax_fr
            cmaxto = smax./vmin_to + abs.(y_to)*vmax_to

            cmax[key] = max.(cmaxfr, cmaxto)
        end
        # create matrix variables
        (Lre,Lim) = variable_mx_hermitian_sqrt_bounds(pm.model, branch_ids, n_cond, cmax; name="CC", prefix="$nw")
    else
        (Lre,Lim) = variable_mx_hermitian(pm.model, branch_ids, n_cond; name="CC", prefix="$nw", lb_diag_zero=true)
    end

    # save reference
    _PMs.var(pm, nw)[:CC_re] = Lre
    _PMs.var(pm, nw)[:CC_im] = Lim
    for c in 1:n_cond
        _PMs.var(pm, nw, c)[:cm] = Dict([(id, Lre[id][c,c]) for id in branch_ids])
    end
end


""
function variable_tp_branch_flow(pm::_PMs.GenericPowerModel{T}; n_cond::Int=3, nw::Int=pm.cnw, bounded = true) where T <: AbstractUBFForm
    # calculate S bound
    branch_arcs = _PMs.ref(pm, nw, :arcs)
    if bounded
        bound = Dict{eltype(branch_arcs), Array{Real,2}}()
        for (l,i,j) in branch_arcs
            vmin = _PMs.ref(pm, nw, :bus, i)["vmin"].values
            vmax = _PMs.ref(pm, nw, :bus, i)["vmax"].values
            smax = _PMs.ref(pm, nw, :branch, l)["rate_a"].values
            cmax = smax./vmin
            bound[(l,i,j)] = vmax*cmax'
            for c in 1:length(smax)
                bound[(l,i,j)][c,c] = smax[c]
            end
        end
        # create matrix variables
        (P,Q) = variable_mx_complex(pm.model, branch_arcs, n_cond, n_cond, bound; name=("P", "Q"), prefix="$nw")
    else
        (P,Q) = variable_mx_complex(pm.model, branch_arcs, n_cond, n_cond; name=("P", "Q"), prefix="$nw")
    end
    # save reference
    _PMs.var(pm, nw)[:P_mx] = P
    _PMs.var(pm, nw)[:Q_mx] = Q
    for c in 1:n_cond
        _PMs.var(pm, nw, c)[:p] = Dict([(id,P[id][c,c]) for id in branch_arcs])
        _PMs.var(pm, nw, c)[:q] = Dict([(id,Q[id][c,c]) for id in branch_arcs])
    end
end


""
function constraint_tp_theta_ref(pm::_PMs.GenericPowerModel{T}, i::Int; nw::Int=pm.cnw) where T <: AbstractUBFForm
    constraint_tp_theta_ref(pm, nw, i)
end


"Defines branch flow model power flow equations"
function constraint_tp_flow_losses(pm::_PMs.GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to) where T <: AbstractUBFForm
    p_to = _PMs.var(pm, n, :P_mx)[t_idx]
    q_to = _PMs.var(pm, n, :Q_mx)[t_idx]

    p_fr = _PMs.var(pm, n, :P_mx)[f_idx]
    q_fr = _PMs.var(pm, n, :Q_mx)[f_idx]

    w_to_re = _PMs.var(pm, n, :W_re)[t_bus]
    w_fr_re = _PMs.var(pm, n, :W_re)[f_bus]

    w_to_im = _PMs.var(pm, n, :W_im)[t_bus]
    w_fr_im = _PMs.var(pm, n, :W_im)[f_bus]

    ccm_re =  _PMs.var(pm, n, :CC_re)[i]
    ccm_im =  _PMs.var(pm, n, :CC_im)[i]

    JuMP.@constraint(pm.model, p_fr + p_to .==  w_fr_re*(g_sh_fr)' + w_fr_im*(b_sh_fr)' + r*ccm_re - x*ccm_im +  w_to_re*(g_sh_to)'  + w_to_im*(b_sh_to)')
    JuMP.@constraint(pm.model, q_fr + q_to .==  w_fr_im*(g_sh_fr)' - w_fr_re*(b_sh_fr)' + x*ccm_re + r*ccm_im +  w_to_im*(g_sh_to)'  - w_to_re*(b_sh_to)')
end


""
function constraint_tp_theta_ref(pm::_PMs.GenericPowerModel{T}, n::Int, i) where T <: AbstractUBFForm
    nconductors = length(_PMs.conductor_ids(pm))

    w_re = _PMs.var(pm, n, :W_re)[i]
    w_im = _PMs.var(pm, n, :W_im)[i]

    alpha = exp(-im*_wrap_to_pi(2 * pi / nconductors ))
    beta = (alpha*ones(nconductors)).^(0:nconductors-1)
    gamma = beta*beta'

    w_re_ref = real(gamma).*w_re[1,1]
    w_im_ref = imag(gamma).*w_re[1,1]
    JuMP.@constraint(pm.model, diag(w_re)[2:nconductors]        .== diag(w_re_ref)[2:nconductors]) # first equality is implied
    JuMP.@constraint(pm.model, _mat2utrivec!(w_re) .== _mat2utrivec!(w_re_ref))
    JuMP.@constraint(pm.model, _mat2utrivec!(w_im) .== _mat2utrivec!(w_im_ref))
end


"Defines voltage drop over a branch, linking from and to side voltage"
function constraint_tp_model_voltage_magnitude_difference(pm::_PMs.GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm) where T <: AbstractUBFForm
    w_fr_re = _PMs.var(pm, n, :W_re)[f_bus]
    w_fr_im = _PMs.var(pm, n, :W_im)[f_bus]

    w_to_re = _PMs.var(pm, n, :W_re)[t_bus]
    w_to_im = _PMs.var(pm, n, :W_im)[t_bus]

    p_fr = _PMs.var(pm, n, :P_mx)[f_idx]
    q_fr = _PMs.var(pm, n, :Q_mx)[f_idx]

    p_s_fr = p_fr - (w_fr_re*(g_sh_fr)' + w_fr_im*(b_sh_fr)')
    q_s_fr = q_fr - (w_fr_im*(g_sh_fr)' - w_fr_re*(b_sh_fr)')

    ccm_re =  _PMs.var(pm, n, :CC_re)[i]
    ccm_im =  _PMs.var(pm, n, :CC_im)[i]

    #KVL over the line:
    JuMP.@constraint(pm.model, diag(w_to_re) .== diag(
    w_fr_re   - p_s_fr  *r' - q_s_fr*x'        - r*p_s_fr'    - x*q_s_fr'
    + r*ccm_re*r' - x     *ccm_im*r' + x*ccm_re *x' + r*ccm_im *x'))
    JuMP.@constraint(pm.model, _mat2utrivec!(w_to_re) .== _mat2utrivec!(
    w_fr_re   - p_s_fr  *r' - q_s_fr*x'        - r*p_s_fr'    - x*q_s_fr'
    + r*ccm_re*r' - x     *ccm_im*r' + x*ccm_re *x' + r*ccm_im *x'))
    JuMP.@constraint(pm.model, _mat2utrivec!(w_to_im) .== _mat2utrivec!(
    w_fr_im   - q_s_fr  *r' + p_s_fr*x'        - x*p_s_fr'    + r*q_s_fr'
    + x*ccm_re*r' + r     *ccm_im*r' - r*ccm_re *x' + x*ccm_im *x'))
end
