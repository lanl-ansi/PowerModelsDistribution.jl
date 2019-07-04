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
    _PMs.var(pm, nw)[:Wr] = Wre
    _PMs.var(pm, nw)[:Wi] = Wim
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
    _PMs.var(pm, nw)[:P] = P
    _PMs.var(pm, nw)[:Q] = Q
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
    p_to = _PMs.var(pm, n, :P)[t_idx]
    q_to = _PMs.var(pm, n, :Q)[t_idx]

    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    w_to_re = _PMs.var(pm, n, :Wr)[t_bus]
    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]

    w_to_im = _PMs.var(pm, n, :Wi)[t_bus]
    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]

    ccm_re =  _PMs.var(pm, n, :CC_re)[i]
    ccm_im =  _PMs.var(pm, n, :CC_im)[i]

    JuMP.@constraint(pm.model, p_fr + p_to .==  w_fr_re*(g_sh_fr)' + w_fr_im*(b_sh_fr)' + r*ccm_re - x*ccm_im +  w_to_re*(g_sh_to)'  + w_to_im*(b_sh_to)')
    JuMP.@constraint(pm.model, q_fr + q_to .==  w_fr_im*(g_sh_fr)' - w_fr_re*(b_sh_fr)' + x*ccm_re + r*ccm_im +  w_to_im*(g_sh_to)'  - w_to_re*(b_sh_to)')
end


""
function constraint_tp_theta_ref(pm::_PMs.GenericPowerModel{T}, n::Int, i) where T <: AbstractUBFForm
    nconductors = length(_PMs.conductor_ids(pm))

    w_re = _PMs.var(pm, n, :Wr)[i]
    w_im = _PMs.var(pm, n, :Wi)[i]

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
    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]

    w_to_re = _PMs.var(pm, n, :Wr)[t_bus]
    w_to_im = _PMs.var(pm, n, :Wi)[t_bus]

    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

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


function variable_tp_generation_power_mx(pm::_PMs.GenericPowerModel; nw=pm.cnw)
    gen_ids = collect(_PMs.ids(pm, nw, :gen))
    # delegate creation of diagonal elements back to PMs as before
    for id in gen_ids, c in _PMs.conductor_ids(pm, nw)
        _PMs.variable_generation(pm, nw=nw, cnd=c)
    end
    # create vectors of the diagonal elements
    ncnds = length(_PMs.conductor_ids(pm, nw))
    diag_re = Dict([(id, [_PMs.var(pm, nw, c, :pg, id) for c in 1:ncnds]) for id in gen_ids])
    diag_im = Dict([(id, [_PMs.var(pm, nw, c, :qg, id) for c in 1:ncnds]) for id in gen_ids])
    # calculate bounds for matrix variable
    bound = Dict{eltype(gen_ids), Array{Real,2}}()
    for (id, gen) in _PMs.ref(pm, nw, :gen)
        bus_id = gen["gen_bus"]
        vmax = _PMs.ref(pm, nw, :bus, bus_id)["vmax"].values
        vmin = _PMs.ref(pm, nw, :bus, bus_id)["vmin"].values
        pmax = max(abs.(gen["pmax"].values), abs.(gen["pmin"].values))
        qmax = max(abs.(gen["qmax"].values), abs.(gen["qmin"].values))
        smax = sqrt.(pmax.^2 + qmax.^2)
        cmax = smax./vmin
        bound[id] = vmax*cmax'
    end
    # create matrix variables, whilst injecting diagonals
    (Pg,Qg) = variable_mx_complex_with_diag(pm.model, gen_ids, ncnds, bound;
        diag_re=diag_re, diag_im=diag_im, name=("Pg", "Qg"), prefix="$nw")
    # save references
    _PMs.var(pm, nw)[:Pg] = Pg
    _PMs.var(pm, nw)[:Qg] = Qg
end


function variable_tp_generation_current_mx(pm::_PMs.GenericPowerModel; nw=pm.cnw)
    gen_ids = collect(_PMs.ids(pm, nw, :gen))
    ncnds = length(_PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(gen_ids), Array{Real,2}}()
    for (id, gen) in _PMs.ref(pm, nw, :gen)
        bus_id = gen["gen_bus"]
        vmax = _PMs.ref(pm, nw, :bus, bus_id)["vmax"].values
        vmin = _PMs.ref(pm, nw, :bus, bus_id)["vmin"].values
        pmax = max(abs.(gen["pmax"].values), abs.(gen["pmin"].values))
        qmax = max(abs.(gen["qmax"].values), abs.(gen["qmin"].values))
        smax = sqrt.(pmax.^2 + qmax.^2)
        cmax = smax./vmin
        bound[id] = cmax*cmax'
    end
    # create matrix variables
    (Lre,Lim) = variable_mx_hermitian(pm.model, gen_ids, ncnds, bound; name="Ld", prefix="$nw")
    # save references
    _PMs.var(pm, nw)[:CCgr] = Lre
    _PMs.var(pm, nw)[:CCgi] = Lim
end


function variable_tp_load_mx(pm::_PMs.GenericPowerModel; nw=pm.cnw)
    load_wye_ids = [id for (id, load) in _PMs.ref(pm, nw, :load) if load["conn"]=="wye"]
    load_del_ids = [id for (id, load) in _PMs.ref(pm, nw, :load) if load["conn"]=="delta"]
    # first create the wye loads; will create keys :Pd, :Qd
    variable_tp_load_power_wye_mx(pm, load_wye_ids)
    # now, create delta loads; will create :Xdre, :Xdim
    variable_tp_load_power_delta(pm, load_del_ids)
    # define :Pd, :Qd for delta loads as lin. transformation of :Xdre and :Wdim
    Td = [1 -1 0; 0 1 -1; -1 0 1]
    for id in load_del_ids
        _PMs.var(pm, nw, :Pd)[id] = _PMs.var(pm, nw, :Xdre, id)*Td
        _PMs.var(pm, nw, :Qd)[id] = _PMs.var(pm, nw, :Xdim, id)*Td
    end
    # both loads need a current variable; bounds adjusted for connection type
    variable_tp_load_current_mx(pm, [load_wye_ids..., load_del_ids...])
end


function variable_tp_load(pm::_PMs.GenericPowerModel{T}; nw=pm.cnw) where T <: AbstractUBFForm
    load_wye_ids = [id for (id, load) in _PMs.ref(pm, nw, :load) if load["conn"]=="wye"]
    load_del_ids = [id for (id, load) in _PMs.ref(pm, nw, :load) if load["conn"]=="delta"]
    # create dictionary for wye loads
    for c in _PMs.conductor_ids(pm)
        variable_load(pm, nw=nw, cnd=c)
    end
    # now, create delta loads; will create :Xdre, :Xdim
    variable_tp_load_power_delta(pm, load_del_ids)
    # define :Pd, :Qd for delta loads as lin. transformation of :Xdre and :Wdim
    Td = [1 -1 0; 0 1 -1; -1 0 1]
    for id in load_del_ids
        Pd = _PMs.var(pm, nw, :Xdre, id)*Td
        Qd = _PMs.var(pm, nw, :Xdim, id)*Td
        for c in _PMs.conductor_ids(pm)
            _PMs.var(pm, nw, c, :pd)[id] = Pd[c,c]
            _PMs.var(pm, nw, c, :qd)[id] = Qd[c,c]
        end
    end
    # only delta loads need a current matrix variable
    variable_tp_load_current_mx(pm, load_del_ids)
end


function variable_tp_load(pm::_PMs.GenericPowerModel{T}; nw=pm.cnw) where T <: LPLinUBFForm
    load_wye_ids = [id for (id, load) in _PMs.ref(pm, nw, :load) if load["conn"]=="wye"]
    load_del_ids = [id for (id, load) in _PMs.ref(pm, nw, :load) if load["conn"]=="delta"]
    # create dictionary for wye loads
    for c in _PMs.conductor_ids(pm)
        variable_load(pm, nw=nw, cnd=c)
    end
    #TODO figure out how to include delta loads in the LPLinUBFForm
end


function variable_tp_load_power_wye_mx(pm::_PMs.GenericPowerModel, load_ids::Array{Int,1}; nw=pm.cnw)
    ncnds = length(_PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(load_ids), Array{Real,2}}()
    for id in load_ids
        load = _PMs.ref(pm, nw, :load, id)
        bus_id = load["load_bus"]
        vmax = _PMs.ref(pm, nw, :bus, bus_id)["vmax"].values
        vmin = _PMs.ref(pm, nw, :bus, bus_id)["vmin"].values
        pmax = abs.(load["pd"].values)
        qmax = abs.(load["qd"].values)
        smax = sqrt.(pmax.^2 + qmax.^2)
        cmax = smax./vmin
        bound[id] = vmax*cmax'
    end
    # create matrix variables
    (Pd,Qd) = variable_mx_complex_with_diag(pm.model, load_ids, ncnds, bound; name=("Pd", "Qd"), prefix="$nw")
    # save references
    _PMs.var(pm, nw)[:Pd] = Pd
    _PMs.var(pm, nw)[:Qd] = Qd
    for c in 1:ncnds
        _PMs.var(pm, nw, c)[:pd] =Dict([(id, Pd[id][c,c]) for id in load_ids])
        _PMs.var(pm, nw, c)[:qd] =Dict([(id, Qd[id][c,c]) for id in load_ids])
    end
end


"""
Creates power matrix variable X for delta windings; this defines both the
wye-side power Sy and the delta-side power Sd through the lin. transformations
Sy = X.Td, Sd = Td.X with Td=[1 -1 0; 0 1 -1; -1 0 1]

See the paper by Zhao et al. for the first convex relaxation of delta transformations.
@INPROCEEDINGS{zhao_optimal_2017,
	author={C. Zhao, E. Dall'Anese and S. Low},
	booktitle={IREP 2017 Bulk Power Systems Dynamics and Control Symposium},
	title={{Optimal Power Flow in Multiphase Radial Networks with Delta Connections}},
	year={2017},
	month={},
    url={https://www.nrel.gov/docs/fy18osti/67852.pdf}
}

See upcoming paper for discussion of bounds. [reference added when accepted]

Note: this does not have the mx suffix because it is needed for both vec and
mat KCL.
"""
function variable_tp_load_power_delta(pm::_PMs.GenericPowerModel, load_ids::Array{Int,1}; nw=pm.cnw, eps=0.1)
    ncnds = length(_PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(load_ids), Array{Real,2}}()
    for id in load_ids
        load = _PMs.ref(pm, nw, :load, id)
        # get voltage LL bounds
        bus_id = load["load_bus"]
        bus = _PMs.ref(pm, nw, :bus, bus_id)
        vdmax, vdmin = _bus_vm_ll_bounds(bus)
        pmax = abs.(load["pd"].values)
        qmax = abs.(load["qd"].values)
        smax = sqrt.(pmax.^2 + qmax.^2)
        cdmax = smax./vdmin
        bound[id] = bus["vmax"].values*cdmax'
    end
    # create matrix variables
    (Xdre,Xdim) = variable_mx_complex(pm.model, load_ids, ncnds, ncnds, bound; name="Xd", prefix="$nw")
    # save references
    _PMs.var(pm, nw)[:Xdre] = Xdre
    _PMs.var(pm, nw)[:Xdim] = Xdim
end


function variable_tp_load_current_mx(pm::_PMs.GenericPowerModel, load_ids::Array{Int,1}; nw=pm.cnw)
    ncnds = length(_PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(load_ids), Array{Real,2}}()
    for (id, load) in _PMs.ref(pm, nw, :load)
        bus_id = load["load_bus"]
        bus = _PMs.ref(pm, nw, :bus, bus_id)
        vmax = bus["vmax"].values
        vmin = bus["vmin"].values
        # this presumes constant power, wye loads!
        @assert(load["model"]=="constant_power")
        #TODO extend to other load models
        pmax = abs.(load["pd"].values)
        qmax = abs.(load["qd"].values)
        smax = sqrt.(pmax.^2 + qmax.^2)
        if load["conn"]=="wye"
            cmax = smax./vmin
        elseif load["conn"]=="delta"
            vdmax, vdmin = _bus_vm_ll_bounds(bus)
            cmax = smax./vdmin
        end
        bound[id] = cmax*cmax'
    end
    # create matrix variables
    (Ldre, Ldim) = variable_mx_hermitian(pm.model, load_ids, ncnds, bound; name="Ld", prefix="$nw")
    # save references
    _PMs.var(pm, nw)[:CCdr] = Ldre
    _PMs.var(pm, nw)[:CCdi] = Ldim
end


function constraint_tp_generation_mx(pm::_PMs.GenericPowerModel{T}, gen_id::Int; nw::Int=pm.cnw) where T <: AbstractUBFForm
    Pg = _PMs.var(pm, nw, :Pg, gen_id)
    Qg = _PMs.var(pm, nw, :Qg, gen_id)
    bus_id = _PMs.ref(pm, nw, :gen, gen_id)["gen_bus"]
    W_re = _PMs.var(pm, nw, :Wr, bus_id)
    W_im = _PMs.var(pm, nw, :Wi, bus_id)
    Lgre = _PMs.var(pm, nw, :CCgr, gen_id)
    Lgim = _PMs.var(pm, nw, :CCgi, gen_id)
    constraint_SWL_psd(pm.model, Pg, Qg, W_re, W_im, Lgre, Lgim)
end


function constraint_tp_load_mx(pm::_PMs.GenericPowerModel{T}, load_id::Int; nw::Int=pm.cnw) where T <: AbstractUBFForm
    # shared variables and parameters
    load = _PMs.ref(pm, nw, :load, load_id)
    pd = load["pd"].values
    qd = load["qd"].values
    bus_id = _PMs.ref(pm, nw, :load, load_id)["load_bus"]
    W_re = _PMs.var(pm, nw, :Wr, bus_id)
    W_im = _PMs.var(pm, nw, :Wi, bus_id)
    Ldre = _PMs.var(pm, nw, :CCdr, load_id)
    Ldim = _PMs.var(pm, nw, :CCdi, load_id)

    @assert(load["model"]=="constant_power")
    if load["conn"]=="wye"
        # set the diagonal values
        Pd = _PMs.var(pm, nw, :Pd, load_id)
        Qd = _PMs.var(pm, nw, :Qd, load_id)
        for c in 1:length(pd)
            Pd[c,c] = pd[c]
            Qd[c,c] = qd[c]
        end
        # link S, W and L
        constraint_SWL_psd(pm.model, Pd, Qd, W_re, W_im, Ldre, Ldim)
    elseif load["conn"]=="delta"
        Xdre = _PMs.var(pm, nw, :Xdre, load_id)
        Xdim = _PMs.var(pm, nw, :Xdim, load_id)
        Td = [1 -1 0; 0 1 -1; -1 0 1]
        JuMP.@constraint(pm.model, LinearAlgebra.diag(Td*Xdre) .== pd)
        JuMP.@constraint(pm.model, LinearAlgebra.diag(Td*Xdim) .== qd)
        constraint_SWL_psd(pm.model, Xdre, Xdim, W_re, W_im, Ldre, Ldim)
    end
end


function constraint_tp_load(pm::_PMs.GenericPowerModel{T}, load_id::Int; nw::Int=pm.cnw) where T <: AbstractUBFForm
    # shared variables and parameters
    load = _PMs.ref(pm, nw, :load, load_id)
    pd = load["pd"].values
    qd = load["qd"].values
    bus_id = _PMs.ref(pm, nw, :load, load_id)["load_bus"]

    @assert(load["model"]=="constant_power")
    if load["conn"]=="wye"
        for c in 1:length(pd)
            _PMs.var(pm, nw, c, :pd)[load_id] = pd[c]
            _PMs.var(pm, nw, c, :qd)[load_id] = qd[c]
        end
    elseif load["conn"]=="delta"
        W_re = _PMs.var(pm, nw, :Wr, bus_id)
        W_im = _PMs.var(pm, nw, :Wi, bus_id)
        Ldre = _PMs.var(pm, nw, :CCdr, load_id)
        Ldim = _PMs.var(pm, nw, :CCdi, load_id)
        Xdre = _PMs.var(pm, nw, :Xdre, load_id)
        Xdim = _PMs.var(pm, nw, :Xdim, load_id)
        Td = [1 -1 0; 0 1 -1; -1 0 1]
        JuMP.@constraint(pm.model, LinearAlgebra.diag(Td*Xdre) .== pd)
        JuMP.@constraint(pm.model, LinearAlgebra.diag(Td*Xdim) .== qd)
        constraint_SWL_psd(pm.model, Xdre, Xdim, W_re, W_im, Ldre, Ldim)
    end
end


function constraint_tp_voltage_psd(pm::_PMs.GenericPowerModel; nw=pm.cnw)
    buses_covered = [i for (l,i,j) in _PMs.ref(pm, nw, :arcs)]
    buses_psd = [i for i in _PMs.ids(pm, nw, :bus) if !(i in buses_covered)]
    for bus_id in buses_psd
        W_re = _PMs.var(pm, nw, :Wr, bus_id)
        W_im = _PMs.var(pm, nw, :Wi, bus_id)
        constraint_M_psd(W_re, W_im)
    end
end


function constraint_SWL_psd(model::JuMP.Model, P, Q, W_re, W_im, L_re, L_im)
    M_re = [W_re P; P' L_re]
    M_im = [W_im Q; -Q' L_im]
    constraint_M_psd(model, M_re, M_im)
end


function constraint_M_psd(model::JuMP.Model, M_re, M_im)
    JuMP.@constraint(model, [M_re -M_im; M_im M_re] in JuMP.PSDCone())
end


""
function constraint_tp_power_balance_mx_shunt(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    if !haskey(_PMs.con(pm, nw), :kcl_P)
        _PMs.con(pm, nw)[:kcl_P] = Dict{Int,Array{JuMP.ConstraintRef,2}}()
    end
    if !haskey(_PMs.con(pm, nw), :kcl_Q)
        _PMs.con(pm, nw)[:kcl_Q] = Dict{Int,Array{JuMP.ConstraintRef,2}}()
    end

    bus = _PMs.ref(pm, nw, :bus, i)
    bus_arcs = _PMs.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = _PMs.ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = _PMs.ref(pm, nw, :bus_gens, i)
    bus_loads = _PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PMs.ref(pm, nw, :bus_shunts, i)

    bus_Gs = Dict(k => LinearAlgebra.diagm(0=>_PMs.ref(pm, nw, :shunt, k, "gs").values) for k in bus_shunts)
    bus_Bs = Dict(k => LinearAlgebra.diagm(0=>_PMs.ref(pm, nw, :shunt, k, "bs").values) for k in bus_shunts)

    constraint_tp_power_balance_mx_shunt(pm, nw, i, bus_arcs, bus_arcs_dc, bus_gens, bus_loads, bus_Gs, bus_Bs)
end


"""
Shunt handling in matrix form:
I = Y.U
S = U.I' = U.(Y.U)' = U.U'.Y' = W.Y'
  = (Wre+j.Wim)(G+jB)' = (Wre+j.Wim)(G'-j.B') = (Wre.G'+Wim.B')+j(-Wre.B'+Wim.G')
P =  Wre.G'+Wim.B'
Q = -Wre.B'+Wim.G'
"""
function constraint_tp_power_balance_mx_shunt(pm::_PMs.GenericPowerModel{T}, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_loads, bus_Gs, bus_Bs) where T <: AbstractUBFForm
    Wre = _PMs.var(pm, n, :Wr, i)
    Wim = _PMs.var(pm, n, :Wi, i)
    P = _PMs.var(pm, n, :P)
    Q = _PMs.var(pm, n, :Q)
    Pg = _PMs.var(pm, n, :Pg)
    Qg = _PMs.var(pm, n, :Qg)
    Pd = _PMs.var(pm, n, :Pd)
    Qd = _PMs.var(pm, n, :Qd)
    # ignore dc for now
    #TODO add DC in matrix version?
    ncnds = size(Wre)[1]
    G = (length(bus_Gs)>0) ? sum(values(bus_Gs)) : zeros(ncnds, ncnds)
    B = (length(bus_Bs)>0) ? sum(values(bus_Bs)) : zeros(ncnds, ncnds)

    # changed the ordering
    # LHS: all variables with generator sign convention
    # RHS: all variables with load sign convention
    _PMs.con(pm, n, :kcl_P)[i] = JuMP.@constraint(pm.model, sum(Pg[g] for g in bus_gens) .== sum(P[a] for a in bus_arcs) + sum(Pd[d] for d in bus_loads) + ( Wre*G'+Wim*B'))
    _PMs.con(pm, n, :kcl_Q)[i] = JuMP.@constraint(pm.model, sum(Qg[g] for g in bus_gens) .== sum(Q[a] for a in bus_arcs) + sum(Qd[d] for d in bus_loads) + (-Wre*B'+Wim*G'))
end
