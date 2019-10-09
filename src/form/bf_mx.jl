import LinearAlgebra: diag, diagm

""
function variable_mc_branch_current(pm::AbstractUBFModels; kwargs...)
    variable_mc_branch_series_current_prod_hermitian(pm; kwargs...)
end


""
function variable_mc_voltage(pm::AbstractUBFModels; kwargs...)
    variable_mc_voltage_prod_hermitian(pm; kwargs...)
end


""
function variable_mc_voltage_prod_hermitian(pm::AbstractUBFModels; n_cond::Int=3, nw::Int=pm.cnw, bounded = true)
    bus_ids = collect(_PMs.ids(pm, nw, :bus))

    if bounded
        # get bounds
        vmax = Dict([(id, _PMs.ref(pm, nw, :bus, id, "vmax").values) for id in bus_ids])
        vmin = Dict([(id, _PMs.ref(pm, nw, :bus, id, "vmin").values) for id in bus_ids])
        # create bounded Hermitian matrix variables
        (Wr,Wi) = variable_mx_hermitian_sqrt_bounds(pm.model, bus_ids, n_cond,
            vmax, vmin; name="W", prefix="$nw")
    else
        # create unbounded Hermitian matrix variables
        (Wr,Wi) = variable_mx_hermitian(pm.model, bus_ids, n_cond; name="W", prefix="$nw", lb_diag_zero=0)
    end

    # save references in dict
    _PMs.var(pm, nw)[:Wr] = Wr
    _PMs.var(pm, nw)[:Wi] = Wi
    for c in 1:n_cond
        _PMs.var(pm, nw, c)[:w] = Dict{Int, Any}([(id, Wr[id][c,c]) for id in bus_ids])
    end
end


""
function variable_mc_branch_series_current_prod_hermitian(pm::AbstractUBFModels; n_cond::Int=3, nw::Int=pm.cnw, bounded = true)
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

            smax = branch["rate_a"].values
            cmaxfr = smax./vmin_fr + abs.(y_fr)*vmax_fr
            cmaxto = smax./vmin_to + abs.(y_to)*vmax_to

            cmax[key] = max.(cmaxfr, cmaxto)
        end
        # create matrix variables
        (Lr,Li) = variable_mx_hermitian_sqrt_bounds(pm.model, branch_ids, n_cond, cmax; name="CC", prefix="$nw")
    else
        (Lr,Li) = variable_mx_hermitian(pm.model, branch_ids, n_cond; name="CC", prefix="$nw", lb_diag_zero=true)
    end

    # save reference
    _PMs.var(pm, nw)[:CCr] = Lr
    _PMs.var(pm, nw)[:CCi] = Li
    for c in 1:n_cond
        _PMs.var(pm, nw, c)[:cm] = Dict([(id, Lr[id][c,c]) for id in branch_ids])
    end
end


""
function variable_mc_branch_flow(pm::AbstractUBFModels; n_cond::Int=3, nw::Int=pm.cnw, bounded = true)
    @assert n_cond<=5

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
function constraint_mc_theta_ref(pm::AbstractUBFModels, i::Int; nw::Int=pm.cnw)
    constraint_mc_theta_ref(pm, nw, i)
end


"Defines branch flow model power flow equations"
function constraint_mc_flow_losses(pm::AbstractUBFModels, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to)
    p_to = _PMs.var(pm, n, :P)[t_idx]
    q_to = _PMs.var(pm, n, :Q)[t_idx]

    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    w_to_re = _PMs.var(pm, n, :Wr)[t_bus]
    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]

    w_to_im = _PMs.var(pm, n, :Wi)[t_bus]
    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]

    ccm_re =  _PMs.var(pm, n, :CCr)[i]
    ccm_im =  _PMs.var(pm, n, :CCi)[i]

    JuMP.@constraint(pm.model, p_fr + p_to .==  w_fr_re*(g_sh_fr)' + w_fr_im*(b_sh_fr)' + r*ccm_re - x*ccm_im +  w_to_re*(g_sh_to)'  + w_to_im*(b_sh_to)')
    JuMP.@constraint(pm.model, q_fr + q_to .==  w_fr_im*(g_sh_fr)' - w_fr_re*(b_sh_fr)' + x*ccm_re + r*ccm_im +  w_to_im*(g_sh_to)'  - w_to_re*(b_sh_to)')
end


""
function constraint_mc_theta_ref(pm::AbstractUBFModels, n::Int, i)
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
function constraint_mc_model_voltage_magnitude_difference(pm::AbstractUBFModels, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
    w_fr_re = _PMs.var(pm, n, :Wr)[f_bus]
    w_fr_im = _PMs.var(pm, n, :Wi)[f_bus]

    w_to_re = _PMs.var(pm, n, :Wr)[t_bus]
    w_to_im = _PMs.var(pm, n, :Wi)[t_bus]

    p_fr = _PMs.var(pm, n, :P)[f_idx]
    q_fr = _PMs.var(pm, n, :Q)[f_idx]

    p_s_fr = p_fr - (w_fr_re*(g_sh_fr)' + w_fr_im*(b_sh_fr)')
    q_s_fr = q_fr - (w_fr_im*(g_sh_fr)' - w_fr_re*(b_sh_fr)')

    ccm_re =  _PMs.var(pm, n, :CCr)[i]
    ccm_im =  _PMs.var(pm, n, :CCi)[i]

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


function variable_mc_generation_mx(pm::AbstractUBFModels; nw=pm.cnw)
    variable_mc_generation_current_mx(pm; nw=nw)
    variable_mc_generation_power_mx(pm; nw=nw)
end


function variable_mc_generation_power_mx(pm::AbstractUBFModels; nw=pm.cnw)
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


function variable_mc_generation_current_mx(pm::AbstractUBFModels; nw=pm.cnw)
    gen_ids = collect(_PMs.ids(pm, nw, :gen))
    ncnds = length(_PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(gen_ids), Array{Real,2}}()
    for (id, gen) in _PMs.ref(pm, nw, :gen)
        bus = _PMs.ref(pm, nw, :bus, gen["gen_bus"])
        cmax = _gen_curr_max(gen, bus)
        bound[id] = cmax*cmax'
    end
    # create matrix variables
    (CCgr,CCgi) = variable_mx_hermitian(pm.model, gen_ids, ncnds, bound; name="CCg", prefix="$nw")
    # save references
    _PMs.var(pm, nw)[:CCgr] = CCgr
    _PMs.var(pm, nw)[:CCgi] = CCgi
end


function variable_mc_load_mx(pm::AbstractUBFModels; nw=pm.cnw)
    load_wye_ids = [id for (id, load) in _PMs.ref(pm, nw, :load) if load["conn"]=="wye"]
    load_del_ids = [id for (id, load) in _PMs.ref(pm, nw, :load) if load["conn"]=="delta"]
    load_cone_ids = [id for (id, load) in _PMs.ref(pm, nw, :load) if _load_needs_cone(load)]
    # first create the wye loads; will create keys :Pd, :Qd
    variable_mc_load_power_wye_mx(pm, load_wye_ids)
    # now, create delta loads; will create :Xdr, :Xdi
    variable_mc_load_power_delta(pm, load_del_ids)
    # define :Pd, :Qd for delta loads as lin. transformation of :Xdr and :Wdim
    Td = [1 -1 0; 0 1 -1; -1 0 1]
    for id in load_del_ids
        _PMs.var(pm, nw, :Pd)[id] = _PMs.var(pm, nw, :Xdr, id)*Td
        _PMs.var(pm, nw, :Qd)[id] = _PMs.var(pm, nw, :Xdi, id)*Td
    end
    # both loads need a current variable; bounds adjusted for connection type
    variable_mc_load_current_mx(pm, [load_wye_ids..., load_del_ids...])
    variable_mc_load_power_vector(pm, load_cone_ids)
end


function variable_mc_load(pm::AbstractUBFModels; nw=pm.cnw)
    load_wye_ids = [id for (id, load) in _PMs.ref(pm, nw, :load) if load["conn"]=="wye"]
    load_del_ids = [id for (id, load) in _PMs.ref(pm, nw, :load) if load["conn"]=="delta"]
    load_cone_ids = [id for (id, load) in _PMs.ref(pm, nw, :load) if _load_needs_cone(load)]
    # create dictionary for wye loads
    for c in _PMs.conductor_ids(pm)
        _PMs.var(pm, nw, c)[:pd] = Dict()
        _PMs.var(pm, nw, c)[:qd] = Dict()
    end
    # now, create delta loads; will create :Xdr, :Xdi
    variable_mc_load_power_delta(pm, load_del_ids)
    variable_mc_load_current_mx(pm, load_del_ids)
    # only delta loads need a current matrix variable
    variable_mc_load_power_vector(pm, load_cone_ids)
end


function variable_mc_load_power_vector(pm::AbstractUBFModels, load_ids::Array{Int,1}; nw=pm.cnw)
    # calculate bounds for all loads
    pmin = Dict()
    pmax = Dict()
    qmin = Dict()
    qmax = Dict()
    for id in load_ids
        load = _PMs.ref(pm, nw, :load, id)
        bus = _PMs.ref(pm, nw, :bus, load["load_bus"])
        pmin[id], pmax[id], qmin[id], qmax[id] = _load_pq_bounds(load, bus)
    end

    # create variables
    ncnds = length(_PMs.conductor_ids(pm, nw))
    for c in 1:ncnds
        pl = JuMP.@variable(pm.model, [id in load_ids], base_name="$(nw)_$(c)_pl",
            lower_bound=pmin[id][c], upper_bound=pmax[id][c]
        )
        ql = JuMP.@variable(pm.model, [id in load_ids], base_name="$(nw)_$(c)_ql",
            lower_bound=qmin[id][c], upper_bound=qmax[id][c]
        )
        # save as dict and not JuMP Array, so constants can be added later
        _PMs.var(pm, nw, c)[:pl] = Dict{Int, Any}([(id, pl[id]) for id in load_ids])
        _PMs.var(pm, nw, c)[:ql] = Dict{Int, Any}([(id, ql[id]) for id in load_ids])
    end
end


function variable_mc_load_power_wye_mx(pm::AbstractUBFModels, load_ids::Array{Int,1}; nw=pm.cnw)
    ncnds = length(_PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(load_ids), Array{Real,2}}()
    for id in load_ids
        load = _PMs.ref(pm, nw, :load, id)
        @assert(load["conn"]=="wye")
        bus = _PMs.ref(pm, nw, :bus, load["load_bus"])
        cmax = _load_curr_max(load, bus)
        bound[id] = bus["vmax"].values*cmax'
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
function variable_mc_load_power_delta(pm::AbstractUBFModels, load_ids::Array{Int,1}; nw=pm.cnw, eps=0.1)
    ncnds = length(_PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(load_ids), Array{Real,2}}()
    for id in load_ids
        load = _PMs.ref(pm, nw, :load, id)
        bus_id = load["load_bus"]
        bus = _PMs.ref(pm, nw, :bus, bus_id)
        cmax = _load_curr_max(load, bus)
        bound[id] = bus["vmax"].values*cmax'
    end
    # create matrix variables
    (Xdre,Xdim) = variable_mx_complex(pm.model, load_ids, ncnds, ncnds, bound; name="Xd", prefix="$nw")
    # save references
    _PMs.var(pm, nw)[:Xdr] = Xdre
    _PMs.var(pm, nw)[:Xdi] = Xdim
end


function variable_mc_load_current_mx(pm::AbstractUBFModels, load_ids::Array{Int,1}; nw=pm.cnw)
    ncnds = length(_PMs.conductor_ids(pm, nw))
    # calculate bounds
    cmin = Dict{eltype(load_ids), Array{Real,1}}()
    cmax = Dict{eltype(load_ids), Array{Real,1}}()
    for (id, load) in _PMs.ref(pm, nw, :load)
        bus_id = load["load_bus"]
        bus = _PMs.ref(pm, nw, :bus, bus_id)
        cmin[id], cmax[id] = _load_curr_mag_bounds(load, bus)
    end
    # create matrix variables
    (CCdr, CCdi) = variable_mx_hermitian_sqrt_bounds(pm.model, load_ids, ncnds, cmax, cmin; name="CCd", prefix="$nw")
    # save references
    _PMs.var(pm, nw)[:CCdr] = CCdr
    _PMs.var(pm, nw)[:CCdi] = CCdi
end


function constraint_mc_generation_mx_SWL(pm::AbstractUBFModels, gen_id::Int; nw::Int=pm.cnw)
    Pg = _PMs.var(pm, nw, :Pg, gen_id)
    Qg = _PMs.var(pm, nw, :Qg, gen_id)
    bus_id = _PMs.ref(pm, nw, :gen, gen_id)["gen_bus"]
    Wr = _PMs.var(pm, nw, :Wr, bus_id)
    Wi = _PMs.var(pm, nw, :Wi, bus_id)
    CCgr = _PMs.var(pm, nw, :CCgr, gen_id)
    CCgi = _PMs.var(pm, nw, :CCgi, gen_id)
    constraint_SWL_psd(pm.model, Pg, Qg, Wr, Wi, CCgr, CCgi)
end


function constraint_mc_load_vector(pm::AbstractUBFModels, load_id::Int; nw=pm.cnw)
    load = _PMs.ref(pm, nw, :load, load_id)
    pd = load["pd"].values
    qd = load["qd"].values
    bus = _PMs.ref(pm, nw, :bus, load["load_bus"])
    ncnds = length(pd)

    if load["model"]=="constant_power"
        for c in 1:ncnds
            _PMs.var(pm, nw, c, :pl)[load_id] = pd[c]
            _PMs.var(pm, nw, c, :ql)[load_id] = qd[c]
            println("qd: $(qd[c])")
        end
    else
        # obtain a reference for w
        Wyr = _PMs.var(pm, nw, :Wr, bus["index"])
        Td = [1 -1 0; 0 1 -1; -1 0 1]
        if load["conn"]=="wye"
            Wr = Wyr
        elseif load["conn"]=="delta"
            Wr = Td*Wyr*Td'
        end
        # create a, α, b and β
        a, α, b, β = _load_expmodel_params(load, bus)
        if load["model"]=="constant_impedance"
            for c in 1:ncnds
                _PMs.var(pm, nw, c, :pl)[load_id] = a[c]*Wr[c,c]
                _PMs.var(pm, nw, c, :ql)[load_id] = b[c]*Wr[c,c]
            end
        else
            vmin, vmax = _load_vbounds(load, bus)
            wmin = vmin.^2
            wmax = vmax.^2
            pmin, pmax, qmin, qmax = _load_pq_bounds(load, bus)
            for c in 1:ncnds
                pl = _PMs.var(pm, nw, c, :pl)[load_id]
                ql = _PMs.var(pm, nw, c, :ql)[load_id]
                constraint_pqw(pm.model, Wr[c,c], pl, a[c], α[c], wmin[c], wmax[c], pmin[c], pmax[c])
                constraint_pqw(pm.model, Wr[c,c], ql, b[c], β[c], wmin[c], wmax[c], qmin[c], qmax[c])
            end
        end
    end
end

function constraint_pqw(model::JuMP.Model, w, p, a::Real, α::Real, wmin::Real, wmax::Real, pmin::Real, pmax::Real)
    if a==0
        JuMP.@constraint(model, p==0)
    else
        # AFFINE BOUNDARY
        if a>0
            l = (1/a)*(pmax-pmin)/(wmax-wmin)*(w-wmin) + pmin/a
        else
            # swap pmin and pmax if a<0, because pmin/a > pmax/a
            l = (1/a)*(pmin-pmax)/(wmax-wmin)*(w-wmin) + pmax/a
        end
        # affine overestimator
        if α>2
            JuMP.@constraint(model, p/a <= l)
        # affine underestimator
        elseif 0<α<2
            JuMP.@constraint(model, p/a >= l)
        end
        # CONE INCLUSIONS
        # constant current case
        # simplifies to a RotatedSecondOrderCone
        @assert(α>=0, "α has to greater than or equal to zero.")
        if α==0
            JuMP.@constraint(model, p==a)
        elseif α==1
            #       p/a <= w^(1/2)
            # <=>   (p/a)^2 <= w
            # <=>   2*(w/2)*1 >= ||p/a||^2_2
            # <=>   (w/2, 1, p/a) ∈ RotatedSecondOrderCone(3)
            JuMP.@constraint(model, [w/2, 1, p/a] in JuMP.RotatedSecondOrderCone())
        elseif 0<α<2
            #       p/a <= w^(α/2)
            # <=>   w^(α/2) >= p/a
            # <=>   (w, 1, p/a) ∈ PowerCone(3)
            JuMP.@constraint(model, [w, 1, p/a] in MathOptInterface.PowerCone(α/2))
        elseif α==2
            JuMP.@constraint(model, p==a*w)
        else # α>2
            #       p/a >= w^(α/2)
            # <=>   (p/a)^(2/α) >= w
            # <=>   (p/a, 1, w) ∈ PowerCone(3)
            JuMP.@constraint(model, [p/a, 1, w] in MathOptInterface.PowerCone(2/α))
        end
    end
end


function constraint_mc_load(pm::AbstractUBFModels, load_id::Int; nw::Int=pm.cnw)
    # shared variables and parameters
    load = _PMs.ref(pm, nw, :load, load_id)
    pd = load["pd"].values
    qd = load["qd"].values
    bus_id = load["load_bus"]
    ncnds = length(pd)

    # take care of voltage-dependency load models
    constraint_mc_load_vector(pm, load_id, nw=nw)
    pl = [_PMs.var(pm, nw, c, :pl, load_id) for c in 1:ncnds]
    ql = [_PMs.var(pm, nw, c, :ql, load_id) for c in 1:ncnds]

    # take care of connections
    if load["conn"]=="wye"
        for c in 1:length(pd)
            _PMs.var(pm, nw, c, :pd)[load_id] = pl[c]
            _PMs.var(pm, nw, c, :qd)[load_id] = ql[c]
        end
    elseif load["conn"]=="delta"
        Wr = _PMs.var(pm, nw, :Wr, bus_id)
        Wi = _PMs.var(pm, nw, :Wi, bus_id)
        CCdr = _PMs.var(pm, nw, :CCdr, load_id)
        CCdi = _PMs.var(pm, nw, :CCdi, load_id)
        Xdr = _PMs.var(pm, nw, :Xdr, load_id)
        Xdi = _PMs.var(pm, nw, :Xdi, load_id)
        Td = [1 -1 0; 0 1 -1; -1 0 1]
        JuMP.@constraint(pm.model, LinearAlgebra.diag(Td*Xdr) .== pl)
        JuMP.@constraint(pm.model, LinearAlgebra.diag(Td*Xdi) .== ql)
        constraint_SWL_psd(pm.model, Xdr, Xdi, Wr, Wi, CCdr, CCdi)
        Pd = Xdr*Td
        Qd = Xdi*Td
        for c in 1:length(pd)
            _PMs.var(pm, nw, c, :pd)[load_id] = Pd[c,c]
            _PMs.var(pm, nw, c, :qd)[load_id] = Qd[c,c]
        end
    end
end



function constraint_mc_load_mx(pm::AbstractUBFModels, load_id::Int; nw::Int=pm.cnw)
    # shared variables and parameters
    load = _PMs.ref(pm, nw, :load, load_id)
    bus_id = load["load_bus"]
    ncnds = length(load["pd"])

    # take care of voltage-dependency load models
    constraint_mc_load_vector(pm, load_id, nw=nw)
    pl = [_PMs.var(pm, nw, c, :pl, load_id) for c in 1:ncnds]
    ql = [_PMs.var(pm, nw, c, :ql, load_id) for c in 1:ncnds]


    Wr = _PMs.var(pm, nw, :Wr, bus_id)
    Wi = _PMs.var(pm, nw, :Wi, bus_id)
    CCdr = _PMs.var(pm, nw, :CCdr, load_id)
    CCdi = _PMs.var(pm, nw, :CCdi, load_id)
    if load["conn"]=="wye"
        # set the diagonal values
        Pd = _PMs.var(pm, nw, :Pd, load_id)
        Qd = _PMs.var(pm, nw, :Qd, load_id)
        for c in 1:ncnds
            Pd[c,c] = pl[c]
            Qd[c,c] = ql[c]
            _PMs.var(pm, nw, c, :pd)[load_id] = pl[c]
            _PMs.var(pm, nw, c, :qd)[load_id] = ql[c]
        end
        # link S, W and L
        #constraint_SWL_psd(pm.model, Pd, Qd, Wr, Wi, CCdr, CCdi)
    elseif load["conn"]=="delta"
        Xdre = _PMs.var(pm, nw, :Xdr, load_id)
        Xdim = _PMs.var(pm, nw, :Xdi, load_id)
        Td = [1 -1 0; 0 1 -1; -1 0 1]
        JuMP.@constraint(pm.model, LinearAlgebra.diag(Td*Xdre) .== pl)
        JuMP.@constraint(pm.model, LinearAlgebra.diag(Td*Xdim) .== ql)
        #constraint_SWL_psd(pm.model, Xdre, Xdim, Wr, Wi, CCdr, CCdi)
        Pd = Xdre*Td
        Qd = Xdim*Td
        _PMs.var(pm, nw, :Pd)[load_id] = Pd
        _PMs.var(pm, nw, :Qd)[load_id] = Qd
    end
end


function constraint_mc_load_mx_SWL(pm::AbstractUBFModels, load_id::Int; nw::Int=pm.cnw)
    # shared variables and parameters
    load = _PMs.ref(pm, nw, :load, load_id)
    bus_id = load["load_bus"]
    ncnds = length(load["pd"])

    Wr = _PMs.var(pm, nw, :Wr, bus_id)
    Wi = _PMs.var(pm, nw, :Wi, bus_id)
    CCdr = _PMs.var(pm, nw, :CCdr, load_id)
    CCdi = _PMs.var(pm, nw, :CCdi, load_id)
    if load["conn"]=="wye"
        # set the diagonal values
        Pd = _PMs.var(pm, nw, :Pd, load_id)
        Qd = _PMs.var(pm, nw, :Qd, load_id)
        # link S, W and L
        constraint_SWL_psd(pm.model, Pd, Qd, Wr, Wi, CCdr, CCdi)
    elseif load["conn"]=="delta"
        Xdre = _PMs.var(pm, nw, :Xdr, load_id)
        Xdim = _PMs.var(pm, nw, :Xdi, load_id)
        constraint_SWL_psd(pm.model, Xdre, Xdim, Wr, Wi, CCdr, CCdi)
    end
end


function constraint_mc_load_mx_SWL_only_delta(pm::AbstractUBFModels, load_id::Int; nw::Int=pm.cnw)
    # shared variables and parameters
    load = _PMs.ref(pm, nw, :load, load_id)
    bus_id = load["load_bus"]
    ncnds = length(load["pd"])

    Wr = _PMs.var(pm, nw, :Wr, bus_id)
    Wi = _PMs.var(pm, nw, :Wi, bus_id)
    CCdr = _PMs.var(pm, nw, :CCdr, load_id)
    CCdi = _PMs.var(pm, nw, :CCdi, load_id)
    if load["conn"]=="delta"
        Xdre = _PMs.var(pm, nw, :Xdr, load_id)
        Xdim = _PMs.var(pm, nw, :Xdi, load_id)
        constraint_SWL_psd(pm.model, Xdre, Xdim, Wr, Wi, CCdr, CCdi)
    end
end


function constraint_mc_voltage_psd(pm::AbstractUBFModels; nw=pm.cnw)
    buses_covered = [i for (l,i,j) in _PMs.ref(pm, nw, :arcs)]
    buses_psd = [i for i in _PMs.ids(pm, nw, :bus) if !(i in buses_covered)]
    for bus_id in buses_psd
        Wr = _PMs.var(pm, nw, :Wr, bus_id)
        Wi = _PMs.var(pm, nw, :Wi, bus_id)
        constraint_M_psd(Wr, Wi)
    end
end


function constraint_SWL_psd(model::JuMP.Model, P, Q, Wr, Wi, L_re, L_im)
    M_re = [Wr P; P' L_re]
    M_im = [Wi Q; -Q' L_im]
    constraint_M_psd(model, M_re, M_im)
end


function constraint_M_psd(model::JuMP.Model, M_re, M_im)
    JuMP.@constraint(model, [M_re -M_im; M_im M_re] in JuMP.PSDCone())
end


""
function constraint_mc_power_balance_shunt(pm::AbstractUBFModels, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    if !haskey(_PMs.con(pm, nw, cnd), :kcl_p)
        _PMs.con(pm, nw, cnd)[:kcl_p] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(_PMs.con(pm, nw, cnd), :kcl_q)
        _PMs.con(pm, nw, cnd)[:kcl_q] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = _PMs.ref(pm, nw, :bus, i)
    bus_arcs = _PMs.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = _PMs.ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = _PMs.ref(pm, nw, :bus_gens, i)
    bus_loads = _PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PMs.ref(pm, nw, :bus_shunts, i)

    bus_gs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "gs", cnd) for k in bus_shunts)
    bus_bs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "bs", cnd) for k in bus_shunts)

    constraint_power_balance_shunt(pm, nw, cnd, i, bus_arcs, bus_arcs_dc, bus_gens, bus_loads, bus_gs, bus_bs)
end


function constraint_mc_power_balance_shunt(pm::AbstractUBFModels, n::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_loads, bus_gs, bus_bs)
    w    = _PMs.var(pm, n, c, :w, i)
    pg   = _PMs.var(pm, n, c, :pg)
    qg   = _PMs.var(pm, n, c, :qg)
    pd   = _PMs.var(pm, n, c, :pd)
    qd   = _PMs.var(pm, n, c, :qd)
    p    = _PMs.var(pm, n, c, :p)
    q    = _PMs.var(pm, n, c, :q)
    p_dc = _PMs.var(pm, n, c, :p_dc)
    q_dc = _PMs.var(pm, n, c, :q_dc)

    JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd[d] for d in bus_loads) - sum(gs for gs in values(bus_gs))*w)
    JuMP.@constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd[d] for d in bus_loads) + sum(bs for bs in values(bus_bs))*w)
end


""
function constraint_mc_power_balance_mx_shunt(pm::AbstractUBFModels, i::Int; nw::Int=pm.cnw)
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

    constraint_mc_power_balance_mx_shunt(pm, nw, i, bus_arcs, bus_arcs_dc, bus_gens, bus_loads, bus_Gs, bus_Bs)
end


"""
Shunt handling in matrix form:
I = Y.U
S = U.I' = U.(Y.U)' = U.U'.Y' = W.Y'
  = (Wr+j.Wi)(G+jB)' = (Wr+j.Wi)(G'-j.B') = (Wr.G'+Wi.B')+j(-Wr.B'+Wi.G')
P =  Wr.G'+Wi.B'
Q = -Wr.B'+Wi.G'
"""
function constraint_mc_power_balance_mx_shunt(pm::AbstractUBFModels, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_loads, bus_Gs, bus_Bs)
    Wr = _PMs.var(pm, n, :Wr, i)
    Wi = _PMs.var(pm, n, :Wi, i)
    P = _PMs.var(pm, n, :P)
    Q = _PMs.var(pm, n, :Q)
    Pg = _PMs.var(pm, n, :Pg)
    Qg = _PMs.var(pm, n, :Qg)
    Pd = _PMs.var(pm, n, :Pd)
    Qd = _PMs.var(pm, n, :Qd)
    # ignore dc for now
    #TODO add DC in matrix version?
    ncnds = size(Wr)[1]
    G = (length(bus_Gs)>0) ? sum(values(bus_Gs)) : zeros(ncnds, ncnds)
    B = (length(bus_Bs)>0) ? sum(values(bus_Bs)) : zeros(ncnds, ncnds)

    # changed the ordering
    # LHS: all variables with generator sign convention
    # RHS: all variables with load sign convention
    _PMs.con(pm, n, :kcl_P)[i] = JuMP.@constraint(pm.model, sum(Pg[g] for g in bus_gens) .== sum(P[a] for a in bus_arcs) + sum(Pd[d] for d in bus_loads) + ( Wr*G'+Wi*B'))
    _PMs.con(pm, n, :kcl_Q)[i] = JuMP.@constraint(pm.model, sum(Qg[g] for g in bus_gens) .== sum(Q[a] for a in bus_arcs) + sum(Qd[d] for d in bus_loads) + (-Wr*B'+Wi*G'))
end
