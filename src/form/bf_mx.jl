import LinearAlgebra: diag


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
    n_diag_el = n_cond
    n_lower_triangle_el = Int((n_cond^2 - n_cond)/2)
    for c in 1:n_diag_el
        PowerModels.variable_voltage_magnitude_sqr(pm, nw=nw, cnd=c, bounded=bounded)
    end

    wmaxdict = Dict{Int64, Any}()
    for i in _PMs.ids(pm, nw, :bus)
        wmax = _PMs.ref(pm, nw, :bus, i, "vmax").values*_PMs.ref(pm, nw, :bus, i, "vmax").values'
        wmaxltri = _mat2ltrivec!(wmax)
        wmaxdict[i] = wmaxltri
    end

    for c in 1:n_lower_triangle_el
        if bounded
            _PMs.var(pm, nw, c)[:wr] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_$(c)_wr",
            lower_bound = -wmaxdict[i][c],
            upper_bound =  wmaxdict[i][c],
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i), "w_start", c, 1.001)
            )
            _PMs.var(pm, nw, c)[:wi] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_$(c)_wi",
            lower_bound = -wmaxdict[i][c],
            upper_bound =  wmaxdict[i][c],
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i), "w_start", c, 1.001)
            )
        else
            _PMs.var(pm, nw, c)[:wr] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_$(c)_wr",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i), "w_start", c, 1.001)
            )
            _PMs.var(pm, nw, c)[:wi] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_$(c)_wi",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i), "w_start", c, 1.001)
            )
        end
    end

    #Store dictionary with matrix variables by bus
    w_re_dict = Dict{Int64, Any}()
    w_im_dict = Dict{Int64, Any}()
    for i in _PMs.ids(pm, nw, :bus)
        w =  [_PMs.var(pm, nw, h, :w,  i) for h in 1:n_diag_el]
        wr = [_PMs.var(pm, nw, h, :wr, i) for h in 1:n_lower_triangle_el]
        wi = [_PMs.var(pm, nw, h, :wi, i) for h in 1:n_lower_triangle_el]

        (w_re, w_im) = _make_hermitian_matrix_variable(w, wr, wi)
        w_re_dict[i] = w_re
        w_im_dict[i] = w_im
    end
    _PMs.var(pm, nw)[:W_re] = w_re_dict
    _PMs.var(pm, nw)[:W_im] = w_im_dict
end


""
function variable_mc_branch_series_current_prod_hermitian(pm::AbstractUBFModels; n_cond::Int=3, nw::Int=pm.cnw, bounded = true)
    branches = _PMs.ref(pm, nw, :branch)
    buses = _PMs.ref(pm, nw, :bus)

    n_diag_el = n_cond
    n_lower_triangle_el = Int((n_cond^2 - n_cond)/2)

    cmax = Dict([(key, zeros(n_cond)) for key in keys(branches)])

    for (key, branch) in branches
        bus_fr = buses[branch["f_bus"]]
        bus_to = buses[branch["t_bus"]]

        vmin_fr = bus_fr["vmin"].values
        vmin_to = bus_to["vmin"].values

        vmax_fr = bus_fr["vmax"].values
        vmax_to = bus_to["vmax"].values

        y_fr_mag = abs.(branch["g_fr"].values + im* branch["b_fr"].values)
        y_to_mag = abs.(branch["g_to"].values + im* branch["b_to"].values)

        smax = branch["rate_a"].values
        cmaxfr = smax./vmin_fr + y_fr_mag*vmax_fr
        cmaxto = smax./vmin_to + y_to_mag*vmax_to

        cmax[key] = max.(cmaxfr, cmaxto)
    end


    for c in 1:n_diag_el
        if bounded
            _PMs.var(pm, nw, c)[:ccm] = JuMP.@variable(pm.model,
            [l in _PMs.ids(pm, nw, :branch)], base_name="$(nw)_$(c)_cm",
            lower_bound = 0,
            upper_bound = (cmax[l][c])^2,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "i_start", c) #TODO shouldn't this be squared?
            )
        else
            _PMs.var(pm, nw, c)[:ccm] = JuMP.@variable(pm.model,
            [l in _PMs.ids(pm, nw, :branch)], base_name="$(nw)_$(c)_cm",
            lower_bound = 0,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "i_start", c)
            )
        end
        # PowerModels.variable_current_magnitude_sqr(pm, nw=nw, cnd=c)
    end

    for c in 1:n_lower_triangle_el
        if bounded
            _PMs.var(pm, nw, c)[:ccmr] = JuMP.@variable(pm.model,
            [l in _PMs.ids(pm, nw, :branch)], base_name="$(nw)_$(c)_ccmr",
            lower_bound = -(cmax[l][c])^2,
            upper_bound = (cmax[l][c])^2,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "i_start", c) #TODO shouldn't this be squared?
            )
            _PMs.var(pm, nw, c)[:ccmi] = JuMP.@variable(pm.model,
            [l in _PMs.ids(pm, nw, :branch)], base_name="$(nw)_$(c)_ccmi",
            lower_bound = -(cmax[l][c])^2,
            upper_bound = (cmax[l][c])^2,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "i_start", c)
            )
        else
            _PMs.var(pm, nw, c)[:ccmr] = JuMP.@variable(pm.model,
            [l in _PMs.ids(pm, nw, :branch)], base_name="$(nw)_$(c)_ccmr",
            # lower_bound = 0,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "i_start", c)
            )
            _PMs.var(pm, nw, c)[:ccmi] = JuMP.@variable(pm.model,
            [l in _PMs.ids(pm, nw, :branch)], base_name="$(nw)_$(c)_ccmi",
            # lower_bound = 0,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "i_start", c)
            )
        end
    end

    #Store dictionary with matrix variables by branch
    ccm_re_dict = Dict{Int64, Any}()
    ccm_im_dict = Dict{Int64, Any}()
    for i in _PMs.ids(pm, nw, :branch)
        ccm =  [_PMs.var(pm, nw, h, :ccm,  i) for h in 1:n_diag_el]
        ccmr = [_PMs.var(pm, nw, h, :ccmr, i) for h in 1:n_lower_triangle_el]
        ccmi = [_PMs.var(pm, nw, h, :ccmi, i) for h in 1:n_lower_triangle_el]

        (ccm_re, ccm_im) = _make_hermitian_matrix_variable(ccm, ccmr, ccmi)
        ccm_re_dict[i] = ccm_re
        ccm_im_dict[i] = ccm_im
    end
    _PMs.var(pm, nw)[:CC_re] = ccm_re_dict
    _PMs.var(pm, nw)[:CC_im] = ccm_im_dict
end


""
function variable_mc_branch_flow(pm::AbstractUBFModels; n_cond::Int=3, nw::Int=pm.cnw, bounded = true)
    n_diag_el = n_cond
    n_lower_triangle_el = Int((n_cond^2 - n_cond)/2)
    @assert n_cond<=5

    for i in 1:n_diag_el
        _PMs.variable_active_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
        _PMs.variable_reactive_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
    end

    for i in 1:n_lower_triangle_el
        variable_mc_lower_triangle_active_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
        variable_mc_lower_triangle_reactive_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
        variable_mc_upper_triangle_active_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
        variable_mc_upper_triangle_reactive_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
    end
    #Store dictionary with matrix variables by arc
    p_mat_dict = Dict{Tuple{Int64,Int64,Int64}, Any}()
    q_mat_dict = Dict{Tuple{Int64,Int64,Int64}, Any}()

    for i in _PMs.ref(pm, nw, :arcs)
        p_d =  [_PMs.var(pm, nw, c, :p,    i) for c in 1:n_diag_el]
        p_ut = [_PMs.var(pm, nw, c, :p_ut, i) for c in 1:n_lower_triangle_el]
        p_lt = [_PMs.var(pm, nw, c, :p_lt, i) for c in 1:n_lower_triangle_el]

        q_d =  [_PMs.var(pm, nw, c, :q,    i) for c in 1:n_diag_el]
        q_ut = [_PMs.var(pm, nw, c, :q_ut, i) for c in 1:n_lower_triangle_el]
        q_lt = [_PMs.var(pm, nw, c, :q_lt, i) for c in 1:n_lower_triangle_el]

        p_mat = _make_full_matrix_variable(p_d, p_lt, p_ut)
        q_mat = _make_full_matrix_variable(q_d, q_lt, q_ut)

        p_mat_dict[i] = p_mat
        q_mat_dict[i] = q_mat
    end
    _PMs.var(pm, nw)[:P_mx] = p_mat_dict
    _PMs.var(pm, nw)[:Q_mx] = q_mat_dict
end


"variable: `p_lt[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_lower_triangle_active_branch_flow(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded = true)
    smaxdict = Dict{Tuple{Int64, Int64, Int64}, Any}()

    for (l,i,j) in _PMs.ref(pm, nw, :arcs)
        cmax = _PMs.ref(pm, nw, :branch, l, "rate_a").values./_PMs.ref(pm, nw, :bus, i, "vmin").values
        smax = _PMs.ref(pm, nw, :bus, i, "vmax").values.*cmax'
        smaxltri = _mat2ltrivec!(smax)
        smaxdict[(l,i,j)] = smaxltri
    end

    if bounded
        _PMs.var(pm, nw, cnd)[:p_lt] = JuMP.@variable(pm.model,
        [(l,i,j) in _PMs.ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_p_lt",
        lower_bound = -smaxdict[(l,i,j)][cnd],
        upper_bound =  smaxdict[(l,i,j)][cnd],
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "p_start", cnd)
        )
    else
        _PMs.var(pm, nw, cnd)[:p_lt] = JuMP.@variable(pm.model,
        [(l,i,j) in _PMs.ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_p_lt",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "p_start", cnd)
        )
    end
end


"variable: `q_lt[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_lower_triangle_reactive_branch_flow(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded = true)
    smaxdict = Dict{Tuple{Int64, Int64, Int64}, Any}()

    for (l,i,j) in _PMs.ref(pm, nw, :arcs)
        cmax = _PMs.ref(pm, nw, :branch, l, "rate_a").values./_PMs.ref(pm, nw, :bus, i, "vmin").values
        smax = _PMs.ref(pm, nw, :bus, i, "vmax").values.*cmax'
        smaxltri = _mat2ltrivec!(smax)
        smaxdict[(l,i,j)] = smaxltri
    end

    if bounded
        _PMs.var(pm, nw, cnd)[:q_lt] = JuMP.@variable(pm.model,
        [(l,i,j) in _PMs.ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_q_lt",
        lower_bound = -smaxdict[(l,i,j)][cnd],
        upper_bound =  smaxdict[(l,i,j)][cnd],
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "q_start", cnd)
        )
    else
        _PMs.var(pm, nw, cnd)[:q_lt] = JuMP.@variable(pm.model,
        [(l,i,j) in _PMs.ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_q_lt",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "q_start", cnd)
        )
    end
end


"variable: `p_ut[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_upper_triangle_active_branch_flow(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded = true)
    smaxdict = Dict{Tuple{Int64, Int64, Int64}, Any}()

    for (l,i,j) in _PMs.ref(pm, nw, :arcs)
        cmax = _PMs.ref(pm, nw, :branch, l, "rate_a").values./_PMs.ref(pm, nw, :bus, i, "vmin").values
        smax = _PMs.ref(pm, nw, :bus, i, "vmax").values.*cmax'
        smaxutri = _mat2utrivec!(smax)
        smaxdict[(l,i,j)] = smaxutri
    end

    if bounded
        _PMs.var(pm, nw, cnd)[:p_ut] = JuMP.@variable(pm.model,
        [(l,i,j) in _PMs.ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_p_ut",
        lower_bound = -smaxdict[(l,i,j)][cnd],
        upper_bound =  smaxdict[(l,i,j)][cnd],
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "p_start", cnd)
        )
    else
        _PMs.var(pm, nw, cnd)[:p_ut] = JuMP.@variable(pm.model,
        [(l,i,j) in _PMs.ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_p_ut",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "p_start", cnd)
        )
    end
end


"variable: `q_ut[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_upper_triangle_reactive_branch_flow(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded = true)
    smaxdict = Dict{Tuple{Int64, Int64, Int64}, Any}()

    for (l,i,j) in _PMs.ref(pm, nw, :arcs)
        cmax = _PMs.ref(pm, nw, :branch, l, "rate_a").values./_PMs.ref(pm, nw, :bus, i, "vmin").values
        smax = _PMs.ref(pm, nw, :bus, i, "vmax").values.*cmax'
        smaxutri = _mat2utrivec!(smax)
        smaxdict[(l,i,j)] = smaxutri
    end
    if bounded
        _PMs.var(pm, nw, cnd)[:q_ut] = JuMP.@variable(pm.model,
        [(l,i,j) in _PMs.ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_q_ut",
        lower_bound = -smaxdict[(l,i,j)][cnd],
        upper_bound =  smaxdict[(l,i,j)][cnd],
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "q_start", cnd)
        )
    else
        _PMs.var(pm, nw, cnd)[:q_ut] = JuMP.@variable(pm.model,
        [(l,i,j) in _PMs.ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_q_ut",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :branch, l), "q_start", cnd)
        )
    end
end


""
function constraint_mc_theta_ref(pm::AbstractUBFModels, i::Int; nw::Int=pm.cnw)
    constraint_mc_theta_ref(pm, nw, i)
end


"Defines branch flow model power flow equations"
function constraint_mc_flow_losses(pm::AbstractUBFModels, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to)
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
function constraint_mc_theta_ref(pm::AbstractUBFModels, n::Int, i)
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
function constraint_mc_model_voltage_magnitude_difference(pm::AbstractUBFModels, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
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
