""
function variable_tp_voltage(pm::_PMs.GenericPowerModel{T}; kwargs...) where T <: LPLinUBFForm
    for cnd in _PMs.conductor_ids(pm)
        _PMs.variable_voltage_magnitude_sqr(pm; cnd=cnd, kwargs...)
    end
end


"nothing to do, variables not used in linearised branch flow model"
function variable_tp_branch_current(pm::_PMs.GenericPowerModel{T}; kwargs...) where T <: LPLinUBFForm
end


"Defines voltage drop over a branch, linking from and to side voltage magnitude"
function constraint_tp_voltage_magnitude_difference(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm) where T <: LPLinUBFForm
    p_fr = [_PMs.var(pm, n, d, :p, f_idx) for d in _PMs.conductor_ids(pm)]
    q_fr = [_PMs.var(pm, n, d, :q, f_idx) for d in _PMs.conductor_ids(pm)]
    w_fr = _PMs.var(pm, n, c, :w, f_bus)
    w_to = _PMs.var(pm, n, c, :w, t_bus)

    np = length(_PMs.conductor_ids(pm))
    rot = _roll([_wrap_to_pi(2*pi/np*(1-d)) for d in _PMs.conductor_ids(pm)], c-1)

    #KVL over the line:
    JuMP.@constraint(pm.model, w_to == w_fr - 2*sum((r[c,d]*cos(rot[d])-x[c,d]*sin(rot[d]))*(p_fr[d] - g_sh_fr*(w_fr/tm^2)) +
                                               (r[c,d]*sin(rot[d])+x[c,d]*cos(rot[d]))*(q_fr[d] + b_sh_fr*(w_fr/tm^2)) for d in _PMs.conductor_ids(pm)) )
end


"do nothing"
function _PMs.constraint_model_current(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr, tm) where T <: LPLinUBFForm
end


"Defines branch flow model power flow loss equations"
function constraint_tp_flow_losses(pm::_PMs.GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to) where T <: LPLinUBFForm
    tm = [1,1,1] #TODO
    for c in _PMs.conductor_ids(pm)
        constraint_tp_flow_losses(pm, n, c, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr[c,c], g_sh_to[c,c], b_sh_fr[c,c], b_sh_to[c,c], tm[c])
    end
end


"Defines branch flow model power flow equations"
function constraint_tp_model_voltage_magnitude_difference(pm::_PMs.GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm) where T <: LPLinUBFForm
    tm = [1,1,1] #TODO
    for c in _PMs.conductor_ids(pm)
        constraint_tp_voltage_magnitude_difference(pm, n, c, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr[c,c], b_sh_fr[c,c], tm[c])
    end
end


"Defines branch flow model power flow equations"
function constraint_tp_flow_losses(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to, tm) where T <: LPLinUBFForm
    p_fr = _PMs.var(pm, n, c, :p, f_idx)
    q_fr = _PMs.var(pm, n, c, :q, f_idx)
    p_to = _PMs.var(pm, n, c, :p, t_idx)
    q_to = _PMs.var(pm, n, c, :q, t_idx)
    w_fr = _PMs.var(pm, n, c, :w, f_bus)
    w_to = _PMs.var(pm, n, c, :w, t_bus)

    JuMP.@constraint(pm.model, p_fr + p_to ==  g_sh_fr*(w_fr/tm^2) +  g_sh_to*w_to)
    JuMP.@constraint(pm.model, q_fr + q_to == -b_sh_fr*(w_fr/tm^2) + -b_sh_to*w_to)
end


"Create voltage variables for branch flow model"
function variable_tp_bus_voltage_on_off(pm::_PMs.GenericPowerModel{T}; kwargs...) where T <: LPLinUBFForm
    for cnd in _PMs.conductor_ids(pm)
        variable_tp_voltage_magnitude_sqr_on_off(pm; cnd=cnd, kwargs...)
    end
end


"nothing to do, no voltage variables"
function constraint_tp_trans_voltage(pm::_PMs.GenericPowerModel{T}, nw::Int, i::Int, f_bus::Int, t_bus::Int, tm::_PMs.MultiConductorVector, Tv_fr, Tv_im, Cv_to) where T <: LPLinUBFForm
end


"nothing to do, this model is symmetric"
function constraint_tp_trans_flow(pm::_PMs.GenericPowerModel{T}, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, tm::_PMs.MultiConductorVector, Ti_fr, Ti_im, Cv_to) where T <: LPLinUBFForm
end


"This is duplicated at PMD level to correctly handle the indexing of the shunts."
function constraint_voltage_angle_difference(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, f_idx, angmin, angmax) where T <: _PMs.AbstractBFForm
    i, f_bus, t_bus = f_idx
    t_idx = (i, t_bus, f_bus)

    branch = _PMs.ref(pm, n, :branch, i)
    tm = branch["tap"][c]
    g_fr = branch["g_fr"][c,c]
    g_to = branch["g_to"][c,c]
    b_fr = branch["b_fr"][c,c]
    b_to = branch["b_to"][c,c]

    tr, ti = _PMs.calc_branch_t(branch)
    tr, ti = tr[c], ti[c]

    r = branch["br_r"][c,c]
    x = branch["br_x"][c,c]

    # getting the variables
    w_fr = _PMs.var(pm, n, c, :w, f_bus)
    p_fr = _PMs.var(pm, n, c, :p, f_idx)
    q_fr = _PMs.var(pm, n, c, :q, f_idx)

    tzr = r*tr + x*ti
    tzi = r*ti - x*tr

    JuMP.@constraint(pm.model,
        tan(angmin)*((tr + tzr*g_fr + tzi*b_fr)*(w_fr/tm^2) - tzr*p_fr + tzi*q_fr)
                 <= ((ti + tzi*g_fr - tzr*b_fr)*(w_fr/tm^2) - tzi*p_fr - tzr*q_fr)
        )
    JuMP.@constraint(pm.model,
        tan(angmax)*((tr + tzr*g_fr + tzi*b_fr)*(w_fr/tm^2) - tzr*p_fr + tzi*q_fr)
                 >= ((ti + tzi*g_fr - tzr*b_fr)*(w_fr/tm^2) - tzi*p_fr - tzr*q_fr)
        )
end
