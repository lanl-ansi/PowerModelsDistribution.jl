"This is duplicated at PMD level to correctly handle the indexing of the shunts."
function constraint_mc_voltage_angle_difference(pm::_PMs.AbstractBFModel, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx
    t_idx = (i, t_bus, f_bus)

    branch = _PMs.ref(pm, n, :branch, i)

    for c in _PMs.conductor_ids(pm; nw=n)
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
        w_fr = _PMs.var(pm, n, :w, f_bus)[c]
        p_fr = _PMs.var(pm, n, :p, f_idx)[c]
        q_fr = _PMs.var(pm, n, :q, f_idx)[c]

        tzr = r*tr + x*ti
        tzi = r*ti - x*tr

        JuMP.@constraint(pm.model,
            tan(angmin[c])*((tr + tzr*g_fr + tzi*b_fr)*(w_fr/tm^2) - tzr*p_fr + tzi*q_fr)
                     <= ((ti + tzi*g_fr - tzr*b_fr)*(w_fr/tm^2) - tzi*p_fr - tzr*q_fr)
            )
        JuMP.@constraint(pm.model,
            tan(angmax[c])*((tr + tzr*g_fr + tzi*b_fr)*(w_fr/tm^2) - tzr*p_fr + tzi*q_fr)
                     >= ((ti + tzi*g_fr - tzr*b_fr)*(w_fr/tm^2) - tzi*p_fr - tzr*q_fr)
            )
    end
end


"Create voltage variables for branch flow model"
function variable_mc_bus_voltage_on_off(pm::LPUBFDiagModel; kwargs...)
    variable_mc_voltage_magnitude_sqr_on_off(pm; kwargs...)
end


"nothing to do, this model is symmetric"
function constraint_mc_trans_yy(pm::LPUBFDiagModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
end
