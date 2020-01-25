


""
function variable_mc_voltage(pm::LPLinUBFModel; kwargs...)
    variable_mc_voltage_magnitude_sqr(pm; kwargs...)
end


"nothing to do, variables not used in linearised branch flow model"
function variable_mc_branch_current(pm::LPLinUBFModel; kwargs...)
end


"Defines voltage drop over a branch, linking from and to side voltage magnitude"
function constraint_mc_model_voltage_magnitude_difference(pm::LPLinUBFModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
    p_fr = _PMs.var(pm, n, :p, f_idx)
    q_fr = _PMs.var(pm, n, :q, f_idx)
    w_fr = _PMs.var(pm, n, :w, f_bus)
    w_to = _PMs.var(pm, n, :w, t_bus)

    for c in _PMs.conductor_ids(pm, n)
        np = length(_PMs.conductor_ids(pm))
        rot = _roll([_wrap_to_pi(2*pi/np*(1-d)) for d in _PMs.conductor_ids(pm)], c-1)

        #KVL over the line:
        JuMP.@constraint(pm.model, w_to[c] == w_fr[c] - 2*sum((r[c,d]*cos(rot[d])-x[c,d]*sin(rot[d]))*(p_fr[d] - g_sh_fr[c]*(w_fr[c]/tm[c]^2)) +
                                                   (r[c,d]*sin(rot[d])+x[c,d]*cos(rot[d]))*(q_fr[d] + b_sh_fr[c]*(w_fr[c]/tm[c]^2)) for d in _PMs.conductor_ids(pm)) )
    end
end


"do nothing"
function constraint_mc_model_current(pm::LPLinUBFModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
end


"Defines branch flow model power flow equations"
function constraint_mc_flow_losses(pm::LPLinUBFModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to, tm)
    p_fr = _PMs.var(pm, n, :p, f_idx)
    q_fr = _PMs.var(pm, n, :q, f_idx)
    p_to = _PMs.var(pm, n, :p, t_idx)
    q_to = _PMs.var(pm, n, :q, t_idx)
    w_fr = _PMs.var(pm, n, :w, f_bus)
    w_to = _PMs.var(pm, n, :w, t_bus)

    JuMP.@constraint(pm.model, p_fr + p_to .==  diag(g_sh_fr).*(w_fr./tm.^2) +  diag(g_sh_to).*w_to)
    JuMP.@constraint(pm.model, q_fr + q_to .== -diag(b_sh_fr).*(w_fr./tm.^2) + -diag(b_sh_to).*w_to)
end


"Create voltage variables for branch flow model"
function variable_mc_bus_voltage_on_off(pm::LPLinUBFModel; kwargs...)
    variable_mc_voltage_magnitude_sqr_on_off(pm; kwargs...)
end


"nothing to do, this model is symmetric"
function constraint_mc_trans(pm::LPLinUBFModel, i::Int; nw::Int=pm.cnw)
end

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
