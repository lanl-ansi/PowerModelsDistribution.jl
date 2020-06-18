"This is duplicated at PowerModelsDistribution level to correctly handle the indexing of the shunts."
function constraint_mc_voltage_angle_difference(pm::_PM.AbstractBFModel, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx
    t_idx = (i, t_bus, f_bus)

    branch = ref(pm, n, :branch, i)

    for c in conductor_ids(pm; nw=n)
        tm = branch["tap"][c]
        g_fr = branch["g_fr"][c,c]
        g_to = branch["g_to"][c,c]
        b_fr = branch["b_fr"][c,c]
        b_to = branch["b_to"][c,c]

        tr, ti = _PM.calc_branch_t(branch)
        tr, ti = tr[c], ti[c]

        r = branch["br_r"][c,c]
        x = branch["br_x"][c,c]

        # getting the variables
        w_fr = var(pm, n, :w, f_bus)[c]
        p_fr = var(pm, n, :p, f_idx)[c]
        q_fr = var(pm, n, :q, f_idx)[c]

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
    variable_mc_bus_voltage_magnitude_sqr_on_off(pm; kwargs...)
end


# TODO: Throw error if tm_fixed is false
"""
Links to and from power and voltages in a wye-wye transformer, assumes tm_fixed is true
w_fr_i=(pol_i*tm_scale*tm_i)^2w_to_i
"""
function constraint_mc_transformer_power_yy(pm::LPUBFDiagModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
    tm = [tm_fixed[c] ? tm_set[c] : var(pm, nw, :tap, trans_id)[c] for c in conductor_ids(pm)]
    nph = length(conductor_ids(pm))

    p_fr = [var(pm, nw, :pt, f_idx)[p] for p in f_cnd]
    p_to = [var(pm, nw, :pt, t_idx)[p] for p in t_cnd]
    q_fr = [var(pm, nw, :qt, f_idx)[p] for p in f_cnd]
    q_to = [var(pm, nw, :qt, t_idx)[p] for p in t_cnd]

    w_fr = var(pm, nw, :w)[f_bus]
    w_to = var(pm, nw, :w)[t_bus]

    for p in 1:nph
        JuMP.@constraint(pm.model, w_fr[p] == (pol*tm_scale*tm[p])^2*w_to[p])
    end

    JuMP.@constraint(pm.model, p_fr + p_to .== 0)
    JuMP.@constraint(pm.model, q_fr + q_to .== 0)
end


# TODO: Throw error if tm_fixed is false
raw"""
Links to and from power and voltages in a delta-wye transformer, assumes tm_fixed is true
3(w_fr_i+w_fr_j)=2(pol_i*tm_scale*tm_i)^2w_to_i \quad \forall (i,j) \in \{(1,2),(2,3),(3,1)\}

2P_fr_i=-(P_to_i+P_to_j)+(Q_to_j-Q_to_i)/\sqrt{3} \quad \forall (i,j) \in \{(1,3),(2,1),(3,2)\}
2Q_fr_i=(P_to_i-P_to_j)/\sqrt{3}-(Q_to_j+Q_to_i)  \quad \forall (i,j) \in \{(1,3),(2,1),(3,2)\}
"""
function constraint_mc_transformer_power_dy(pm::LPUBFDiagModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
    tm = [tm_fixed[c] ? tm_set[c] : var(pm, nw, :tap, trans_id)[c] for c in conductor_ids(pm)]
    nph = length(conductor_ids(pm))

    p_fr = [var(pm, nw, :pt, f_idx)[p] for p in f_cnd]
    p_to = [var(pm, nw, :pt, t_idx)[p] for p in t_cnd]
    q_fr = [var(pm, nw, :qt, f_idx)[p] for p in f_cnd]
    q_to = [var(pm, nw, :qt, t_idx)[p] for p in t_cnd]

    w_fr = var(pm, nw, :w)[f_bus]
    w_to = var(pm, nw, :w)[t_bus]

    for p in 1:nph
        # rotate by 1 to get 'previous' phase
        # e.g., for nph=3: 1->3, 2->1, 3->2
        q = (p-1+1)%nph+1
	    JuMP.@constraint(pm.model, 3.0*(w_fr[p] + w_fr[q]) == 2.0*(pol*tm_scale*tm[p])^2*w_to[p])
    end

    for p in 1:nph
        # rotate by nph-1 to get 'previous' phase
        # e.g., for nph=3: 1->3, 2->1, 3->2
        q = (p-1+nph-1)%nph+1
	    JuMP.@constraint(pm.model, 2*p_fr[p] == -(p_to[p]+p_to[q])+(q_to[q]-q_to[p])/sqrt(3.0))
        JuMP.@constraint(pm.model, 2*q_fr[p] == (p_to[p]-p_to[q])/sqrt(3.0)-(q_to[q]+q_to[p]))
    end
end
