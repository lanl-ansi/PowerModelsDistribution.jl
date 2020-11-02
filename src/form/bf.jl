"This is duplicated at PowerModelsDistribution level to correctly handle the indexing of the shunts."
function constraint_mc_voltage_angle_difference(pm::_PM.AbstractBFModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})
    i, f_bus, t_bus = f_idx
    t_idx = (i, t_bus, f_bus)

    branch = ref(pm, nw, :branch, i)

    for (idx, (fc, tc)) in enumerate(zip(branch["f_connections"], branch["t_connections"]))
        g_fr = branch["g_fr"][idx,idx]
        g_to = branch["g_to"][idx,idx]
        b_fr = branch["b_fr"][idx,idx]
        b_to = branch["b_to"][idx,idx]

        r = branch["br_r"][idx,idx]
        x = branch["br_x"][idx,idx]

        # getting the variables
        w_fr = var(pm, nw, :w, f_bus)[fc]
        p_fr = var(pm, nw, :p, f_idx)[fc]
        q_fr = var(pm, nw, :q, f_idx)[fc]

        JuMP.@constraint(pm.model,
            tan(angmin[idx])*((1 + r*g_fr - x*b_fr)*(w_fr) - r*p_fr - x*q_fr)
                     <= ((-x*g_fr - r*b_fr)*(w_fr) + x*p_fr - r*q_fr)
            )
        JuMP.@constraint(pm.model,
            tan(angmax[idx])*((1 + r*g_fr - x*b_fr)*(w_fr) - r*p_fr - x*q_fr)
                     >= ((-x*g_fr - r*b_fr)*(w_fr) + x*p_fr - r*q_fr)
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
function constraint_mc_transformer_power_yy(pm::LPUBFDiagModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[fc] for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))]

    nph = length(tm_set)

    p_fr = [var(pm, nw, :pt, f_idx)[p] for p in f_connections]
    p_to = [var(pm, nw, :pt, t_idx)[p] for p in t_connections]
    q_fr = [var(pm, nw, :qt, f_idx)[p] for p in f_connections]
    q_to = [var(pm, nw, :qt, t_idx)[p] for p in t_connections]

    w_fr = var(pm, nw, :w)[f_bus]
    w_to = var(pm, nw, :w)[t_bus]

    for (idx, (fc, tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, w_fr[fc] == (pol*tm_scale*tm[idx])^2*w_to[tc])
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
function constraint_mc_transformer_power_dy(pm::LPUBFDiagModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[fc] for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))]
    nph = length(tm_set)

    p_fr = [var(pm, nw, :pt, f_idx)[c] for c in f_connections]
    p_to = [var(pm, nw, :pt, t_idx)[c] for c in t_connections]
    q_fr = [var(pm, nw, :qt, f_idx)[c] for c in f_connections]
    q_to = [var(pm, nw, :qt, t_idx)[c] for c in t_connections]

    w_fr = var(pm, nw, :w)[f_bus]
    w_to = var(pm, nw, :w)[t_bus]

    for (idx,(fc, tc)) in enumerate(zip(f_connections,t_connections))
        # rotate by 1 to get 'previous' phase
        # e.g., for nph=3: 1->3, 2->1, 3->2
        jdx = (idx-1+1)%nph+1
        fd = f_connections[jdx]
	    JuMP.@constraint(pm.model, 3.0*(w_fr[fc] + w_fr[fd]) == 2.0*(pol*tm_scale*tm[idx])^2*w_to[tc])
    end

    for (idx,(fc, tc)) in enumerate(zip(f_connections,t_connections))
        # rotate by nph-1 to get 'previous' phase
        # e.g., for nph=3: 1->3, 2->1, 3->2
        jdx = (idx-1+nph-1)%nph+1
        fd = f_connections[jdx]
        td = t_connections[jdx]
	    JuMP.@constraint(pm.model, 2*p_fr[fc] == -(p_to[tc]+p_to[td])+(q_to[td]-q_to[tc])/sqrt(3.0))
        JuMP.@constraint(pm.model, 2*q_fr[fc] == (p_to[tc]-p_to[td])/sqrt(3.0)-(q_to[td]+q_to[tc]))
    end
end
