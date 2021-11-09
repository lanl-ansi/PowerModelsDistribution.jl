"This is duplicated at PowerModelsDistribution level to correctly handle the indexing of the shunts."
function constraint_mc_voltage_angle_difference(pm::AbstractUBFModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})
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


"""
Links to and from power and voltages in a wye-wye transformer, assumes tm_fixed is true

```math
w_fr_i=(pol_i*tm_scale*tm_i)^2w_to_i
```
"""
function constraint_mc_transformer_power_yy(pm::LPUBFDiagModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))]
    transformer = ref(pm, nw, :transformer, trans_id)

    p_fr = [var(pm, nw, :pt, f_idx)[p] for p in f_connections]
    p_to = [var(pm, nw, :pt, t_idx)[p] for p in t_connections]
    q_fr = [var(pm, nw, :qt, f_idx)[p] for p in f_connections]
    q_to = [var(pm, nw, :qt, t_idx)[p] for p in t_connections]

    w_fr = var(pm, nw, :w)[f_bus]
    w_to = var(pm, nw, :w)[t_bus]

    tmsqr = [tm_fixed[i] ? tm[i]^2 : JuMP.@variable(pm.model, base_name="$(nw)_tmsqr_$(trans_id)_$(f_connections[i])", start=JuMP.start_value(tm[i])^2, lower_bound=JuMP.has_lower_bound(tm[i]) ? JuMP.lower_bound(tm[i])^2 : 0.0, upper_bound=JuMP.has_upper_bound(tm[i]) ? JuMP.upper_bound(tm[i])^2 : 10.0^2) for i in 1:length(tm)]

    for (idx, (fc, tc)) in enumerate(zip(f_connections, t_connections))
        if tm_fixed[idx]
            JuMP.@constraint(pm.model, w_fr[fc] == (pol*tm_scale*tm[idx])^2*w_to[tc])
        else
            PolyhedralRelaxations.construct_univariate_relaxation!(pm.model, x->x^2, tm[idx], tmsqr[idx], [JuMP.has_lower_bound(tm[idx]) ? JuMP.lower_bound(tm[idx]) : 0.0, JuMP.has_upper_bound(tm[idx]) ? JuMP.upper_bound(tm[idx]) : 10.0], false)

            tmsqr_w_to = JuMP.@variable(pm.model, base_name="$(nw)_tmsqr_w_to_$(trans_id)_$(t_bus)_$(tc)")
            PolyhedralRelaxations.construct_bilinear_relaxation!(pm.model, tmsqr[idx], w_to[tc], tmsqr_w_to, [JuMP.lower_bound(tmsqr[idx]), JuMP.upper_bound(tmsqr[idx])], [JuMP.has_lower_bound(w_to[tc]) ? JuMP.lower_bound(w_to[tc]) : 0.0^2, JuMP.has_upper_bound(w_to[tc]) ? JuMP.upper_bound(w_to[tc]) : 10.0^2])

            JuMP.@constraint(pm.model, w_fr[fc] == (pol*tm_scale)^2*tmsqr_w_to)

            # with regcontrol
            if haskey(transformer,"controls")
                v_ref = transformer["controls"]["vreg"][idx]
                δ = transformer["controls"]["band"][idx]
                r = transformer["controls"]["r"][idx]
                x = transformer["controls"]["x"][idx]

                # linearized voltage squared: w_drop = (2⋅r⋅p+2⋅x⋅q)
                w_drop = JuMP.@expression(pm.model, 2*r*p_to[idx] + 2*x*q_to[idx])

                # (v_ref-δ)^2 ≤ w_fr-w_drop ≤ (v_ref+δ)^2
                # w_fr/1.1^2 ≤ w_to ≤ w_fr/0.9^2
                JuMP.@constraint(pm.model, w_fr[fc] - w_drop ≥ (v_ref - δ)^2)
                JuMP.@constraint(pm.model, w_fr[fc] - w_drop ≤ (v_ref + δ)^2)
                JuMP.@constraint(pm.model, w_fr[fc]/1.1^2 ≤ w_to[tc])
                JuMP.@constraint(pm.model, w_fr[fc]/0.9^2 ≥ w_to[tc])
            end
        end
    end

    JuMP.@constraint(pm.model, p_fr + p_to .== 0)
    JuMP.@constraint(pm.model, q_fr + q_to .== 0)
end


raw"""
Links to and from power and voltages in a delta-wye transformer, assumes tm_fixed is true

```math
\begin{align}
3(w_fr_i+w_fr_j)=2(pol_i*tm_scale*tm_i)^2w_to_i & \quad \forall (i,j) \in \{(1,2),(2,3),(3,1)\} \\
2P_fr_i=-(P_to_i+P_to_j)+(Q_to_j-Q_to_i)/\sqrt{3} & \quad \forall (i,j) \in \{(1,3),(2,1),(3,2)\} \\
2Q_fr_i=(P_to_i-P_to_j)/\sqrt{3}-(Q_to_j+Q_to_i)  & \quad \forall (i,j) \in \{(1,3),(2,1),(3,2)\}
\end{align}
````
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


"Neglects the active and reactive loss terms associated with the squared current magnitude."
function constraint_storage_losses(pm::AbstractUBFAModel, n::Int, i, bus, r, x, p_loss, q_loss; conductors=[1])
    ps = var(pm, n, :ps, i)
    qs = var(pm, n, :qs, i)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    qsc = var(pm, n, :qsc, i)

    JuMP.@constraint(pm.model,
        sum(ps[c] for c in conductors) + (sd - sc)
        ==
        p_loss
    )

    JuMP.@constraint(pm.model,
        sum(qs[c] for c in conductors)
        ==
        qsc + q_loss
    )
end
