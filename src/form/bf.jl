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
function variable_mc_bus_voltage_on_off(pm::LPUBFDiagModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_bus_voltage_magnitude_sqr_on_off(pm; nw=nw, bounded=bounded, report=report)
end


"""
Links to and from power and voltages in a wye-wye transformer, assumes tm_fixed is true

```math
w_fr_i=(pol_i*tm_scale*tm_i)^2w_to_i
```
"""
function constraint_mc_transformer_power_yy(pm::AbstractUBFModels, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))]
    transformer = ref(pm, nw, :transformer, trans_id)

    p_fr = [var(pm, nw, :pt, f_idx)[p] for p in f_connections]
    p_to = [var(pm, nw, :pt, t_idx)[p] for p in t_connections]
    q_fr = [var(pm, nw, :qt, f_idx)[p] for p in f_connections]
    q_to = [var(pm, nw, :qt, t_idx)[p] for p in t_connections]

    w_fr = var(pm, nw, :w)[f_bus]
    w_to = var(pm, nw, :w)[t_bus]

    tmsqr = [tm_fixed[i] ? tm[i]^2 : JuMP.@variable(pm.model, base_name="$(nw)_tmsqr_$(trans_id)_$(f_connections[i])", start=JuMP.start_value(tm[i])^2, lower_bound=JuMP.has_lower_bound(tm[i]) ? JuMP.lower_bound(tm[i])^2 : 0.9^2, upper_bound=JuMP.has_upper_bound(tm[i]) ? JuMP.upper_bound(tm[i])^2 : 1.1^2) for i in 1:length(tm)]

    for (idx, (fc, tc)) in enumerate(zip(f_connections, t_connections))
        if tm_fixed[idx]
            JuMP.@constraint(pm.model, w_fr[fc] == (pol*tm_scale*tm[idx])^2*w_to[tc])
        else
            PolyhedralRelaxations.construct_univariate_relaxation!(pm.model, x->x^2, tm[idx], tmsqr[idx], [JuMP.has_lower_bound(tm[idx]) ? JuMP.lower_bound(tm[idx]) : 0.9, JuMP.has_upper_bound(tm[idx]) ? JuMP.upper_bound(tm[idx]) : 1.1], false)

            tmsqr_w_to = JuMP.@variable(pm.model, base_name="$(nw)_tmsqr_w_to_$(trans_id)_$(t_bus)_$(tc)")
            PolyhedralRelaxations.construct_bilinear_relaxation!(pm.model, tmsqr[idx], w_to[tc], tmsqr_w_to, [JuMP.lower_bound(tmsqr[idx]), JuMP.upper_bound(tmsqr[idx])], [JuMP.has_lower_bound(w_to[tc]) ? JuMP.lower_bound(w_to[tc]) : 0.9^2, JuMP.has_upper_bound(w_to[tc]) ? JuMP.upper_bound(w_to[tc]) : 1.1^2])

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


"This function adds all constraints required to model a two-winding, delta-wye connected transformer."
function constraint_mc_transformer_power_dy(pm::AbstractUBFModels, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    nph = length(tm_set)
    @assert length(f_connections) == length(t_connections) && nph == 3 "only phases == 3 dy transformers are currently supported"

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[fc] for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))]
    @assert all(isa.(tm, Real)) "This formulation does not support variable taps"

    vr_fr = [JuMP.@variable(pm.model, base_name="$(nw)_vr_$(f_bus)_$(c)", start=0.0, lower_bound=-sqrt(JuMP.upper_bound(var(pm, nw, :w, f_bus)[c])), upper_bound=sqrt(JuMP.upper_bound(var(pm, nw, :w, f_bus)[c]))) for c in f_connections]
    vi_fr = [JuMP.@variable(pm.model, base_name="$(nw)_vi_$(f_bus)_$(c)", start=0.0, lower_bound=-sqrt(JuMP.upper_bound(var(pm, nw, :w, f_bus)[c])), upper_bound=sqrt(JuMP.upper_bound(var(pm, nw, :w, f_bus)[c]))) for c in f_connections]
    vr_to = [JuMP.@variable(pm.model, base_name="$(nw)_vr_$(t_bus)_$(c)", start=0.0, lower_bound=-sqrt(JuMP.upper_bound(var(pm, nw, :w, t_bus)[c])), upper_bound=sqrt(JuMP.upper_bound(var(pm, nw, :w, t_bus)[c]))) for c in t_connections]
    vi_to = [JuMP.@variable(pm.model, base_name="$(nw)_vi_$(t_bus)_$(c)", start=0.0, lower_bound=-sqrt(JuMP.upper_bound(var(pm, nw, :w, t_bus)[c])), upper_bound=sqrt(JuMP.upper_bound(var(pm, nw, :w, t_bus)[c]))) for c in t_connections]

    w_fr = [var(pm, nw, :w, f_bus)[c] for c in f_connections]
    w_to = [var(pm, nw, :w, t_bus)[c] for c in t_connections]

    secant_vr_fr = [JuMP.@variable(pm.model, base_name="$(nw)_secant_vr_$(f_bus)_$(c)", start=0.0) for c in f_connections]
    secant_vi_fr = [JuMP.@variable(pm.model, base_name="$(nw)_secant_vi_$(f_bus)_$(c)", start=0.0) for c in f_connections]
    secant_vr_to = [JuMP.@variable(pm.model, base_name="$(nw)_secant_vr_$(t_bus)_$(c)", start=0.0) for c in t_connections]
    secant_vi_to = [JuMP.@variable(pm.model, base_name="$(nw)_secant_vi_$(t_bus)_$(c)", start=0.0) for c in t_connections]

    JuMP.@constraint(pm.model, w_fr .>= vr_fr.^2 + vi_fr.^2)
    _IM.relaxation_sqr.(pm.model, vr_fr, secant_vr_fr)
    _IM.relaxation_sqr.(pm.model, vi_fr, secant_vi_fr)
    JuMP.@constraint(pm.model, w_fr .<= secant_vr_fr + secant_vi_fr)

    JuMP.@constraint(pm.model, w_to .>= vr_to.^2 + vi_to.^2)
    _IM.relaxation_sqr.(pm.model, vr_to, secant_vr_to)
    _IM.relaxation_sqr.(pm.model, vi_to, secant_vi_to)
    JuMP.@constraint(pm.model, w_to .<= secant_vr_to + secant_vi_to)

    p_fr = var(pm, nw, :pt, f_idx)
    p_to = var(pm, nw, :pt, t_idx)
    q_fr = var(pm, nw, :qt, f_idx)
    q_to = var(pm, nw, :qt, t_idx)

    c = (pol*tm_scale)*tm

    lift_w_to_ir_to = [JuMP.@variable(pm.model, base_name="$(nw)_lift_w_ir_$(t_bus)_$(c)", start=0.0) for c in t_connections]
    lift_w_to_ii_to = [JuMP.@variable(pm.model, base_name="$(nw)_lift_w_ii_$(t_bus)_$(c)", start=0.0) for c in t_connections]
    lift_p_to_vr_to = [JuMP.@variable(pm.model, base_name="$(nw)_lift_p_vr_$(t_bus)_$(c)", start=0.0) for c in t_connections]
    lift_q_to_vi_to = [JuMP.@variable(pm.model, base_name="$(nw)_lift_p_vi_$(t_bus)_$(c)", start=0.0) for c in t_connections]

    c_fr_ub, c_to_ub = _calc_transformer_current_max_frto(ref(pm, nw, :transformer, trans_id), ref(pm, nw, :bus, f_bus), ref(pm, nw, :bus, t_bus))

    cr_to = [JuMP.@variable(pm.model, base_name="$(nw)_cr_$(t_bus)_$(c)", start=0.0, lower_bound=-c_to_ub[idx], upper_bound=c_to_ub[idx]) for (idx,c) in enumerate(t_connections)]
    ci_to = [JuMP.@variable(pm.model, base_name="$(nw)_ci_$(t_bus)_$(c)", start=0.0, lower_bound=-c_to_ub[idx], upper_bound=c_to_ub[idx]) for (idx,c) in enumerate(t_connections)]

    JuMP.@constraint(pm.model, lift_w_to_ir_to .== (lift_p_to_vr_to + lift_q_to_vi_to) ./ c)
    JuMP.@constraint(pm.model, lift_w_to_ii_to .== (lift_p_to_vr_to - lift_q_to_vi_to) ./ c)

    _IM.relaxation_product.(pm.model, w_to, cr_to, lift_w_to_ir_to)
    _IM.relaxation_product.(pm.model, w_to, ci_to, lift_w_to_ii_to)

    _IM.relaxation_product.(pm.model, p_to, vr_to, lift_p_to_vr_to)
    _IM.relaxation_product.(pm.model, q_to, vi_to, lift_q_to_vi_to)

    M = _get_delta_transformation_matrix(nph)

    # introduce auxialiary variable vd = Md*v_fr
    vrd = M*vr_fr
    vid = M*vi_fr

    JuMP.@constraint(pm.model, vrd .== (pol*tm_scale)*tm.*vr_to)
    JuMP.@constraint(pm.model, vid .== (pol*tm_scale)*tm.*vi_to)

    lift_vr_fr_ir_to = Array{Any,2}(undef, nph, nph)
    lift_vi_fr_ii_to = Array{Any,2}(undef, nph, nph)
    lift_vr_fr_ii_to = Array{Any,2}(undef, nph, nph)
    lift_vi_fr_ir_to = Array{Any,2}(undef, nph, nph)
    for (idx, fc) in enumerate(f_connections)
        for (jdx, tc) in enumerate(t_connections)
            if jdx >= idx
                lift_vr_fr_ir_to[idx,jdx] = JuMP.@variable(pm.model, base_name="$(nw)_lift_vr_fr_ir_to_$(fc)_$(tc)", start=0.0)
                lift_vi_fr_ii_to[idx,jdx] = JuMP.@variable(pm.model, base_name="$(nw)_lift_vi_fr_ii_to_$(fc)_$(tc)", start=0.0)
                lift_vi_fr_ii_to[idx,jdx] = JuMP.@variable(pm.model, base_name="$(nw)_lift_vi_fr_ii_to_$(fc)_$(tc)", start=0.0)
                lift_vr_fr_ii_to[idx,jdx] = JuMP.@variable(pm.model, base_name="$(nw)_lift_vr_fr_ii_to_$(fc)_$(tc)", start=0.0)
                lift_vi_fr_ir_to[idx,jdx] = JuMP.@variable(pm.model, base_name="$(nw)_lift_vi_fr_ir_to_$(fc)_$(tc)", start=0.0)

                if idx != jdx
                    lift_vr_fr_ir_to[jdx,idx] = lift_vr_fr_ir_to[idx,jdx]
                    lift_vi_fr_ii_to[jdx,idx] = lift_vi_fr_ii_to[idx,jdx]
                    lift_vr_fr_ii_to[jdx,idx] = lift_vr_fr_ii_to[idx,jdx]
                    lift_vi_fr_ir_to[jdx,idx] = lift_vi_fr_ir_to[idx,jdx]
                end
            end
        end
    end

    for (idx,(fc,tc)) in enumerate(zip(f_connections, t_connections))
        jdx = (idx-1+nph-1)%nph+1

        JuMP.@constraint(pm.model, p_fr[idx] == lift_vr_fr_ir_to[idx,jdx] - lift_vr_fr_ir_to[idx,idx] + lift_vi_fr_ii_to[idx,jdx] - lift_vi_fr_ii_to[idx,idx])
        JuMP.@constraint(pm.model, q_fr[idx] == lift_vi_fr_ir_to[idx,jdx] - lift_vi_fr_ir_to[idx,idx] - lift_vr_fr_ii_to[idx,jdx] - lift_vr_fr_ii_to[idx,idx])

        _IM.relaxation_product(pm.model, vr_fr[idx], cr_to[jdx], lift_vr_fr_ir_to[idx,jdx])
        _IM.relaxation_product(pm.model, vr_fr[idx], cr_to[idx], lift_vr_fr_ir_to[idx,idx])
        _IM.relaxation_product(pm.model, vi_fr[idx], ci_to[jdx], lift_vi_fr_ii_to[idx,jdx])
        _IM.relaxation_product(pm.model, vi_fr[idx], ci_to[idx], lift_vi_fr_ii_to[idx,idx])
        _IM.relaxation_product(pm.model, vi_fr[idx], cr_to[jdx], lift_vi_fr_ir_to[idx,jdx])
        _IM.relaxation_product(pm.model, vi_fr[idx], cr_to[idx], lift_vi_fr_ir_to[idx,idx])
        _IM.relaxation_product(pm.model, vr_fr[idx], ci_to[jdx], lift_vr_fr_ii_to[idx,jdx])
        _IM.relaxation_product(pm.model, vr_fr[idx], ci_to[idx], lift_vr_fr_ii_to[idx,idx])
    end
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
function constraint_mc_storage_losses(pm::AbstractUBFAModel, n::Int, i::Int, bus::Int, connections::Vector{Int}, r::Real, x::Real, p_loss::Real, q_loss::Real)
    ps = var(pm, n, :ps, i)
    qs = var(pm, n, :qs, i)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    qsc = var(pm, n, :qsc, i)

    JuMP.@constraint(pm.model,
        sum(ps[c] for c in connections) + (sd - sc)
        ==
        p_loss
    )

    JuMP.@constraint(pm.model,
        sum(qs[c] for c in connections)
        ==
        qsc + q_loss
    )
end
