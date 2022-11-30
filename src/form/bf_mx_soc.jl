"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SOCUBFModels, nw::Int, i::Int, f_bus::Int, f_idx::Tuple{Int,Int,Int}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real})
    branch = ref(pm, nw, :branch, f_idx[1])
    f_connections = branch["f_connections"]
    t_connections = branch["t_connections"]
    bus = ref(pm, nw, :bus, f_bus)
    terminals = bus["terminals"]

    p_fr = var(pm, nw, :P)[f_idx]
    q_fr = var(pm, nw, :Q)[f_idx]

    w_fr_re = var(pm, nw, :Wr, f_bus)[[findfirst(isequal(fc), terminals) for fc in f_connections],[findfirst(isequal(tc), terminals) for tc in t_connections]]
    w_fr_im = var(pm, nw, :Wi, f_bus)[[findfirst(isequal(fc), terminals) for fc in f_connections],[findfirst(isequal(tc), terminals) for tc in t_connections]]

    ccm_re =  var(pm, nw, :CCr)[i]
    ccm_im =  var(pm, nw, :CCi)[i]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

    mat_real = [
    w_fr_re     p_s_fr  ;
    p_s_fr'    ccm_re  ;
    ]

    mat_imag = [
    w_fr_im     q_s_fr  ;
    -q_s_fr'    ccm_im  ;
    ]

    relaxation_psd_to_soc(pm.model, mat_real, mat_imag, complex=true)

    # code below useful for debugging: valid inequality equired to make the SOC-NLP formulation more accurate
    # (l,i,j) = f_idx
    # t_idx = (l,j,i)
    # p_to = var(pm, n, :P)[t_idx]
    # total losses are positive when g_fr, g_to and r are positive
    # not guaranteed for individual phases though when matrix obtained through Kron's reduction
    # JuMP.@constraint(pm.model, tr(p_fr) + tr(p_to) >= 0)
end


"Add explicit PSD-ness of W for nodes where it is not implied"
function constraint_mc_voltage_psd(pm::SOCUBFModels, nw::Int, i::Int)
    Wr = var(pm, nw, :Wr)[i]
    Wi = var(pm, nw, :Wi)[i]

    relaxation_psd_to_soc(pm.model, Wr, Wi)
end


"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SOCConicUBFModel, nw::Int, i::Int, f_bus::Int, f_idx::Tuple{Int,Int,Int}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real})
    branch = ref(pm, nw, :branch, f_idx[1])
    f_connections = branch["f_connections"]
    t_connections = branch["t_connections"]
    bus = ref(pm, nw, :bus, f_bus)
    terminals = bus["terminals"]

    p_fr = var(pm, nw, :P)[f_idx]
    q_fr = var(pm, nw, :Q)[f_idx]

    w_fr_re = var(pm, nw, :Wr, f_bus)[[findfirst(isequal(fc), terminals) for fc in f_connections],[findfirst(isequal(tc), terminals) for tc in t_connections]]
    w_fr_im = var(pm, nw, :Wi, f_bus)[[findfirst(isequal(fc), terminals) for fc in f_connections],[findfirst(isequal(tc), terminals) for tc in t_connections]]

    ccm_re =  var(pm, nw, :CCr)[i]
    ccm_im =  var(pm, nw, :CCi)[i]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

    mat_real = [
    w_fr_re     p_s_fr  ;
    p_s_fr'    ccm_re  ;
    ]

    mat_imag = [
    w_fr_im     q_s_fr  ;
    -q_s_fr'    ccm_im  ;
    ]

    relaxation_psd_to_soc_conic(pm.model, mat_real, mat_imag, complex=true)
end


"Add explicit PSD-ness of W for nodes where it is not implied"
function constraint_mc_voltage_psd(pm::SOCConicUBFModel, nw::Int, i::Int)
    Wr = var(pm, nw, :Wr)[i]
    Wi = var(pm, nw, :Wi)[i]

    relaxation_psd_to_soc_conic(pm.model, Wr, Wi)
end


@doc raw"""
    constraint_mc_transformer_power_yy(pm::SOCUBFModels, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)

Constraints to model a two-winding, wye-wye connected transformer.

```math
\begin{align}
    & {W}_{fr} = {T}_{m}{T}_{m}^{H} {W}_{to}      \\
    & {s}_{fr} + {s}_{to} = 0 \\
    & \sum_{j=1}^3 S_{fr}^{ij} = 0, ~\sum_{j=1}^3 S_{to}^{ij} = 0, ~\forall i \in \{1,2,3\}
\end{align}
```
"""
function constraint_mc_transformer_power_yy(pm::SOCUBFModels, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    transformer = ref(pm, nw, :transformer, trans_id)
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))]

    Wr_fr = var(pm, nw, :Wr, f_bus)
    Wr_to = var(pm, nw, :Wr, t_bus)
    Wi_fr = var(pm, nw, :Wi, f_bus)
    Wi_to = var(pm, nw, :Wi, t_bus)

    P_fr = var(pm, nw, :Pt, f_idx)
    P_to = var(pm, nw, :Pt, t_idx)
    Q_fr = var(pm, nw, :Qt, f_idx)
    Q_to = var(pm, nw, :Qt, t_idx)
    p_fr = [var(pm, nw, :pt, f_idx)[c] for c in f_connections]
    p_to = [var(pm, nw, :pt, t_idx)[c] for c in t_connections]
    q_fr = [var(pm, nw, :qt, f_idx)[c] for c in f_connections]
    q_to = [var(pm, nw, :qt, t_idx)[c] for c in t_connections]

    for (f_idx,fc) in enumerate(f_connections)
        for (t_idx,tc) in enumerate(t_connections)
            if tm_fixed[t_idx]
                JuMP.@constraint(pm.model, Wr_fr[fc,tc] == (pol*tm_scale)^2*tm[f_idx]*tm[t_idx]*Wr_to[fc,tc])
                JuMP.@constraint(pm.model, Wi_fr[fc,tc] == (pol*tm_scale)^2*tm[f_idx]*tm[t_idx]*Wi_to[fc,tc])
            else
                PolyhedralRelaxations.construct_univariate_relaxation!(pm.model, x->x^2, tm[idx], tmsqr[idx], [JuMP.has_lower_bound(tm[idx]) ? JuMP.lower_bound(tm[idx]) : 0.9, JuMP.has_upper_bound(tm[idx]) ? JuMP.upper_bound(tm[idx]) : 1.1], false)
                tmsqr_Wr = JuMP.@variable(pm.model, base_name="$(nw)_tmsqr_Wr_to_$(trans_id)_$(f_bus)_$(fc)_$(t_bus)_$(tc)")
                tmsqr_Wi = JuMP.@variable(pm.model, base_name="$(nw)_tmsqr_Wi_to_$(trans_id)_$(f_bus)_$(fc)_$(t_bus)_$(tc)")
                PolyhedralRelaxations.construct_bilinear_relaxation!(pm.model, tmsqr[idx], Wr_to[fc,tc], tmsqr_Wr, [JuMP.lower_bound(tmsqr[idx]), JuMP.upper_bound(tmsqr[idx])], [JuMP.has_lower_bound(Wr_to[fc,tc]) ? JuMP.lower_bound(Wr_to[fc,tc]) : -(1.1^2), JuMP.has_upper_bound(Wr_to[fc,tc]) ? JuMP.upper_bound(Wr_to[fc,tc]) : 1.1^2])
                PolyhedralRelaxations.construct_bilinear_relaxation!(pm.model, tmsqr[idx], Wi_to[fc,tc], tmsqr_Wi, [JuMP.lower_bound(tmsqr[idx]), JuMP.upper_bound(tmsqr[idx])], [JuMP.has_lower_bound(Wi_to[fc,tc]) ? JuMP.lower_bound(Wi_to[fc,tc]) : -(1.1^2), JuMP.has_upper_bound(Wi_to[fc,tc]) ? JuMP.upper_bound(Wi_to[fc,tc]) : 1.1^2])

                JuMP.@constraint(pm.model, Wr_fr[fc,tc] == (pol*tm_scale)^2*tmsqr_Wr_to)
                JuMP.@constraint(pm.model, Wi_fr[fc,tc] == (pol*tm_scale)^2*tmsqr_Wi_to)

                # with regcontrol
                if haskey(transformer,"controls")
                    # voltage squared ignoring losses: w_drop = (2⋅r⋅p+2⋅x⋅q)
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
        if haskey(transformer,"controls")  # regcontrol settings
            w_fr = var(pm, nw, :w)[f_bus]
            w_to = var(pm, nw, :w)[t_bus]
            v_ref = transformer["controls"]["vreg"][f_idx]
            δ = transformer["controls"]["band"][f_idx]
            r = transformer["controls"]["r"][f_idx]
            x = transformer["controls"]["x"][f_idx]
        end

        # conservation of current
        JuMP.@constraint(pm.model, sum(P_fr[f_idx,:]) == 0)
        JuMP.@constraint(pm.model, sum(Q_fr[f_idx,:]) == 0)
        JuMP.@constraint(pm.model, sum(P_to[f_idx,:]) == 0)
        JuMP.@constraint(pm.model, sum(Q_to[f_idx,:]) == 0)
    end

    JuMP.@constraint(pm.model, p_fr + p_to .== 0)
    JuMP.@constraint(pm.model, q_fr + q_to .== 0)
end

@doc raw"""
    constraint_mc_transformer_power_dy(pm::SOCUBFModels, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)

Constraints to model a two-winding, delta-wye connected transformer.

```math
\begin{align}
    &{W}_{fr}^{ij}-{W}_{fr}^{ik}-{W}_{fr}^{lj}+{W}_{fr}^{lk} = t_m^2{W}_{to}^{ij} ~\forall i,j \in \{1,2,3\}~ \text{and}~ k,l \in \{2,3,1\}   \\
    &{S}_{fr} = X_tT_t \\
    &{S}_{fr}^\Delta = T_tX_t \\
    & {s}_{fr}^\Delta + {s}_{to} = 0\\
    & {M}_{\Delta} =
    \begin{bmatrix}
    {W}_{fr} & {X}_{t} \\
     {X}_{t}^{\text{H}} &  {L}_{\Delta}
    \end{bmatrix} \succeq 0
\end{align}
```
"""
function constraint_mc_transformer_power_dy(pm::SOCUBFModels, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    transformer = ref(pm, nw, :transformer, trans_id)
    bus_id = transformer["f_bus"]
    bus = ref(pm, nw, :bus, bus_id)

    nph = length(tm_set)
    @assert length(f_connections) == length(t_connections) && nph == 3 "only phases == 3 dy transformers are currently supported"
    next = Dict(c=>f_connections[idx%nph+1] for (idx,c) in enumerate(f_connections))

    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for (idx,(fc,tc)) in enumerate(zip(f_connections,t_connections))]

    Wr_fr = var(pm, nw, :Wr, f_bus)
    Wr_to = var(pm, nw, :Wr, t_bus)
    Wi_fr = var(pm, nw, :Wi, f_bus)
    Wi_to = var(pm, nw, :Wi, t_bus)
    Xtr = var(pm, nw, :Xtr, trans_id)
    Xti = var(pm, nw, :Xti, trans_id)
    CCtr = var(pm, nw, :CCtr, trans_id)
    CCti = var(pm, nw, :CCti, trans_id)

    P_fr = var(pm, nw, :Pt, f_idx)
    P_to = var(pm, nw, :Pt, t_idx)
    Q_fr = var(pm, nw, :Qt, f_idx)
    Q_to = var(pm, nw, :Qt, t_idx)
    p_fr = [var(pm, nw, :pt, f_idx)[c] for c in f_connections]
    p_to = [var(pm, nw, :pt, t_idx)[c] for c in t_connections]
    q_fr = [var(pm, nw, :qt, f_idx)[c] for c in f_connections]
    q_to = [var(pm, nw, :qt, t_idx)[c] for c in t_connections]

    for (f_idx,fc) in enumerate(f_connections)
        for (t_idx,tc) in enumerate(t_connections)
            JuMP.@constraint(pm.model, Wr_fr[fc,tc]-Wr_fr[fc,next[tc]]-Wr_fr[next[fc],tc]+Wr_fr[next[fc],next[tc]] == (pol*tm_scale)^2*tm[f_idx]*tm[t_idx]*Wr_to[fc,tc])
            JuMP.@constraint(pm.model, Wi_fr[fc,tc]-Wi_fr[fc,next[tc]]-Wi_fr[next[fc],tc]+Wi_fr[next[fc],next[tc]] == (pol*tm_scale)^2*tm[f_idx]*tm[t_idx]*Wi_to[fc,tc])
        end

        # conservation of current
        JuMP.@constraint(pm.model, sum(P_fr[f_idx,:]) == 0)
        JuMP.@constraint(pm.model, sum(Q_fr[f_idx,:]) == 0)
        JuMP.@constraint(pm.model, sum(P_to[f_idx,:]) == 0)
        JuMP.@constraint(pm.model, sum(Q_to[f_idx,:]) == 0)
    end

    Tt = [1 -1 0; 0 1 -1; -1 0 1]  # TODO
    constraint_SWL_psd(pm.model, Xtr, Xti, Wr_fr, Wi_fr, CCtr, CCti) # link W, CCt and Xt
    # define powers as affine transformations of X
    JuMP.@constraint(pm.model, P_fr .== Xtr*Tt)
    JuMP.@constraint(pm.model, Q_fr .== Xti*Tt)
    JuMP.@constraint(pm.model, LinearAlgebra.diag(Tt*Xtr) + p_to .== 0)
    JuMP.@constraint(pm.model, LinearAlgebra.diag(Tt*Xti) + q_to .== 0)
end
