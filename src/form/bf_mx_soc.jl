"""
    variable_mc_generator_power(pm::SOCUBFModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

The variable creation for generators in SOC branch flow model.
Delta generators always need an auxilary power variable (X) and current squared variable (CC) similar to delta loads.
Wye generators however, don't need any variables.
"""
function variable_mc_generator_power(pm::SOCUBFModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_generator_power_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_generator_power_imaginary(pm; nw=nw, bounded=bounded, report=report)
    # create auxilary variables for delta generators
    gen_del_ids = [id for (id, gen) in ref(pm, nw, :gen) if gen["configuration"]==DELTA]
    variable_mc_generator_power_delta_aux(pm, gen_del_ids; nw=nw)
    bounded && variable_mc_generator_current(pm, gen_del_ids; nw=nw, bounded=bounded, report=report)
end


"""
    variable_mc_generator_current(pm::SOCUBFModels, gen_ids::Vector{Int}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

For the SOC branch-flow formulation, the delta-generator needs an explicit current variable.
"""
function variable_mc_generator_current(pm::SOCUBFModels, gen_ids::Vector{Int}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    @assert(bounded)

    connections = Dict{Int,Vector{Int}}(i => gen["connections"] for (i,gen) in ref(pm, nw, :gen))
    # calculate bounds
    bound = Dict{eltype(gen_ids), Matrix{Real}}()
    for (id, gen) in ref(pm, nw, :gen)
        bus = ref(pm, nw, :bus, gen["gen_bus"])
        cmax = _calc_gen_current_max(gen, bus)
        bound[id] = cmax*cmax'
    end
    # create matrix variables
    (CCgr,CCgi) = variable_mx_hermitian(pm.model, gen_ids, connections; symm_bound=bound, name="CCg", prefix="$nw")
    # save references
    var(pm, nw)[:CCgr] = CCgr
    var(pm, nw)[:CCgi] = CCgi

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :CCgr, gen_ids, CCgr)
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :CCgi, gen_ids, CCgi)
end


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
    & {W}_{fr} = {T}_{m}{T}_{m}^{H} {W}_{to}  \\
    & {s}_{fr} + {s}_{to} = 0
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

    p_fr = [var(pm, nw, :pt, f_idx)[c] for c in f_connections]
    p_to = [var(pm, nw, :pt, t_idx)[c] for c in t_connections]
    q_fr = [var(pm, nw, :qt, f_idx)[c] for c in f_connections]
    q_to = [var(pm, nw, :qt, t_idx)[c] for c in t_connections]

    for (f_idx,fc) in enumerate(f_connections)
        if haskey(transformer,"controls")  # regcontrol settings
            w_fr = var(pm, nw, :w)[f_bus]
            w_to = var(pm, nw, :w)[t_bus]
            v_ref = transformer["controls"]["vreg"][f_idx]
            δ = transformer["controls"]["band"][f_idx]
            r = transformer["controls"]["r"][f_idx]
            x = transformer["controls"]["x"][f_idx]
        end

        for (t_idx,tc) in enumerate(t_connections)
            if tm_fixed[t_idx]
                JuMP.@constraint(pm.model, Wr_fr[f_idx,t_idx] == (pol*tm_scale)^2*tm[f_idx]*tm[t_idx]*Wr_to[f_idx,t_idx])
                JuMP.@constraint(pm.model, Wi_fr[f_idx,t_idx] == (pol*tm_scale)^2*tm[f_idx]*tm[t_idx]*Wi_to[f_idx,t_idx])
            else
                PolyhedralRelaxations.construct_univariate_relaxation!(pm.model, x->x^2, tm[idx], tmsqr[idx], [JuMP.has_lower_bound(tm[idx]) ? JuMP.lower_bound(tm[idx]) : 0.9, JuMP.has_upper_bound(tm[idx]) ? JuMP.upper_bound(tm[idx]) : 1.1], false)
                tmsqr_Wr = JuMP.@variable(pm.model, base_name="$(nw)_tmsqr_Wr_to_$(trans_id)_$(f_bus)_$(fc)_$(t_bus)_$(tc)")
                tmsqr_Wi = JuMP.@variable(pm.model, base_name="$(nw)_tmsqr_Wi_to_$(trans_id)_$(f_bus)_$(fc)_$(t_bus)_$(tc)")
                PolyhedralRelaxations.construct_bilinear_relaxation!(pm.model, tmsqr[idx], Wr_to[f_idx,t_idx], tmsqr_Wr, [JuMP.lower_bound(tmsqr[idx]), JuMP.upper_bound(tmsqr[idx])], [JuMP.has_lower_bound(Wr_to[f_idx,t_idx]) ? JuMP.lower_bound(Wr_to[f_idx,t_idx]) : -(1.1^2), JuMP.has_upper_bound(Wr_to[f_idx,t_idx]) ? JuMP.upper_bound(Wr_to[f_idx,t_idx]) : 1.1^2])
                PolyhedralRelaxations.construct_bilinear_relaxation!(pm.model, tmsqr[idx], Wi_to[f_idx,t_idx], tmsqr_Wi, [JuMP.lower_bound(tmsqr[idx]), JuMP.upper_bound(tmsqr[idx])], [JuMP.has_lower_bound(Wi_to[f_idx,t_idx]) ? JuMP.lower_bound(Wi_to[f_idx,t_idx]) : -(1.1^2), JuMP.has_upper_bound(Wi_to[f_idx,t_idx]) ? JuMP.upper_bound(Wi_to[f_idx,t_idx]) : 1.1^2])

                JuMP.@constraint(pm.model, Wr_fr[f_idx,t_idx] == (pol*tm_scale)^2*tmsqr_Wr_to)
                JuMP.@constraint(pm.model, Wi_fr[f_idx,t_idx] == (pol*tm_scale)^2*tmsqr_Wi_to)

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
    Q_fr = var(pm, nw, :Qt, f_idx)
    p_to = [var(pm, nw, :pt, t_idx)[c] for c in t_connections]
    q_to = [var(pm, nw, :qt, t_idx)[c] for c in t_connections]

    for (f_idx,fc) in enumerate(f_connections)
        for (t_idx,tc) in enumerate(t_connections)
            JuMP.@constraint(pm.model, Wr_fr[fc,tc]-Wr_fr[fc,next[tc]]-Wr_fr[next[fc],tc]+Wr_fr[next[fc],next[tc]] == (pol*tm_scale)^2*tm[f_idx]*tm[t_idx]*Wr_to[fc,tc])
            JuMP.@constraint(pm.model, Wi_fr[fc,tc]-Wi_fr[fc,next[tc]]-Wi_fr[next[fc],tc]+Wi_fr[next[fc],next[tc]] == (pol*tm_scale)^2*tm[f_idx]*tm[t_idx]*Wi_to[fc,tc])
        end
    end

    Tt = [1 -1 0; 0 1 -1; -1 0 1]  # TODO
    constraint_SWL_psd(pm.model, Xtr, Xti, Wr_fr, Wi_fr, CCtr, CCti) # link W, CCt and Xt
    # define powers as affine transformations of X
    JuMP.@constraint(pm.model, P_fr .== Xtr*Tt)
    JuMP.@constraint(pm.model, Q_fr .== Xti*Tt)
    JuMP.@constraint(pm.model, LinearAlgebra.diag(Tt*Xtr) + p_to .== 0)
    JuMP.@constraint(pm.model, LinearAlgebra.diag(Tt*Xti) + q_to .== 0)
end


@doc raw"""
    constraint_mc_generator_power_delta(pm::SOCUBFModels, nw::Int, gen_id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)

Adds constraints for delta-connected generators similar to delta-connected loads (using auxilary variable X).

```math
\begin{align}
&\text{Three-phase delta transformation matrix: }  T^\Delta = \begin{bmatrix}\;\;\;1 & -1 & \;\;0\\ \;\;\;0 & \;\;\;1 & -1\\ -1 & \;\;\;0 & \;\;\;1\end{bmatrix} \\
&\text{Single-phase delta transformation matrix (triple nodes): }  T^\Delta = \begin{bmatrix}\;1 & -1 \end{bmatrix} \\
&\text{Line-neutral generation power: }  S_{bus} = diag(T^\Delta X_g) \\
&\text{Line-line generation power: }  S^\Delta = diag(X_g T^\Delta)
\end{align}
```
"""
function constraint_mc_generator_power_delta(pm::SOCUBFModels, nw::Int, gen_id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    gen = ref(pm, nw, :gen, gen_id)
    bus_id = gen["gen_bus"]
    bus = ref(pm, nw, :bus, bus_id)
    terminals = bus["terminals"]
    is_triplex = length(connections)<3
    conn_bus = is_triplex ? bus["terminals"] : connections

    pg = var(pm, nw, :pg, gen_id)
    qg = var(pm, nw, :qg, gen_id)
    Xgr = var(pm, nw, :Xgr, gen_id)
    Xgi = var(pm, nw, :Xgi, gen_id)
    CCgr = var(pm, nw, :CCgr, gen_id)
    CCgi = var(pm, nw, :CCgi, gen_id)
    Wr = var(pm, nw, :Wr, bus_id)[[findfirst(isequal(c), terminals) for c in conn_bus],[findfirst(isequal(c), terminals) for c in conn_bus]]
    Wi = var(pm, nw, :Wi, bus_id)[[findfirst(isequal(c), terminals) for c in conn_bus],[findfirst(isequal(c), terminals) for c in conn_bus]]

    Tg = is_triplex ? [1 -1] : [1 -1 0; 0 1 -1; -1 0 1]  # TODO
    constraint_SWL_psd(pm.model, Xgr, Xgi, Wr, Wi, CCgr, CCgi)
    # define pg/qg and pg_bus/qg_bus as affine transformations of X
    JuMP.@constraint(pm.model, pg .== LinearAlgebra.diag(Tg*Xgr))
    JuMP.@constraint(pm.model, qg .== LinearAlgebra.diag(Tg*Xgi))

    pg_bus = LinearAlgebra.diag(Xgr*Tg)
    qg_bus = LinearAlgebra.diag(Xgi*Tg)
    pg_bus = JuMP.Containers.DenseAxisArray(pg_bus, conn_bus)
    qg_bus = JuMP.Containers.DenseAxisArray(qg_bus, conn_bus)
    var(pm, nw, :pg_bus)[gen_id] = pg_bus
    var(pm, nw, :qg_bus)[gen_id] = qg_bus

    if report
        sol(pm, nw, :gen, gen_id)[:pg_bus] = pg_bus
        sol(pm, nw, :gen, gen_id)[:qg_bus] = qg_bus
    end
end
