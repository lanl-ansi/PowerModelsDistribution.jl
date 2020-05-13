import LinearAlgebra: diag, diagm


""
function variable_mc_branch_current(pm::AbstractUBFModels; kwargs...)
    constraint_mc_branch_current_series_product_hermitian(pm; kwargs...)
end


""
function variable_mc_bus_voltage(pm::AbstractUBFModels; kwargs...)
    variable_mc_bus_voltage_prod_hermitian(pm; kwargs...)

    nw = get(kwargs, :nw, pm.cnw)
    allbuses = Set(ids(pm, nw, :bus))
    startingbuses = Set(i for (l,i,j)  in ref(pm, nw, :arcs_from))
    leafnodes = setdiff(allbuses, startingbuses)
    for i in leafnodes
        constraint_mc_voltage_psd(pm, nw, i)
    end
end


""
function variable_mc_bus_voltage_prod_hermitian(pm::AbstractUBFModels; n_cond::Int=3, nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    bus_ids = collect(ids(pm, nw, :bus))

    if bounded
        # get bounds
        vmax = Dict([(id, ref(pm, nw, :bus, id, "vmax")) for id in bus_ids])
        vmin = Dict([(id, ref(pm, nw, :bus, id, "vmin")) for id in bus_ids])
        # create bounded Hermitian matrix variables
        (Wr,Wi) = variable_mx_hermitian(pm.model, bus_ids, n_cond;
            sqrt_upper_bound=vmax, sqrt_lower_bound=vmin, name="W", prefix="$nw")
    else
        # create unbounded Hermitian matrix variables
        (Wr,Wi) = variable_mx_hermitian(pm.model, bus_ids, n_cond;
            set_lower_bound_diag_to_zero=true, name="W", prefix="$nw")
    end
    v_start = exp.((im*2*pi/3).*[0; -1; 1]) #TODO this should be made more generic eventually
    W_start = v_start*v_start'
    for (id,_) in Wr
        for i in 1:3
            for j in 1:i
                JuMP.set_start_value(Wr[id][i,j], real.(W_start)[i,j])
                if j<i
                    Wi_ij = collect(keys(Wi[id][i,j].terms))[1]
                    JuMP.set_start_value(Wi_ij, imag.(W_start)[i,j])
                end
            end
        end
    end

    # save references in dict
    var(pm, nw)[:Wr] = Wr
    var(pm, nw)[:Wi] = Wi
    # maintain compatibility
    var(pm, nw)[:w] = Dict{Int, Any}([(id, diag(Wr[id])) for id in bus_ids])

    report && _IM.sol_component_value(pm, nw, :bus, :Wr, ids(pm, nw, :bus), Wr)
    report && _IM.sol_component_value(pm, nw, :bus, :Wi, ids(pm, nw, :bus), Wi)
end


""
function constraint_mc_branch_current_series_product_hermitian(pm::AbstractUBFModels; n_cond::Int=3, nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    branches = ref(pm, nw, :branch)
    buses = ref(pm, nw, :bus)

    branch_ids = collect(keys(branches))

    if bounded
        # calculate max series current for each branch
        cmax = Dict{eltype(branch_ids), Array{Real,1}}()
        for (key, branch) in branches
            bus_fr = buses[branch["f_bus"]]
            bus_to = buses[branch["t_bus"]]
            cmax[key] = _calc_branch_series_current_max(branch, bus_fr, bus_to)
        end
        # create matrix variables
        (Lr,Li) = variable_mx_hermitian(pm.model, branch_ids, n_cond;
            sqrt_upper_bound=cmax, set_lower_bound_diag_to_zero=true,
            name="CC", prefix="$nw")
    else
        (Lr,Li) = variable_mx_hermitian(pm.model, branch_ids, n_cond;
            set_lower_bound_diag_to_zero=true, name="CC", prefix="$nw")
    end

    for (id,L) in Lr
        JuMP.set_start_value.(LinearAlgebra.diag(Lr[id]), 0.01)
    end

    # save reference
    var(pm, nw)[:CCr] = Lr
    var(pm, nw)[:CCi] = Li
    var(pm, nw)[:cm] = Dict([(id, diag(Lr[id])) for id in branch_ids])

    report && _IM.sol_component_value(pm, nw, :branch, :CCr, ids(pm, nw, :branch), Lr)
    report && _IM.sol_component_value(pm, nw, :branch, :CCi, ids(pm, nw, :branch), Li)
end


""
function variable_mc_branch_power(pm::AbstractUBFModels; n_cond::Int=3, nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    @assert n_cond<=5

    # calculate S bound
    branch_arcs = ref(pm, nw, :arcs)
    if bounded
        bound = Dict{eltype(branch_arcs), Array{Real,2}}()
        for (br, branch) in ref(pm, nw, :branch)
            bus_fr = ref(pm, nw, :bus, branch["f_bus"])
            bus_to = ref(pm, nw, :bus, branch["t_bus"])

            smax_fr = _calc_branch_power_max(branch, bus_fr)
            smax_to = _calc_branch_power_max(branch, bus_to)
            cmax_fr, cmax_to = _calc_branch_current_max_frto(branch, bus_fr, bus_to)

            tuple_fr = (br, bus_fr["index"], bus_to["index"])
            tuple_to = (br, bus_to["index"], bus_fr["index"])

            bound[tuple_fr] = bus_fr["vmax"].*cmax_fr'
            bound[tuple_to] = bus_to["vmax"].*cmax_to'

            for c in 1:length(smax_fr)
                bound[tuple_fr][c,c] = smax_fr[c]
                bound[tuple_to][c,c] = smax_to[c]
            end
        end
        # create matrix variables
        (P,Q) = variable_mx_complex(pm.model, branch_arcs, n_cond, n_cond;
            symm_bound=bound, name=("P", "Q"), prefix="$nw")
    else
        (P,Q) = variable_mx_complex(pm.model, branch_arcs, n_cond, n_cond;
            name=("P", "Q"), prefix="$nw")
    end
    # save reference
    var(pm, nw)[:P] = P
    var(pm, nw)[:Q] = Q

    var(pm, nw)[:p] = Dict([(id,diag(P[id])) for id in branch_arcs])
    var(pm, nw)[:q] = Dict([(id,diag(Q[id])) for id in branch_arcs])

    report && _IM.sol_component_value_edge(pm, nw, :branch, :Pf, :Pt, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), P)
    report && _IM.sol_component_value_edge(pm, nw, :branch, :Qf, :Qt, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), Q)
end


"Defines branch flow model power flow equations"
function constraint_mc_power_losses(pm::AbstractUBFModels, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to, tm)
    P_to = var(pm, n, :P)[t_idx]
    Q_to = var(pm, n, :Q)[t_idx]

    P_fr = var(pm, n, :P)[f_idx]
    Q_fr = var(pm, n, :Q)[f_idx]

    Wr_to = var(pm, n, :Wr)[t_bus]
    Wr_fr = var(pm, n, :Wr)[f_bus]

    Wi_to = var(pm, n, :Wi)[t_bus]
    Wi_fr = var(pm, n, :Wi)[f_bus]

    CCr =  var(pm, n, :CCr)[i]
    CCi =  var(pm, n, :CCi)[i]

    JuMP.@constraint(pm.model, P_fr + P_to .==  Wr_fr*(g_sh_fr)' + Wi_fr*(b_sh_fr)' + r*CCr - x*CCi +  Wr_to*(g_sh_to)'  + Wi_to*(b_sh_to)')
    JuMP.@constraint(pm.model, Q_fr + Q_to .==  Wi_fr*(g_sh_fr)' - Wr_fr*(b_sh_fr)' + x*CCr + r*CCi +  Wi_to*(g_sh_to)'  - Wr_to*(b_sh_to)')
end


""
function constraint_mc_theta_ref(pm::AbstractUBFModels, n::Int, i::Int, va_ref)
    nconductors = length(conductor_ids(pm))

    Wr = var(pm, n, :Wr)[i]
    Wi = var(pm, n, :Wi)[i]

    beta = exp.(im.*va_ref)
    gamma = beta*beta'

    Wr_ref = real(gamma).*Wr[1,1]
    Wi_ref = imag(gamma).*Wr[1,1]
    JuMP.@constraint(pm.model, diag(Wr)[2:nconductors]        .== diag(Wr_ref)[2:nconductors]) # first equality is implied
    JuMP.@constraint(pm.model, _mat2utrivec!(Wr) .== _mat2utrivec!(Wr_ref))
    JuMP.@constraint(pm.model, _mat2utrivec!(Wi) .== _mat2utrivec!(Wi_ref))
end


"Defines voltage drop over a branch, linking from and to side voltage"
function constraint_mc_model_voltage_magnitude_difference(pm::AbstractUBFModels, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
    Wr_fr = var(pm, n, :Wr)[f_bus]
    Wi_fr = var(pm, n, :Wi)[f_bus]

    Wr_to = var(pm, n, :Wr)[t_bus]
    Wi_to = var(pm, n, :Wi)[t_bus]

    p_fr = var(pm, n, :P)[f_idx]
    q_fr = var(pm, n, :Q)[f_idx]

    p_s_fr = p_fr - (Wr_fr*(g_sh_fr)' + Wi_fr*(b_sh_fr)')
    q_s_fr = q_fr - (Wi_fr*(g_sh_fr)' - Wr_fr*(b_sh_fr)')

    CCr =  var(pm, n, :CCr)[i]
    CCi =  var(pm, n, :CCi)[i]

    #KVL over the line:
    JuMP.@constraint(pm.model, diag(Wr_to) .== diag(
    Wr_fr   - p_s_fr  *r' - q_s_fr*x'        - r*p_s_fr'    - x*q_s_fr'
    + r*CCr*r' - x     *CCi*r' + x*CCr *x' + r*CCi *x'))
    JuMP.@constraint(pm.model, _mat2utrivec!(Wr_to) .== _mat2utrivec!(
    Wr_fr   - p_s_fr  *r' - q_s_fr*x'        - r*p_s_fr'    - x*q_s_fr'
    + r*CCr*r' - x     *CCi*r' + x*CCr *x' + r*CCi *x'))
    JuMP.@constraint(pm.model, _mat2utrivec!(Wi_to) .== _mat2utrivec!(
    Wi_fr   - q_s_fr  *r' + p_s_fr*x'        - x*p_s_fr'    + r*q_s_fr'
    + x*CCr*r' + r     *CCi*r' - r*CCr *x' + x*CCi *x'))
end


"""
For the matrix KCL formulation, the generator needs an explicit current and
power variable.
"""
function variable_mc_gen_power_setpoint(pm::SDPUBFKCLMXModel; kwargs...)
    variable_mc_gen_current(pm; kwargs...)
    variable_mc_gen_power(pm; kwargs...)
end


"""
For the matrix KCL formulation, the generator needs an explicit power
variable.
"""
function variable_mc_gen_power(pm::SDPUBFKCLMXModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    @assert(bounded)

    gen_ids = collect(ids(pm, nw, :gen))
    ncnds = length(conductor_ids(pm, nw))

    # calculate bounds for matrix variable
    Pg_min = Dict{eltype(gen_ids), Array{Real,2}}()
    Pg_max = Dict{eltype(gen_ids), Array{Real,2}}()
    Qg_min = Dict{eltype(gen_ids), Array{Real,2}}()
    Qg_max = Dict{eltype(gen_ids), Array{Real,2}}()
    for (id, gen) in ref(pm, nw, :gen)
        bus = ref(pm, nw, :bus, gen["gen_bus"])
        vmax = get(bus, "vmax", fill(Inf, 3))
        cmax = _calc_gen_current_max(gen, bus)
        S_bound = _mat_mult_rm_nan(vmax, cmax')

        Pg_min[id] = Qg_min[id] = -S_bound
        Pg_max[id] = Qg_max[id] =  S_bound

        pmin = get(gen, "pmin", fill(-Inf, 3))
        pmax = get(gen, "pmax", fill( Inf, 3))
        qmin = get(gen, "qmin", fill(-Inf, 3))
        qmax = get(gen, "qmax", fill( Inf, 3))

        for c in 1:3
            Pg_min[id][c,c] = max(pmin[c], Pg_min[id][c,c])
            Pg_max[id][c,c] = min(pmax[c], Pg_max[id][c,c])
            Qg_min[id][c,c] = max(qmin[c], Qg_min[id][c,c])
            Qg_max[id][c,c] = min(qmax[c], Qg_max[id][c,c])
        end
    end

    # create matrix variables
    Pg = variable_mx_real(pm.model, gen_ids, ncnds, ncnds;
        lower_bound=Pg_min, upper_bound=Pg_max,
        name="Pg", prefix="$nw")
    Qg = variable_mx_real(pm.model, gen_ids, ncnds, ncnds;
        lower_bound=Qg_min, upper_bound=Qg_max,
        name="Qg", prefix="$nw")

    # save references
    var(pm, nw)[:Pg_bus] = Pg
    var(pm, nw)[:Qg_bus] = Qg
    var(pm, nw)[:pg] = Dict{Int, Any}([(id, diag(Pg[id])) for id in gen_ids])
    var(pm, nw)[:qg] = Dict{Int, Any}([(id, diag(Qg[id])) for id in gen_ids])

    report && _IM.sol_component_value(pm, nw, :gen, :Pg_bus, ids(pm, nw, :gen), Pg)
    report && _IM.sol_component_value(pm, nw, :gen, :Qg_bus, ids(pm, nw, :gen), Qg)
end


"""
For the matrix KCL formulation, the generator needs an explicit current
variable.
"""
function variable_mc_gen_current(pm::AbstractUBFModels; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    @assert(bounded)

    gen_ids = collect(ids(pm, nw, :gen))
    ncnds = length(conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(gen_ids), Array{Real,2}}()
    for (id, gen) in ref(pm, nw, :gen)
        bus = ref(pm, nw, :bus, gen["gen_bus"])
        cmax = _calc_gen_current_max(gen, bus)
        bound[id] = cmax*cmax'
    end
    # create matrix variables
    (CCgr,CCgi) = variable_mx_hermitian(pm.model, gen_ids, ncnds; symm_bound=bound, name="CCg", prefix="$nw")
    # save references
    var(pm, nw)[:CCgr] = CCgr
    var(pm, nw)[:CCgi] = CCgi

    report && _IM.sol_component_value(pm, nw, :gen, :CCgr, ids(pm, nw, :gen), CCgr)
    report && _IM.sol_component_value(pm, nw, :gen, :CCgi, ids(pm, nw, :gen), CCgi)
end


"""
The variable creation for the loads is rather complicated because Expressions
are used wherever possible instead of explicit variables.
Delta loads always need a current variable and auxilary power variable (X), and
all other load model variables are then linear transformations of these
(linear Expressions).
Wye loads however, don't need any variables when the load is modelled as
constant power or constant impedance. In all other cases (e.g. when a cone is
used to constrain the power), variables need to be created.
"""
function variable_mc_load_setpoint(pm::AbstractUBFModels; nw=pm.cnw)
    load_wye_ids = [id for (id, load) in ref(pm, nw, :load) if load["configuration"]==WYE]
    load_del_ids = [id for (id, load) in ref(pm, nw, :load) if load["configuration"]==DELTA]
    load_cone_ids = [id for (id, load) in ref(pm, nw, :load) if _check_load_needs_cone(load)]
    # create dictionaries
    var(pm, nw)[:pd_bus] = Dict()
    var(pm, nw)[:qd_bus] = Dict()
    var(pm, nw)[:pd] = Dict()
    var(pm, nw)[:qd] = Dict()
    # now, create auxilary power variable X for delta loads
    variable_mc_load_power_delta_aux(pm, load_del_ids)
    # only delta loads need a current variable
    variable_mc_load_current(pm, load_del_ids)
    # for wye loads with a cone inclusion constraint, we need to create a variable
    variable_mc_load_power(pm, intersect(load_wye_ids, load_cone_ids))

end


"""
The variable creation for the loads is rather complicated because Expressions
are used wherever possible instead of explicit variables.
All loads need a current variable; for wye loads, this variable will be in the
wye reference frame whilst for delta currents it will be in the delta reference
frame.
All loads need variables for the off-diagonals of the nodal power variables. In
some cases, the diagonals elements can be created as Expressions.
Delta loads only need a current variable and auxilary power variable (X), and
all other load model variables are then linear transformations of these
(linear Expressions).
"""
function variable_mc_load_setpoint(pm::SDPUBFKCLMXModel; nw=pm.cnw)
    load_wye_ids = [id for (id, load) in ref(pm, nw, :load) if load["configuration"]==WYE]
    load_del_ids = [id for (id, load) in ref(pm, nw, :load) if load["configuration"]==DELTA]
    load_cone_ids = [id for (id, load) in ref(pm, nw, :load) if _check_load_needs_cone(load)]
    # create dictionaries
    var(pm, nw)[:Pd_bus] = Dict{Int, Any}()
    var(pm, nw)[:Qd_bus] = Dict{Int, Any}()
    var(pm, nw)[:pd] = Dict{Int, Any}()
    var(pm, nw)[:qd] = Dict{Int, Any}()    # now, create auxilary power variable X for delta loads
    variable_mc_load_power_delta_aux(pm, load_del_ids)
    # all loads need a current variable now
    variable_mc_load_current(pm, collect(ids(pm, nw, :load)))
    # for all wye-connected loads, we need variables for the off-diagonals of Pd/Qd
    variable_mc_load_power_bus(pm, load_wye_ids)
    # for wye loads with a cone inclusion constraint, we need to create a variable for the diagonal
    variable_mc_load_power(pm, intersect(load_wye_ids, load_cone_ids))
end


"""
These variables reflect the power consumed by the load, NOT the power injected
into the bus nodes; these variables only coincide for wye-connected loads
with a grounded neutral.
"""
function variable_mc_load_power(pm::AbstractUBFModels, load_ids::Array{Int,1}; nw=pm.cnw, bounded::Bool=true, report::Bool=true)
    @assert(bounded)
    # calculate bounds for all loads
    pmin = Dict()
    pmax = Dict()
    qmin = Dict()
    qmax = Dict()
    for id in load_ids
        load = ref(pm, nw, :load, id)
        bus = ref(pm, nw, :bus, load["load_bus"])
        pmin[id], pmax[id], qmin[id], qmax[id] = _calc_load_pq_bounds(load, bus)
    end

    # create variables
    ncnds = length(conductor_ids(pm, nw))

    pd = Dict(i => JuMP.@variable(pm.model,
        [c in 1:ncnds], base_name="$(nw)_pd_$(i)"
        ) for i in load_ids
    )
    qd = Dict(i => JuMP.@variable(pm.model,
        [c in 1:ncnds], base_name="$(nw)_qd_$(i)"
        ) for i in load_ids
    )

    if bounded
        for i in load_ids
            load = ref(pm, nw, :load, i)
            bus = ref(pm, nw, :bus, load["load_bus"])
            pmin, pmax, qmin, qmax = _calc_load_pq_bounds(load, bus)
            set_lower_bound.(pd[i], pmin)
            set_upper_bound.(pd[i], pmax)
            set_lower_bound.(qd[i], qmin)
            set_upper_bound.(qd[i], qmax)
        end
    end

    #store in dict, but do not overwrite
    for i in load_ids
        var(pm, nw)[:pd][i] = pd[i]
        var(pm, nw)[:qd][i] = qd[i]
    end

    report && _IM.sol_component_value(pm, nw, :load, :pd, load_ids, pd)
    report && _IM.sol_component_value(pm, nw, :load, :qd, load_ids, qd)
end


"""
The bus qualifier denotes that this is the power withdrawn at the bus; Only for
grounded wye-connected loads, this is the same as the power consumed by the
multi-phase load. The off-diagonals only need to be created for the matrix KCL
formulation.
"""
function variable_mc_load_power_bus(pm::SDPUBFKCLMXModel, load_ids::Array{Int,1}; nw=pm.cnw, bounded::Bool=true, report::Bool=true)
    @assert(bounded)
    ncnds = length(conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(load_ids), Array{Real,2}}()
    for id in load_ids
        load = ref(pm, nw, :load, id)
        @assert(load["configuration"]==WYE)
        bus = ref(pm, nw, :bus, load["load_bus"])
        cmax = _calc_load_current_max(load, bus)
        bound[id] = bus["vmax"]*cmax'
    end
    # create matrix variables
    (Pd_bus,Qd_bus) = variable_mx_complex_with_diag(pm.model, load_ids, ncnds; symm_bound=bound, name=("Pd_bus", "Qd_bus"), prefix="$nw")
    for id in load_ids
        var(pm, nw, :Pd_bus)[id] = Pd_bus[id]
        var(pm, nw, :Qd_bus)[id] = Qd_bus[id]
    end

    report && _IM.sol_component_value(pm, nw, :load, :Pd_bus, load_ids, Pd_bus)
    report && _IM.sol_component_value(pm, nw, :load, :Qd_bus, load_ids, Qd_bus)
end


"""
Creates power matrix variable X for delta windings; this defines both the
wye-side power Sy and the delta-side power Sd through the lin. transformations
Sy = X.Td, Sd = Td.X with Td=[1 -1 0; 0 1 -1; -1 0 1]

See the paper by Zhao et al. for the first convex relaxation of delta transformations.
@INPROCEEDINGS{zhao_optimal_2017,
	author={C. Zhao, E. Dall'Anese and S. Low},
	booktitle={IREP 2017 Bulk Power Systems Dynamics and Control Symposium},
	title={{Optimal Power Flow in Multiphase Radial Networks with Delta Connections}},
	year={2017},
	month={},
    url={https://www.nrel.gov/docs/fy18osti/67852.pdf}
}

See upcoming paper for discussion of bounds. [reference added when accepted]
"""
function variable_mc_load_power_delta_aux(pm::AbstractUBFModels, load_ids::Array{Int,1}; nw=pm.cnw, eps=0.1, bounded::Bool=true, report::Bool=true)
    @assert(bounded)
    ncnds = length(conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(load_ids), Array{Real,2}}()
    for id in load_ids
        load = ref(pm, nw, :load, id)
        bus_id = load["load_bus"]
        bus = ref(pm, nw, :bus, bus_id)
        cmax = _calc_load_current_max(load, bus)
        bound[id] = bus["vmax"]*cmax'
    end
    # create matrix variables
    (Xdr,Xdi) = variable_mx_complex(pm.model, load_ids, ncnds, ncnds;
        symm_bound=bound, name="Xd", prefix="$nw")
    # save references
    var(pm, nw)[:Xdr] = Xdr
    var(pm, nw)[:Xdi] = Xdi

    report && _IM.sol_component_value(pm, nw, :load, :Xdr, load_ids, Xdr)
    report && _IM.sol_component_value(pm, nw, :load, :Xdi, load_ids, Xdi)
end


"""
All loads need a current variable; for wye loads, this variable will be in the
wye reference frame whilst for delta currents it will be in the delta reference
frame.
"""
function variable_mc_load_current(pm::AbstractUBFModels, load_ids::Array{Int,1}; nw=pm.cnw, bounded::Bool=true, report::Bool=true)
    @assert(bounded)

    ncnds = length(conductor_ids(pm, nw))
    # calculate bounds
    cmin = Dict{eltype(load_ids), Array{Real,1}}()
    cmax = Dict{eltype(load_ids), Array{Real,1}}()
    for (id, load) in ref(pm, nw, :load)
        bus_id = load["load_bus"]
        bus = ref(pm, nw, :bus, bus_id)
        cmin[id], cmax[id] = _calc_load_current_magnitude_bounds(load, bus)
    end
    # create matrix variables
    (CCdr, CCdi) = variable_mx_hermitian(pm.model, load_ids, ncnds;
        sqrt_upper_bound=cmax, sqrt_lower_bound=cmin, name="CCd", prefix="$nw")
    # save references
    var(pm, nw)[:CCdr] = CCdr
    var(pm, nw)[:CCdi] = CCdi

    report && _IM.sol_component_value(pm, nw, :load, :CCdr, load_ids, CCdr)
    report && _IM.sol_component_value(pm, nw, :load, :CCdi, load_ids, CCdi)
end


"""
Link the current and power withdrawn by a generator at the bus through a PSD
constraint. The rank-1 constraint is dropped in this formulation.
"""
function constraint_mc_gen_setpoint(pm::SDPUBFKCLMXModel, gen_id::Int; nw::Int=pm.cnw)
    Pg = var(pm, nw, :Pg_bus, gen_id)
    Qg = var(pm, nw, :Qg_bus, gen_id)
    bus_id = ref(pm, nw, :gen, gen_id)["gen_bus"]
    Wr = var(pm, nw, :Wr, bus_id)
    Wi = var(pm, nw, :Wi, bus_id)
    CCgr = var(pm, nw, :CCgr, gen_id)
    CCgi = var(pm, nw, :CCgi, gen_id)
    constraint_SWL_psd(pm.model, Pg, Qg, Wr, Wi, CCgr, CCgi)
end


"""
Creates the constraints modelling the (relaxed) voltage-dependency of the
power consumed in each phase, s=p+jq. This is completely symmetrical for p and
q, with appropriate substitutions of the variables and parameters:
p->q, a->b, alpha->beta, pmin->qmin, pmax->qmax
"""
function constraint_pqw(model::JuMP.Model, w, p, a::Real, alpha::Real, wmin::Real, wmax::Real, pmin::Real, pmax::Real)
    if a==0
        JuMP.@constraint(model, p==0)
    else
        @assert(alpha>=0, "alpha has to greater than or equal to zero.")

        # CONSTANT POWER
        if alpha==0
            JuMP.@constraint(model, p==a)

        # CONSTANT IMPEDANCE
        elseif alpha==2
            JuMP.@constraint(model, p==a*w)

        # CONE INCLUSIONS
        else
            # cone inclusions have an affine over/under estimator
            # boundary
            if a>0
                l = (1/a)*(pmax-pmin)/(wmax-wmin)*(w-wmin) + pmin/a
            else
                # swap pmin and pmax if a<0, because pmin/a > pmax/a
                l = (1/a)*(pmin-pmax)/(wmax-wmin)*(w-wmin) + pmax/a
            end
            # affine overestimator
            if alpha>2
                JuMP.@constraint(model, p/a <= l)
            # affine underestimator
            elseif 0<alpha<2
                JuMP.@constraint(model, p/a >= l)
            end

            # constant current case, simplifies to a RotatedSecondOrderCone
            if alpha==1
                #       p/a <= w^(1/2)
                # <=>   (p/a)^2 <= w
                # <=>   2*(w/2)*1 >= ||p/a||^2_2
                # <=>   (w/2, 1, p/a) ∈ RotatedSecondOrderCone(3)
                JuMP.@constraint(model, [w/2, 1, p/a] in JuMP.RotatedSecondOrderCone())
            # general power cone
            elseif 0<alpha<2
                #       p/a <= w^(alpha/2)
                # <=>   w^(alpha/2) >= p/a
                # <=>   (w, 1, p/a) ∈ PowerCone(3)
                JuMP.@constraint(model, [w, 1, p/a] in MathOptInterface.PowerCone(alpha/2))
            # general power cone
            else # alpha>2
                #       p/a >= w^(alpha/2)
                # <=>   (p/a)^(2/alpha) >= w
                # <=>   (p/a, 1, w) ∈ PowerCone(3)
                JuMP.@constraint(model, [p/a, 1, w] in MathOptInterface.PowerCone(2/alpha))
            end
        end
    end
end


"""
Creates the constraints modelling the (relaxed) voltage-dependent loads.
"""
function constraint_mc_load_setpoint(pm::AbstractUBFModels, load_id::Int; nw::Int=pm.cnw, report::Bool=true)
    # shared variables and parameters
    load = ref(pm, nw, :load, load_id)
    pd0 = load["pd"]
    qd0 = load["qd"]
    bus_id = load["load_bus"]
    bus = ref(pm, nw, :bus, bus_id)
    ncnds = length(pd0)

    # calculate load params
    a, alpha, b, beta = _load_expmodel_params(load, bus)
    vmin, vmax = _calc_load_vbounds(load, bus)
    wmin = vmin.^2
    wmax = vmax.^2
    pmin, pmax, qmin, qmax = _calc_load_pq_bounds(load, bus)

    # take care of connections
    if load["configuration"]==WYE
        if load["model"]==POWER
            var(pm, nw, :pd)[load_id] = pd0
            var(pm, nw, :qd)[load_id] = qd0
        elseif load["model"]==IMPEDANCE
            w = var(pm, nw, :w)[bus_id]
            var(pm, nw, :pd)[load_id] = a.*w
            var(pm, nw, :qd)[load_id] = b.*w
        # in this case, :pd has a JuMP variable
        else
            Wr = var(pm, nw, :Wr, bus_id)
            pd = var(pm, nw, :pd)[load_id]
            qd = var(pm, nw, :qd)[load_id]
            for c in 1:ncnds
                constraint_pqw(pm.model, Wr[c,c], pd[c], a[c], alpha[c], wmin[c], wmax[c], pmin[c], pmax[c])
                constraint_pqw(pm.model, Wr[c,c], qd[c], b[c], beta[c], wmin[c], wmax[c], qmin[c], qmax[c])
            end
        end
        # :pd_bus is identical to :pd now
        var(pm, nw, :pd_bus)[load_id] = var(pm, nw, :pd)[load_id]
        var(pm, nw, :qd_bus)[load_id] = var(pm, nw, :qd)[load_id]

        ## reporting
        if report
            sol(pm, nw, :load, load_id)[:pd] = var(pm, nw, :pd)[load_id]
            sol(pm, nw, :load, load_id)[:qd] = var(pm, nw, :qd)[load_id]
            sol(pm, nw, :load, load_id)[:pd_bus] = var(pm, nw, :pd_bus)[load_id]
            sol(pm, nw, :load, load_id)[:qd_bus] = var(pm, nw, :qd_bus)[load_id]
        end
    elseif load["configuration"]==DELTA
        # link Wy, CCd and X
        Wr = var(pm, nw, :Wr, bus_id)
        Wi = var(pm, nw, :Wi, bus_id)
        CCdr = var(pm, nw, :CCdr, load_id)
        CCdi = var(pm, nw, :CCdi, load_id)
        Xdr = var(pm, nw, :Xdr, load_id)
        Xdi = var(pm, nw, :Xdi, load_id)
        Td = [1 -1 0; 0 1 -1; -1 0 1]
        constraint_SWL_psd(pm.model, Xdr, Xdi, Wr, Wi, CCdr, CCdi)
        # define pd/qd and pd_bus/qd_bus as affine transformations of X
        pd_bus = LinearAlgebra.diag(Xdr*Td)
        qd_bus = LinearAlgebra.diag(Xdi*Td)
        pd = LinearAlgebra.diag(Td*Xdr)
        qd = LinearAlgebra.diag(Td*Xdi)

        var(pm, nw, :pd_bus)[load_id] = pd_bus
        var(pm, nw, :qd_bus)[load_id] = qd_bus
        var(pm, nw, :pd)[load_id] = pd
        var(pm, nw, :qd)[load_id] = qd

        # |Vd|^2 is a linear transformation of Wr
        wd = LinearAlgebra.diag(Td*Wr*Td')
        if load["model"]==POWER
            for c in 1:ncnds
                JuMP.@constraint(pm.model, pd[c]==pd0[c])
                JuMP.@constraint(pm.model, qd[c]==qd0[c])
            end
        elseif load["model"]==IMPEDANCE
            for c in 1:ncnds
                JuMP.@constraint(pm.model, pd[c]==a[c]*wd[c])
                JuMP.@constraint(pm.model, qd[c]==b[c]*wd[c])
            end
        else
            for c in 1:ncnds
                constraint_pqw(pm.model, wd[c], pd[c], a[c], alpha[c], wmin[c], wmax[c], pmin[c], pmax[c])
                constraint_pqw(pm.model, wd[c], qd[c], b[c], beta[c],  wmin[c], wmax[c], qmin[c], qmax[c])
            end
            end

        ## reporting; for delta these are not available as saved variables!
        if report
            sol(pm, nw, :load, load_id)[:pd] = pd
            sol(pm, nw, :load, load_id)[:qd] = qd
            sol(pm, nw, :load, load_id)[:pd_bus] = pd_bus
            sol(pm, nw, :load, load_id)[:qd_bus] = qd_bus
        end
    end
end


"""
Creates the constraints modelling the (relaxed) voltage-dependent loads for
the matrix KCL formulation.
"""
function constraint_mc_load_setpoint(pm::SDPUBFKCLMXModel, load_id::Int; nw::Int=pm.cnw, report::Bool=true)
    # shared variables and parameters
    load = ref(pm, nw, :load, load_id)
    pd0 = load["pd"]
    qd0 = load["qd"]
    bus_id = load["load_bus"]
    bus = ref(pm, nw, :bus, bus_id)
    ncnds = length(pd0)

    # calculate load params
    a, alpha, b, beta = _load_expmodel_params(load, bus)
    vmin, vmax = _calc_load_vbounds(load, bus)
    wmin = vmin.^2
    wmax = vmax.^2
    pmin, pmax, qmin, qmax = _calc_load_pq_bounds(load, bus)

    # take care of connections
    Wr = var(pm, nw, :Wr, bus_id)
    Wi = var(pm, nw, :Wi, bus_id)
    CCdr = var(pm, nw, :CCdr, load_id)
    CCdi = var(pm, nw, :CCdi, load_id)

    if load["configuration"]==WYE
        if load["model"]==POWER
            var(pm, nw, :pd)[load_id] = pd0
            var(pm, nw, :qd)[load_id] = qd0
        elseif load["model"]==IMPEDANCE
            w = var(pm, nw, :w, bus_id)
            # for c in 1:ncnds
            var(pm, nw, :pd)[load_id] = a.*w
            var(pm, nw, :qd)[load_id] = b.*w
            # end
        # in this case, :pd has a JuMP variable
        else
            pd = var(pm, nw, :pd)[load_id]
            qd = var(pm, nw, :qd)[load_id]
            for c in 1:ncnds
                constraint_pqw(pm.model, Wr[c,c], pd[c], a[c], alpha[c], wmin[c], wmax[c], pmin[c], pmax[c])
                constraint_pqw(pm.model, Wr[c,c], qd[c], b[c], beta[c], wmin[c], wmax[c], qmin[c], qmax[c])
            end
        end
        # diagonal of :Pd is identical to :pd now
        Pd_bus = var(pm, nw, :Pd_bus)[load_id]
        Qd_bus = var(pm, nw, :Qd_bus)[load_id]
        for c in 1:ncnds
            Pd_bus[c,c] = var(pm, nw, :pd)[load_id][c]
            Qd_bus[c,c] = var(pm, nw, :qd)[load_id][c]
        end

    elseif load["configuration"]==DELTA
        # link Wy, CCd and X
        Xdr = var(pm, nw, :Xdr, load_id)
        Xdi = var(pm, nw, :Xdi, load_id)
        Td = [1 -1 0; 0 1 -1; -1 0 1]
        constraint_SWL_psd(pm.model, Xdr, Xdi, Wr, Wi, CCdr, CCdi)
        # define pd_bus/qd_bus and pd/qd as affine transformations of X
        Pd_bus = Xdr*Td
        Qd_bus = Xdi*Td
        pd = LinearAlgebra.diag(Td*Xdr)
        qd = LinearAlgebra.diag(Td*Xdi)

        var(pm, nw, :Pd_bus)[load_id] = Pd_bus
        var(pm, nw, :Qd_bus)[load_id] = Qd_bus
        var(pm, nw, :pd)[load_id] = pd
        var(pm, nw, :qd)[load_id] = qd

        # |Vd|^2 is a linear transformation of Wr
        wd = LinearAlgebra.diag(Td*Wr*Td')
        if load["model"]==POWER
            for c in 1:ncnds
                JuMP.@constraint(pm.model, pd[c]==pd0[c])
                JuMP.@constraint(pm.model, qd[c]==qd0[c])
            end
        elseif load["model"]==IMPEDANCE
            for c in 1:ncnds
                JuMP.@constraint(pm.model, pd[c]==a[c]*wd[c])
                JuMP.@constraint(pm.model, qd[c]==b[c]*wd[c])
            end
        else
            for c in 1:ncnds
                constraint_pqw(pm.model, wd[c], pd[c], a[c], alpha[c], wmin[c], wmax[c], pmin[c], pmax[c])
                constraint_pqw(pm.model, wd[c], qd[c], b[c], beta[c], wmin[c], wmax[c], qmin[c], qmax[c])
            end
        end
    end
end


"""
Take a multi-conductor voltage variable V and a current variable I. The
associated power is then defined as S = V*I^H
Define the lifted variables as W and L as
W = V*V^H, L = I*I^H
Then, it is equally valid that
[W S; S^H L] ∈ PSDCone, rank([W S; S^H L])=1
This function adds this PSD constraint for the rectangular coordinates of S, W
and L.
"""
function constraint_SWL_psd(model::JuMP.Model, P, Q, Wr, Wi, L_re, L_im)
    M_re = [Wr P; P' L_re]
    M_im = [Wi Q; -Q' L_im]
    constraint_M_psd(model, M_re, M_im)
end


"""
For rectangular coordinates of a complex matrix M=M_re+im*M_im,
this function applies constraints equivalent to requiring that M itself is PSD.
"""
function constraint_M_psd(model::JuMP.Model, M_re, M_im)
    JuMP.@constraint(model, [M_re -M_im; M_im M_re] in JuMP.PSDCone())
end


"""
Shunt handling in matrix form:
I = Y.U
S = U.I' = U.(Y.U)' = U.U'.Y' = W.Y'
  = (Wr+j.Wi)(G+jB)' = (Wr+j.Wi)(G'-j.B') = (Wr.G'+Wi.B')+j(-Wr.B'+Wi.G')
P =  Wr.G'+Wi.B'
Q = -Wr.B'+Wi.G'
"""
function constraint_mc_load_power_balance(pm::KCLMXModels, n::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_gs, bus_bs)
    Wr = var(pm, n, :Wr, i)
    Wi = var(pm, n, :Wi, i)

    P = get(var(pm, n), :P, Dict()); _PM._check_var_keys(P, bus_arcs, "active power", "branch")
    Q = get(var(pm, n), :Q, Dict()); _PM._check_var_keys(Q, bus_arcs, "reactive power", "branch")
    Pg = get(var(pm, n), :Pg_bus, Dict()); _PM._check_var_keys(Pg, bus_gens, "active power", "generator")
    Qg = get(var(pm, n), :Qg_bus, Dict()); _PM._check_var_keys(Qg, bus_gens, "reactive power", "generator")
    Pd = get(var(pm, n), :Pd_bus, Dict()); _PM._check_var_keys(Pd, bus_loads, "active power", "load")
    Qd = get(var(pm, n), :Qd_bus, Dict()); _PM._check_var_keys(Qd, bus_loads, "reactive power", "load")

    cnds = conductor_ids(pm; nw=n)
    ncnds = length(cnds)

    Gt = isempty(bus_gs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_gs))
    Bt = isempty(bus_bs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_bs))

    cstr_p = JuMP.@constraint(pm.model,
          sum(P[a] for a in bus_arcs)
        # + sum(Psw[a_sw] for a_sw in bus_arcs_sw)
        # + sum(Pt[a_trans] for a_trans in bus_arcs_trans)
        .==
          sum(Pg[g] for g in bus_gens)
        # - sum(ps[s] for s in bus_storage)
        - sum(Pd[d] for d in bus_loads)
        - diag(Wr*Gt'+Wi*Bt')
    )

    cstr_q = JuMP.@constraint(pm.model,
          sum(Q[a] for a in bus_arcs)
        # + sum(diag(Qsw[a_sw]) for a_sw in bus_arcs_sw)
        # + sum(diag(Qt[a_trans]) for a_trans in bus_arcs_trans)
        .==
          sum(Qg[g] for g in bus_gens)
        # - sum(qs[s] for s in bus_storage)
        - sum(Qd[d] for d in bus_loads)
        - diag(-Wr*Bt'+Wi*Gt')
    )

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end