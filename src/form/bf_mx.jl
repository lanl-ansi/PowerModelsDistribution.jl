import LinearAlgebra: diag, diagm


""
function variable_mc_branch_current(pm::AbstractUBFModels; kwargs...)
    constraint_mc_branch_current_series_product_hermitian(pm; kwargs...)
end


""
function variable_mc_bus_voltage(pm::AbstractUBFModels; kwargs...)
    variable_mc_bus_voltage_prod_hermitian(pm; kwargs...)

    nw = get(kwargs, :nw, nw_id_default)
    allbuses = Set(ids(pm, nw, :bus))
    startingbuses = Set(i for (l,i,j)  in ref(pm, nw, :arcs_branch_from))
    leafnodes = setdiff(allbuses, startingbuses)
    for i in leafnodes
        constraint_mc_voltage_psd(pm, nw, i)
    end
end


""
function variable_mc_bus_voltage_prod_hermitian(pm::AbstractUBFModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    bus_ids = collect(ids(pm, nw, :bus))
    terminals = Dict{Int,Vector{Int}}(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))

    if bounded
        # get bounds
        vmax = Dict([(id, ref(pm, nw, :bus, id, "vmax")) for id in bus_ids])
        vmin = Dict([(id, ref(pm, nw, :bus, id, "vmin")) for id in bus_ids])
        # create bounded Hermitian matrix variables
        (Wr,Wi) = variable_mx_hermitian(pm.model, bus_ids, terminals; sqrt_upper_bound=vmax, sqrt_lower_bound=vmin, name="W", prefix="$nw")
    else
        # create unbounded Hermitian matrix variables
        (Wr,Wi) = variable_mx_hermitian(pm.model, bus_ids, terminals; set_lower_bound_diag_to_zero=true, name="W", prefix="$nw")
    end
    v_start = exp.((im*2*pi/3).*[0; -1; 1]) #TODO this should be made more generic eventually
    W_start = v_start*v_start'
    for (id,_) in Wr
        for (i,t) in enumerate(terminals[id])
            for (j,u) in enumerate(terminals[id][1:i])
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :Wr, ids(pm, nw, :bus), Wr)
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :Wi, ids(pm, nw, :bus), Wi)
end


""
function constraint_mc_branch_current_series_product_hermitian(pm::AbstractUBFModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    branches = ref(pm, nw, :branch)
    buses = ref(pm, nw, :bus)

    branch_ids = collect(keys(branches))
    connections = Dict{Int,Vector{Int}}(l => br["f_connections"] for (l,br) in ref(pm, nw, :branch))

    if bounded
        # calculate max series current for each branch
        cmax = Dict{eltype(branch_ids), Vector{Real}}()
        for (key, branch) in branches
            bus_fr = buses[branch["f_bus"]]
            bus_to = buses[branch["t_bus"]]
            cmax[key] = _calc_branch_series_current_max(branch, bus_fr, bus_to)
        end
        # create matrix variables
        (Lr,Li) = variable_mx_hermitian(pm.model, branch_ids, connections; sqrt_upper_bound=cmax, set_lower_bound_diag_to_zero=true, name="CC", prefix="$nw")
    else
        (Lr,Li) = variable_mx_hermitian(pm.model, branch_ids, connections; set_lower_bound_diag_to_zero=true, name="CC", prefix="$nw")
    end

    for (id,L) in Lr
        JuMP.set_start_value.(LinearAlgebra.diag(Lr[id]), 0.01)
    end

    # save reference
    var(pm, nw)[:CCr] = Lr
    var(pm, nw)[:CCi] = Li
    var(pm, nw)[:cm] = Dict([(id, diag(Lr[id])) for id in branch_ids])

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch, :CCr, ids(pm, nw, :branch), Lr)
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch, :CCi, ids(pm, nw, :branch), Li)
end


""
function variable_mc_branch_power(pm::AbstractUBFModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # calculate S bound
    branch_arcs = ref(pm, nw, :arcs_branch)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)

    if bounded
        bound = Dict{eltype(branch_arcs), Matrix{Real}}()
        for (br, branch) in ref(pm, nw, :branch)
            bus_fr = ref(pm, nw, :bus, branch["f_bus"])
            bus_to = ref(pm, nw, :bus, branch["t_bus"])

            smax_fr = _calc_branch_power_max(branch, bus_fr)
            smax_to = _calc_branch_power_max(branch, bus_to)
            cmax_fr, cmax_to = _calc_branch_current_max_frto(branch, bus_fr, bus_to)

            tuple_fr = (br, bus_fr["index"], bus_to["index"])
            tuple_to = (br, bus_to["index"], bus_fr["index"])

            bound[tuple_fr] = bus_fr["vmax"][[findfirst(isequal(c), bus_fr["terminals"]) for c in branch["f_connections"]]].*cmax_fr'
            bound[tuple_to] = bus_to["vmax"][[findfirst(isequal(c), bus_to["terminals"]) for c in branch["t_connections"]]].*cmax_to'

            for (idx, (fc,tc)) in enumerate(zip(branch["f_connections"], branch["t_connections"]))
                bound[tuple_fr][idx,idx] = smax_fr[idx]
                bound[tuple_to][idx,idx] = smax_to[idx]
            end
        end
        # create matrix variables
        (P,Q) = variable_mx_complex(pm.model, branch_arcs, connections, connections; symm_bound=bound, name=("P", "Q"), prefix="$nw")
    else
        (P,Q) = variable_mx_complex(pm.model, branch_arcs, connections, connections; name=("P", "Q"), prefix="$nw")
    end
    # save reference
    var(pm, nw)[:P] = P
    var(pm, nw)[:Q] = Q

    var(pm, nw)[:p] = Dict([(id,diag(P[id])) for id in branch_arcs])
    var(pm, nw)[:q] = Dict([(id,diag(Q[id])) for id in branch_arcs])

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch, :Pf, :Pt, ref(pm, nw, :arcs_branch_from), ref(pm, nw, :arcs_branch_to), P)
    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch, :Qf, :Qt, ref(pm, nw, :arcs_branch_from), ref(pm, nw, :arcs_branch_to), Q)
end


"Defines branch flow model power flow equations"
function constraint_mc_power_losses(pm::AbstractUBFModels, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, r::Matrix{<:Real}, x::Matrix{<:Real}, g_sh_fr::Matrix{<:Real}, g_sh_to::Matrix{<:Real}, b_sh_fr::Matrix{<:Real}, b_sh_to::Matrix{<:Real})
    fr_bus_terminals = ref(pm, nw, :bus, f_bus, "terminals")
    to_bus_terminals = ref(pm, nw, :bus, t_bus, "terminals")

    f_connections = ref(pm, nw, :branch, i, "f_connections")
    t_connections = ref(pm, nw, :branch, i, "t_connections")

    P_to = var(pm, nw, :P)[t_idx]
    Q_to = var(pm, nw, :Q)[t_idx]

    P_fr = var(pm, nw, :P)[f_idx]
    Q_fr = var(pm, nw, :Q)[f_idx]

    Wr_to = var(pm, nw, :Wr)[t_bus][[findfirst(isequal(fc), to_bus_terminals) for fc in f_connections],[findfirst(isequal(tc), to_bus_terminals) for tc in t_connections]]
    Wr_fr = var(pm, nw, :Wr)[f_bus][[findfirst(isequal(fc), fr_bus_terminals) for fc in f_connections],[findfirst(isequal(tc), fr_bus_terminals) for tc in t_connections]]

    Wi_to = var(pm, nw, :Wi)[t_bus][[findfirst(isequal(fc), to_bus_terminals) for fc in f_connections],[findfirst(isequal(tc), to_bus_terminals) for tc in t_connections]]
    Wi_fr = var(pm, nw, :Wi)[f_bus][[findfirst(isequal(fc), fr_bus_terminals) for fc in f_connections],[findfirst(isequal(tc), fr_bus_terminals) for tc in t_connections]]

    CCr =  var(pm, nw, :CCr)[i]
    CCi =  var(pm, nw, :CCi)[i]

    JuMP.@constraint(pm.model, P_fr + P_to .==  Wr_fr*(g_sh_fr)' + Wi_fr*(b_sh_fr)' + r*CCr - x*CCi +  Wr_to*(g_sh_to)'  + Wi_to*(b_sh_to)')
    JuMP.@constraint(pm.model, Q_fr + Q_to .==  Wi_fr*(g_sh_fr)' - Wr_fr*(b_sh_fr)' + x*CCr + r*CCi +  Wi_to*(g_sh_to)'  - Wr_to*(b_sh_to)')
end


""
function constraint_mc_theta_ref(pm::AbstractUBFModels, nw::Int, i::Int, va_ref::Vector{<:Real})
    nconductors = length(va_ref)

    Wr = var(pm, nw, :Wr, i)
    Wi = var(pm, nw, :Wi, i)

    beta = exp.(im.*va_ref)
    gamma = beta*beta'

    Wr_ref = real(gamma).*Wr[1,1]
    Wi_ref = imag(gamma).*Wi[1,1]
    JuMP.@constraint(pm.model, diag(Wr)[2:nconductors]        .== diag(Wr_ref)[2:nconductors]) # first equality is implied
    JuMP.@constraint(pm.model, _mat2utrivec!(Wr) .== _mat2utrivec!(Wr_ref))
    JuMP.@constraint(pm.model, _mat2utrivec!(Wi) .== _mat2utrivec!(Wi_ref))
end


"Defines voltage drop over a branch, linking from and to side voltage"
function constraint_mc_model_voltage_magnitude_difference(pm::AbstractUBFModels, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, r::Matrix{<:Real}, x::Matrix{<:Real}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real})
    fr_bus_terminals = ref(pm, nw, :bus, f_bus, "terminals")
    to_bus_terminals = ref(pm, nw, :bus, t_bus, "terminals")

    f_connections = ref(pm, nw, :branch, i, "f_connections")
    t_connections = ref(pm, nw, :branch, i, "t_connections")

    Wr_fr = var(pm, nw, :Wr)[f_bus][[findfirst(isequal(fc), fr_bus_terminals) for fc in f_connections],[findfirst(isequal(tc), fr_bus_terminals) for tc in t_connections]]
    Wi_fr = var(pm, nw, :Wi)[f_bus][[findfirst(isequal(fc), fr_bus_terminals) for fc in f_connections],[findfirst(isequal(tc), fr_bus_terminals) for tc in t_connections]]

    Wr_to = var(pm, nw, :Wr)[t_bus][[findfirst(isequal(fc), to_bus_terminals) for fc in f_connections],[findfirst(isequal(tc), to_bus_terminals) for tc in t_connections]]
    Wi_to = var(pm, nw, :Wi)[t_bus][[findfirst(isequal(fc), to_bus_terminals) for fc in f_connections],[findfirst(isequal(tc), to_bus_terminals) for tc in t_connections]]

    p_fr = var(pm, nw, :P)[f_idx]
    q_fr = var(pm, nw, :Q)[f_idx]

    p_s_fr = p_fr - (Wr_fr*(g_sh_fr)' + Wi_fr*(b_sh_fr)')
    q_s_fr = q_fr - (Wi_fr*(g_sh_fr)' - Wr_fr*(b_sh_fr)')

    CCr =  var(pm, nw, :CCr)[i]
    CCi =  var(pm, nw, :CCi)[i]

    #KVL over the line:
    JuMP.@constraint(pm.model,          diag(Wr_to) .==          diag(Wr_fr - p_s_fr *r' - q_s_fr*x' - r*p_s_fr' - x*q_s_fr' + r*CCr*r' - x*CCi*r' + x*CCr*x' + r*CCi*x'))
    JuMP.@constraint(pm.model, _mat2utrivec!(Wr_to) .== _mat2utrivec!(Wr_fr - p_s_fr *r' - q_s_fr*x' - r*p_s_fr' - x*q_s_fr' + r*CCr*r' - x*CCi*r' + x*CCr*x' + r*CCi*x'))
    JuMP.@constraint(pm.model, _mat2utrivec!(Wi_to) .== _mat2utrivec!(Wi_fr - q_s_fr *r' + p_s_fr*x' - x*p_s_fr' + r*q_s_fr' + x*CCr*r' + r*CCi*r' - r*CCr*x' + x*CCi*x'))
end


"""
For the matrix KCL formulation, the generator needs an explicit current and
power variable.
"""
function variable_mc_generator_power(pm::SDPUBFKCLMXModel; kwargs...)
    variable_mc_generator_current(pm; kwargs...)
    variable_mc_generator_power_mx(pm; kwargs...)
end


"""
For the matrix KCL formulation, the generator needs an explicit power
variable.
"""
function variable_mc_generator_power_mx(pm::SDPUBFKCLMXModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    @assert(bounded)

    gen_ids = collect(ids(pm, nw, :gen))
    connections = Dict{Int,Vector{Int}}(id => gen["connections"] for (id,gen) in ref(pm, nw, :gen))

    # calculate bounds for matrix variable
    Pg_min = Dict{eltype(gen_ids), Matrix{Real}}()
    Pg_max = Dict{eltype(gen_ids), Matrix{Real}}()
    Qg_min = Dict{eltype(gen_ids), Matrix{Real}}()
    Qg_max = Dict{eltype(gen_ids), Matrix{Real}}()
    for (id, gen) in ref(pm, nw, :gen)
        ncnds = length(connections[id])

        bus = ref(pm, nw, :bus, gen["gen_bus"])
        vmax = haskey(bus, "vmax") ? bus["vmax"][[findfirst(isequal(c), bus["terminals"]) for c in connections[id]]] : fill(Inf, ncnds)
        cmax = _calc_gen_current_max(gen, bus)
        S_bound = _mat_mult_rm_nan(vmax, cmax')

        Pg_min[id] = Qg_min[id] = -S_bound
        Pg_max[id] = Qg_max[id] =  S_bound

        pmin = get(gen, "pmin", fill(-Inf, ncnds))
        pmax = get(gen, "pmax", fill( Inf, ncnds))
        qmin = get(gen, "qmin", fill(-Inf, ncnds))
        qmax = get(gen, "qmax", fill( Inf, ncnds))

        for (idx,c) in enumerate(connections[id])
            Pg_min[id][idx,idx] = max(pmin[idx], Pg_min[id][idx,idx])
            Pg_max[id][idx,idx] = min(pmax[idx], Pg_max[id][idx,idx])
            Qg_min[id][idx,idx] = max(qmin[idx], Qg_min[id][idx,idx])
            Qg_max[id][idx,idx] = min(qmax[idx], Qg_max[id][idx,idx])
        end
    end

    # create matrix variables
    Pg = variable_mx_real(pm.model, gen_ids, connections, connections; lower_bound=Pg_min, upper_bound=Pg_max, name="Pg", prefix="$nw")
    Qg = variable_mx_real(pm.model, gen_ids, connections, connections; lower_bound=Qg_min, upper_bound=Qg_max, name="Qg", prefix="$nw")

    # save references
    var(pm, nw)[:Pg_bus] = Pg
    var(pm, nw)[:Qg_bus] = Qg
    var(pm, nw)[:pg] = Dict{Int, Any}([(id, diag(Pg[id])) for id in gen_ids])
    var(pm, nw)[:qg] = Dict{Int, Any}([(id, diag(Qg[id])) for id in gen_ids])

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :Pg_bus, ids(pm, nw, :gen), Pg)
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :Qg_bus, ids(pm, nw, :gen), Qg)
end


"""
For the matrix KCL formulation, the generator needs an explicit current
variable.
"""
function variable_mc_generator_current(pm::AbstractUBFModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    @assert(bounded)

    gen_ids = collect(ids(pm, nw, :gen))
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

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :CCgr, ids(pm, nw, :gen), CCgr)
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :CCgi, ids(pm, nw, :gen), CCgi)
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
function variable_mc_load_power(pm::AbstractUBFModels; nw=nw_id_default)
    load_wye_ids = [id for (id, load) in ref(pm, nw, :load) if load["configuration"]==WYE]
    load_del_ids = [id for (id, load) in ref(pm, nw, :load) if load["configuration"]==DELTA]
    load_cone_ids = [id for (id, load) in ref(pm, nw, :load) if _check_load_needs_cone(load)]
    # create dictionaries
    var(pm, nw)[:pd_bus] = Dict()
    var(pm, nw)[:qd_bus] = Dict()
    var(pm, nw)[:pd] = Dict()
    var(pm, nw)[:qd] = Dict()
    # now, create auxilary power variable X for delta loads
    variable_mc_load_power_delta_aux(pm, load_del_ids; nw=nw)
    # only delta loads need a current variable
    variable_mc_load_current(pm, load_del_ids; nw=nw)
    # for wye loads with a cone inclusion constraint, we need to create a variable
    variable_mc_load_power(pm, intersect(load_wye_ids, load_cone_ids); nw=nw)

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
function variable_mc_load_power(pm::SDPUBFKCLMXModel; nw::Int=nw_id_default)
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
function variable_mc_load_power(pm::AbstractUBFModels, load_ids::Vector{Int}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
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
    connections = Dict(i => load["connections"] for (i,load) in ref(pm, nw, :load))

    pd = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_pd_$(i)"
        ) for i in load_ids
    )
    qd = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_qd_$(i)"
        ) for i in load_ids
    )

    if bounded
        for i in load_ids
            load = ref(pm, nw, :load, i)
            bus = ref(pm, nw, :bus, load["load_bus"])
            pmin, pmax, qmin, qmax = _calc_load_pq_bounds(load, bus)
            for (idx,c) in enumerate(connections[i])
                set_lower_bound(pd[i][c], pmin[idx])
                set_upper_bound(pd[i][c], pmax[idx])
                set_lower_bound(qd[i][c], qmin[idx])
                set_upper_bound(qd[i][c], qmax[idx])
            end
        end
    end

    #store in dict, but do not overwrite
    for i in load_ids
        var(pm, nw)[:pd][i] = pd[i]
        var(pm, nw)[:qd][i] = qd[i]
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :pd, load_ids, pd)
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :qd, load_ids, qd)
end


"""
The bus qualifier denotes that this is the power withdrawn at the bus; Only for
grounded wye-connected loads, this is the same as the power consumed by the
multi-phase load. The off-diagonals only need to be created for the matrix KCL
formulation.
"""
function variable_mc_load_power_bus(pm::SDPUBFKCLMXModel, load_ids::Vector{Int}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    @assert(bounded)
    connections = Dict{Int,Vector{Int}}(i => load["connections"] for (i,load) in ref(pm, nw, :load))
    # calculate bounds
    bound = Dict{eltype(load_ids), Matrix{Real}}()
    for id in load_ids
        load = ref(pm, nw, :load, id)
        @assert(load["configuration"]==WYE)
        bus = ref(pm, nw, :bus, load["load_bus"])
        cmax = _calc_load_current_max(load, bus)
        bound[id] = bus["vmax"][[findfirst(isequal(c), bus["terminals"]) for c in connections[id]]]*cmax'
    end
    # create matrix variables
    (Pd_bus,Qd_bus) = variable_mx_complex_with_diag(pm.model, load_ids, connections; symm_bound=bound, name=("Pd_bus", "Qd_bus"), prefix="$nw")
    for id in load_ids
        var(pm, nw, :Pd_bus)[id] = Pd_bus[id]
        var(pm, nw, :Qd_bus)[id] = Qd_bus[id]
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :Pd_bus, load_ids, Pd_bus)
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :Qd_bus, load_ids, Qd_bus)
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
function variable_mc_load_power_delta_aux(pm::AbstractUBFModels, load_ids::Vector{Int}; nw::Int=nw_id_default, eps::Real=0.1, bounded::Bool=true, report::Bool=true)
    @assert(bounded)
    connections = Dict{Int,Vector{Int}}(id => load["connections"] for (id,load) in ref(pm, nw, :load))
    # calculate bounds
    bound = Dict{eltype(load_ids), Matrix{Real}}()
    for id in load_ids
        load = ref(pm, nw, :load, id)
        bus_id = load["load_bus"]
        bus = ref(pm, nw, :bus, bus_id)
        cmax = _calc_load_current_max(load, bus)
        bound[id] = bus["vmax"][[findfirst(isequal(c), bus["terminals"]) for c in connections[id]]]*cmax'
    end
    # create matrix variables
    (Xdr,Xdi) = variable_mx_complex(pm.model, load_ids, connections, connections; symm_bound=bound, name="Xd", prefix="$nw")
    # save references
    var(pm, nw)[:Xdr] = Xdr
    var(pm, nw)[:Xdi] = Xdi

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :Xdr, load_ids, Xdr)
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :Xdi, load_ids, Xdi)
end


"""
All loads need a current variable; for wye loads, this variable will be in the
wye reference frame whilst for delta currents it will be in the delta reference
frame.
"""
function variable_mc_load_current(pm::AbstractUBFModels, load_ids::Vector{Int}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    @assert(bounded)

    connections = Dict{Int,Vector{Int}}(i => load["connections"] for (i,load) in ref(pm, nw, :load))

    # calculate bounds
    cmin = Dict{eltype(load_ids), Vector{Real}}()
    cmax = Dict{eltype(load_ids), Vector{Real}}()
    for (id, load) in ref(pm, nw, :load)
        bus_id = load["load_bus"]
        bus = ref(pm, nw, :bus, bus_id)
        cmin[id], cmax[id] = _calc_load_current_magnitude_bounds(load, bus)
    end
    # create matrix variables
    (CCdr, CCdi) = variable_mx_hermitian(pm.model, load_ids, connections; sqrt_upper_bound=cmax, sqrt_lower_bound=cmin, name="CCd", prefix="$nw")
    # save references
    var(pm, nw)[:CCdr] = CCdr
    var(pm, nw)[:CCdi] = CCdi

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :CCdr, load_ids, CCdr)
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :CCdi, load_ids, CCdi)
end


"""
Link the current and power withdrawn by a generator at the bus through a PSD
constraint. The rank-1 constraint is dropped in this formulation.
"""
function constraint_mc_generator_power(pm::SDPUBFKCLMXModel, gen_id::Int; nw::Int=nw_id_default)
    bus_id = ref(pm, nw, :gen, gen_id, "gen_bus")
    connections = ref(pm, nw, :gen, gen_id, "connections")
    terminals = ref(pm, nw, :bus, bus_id, "terminals")
    Pg = var(pm, nw, :Pg_bus, gen_id)
    Qg = var(pm, nw, :Qg_bus, gen_id)
    Wr = var(pm, nw, :Wr, bus_id)[[findfirst(isequal(c), terminals) for c in connections],[findfirst(isequal(c), terminals) for c in connections]]
    Wi = var(pm, nw, :Wi, bus_id)[[findfirst(isequal(c), terminals) for c in connections],[findfirst(isequal(c), terminals) for c in connections]]
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
function constraint_pqw(model::JuMP.Model, w::JuMP.VariableRef, p::JuMP.VariableRef, a::Real, alpha::Real, wmin::Real, wmax::Real, pmin::Real, pmax::Real)
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
function constraint_mc_load_power(pm::AbstractUBFModels, load_id::Int; nw::Int=nw_id_default, report::Bool=true)
    # shared variables and parameters
    load = ref(pm, nw, :load, load_id)
    connections = load["connections"]
    pd0 = load["pd"]
    qd0 = load["qd"]
    bus_id = load["load_bus"]
    bus = ref(pm, nw, :bus, bus_id)
    terminals = bus["terminals"]

    # calculate load params
    a, alpha, b, beta = _load_expmodel_params(load, bus)
    vmin, vmax = _calc_load_vbounds(load, bus)
    wmin = vmin.^2
    wmax = vmax.^2
    pmin, pmax, qmin, qmax = _calc_load_pq_bounds(load, bus)

    # take care of connections
    if load["configuration"]==WYE
        if load["model"]==POWER
            var(pm, nw, :pd)[load_id] = JuMP.Containers.DenseAxisArray(pd0, connections)
            var(pm, nw, :qd)[load_id] = JuMP.Containers.DenseAxisArray(qd0, connections)
        elseif load["model"]==IMPEDANCE
            w = var(pm, nw, :w)[bus_id][[findfirst(isequal(c), terminals) for c in connections]]
            var(pm, nw, :pd)[load_id] = a.*w
            var(pm, nw, :qd)[load_id] = b.*w
        # in this case, :pd has a JuMP variable
        else
            Wr = var(pm, nw, :Wr, bus_id)[[findfirst(isequal(c), terminals) for c in connections],[findfirst(isequal(c), terminals) for c in connections]]
            pd = var(pm, nw, :pd)[load_id]
            qd = var(pm, nw, :qd)[load_id]
            for (idx, c) in enumerate(load["connections"])
                constraint_pqw(pm.model, Wr[idx,idx], pd[idx], a[idx], alpha[idx], wmin[idx], wmax[idx], pmin[idx], pmax[idx])
                constraint_pqw(pm.model, Wr[idx,idx], qd[idx], b[idx], beta[idx], wmin[idx], wmax[idx], qmin[idx], qmax[idx])
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
        Wr = var(pm, nw, :Wr, bus_id)[[findfirst(isequal(c), terminals) for c in connections],[findfirst(isequal(c), terminals) for c in connections]]
        Wi = var(pm, nw, :Wi, bus_id)[[findfirst(isequal(c), terminals) for c in connections],[findfirst(isequal(c), terminals) for c in connections]]
        CCdr = var(pm, nw, :CCdr, load_id)
        CCdi = var(pm, nw, :CCdi, load_id)
        Xdr = var(pm, nw, :Xdr, load_id)
        Xdi = var(pm, nw, :Xdi, load_id)
        Td = [1 -1 0; 0 1 -1; -1 0 1]  # TODO
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
            for (idx, c) in enumerate(connections)
                JuMP.@constraint(pm.model, pd[idx]==pd0[idx])
                JuMP.@constraint(pm.model, qd[idx]==qd0[idx])
            end
        elseif load["model"]==IMPEDANCE
            for (idx,idx) in enumerate(connections)
                JuMP.@constraint(pm.model, pd[idx]==a[idx]*wd[idx])
                JuMP.@constraint(pm.model, qd[idx]==b[idx]*wd[idx])
            end
        else
            for (idx,c) in enumerate(connections)
                constraint_pqw(pm.model, wd[idx], pd[idx], a[idx], alpha[idx], wmin[idx], wmax[idx], pmin[idx], pmax[idx])
                constraint_pqw(pm.model, wd[idx], qd[idx], b[idx], beta[idx],  wmin[idx], wmax[idx], qmin[idx], qmax[idx])
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
function constraint_mc_load_power(pm::SDPUBFKCLMXModel, load_id::Int; nw::Int=nw_id_default, report::Bool=true)
    # shared variables and parameters
    load = ref(pm, nw, :load, load_id)
    connections = load["connections"]
    pd0 = load["pd"]
    qd0 = load["qd"]
    bus_id = load["load_bus"]
    bus = ref(pm, nw, :bus, bus_id)
    terminals = bus["terminals"]

    # calculate load params
    a, alpha, b, beta = _load_expmodel_params(load, bus)
    vmin, vmax = _calc_load_vbounds(load, bus)
    wmin = vmin.^2
    wmax = vmax.^2
    pmin, pmax, qmin, qmax = _calc_load_pq_bounds(load, bus)

    # take care of connections
    Wr = var(pm, nw, :Wr, bus_id)[[findfirst(isequal(c), terminals) for c in connections],[findfirst(isequal(c), terminals) for c in connections]]
    Wi = var(pm, nw, :Wi, bus_id)[[findfirst(isequal(c), terminals) for c in connections],[findfirst(isequal(c), terminals) for c in connections]]
    CCdr = var(pm, nw, :CCdr, load_id)
    CCdi = var(pm, nw, :CCdi, load_id)

    if load["configuration"]==WYE
        if load["model"]==POWER
            var(pm, nw, :pd)[load_id] = JuMP.Containers.DenseAxisArray(pd0, connections)
            var(pm, nw, :qd)[load_id] = JuMP.Containers.DenseAxisArray(qd0, connections)
        elseif load["model"]==IMPEDANCE
            w = var(pm, nw, :w, bus_id)[[findfirst(isequal(c), terminals) for c in connections]]
            var(pm, nw, :pd)[load_id] = a.*w
            var(pm, nw, :qd)[load_id] = b.*w
        else
            pd = var(pm, nw, :pd)[load_id]
            qd = var(pm, nw, :qd)[load_id]
            for (idx, c) in enumerate(connections)
                constraint_pqw(pm.model, Wr[idx,idx], pd[idx], a[idx], alpha[idx], wmin[idx], wmax[idx], pmin[idx], pmax[idx])
                constraint_pqw(pm.model, Wr[idx,idx], qd[idx], b[idx], beta[idx], wmin[idx], wmax[idx], qmin[idx], qmax[idx])
            end
        end
        # diagonal of :Pd is identical to :pd now
        Pd_bus = var(pm, nw, :Pd_bus)[load_id]
        Qd_bus = var(pm, nw, :Qd_bus)[load_id]
        for (idx,c) in enumerate(connections)
            Pd_bus[idx,idx] = var(pm, nw, :pd)[load_id][c]
            Qd_bus[idx,idx] = var(pm, nw, :qd)[load_id][c]
        end

    elseif load["configuration"]==DELTA
        # link Wy, CCd and X
        Xdr = var(pm, nw, :Xdr, load_id)
        Xdi = var(pm, nw, :Xdi, load_id)
        Td = [1 -1 0; 0 1 -1; -1 0 1]  # TODO
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
            for (idx,c) in enumerate(connections)
                JuMP.@constraint(pm.model, pd[idx]==pd0[idx])
                JuMP.@constraint(pm.model, qd[idx]==qd0[idx])
            end
        elseif load["model"]==IMPEDANCE
            for (idx,c) in enumerate(connections)
                JuMP.@constraint(pm.model, pd[idx]==a[idx]*wd[idx])
                JuMP.@constraint(pm.model, qd[idx]==b[idx]*wd[idx])
            end
        else
            for (idx,c) in enumerate(connections)
                constraint_pqw(pm.model, wd[idx], pd[idx], a[idx], alpha[idx], wmin[idx], wmax[idx], pmin[idx], pmax[idx])
                constraint_pqw(pm.model, wd[idx], qd[idx], b[idx], beta[idx], wmin[idx], wmax[idx], qmin[idx], qmax[idx])
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
function constraint_mc_power_balance(pm::KCLMXModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    Wr = var(pm, nw, :Wr, i)
    Wi = var(pm, nw, :Wi, i)

    P = get(var(pm, nw), :P, Dict()); _check_var_keys(P, bus_arcs, "active power", "branch")
    Q = get(var(pm, nw), :Q, Dict()); _check_var_keys(Q, bus_arcs, "reactive power", "branch")
    Pg = get(var(pm, nw), :Pg_bus, Dict()); _check_var_keys(Pg, bus_gens, "active power", "generator")
    Qg = get(var(pm, nw), :Qg_bus, Dict()); _check_var_keys(Qg, bus_gens, "reactive power", "generator")
    Pd = get(var(pm, nw), :Pd_bus, Dict()); _check_var_keys(Pd, bus_loads, "active power", "load")
    Qd = get(var(pm, nw), :Qd_bus, Dict()); _check_var_keys(Qd, bus_loads, "reactive power", "load")

    Gt, Bt = _build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx, t) in ungrounded_terminals
        for (jdx, u) in ungrounded_terminals
            cp = JuMP.@constraint(pm.model,
                  sum(       P[a][findfirst(isequal(t), conns),findfirst(isequal(u), conns)] for (a, conns) in bus_arcs if t in conns && u in conns)
                # + sum(  Psw[a_sw][findfirst(isequal(t), conns),findfirst(isequal(u), conns)] for (a_sw, conns) in bus_arcs_sw if t in conns && u in conns)
                # + sum(Pt[a_trans][findfirst(isequal(t), conns),findfirst(isequal(u), conns)] for (a_trans, conns) in bus_arcs_trans if t in conns && u in conns)
                ==
                  sum(      Pg[g][findfirst(isequal(t), conns),findfirst(isequal(u), conns)] for (g, conns) in bus_gens if t in conns && u in conns)
                - sum(      Pd[d][findfirst(isequal(t), conns),findfirst(isequal(u), conns)] for (d, conns) in bus_loads if t in conns && u in conns)
                - diag(Wr*Gt'+Wi*Bt')[idx]
            )
            push!(cstr_p, cp)

            cq = JuMP.@constraint(pm.model,
                  sum(       Q[a][findfirst(isequal(t), conns),findfirst(isequal(u), conns)] for (a, conns) in bus_arcs if t in conns && u in conns)
                # + sum(  Qsw[a_sw][findfirst(isequal(t), conns),findfirst(isequal(u), conns)] for (a_sw, conns) in bus_arcs_sw if t in conns && u in conns)
                # + sum(Qt[a_trans][findfirst(isequal(t), conns),findfirst(isequal(u), conns)] for (a_trans, conns) in bus_arcs_trans if t in conns && u in conns)
                ==
                  sum(      Qg[g][findfirst(isequal(t), conns),findfirst(isequal(u), conns)] for (g, conns) in bus_gens if t in conns && u in conns)
                - sum(      Qd[d][findfirst(isequal(t), conns),findfirst(isequal(u), conns)] for (d, conns) in bus_loads if t in conns && u in conns)
                - diag(-Wr*Bt'+Wi*Gt')[idx]
            )
            push!(cstr_q, cq)
        end
    end

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = cstr_q

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end
