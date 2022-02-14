"Power Flow Problem"
function solve_mc_pf(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mc_pf; kwargs...)
end


"Constructor for Power Flow Problem"
function build_mc_pf(pm::AbstractUnbalancedPowerModel)
    variable_mc_bus_voltage(pm; bounded=false)
    variable_mc_branch_power(pm; bounded=false)
    variable_mc_switch_power(pm; bounded=false)
    variable_mc_transformer_power(pm; bounded=false)
    variable_mc_generator_power(pm; bounded=false)
    variable_mc_load_power(pm; bounded=false)
    variable_mc_storage_power(pm; bounded=false)

    constraint_mc_model_voltage(pm)

    for (i,bus) in ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3

        constraint_mc_theta_ref(pm, i)
        constraint_mc_voltage_magnitude_only(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for (i,bus) in ref(pm, :bus)
        constraint_mc_power_balance(pm, i)

        # PV Bus Constraints
        if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            constraint_mc_voltage_magnitude_only(pm, i)
            for j in ref(pm, :bus_gens, i)
                constraint_mc_gen_power_setpoint_real(pm, j)
            end
        end
    end

    for i in ids(pm, :storage)
        constraint_storage_state(pm, i)
        constraint_storage_complementarity_nl(pm, i)
        constraint_mc_storage_losses(pm, i)
        constraint_mc_storage_thermal_limit(pm, i)
        constraint_mc_storage_power_setpoint_real(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end
end


"Constructor for Power Flow in current-voltage variable space"
function build_mc_pf(pm::AbstractUnbalancedIVRModel)
    # Variables
    variable_mc_bus_voltage(pm, bounded = false)
    variable_mc_branch_current(pm, bounded = false)
    variable_mc_switch_current(pm, bounded=false)
    variable_mc_transformer_current(pm, bounded = false)
    variable_mc_generator_current(pm, bounded = false)
    variable_mc_load_current(pm, bounded = false)

    # Constraints
    for (i,bus) in ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3
        constraint_mc_theta_ref(pm, i)
        constraint_mc_voltage_magnitude_only(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for (i,bus) in ref(pm, :bus)
        constraint_mc_current_balance(pm, i)

        # PV Bus Constraints
        if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2
            constraint_mc_voltage_magnitude_only(pm, i)
            for j in ref(pm, :bus_gens, i)
                constraint_mc_gen_power_setpoint_real(pm, j)
            end
        end
    end

    for i in ids(pm, :branch)
        constraint_mc_current_from(pm, i)
        constraint_mc_current_to(pm, i)

        constraint_mc_bus_voltage_drop(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end
end


"Constructor for Branch Flow Power Flow"
function build_mc_pf(pm::AbstractUBFModels)
    # Variables
    variable_mc_bus_voltage(pm; bounded=false)
    variable_mc_branch_current(pm)
    variable_mc_branch_power(pm)
    variable_mc_switch_power(pm)
    variable_mc_transformer_power(pm; bounded=false)
    variable_mc_generator_power(pm; bounded=false)
    variable_mc_load_power(pm)
    variable_mc_storage_power(pm; bounded=false)

    # Constraints
    constraint_mc_model_current(pm)

    for (i,bus) in ref(pm, :ref_buses)
        if !(typeof(pm)<:LPUBFDiagPowerModel)
            constraint_mc_theta_ref(pm, i)
        end

        @assert bus["bus_type"] == 3
        constraint_mc_voltage_magnitude_only(pm, i)
    end

    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
    end

    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for (i,bus) in ref(pm, :bus)
        constraint_mc_power_balance(pm, i)

        # PV Bus Constraints
        if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            constraint_mc_voltage_magnitude_only(pm, i)
            for j in ref(pm, :bus_gens, i)
                constraint_mc_gen_power_setpoint_real(pm, j)
            end
        end
    end

    for i in ids(pm, :storage)
        constraint_storage_state(pm, i)
        constraint_storage_complementarity_nl(pm, i)
        constraint_mc_storage_losses(pm, i)
        constraint_mc_storage_thermal_limit(pm, i)
        constraint_mc_storage_power_setpoint_real(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_power_losses(pm, i)
        constraint_mc_model_voltage_magnitude_difference(pm, i)
        constraint_mc_voltage_angle_difference(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end
end


## Native power flow solver
"""
internal data required used solving an ac unbalanced power flow
the primary use of this data structure is to prevent re-allocation of memory
between successive power flow solves
* `data_math` -- a PowerModelsDistribution mathematical data dictionary
* `ntype` -- node type: @enum NodeType
* `cc_ns_func_pairs` -- a mapping of nodes to compensation current evaluation functions 
* `indexed_nodes` -- all nodes in the network that are not grounded
* `node_to_idx` -- a mapping of all nodes (bus, terminal) to a unique value
* `fnode_to_idx` -- a mapping of nodes with fixed voltage (bus, terminal) to a unique value
* `vnode_to_idx` -- a mapping of nodes with variable voltage (bus, terminal) to a unique value
* `fixed_nodes` -- all indexed nodes with a known voltage reference
* `Uf` -- vector of known voltages
* `Yf` --  admittance matrix partition: Y[node_other_idx, node_fixed_idx]
* `Yv_LU` -- LU factorization of admittance matrix partition: Y[node_other_idx, node_other_idx]
The postfix `_idx` indicates the admittance matrix indexing convention.
"""
mutable struct PowerFlowData
    data_math::Dict
    ntype::Dict
    cc_ns_func_pairs::Vector{Tuple}
    indexed_nodes::Vector
    node_to_idx::Dict
    fnode_to_idx::Dict
    vnode_to_idx::Dict
    fixed_nodes::Vector
    Uf::Vector
    Yf::Matrix
    Yv_LU::Any
end

"""
    PowerFlowData(
      data_math::Dict,
      v_start::Dict
    )

Constructor for PowerFlowData struct that requires mathematical model and v_start.

v_start assigns the initialisation voltages to appropriate bus terminals.
"""
function PowerFlowData(data_math::Dict{String,<:Any}, v_start::Dict{<:Any,<:Any})
    @assert ismath(data_math) "The model is not a mathematical model "
    ntype = Dict{Any, NodeType}()
    for (_, bus) in data_math["bus"]
        id = bus["index"]
        for (i,t) in enumerate(bus["terminals"])
            if bus["grounded"][i]
                ntype[(id,t)] = GROUNDED
            else
                if haskey(bus, "vm")
                    ntype[(id,t)] = FIXED
                else
                    ntype[(id,t)] = VARIABLE
                end
            end
        end
    end

    cc_ns_func_pairs = []
    ns_yprim = []
    virtual_count = 0
    comp_status_prop = Dict(x=>"status" for x in ["load", "shunt", "storage", "transformer", "switch"])
    comp_status_prop["branch"] = "br_status"
    comp_status_prop["gen"]    = "gen_status"
    for (comp_type, comp_interface) in _CPF_COMPONENT_INTERFACES
        for (id, comp) in data_math[comp_type]
            if comp[comp_status_prop[comp_type]]==1
                (bts, nr_vns, y_prim, cc, ccd) = comp_interface(comp, v_start)

                ungr_btidx  = [ntype[bt]!=GROUNDED for bt in bts]
                ungr_filter = [ungr_btidx..., fill(true, nr_vns)...]

                vns = collect(virtual_count+1:virtual_count+nr_vns)
                ns = [bts..., vns...]

                push!(ns_yprim, (ns[ungr_filter], y_prim[ungr_filter, ungr_filter]))
                if cc != nothing && ccd != nothing
                    push!(cc_ns_func_pairs, (ns, cc, ccd, comp_type, id))
                end

                virtual_count += nr_vns
            end
        end
    end

    for i in 1:virtual_count
        ntype[i] = VIRTUAL
    end

    indexed_nodes = [keys(filter(x->x.second!=GROUNDED,  ntype))...]
    node_to_idx = Dict(bt=>i for (i, bt) in enumerate(indexed_nodes))

    N = length(indexed_nodes)
    Y = spzeros(Complex{Float64}, N, N)
    for (ns, y_prim) in ns_yprim
        inds = [node_to_idx[n] for n in ns]
        Y[inds,inds] .+= y_prim
    end

    node_fixed_idx = [i for i in 1:N if ntype[indexed_nodes[i]]==FIXED]
    node_other_idx = setdiff(1:N, node_fixed_idx)

    fnode_to_idx = Dict(node=>idx for (idx,node) in enumerate(indexed_nodes[node_fixed_idx]))
    vnode_to_idx = Dict(node=>idx for (idx,node) in enumerate(indexed_nodes[node_other_idx]))

    Yf = Y[node_other_idx, node_fixed_idx]
    Yv = Y[node_other_idx, node_other_idx]

    fixed_nodes = [indexed_nodes[i] for i in 1:N if ntype[indexed_nodes[i]]==FIXED]
    Uf = fill(NaN+im*NaN, length(fixed_nodes))

    for (i,(b,t)) in enumerate(fixed_nodes)
        vm_t = Dict(data_math["bus"]["$b"]["terminals"].=>data_math["bus"]["$b"]["vm"])
        va_t = Dict(data_math["bus"]["$b"]["terminals"].=>data_math["bus"]["$b"]["va"])
        Uf[i] = vm_t[t]*exp(im*va_t[t])
    end

    Yv_LU = factorize(Yv)

    return PowerFlowData(data_math, ntype, cc_ns_func_pairs, indexed_nodes, node_to_idx, fnode_to_idx, vnode_to_idx, fixed_nodes, Uf, Yf, Yv_LU)
end


"""
    _get_v(
      pfd::struct,
      Vp::Vector,
      n::Vector
    )

Calculates the voltage from PowerFlowData struct.
"""
function _get_v(pfd::PowerFlowData, Vp::Vector{Complex{Float64}}, n::Tuple{Int64, Int64})
    if pfd.ntype[n] == GROUNDED
        return 0.0+im*0.0
    elseif pfd.ntype[n] == FIXED
        return pfd.Uf[pfd.fnode_to_idx[n]]
    else
        return Vp[pfd.vnode_to_idx[n]]
    end
end


"""
    compute_pf(
      data_math::Dict,
      v_start::Dict,
      max_iter::Int,
      stat_tol::Float,
      verbose::Bool
    )

Computes native power flow and outputs the result dict.
"""
function compute_pf(data_math::Dict{String, Any}; v_start::Union{Dict{<:Any,<:Any},Missing}=missing, max_iter::Int=100, stat_tol::Float64=1E-8, verbose::Bool=false)
    if !ismultinetwork(data_math)
        nw_dm = Dict("0"=>data_math)
    else
        nw_dm = data_math["nw"]
    end

    sol = Dict{String, Any}()
    sol["nw"] = Dict{String, Any}()
    time_build = Dict{String, Any}()
    time_solve = Dict{String, Any}()
    time_post = Dict{String, Any}()
    status = Dict{String, Any}()
    its = Dict{String, Any}()
    stat = Dict{String, Any}()

    for (nw, dm) in nw_dm
        if ismissing(v_start)
            # if epsilon is not zero, voltage initialisation is not correct which then leads to unnecessary iterations
            add_start_voltage!(dm, coordinates=:rectangular, epsilon=0)
            v_start = _bts_to_start_voltage(dm)
        end

        time_build[nw] = @elapsed pfd = PowerFlowData(dm, v_start)

        time_solve[nw] = @elapsed (Uv, status[nw], its[nw], stat[nw]) = _compute_Uv(pfd, max_iter=max_iter, stat_tol=stat_tol)

        time_post[nw] = @elapsed sol["nw"][nw] = build_solution(pfd, Uv)
    end

    res = Dict{String, Any}()
    if !ismultinetwork(data_math)
        res["time_build"] = time_build["0"]
        res["time_solve"] = time_solve["0"]
        res["time_post"] = time_post["0"]
        res["solution"] = sol["nw"]["0"]
        res["termination_status"] = status["0"]
        res["iterations"] = its["0"]
        res["stationarity"] = stat["0"]
        res["time_total"] = time_build["0"]+time_solve["0"]+time_post["0"]
    else
        res["time_build"] = time_build
        res["time_solve"] = time_solve
        res["time_post"] = time_post
        res["solution"] = sol
        res["termination_status"] = status
        res["iterations"] = its
        res["stationarity"] = stat
        res["time_total"] = sum(values(time_build))+sum(values(time_solve))+sum(values(time_post))
        res["all_converged"] = all(x==CONVERGED for (_, x) in status)
    end

    return res
end


"""
    compute_pf(
      pdf::PowerFlowData,
      max_iter::Int,
      stat_tol::Float,
      verbose::Bool
    )

Computes native power flow and requires PowerFlowData.
"""
function compute_pf(pfd::PowerFlowData; max_iter::Int=100, stat_tol::Float64=1E-8, verbose::Bool=false)
    time = @elapsed (Uv, status, its, stat) = _compute_Uv(pfd, max_iter=max_iter, stat_tol=stat_tol)
    return build_result(pfd, Uv, status, its, time, stat)
end


"""
    _compute_Uv(
      pfd::PowerFlowData,
      max_iter::Int,
      stat_tol::Float,
      verbose::Bool
    )
Computes a nonlinear AC power flow in rectangular coordinates based on the admittance
matrix of the network data using the fixed-point current injection method.  
Returns a solution data structure in PowerModelsDistribution Dict format.
"""
function _compute_Uv(pfd::PowerFlowData; max_iter::Int=100, stat_tol::Float64=1E-8, verbose::Bool=false)

    Nv = length(pfd.indexed_nodes)-length(pfd.fixed_nodes)

    Uv0 = pfd.Yv_LU\(-pfd.Yf*pfd.Uf)
    Uv = Uv0

    for it in 1:max_iter
        Iv = zeros(Complex{Float64}, Nv)
        for (ns, cc_func, ccd_func, comp_type, id) in pfd.cc_ns_func_pairs
            Un = [_get_v(pfd, Uv, n) for n in ns]
            Icc = cc_func(Un)
            for (i,n) in enumerate(ns)
                if pfd.ntype[n]==VARIABLE
                    Iv[pfd.vnode_to_idx[n]] += Icc[i]
                end
            end
        end

        Uv_next = pfd.Yv_LU\(Iv.-pfd.Yf*pfd.Uf)

        change = maximum(abs.(Uv .- Uv_next))
        if change <= stat_tol
            return (Uv_next, CONVERGED, it, change)
        elseif it==max_iter
            return (Uv_next, ITERATION_LIMIT, it, change)
        end

        Uv = Uv_next
    end
end


"""
    build_result(
      pfd::PowerFlowData,
      Uv::Vector,
      status::PFTerminationStatus,
      its::Int,
      time::Real,
      stationarity::Real,
      verbose::Bool
    )

Builds the result dict from the solution dict.
"""
function build_result(pfd::PowerFlowData, Uv::Vector{Complex{Float64}}, status::PFTerminationStatus, its::Int, time::Real, stationarity::Real; verbose::Bool=false)
    result = Dict{String, Any}()
    result["termination_status"] = status
    result["solution"] = build_solution(pfd, Uv)
    result["solve_time"] = time
    result["iterations"] = its
    result["stationarity"] = stationarity
    return result
end


"""
    build_solution(
      pfd::PowerFlowData,
      Uv::Vector
    )

Builds the solution dict.
"""
function build_solution(pfd::PowerFlowData, Uv::Vector{Complex{Float64}})
    solution = Dict{String, Any}("bus"=>Dict{String, Any}(), "load"=>Dict{String, Any}(), "gen"=>Dict{String, Any}())
    for (id, bus) in pfd.data_math["bus"]
        ind = bus["index"]
        solution["bus"][id] = Dict{String, Any}()
        v = Dict(t=>_get_v(pfd, Uv, (ind,t)) for t in bus["terminals"])
        solution["bus"][id]["vm"] = Dict("$t"=>abs.(v[t]) for t in bus["terminals"])
        solution["bus"][id]["va"] = Dict("$t"=>angle.(v[t]) for t in bus["terminals"])
    end
    for (ns, cc_func, ccd_func, comp_type, id) in pfd.cc_ns_func_pairs
        v_bt = [solution["bus"]["$busid"]["vm"]["$t"] for (busid, t) in sort(ns)]
        if comp_type == "load"
            solution["load"][id] = ccd_func(v_bt)
        elseif comp_type == "gen"
            solution["gen"][id] = ccd_func(v_bt)
        end
    end
    solution["settings"] = deepcopy(pfd.data_math["settings"])
    solution["per_unit"] = deepcopy(pfd.data_math["per_unit"])

    return solution
end


# COMPONENT INTERFACES
"""
    _cpf_branch_interface(
      branch::Dict,
      v_start::Dict
    )

Branch component interface outputs branch primitive Y matrix.
"""
function _cpf_branch_interface(branch::Dict{String,<:Any}, v_start::Dict{<:Any,<:Any})
    br_size = size(branch["br_r"],1)
    @assert br_size <= 4  "Line impedance matrix should be up to 4x4, but is $(br_size)x$(br_size)"
    Ys = inv(branch["br_r"]+im*branch["br_x"])
    Yfr = branch["g_fr"]+im*branch["b_fr"]
    Yto = branch["g_to"]+im*branch["b_to"]
    y_prim = [(Ys.+Yfr) -Ys; -Ys (Ys.+Yto)]
    bts = [[(branch["f_bus"], t) for t in branch["f_connections"]]..., [(branch["t_bus"], t) for t in branch["t_connections"]]...]
    return  bts, 0, y_prim, nothing, nothing
end


"""
    _cpf_shunt_interface(
      shunt::Dict,
      v_start::Dict
    )

Shunt component interface outputs shunt primitive Y matrix.
"""
function _cpf_shunt_interface(shunt::Dict{String,<:Any}, v_start::Dict{<:Any,<:Any})
    Y = shunt["gs"]+im*shunt["bs"]
    bts = [(shunt["shunt_bus"], t) for t in shunt["connections"]]
    return  bts, 0, Y, nothing, nothing
end


"""
    _cpf_transformer_interface(
      tr::Dict,
      v_start::Dict
    )

Transformer component interface outputs transformer primitive Y matrix.
"""
function _cpf_transformer_interface(tr::Dict{String,<:Any}, v_start::Dict{<:Any,<:Any})
    f_ns = [(tr["f_bus"], t) for t in tr["f_connections"]]
    t_ns = [(tr["t_bus"], t) for t in tr["t_connections"]]
    ts = tr["tm_set"]*tr["tm_nom"]*tr["polarity"]
    if tr["configuration"]==WYE
        npairs_fr = [(f_ns[i], f_ns[end]) for i in 1:length(f_ns)-1]
        npairs_to = [(t_ns[i], t_ns[end]) for i in 1:length(t_ns)-1]
        bts, nr_vns, Y = _compose_yprim_banked_ideal_transformers(ts, npairs_fr, npairs_to)
    else
        @assert length(f_ns)==3
        npairs_fr = [(f_ns[1], f_ns[2]), (f_ns[2], f_ns[3]), (f_ns[3], f_ns[1])]
        npairs_to = [(t_ns[i], t_ns[end]) for i in 1:length(t_ns)-1]
        bts, nr_vns, Y = _compose_yprim_banked_ideal_transformers(ts, npairs_fr, npairs_to)
    end
    return  bts, nr_vns, Y, nothing, nothing
end


"""
    _compose_yprim_banked_ideal_transformers(
      ts::Vector,
      npairs_fr::Tuple,
      npairs_to::Tuple,
      ppm::Float
    )

Modifies ideal transformers to avoid singularity error.
"""
function _compose_yprim_banked_ideal_transformers(ts::Vector{Float64}, npairs_fr::Vector{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}}}, npairs_to::Vector{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}}}; ppm=1)
    y_prim_ind = Dict()

    nph = length(ts)

    ns_fr = unique([[a for (a,b) in npairs_fr]..., [b for (a,b) in npairs_fr]...])
    ns_to = unique([[a for (a,b) in npairs_to]..., [b for (a,b) in npairs_to]...])

    ns = [unique([ns_fr..., ns_to...])..., 1:nph...]
    n_to_idx = Dict(n=>idx for (idx,n) in enumerate(ns))
    Y = spzeros(Float64, length(ns), length(ns))

    for (i,t) in enumerate(ts)
        ns_i = [npairs_fr[i]..., npairs_to[i]..., i]
        inds = [n_to_idx[n] for n in ns_i]
        Y[inds[1:4],inds[5]] .= [1/t, -1/t, -1, 1]
        Y[inds[5],inds[1:4]] .= [1/t, -1/t, -1, 1]
        for k in 1:5
            Y[k,k] += ppm*1E-6
        end
    end

    return ns[1:end-nph], nph, Y
end


"""
    _cpf_load_interface(
      load::Dict,
      v_start::Dict
    )

Load component interface outputs load primitive Y matrix.
"""
function _cpf_load_interface(load::Dict{String,<:Any}, v_start::Dict{<:Any,<:Any})
    bts = [(load["load_bus"], t) for t in load["connections"]]
    v0_bt = [v_start[bt] for bt in bts]
    wires = length(bts)

    conf = load["configuration"]

    if conf==WYE || length(v0_bt)==2
        vd0 = v0_bt[1:end-1] .- v0_bt[end]
        if load["model"]==IMPEDANCE
            g =  load["pd"]./load["vnom_kv"]^2
            b = -load["qd"]./load["vnom_kv"]^2
            y = g+im*b
            if wires == 3
                y_prim = diagm(y)
            elseif wires == 4
                y_prim = [diagm(y) -y; -transpose(y) sum(y)]
            end
            cc_func = nothing
            ccd_func = function(v_bt)
                return y_prim * v_bt
            end
        else
            sd0 = load["pd"]+im*load["qd"]
            c0 = conj.(sd0./vd0)
            y0 = c0./vd0
            y_prim = [diagm(y0) -y0; -transpose(y0) sum(y0)]
            if load["model"]==POWER
                cc_func = function(v_bt)
                    sd = load["pd"]+im*load["qd"]
                    vd = v_bt[1:end-1].-v_bt[end]
                    cd = conj.(sd./vd)
                    cd_bus = [cd..., -sum(cd)]
                    return -(cd_bus .- y_prim*v_bt)
                end
                ccd_func = function(v_bt)
                    return y_prim * v_bt .+ cc_func(v_bt)
                end
            elseif load["model"]==CURRENT
                cc_func = function(v_bt)
                    vd = v_bt[1:end-1].-v_bt[end]
                    pd = load["pd"].*(abs.(vd)./load["vnom_kv"])
                    qd = load["qd"].*(abs.(vd)./load["vnom_kv"])
                    sd = pd+im*qd
                    cd = conj.(sd./vd)
                    cd_bus = [cd..., -sum(cd)]
                    return -(cd_bus .- y_prim*v_bt)
                end
                ccd_func = function(v_bt)
                    return y_prim * v_bt .+ cc_func(v_bt)
                end
            elseif load["model"]==EXPONENTIAL
                cc_func = function(v_bt)
                    vd = v_bt[1:end-1].-v_bt[end]
                    pd = load["pd"].*(abs.(vd)./load["vnom_kv"]).^load["alpha"]
                    qd = load["qd"].*(abs.(vd)./load["vnom_kv"]).^load["beta"]
                    sd = pd+im*qd
                    cd = conj.(sd./vd)
                    cd_bus = [cd..., -sum(cd)]
                    return -(cd_bus .- y_prim*v_bt)
                end
                ccd_func = function(v_bt)
                    return y_prim * v_bt .+ cc_func(v_bt)
                end
            end
        end
    elseif conf==DELTA
        Md = [1 -1 0; 0 1 -1; -1 0 1]
        vd0 = Md*v0_bt
        if load["model"]==IMPEDANCE
            g =  load["pd"]./load["vnom_kv"]^2
            b = -load["qd"]./load["vnom_kv"]^2
            y = g+im*b
            y_prim = Md'*diagm(0=>y)*Md
            cc_func = nothing
            ccd_func = function(v_bt)
                return y_prim * v_bt
            end
        else
            sd0 = load["pd"]+im*load["qd"]
            c0 = conj.(sd0./vd0)
            y0 = c0./vd0
            y_prim = Md'*diagm(0=>y0)*Md
            if load["model"]==POWER
                cc_func = function(v_bt)
                    sd = load["pd"]+im*load["qd"]
                    vd = Md*v_bt
                    cd = conj.(sd./vd)
                    cd_bus = Md'*cd
                    return -(cd_bus .- y_prim*v_bt)
                end
                ccd_func = function(v_bt)
                    return y_prim * v_bt .+ cc_func(v_bt)
                end
            elseif load["model"]==CURRENT
                cc_func = function(v_bt)
                    vd = Md*v_bt
                    pd = load["pd"].*(abs.(vd)./load["vnom_kv"])
                    qd = load["qd"].*(abs.(vd)./load["vnom_kv"])
                    sd = pd+im*qd
                    cd = conj.(sd./vd)
                    cd_bus = Md'*cd
                    return -(cd_bus .- y_prim*v_bt)
                end
                ccd_func = function(v_bt)
                    return y_prim * v_bt .+ cc_func(v_bt)
                end
            elseif load["model"]==EXPONENTIAL
                cc_func = function(v_bt)
                    vd = Md*v_bt
                    pd = load["pd"].*(abs.(vd)./load["vnom_kv"]).^load["alpha"]
                    qd = load["qd"].*(abs.(vd)./load["vnom_kv"]).^load["beta"]
                    sd = pd+im*qd
                    cd = conj.(sd./vd)
                    cd_bus = Md'*cd
                    return -(cd_bus .- y_prim*v_bt)
                end
                ccd_func = function(v_bt)
                    return y_prim * v_bt .+ cc_func(v_bt)
                end
            end
        end
    end

    return bts, 0, y_prim, cc_func, ccd_func
end


"""
    _cpf_generator_interface(
      gen::Dict,
      v_start::Dict
    )

Generator component interface outputs generator primitive Y matrix.
"""
function _cpf_generator_interface(gen::Dict{String,<:Any}, v_start::Dict{<:Any,<:Any})
    bts = [(gen["gen_bus"], t) for t in gen["connections"]]
    wires = length(bts)
    v0_bt = [v_start[bt] for bt in bts]

    conf = gen["configuration"]

    if conf==WYE || length(v0_bt)==2
        if wires == 3
            vd0 = v0_bt
        elseif wires == 4 || wires == 2
            vd0 = v0_bt[1:end-1] .- v0_bt[end]
        end
        sg0 = gen["pg"]+im*gen["qg"]
        sd0 = -sg0
        c0 = conj.(sd0./vd0)
        y0 = c0./vd0
        if wires == 3
            y_prim = diagm(y0)
        elseif wires == 4 || wires == 2
            y_prim = [diagm(y0) -y0; -transpose(y0) sum(y0)]
        end
        
        cc_func = function(v_bt)
            sg = gen["pg"]+im*gen["qg"]
            sd = -sg

            if wires == 3
                vd = v_bt
                cd = conj.(sd./vd)
                cd_bus = cd
                return -(cd_bus .- y_prim*v_bt)
            elseif wires == 4 || wires == 2
                vd = v_bt[1:end-1].-v_bt[end]
                cd = conj.(sd./vd)
                cd_bus = [cd..., -sum(cd)]
                return -(cd_bus .- y_prim*v_bt)
            end
        end
        ccd_func = function(v_bt)
            return y_prim * v_bt .+ cc_func(v_bt)
        end

    elseif conf==DELTA
        Md = [1 -1 0; 0 1 -1; -1 0 1]
        vd0 = Md*v0_bt
        sg0 = gen["pg"]+im*gen["qg"]
        sd0 = -sg0
        c0 = conj.(sd0./vd0)
        y0 = c0./vd0
        y_prim = Md'*diagm(0=>y0)*Md
        cc_func = function(v_bt)
            sg = gen["pg"]+im*gen["qg"]
            sd = -sg
            vd = Md*v_bt
            cd = conj.(sd./vd)
            cd_bus = Md'*cd
            return -(cd_bus .- y_prim*v_bt)
        end
        ccd_func = function(v_bt)
            return y_prim * v_bt .+ cc_func(v_bt)
        end
    end

    return bts, 0, y_prim, cc_func, ccd_func
end


"""
    _cpf_switch_interface(
      switch::Dict,
      v_start::Dict
    )

Branch component interface outputs branch primitive Y matrix.
"""
function _cpf_switch_interface(switch::Dict{String,<:Any}, v_start::Dict{<:Any,<:Any}; small_impedance=1E-8)
    len = length(switch["f_connections"])
    Ys = LinearAlgebra.I(len) .* 1/small_impedance
    y_prim = [Ys -Ys; -Ys Ys]
    bts = [[(switch["f_bus"], t) for t in switch["f_connections"]]..., [(switch["t_bus"], t) for t in switch["t_connections"]]...]
    return  bts, 0, y_prim, nothing, nothing
end


"""
A mapping of supported component types to their functional interfaces.
"""
const _CPF_COMPONENT_INTERFACES = Dict(
    "load"   => _cpf_load_interface,
    "gen"   => _cpf_generator_interface,
    "branch" => _cpf_branch_interface,
    "transformer" => _cpf_transformer_interface,
    "shunt" => _cpf_shunt_interface,
    "switch" => _cpf_switch_interface,
)


"""
    _bts_to_start_voltage(
        dm::Dict
    )

Assigns the initialisation voltages to appropriate bus terminals.
"""
function _bts_to_start_voltage(dm::Dict{String,<:Any})
    v_start = Dict()
    for (i,bus) in dm["bus"]
        for (t, terminal) in enumerate(bus["terminals"])
            v_start[(bus["index"],terminal)] = bus["vr_start"][t] + im* bus["vi_start"][t]
        end
    end
    return v_start
end