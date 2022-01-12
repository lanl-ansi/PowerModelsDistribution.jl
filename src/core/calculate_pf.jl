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
    @assert ismath(data_math) """The model is not a mathematical model """
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
    comp_status_prop = Dict(x=>"status" for x in ["load", "shunt", "storage", "transformer"])
    comp_status_prop["branch"] = "br_status"
    comp_status_prop["gen"]    = "gen_status"
    for (comp_type, comp_interface) in _CPF_COMPONENT_INTERFACES
        for (id, comp) in data_math[comp_type]
            if comp[comp_status_prop[comp_type]]==1
                (bts, nr_vns, y_prim, cc) = comp_interface(comp, v_start)

                ungr_btidx  = [ntype[bt]!=GROUNDED for bt in bts]
                ungr_filter = [ungr_btidx..., fill(true, nr_vns)...]

                vns = collect(virtual_count+1:virtual_count+nr_vns)
                ns = [bts..., vns...]

                push!(ns_yprim, (ns[ungr_filter], y_prim[ungr_filter, ungr_filter]))
                if cc != nothing
                    push!(cc_ns_func_pairs, (ns, cc))
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
function compute_pf(data_math::Dict{String, Any}; v_start::Union{Dict{<:Any,<:Any},Missing}=missing, max_iter::Int=100, stat_tol::Float64=1E-8, verbose::Bool=true)
    if !ismultinetwork(data_math)
        nw_dm = Dict("0"=>data_math)
    else
        nw_dm = data_math["nw"]
    end

    sol = Dict{String, Any}()
    time_build = Dict{String, Any}()
    time_solve = Dict{String, Any}()
    time_post = Dict{String, Any}()
    Yf_matrix = Dict{String, Any}()
    Yv_matrix = Dict{String, Any}()
    status = Dict{String, Any}()
    its = Dict{String, Any}()
    stat = Dict{String, Any}()

    for (nw, dm) in nw_dm
        if ismissing(v_start)
            add_start_vrvi!(dm; epsilon=0)  # if epsilon is not zero, voltage initialisation is not correct which then leads to unnecessary iterations
            v_start = _bts_to_start_voltage(dm)
        end


        time_build[nw] = @elapsed pfd = PowerFlowData(dm, v_start)

        time_solve[nw] = @elapsed (Uv, status[nw], its[nw], stat[nw]) = _compute_Uv(pfd, max_iter=max_iter, stat_tol=stat_tol)

        time_post[nw] = @elapsed sol[nw] = build_solution(pfd, Uv)

        Yf_matrix[nw] = pfd.Yf
        Yv_matrix[nw] = pfd.Yv_LU
    end

    res = Dict{String, Any}()
    if !ismultinetwork(data_math)
        res["time_build"] = time_build["0"]
        res["time_solve"] = time_solve["0"]
        res["time_post"] = time_post["0"]
        res["solution"] = sol["0"]
        res["termination_status"] = status["0"]
        res["iterations"] = its["0"]
        res["stationarity"] = stat["0"]
        res["time_total"] = time_build["0"]+time_solve["0"]+time_post["0"]
        res["Yf"] = Yf_matrix["0"]
        res["Yv"] = Yv_matrix["0"]
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
        res["Yf"] = Yf_matrix
        res["Yv"] = Yv_matrix
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

Constructor for compute_pf that requires PowerFlowData.
"""
function compute_pf(pfd::PowerFlowData; max_iter::Int=100, stat_tol::Float64=1E-8, verbose::Bool=true)
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

Computes variable voltages until convergence or reaching max_iter.
"""
function _compute_Uv(pfd::PowerFlowData; max_iter::Int=100, stat_tol::Float64=1E-8, verbose::Bool=true)

    Nv = length(pfd.indexed_nodes)-length(pfd.fixed_nodes)

    Uv0 = pfd.Yv_LU\(-pfd.Yf*pfd.Uf)
    Uv = Uv0

    for it in 1:max_iter
        Iv = zeros(Complex{Float64}, Nv)
        for (ns, cc_func) in pfd.cc_ns_func_pairs
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
function build_result(pfd::PowerFlowData, Uv::Vector{Complex{Float64}}, status::PFTerminationStatus, its::Int, time::Real, stationarity::Real; verbose::Bool=true)
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
    solution = Dict{String, Any}("bus"=>Dict{String, Any}())
    for (id, bus) in pfd.data_math["bus"]
        ind = bus["index"]
        solution["bus"][id] = Dict{String, Any}()
        v = Dict(t=>_get_v(pfd, Uv, (ind,t)) for t in bus["terminals"])
        solution["bus"][id]["vm"] = Dict("$t"=>abs.(v[t]) for t in bus["terminals"])
        solution["bus"][id]["va"] = Dict("$t"=>angle.(v[t]) for t in bus["terminals"])
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
    Ys = inv(branch["br_r"]+im*branch["br_x"])
    Yfr = branch["g_fr"]+im*branch["b_fr"]
    Yto = branch["g_to"]+im*branch["b_to"]
    y_prim = [(Ys.+Yfr) -Ys; -Ys (Ys.+Yto)]
    bts = [[(branch["f_bus"], t) for t in branch["f_connections"]]..., [(branch["t_bus"], t) for t in branch["t_connections"]]...]
    return  bts, 0, y_prim, nothing
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
    return  bts, 0, Y, nothing
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
    return  bts, nr_vns, Y, nothing
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
            end
        end
    end

    return bts, 0, y_prim, cc_func
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
    end

    return bts, 0, y_prim, cc_func
end


const _CPF_COMPONENT_INTERFACES = Dict(
    "load"   => _cpf_load_interface,
    "gen"   => _cpf_generator_interface,
    "branch" => _cpf_branch_interface,
    "transformer" => _cpf_transformer_interface,
    "shunt" => _cpf_shunt_interface,
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