
"""
calc_start_voltage(
    data_math::Dict{String,Any};
    max_iter=Inf,
    epsilon::Number=1E-3
)::Dict{Tuple{Int,Any},Union{Complex,Missing}}

Calculate no-load starting values for all bus-terminals pairs.
"""
function calc_start_voltage(
    data_math::Dict{String,Any};
    max_iter=Inf,
    epsilon::Number=1E-3,
    explicit_neutral=true
    )::Dict{Tuple{Int,Any},Union{Complex,Missing}}

    if haskey(data_math, "multinetwork")
        @assert !data_math["multinetwork"] "This method should be called on individual networks."
    end
    if haskey(data_math, "data_model")
        @assert data_math["data_model"]==MATHEMATICAL
    end

    node_links = Dict(vcat([[((bus["index"],t), []) for t in bus["terminals"]] for (_, bus) in data_math["bus"]]...))
    for (id, comp) in [data_math["branch"]..., data_math["switch"]...]
        f_bus   = comp["f_bus"]
        f_conns = comp["f_connections"]
        t_bus   = comp["t_bus"]
        t_conns = comp["t_connections"]

        arcs = []
        append!(arcs, [((f_bus,f_conns[k]), (t_bus,t_conns[k]), 1.0) for k in 1:length(f_conns)])
        append!(arcs, [((t_bus,t_conns[k]), (f_bus,f_conns[k]), 1.0) for k in 1:length(f_conns)])
        for (fr, to, scale) in arcs
            push!(node_links[fr], (to, scale))
        end
    end

    # initialize v_start for all bts to missing
    v_start = Dict{Tuple{Int,Any},Union{Complex,Missing}}((bus["index"],t)=>missing for (_, bus) in data_math["bus"] for t in bus["terminals"])

    for (_, bus) in data_math["bus"]
        # set fixed nodes
        if haskey(bus, "vm") && haskey(bus, "va")
            for (i,t) in enumerate(bus["terminals"])
                v_start[(bus["index"],t)] = bus["vm"][i]*exp(im*bus["va"][i])
            end

        end
        # set grounded nodes to zero
        for t in bus["terminals"][bus["grounded"]]
            v_start[(bus["index"],t)] = 0.0+im*0.0
        end
    end

    progress = round(sum(ismissing.(values(v_start)))/length(v_start)*100, digits=2)
    @debug "it. 0:\t$(100-progress)% specified at start"

    # propogate within zones and then between them across transformers as long as progress is made
    # a 'zone' in this context is a collection of buses which are galvanically connected
    # zones are connected through transformers
    count = 0
    stack = collect(keys(filter(kv->!ismissing(kv[2]), v_start)))
    while !isempty(stack) && count<max_iter
        count += 1

        # propogate all nodes in stack through switches/branches
        while !isempty(stack)
            node = pop!(stack)
            for (to, scale) in node_links[node]
                if ismissing(v_start[to])
                    v_start[to] = v_start[node]*scale
                    push!(stack, to)
                end
            end
        end

        # cross zones through transformers
        for (_, tr) in data_math["transformer"]
            f_bus = tr["f_bus"]; f_conns = tr["f_connections"]
            t_bus = tr["t_bus"]; t_conns = tr["t_connections"]
            tm_scale = calculate_tm_scale(tr, data_math["bus"]["$f_bus"], data_math["bus"]["$t_bus"])
            scale = (tm_scale*tr["polarity"]).*tr["tm_set"]
            if tr["configuration"]==WYE || (tr["configuration"]==DELTA && length(tr["tm_set"])==1)
                v_fr = Array{Union{Complex,Missing}}([v_start[(f_bus, t)] for t in f_conns])
                v_to = Array{Union{Complex,Missing}}([v_start[(t_bus, t)] for t in t_conns])
                # forward propagation
                if all((!).(ismissing.(v_fr))) && any(ismissing.(v_to))
                    if explicit_neutral
                        N = length(v_fr)-1
                        Mpn = [LinearAlgebra.diagm(0=>ones(N)) fill(-1.0, N, 1)]
                    else
                        N = length(v_fr)
                        Mpn = LinearAlgebra.diagm(0=>ones(N)) 
                    end
                    v_fr_pn = Mpn*v_fr
                    if all(ismissing.(v_to))
                        anchor_ind = length(v_to)
                        anchor_val = 0.0*im
                    else
                        (anchor_ind, anchor_val) = [(i, v) for (i, v) in enumerate(v_to) if !ismissing(v)][1]
                    end
                    if explicit_neutral
                        v_to_prop = inv([Mpn; [i==anchor_ind ? 1 : 0 for i in 1:N+1]'])*[v_fr_pn./scale..., 0]
                    else
                        v_to_prop = [v_fr_pn./scale...]
                    end

                    for (i,t) in enumerate(t_conns)
                        if ismissing(v_start[(t_bus,t)])
                            v_start[(t_bus,t)] = v_to_prop[i]
                            push!(stack, (t_bus,t))
                        end
                    end
                end
                # backward propagation
                if all((!).(ismissing.(v_to))) && any(ismissing.(v_fr))
                    if explicit_neutral
                        N = length(v_fr)-1
                        Mpn = [LinearAlgebra.diagm(0=>ones(N)) fill(-1.0, N, 1)]
                    else
                        N = length(v_fr)
                        Mpn = LinearAlgebra.diagm(0=>ones(N)) 
                    end
                    v_to_pn = Mpn*v_to
                    if all(ismissing.(v_fr))
                        anchor_ind = length(v_fr)
                        anchor_val = 0.0*im
                    else
                        (anchor_ind, anchor_val) = [(i, v) for (i, v) in enumerate(v_fr) if !ismissing(v)][1]
                    end
                    if explicit_neutral
                        v_fr_prop = inv([Mpn; [i==anchor_ind ? 1 : 0 for i in 1:N+1]'])*[scale.*v_to_pn..., 0]
                    else
                        v_fr_prop = [scale.*v_to_pn...]
                    end

                    for (i,t) in enumerate(f_conns)
                        if ismissing(v_start[(f_bus,t)])
                            v_start[(f_bus,t)] = v_fr_prop[i]
                            push!(stack, (f_bus,t))
                        end
                    end
                end
            elseif tr["configuration"]==DELTA
                v_fr = Array{Union{Complex,Missing}}([v_start[(f_bus, t)] for t in f_conns])
                v_to = Array{Union{Complex,Missing}}([v_start[(t_bus, t)] for t in t_conns])
                # forward propagation
                if all((!).(ismissing.(v_fr))) && any(ismissing.(v_to))
                    if ismissing(v_to[end])
                        v_to[end] = 0.0*im
                    end
                    v_to_n = v_to[end]
                    M = _get_delta_transformation_matrix(length(v_fr))
                    v_to = [((M*v_fr)./scale.+v_to_n)..., v_to_n]
                    for (i,t) in enumerate(t_conns)
                        if ismissing(v_start[(t_bus,t)])
                            v_start[(t_bus,t)] = v_to[i]
                            push!(stack, (t_bus,t))
                        end
                    end
                end
                # backward propagation
                if all((!).(ismissing.(v_to))) && any(ismissing.(v_fr))
                    v_to = [v_to..., 0]
                    v_to_p = v_to[1:end-1]; v_to_n = v_to[end]
                    M = _get_delta_transformation_matrix(length(v_fr))
                    Mp = [M[1:end-1,:]; fill(1.0, 1, length(v_fr))]
                    v_to_pn_scaled = (v_to_p.-v_to_n).*scale
                    v_fr = inv(Mp)*[v_to_pn_scaled[1:end-1]..., 0.0]
                    for (i,t) in enumerate(f_conns)
                        if ismissing(v_start[(f_bus,t)])
                            v_start[(f_bus,t)] = v_fr[i]
                            push!(stack, (f_bus,t))
                        end
                    end
                end
            end
        end
        progress = round(sum(ismissing.(values(v_start)))/length(v_start)*100, digits=2)
        @debug "it. $count:\t$progress% left to initialize at end"
    end

    # increment non-grounded zero values with epsilon
    for (k,v) in v_start
        (b,t) = k
        data_bus = data_math["bus"][string(b)]
        grounded_t = Dict(data_bus["terminals"].=>data_bus["grounded"])
        if !ismissing(v) && iszero(v) && !grounded_t[t]
            v_start[k] = convert(Complex, epsilon)
        end
    end

    return v_start
end


"""
    add_start_voltage!(
        data_math::Dict{String,Any};
        coordinates=:rectangular,
        uniform_v_start=missing,
        vr_default=0.0,
        vi_default=0.0,
        vm_default=0.0,
        va_default=0.0,
        epsilon::Number=1E-3,
    )::Dict{String,Any}

Adds start values for the voltage to the buses.
For a multinetwork data model, you can calculate the start voltages for a representative network through 'calc_start_voltage',
and pass the result as 'uniform_v_start' to use the same values for all networks and avoid recalculating it for each network.
The argument 'epsilon' controls the offset added to ungrounded terminals which would otherwise be set to zero.
"""
function add_start_voltage!(
    data_math::Dict{String,Any};
    coordinates=:rectangular,
    uniform_v_start=missing,
    vr_default=0.0,
    vi_default=0.0,
    vm_default=0.0,
    va_default=0.0,
    epsilon::Number=1E-3,
    explicit_neutral=true
    )::Dict{String,Any}

    @assert data_math["data_model"]==MATHEMATICAL
    @assert coordinates in [:polar, :rectangular] "Legal values for the 'coordinates' argument are [:polar,:rectangular], not :$coordinates."

    is_mn = haskey(data_math, "multinetwork")

    for (nw,dm) in (is_mn ? data_math["nw"] : [("", data_math)])
        if ismissing(uniform_v_start)
            v_start = calc_start_voltage(dm, epsilon=epsilon, explicit_neutral=explicit_neutral)
        else
            v_start = uniform_v_start
        end
        for (_, bus) in dm["bus"]
            index = bus["index"]
            if coordinates==:rectangular
                bus["vr_start"] = [ismissing(v_start[(index, t)]) ? vr_default : real(v_start[(index, t)]) for t in bus["terminals"]]
                bus["vi_start"] = [ismissing(v_start[(index, t)]) ? vi_default : imag(v_start[(index, t)]) for t in bus["terminals"]]
            elseif coordinates==:polar
                bus["vm_start"] = [ismissing(v_start[(index, t)]) ? vm_default : abs(v_start[(index, t)]) for t in bus["terminals"]]
                bus["va_start"] = [ismissing(v_start[(index, t)]) ? va_default : angle(v_start[(index, t)]) for t in bus["terminals"]]
            end
        end
    end

    return data_math
end


"""
add_start_vrvi!(data_math::Dict{String,Any}; kwargs...)

Short-hand for [`add_start_voltage`](@ref add_start_voltage) with rectangular coordinates (coordinates=:rectangular).
"""
add_start_vrvi!(data_math::Dict{String,Any}; kwargs...) = add_start_voltage!(data_math, coordinates=:rectangular, kwargs...)


"""
add_start_vmva!(data_math::Dict{String,Any}; kwargs...)

Short-hand for [`add_start_voltage`](@ref add_start_voltage) with polar coordinates (coordinates=:polar).
"""
add_start_vmva!(data_math::Dict{String,Any}; kwargs...) = add_start_voltage!(data_math, coordinates=:polar, kwargs...)
