
""
function scale(dict, key, scale)
    if haskey(dict, key)
        dict[key] *= scale
    end
end


""
function _get_next_index(last_index, presets)
    new_index = last_index+1
    while new_index in presets
        new_index += 1
    end
    return new_index
end


""
function solution_identify!(solution, data_model; id_prop="id")
    for comp_type in keys(solution)
        if isa(solution[comp_type], Dict)
            comp_dict = Dict{Any, Any}()
            for (ind, comp) in solution[comp_type]
                id = data_model[comp_type][ind][id_prop]
                comp_dict[id] = comp
            end
            solution[comp_type] = comp_dict
        end
    end

    return solution
end


""
function add_solution!(solution, comp_type, id, data)
    if !haskey(solution, comp_type)
        solution[comp_type] = Dict()
    end

    if !haskey(solution[comp_type], id)
        solution[comp_type][id] = Dict{String, Any}()
    end

    for (key, prop) in data
        solution[comp_type][id][key] = prop
    end
end


""
function delete_solution!(solution, comp_type, id, props)
    if haskey(solution, comp_type)
        if haskey(solution[comp_type], id)
            for prop in props
                delete!(solution[comp_type][id], prop)
            end
        end
    end
end


""
function delete_solution!(solution, comp_type, id)
    if haskey(solution, comp_type)
        if haskey(solution[comp_type], id)
            delete!(solution[comp_type], id)
        end
    end
end


""
function _get_new_ground(terminals)
    if isa(terminals, Vector{Int})
        return maximum(terminals)+1
    else
        nrs = [parse(Int, x[1]) for x in [match(r"n([1-9]{1}[0-9]*)", string(t)) for t in terminals] if x!=nothing]
        new = isempty(nrs) ? 1 : maximum(nrs)+1
        if isa(terminals, Vector{Symbol})
            return Symbol("g$new")
        else
            return "g$new"
        end
    end
end


""
function _get_ground!(bus::Dict{String,<:Any})
    # find perfect groundings (true ground)
    grounded_perfect = []
    for i in 1:length(bus["grounded"])
        if bus["rg"][i]==0 && bus["xg"][i]==0
            push!(grounded_perfect, bus["grounded"][i])
        end
    end

    if !isempty(grounded_perfect)
        return grounded_perfect[1]
    else
        g = _get_new_ground(bus["terminals"])
        push!(bus["terminals"], g)
        push!(bus["rg"], 0.0)
        push!(bus["xg"], 0.0)
        return g
    end
end


"""
Converts a set of short-circuit tests to an equivalent reactance network.
Reference:
R. C. Dugan, “A perspective on transformer modeling for distribution system analysis,”
in 2003 IEEE Power Engineering Society General Meeting (IEEE Cat. No.03CH37491), 2003, vol. 1, pp. 114-119 Vol. 1.
"""
function _sc2br_impedance(Zsc::Dict{Tuple{Int,Int},Complex{Float64}})::Dict{Tuple{Int,Int},Complex}
    N = maximum([maximum(k) for k in keys(Zsc)])
    # check whether no keys are missing
    # Zsc should contain tupples for upper triangle of NxN
    for i in 1:N
        for j in i+1:N
            if !haskey(Zsc, (i,j))
                if haskey(Zsc, (j,i))
                    # Zsc is symmetric; use value of lower triangle if defined
                    Zsc[(i,j)] =  Zsc[(j,i)]
                else
                    Memento.error(_LOGGER, "Short-circuit impedance between winding $i and $j is missing.")
                end
            end
        end
    end

    # if all zero, return all zeros
    if all(values(Zsc).==0.0)
        return Zsc
    end

    # make Zb
    Zb = zeros(Complex{Float64}, N-1,N-1)
    for i in 1:N-1
        Zb[i,i] = Zsc[(1,i+1)]
    end
    for i in 1:N-1
        for j in 1:i-1
            Zb[i,j] = (Zb[i,i]+Zb[j,j]-Zsc[(j+1,i+1)])/2
            Zb[j,i] = Zb[i,j]
        end
    end
    # get Ybus
    Y = LinearAlgebra.pinv(Zb)
    Y = [-Y*ones(N-1) Y]
    Y = [-ones(1,N-1)*Y; Y]
    # extract elements
    Zbr = Dict{Tuple{Int,Int},Complex}()
    for k in keys(Zsc)
        Zbr[k] = (abs(Y[k...])==0) ? Inf : -1/Y[k...]
    end
    return Zbr
end


"loss model builder for transformer decomposition"
function _build_loss_model!(data_math::Dict{String,<:Any}, transformer_name::String, to_map::Vector{String}, r_s::Vector{Float64}, zsc::Dict{Tuple{Int,Int},Complex{Float64}}, ysh::Complex{Float64}; nphases::Int=3)::Vector{Int}
    # precompute the minimal set of buses and lines
    N = length(r_s)
    tr_t_bus = collect(1:N)
    buses = Set(1:2*N)

    zbr = _sc2br_impedance(zsc)

    edges = [[[i,i+N] for i in 1:N]..., [[i+N,j+N] for (i,j) in keys(zbr)]...]
    lines = Dict(enumerate(edges))

    z = Dict(enumerate([r_s..., values(zbr)...]))

    shunts = Dict(2=>ysh)

    # remove Inf lines
    for (l,edge) in lines
        if real(z[l])==Inf || imag(z[l])==Inf
            delete!(lines, l)
            delete!(z, l)
        end
    end

    # merge short circuits
    stack = Set(keys(lines))

    while !isempty(stack)
        l = pop!(stack)
        if z[l] == 0
            (i,j) = lines[l]

            # remove line
            delete!(lines, l)

            # remove  bus j
            delete!(buses, j)

            # update lines
            for (k,(edge)) in lines
                if edge[1] == j
                    edge[1] = i
                end
                if edge[2] == j
                    edge[2] = i
                end
                if edge[1]==edge[2]
                    delete!(lines, k)
                    delete!(stack, k)
                end
            end

            # move shunts
            if haskey(shunts, j)
                if haskey(shunts, i)
                    shunts[i] += shunts[j]
                else
                    shunts[i] = shunts[j]
                end
            end

            # update transformer buses
            for w in 1:N
                if tr_t_bus[w] == j
                    tr_t_bus[w] = i
                end
            end
        end
    end

    bus_ids = Dict{Int,Int}()
    for bus in buses
        bus_obj = Dict{String,Any}(
            "name" => "_virtual_bus.transformer.$(transformer_name)_$(bus)",
            "bus_i" => length(data_math["bus"])+1,
            "vm" => fill(1.0, nphases),
            "va" => fill(0.0, nphases),
            "vmin" => fill(0.0, nphases),
            "vmax" => fill(Inf, nphases),
            "base_kv" => 1.0,
            "bus_type" => 1,
            "status" => 1,
            "index" => length(data_math["bus"])+1,
        )

        data_math["bus"]["$(bus_obj["index"])"] = bus_obj

        bus_ids[bus] = bus_obj["bus_i"]

        push!(to_map, "bus.$(bus_obj["index"])")
    end

    for (l,(i,j)) in lines
        # merge the shunts into the shunts of the pi model of the line
        g_fr = b_fr = g_to = b_to = 0

        if haskey(shunts, i)
            g_fr = real(shunts[i])
            b_fr = imag(shunts[i])
            delete!(shunts, i)
        end

        if haskey(shunts, j)
            g_to = real(shunts[j])
            b_to = imag(shunts[j])
            delete!(shunts, j)
        end

        branch_obj = Dict{String,Any}(
            "name" => "_virtual_branch.transformer.$(transformer_name)_$(l)",
            "source_id" => "_virtual_branch.transformer.$(transformer_name)_$(l)",
            "index" => length(data_math["branch"])+1,
            "br_status"=>1,
            "f_bus"=>bus_ids[i],
            "t_bus"=>bus_ids[j],
            "f_connections"=>collect(1:nphases),
            "t_connections"=>collect(1:nphases),
            "br_r" => diagm(0=>fill(real(z[l]), nphases)),
            "br_x" => diagm(0=>fill(imag(z[l]), nphases)),
            "g_fr" => diagm(0=>fill(g_fr, nphases)),
            "b_fr" => diagm(0=>fill(b_fr, nphases)),
            "g_to" => diagm(0=>fill(g_to, nphases)),
            "b_to" => diagm(0=>fill(b_to, nphases)),
            "angmin" => fill(-60.0, nphases),
            "angmax" => fill( 60.0, nphases),
            "shift" => zeros(nphases),
            "tap" => ones(nphases),
            "switch" => false,
            "transformer" => false,
        )

        data_math["branch"]["$(branch_obj["index"])"] = branch_obj

        push!(to_map, "branch.$(branch_obj["index"])")
    end

    return Vector{Int}([bus_ids[bus] for bus in tr_t_bus])
end


""
function _kron_reduce_branch!(object::Dict{String,<:Any}, Zs_keys::Vector{String}, Ys_keys::Vector{String}, terminals::Vector{Int}, neutral::Int)::Vector{Int}
    Zs = Vector{Matrix}([object[k] for k in Zs_keys])
    Ys = Vector{Matrix}([object[k] for k in Ys_keys])
    Zs_kr, Ys_kr, terminals_kr = _kron_reduce_branch(Zs, Ys, terminals, neutral)

    for (i,k) in enumerate(Zs_keys)
        object[k] = Zs_kr[i]
    end

    for (i,k) in enumerate(Ys_keys)
        object[k] = Ys_kr[i]
    end

    return _get_idxs(terminals, terminals_kr)
end


""
function _kron_reduce_branch(Zs::Vector{Matrix}, Ys::Vector{Matrix}, terminals::Vector{Int}, neutral::Int)::Tuple{Vector{Matrix}, Vector{Matrix}, Vector{Int}}
    Zs_kr = Vector{Matrix}([deepcopy(Z) for Z in Zs])
    Ys_kr = Vector{Matrix}([deepcopy(Y) for Y in Ys])
    terminals_kr = deepcopy(terminals)

    while neutral in terminals_kr
        n = _get_ilocs(terminals_kr, neutral)[1]
        P = setdiff(collect(1:length(terminals_kr)), n)

        if all(size(Z) == (length(terminals_kr), length(terminals_kr)) for Z in Zs_kr)
            Zs_kr = Vector{Matrix}([Z[P,P]-(1/Z[n,n])*Z[P,[n]]*Z[[n],P] for Z in Zs_kr])
            Ys_kr = Vector{Matrix}([Y[P,P] for Y in Ys_kr])
        end

        terminals_kr = terminals_kr[P]
    end

    return Zs_kr, Ys_kr, terminals_kr
end


""
function _pad_properties!(object::Dict{<:Any,<:Any}, properties::Vector{String}, connections::Vector{Int}, phases::Vector{Int}; pad_value::Real=0.0)
    @assert(all(c in phases for c in connections))
    inds = _get_idxs(phases, connections)

    for property in properties
        if haskey(object, property)
            if isa(object[property], Vector)
                tmp = fill(pad_value, length(phases))
                tmp[inds] = object[property]
                object[property] = tmp
            elseif isa(object[property], Matrix)
                tmp = fill(pad_value, length(phases), length(phases))
                tmp[inds, inds] = object[property]
                object[property] = tmp
            end
        end
    end
end


""
function _pad_properties_delta!(object::Dict{<:Any,<:Any}, properties::Vector{String}, connections::Vector{Int}, phases::Vector{Int}; invert::Bool=false)
    @assert(all(c in phases for c in connections))
    @assert(length(connections) in [2, 3], "A delta configuration has to have at least 2 or 3 connections!")
    @assert(length(phases)==3, "Padding only possible to a |phases|==3!")

    for property in properties
        val = object[property]
        val_length = length(connections)==2 ? 1 : length(connections)
        @assert(isa(val, Vector) && length(val)==val_length)

        # build tmp
        tmp = Dict()
        sign = invert ? -1 : 1
        if val_length==1
            tmp[(connections[1], connections[2])] =      val[1]
            tmp[(connections[2], connections[1])] = sign*val[1]
        else
            tmp[(connections[1], connections[2])] =      val[1]
            tmp[(connections[2], connections[3])] =      val[2]
            tmp[(connections[3], connections[1])] =      val[3]
        end
        merge!(tmp, Dict((k[2], k[1])=>sign*v for (k,v) in tmp))
        get_val(x,y) = haskey(tmp, (x,y)) ? tmp[(x,y)] : 0.0

        object[property] = [get_val(phases[1], phases[2]), get_val(phases[2], phases[3]), get_val(phases[3], phases[1])]
    end
end


""
function _apply_filter!(obj, properties, filter)
    for property in properties
        if haskey(obj, property)
            if isa(obj[property], Vector)
                obj[property] = obj[property][filter]
            elseif isa(obj[property], Matrix)
                obj[property] = obj[property][filter, filter]
            else
                Memento.error(_LOGGER, "The property $property is not a Vector or a Matrix!")
            end
        end
    end
end


"""
Given a set of addmittances 'y' connected from the conductors 'f_cnds' to the
conductors 't_cnds', this method will return a list of conductors 'cnd' and a
matrix 'Y', which will satisfy I[cnds] = Y*V[cnds].
"""
function calc_shunt(f_cnds::Vector{Int}, t_cnds::Vector{Int}, y::Vector{T})::Tuple{Vector{Int}, Matrix{T}} where T <: Number
    cnds = unique([f_cnds..., t_cnds...])
    e(f,t) = reshape([c==f ? 1 : c==t ? -1 : 0 for c in cnds], length(cnds), 1)
    Y = sum([e(f_cnds[i], t_cnds[i])*y[i]*e(f_cnds[i], t_cnds[i])' for i in 1:length(y)])
    return (cnds, Y)
end



"""
Given a set of terminals 'cnds' with associated shunt addmittance 'Y', this
method will calculate the reduced addmittance matrix if terminal 'ground' is
grounded.
"""
function _calc_ground_shunt_admittance_matrix(cnds::Vector{Int}, Y::Matrix{T}, ground::Int)::Tuple{Vector{Int}, Matrix{T}} where T <: Number
    # TODO add types
    if ground in cnds
        cndsr = setdiff(cnds, ground)
        cndsr_inds = _get_idxs(cnds, cndsr)
        Yr = Y[cndsr_inds, cndsr_inds]
        return (cndsr, Yr)
    else
        return cnds, Y
    end
end
