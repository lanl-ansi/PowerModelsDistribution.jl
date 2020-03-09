
#import PowerModelsDistribution
#get = PowerModelsDistribution.get

function scale(dict, key, scale)
    if haskey(dict, key)
        dict[key] *= scale
    end
end


function add_virtual!(data_model, comp_type, comp)
    if !haskey(data_model, comp_type)
        data_model[comp_type] = Dict{Any, Any}()
    end
    comp_dict = data_model[comp_type]
    virtual_ids = [parse(Int, x[1]) for x in [match(r"_virtual_([1-9]{1}[0-9]*)", id) for id in keys(comp_dict) if isa(id, AbstractString)] if x!=nothing]
    if isempty(virtual_ids)
        id = "_virtual_1"
    else
        id = "_virtual_$(maximum(virtual_ids)+1)"
    end
    comp["id"] = id
    comp_dict[id] = comp
    return comp
end

add_virtual_get_id!(data_model, comp_type, comp) = add_virtual!(data_model, comp_type, comp)["id"]

function delete_component!(data_model, comp_type, comp::Dict)
    delete!(data_model[comp_type], comp["id"])
    if isempty(data_model[comp_type])
        delete!(data_model, comp_type)
    end
end

function delete_component!(data_model, comp_type, id::Any)
    delete!(data_model[comp_type], id)
    if isempty(data_model[comp_type])
        delete!(data_model, comp_type)
    end
end

function add_mappings!(data_model::Dict{String, Any}, mapping_type::String, mappings::Vector)
    if !haskey(data_model, "mappings")
        data_model["mappings"] = []
    end

    append!(data_model["mappings"], [(mapping_type, mapping) for mapping in mappings])
end


function _get_next_index(last_index, presets)
    new_index = last_index+1
    while new_index in presets
        new_index += 1
    end
    return new_index
end


function data_model_index!(data_model; components=["bus", "line", "shunt", "generator", "load", "transformer_2wa"], index_presets=Dict())
    comp_id2ind = Dict()

    # bus should be the first component, because we want to
    for comp_type in components
        comp_dict = Dict{String, Any}()

        if !haskey(index_presets, comp_type)
            for (i,(id,comp)) in enumerate(data_model[comp_type])
                @assert(!haskey(comp, "index"), "$comp_type $id: component already has an index.")
                comp["index"] = i
                comp["id"] = id
                comp_dict["$i"] = comp
            end
        else
            last_index = 0

            for (id, comp) in data_model[comp_type]
                @assert(!haskey(comp, "index"), "$comp_type $id: component already has an index.")
                if haskey(index_presets[comp_type], id)
                    comp["index"] = index_presets[comp_type][id]
                else
                    comp["index"] = _get_next_index(last_index, values(index_presets[comp_type]))
                    last_index = comp["index"]
                end

                comp["id"] = id
                comp_dict["$(comp["index"])"] = comp
            end
        end

        data_model[comp_type] = comp_dict
        comp_id2ind[comp_type] = Dict(comp["id"]=>comp["index"] for comp in values(comp_dict))
    end

    # update bus references
    for comp_type in components
        for (_, comp) in data_model[comp_type]
            for bus_key in ["f_bus", "t_bus", "bus"]
                if haskey(comp, bus_key)
                    comp[bus_key] = comp_id2ind["bus"][comp[bus_key]]
                end
            end
        end
    end

    return data_model
end


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


function delete_solution!(solution, comp_type, id, props)
    if haskey(solution, comp_type)
        if haskey(solution[comp_type], id)
            for prop in props
                delete!(solution[comp_type][id], prop)
            end
        end
    end
end


function delete_solution!(solution, comp_type, id)
    if haskey(solution, comp_type)
        if haskey(solution[comp_type], id)
            delete!(solution[comp_type], id)
        end
    end
end

##
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


function _get_ground!(bus)
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
function _sc2br_impedance(Zsc)
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
    Zbr = Dict()
    for k in keys(Zsc)
        Zbr[k] = (abs(Y[k...])==0) ? Inf : -1/Y[k...]
    end
    return Zbr
end




""
function _build_loss_model!(data_math::Dict{String,<:Any}, transformer_name::String, to_map::Vector{String}, r_s::Vector{Float64}, zsc::Dict{Tuple{Int,Int},Complex{Float64}}, ysh::Complex{Float64}; nphases::Int=3)
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

    bus_ids = Dict()
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
            "tap" => ones(nphases)
        )

        data_math["branch"]["$(branch_obj["index"])"] = branch_obj

        push!(to_map, "branch.$(branch_obj["index"])")
    end

    return [bus_ids[bus] for bus in tr_t_bus]
end


""
function _kron_reduce_branch!(obj, Zs_keys, Ys_keys, terminals, neutral)
    Zs = [obj[k] for k in Zs_keys]
    Ys = [obj[k] for k in Ys_keys]
    Zs_kr, Ys_kr, terminals_kr = _kron_reduce_branch(Zs, Ys, terminals, neutral)

    for (i,k) in enumerate(Zs_keys)
        obj[k] = Zs_kr[i]
    end

    for (i,k) in enumerate(Ys_keys)
        obj[k] = Ys_kr[i]
    end

    return _get_idxs(terminals, terminals_kr)
end


""
function _kron_reduce_branch(Zs, Ys, terminals, neutral)
    Zs_kr = [deepcopy(Z) for Z in Zs]
    Ys_kr = [deepcopy(Y) for Y in Ys]
    terminals_kr = deepcopy(terminals)

    while neutral in terminals_kr
        n = _get_ilocs(terminals_kr, neutral)[1]
        P = setdiff(collect(1:length(terminals_kr)), n)

        Zs_kr = [Z[P,P]-(1/Z[n,n])*Z[P,[n]]*Z[[n],P] for Z in Zs_kr]
        Ys_kr = [Y[P,P] for Y in Ys_kr]

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
function _pad_properties_delta!(object::Dict{<:Any,<:Any}, properties::Vector{String}, connections::Vector{Int}, phases::Vector{Int}; invert=false)
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


function _apply_filter!(obj, properties, filter)
    for property in properties
        if haskey(obj, property)
            if isa(obj[property], Vector)
                obj[property] = obj[property][filter]
            elseif isa(obj[property], Matrix)
                obj[property] = obj[property][filter, filter]
            else
                error("The property $property is not a Vector or a Matrix!")
            end
        end
    end
end


"""
Given a set of addmittances 'y' connected from the conductors 'f_cnds' to the
conductors 't_cnds', this method will return a list of conductors 'cnd' and a
matrix 'Y', which will satisfy I[cnds] = Y*V[cnds].
"""
function calc_shunt(f_cnds::Vector{Int}, t_cnds::Vector{Int}, y::Vector{Float64})
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
function _calc_ground_shunt_admittance_matrix(cnds::Vector{Int}, Y::Matrix{Float64}, ground::Int)
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


""
function _rm_floating_cnd(cnds::Vector{Int}, Y::Matrix{Float64}, f::Int)
    P = setdiff(cnds, f)
    f_inds = _get_idxs(cnds, [f])
    P_inds = _get_idxs(cnds, P)
    Yrm = Y[P_inds,P_inds]-(1/Y[f_inds,f_inds][1])*Y[P_inds,f_inds]*Y[f_inds,P_inds]
    return (P,Yrm)
end


""
function _expand_linecode!(data_model)
    # expand line codes
    for (id, line) in data_model["line"]
        if haskey(line, "linecode")
            linecode = data_model["linecode"][line["linecode"]]
            for key in ["rs", "xs", "g_fr", "g_to", "b_fr", "b_to"]
                line[key] = line["length"]*linecode[key]
            end
            delete!(line, "linecode")
            delete!(line, "length")
        end
    end
    delete!(data_model, "linecode")
end


""
function _lossy_ground_to_shunt!(data_model)
    mappings = []

    if haskey(data_model, "bus")
        for (id, bus) in data_model["bus"]
            grounding_lossy_inds = [i for (i,t) in enumerate(bus["grounded"]) if bus["rg"][i]!=0 || bus["xg"][i]!=0]
            grounding_lossy = bus["grounded"][grounding_lossy_inds]
            grounding_perfect = bus["grounded"][setdiff(1:length(bus["grounded"]), grounding_lossy_inds)]

            if !isempty(grounding_lossy)
                zg = bus["rg"][grounding_lossy_inds].+im*bus["xg"][grounding_lossy_inds]
                Y_sh = diagm(0=>inv.(zg)) # diagonal matrix, so matrix inverse is element-wise inverse
                add_virtual!(data_model, "shunt", create_shunt(bus=bus["id"], connections=grounding_lossy,
                    g_sh=real.(Y_sh), b_sh=imag.(Y_sh)
                ))
            end
        end
    end
    return mappings
end


""
function _load_to_shunt!(data_model)
    mappings = []
    if haskey(data_model, "load")
        for (id, load) in data_model["load"]
            if load["model"]=="constant_impedance"
                b = load["qd_ref"]./load["vnom"].^2*1E3
                g = load["pd_ref"]./load["vnom"].^2*1E3
                y = b.+im*g
                N = length(b)

                if load["configuration"]=="delta"
                    # create delta transformation matrix Md
                    Md = diagm(0=>ones(N), 1=>-ones(N-1))
                    Md[N,1] = -1
                    Y = Md'*diagm(0=>y)*Md

                else # load["configuration"]=="wye"
                    Y_fr = diagm(0=>y)
                    # B = [[b]; -1'*[b]]*[I -1]
                    Y = vcat(Y_fr, -ones(N)'*Y_fr)*hcat(diagm(0=>ones(N)),  -ones(N))
                end

                shunt = add_virtual!(data_model, "shunt", create_shunt(bus=load["bus"], connections=load["connections"], b_sh=imag.(Y), g_sh=real.(Y)))

                delete_component!(data_model, "load", load)

                push!(mappings, Dict(
                    "load" => load,
                    "shunt_id" => shunt["id"],
                ))
            end
        end
    end

    return mappings
end


""
function _capacitor_to_shunt!(data_model)
    mappings = []

    if haskey(data_model, "capacitor")
        for (id, cap) in data_model["capacitor"]
            b = cap["qd_ref"]./cap["vnom"]^2*1E3
            N = length(b)

            if cap["configuration"]=="delta"
                # create delta transformation matrix Md
                Md = diagm(0=>ones(N), 1=>-ones(N-1))
                Md[N,1] = -1
                B = Md'*diagm(0=>b)*Md

            elseif cap["configuration"]=="wye-grounded"
                B = diagm(0=>b)

            elseif cap["configuration"]=="wye-floating"
                # this is a floating wye-segment
                # B = [b]*(I-1/(b'*1)*[b';...;b'])
                B = diagm(0=>b)*(diagm(0=>ones(N)) - 1/sum(b)*repeat(b',N,1))

            else # cap["configuration"]=="wye"
                B_fr = diagm(0=>b)
                # B = [[b]; -1'*[b]]*[I -1]
                B = vcat(B_fr, -ones(N)'*B_fr)*hcat(diagm(0=>ones(N)),  -ones(N))
            end

            shunt = create_shunt(NaN, cap["bus"], cap["connections"], b_sh=B)
            add_virtual!(data_model, "shunt", shunt)
            delete_component!(data_model, "capacitor", cap)

            push!(mappings, Dict(
                "capacitor" => cap,
                "shunt_id" => shunt["id"],
            ))
        end
    end

    return mappings
end
