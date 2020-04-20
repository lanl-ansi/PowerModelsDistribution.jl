"initializes the base math object of any type, and copies any one-to-one mappings"
function _init_math_obj(obj_type::String, eng_id::Any, eng_obj::Dict{String,<:Any}, index::Int)::Dict{String,Any}
    math_obj = Dict{String,Any}(
        "name" => "$eng_id"
    )

    for key in _1to1_maps[obj_type]
        if haskey(eng_obj, key)
            math_obj[key] = eng_obj[key]
        end
    end

    math_obj["index"] = index

    return math_obj
end


"initializes the base components that are expected by powermodelsdistribution in the mathematical model"
function _init_base_components!(data_math::Dict{String,<:Any})
    for key in ["bus", "load", "shunt", "gen", "branch", "switch", "transformer", "storage", "dcline"]
        if !haskey(data_math, key)
            data_math[key] = Dict{String,Any}()
        end
    end
end


"Initializes the lookup table"
function _init_lookup!(data_math::Dict{String,<:Any})
    for key in keys(_1to1_maps)
        if !haskey(data_math["lookup"], key)
            data_math["lookup"][key] = Dict{Any,Int}()
        end
    end
end


"function for applying a scale to a paramter"
function _scale(dict::Dict{String,<:Any}, key::String, scale::Real)
    if haskey(dict, key)
        dict[key] *= scale
    end
end


"_get_ground helper function"
function _get_new_ground(terminals::Vector{<:Any})
    if isa(terminals, Vector{Int})
        return maximum(terminals)+1
    else
        nrs = [parse(Int, x[1]) for x in [match(r"n([1-9]{1}[0-9]*)", string(t)) for t in terminals] if x !== nothing]
        new = isempty(nrs) ? 1 : maximum(nrs)+1
        if isa(terminals, Vector{Symbol})
            return Symbol("g$new")
        else
            return "g$new"
        end
    end
end


"gets the grounding information for a bus"
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
function _build_loss_model!(data_math::Dict{String,<:Any}, transformer_name::String, to_map::Vector{String}, r_s::Vector{Float64}, zsc::Dict{Tuple{Int,Int},Complex{Float64}}, ysh::Complex{Float64}; nphases::Int=3, kron_reduced=false)::Vector{Int}
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
            "vmin" => fill(0.0, nphases),
            "vmax" => fill(Inf, nphases),
            "grounded" => fill(false, nphases),
            "base_kv" => 1.0,
            "bus_type" => 1,
            "status" => 1,
            "index" => length(data_math["bus"])+1,
        )

        if !kron_reduced
            if bus in tr_t_bus
                bus_obj["terminals"] = collect(1:nphases+1)
                bus_obj["vmin"] = fill(0.0, nphases+1)
                bus_obj["vmax"] = fill(Inf, nphases+1)
                bus_obj["grounded"] = [fill(false, nphases)..., true]
                bus_obj["rg"] = [0.0]
                bus_obj["xg"] = [0.0]
            else
                bus_obj["terminals"] = collect(1:nphases)
                bus_obj["vmin"] = fill(0.0, nphases)
                bus_obj["vmax"] = fill(Inf, nphases)
            end
        end

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


"performs kron reduction on branch"
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


"get locations of terminal in connections list"
function _get_ilocs(vec::Vector{<:Any}, loc::Any)::Vector{Int}
    return collect(1:length(vec))[vec.==loc]
end


"performs kron reduction on branch - helper function"
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


"pads properties to have the total number of conductors for the whole system"
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


"pads properties to have the total number of conductors for the whole system - delta connection variant"
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


"Filters out values of a vector or matrix for certain properties"
function _apply_filter!(obj::Dict{String,<:Any}, properties::Vector{String}, filter::Union{Array,BitArray})
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
function _calc_shunt(f_cnds::Vector{Int}, t_cnds::Vector{Int}, y)::Tuple{Vector{Int}, Matrix{Real}}
    #TODO fix y::Type
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
    if ground in cnds
        cndsr = setdiff(cnds, ground)
        cndsr_inds = _get_idxs(cnds, cndsr)
        Yr = Y[cndsr_inds, cndsr_inds]
        return (cndsr, Yr)
    else
        return cnds, Y
    end
end


"initialization actions for unmapping"
function _init_unmap_eng_obj!(data_eng::Dict{String,<:Any}, eng_obj_type::String, map::Dict{Symbol,Any})::Dict{String,Any}
    if !haskey(data_eng, eng_obj_type)
        data_eng[eng_obj_type] = Dict{Any,Any}()
    end

    eng_obj = Dict{String,Any}()

    return eng_obj
end


"returns component from the mathematical data model"
function _get_math_obj(data_math::Dict{String,<:Any}, to_id::String)::Dict{String,Any}
    math_type, math_id = split(to_id, '.')
    return haskey(data_math, math_type) && haskey(data_math[math_type], math_id) ? data_math[math_type][math_id] : Dict{String,Any}()
end


"convert cost model names"
function _add_gen_cost_model!(math_obj::Dict{String,<:Any}, eng_obj::Dict{String,<:Any})
    math_obj["model"] = get(eng_obj, "cost_pg_model", 2)
    math_obj["startup"] = 0.0
    math_obj["shutdown"] = 0.0
    math_obj["cost"] = get(eng_obj, "cost_pg_parameters", [0.0, 1.0, 0.0])
    math_obj["ncost"] = length(math_obj["cost"])
end


"applies a xfmrcode to a transformer in preparation for converting to mathematical model"
function _apply_xfmrcode!(eng_obj::Dict{String,<:Any}, data_eng::Dict{String,<:Any})
    if haskey(eng_obj, "xfmrcode") && haskey(data_eng, "xfmrcode") && haskey(data_eng["xfmrcode"], eng_obj["xfmrcode"])
        xfmrcode = data_eng["xfmrcode"][eng_obj["xfmrcode"]]

        for (k, v) in xfmrcode
            if !haskey(eng_obj, k)
                eng_obj[k] = v
            elseif haskey(eng_obj, k) && k in ["vnom", "snom", "tm", "rs"]
                for (w, vw) in enumerate(eng_obj[k])
                    if ismissing(vw)
                        eng_obj[k][w] = v[w]
                    end
                end
            end
        end
    end
end


"applies a linecode to a line in preparation for converting to mathematical model"
function _apply_linecode!(eng_obj::Dict{String,<:Any}, data_eng::Dict{String,<:Any})
    if haskey(eng_obj, "linecode") && haskey(data_eng, "linecode") && haskey(data_eng["linecode"], eng_obj["linecode"])
        linecode = data_eng["linecode"][eng_obj["linecode"]]

        for property in ["rs", "xs", "g_fr", "g_to", "b_fr", "b_to"]
            if !haskey(eng_obj, property) && haskey(linecode, property)
                eng_obj[property] = linecode[property]
            end
        end
    end
end


"parses {}_times_series parameters into the expected InfrastructureModels format"
function _parse_time_series_parameter!(data_math::Dict{String,<:Any}, time_series::Dict{String,<:Any}, to_component_type::String, to_component_id::String, fr_parameter::Any, to_parameter::String, conversion_function::Function)
    if !haskey(data_math, "time_series")
        data_math["time_series"] => Dict{String,Any}()
    end

    if !haskey(data_math["time_series"], "time")
        data_math["time_series"]["time_point"] = time_series["time"]
    else
        if time_series["time"] == data_math["time_series"]["time_point"]
            Memento.warn(_LOGGER, "Time series data doesn't match between different objects, aborting")
        end
    end

    data_math["time_series"]["num_steps"] = length(time_series["time"])

    if !haskey(data_math["time_series"], to_component_type)
        data_math["time_series"][to_component_type] = Dict{String,Any}()
    end

    if !haskey(data_math["time_series"][to_component_type], to_component_id)
        data_math["time_series"][to_component_type][to_component_id] = Dict{String,Any}()
    end

    if time_series["replace"]
        data_math["time_series"][to_component_type][to_component_id][to_parameter] = time_series["values"]
    else
        data_math["time_series"][to_component_type][to_component_id][to_parameter] = fr_parameter .* time_series["values"]
    end
end


""
function _no_conversion(data_eng::Dict{String,<:Any}, eng_obj::Dict{String,<:Any}, key::String)
    eng_obj[key]
end


""
function _impedance_conversion(data_eng::Dict{String,<:Any}, eng_obj::Dict{String,<:Any}, key::String)
    eng_obj[key] .* get(eng_obj, "length", 1.0)
end


""
function _admittance_conversion(data_eng::Dict{String,<:Any}, eng_obj::Dict{String,<:Any}, key::String)
    2.0 .* pi .* data_eng["settings"]["base_frequency"] .* eng_obj[key] .* get(eng_obj, "length", 1.0) ./ 1e9
end


""
function _bus_type_conversion(data_eng::Dict{String,<:Any}, eng_obj::Dict{String,<:Any}, key::String)
    eng_obj[key] == 0 ? 4 : 1
end


""
function _vnom_conversion(data_eng::Dict{String,<:Any}, eng_obj::Dict{String,<:Any}, key::String)
    eng_obj[key] ./ data_eng["settings"]["vbase"]
end


""
function _angle_shift_conversion(data_eng::Dict{String,<:Any}, eng_obj::Dict{String,<:Any}, key::String)
    _wrap_to_180([120.0 * i for i in 1:3] .+ eng_obj[key])
end


"lossy grounding to perfect grounding and shunts"
function _convert_grounding(terminals, grounded, rg, xg)
    grouped = Dict(t=>[] for t in unique(grounded))
    for (i,t) in enumerate(grounded)
        push!(grouped[t], rg[i]+im*xg[i])
    end
    t_lookup = Dict(t=>i for (i,t) in enumerate(terminals))
    grounded_lossless = fill(false, length(terminals))
    shunts = []
    for (t, zgs) in grouped
        if any(iszero.(zgs))
            grounded_lossless[t_lookup[t]] = true
        else
            ygs = 1 ./zgs
            yg = sum(ygs)
            push!(shunts, ([t], [yg]))
        end
    end
    return grounded_lossless, shunts
end


"slices branches based on connected terminals"
function _slice_branches!(data_math::Dict{String,<:Any})
    for (_, branch) in data_math["branch"]
        if haskey(branch, "f_connections")
            N = length(branch["f_connections"])
            for prop in ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to"]
                branch[prop] = branch[prop][1:N,1:N]
            end
        end
    end
end
