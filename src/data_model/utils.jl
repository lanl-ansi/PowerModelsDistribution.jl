const _pmd_eng_global_keys = Set{String}([
    "files", "name", "data_model", "dss_options"
])


"Hacky helper function to transform single-conductor network data, from, e.g., matpower/psse, into multi-conductor data"
function make_multiconductor!(data::Dict{String,<:Any}, conductors::Int)
    if ismultinetwork(data)
        for (i,nw_data) in data["nw"]
            _make_multiconductor!(nw_data, conductors)
        end
    else
         _make_multiconductor!(data, conductors)
    end
end


"Hacky helper function to transform single-conductor network data, from, e.g., matpower/psse, into multi-conductor data"
function _make_multiconductor!(data::Dict{String,<:Any}, conductors::Real)
    if haskey(data, "conductor_ids")
        @warn "skipping network that is already multiconductor"
        return
    end

    data["conductor_ids"] = collect(1:conductors)
    data["data_model"] = MATHEMATICAL
    data["settings"] = get(
        data,
        "settings",
        Dict{String,Any}(
            "sbase_default" => get(data, "baseMVA", 1e6)
        )
    )

    for (key, item) in data
        if isa(item, Dict{String,Any})
            for (item_id, item_data) in item
                if isa(item_data, Dict{String,Any})
                    item_ref_data = Dict{String,Any}()
                    for (param, value) in item_data
                        if param in _conductorless
                            item_ref_data[param] = value
                        else
                            if param in _conductor_matrix
                                item_ref_data[param] = LinearAlgebra.diagm(0=>fill(value, conductors))
                            else
                                item_ref_data[param] = fill(value, conductors)
                            end
                        end
                    end
                    item[item_id] = item_ref_data
                end
            end
        else
            #root non-dict items
        end
    end

    for (_, load) in data["load"]
        load["model"] = POWER
        load["configuration"] = WYE
    end

    for (_, gen) in data["gen"]
        gen["configuration"] = WYE
    end

    for type in ["load", "gen", "storage", "shunt"]
        if haskey(data, type)
            for (_,obj) in data[type]
                obj["connections"] = collect(1:conductors)
            end
        end
    end

    for type in ["branch", "transformer", "switch"]
        if haskey(data, type)
            for (_,obj) in data[type]
                obj["f_connections"] = collect(1:conductors)
                obj["t_connections"] = collect(1:conductors)
            end
        end
    end

    for (_,bus) in data["bus"]
        bus["terminals"] = collect(1:conductors)
        bus["grounded"] = fill(false, conductors)
    end
end




"initializes the base math object of any type, and copies any one-to-one mappings"
function _init_math_obj(obj_type::String, eng_id::Any, eng_obj::Dict{String,<:Any}, index::Int)::Dict{String,Any}
    math_obj = Dict{String,Any}(
        "name" => "$eng_id",
        "source_id" => "$obj_type.$eng_id"
    )

    for key in _1to1_maps[obj_type]
        if haskey(eng_obj, key)
            if key in ["status", "dispatchable"]
                math_obj[key] = Int(eng_obj[key])
            else
                math_obj[key] = eng_obj[key]
            end
        end
    end

    math_obj["index"] = index

    return math_obj
end


"initializes the base components that are expected by powermodelsdistribution in the mathematical model"
function _init_base_components!(data_math::Dict{String,<:Any})
    for key in pmd_math_asset_types
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
                    error("Short-circuit impedance between winding $i and $j is missing.")
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
function _build_loss_model!(
    data_math::Dict{String,<:Any},
    transformer_name::Any,
    to_map::Vector{String},
    r_s::Vector{Float64},
    zsc::Dict{Tuple{Int,Int},Complex{Float64}},
    ysh::Complex{Float64};
    nphases::Int=3
    )::Vector{Int}

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
            "terminals" => collect(1:nphases),
            "grounded" => fill(false, nphases),
            "base_kv" => 1.0,
            "bus_type" => 1,
            "status" => 1,
            "source_id" => "transformer.$(transformer_name)",
            "index" => length(data_math["bus"])+1,
        )

        if !get(data_math, "is_kron_reduced", false)
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
function _pad_properties!(object::Dict{String,<:Any}, properties::Vector{String}, connections::Vector{Int}, phases::Vector{Int}; pad_value::Real=0.0)
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


"pads properties to have the total number of conductors for the whole system (transformer winding variant)"
function _pad_properties!(object::Dict{String,<:Any}, properties::Vector{String}, wdg::Int, connections::Vector{Int}, phases::Vector{Int}; pad_value::Real=0.0)
    @assert(all(c in phases for c in connections))
    inds = _get_idxs(phases, connections)

    for property in properties
        if haskey(object, property)
            if isa(object[property][wdg], Vector)
                tmp = fill(pad_value, length(phases))
                tmp[inds] = object[property][wdg]
                object[property][wdg] = tmp
            elseif isa(object[property][wdg], Matrix)
                tmp = fill(pad_value, length(phases), length(phases))
                tmp[inds, inds] = object[property][wdg]
                object[property][wdg] = tmp
            end
        end
    end
end


"pads properties to have the total number of conductors for the whole system - delta connection variant"
function _pad_properties_delta!(object::Dict{String,<:Any}, properties::Vector{String}, connections::Vector{Int}, phases::Vector{Int}; invert::Bool=false)
    @assert(all(c in phases for c in connections))
    @assert(length(connections) in [2, 3], "A delta configuration has to have at least 2 or 3 connections!")
    @assert(length(phases)==3, "Padding only possible to a |phases|==3!")

    for property in properties
        if haskey(object, property)
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
end


"pads properties to have the total number of conductors for the whole system - delta connection variant"
function _pad_properties_delta!(object::Dict{String,<:Any}, properties::Vector{String}, connections::Vector{Int}, wdg::Int, phases::Vector{Int}; invert::Bool=false)
    @assert(all(c in phases for c in connections))
    @assert(length(connections) in [2, 3], "A delta configuration has to have at least 2 or 3 connections!")
    @assert(length(phases)==3, "Padding only possible to a |phases|==3!")

    for property in properties
        if haskey(object, property)
            val = object[property][wdg]
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

            object[property][wdg] = [get_val(phases[1], phases[2]), get_val(phases[2], phases[3]), get_val(phases[3], phases[1])]
        end
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
                error("The property $property is not a Vector or a Matrix!")
            end
        end
    end
end


"Filters out values of a vector or matrix for certain properties (transformer winding variant)"
function _apply_filter!(obj::Dict{String,<:Any}, properties::Vector{String}, wdg::Int, filter::Union{Array,BitArray})
    for property in properties
        if haskey(obj, property)
            if isa(obj[property], Vector)
                obj[property][wdg] = obj[property][wdg][filter]
            elseif isa(obj[property], Matrix)
                obj[property][wdg] = obj[property][wdg][filter, filter]
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
function _calc_shunt(f_cnds::Vector{Int}, t_cnds::Vector{Int}, y::Vector{<:Union{Real,Vector{<:Real}}})::Tuple{Vector{Int}, Matrix{Real}}
    cnds = unique([f_cnds..., t_cnds...])
    e(f,t) = reshape([c==f ? 1 : c==t ? -1 : 0 for c in cnds], length(cnds), 1)
    Y = sum([e(f_cnds[i], t_cnds[i])*y[i]*e(f_cnds[i], t_cnds[i])' for (i,_y) in enumerate(y)])
    return (cnds, Y)
end


"""
Given a set of terminals 'cnds' with associated shunt admittance 'Y', this
method will calculate the reduced admittance matrix if terminal 'ground' is
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
function _init_unmap_eng_obj!(data_eng::Dict{String,<:Any}, eng_obj_type::String, map::Dict{String,<:Any})::Dict{String,Any}
    if !haskey(data_eng, eng_obj_type)
        data_eng[eng_obj_type] = Dict{String,Any}()
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
                eng_obj[k] = deepcopy(v)
            elseif haskey(eng_obj, k) && k in ["vm_nom", "sm_nom", "tm_set", "rw"]
                for (w, vw) in enumerate(eng_obj[k])
                    if ismissing(vw)
                        eng_obj[k][w] = deepcopy(v[w])
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
                eng_obj[property] = deepcopy(linecode[property])
            end
        end
    end
end


"converts impendance in Ohm/m by multiplying by length"
function _impedance_conversion(data_eng::Dict{String,<:Any}, eng_obj::Dict{String,<:Any}, key::String)
    eng_obj[key] .* get(eng_obj, "length", 1.0)
end


"converts admittance by multiplying by 2πωl"
function _admittance_conversion(data_eng::Dict{String,<:Any}, eng_obj::Dict{String,<:Any}, key::String)
    2.0 .* pi .* data_eng["settings"]["base_frequency"] .* eng_obj[key] .* get(eng_obj, "length", 1.0) ./ 1e9
end


"converts Int(Status) Enums into bus_type"
function _bus_type_conversion(data_eng::Dict{String,<:Any}, eng_obj::Dict{String,<:Any}, key::String)
    eng_obj[key] == 0 ? 4 : 1
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


"transformations might have introduced buses with four-terminals; crop here"
function _kron_reduce_buses!(data_math)
    for (_, bus) in data_math["bus"]
        for prop in ["vm", "va", "vmax", "vmin"]
            if haskey(bus, prop) && length(bus[prop])>3
                bus[prop] = bus[prop][1:3]
            end
        end
    end
end


"generate a new, unique terminal"
_new_terminal(terms) = maximum([terms[isa.(terms, Int)]..., 3])+1


"get a grounded terminal from a bus; if not present, create one"
function _get_ground_math!(bus; exclude_terminals=[])
    tgs = setdiff(bus["terminals"][bus["grounded"]], exclude_terminals)
    if !isempty(tgs)
        return tgs[1]
    else
        n = _new_terminal([bus["terminals"]])
        push!(bus["terminals"], n)
        push!(bus["grounded"], true)
        return n
    end
end


"finds maximal set of ungrounded phases"
function _get_complete_conductor_set(data::Dict{String,<:Any})
    conductors = Set([])
    for (_, obj) in data["bus"]
        for t in obj["terminals"]
            push!(conductors, t)
        end
    end

    return sort([c for c in conductors])
end


"checks if data structures are equivalent, and if not, will enumerate the differences"
function _check_equal(data1::Dict{String,<:Any}, data2::Dict{String,<:Any}; context::String="", ignore::Vector{String}=Vector{String}(["connections", "f_connections", "t_connections", "terminals"]))
    lines = []
    for (i, obj) in data1
        if !haskey(data2, i)
            push!(lines, "  $i missing in $data2")
        else
            if obj != data2[i]
                if isa(obj, Dict)
                    push!(lines, _check_equal(obj, data2[i]; context="  $context $i"))
                else
                    if !(i in ignore)
                        push!(lines, "  $context $i: $obj != $(data2[i])")
                    end
                end
            end
        end
    end

    return lines
end


"adds conductors to connections during padding process"
function _pad_connections!(eng_obj::Dict{String,<:Any}, connection_key::String, conductors::Union{Vector{Int},Vector{String}})
    for cond in conductors
        if !(cond in eng_obj[connection_key])
            push!(eng_obj[connection_key], cond)
        end
    end
end


"adds conductors to connections during padding process, transformer winding variant"
function _pad_connections!(eng_obj::Dict{String,<:Any}, connection_key::String, wdg::Int, conductors::Union{Vector{Int},Vector{String}})
    for cond in conductors
        if !(cond in eng_obj[connection_key][wdg])
            push!(eng_obj[connection_key][wdg], cond)
        end
    end
end


"helper function to map non integer conductor ids into integers"
function _map_conductor_ids!(data_math::Dict{String,<:Any})
    if all(typeof(c) <: Int for c in data_math["conductor_ids"])
        cnd_map = Dict{Any,Int}(c => c for c in data_math["conductor_ids"])
    else
        cnd_map = Dict{Any,Int}(c => idx for (idx, c) in enumerate(data_math["conductor_ids"]))
    end

    data_math["conductor_ids"] = Vector{Int}([cnd_map[c] for c in data_math["conductor_ids"]])

    for type in ["branch", "switch", "transformer"]
        if haskey(data_math, type)
            for (_,obj) in data_math[type]
                obj["f_connections"] = Vector{Int}([cnd_map[c] for c in obj["f_connections"]])
                obj["t_connections"] = Vector{Int}([cnd_map[c] for c in obj["t_connections"]])
            end
        end
    end

    for type in ["load", "shunt", "gen", "storage"]
        if haskey(data_math, type)
            for (_,obj) in data_math[type]
                obj["connections"] = Vector{Int}([cnd_map[c] for c in obj["connections"]])
            end
        end
    end

    for (_,bus) in data_math["bus"]
        bus["terminals"] = Vector{Int}([cnd_map[t] for t in bus["terminals"]])
    end
end
