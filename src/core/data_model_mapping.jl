import LinearAlgebra

# MAP DATA MODEL DOWN

function data_model_map!(data_model)

    !haskey(data_model, "mappings")

    # needs to happen before _expand_linecode, as it might contain a linecode for the internal impedance
    add_mappings!(data_model, "decompose_voltage_source", _decompose_voltage_source!(data_model))
    _expand_linecode!(data_model)
    # creates shunt of 4x4; disabled for now (incompatible 3-wire kron-reduced)
    #add_mappings!(data_model, "load_to_shunt", _load_to_shunt!(data_model))
    add_mappings!(data_model, "capacitor_to_shunt", _capacitor_to_shunt!(data_model))
    add_mappings!(data_model, "decompose_transformer_nw", _decompose_transformer_nw!(data_model))
    add_mappings!(data_model, "_lossy_ground_to_shunt", _lossy_ground_to_shunt!(data_model))

    # add low level component types if not present yet
    for comp_type in ["load", "generator", "bus", "line", "shunt", "transformer_2wa", "storage", "switch"]
        if !haskey(data_model, comp_type)
            data_model[comp_type] = Dict{String, Any}()
        end
    end
    return data_model
end


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


function _lossy_ground_to_shunt!(data_model)
    mappings = []

    if haskey(data_model, "bus")
        for (id, bus) in data_model["bus"]
            grounding_lossy_inds = [i for (i,t) in enumerate(bus["grounded"]) if bus["rg"][i]!=0 || bus["xg"][i]!=0]
            grounding_lossy = bus["grounded"][grounding_lossy_inds]
            grounding_perfect = bus["grounded"][setdiff(1:length(bus["grounded"]), grounding_lossy_inds)]

            if !isempty(grounding_lossy)
                zg = bus["rg"][grounding_lossy_inds].+im*bus["xg"][grounding_lossy_inds]
                Y_sh = LinearAlgebra.diagm(0=>inv.(zg)) # diagonal matrix, so matrix inverse is element-wise inverse
                add_virtual!(data_model, "shunt", create_shunt(bus=bus["id"], connections=grounding_lossy,
                    g_sh=real.(Y_sh), b_sh=imag.(Y_sh)
                ))
            end
        end
    end
    return mappings
end


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
                    Md = LinearAlgebra.diagm(0=>ones(N), 1=>-ones(N-1))
                    Md[N,1] = -1
                    Y = Md'*LinearAlgebra.diagm(0=>y)*Md

                else # load["configuration"]=="wye"
                    Y_fr = LinearAlgebra.diagm(0=>y)
                    # B = [[b]; -1'*[b]]*[I -1]
                    Y = vcat(Y_fr, -ones(N)'*Y_fr)*hcat(LinearAlgebra.diagm(0=>ones(N)),  -ones(N))
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


function _capacitor_to_shunt!(data_model)
    mappings = []

    if haskey(data_model, "capacitor")
        for (id, cap) in data_model["capacitor"]
            b = cap["qd_ref"]./cap["vnom"]^2*1E3
            N = length(b)

            if cap["configuration"]=="delta"
                # create delta transformation matrix Md
                Md = LinearAlgebra.diagm(0=>ones(N), 1=>-ones(N-1))
                Md[N,1] = -1
                B = Md'*LinearAlgebra.diagm(0=>b)*Md

            elseif cap["configuration"]=="wye-grounded"
                B = LinearAlgebra.diagm(0=>b)

            elseif cap["configuration"]=="wye-floating"
                # this is a floating wye-segment
                # B = [b]*(I-1/(b'*1)*[b';...;b'])
                B = LinearAlgebra.diagm(0=>b)*(LinearAlgebra.diagm(0=>ones(N)) - 1/sum(b)*repeat(b',N,1))

            else # cap["configuration"]=="wye"
                B_fr = LinearAlgebra.diagm(0=>b)
                # B = [[b]; -1'*[b]]*[I -1]
                B = vcat(B_fr, -ones(N)'*B_fr)*hcat(LinearAlgebra.diagm(0=>ones(N)),  -ones(N))
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


function _decompose_voltage_source!(data_model)
    mappings = []

    if haskey(data_model, "voltage_source")
        for (id, vs) in data_model["voltage_source"]

            bus = data_model["bus"][vs["bus"]]

            line_kwargs = Dict(Symbol(prop)=>vs[prop] for prop in ["rs", "xs", "g_fr", "b_fr", "g_to", "b_to", "linecode", "length"] if haskey(vs, prop))
            lossy = !isempty(line_kwargs)

            # if any loss parameters (or linecode) were supplied, then create a line and internal bus
            if lossy
                sourcebus = add_virtual!(data_model, "bus", create_bus(terminals=deepcopy(vs["connections"])))

                line = add_virtual!(data_model, "line", create_line(;
                    f_bus=sourcebus["id"], f_connections=vs["connections"], t_bus=bus["id"], t_connections=vs["connections"],
                    line_kwargs...
                ))
            else
                sourcebus = bus
            end

            ground = _get_ground!(sourcebus)
            gen = create_generator(bus=sourcebus["id"], connections=[vs["connections"]..., ground])

            for prop in ["pg_max", "pg_min", "qg_max", "qg_min"]
                if haskey(vs, prop)
                    gen[prop] = vs[prop]
                end
            end

            add_virtual!(data_model, "generator", gen)

            conns = vs["connections"]
            terminals = bus["terminals"]

            tmp = Dict(enumerate(conns))
            sourcebus["vm"] = sourcebus["vmax"] = sourcebus["vmin"] = [haskey(tmp, t) ? vs["vm"][tmp[t]] : NaN for t in terminals]
            sourcebus["va"] = [haskey(tmp, t) ? vs["va"][tmp[t]] : NaN for t in terminals]
            sourcebus["bus_type"] = 3

            delete_component!(data_model, "voltage_source", vs["id"])
            push!(mappings, Dict("voltage_source"=>vs, "gen_id"=>gen["id"],
                "vbus_id"  => lossy ? sourcebus["id"] : nothing,
                "vline_id" => lossy ? line["id"]      : nothing,
            ))
        end
    end

    return mappings
end


"""

    function decompose_transformer_nw_lossy!(data_model)

Replaces complex transformers with a composition of ideal transformers and lines
which model losses. New buses (virtual, no physical meaning) are added.
"""
function _decompose_transformer_nw!(data_model)
    mappings = []

    if haskey(data_model, "transformer_nw")
        for (tr_id, trans) in data_model["transformer_nw"]

            vnom = trans["vnom"]*data_model["settings"]["v_var_scalar"]
            snom = trans["snom"]*data_model["settings"]["v_var_scalar"]

            nrw = length(trans["bus"])

            # calculate zbase in which the data is specified, and convert to SI
            zbase = (vnom.^2)./snom
            # x_sc is specified with respect to first winding
            x_sc = trans["xsc"].*zbase[1]
            # rs is specified with respect to each winding
            r_s = trans["rs"].*zbase

            g_sh =  (trans["noloadloss"]*snom[1])/vnom[1]^2
            b_sh = -(trans["imag"]*snom[1])/vnom[1]^2

            # data is measured externally, but we now refer it to the internal side
            ratios = vnom/data_model["settings"]["v_var_scalar"]
            x_sc = x_sc./ratios[1]^2
            r_s = r_s./ratios.^2
            g_sh = g_sh*ratios[1]^2
            b_sh = b_sh*ratios[1]^2

            # convert x_sc from list of upper triangle elements to an explicit dict
            y_sh = g_sh + im*b_sh
            z_sc = Dict([(key, im*x_sc[i]) for (i,key) in enumerate([(i,j) for i in 1:nrw for j in i+1:nrw])])

            vbuses, vlines, trans_t_bus_w = _build_loss_model!(data_model, r_s, z_sc, y_sh)

            trans_ids_w = Array{String, 1}(undef, nrw)
            for w in 1:nrw
                # 2-WINDING TRANSFORMER
                # make virtual bus and mark it for reduction
                tm_nom = trans["configuration"][w]=="delta" ? trans["vnom"][w]*sqrt(3) : trans["vnom"][w]
                trans_2wa = add_virtual!(data_model, "transformer_2wa", Dict(
                    "f_bus"         => trans["bus"][w],
                    "t_bus"         => trans_t_bus_w[w],
                    "tm_nom"        => tm_nom,
                    "f_connections" => trans["connections"][w],
                    "t_connections" => collect(1:4),
                    "configuration" => trans["configuration"][w],
                    "polarity"      => trans["polarity"][w],
                    "tm"            => trans["tm"][w],
                    "fixed"         => trans["fixed"][w],
                ))

                for prop in ["tm_min", "tm_max", "tm_step"]
                    if haskey(trans, prop)
                        trans_2wa[prop] = trans[prop][w]
                    end
                end

                trans_ids_w[w] = trans_2wa["id"]
            end

            delete_component!(data_model, "transformer_nw", trans)

            push!(mappings, Dict(
                "trans"=>trans,
                "trans_2wa"=>trans_ids_w,
                "vlines"=>vlines,
                "vbuses"=>vbuses,
            ))
        end
    end

    return mappings
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


function _build_loss_model!(data_model, r_s, zsc, ysh; n_phases=3)
    # precompute the minimal set of buses and lines
    N = length(r_s)
    tr_t_bus = collect(1:N)
    buses = Set(1:2*N)
    edges = [[[i,i+N] for i in 1:N]..., [[i+N,j+N] for (i,j) in keys(zsc)]...]
    lines = Dict(enumerate(edges))
    z = Dict(enumerate([r_s..., values(zsc)...]))
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
                if edge[1]==j
                    edge[1] = i
                end
                if edge[2]==j
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
                if tr_t_bus[w]==j
                    tr_t_bus[w] = i
                end
            end
        end
    end

    bus_ids = Dict()
    for bus in buses
        bus_ids[bus] = add_virtual_get_id!(data_model, "bus", create_bus(id=""))
    end
    line_ids = Dict()
    for (l,(i,j)) in lines
        # merge the shunts into the shunts of the pi model of the line
        g_fr = b_fr = g_to = b_to = 0
        if haskey(shunts, i)
            g_fr = real(shunts[i])
            b_fr = imag(shunts[i])
            delete!(shunts, i)
        end
        if haskey(shunts, j)
            g_fr = real(shunts[j])
            b_fr = imag(shunts[j])
            delete!(shunts, j)
        end
        line_ids[l] = add_virtual_get_id!(data_model, "line", Dict(
            "status"=>1,
            "f_bus"=>bus_ids[i], "t_bus"=>bus_ids[j],
            "f_connections"=>collect(1:n_phases),
            "t_connections"=>collect(1:n_phases),
            "rs"=>LinearAlgebra.diagm(0=>fill(real(z[l]), n_phases)),
            "xs"=>LinearAlgebra.diagm(0=>fill(imag(z[l]), n_phases)),
            "g_fr"=>LinearAlgebra.diagm(0=>fill(g_fr, n_phases)),
            "b_fr"=>LinearAlgebra.diagm(0=>fill(b_fr, n_phases)),
            "g_to"=>LinearAlgebra.diagm(0=>fill(g_to, n_phases)),
            "b_to"=>LinearAlgebra.diagm(0=>fill(b_to, n_phases)),
        ))
    end

    return bus_ids, line_ids, [bus_ids[bus] for bus in tr_t_bus]
end

function _alias!(dict, fr, to)
    if haskey(dict, fr)
        dict[to] = dict[fr]
    end
end

function _pad_props!(comp, keys, phases_comp, phases_all)
    pos = Dict((x,i) for (i,x) in enumerate(phases_all))
    inds = [pos[x] for x in phases_comp]
    for prop in keys
        if haskey(comp, prop)
            if isa(comp[prop], Vector)
                tmp = zeros(length(phases_all))
                tmp[inds] = comp[prop]
                comp[prop] = tmp
            elseif isa(comp[prop], Matrix)
                tmp = zeros(length(phases_all), length(phases_all))
                tmp[inds, inds] = comp[prop]
                comp[prop] = tmp
            else
                error("Property is not a vector or matrix!")
            end
        end
    end
end

function data_model_make_compatible_v8!(data_model; phases=[1, 2, 3], neutral=4)
    data_model["conductors"] = 3
    data_model["buspairs"] = nothing
    for (_, bus) in data_model["bus"]
        bus["bus_i"] = bus["index"]
        terminals = bus["terminals"]
        @assert(all(t in [phases..., neutral] for t in terminals))
        for prop in ["vm", "va", "vmin", "vmax"]
            if haskey(bus, prop)
                if length(bus[prop])==4
                    val = bus[prop]
                    bus[prop] = val[terminals.!=neutral]
                end
            end
        end
    end

    for (_, load) in data_model["load"]
        # remove neutral
        if load["configuration"]=="wye"
            bus = data_model["bus"][string(load["bus"])]
            @assert(length(bus["grounded"])==1 && bus["grounded"][1]==load["connections"][end])
            load["connections"] = load["connections"][1:end-1]
            _pad_props!(load, ["pd", "qd"], load["connections"], phases)
        else
            # three-phase loads can only be delta-connected
            #@assert(all(load["connections"].==phases))
        end
        _alias!(load, "bus", "load_bus")
    end

    data_model["gen"] = data_model["generator"]

    # has to be three-phase
    for (_, gen) in data_model["gen"]
        if gen["configuration"]=="wye"
            @assert(all(gen["connections"].==[phases..., neutral]))
        else
            @assert(all(gen["connections"].==phases))
        end

        _alias!(gen, "status", "gen_status")
        _alias!(gen, "bus", "gen_bus")
        _alias!(gen, "pg_min", "pmin")
        _alias!(gen, "qg_min", "qmin")
        _alias!(gen, "pg_max", "pmax")
        _alias!(gen, "qg_max", "qmax")
        _alias!(gen, "configuration", "conn")

        gen["model"] = 2
    end

    data_model["branch"] = data_model["line"]
    for (_, br) in data_model["branch"]
        @assert(all(x in phases for x in br["f_connections"]))
        @assert(all(x in phases for x in br["t_connections"]))
        @assert(all(br["f_connections"].==br["t_connections"]))

        _pad_props!(br, ["rs", "xs", "g_fr", "g_to", "b_fr", "b_to", "s_rating", "c_rating"], br["f_connections"], phases)

        # rename
        _alias!(br, "status", "br_status")
        _alias!(br, "rs", "br_r")
        _alias!(br, "xs", "br_x")

        br["tap"] = 1.0
        br["shift"] = 0

        if !haskey(br, "angmin")
            N = size(br["br_r"])[1]
            br["angmin"] = fill(-pi/2, N)
            br["angmax"] = fill(pi/2, N)
        end
    end

    for (_, shunt) in data_model["shunt"]
        @assert(all(x in phases for x in shunt["connections"]))
        _pad_props!(shunt, ["g_sh", "b_sh"], shunt["connections"], phases)
        _alias!(shunt, "bus", "shunt_bus")
        _alias!(shunt, "g_sh", "gs")
        _alias!(shunt, "b_sh", "bs")
    end

    data_model["dcline"] = Dict()
    data_model["transformer"] = data_model["transformer_2wa"]

    data_model["per_unit"] = true
    data_model["baseMVA"] = data_model["settings"]["sbase"]*data_model["settings"]["v_var_scalar"]/1E6
    data_model["name"] = "IDC"


    return data_model
end

# MAP SOLUTION UP

function solution_unmap!(solution::Dict, data_model::Dict)
    for i in length(data_model["mappings"]):-1:1
        (name, data) = data_model["mappings"][i]

        if name=="decompose_transformer_nw"
            for bus_id in values(data["vbuses"])
                delete!(solution["bus"], bus_id)
            end

            for line_id in values(data["vlines"])
                delete!(solution["branch"], line_id)
            end

            pt = [solution["transformer"][tr_id]["pf"] for tr_id in data["trans_2wa"]]
            qt = [solution["transformer"][tr_id]["qf"] for tr_id in data["trans_2wa"]]
            for tr_id in data["trans_2wa"]
                delete!(solution["transformer"], tr_id)
            end

            add_solution!(solution, "transformer_nw", data["trans"]["id"], Dict("pt"=>pt, "qt"=>qt))
        elseif name=="capacitor_to_shunt"
            # shunt has no solutions defined
            delete_solution!(solution, "shunt", data["shunt_id"])
            add_solution!(solution, "capacitor", data["capacitor"]["id"], Dict())
        elseif name=="load_to_shunt"
            # shunt has no solutions, but a load should have!
            delete!(solution, "shunt", data["shunt_id"])
            add_solution!(solution, "load", data["load"]["id"], Dict())
        elseif name=="decompose_voltage_source"
            gen = solution["gen"][data["gen_id"]]
            delete_solution!(solution, "gen", data["gen_id"])
            add_solution!(solution, "voltage_source", data["voltage_source"]["id"], Dict("pg"=>gen["pg"], "qg"=>gen["qg"]))

        end
    end

    # remove component dicts if empty
    for (comp_type, comp_dict) in solution
        if isa(comp_dict, Dict) && isempty(comp_dict)
            delete!(solution, comp_type)
        end
    end
end
