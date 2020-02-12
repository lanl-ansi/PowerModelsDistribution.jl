import LinearAlgebra


function map_down_data_model(data_model_user)
    data_model_base = deepcopy(data_model_user)
    data_model_base["source_data_model"] = data_model_user

    _expand_linecode!(data_model_base)
    _expand_voltage_zone!(data_model_base)
    _load_to_shunt!(data_model_base)
    _capacitor_to_shunt!(data_model_base)
    _decompose_transformer_nw_lossy!(data_model_base)

    # add low level component types if not present yet
    for comp_type in ["load", "generator", "bus", "line", "shunt", "transformer_2w_ideal"]
        if !haskey(data_model_base, comp_type)
            data_model_base[comp_type] = Dict{String, Any}()
        end
    end
    return data_model_base
end


function _expand_voltage_zone!(data_model)
    # expand line codes
    for (id, bus) in data_model["bus"]
        if haskey(bus, "voltage_zone")
            voltage_zone = data_model["voltage_zone"][bus["voltage_zone"]]
            bus["vnom"] = voltage_zone["vnom"]
            for key in ["vm_ln_min", "vm_ln_max", "vm_lg_min", "vm_lg_max", "vm_ng_min", "vm_ng_max", "vm_ll_min", "vm_ll_max"]
                if haskey(voltage_zone, key)
                    bus[key] = voltage_zone[key]
                end
            end
            delete!(bus, "voltage_zone")
        end
    end
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
end


function _load_to_shunt!(data_model)
    if haskey(data_model, "load")
        for (id, load) in data_model["load"]
            if load["model"]=="constant_impedance"
                b = load["qd"]./load["vm_nom"].^2*1E3
                g = load["pd"]./load["vm_nom"].^2*1E3
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

                shunt = create_shunt(NaN, load["bus"], load["terminals"], b_sh=imag.(Y), g_sh=real.(Y))
                add_component!(data_model, "shunt", shunt)
                delete_component!(data_model, "load", load)

                add_mapping!(data_model, "load_to_shunt", Dict(
                    "load" => data_model["source_data_model"]["load"][id],
                    "shunt" => shunt,
                ))
            end
        end
    end
end


function _capacitor_to_shunt!(data_model)
    if haskey(data_model, "capacitor")
        for (id, cap) in data_model["capacitor"]
            b = cap["qd_ref"]./cap["vm_nom"]^2*1E3
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

            shunt = create_shunt(NaN, cap["bus"], cap["terminals"], b_sh=B)
            add_component!(data_model, "shunt", shunt)
            delete_component!(data_model, "capacitor", cap)

            add_mapping!(data_model, "capacitor_to_shunt", Dict(
                "capacitor" => data_model["source_data_model"]["capacitor"][id],
                "shunt" => shunt,
            ))
        end
    end
end


# test
"""

    function decompose_transformer_nw_lossy!(data_model)

Replaces complex transformers with a composition of ideal transformers and lines
which model losses. New buses (virtual, no physical meaning) are added.
"""
function _decompose_transformer_nw_lossy!(data_model)
    for (tr_id, trans) in data_model["transformer_nw_lossy"]

        vnom = trans["vnom"]*data_model["v_var_scalar"]
        snom = trans["snom"]*data_model["v_var_scalar"]

        nrw = length(trans["buses"])

        # calculate zbase in which the data is specified, and convert to SI
        zbase = (vnom.^2)./snom
        # x_sc is specified with respect to first winding
        x_sc = trans["xsc"].*zbase[1]
        # rs is specified with respect to each winding
        r_s = trans["rs"].*zbase

        g_sh = (trans["noloadloss"]*snom[1]/3)/vnom[1]^2
        b_sh = (trans["imag"]*snom[1]/3)/vnom[1]^2

        # data is measured externally, but we now refer it to the internal side
        ratios = vnom/1E3
        x_sc = x_sc./ratios[1]^2
        r_s = r_s./ratios.^2
        g_sh = g_sh*ratios[1]^2
        b_sh = b_sh*ratios[1]^2

        # convert x_sc from list of upper triangle elements to an explicit dict
        y_sh = g_sh + im*b_sh
        z_sc = Dict([(key, im*x_sc[i]) for (i,key) in enumerate([(i,j) for i in 1:nrw for j in i+1:nrw])])

        vbuses, vlines, trans_t_bus_w = _build_loss_model!(data_model, r_s, z_sc, y_sh)

        trans_w = Array{Dict, 1}(undef, nrw)
        for w in 1:nrw
            # 2-WINDING TRANSFORMER
            # make virtual bus and mark it for reduction
            tm_nom = trans["configuration"][w]=="delta" ? trans["vnom"][w]*sqrt(3) : trans["vnom"][w]
            trans_w[w] = create_transformer_2w_ideal(NaN,
                trans["buses"][w], trans_t_bus_w[w], tm_nom,
                f_terminals     = trans["terminals"][w],
                t_terminals     = collect(1:4),
                configuration   = trans["configuration"][w],
                polarity        = trans["polarity"][w],
                #tm_set         = trans["tm_set"][w],
                tm_fix          = trans["tm_fix"][w],
                #tm_max         = trans["tm_max"][w],
                #tm_min         = trans["tm_min"][w],
                #tm_step        = trans["tm_step"][w],
            )

            add_component!(data_model, "transformer_2w_ideal", trans_w[w])
        end

        delete_component!(data_model, "transformer_nw_lossy", trans)

        add_mapping!(data_model, "transformer_decomposition", Dict(
            "trans"=>data_model["source_data_model"]["transformer_nw_lossy"][tr_id],
            "trans_w"=>trans_w,
            "vlines"=>vlines,
            "vbuses"=>vbuses,
        ))
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
        bus_ids[bus] = add_component!(data_model, "bus", create_bus(""))
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
        line_ids[l] = add_component!(data_model, "line", create_line("",
            bus_ids[i], bus_ids[j], n_phases,
            f_terminals = collect(1:n_phases),
            t_terminals = collect(1:n_phases),
            rs = LinearAlgebra.diagm(0=>fill(real(z[l]), n_phases)),
            xs = LinearAlgebra.diagm(0=>fill(imag(z[l]), n_phases)),
            g_fr = LinearAlgebra.diagm(0=>fill(g_fr, n_phases)),
            b_fr = LinearAlgebra.diagm(0=>fill(b_fr, n_phases)),
            g_to = LinearAlgebra.diagm(0=>fill(g_to, n_phases)),
            b_to = LinearAlgebra.diagm(0=>fill(b_to, n_phases)),
        ))
    end

    return bus_ids, line_ids, [bus_ids[bus] for bus in tr_t_bus]
end
