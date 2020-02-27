# OpenDSS parser
import LinearAlgebra: isdiag, diag, pinv


const _exclude_duplicate_check = ["options", "filename", "circuit"]
const _dss_edge_components = ["line", "transformer", "reactor"]
const _dss_supported_components = ["line", "linecode", "load", "generator", "capacitor", "reactor", "circuit", "transformer", "pvsystem", "storage", "loadshape"]
const _dss_option_dtypes = Dict{String,Type}("defaultbasefreq" => Float64, "voltagebases" => Float64)


"Discovers all of the buses (not separately defined in OpenDSS), from \"lines\""
function _discover_buses(data_dss::Dict{String,<:Any})::Array
    bus_names = []
    buses = []
    for obj_type in _dss_edge_components
        if haskey(data_dss, obj_type)
            dss_objs = data_dss[obj_type]
            for dss_obj in values(dss_objs)
                if obj_type == "transformer"
                    dss_obj = _create_transformer(dss_obj["name"]; _to_sym_keys(dss_obj)...)
                    for bus in dss_obj["buses"]
                        name, nodes = _parse_busname(bus)
                        if !(name in bus_names)
                            push!(bus_names, name)
                            push!(buses, (name, nodes))
                        end
                    end
                elseif haskey(dss_obj, "bus2")
                    for key in ["bus1", "bus2"]
                        name, nodes = _parse_busname(dss_obj[key])
                        if !(name in bus_names)
                            push!(bus_names, name)
                            push!(buses, (name, nodes))
                        end
                    end
                end
            end
        end
    end
    if length(buses) == 0
        Memento.error(_LOGGER, "data_dss has no edge components!")
    else
        return buses
    end
end


"Parses buscoords [lon,lat] (if present) into their respective buses"
function _dss2eng_buscoords!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any})
    for (name, coords) in get(data_dss, "buscoords", Dict{String,Any}())
        if haskey(data_eng["bus"], name)
            bus = data_eng["bus"][name]
            bus["lon"] = coords["x"]
            bus["lat"] = coords["y"]
        end
    end
end


"""
    _dss2pmd_bus!(data_eng, data_dss)

Adds PowerModels-style buses to `data_eng` from `data_dss`.
"""
function _dss2eng_bus!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool=false)
    buses = _discover_buses(data_dss)
    for (n, (bus, nodes)) in enumerate(buses)

        @assert(!(length(bus)>=8 && bus[1:8]=="_virtual"), "Bus $bus: identifiers should not start with _virtual.")

        add_bus!(data_eng, id=bus, status=1, bus_type=1)
    end
end


""
function _dss2eng_sourcebus!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool=false)
    # create virtual sourcebus
    circuit = _create_vsource(get(data_dss["circuit"], "bus1", "sourcebus"), data_dss["circuit"]["name"]; _to_sym_keys(data_dss["circuit"])...)

    nodes = Array{Bool}([1 1 1 0])
    ph1_ang = circuit["angle"]
    vm_pu = circuit["pu"]

    phases = circuit["phases"]
    vnom = data_eng["settings"]["set_vbase_val"]

    vm = fill(vm_pu, 3)*vnom
    va = _wrap_to_pi.([-2*pi/phases*(i-1)+deg2rad(ph1_ang) for i in 1:phases])

    add_voltage_source!(data_eng, id="source", bus="sourcebus", connections=collect(1:phases), vm=vm, va=va, rs=circuit["rmatrix"], xs=circuit["xmatrix"])
end


""
function _dss2eng_vsource!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool=false)
end


""
function _dss2eng_loadshape!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool=false)
    for (name, dss_obj) in get(data_dss, "loadshape", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "loadshape")
        defaults = _apply_ordered_properties(_create_loadshape(name; _to_sym_keys(dss_obj)...), dss_obj)

        eng_obj = Dict{String,Any}()

        eng_obj["hour"] = defaults["hour"]
        eng_obj["pmult"] = defaults["pmult"]
        eng_obj["qmult"] = defaults["qmult"]
        eng_obj["use_actual"] = defaults["useactual"]

        if !haskey(data_eng, "loadshape")
            data_eng["loadshape"] = Dict{String,Any}()
        end

        if import_all
            _import_all!(eng_obj, defaults, dss_obj["prop_order"])
        end

        data_eng["loadshape"][name] = eng_obj
    end
end


"Adds loads to `data_eng` from `data_dss`"
function _dss2eng_load!(data_eng::Dict, data_dss::Dict, import_all::Bool, ground_terminal::Int=4)
    for (name, dss_obj) in get(data_dss, "load", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "load")
        defaults = _apply_ordered_properties(_create_load(dss_obj["bus1"], dss_obj["name"]; _to_sym_keys(dss_obj)...), dss_obj)

        # parse the model
        model = defaults["model"]
        # some info on OpenDSS load models
        ##################################
        # Constant can still be scaled by other settings, fixed cannot
        # Note that in the current feature set, fixed therefore equals constant
        # 1: Constant P and Q, default
        if model == 2
        # 2: Constant Z
        elseif model == 3
        # 3: Constant P and quadratic Q
            Memento.warn(_LOGGER, "$load_name: load model 3 not supported. Treating as model 1.")
            model = 1
        elseif model == 4
        # 4: Exponential
            Memento.warn(_LOGGER, "$load_name: load model 4 not supported. Treating as model 1.")
            model = 1
        elseif model == 5
        # 5: Constant I
            #warn(_LOGGER, "$name: load model 5 not supported. Treating as model 1.")
            #model = 1
        elseif model == 6
        # 6: Constant P and fixed Q
            Memento.warn(_LOGGER, "$load_name: load model 6 identical to model 1 in current feature set. Treating as model 1.")
            model = 1
        elseif model == 7
        # 7: Constant P and quadratic Q (i.e., fixed reactance)
            Memento.warn(_LOGGER, "$load_name: load model 7 not supported. Treating as model 1.")
            model = 1
        elseif model == 8
        # 8: ZIP
            Memento.warn(_LOGGER, "$load_name: load model 8 not supported. Treating as model 1.")
            model = 1
        end
        # save adjusted model type to dict, human-readable
        model_int2str = Dict(1=>"constant_power", 2=>"constant_impedance", 5=>"constant_current")
        model = model_int2str[model]

        nphases = defaults["phases"]
        conf = defaults["conn"]


        # connections
        bus = _parse_busname(defaults["bus1"])[1]

        connections_default = conf=="wye" ? [collect(1:nphases)..., 0] : collect(1:nphases)
        connections = _get_conductors_ordered_dm(defaults["bus1"], default=connections_default, check_length=false)
        # if wye connected and neutral not specified, append ground
        if conf=="wye" && length(connections)==nphases
            connections = [connections..., 0]
        end

        # now we can create the load; if you do not have the correct model,
        # pd/qd fields will be populated by default (should not happen for constant current/impedance)
        eng_obj = add_load!(data_eng, id=name, model=model, connections=connections, bus=bus, configuration=conf)

        # if the ground is used directly, register load
        if 0 in connections
            if !haskey(data_eng["bus"][bus], "awaiting_ground")
                data_eng["bus"][bus]["awaiting_ground"] = []
            end
            push!(data_eng["bus"][bus]["awaiting_ground"], eng_obj)
        end

        kv = defaults["kv"]
        if conf=="wye" && nphases in [2, 3]
            kv = kv/sqrt(3)
        end

        if model=="constant_power"
            eng_obj["pd"] = fill(defaults["kw"]/nphases, nphases)
            eng_obj["qd"] = fill(defaults["kvar"]/nphases, nphases)
        else
            eng_obj["pd_ref"] = fill(defaults["kw"]/nphases, nphases)
            eng_obj["qd_ref"] = fill(defaults["kvar"]/nphases, nphases)
            eng_obj["vnom"] = kv
        end

        eng_obj["status"] = convert(Int, defaults["enabled"])
        eng_obj["source_id"] = "load.$name"

        if import_all
            _import_all!(eng_obj, defaults, dss_obj["prop_order"])
        end
    end
end


"Adds capacitors to `data_eng` from `data_dss`"
function _dss2eng_capacitor!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "capacitor", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "capacitor")
        defaults = _apply_ordered_properties(_create_capacitor(dss_obj["bus1"], dss_obj["name"]; _to_sym_keys(dss_obj)...), dss_obj)

        nphases = defaults["phases"]

        dyz_map = Dict("wye"=>"wye", "delta"=>"delta", "ll"=>"delta", "ln"=>"wye")
        conn = dyz_map[defaults["conn"]]

        bus_name = _parse_busname(defaults["bus1"])[1]
        bus2_name = _parse_busname(defaults["bus2"])[1]
        if bus_name!=bus2_name
            Memento.error("Capacitor $(name): bus1 and bus2 should connect to the same bus.")
        end

        f_terminals = _get_conductors_ordered_dm(defaults["bus1"], default=collect(1:nphases))
        if conn=="wye"
            t_terminals = _get_conductors_ordered_dm(defaults["bus2"], default=fill(0,nphases))
        else
            # if delta connected, ignore bus2 and generate t_terminals such that
            # it corresponds to a delta winding
            t_terminals = [f_terminals[2:end]..., f_terminals[1]]
        end

        # 'kv' is specified as phase-to-phase for phases=2/3 (unsure for 4 and more)
        #TODO figure out for more than 3 phases
        vnom_ln = defaults["kv"]
        if defaults["phases"] in [2,3]
            vnom_ln = vnom_ln/sqrt(3)
        end
        # 'kvar' is specified for all phases at once; we want the per-phase one, in MVar
        qnom = (defaults["kvar"]/1E3)/nphases
        b = fill(qnom/vnom_ln^2, nphases)

        # convert to a shunt matrix
        terminals, B = calc_shunt(f_terminals, t_terminals, b)

        # if one terminal is ground (0), reduce shunt addmittance matrix
        terminals, B = _calc_ground_shunt_admittance_matrix(terminals, B, 0)

        eng_obj = add_shunt!(data_eng, id=name, status=convert(Int, defaults["enabled"]), bus=bus_name, connections=terminals, g_sh=fill(0.0, size(B)...), b_sh=B)

        if import_all
            _import_all!(eng_obj, defaults, dss_obj["prop_order"])
        end

        if !haskey(data_eng, "capacitor")
            data_eng["capacitor"] = Dict{String,Any}()
        end

        data_eng["capacitor"][name] = eng_obj
    end
end


"Adds shunt reactors to `data_eng` from `data_dss`"
function _dss2eng_shunt_reactor!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    #TODO revisit this in the future
    # for (name, dss_obj) in get(data_dss, "reactor", Dict{String,Any}())
    #     if !haskey(dss_obj, "bus2")
    #         _apply_like!(dss_obj, data_dss, "reactor")
    #         defaults = _apply_ordered_properties(_create_reactor(dss_obj["bus1"], dss_obj["name"]; _to_sym_keys(dss_obj)...), dss_obj)

    #         eng_obj = Dict{String,Any}()

    #         nphases = defaults["phases"]
    #         name, nodes = _parse_busname(defaults["bus1"])

    #         Zbase = (data_eng["basekv"] / sqrt(3.0))^2 * nphases / data_eng["baseMVA"]  # Use single-phase base impedance for each phase
    #         Gcap = Zbase * sum(defaults["kvar"]) / (nphases * 1e3 * (data_eng["basekv"] / sqrt(3.0))^2)

    #         eng_obj["shunt_bus"] = find_bus(name, data_eng)
    #         eng_obj["name"] = name
    #         eng_obj["gs"] = _parse_array(0.0, nodes, nphases)  # TODO:
    #         eng_obj["bs"] = _parse_array(Gcap, nodes, nconductors)
    #         eng_obj["status"] = convert(Int, defaults["enabled"])
    #         eng_obj["index"] = length(data_eng["shunt"]) + 1

    #         eng_obj["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
    #         eng_obj["source_id"] = "reactor.$(name)"

    #         if import_all
    #             _import_all!(eng_obj, defaults, dss_obj["prop_order"])
    #         end

    #         if !haskey(data_eng, "shunt_reactor")
    #             data_eng["shunt_reactor"] = Dict{String,Any}()
    #         end

    #         data_eng["shunt_reactor"][name] = eng_obj
    #     end
    # end
end


"""
Given a vector and a list of elements to find, this method will return a list
of the positions of the elements in that vector.
"""
function _get_idxs(vec::Array{<:Any, 1}, els::Array{<:Any, 1})
    ret = Array{Int, 1}(undef, length(els))
    for (i,f) in enumerate(els)
        for (j,l) in enumerate(vec)
            if f==l
                ret[i] = j
            end
        end
    end
    return ret
end


"""
Given a set of addmittances 'y' connected from the conductors 'f_cnds' to the
conductors 't_cnds', this method will return a list of conductors 'cnd' and a
matrix 'Y', which will satisfy I[cnds] = Y*V[cnds].
"""
function calc_shunt(f_cnds, t_cnds, y)
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
function _calc_ground_shunt_admittance_matrix(cnds, Y, ground)
    if ground in cnds
        cndsr = setdiff(cnds, ground)
        cndsr_inds = _get_idxs(cnds, cndsr)
        Yr = Y[cndsr_inds, cndsr_inds]
        return (cndsr, Yr)
    else
        return cnds, Y
    end
end


function _rm_floating_cnd(cnds, Y, f)
    P = setdiff(cnds, f)
    f_inds = _get_idxs(cnds, [f])
    P_inds = _get_idxs(cnds, P)
    Yrm = Y[P_inds,P_inds]-(1/Y[f_inds,f_inds][1])*Y[P_inds,f_inds]*Y[f_inds,P_inds]
    return (P,Yrm)
end


"Adds generators to `data_eng` from `data_dss`"
function _dss2eng_generator!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "generator", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "generator")
        defaults = _apply_ordered_properties(_create_generator(dss_obj["bus1"], dss_obj["name"]; _to_sym_keys(dss_obj)...), dss_obj)

        eng_obj = Dict{String,Any}()

        eng_obj["phases"] = defaults["phases"]

        eng_obj["bus"] = _parse_busname(defaults["bus1"])[1]

        eng_obj["pg"] = defaults["kw"]
        eng_obj["qg"] = defaults["kvar"]
        eng_obj["vg"] = defaults["kv"]

        eng_obj["control_model"] = defaults["model"]

        # if PV generator mode convert attached bus to PV bus
        # if eng_obj["control_model"] == 3
        #     data_eng["bus"][eng_obj["bus"]]["bus_type"] = 2
        # end

        eng_obj["model"] = 2
        eng_obj["startup"] = 0.0
        eng_obj["shutdown"] = 0.0
        eng_obj["ncost"] = 3
        eng_obj["cost"] = [0.0, 1.0, 0.0]

        eng_obj["status"] = convert(Int, defaults["enabled"])

        eng_obj["source_id"] = "generator.$(name)"

        if import_all
            _import_all!(eng_obj, defaults, dss_obj["prop_order"])
        end

        if !haskey(data_eng, "generator")
            data_eng["generator"] = Dict{String,Any}()
        end

        data_eng["generator"][name] = eng_obj
    end
end


function _dss2eng_linecode!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "linecode", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "linecode")

        dss_obj["units"] = get(dss_obj, "units", "none")
        dss_obj["circuit_basefreq"] = data_eng["settings"]["basefreq"]

        defaults = _apply_ordered_properties(_create_linecode(name; _to_sym_keys(dss_obj)...), dss_obj)

        eng_obj = Dict{String,Any}()

        nphases = size(defaults["rmatrix"])[1]

        eng_obj["rmatrix"] = reshape(defaults["rmatrix"], nphases, nphases)
        eng_obj["xmatrix"] = reshape(defaults["xmatrix"], nphases, nphases)
        eng_obj["cmatrix"] = reshape(defaults["cmatrix"], nphases, nphases)

        if !haskey(data_eng, "linecode")
            data_eng["linecode"] = Dict{String,Any}()
        end

        data_eng["linecode"][name] = eng_obj
    end
end


"""
    _dss2pmd_line!(data_eng, data_dss, import_all)

Adds PowerModels-style lines to `data_eng` from `data_dss`.
"""
function _dss2eng_line!(data_eng::Dict, data_dss::Dict, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "line", Dict())
        _apply_like!(dss_obj, data_dss, "line")

        if haskey(dss_obj, "basefreq") && dss_obj["basefreq"] != data_eng["settings"]["basefreq"]
            Memento.warn(_LOGGER, "basefreq=$(dss_obj["basefreq"]) on line $(dss_obj["name"]) does not match circuit basefreq=$(data_eng["settings"]["basefreq"])")
            dss_obj["circuit_basefreq"] = data_eng["settings"]["basefreq"]
        end

        defaults = _apply_ordered_properties(_create_line(dss_obj["bus1"], dss_obj["bus2"], dss_obj["name"]; _to_sym_keys(dss_obj)...), dss_obj)

        eng_obj = Dict{String,Any}()

        eng_obj["f_bus"] = _parse_busname(defaults["bus1"])[1]
        eng_obj["t_bus"] = _parse_busname(defaults["bus2"])[1]

        eng_obj["length"] = defaults["length"]

        # TODO nphases not being read correctly from linecode; infer indirectly instead
        nphases = size(defaults["rmatrix"])[1]
        eng_obj["n_conductors"] = nphases

        if haskey(dss_obj, "linecode")
            eng_obj["linecode"] = dss_obj["linecode"]
        end

        for key in ["rmatrix", "xmatrix", "cmatrix"]
            if haskey(dss_obj, key)
                eng_obj[key] = reshape(defaults[key], nphases, nphases)
            end
        end

        eng_obj["f_connections"] = _get_conductors_ordered_dm(defaults["bus1"], default=collect(1:nphases))
        eng_obj["t_connections"] = _get_conductors_ordered_dm(defaults["bus2"], default=collect(1:nphases))


        # TODO calculate these in the conversion/mapping to mathematical model
        # eng_obj["rs"] = rmatrix * defaults["length"]
        # eng_obj["xs"] = xmatrix * defaults["length"]

        # eng_obj["g_fr"] = fill(0.0, nphases, nphases)
        # eng_obj["g_to"] = fill(0.0, nphases, nphases)

        # eng_obj["b_fr"] = (2.0 * pi * data_eng["settings"]["basefreq"] * cmatrix * defaults["length"] / 1e9) / 2.0
        # eng_obj["b_to"] = (2.0 * pi * data_eng["settings"]["basefreq"] * cmatrix * defaults["length"] / 1e9) / 2.0

        eng_obj["status"] = convert(Int, defaults["enabled"])

        eng_obj["source_id"] = "line.$(name)"

        if import_all
            _import_all!(eng_obj, defaults, dss_obj["prop_order"])
        end

        if defaults["switch"]
            if !haskey(data_eng, "switch")
                data_eng["switch"] = Dict{String,Any}()
            end

            data_eng["switch"][name] = eng_obj
        else
            if !haskey(data_eng, "line")
                data_eng["line"] = Dict{String, Any}()
            end

            data_eng["line"][name] = eng_obj
        end
    end
end


"""
    _dss2pmd_transformer!(data_eng, data_dss, import_all)

Adds PMD-style transformers to `data_eng` from `data_dss`.
"""
function _dss2eng_transformer!(data_eng::Dict, data_dss::Dict, import_all::Bool)
   if !haskey(data_eng, "transformer_nw")
        data_eng["transformer_nw"] = Dict{String,Any}()
    end

    for (name, dss_obj) in get(data_dss, "transformer", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "transformer")
        defaults = _apply_ordered_properties(_create_transformer(dss_obj["name"]; _to_sym_keys(dss_obj)...), dss_obj)

        nphases = defaults["phases"]
        nrw = defaults["windings"]

        dyz_map = Dict("wye"=>"wye", "delta"=>"delta", "ll"=>"delta", "ln"=>"wye")
        confs = [dyz_map[x] for x in defaults["conns"]]

        # test if this transformer conforms with limitations
        if nphases<3 && "delta" in confs
            Memento.error("Transformers with delta windings should have at least 3 phases to be well-defined.")
        end
        if nrw>3
            # All of the code is compatible with any number of windings,
            # except for the parsing of the loss model (the pair-wise reactance)
            Memento.error(_LOGGER, "For now parsing of xscarray is not supported. At most 3 windings are allowed, not $nrw.")
        end

        eng_obj = Dict{String, Any}()
        eng_obj["bus"] = Array{String, 1}(undef, nrw)
        eng_obj["connections"] = Array{Array{Int, 1}, 1}(undef, nrw)
        eng_obj["tm"] = Array{Array{Real, 1}, 1}(undef, nrw)
        eng_obj["fixed"] = [[true for i in 1:nphases] for j in 1:nrw]
        eng_obj["vnom"] = [defaults["kvs"][w] for w in 1:nrw]
        eng_obj["snom"] = [defaults["kvas"][w] for w in 1:nrw]
        eng_obj["connections"] = Array{Array{Int, 1}, 1}(undef, nrw)
        eng_obj["configuration"] = Array{String, 1}(undef, nrw)
        eng_obj["polarity"] = Array{Int, 1}(undef, nrw)

        for w in 1:nrw
            eng_obj["bus"][w] = _parse_busname(defaults["buses"][w])[1]

            conn = dyz_map[defaults["conns"][w]]
            eng_obj["configuration"][w] = conn

            terminals_default = conn=="wye" ? [1:nphases..., 0] : collect(1:nphases)
            terminals_w = _get_conductors_ordered_dm(defaults["buses"][w], default=terminals_default)
            eng_obj["connections"][w] = terminals_w
            if 0 in terminals_w
                bus = eng_obj["bus"][w]
                if !haskey(data_eng["bus"][bus], "awaiting_ground")
                    data_eng["bus"][bus]["awaiting_ground"] = []
                end
                push!(data_eng["bus"][bus]["awaiting_ground"], eng_obj)
            end
            eng_obj["polarity"][w] = 1
            eng_obj["tm"][w] = fill(defaults["taps"][w], nphases)
        end

        #eng_obj["source_id"] = "transformer.$(name)"
        if !isempty(defaults["bank"])
            eng_obj["bank"] = defaults["bank"]
        end

        # loss model (converted to SI units, referred to secondary)
        eng_obj["rs"] = [defaults["%rs"][w]/100 for w in 1:nrw]
        eng_obj["noloadloss"] = defaults["%noloadloss"]/100
        eng_obj["imag"] = defaults["%imag"]/100
        if nrw==2
            eng_obj["xsc"] = [defaults["xhl"]]/100
        elseif nrw==3
            eng_obj["xsc"] = [defaults[x] for x in ["xhl", "xht", "xlt"]]/100
        end

        add_virtual!(data_eng, "transformer_nw", create_transformer_nw(;
            Dict(Symbol.(keys(eng_obj)).=>values(eng_obj))...
        ))
    end
end


"""
    _dss2pmd_reactor!(data_eng, data_dss, import_all)

Adds PowerModels-style branch components based on DSS reactors to `data_eng` from `data_dss`
"""
function _dss2eng_reactor!(data_eng::Dict, data_dss::Dict, import_all::Bool)
    if !haskey(data_eng, "branch")
        data_eng["branch"] = []
    end

    if haskey(data_dss, "reactor")
        Memento.warn(_LOGGER, "reactors as constant impedance elements is not yet supported, treating like line")
        for (name, reactor) in data_dss["reactor"]
            if haskey(reactor, "bus2")
                _apply_like!(reactor, data_dss, "reactor")
                defaults = _apply_ordered_properties(_create_reactor(reactor["bus1"], reactor["name"], reactor["bus2"]; _to_sym_keys(reactor)...), reactor)

                eng_obj = Dict{String,Any}()

                nconductors = data_eng["conductors"]

                f_bus, nodes = _parse_busname(defaults["bus1"])
                t_bus = _parse_busname(defaults["bus2"])[1]

                eng_obj["name"] = name
                eng_obj["f_bus"] = find_bus(f_bus, data_eng)
                eng_obj["t_bus"] = find_bus(t_bus, data_eng)

                eng_obj["br_r"] = _PMs.MultiConductorMatrix(_parse_matrix(diagm(0 => fill(0.2, nconductors)), nodes, nconductors))
                eng_obj["br_x"] = _PMs.MultiConductorMatrix(_parse_matrix(zeros(nconductors, nconductors), nodes, nconductors))

                eng_obj["g_fr"] = _parse_array(0.0, nodes, nconductors)
                eng_obj["g_to"] = _parse_array(0.0, nodes, nconductors)
                eng_obj["b_fr"] = _parse_array(0.0, nodes, nconductors)
                eng_obj["b_to"] = _parse_array(0.0, nodes, nconductors)

                for key in ["g_fr", "g_to", "b_fr", "b_to"]
                    eng_obj[key] = _PMs.MultiConductorMatrix(LinearAlgebra.diagm(0=>eng_obj[key].values))
                end

                eng_obj["c_rating_a"] = _parse_array(defaults["normamps"], nodes, nconductors)
                eng_obj["c_rating_b"] = _parse_array(defaults["emergamps"], nodes, nconductors)
                eng_obj["c_rating_c"] = _parse_array(defaults["emergamps"], nodes, nconductors)

                eng_obj["tap"] = _parse_array(1.0, nodes, nconductors, NaN)
                eng_obj["shift"] = _parse_array(0.0, nodes, nconductors)

                eng_obj["br_status"] = convert(Int, defaults["enabled"])

                eng_obj["angmin"] = _parse_array(-60.0, nodes, nconductors, -60.0)
                eng_obj["angmax"] = _parse_array( 60.0, nodes, nconductors,  60.0)

                eng_obj["transformer"] = true

                eng_obj["index"] = length(data_eng["branch"]) + 1

                nodes = .+([_parse_busname(defaults[n])[2] for n in ["bus1", "bus2"]]...)
                eng_obj["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
                eng_obj["source_id"] = "reactor.$(name)"

                if import_all
                    _import_all!(eng_obj, defaults, dss_obj["prop_order"])
                end

                push!(data_eng["branch"], eng_obj)
            end
        end
    end
end


"""
    _dss2pmd_pvsystem!(data_eng, data_dss)

Adds PowerModels-style pvsystems to `data_eng` from `data_dss`.
"""
function _dss2eng_pvsystem!(data_eng::Dict, data_dss::Dict, import_all::Bool)

    for (name, dss_obj) in get(data_dss, "pvsystem", Dict{String,Any}())
        Memento.warn(_LOGGER, "Converting PVSystem \"$(dss_obj["name"])\" into generator with limits determined by OpenDSS property 'kVA'")

        _apply_like!(dss_obj, data_dss, "pvsystem")
        defaults = _apply_ordered_properties(_create_pvsystem(dss_obj["bus1"], dss_obj["name"]; _to_sym_keys(dss_obj)...), dss_obj)

        eng_obj = Dict{String,Any}()

        eng_obj["name"] = name
        eng_obj["bus"] = _parse_busname(defaults["bus1"])[1]

        eng_obj["phases"] = defaults["phases"]

        eng_obj["ppv"] = defaults["kva"]
        eng_obj["qpv"] = defaults["kva"]
        eng_obj["vpv"] = defaults["kv"]

        eng_obj["model"] = 2
        eng_obj["startup"] = 0.0
        eng_obj["shutdown"] = 0.0
        eng_obj["ncost"] = 3
        eng_obj["cost"] = [0.0, 1.0, 0.0]

        eng_obj["status"] = convert(Int, defaults["enabled"])

        eng_obj["source_id"] = "pvsystem.$(name)"

        if import_all
            _import_all!(eng_obj, defaults, dss_obj["prop_order"])
        end

        if !haskey(data_eng, "pvsystem")
            data_eng["pvsystem"] = Dict{String,Any}()
        end

        data_eng["pvsystem"][name] = eng_obj
    end
end


"""
    _dss2pmd_storage!(data_eng, data_dss, import_all)

Adds PowerModels-style storage to `data_eng` from `data_dss`
"""
function _dss2eng_storage!(data_eng::Dict, data_dss::Dict, import_all::Bool)

    for (name, dss_obj) in get(data_dss, "storage", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "storage")
        defaults = _apply_ordered_properties(_create_storage(dss_obj["bus1"], dss_obj["name"]; _to_sym_keys(dss_obj)...), dss_obj)

        eng_obj = Dict{String,Any}()

        eng_obj["phases"] = defaults["phases"]

        eng_obj["name"] = name
        eng_obj["bus"] = _parse_busname(defaults["bus1"])[1]

        eng_obj["energy"] = defaults["kwhstored"]
        eng_obj["energy_rating"] = defaults["kwhrated"]
        eng_obj["charge_rating"] = defaults["%charge"] * defaults["kwrated"] / 100.0
        eng_obj["discharge_rating"] = defaults["%discharge"] * defaults["kwrated"] / 100.0
        eng_obj["charge_efficiency"] = defaults["%effcharge"] / 100.0
        eng_obj["discharge_efficiency"] = defaults["%effdischarge"] / 100.0
        eng_obj["thermal_rating"] = defaults["kva"]
        eng_obj["qmin"] = -defaults["kvar"]
        eng_obj["qmax"] =  defaults["kvar"]
        eng_obj["r"] = defaults["%r"] / 100.0
        eng_obj["x"] = defaults["%x"] / 100.0
        eng_obj["p_loss"] = defaults["%idlingkw"] * defaults["kwrated"]
        eng_obj["q_loss"] = defaults["%idlingkvar"] * defaults["kvar"]

        eng_obj["status"] = convert(Int, defaults["enabled"])

        eng_obj["ps"] = 0.0
        eng_obj["qs"] = 0.0

        eng_obj["source_id"] = "storage.$(name)"

        if import_all
            _import_all(eng_obj, defaults, dss_obj["prop_order"])
        end

        if !haskey(data_eng, "storage")
            data_eng["storage"] = Dict{String,Any}()
        end

        data_eng["storage"][name] = eng_obj
    end
end


"This function appends a component to a component dictionary of a pmd data model"
function _push_dict_ret_key!(dict::Dict{String, Any}, v::Dict{String, Any}; assume_no_gaps=false)
    if isempty(dict)
        k = 1
    elseif assume_no_gaps
        k = length(keys(dict))+1
    else
        k = maximum([parse(Int, x) for x  in keys(dict)])+1
    end

    dict[string(k)] = v
    v["index"] = k
    return k
end


"""
    _where_is_comp(data, comp_id)

Finds existing component of id `comp_id` in array of `data` and returns index.
Assumes all components in `data` are unique.
"""
function _where_is_comp(data::Array, comp_id::AbstractString)::Int
    for (i, e) in enumerate(data)
        if e["name"] == comp_id
            return i
        end
    end
    return 0
end


"""
    _correct_duplicate_components!(data_dss)

Finds duplicate components in `data_dss` and merges up, meaning that older
data (lower indices) is always overwritten by newer data (higher indices).
"""
function _correct_duplicate_components!(data_dss::Dict)
    out = Dict{String,Array}()
    for (k, v) in data_dss
        if !(k in _exclude_duplicate_check)
            out[k] = []
            for comp in v
                if isa(comp, Dict)
                    idx = _where_is_comp(out[k], comp["name"])
                    if idx > 0
                        merge!(out[k][idx], comp)
                    else
                        push!(out[k], comp)
                    end
                end
            end
        end
    end
    merge!(data_dss, out)
end


"Creates a virtual branch between the `virtual_sourcebus` and `sourcebus` with the impedance given by `circuit`"
function _create_sourcebus_vbranch_dm!(data_eng::Dict, circuit::Dict)
    #TODO convert to pu
    rs = circuit["rmatrix"]
    xs = circuit["xmatrix"]

    N = size(rs)[1]

    add_line!(data_eng, id="_virtual_source_imp",
        f_bus="_virtual_sourcebus", t_bus="sourcebus",
        f_connections=collect(1:N), t_connections=collect(1:N),
        rs=rs, xs=xs
    )
    #vbranch = _create_vbranch!(data_eng, sourcebus, vsourcebus; name="sourcebus_vbranch", br_r=br_r, br_x=br_x)
end


"Combines transformers with 'bank' keyword into a single transformer"
function _bank_transformers!(data_eng::Dict)
    transformer_names = Dict(trans["name"] => n for (n, trans) in get(data_eng, "transformer_comp", Dict()))
    bankable_transformers = [trans for trans in values(get(data_eng, "transformer_comp", Dict())) if haskey(trans, "bank")]
    banked_transformers = Dict()
    for transformer in bankable_transformers
        bank = transformer["bank"]

        if !(bank in keys(banked_transformers))
            n = length(data_eng["transformer_comp"])+length(banked_transformers)+1

            banked_transformers[bank] = deepcopy(transformer)
            banked_transformers[bank]["name"] = deepcopy(transformer["bank"])
            banked_transformers[bank]["source_id"] = "transformer.$(transformer["bank"])"
            banked_transformers[bank]["index"] = n
            # set impedances / admittances to zero; only the specified phases should be non-zero
            for key in ["rs", "xs", "bsh", "gsh"]
                inds = key=="xs" ? keys(banked_transformers[bank][key]) : 1:length(banked_transformers[bank][key])
                for w in inds
                    banked_transformers[bank][key][w] *= 0
                end
            end
            delete!(banked_transformers[bank], "bank")
        end

        banked_transformer = banked_transformers[bank]
        for phase in transformer["active_phases"]
            push!(banked_transformer["active_phases"], phase)
            for (k, v) in banked_transformer
                if isa(v, _PMs.MultiConductorVector)
                    banked_transformer[k][phase] = deepcopy(transformer[k][phase])
                elseif isa(v, _PMs.MultiConductorMatrix)
                    banked_transformer[k][phase, :] .= deepcopy(transformer[k][phase, :])
                elseif isa(v, Array) && eltype(v) <: _PMs.MultiConductorVector
                    # most properties are arrays (indexed over the windings)
                    for w in 1:length(v)
                        banked_transformer[k][w][phase] = deepcopy(transformer[k][w][phase])
                    end
                elseif isa(v, Array) && eltype(v) <: _PMs.MultiConductorMatrix
                    # most properties are arrays (indexed over the windings)
                    for w in 1:length(v)
                        banked_transformer[k][w][phase, :] .= deepcopy(transformer[k][w][phase, :])
                    end
                elseif k=="xs"
                    # xs is a Dictionary indexed over pairs of windings
                    for w in keys(v)
                        banked_transformer[k][w][phase, :] .= deepcopy(transformer[k][w][phase, :])
                    end
                end
            end
        end
    end

    for transformer in bankable_transformers
        delete!(data_eng["transformer_comp"], transformer_names[transformer["name"]])
    end

    for transformer in values(banked_transformers)
        data_eng["transformer_comp"]["$(transformer["index"])"] = deepcopy(transformer)
    end
end


"""
    parse_options(options)

Parses options defined with the `set` command in OpenDSS.
"""
function parse_options(options)
    out = Dict{String,Any}()

    for (option, dtype) in _dss_option_dtypes
        if haskey(options, option)
            value = options[option]
            if _isa_array(value)
                out[option] = _parse_array(dtype, value)
            elseif _isa_matrix(value)
                out[option] = _parse_matrix(dtype, value)
            else
                out[option] = parse(dtype, value)
            end
        end
    end
    return out
end


"Parses a Dict resulting from the parsing of a DSS file into a PowerModels usable format"
function parse_opendss_dm(data_dss::Dict{String,<:Any}; import_all::Bool=false, bank_transformers::Bool=true)::Dict{String,Any}
    data_eng = create_data_model()

    parse_dss_with_dtypes!(data_dss, _dss_supported_components)

    data_eng["dss_options"] = parse_options(get(data_dss, "options", Dict{String,Any}()))

    if haskey(data_dss, "circuit")
        circuit = data_dss["circuit"]
        defaults = _create_vsource(get(circuit, "bus1", "sourcebus"), circuit["name"]; _to_sym_keys(circuit)...)

        data_eng["name"] = circuit["name"]
        data_eng["sourcebus"] = defaults["bus1"]
        data_eng["data_model"] = "engineering"

        # TODO rename fields
        data_eng["settings"]["kv_kvar_scalar"] = 1
        data_eng["settings"]["set_vbase_val"] = (defaults["basekv"]/sqrt(3))/data_eng["settings"]["v_var_scalar"]
        data_eng["settings"]["set_vbase_bus"] = data_eng["sourcebus"]
        data_eng["settings"]["set_sbase_val"] = defaults["basemva"]*1e3/data_eng["settings"]["v_var_scalar"]
        data_eng["settings"]["basefreq"] = get(get(data_dss, "options", Dict()), "defaultbasefreq", 60.0)

        data_eng["files"] = data_dss["filename"]
    else
        Memento.error(_LOGGER, "Circuit not defined, not a valid circuit!")
    end

    _dss2eng_bus!(data_eng, data_dss, import_all)
    _dss2eng_buscoords!(data_eng, data_dss)

    _dss2eng_linecode!(data_eng, data_dss, import_all)
    _dss2eng_line!(data_eng, data_dss, import_all)

    # _dss2eng_xfrmcode!(data_eng, data_dss, import_all)
    _dss2eng_transformer!(data_eng, data_dss, import_all)

    #_dss2eng_line_reactor!(data_eng, data_dss, import_all)

    _dss2eng_loadshape!(data_eng, data_dss, import_all)
    _dss2eng_load!(data_eng, data_dss, import_all)

    _dss2eng_capacitor!(data_eng, data_dss, import_all)
    _dss2eng_shunt_reactor!(data_eng, data_dss, import_all)

    _dss2eng_generator!(data_eng, data_dss, import_all)

    _dss2eng_pvsystem!(data_eng, data_dss, import_all)

    _dss2eng_storage!(data_eng, data_dss, import_all)

    if bank_transformers
        _bank_transformers!(data_eng)
    end

    _discover_terminals!(data_eng)


    return data_eng
end


""
function _discover_terminals!(data_eng::Dict{String,<:Any})
    terminals = Dict{String, Set{Int}}([(bus["id"], Set{Int}()) for (_,bus) in data_eng["bus"]])

    for (_,dss_obj) in data_eng["line"]
        # ignore 0 terminal
        push!(terminals[dss_obj["f_bus"]], setdiff(dss_obj["f_connections"], [0])...)
        push!(terminals[dss_obj["t_bus"]], setdiff(dss_obj["t_connections"], [0])...)
    end

    if haskey(data_eng, "transformer_nw")
        for (_,tr) in data_eng["transformer_nw"]
            for w in 1:length(tr["bus"])
                # ignore 0 terminal
                push!(terminals[tr["bus"][w]], setdiff(tr["connections"][w], [0])...)
            end
        end
    end

    for (id, bus) in data_eng["bus"]
        data_eng["bus"][id]["terminals"] = sort(collect(terminals[id]))
    end

    # identify neutrals and propagate along cables
    bus_neutral = _find_neutrals(data_eng)

    for (id,bus) in data_eng["bus"]
        if haskey(bus, "awaiting_ground") || haskey(bus_neutral, id)
            # this bus will need a neutral
            if haskey(bus_neutral, id)
                neutral = bus_neutral[id]
            else
                neutral = maximum(bus["terminals"])+1
                push!(bus["terminals"], neutral)
            end
            bus["neutral"] = neutral
            if haskey(bus, "awaiting_ground")
                bus["grounded"] = [neutral]
                bus["rg"] = [0.0]
                bus["xg"] = [0.0]
                for comp in bus["awaiting_ground"]
                    if eltype(comp["connections"])<:Array
                        for w in 1:length(comp["connections"])
                            if comp["bus"][w]==id
                                comp["connections"][w] .+= (comp["connections"][w].==0)*neutral
                            end
                        end
                    else
                        comp["connections"] .+= (comp["connections"].==0)*neutral
                    end
                    # @show comp["connections"]
                end
                #delete!(bus, "awaiting_ground")
            end
        end
        phases = haskey(bus, "neutral") ? setdiff(bus["terminals"], bus["neutral"]) : bus["terminals"]
        bus["phases"] = phases
    end
end


""
function _find_neutrals(data_eng)
    vertices = [(id, t) for (id, bus) in data_eng["bus"] for t in bus["terminals"]]
    neutrals = []
    edges = Set([((dss_obj["f_bus"], dss_obj["f_connections"][c]),(dss_obj["t_bus"], dss_obj["t_connections"][c])) for (id, dss_obj) in data_eng["line"] for c in 1:length(dss_obj["f_connections"])])

    bus_neutrals = [(id,bus["neutral"]) for (id,bus) in data_eng["bus"] if haskey(bus, "neutral")]
    trans_neutrals = []
    for (_, tr) in data_eng["transformer_nw"]
        for w in 1:length(tr["connections"])
            if tr["configuration"][w] == "wye"
                push!(trans_neutrals, (tr["bus"][w], tr["connections"][w][end]))
            end
        end
    end
    load_neutrals = [(dss_obj["bus"],dss_obj["connections"][end]) for (_,dss_obj) in data_eng["load"] if dss_obj["configuration"]=="wye"]
    neutrals = Set(vcat(bus_neutrals, trans_neutrals, load_neutrals))
    neutrals = Set([(bus,t) for (bus,t) in neutrals if t!=0])
    stack = deepcopy(neutrals)
    while !isempty(stack)
        vertex = pop!(stack)
        candidates_t = [((f,t), t) for (f,t) in edges if f==vertex]
        candidates_f = [((f,t), f) for (f,t) in edges if t==vertex]
        for (edge,next) in [candidates_t..., candidates_f...]
            delete!(edges, edge)
            push!(stack, next)
            push!(neutrals, next)
        end
    end
    bus_neutral = Dict{String, Int}()
    for (bus,t) in neutrals
        bus_neutral[bus] = t
    end
    return bus_neutral
end


"Parses a DSS file into a PowerModels usable format"
function parse_opendss_dm(io::IOStream; import_all::Bool=false, bank_transformers::Bool=true)::Dict
    data_dss = parse_dss(io)

    return parse_opendss_dm(data_dss; import_all=import_all)
end


"Returns an ordered list of defined conductors. If ground=false, will omit any `0`"
function _get_conductors_ordered_dm(busname::AbstractString; default=[], check_length=true)::Array
    parts = split(busname, '.'; limit=2)
    ret = []
    if length(parts)==2
        conds_str = split(parts[2], '.')
        ret = [parse(Int, i) for i in conds_str]
    else
        return default
    end

    if check_length && length(default)!=length(ret)
        Memento.error("An incorrect number of nodes was specified on $(parts[1]); |$(parts[2])|!=$(length(default)).")
    end
    return ret
end


""
function _import_all!(component::Dict{String,<:Any}, defaults::Dict{String,<:Any}, prop_order::Array{<:AbstractString,1})
    component["dss"] = Dict{String,Any}((key, defaults[key]) for key in prop_order)
end
