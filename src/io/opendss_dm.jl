# OpenDSS parser
import LinearAlgebra: isdiag, diag, pinv


"Structure representing OpenDSS `dss_source_id` giving the type of the component `dss_type`, its name `dss_name`, and the active phases `active_phases`"
struct DSSSourceId
    dss_type::AbstractString
    dss_name::AbstractString
    active_phases::Set{Int}
end


"Parses a component's OpenDSS source information into the `dss_source_id` struct"
function _parse_dss_source_id(component::Dict)::DSSSourceId
    dss_type, dss_name = split(component["source_id"], '.')
    return DSSSourceId(dss_type, dss_name, Set(component["active_phases"]))
end


"returns the linecode with name `id`"
function _get_linecode(dss_data::Dict, id::AbstractString)
    if haskey(dss_data, "linecode")
        for item in dss_data["linecode"]
            if item["name"] == id
                return item
            end
        end
    end
    return Dict{String,Any}()
end


"""
    _discover_buses(dss_data)

Discovers all of the buses (not separately defined in OpenDSS), from "lines".
"""
function _discover_buses(dss_data::Dict)::Array
    bus_names = []
    buses = []
    for compType in ["line", "transformer", "reactor"]
        if haskey(dss_data, compType)
            compList = dss_data[compType]
            for compObj in compList
                if compType == "transformer"
                    compObj = _create_transformer(compObj["name"]; _to_sym_keys(compObj)...)
                    for bus in compObj["buses"]
                        name, nodes = _parse_busname(bus)
                        if !(name in bus_names)
                            push!(bus_names, name)
                            push!(buses, (name, nodes))
                        end
                    end
                elseif haskey(compObj, "bus2")
                    for key in ["bus1", "bus2"]
                        name, nodes = _parse_busname(compObj[key])
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
        Memento.error(_LOGGER, "dss_data has no lines!")
    else
        return buses
    end
end


"""
    _dss2pmd_bus!(pmd_data, dss_data)

Adds PowerModels-style buses to `pmd_data` from `dss_data`.
"""
function _dss2pmd_bus_dm!(pmd_data::Dict, dss_data::Dict, import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1)

    buses = _discover_buses(dss_data)
    for (n, (bus, nodes)) in enumerate(buses)

        @assert(!(length(bus)>=8 && bus[1:8]=="_virtual"), "Bus $bus: identifiers should not start with _virtual.")

        add_bus!(pmd_data, id=bus, status=1, bus_type=1)
    end

    # create virtual sourcebus
    circuit = _create_vsource(get(dss_data["circuit"][1], "bus1", "sourcebus"), dss_data["circuit"][1]["name"]; _to_sym_keys(dss_data["circuit"][1])...)

    nodes = Array{Bool}([1 1 1 0])
    ph1_ang = circuit["angle"]
    vm_pu = circuit["pu"]
    vmi = circuit["pu"] - circuit["pu"] / (circuit["mvasc3"] / circuit["basemva"])
    vma = circuit["pu"] + circuit["pu"] / (circuit["mvasc3"] / circuit["basemva"])

    phases = circuit["phases"]
    vnom = pmd_data["settings"]["set_vbase_val"]

    vm = fill(vm_pu, 3)*vnom
    va = _wrap_to_pi.([-2*pi/phases*(i-1)+deg2rad(ph1_ang) for i in 1:phases])

    add_voltage_source!(pmd_data, id="source", bus="sourcebus", connections=collect(1:phases),
        vm=vm, va=va,
        rs=circuit["rmatrix"], xs=circuit["xmatrix"],
    )
end


"""
    find_component(pmd_data, name, compType)

Returns the component of `compType` with `name` from `data` of type
Dict{String,Array}.
"""
function find_component(data::Dict, name::AbstractString, compType::AbstractString)::Dict
    for comp in values(data[compType])
        if comp["name"] == name
            return comp
        end
    end
    Memento.warn(_LOGGER, "Could not find $compType \"$name\"")
    return Dict{String,Any}()
end


"""
    find_bus(busname, pmd_data)

Finds the index number of the bus in existing data from the given `busname`.
"""
function find_bus(busname::AbstractString, pmd_data::Dict)
    bus = find_component(pmd_data, busname, "bus")
    if haskey(bus, "bus_i")
        return bus["bus_i"]
    else
        Memento.error(_LOGGER, "cannot find connected bus with id \"$busname\"")
    end
end


"""
    _dss2pmd_load!(pmd_data, dss_data, import_all)

Adds PowerModels-style loads to `pmd_data` from `dss_data`.
"""
function _dss2pmd_load_dm!(pmd_data::Dict, dss_data::Dict, import_all::Bool, ground_terminal::Int=4)

    for load in get(dss_data, "load", [])
        _apply_like!(load, dss_data, "load")
        defaults = _apply_ordered_properties(_create_load(load["bus1"], load["name"]; _to_sym_keys(load)...), load)

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
        loadDict = add_load!(pmd_data, id=defaults["name"], model=model, connections=connections, bus=bus, configuration=conf)

        # if the ground is used directly, register load
        if 0 in connections
            if !haskey(pmd_data["bus"][bus], "awaiting_ground")
                pmd_data["bus"][bus]["awaiting_ground"] = []
            end
            push!(pmd_data["bus"][bus]["awaiting_ground"], loadDict)
        end

        kv = defaults["kv"]
        if conf=="wye" && nphases in [2, 3]
            kv = kv/sqrt(3)
        end

        if model=="constant_power"
            loadDict["pd"] = fill(defaults["kw"]/nphases, nphases)
            loadDict["qd"] = fill(defaults["kvar"]/nphases, nphases)
        else
            loadDict["pd_ref"] = fill(defaults["kw"]/nphases, nphases)
            loadDict["qd_ref"] = fill(defaults["kvar"]/nphases, nphases)
            loadDict["vnom"] = kv
        end

        #loadDict["status"] = convert(Int, defaults["enabled"])

        #loadDict["source_id"] = "load.$load_name"

        used = ["phases", "bus1", "name"]
        _PMs._import_remaining!(loadDict, defaults, import_all; exclude=used)
    end
end


"""
    _dss2pmd_shunt!(pmd_data, dss_data, import_all)

Adds PowerModels-style shunts to `pmd_data` from `dss_data`.
"""
function _dss2pmd_shunt_dm!(pmd_data::Dict, dss_data::Dict, import_all::Bool)

    for shunt in get(dss_data, "capacitor", [])
        _apply_like!(shunt, dss_data, "capacitor")
        defaults = _apply_ordered_properties(_create_capacitor(shunt["bus1"], shunt["name"]; _to_sym_keys(shunt)...), shunt)

        nphases = defaults["phases"]

        dyz_map = Dict("wye"=>"wye", "delta"=>"delta", "ll"=>"delta", "ln"=>"wye")
        conn = dyz_map[defaults["conn"]]

        bus_name = _parse_busname(defaults["bus1"])[1]
        bus2_name = _parse_busname(defaults["bus2"])[1]
        if bus_name!=bus2_name
            Memento.error("Capacitor $(defaults["name"]): bus1 and bus2 should connect to the same bus.")
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
        terminals, B = ground_shunt(terminals, B, 0)

        shuntDict = add_shunt!(pmd_data, id=defaults["name"], status=convert(Int, defaults["enabled"]),
            bus=bus_name, connections=terminals,
            g_sh=fill(0.0, size(B)...), b_sh=B
        )

        used = ["bus1", "phases", "name"]
        _PMs._import_remaining!(shuntDict, defaults, import_all; exclude=used)

        pmd_data["shunt"][shuntDict["id"]] = shuntDict
    end

    #TODO revisit this in the future
    # for shunt in get(dss_data, "reactor", [])
    #     if !haskey(shunt, "bus2")
    #         _apply_like!(shunt, dss_data, "reactor")
    #         defaults = _apply_ordered_properties(_create_reactor(shunt["bus1"], shunt["name"]; _to_sym_keys(shunt)...), shunt)
    #
    #         shuntDict = Dict{String,Any}()
    #
    #         nconductors = pmd_data["conductors"]
    #         name, nodes = _parse_busname(defaults["bus1"])
    #
    #         Zbase = (pmd_data["basekv"] / sqrt(3.0))^2 * nconductors / pmd_data["baseMVA"]  # Use single-phase base impedance for each phase
    #         Gcap = Zbase * sum(defaults["kvar"]) / (nconductors * 1e3 * (pmd_data["basekv"] / sqrt(3.0))^2)
    #
    #         shuntDict["shunt_bus"] = find_bus(name, pmd_data)
    #         shuntDict["name"] = defaults["name"]
    #         shuntDict["gs"] = _PMs.MultiConductorVector(_parse_array(0.0, nodes, nconductors))  # TODO:
    #         shuntDict["bs"] = _PMs.MultiConductorVector(_parse_array(Gcap, nodes, nconductors))
    #         shuntDict["status"] = convert(Int, defaults["enabled"])
    #         shuntDict["index"] = length(pmd_data["shunt"]) + 1
    #
    #         shuntDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
    #         shuntDict["source_id"] = "reactor.$(defaults["name"])"
    #
    #         used = ["bus1", "phases", "name"]
    #         _PMs._import_remaining!(shuntDict, defaults, import_all; exclude=used)
    #
    #         push!(pmd_data["shunt"], shuntDict)
    #     end
    # end
end


"""
Given a vector and a list of elements to find, this method will return a list
of the positions of the elements in that vector.
"""
function get_inds(vec::Array{<:Any, 1}, els::Array{<:Any, 1})
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
function ground_shunt(cnds, Y, ground)
    if ground in cnds
        cndsr = setdiff(cnds, ground)
        cndsr_inds = get_inds(cnds, cndsr)
        Yr = Y[cndsr_inds, cndsr_inds]
        return (cndsr, Yr)
    else
        return cnds, Y
    end
end


function rm_floating_cnd(cnds, Y, f)
    P = setdiff(cnds, f)
    f_inds = get_inds(cnds, [f])
    P_inds = get_inds(cnds, P)
    Yrm = Y[P_inds,P_inds]-(1/Y[f_inds,f_inds][1])*Y[P_inds,f_inds]*Y[f_inds,P_inds]
    return (P,Yrm)
end


"""
    _dss2pmd_gen!(pmd_data, dss_data, import_all)

Adds PowerModels-style generators to `pmd_data` from `dss_data`.
"""
function _dss2pmd_gen_dm!(pmd_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(pmd_data, "gen")
        pmd_data["gen"] = []
    end

    # # sourcebus generator (created by circuit)
    # circuit = dss_data["circuit"][1]
    # defaults = _create_vsource(get(circuit, "bus1", "sourcebus"), "sourcegen"; _to_sym_keys(circuit)...)
    #
    # genDict = Dict{String,Any}()
    #
    # nconductors = pmd_data["conductors"]
    # name, nodes = _parse_busname(defaults["bus1"])
    #
    # genDict["gen_bus"] = find_bus("virtual_sourcebus", pmd_data)
    # genDict["name"] = defaults["name"]
    # genDict["gen_status"] = convert(Int, defaults["enabled"])
    #
    # # TODO: populate with VSOURCE properties
    # genDict["pg"] = _PMs.MultiConductorVector(_parse_array( 0.0, nodes, nconductors))
    # genDict["qg"] = _PMs.MultiConductorVector(_parse_array( 0.0, nodes, nconductors))
    #
    # genDict["qmin"] = _PMs.MultiConductorVector(_parse_array(-NaN, nodes, nconductors))
    # genDict["qmax"] = _PMs.MultiConductorVector(_parse_array( NaN, nodes, nconductors))
    #
    # genDict["pmin"] = _PMs.MultiConductorVector(_parse_array(-NaN, nodes, nconductors))
    # genDict["pmax"] = _PMs.MultiConductorVector(_parse_array( NaN, nodes, nconductors))
    #
    # genDict["model"] = 2
    # genDict["startup"] = 0.0
    # genDict["shutdown"] = 0.0
    # genDict["ncost"] = 3
    # genDict["cost"] = [0.0, 1.0, 0.0]
    #
    # genDict["index"] = length(pmd_data["gen"]) + 1
    #
    # genDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
    # genDict["source_id"] = "vsource.$(defaults["name"])"
    #
    # used = ["name", "phases", "bus1"]
    # _PMs._import_remaining!(genDict, defaults, import_all; exclude=used)
    #
    # push!(pmd_data["gen"], genDict)


    for gen in get(dss_data, "generator", [])
        _apply_like!(gen, dss_data, "generator")
        defaults = _apply_ordered_properties(_create_generator(gen["bus1"], gen["name"]; _to_sym_keys(gen)...), gen)

        genDict = Dict{String,Any}()

        nconductors = pmd_data["conductors"]
        name, nodes = _parse_busname(defaults["bus1"])

        genDict["gen_bus"] = find_bus(name, pmd_data)
        genDict["name"] = defaults["name"]
        genDict["gen_status"] = convert(Int, defaults["enabled"])
        genDict["pg"] = _PMs.MultiConductorVector(_parse_array(defaults["kw"] / (1e3 * nconductors), nodes, nconductors))
        genDict["qg"] = _PMs.MultiConductorVector(_parse_array(defaults["kvar"] / (1e3 * nconductors), nodes, nconductors))
        genDict["vg"] = _PMs.MultiConductorVector(_parse_array(defaults["kv"] / pmd_data["basekv"], nodes, nconductors))

        genDict["qmin"] = _PMs.MultiConductorVector(_parse_array(defaults["minkvar"] / (1e3 * nconductors), nodes, nconductors))
        genDict["qmax"] = _PMs.MultiConductorVector(_parse_array(defaults["maxkvar"] / (1e3 * nconductors), nodes, nconductors))

        genDict["apf"] = _PMs.MultiConductorVector(_parse_array(0.0, nodes, nconductors))

        genDict["pmax"] = genDict["pg"]  # Assumes generator is at rated power
        genDict["pmin"] = 0.0 * genDict["pg"]  # 0% of pmax

        genDict["pc1"] = genDict["pmax"]
        genDict["pc2"] = genDict["pmin"]
        genDict["qc1min"] = genDict["qmin"]
        genDict["qc1max"] = genDict["qmax"]
        genDict["qc2min"] = genDict["qmin"]
        genDict["qc2max"] = genDict["qmax"]

        # For distributed generation ramp rates are not usually an issue
        # and they are not supported in OpenDSS
        genDict["ramp_agc"] = genDict["pmax"]

        genDict["ramp_q"] = _PMs.MultiConductorVector(_parse_array(max.(abs.(genDict["qmin"].values), abs.(genDict["qmax"].values)), nodes, nconductors))
        genDict["ramp_10"] = genDict["pmax"]
        genDict["ramp_30"] = genDict["pmax"]

        genDict["control_model"] = defaults["model"]

        # if PV generator mode convert attached bus to PV bus
        if genDict["control_model"] == 3
            pmd_data["bus"][genDict["gen_bus"]]["bus_type"] = 2
        end

        genDict["model"] = 2
        genDict["startup"] = 0.0
        genDict["shutdown"] = 0.0
        genDict["ncost"] = 3
        genDict["cost"] = [0.0, 1.0, 0.0]

        genDict["index"] = length(pmd_data["gen"]) + 1

        genDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
        genDict["source_id"] = "generator.$(defaults["name"])"

        used = ["name", "phases", "bus1"]
        _PMs._import_remaining!(genDict, defaults, import_all; exclude=used)

        push!(pmd_data["gen"], genDict)
    end

    for pv in get(dss_data, "pvsystem", [])
        Memento.warn(_LOGGER, "Converting PVSystem \"$(pv["name"])\" into generator with limits determined by OpenDSS property 'kVA'")

        _apply_like!(pv, dss_data, "pvsystem")
        defaults = _apply_ordered_properties(_create_pvsystem(pv["bus1"], pv["name"]; _to_sym_keys(pv)...), pv)

        pvDict = Dict{String,Any}()

        nconductors = pmd_data["conductors"]
        name, nodes = _parse_busname(defaults["bus1"])

        pvDict["name"] = defaults["name"]
        pvDict["gen_bus"] = find_bus(name, pmd_data)

        pvDict["pg"] = _PMs.MultiConductorVector(_parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors))
        pvDict["qg"] = _PMs.MultiConductorVector(_parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors))
        pvDict["vg"] = _PMs.MultiConductorVector(_parse_array(defaults["kv"] / pmd_data["basekv"], nodes, nconductors))

        pvDict["pmin"] = _PMs.MultiConductorVector(_parse_array(0.0, nodes, nconductors))
        pvDict["pmax"] = _PMs.MultiConductorVector(_parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors))

        pvDict["qmin"] = -_PMs.MultiConductorVector(_parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors))
        pvDict["qmax"] =  _PMs.MultiConductorVector(_parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors))

        pvDict["gen_status"] = convert(Int, defaults["enabled"])

        pvDict["model"] = 2
        pvDict["startup"] = 0.0
        pvDict["shutdown"] = 0.0
        pvDict["ncost"] = 3
        pvDict["cost"] = [0.0, 1.0, 0.0]

        pvDict["index"] = length(pmd_data["gen"]) + 1

        pvDict["active_phases"] = [nodes[n] > 0 ? 1 : 0 for n in 1:nconductors]
        pvDict["source_id"] = "pvsystem.$(defaults["name"])"

        used = ["name", "phases", "bus1"]
        _PMs._import_remaining!(pvDict, defaults, import_all; exclude=used)

        push!(pmd_data["gen"], pvDict)
    end
end


"""
    _dss2pmd_line!(pmd_data, dss_data, import_all)

Adds PowerModels-style lines to `pmd_data` from `dss_data`.
"""
function _dss2pmd_line_dm!(pmd_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(pmd_data, "branch")
        pmd_data["line"] = Dict{String, Any}()
    end

    #nconductors = pmd_data["conductors"]

    for line in get(dss_data, "line", [])
        _apply_like!(line, dss_data, "line")

        if haskey(line, "linecode")
            linecode = deepcopy(_get_linecode(dss_data, get(line, "linecode", "")))
            if haskey(linecode, "like")
                linecode = merge(find_component(dss_data, linecode["like"], "linecode"), linecode)
            end

            linecode["units"] = get(line, "units", "none") == "none" ? "none" : get(linecode, "units", "none")
            linecode["circuit_basefreq"] = pmd_data["settings"]["basefreq"]

            linecode = _create_linecode(get(linecode, "name", ""); _to_sym_keys(linecode)...)
            delete!(linecode, "name")
        else
            linecode = Dict{String,Any}()
        end

        if haskey(line, "basefreq") && line["basefreq"] != pmd_data["settings"]["basefreq"]
            Memento.warn(_LOGGER, "basefreq=$(line["basefreq"]) on line $(line["name"]) does not match circuit basefreq=$(pmd_data["settings"]["basefreq"])")
            line["circuit_basefreq"] = pmd_data["settings"]["basefreq"]
        end

        defaults = _apply_ordered_properties(_create_line(line["bus1"], line["bus2"], line["name"]; _to_sym_keys(line)...), line; code_dict=linecode)

        lineDict = Dict{String,Any}()

        lineDict["id"] = defaults["name"]
        lineDict["f_bus"] = _parse_busname(defaults["bus1"])[1]
        lineDict["t_bus"] = _parse_busname(defaults["bus2"])[1]

        #lineDict["length"] = defaults["length"]

        #TODO nphases not being read correctly from linecode; infer indirectly instead
        nphases = size(defaults["rmatrix"])[1]
        lineDict["n_conductors"] = nphases

        #TODO fix this in a cleaner way
        # ensure that this actually is a matrix and not a vector for 1x1 data
        rmatrix = reshape(defaults["rmatrix"], nphases, nphases)
        xmatrix = reshape(defaults["xmatrix"], nphases, nphases)
        cmatrix = reshape(defaults["cmatrix"], nphases, nphases)

        lineDict["f_connections"] = _get_conductors_ordered_dm(defaults["bus1"], default=collect(1:nphases))
        lineDict["t_connections"] = _get_conductors_ordered_dm(defaults["bus2"], default=collect(1:nphases))

        lineDict["rs"] = rmatrix * defaults["length"]
        lineDict["xs"] = xmatrix * defaults["length"]

        lineDict["g_fr"] = fill(0.0, nphases, nphases)
        lineDict["g_to"] = fill(0.0, nphases, nphases)

        lineDict["b_fr"] = (2.0 * pi * pmd_data["settings"]["basefreq"] * cmatrix * defaults["length"] / 1e9) / 2.0
        lineDict["b_to"] = (2.0 * pi * pmd_data["settings"]["basefreq"] * cmatrix * defaults["length"] / 1e9) / 2.0

        #lineDict["c_rating_a"] = fill(defaults["normamps"], nphases)
        #lineDict["c_rating_b"] = fill(defaults["emergamps"], nphases)
        #lineDict["c_rating_c"] = fill(defaults["emergamps"], nphases)

        lineDict["status"] = convert(Int, defaults["enabled"])

        #lineDict["source_id"] = "line.$(defaults["name"])"

        #used = ["name", "bus1", "bus2", "rmatrix", "xmatrix"]
        #_PMs._import_remaining!(lineDict, defaults, import_all; exclude=used)

        pmd_data["line"][lineDict["id"]] = lineDict
    end
end


"""
    _dss2pmd_transformer!(pmd_data, dss_data, import_all)

Adds ThreePhasePowerModels-style transformers to `pmd_data` from `dss_data`.
"""
function _dss2pmd_transformer_dm!(pmd_data::Dict, dss_data::Dict, import_all::Bool)
   if !haskey(pmd_data, "transformer_nw")
        pmd_data["transformer_nw"] = Dict{String,Any}()
    end

    for transformer in get(dss_data, "transformer", [])
        _apply_like!(transformer, dss_data, "transformer")
        defaults = _apply_ordered_properties(_create_transformer(transformer["name"]; _to_sym_keys(transformer)...), transformer)

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

        transDict = Dict{String, Any}()
        transDict["id"] = defaults["name"]
        transDict["bus"] = Array{String, 1}(undef, nrw)
        transDict["connections"] = Array{Array{Int, 1}, 1}(undef, nrw)
        transDict["tm"] = Array{Array{Real, 1}, 1}(undef, nrw)
        transDict["fixed"] = [[true for i in 1:nphases] for j in 1:nrw]
        transDict["vnom"] = [defaults["kvs"][w] for w in 1:nrw]
        transDict["snom"] = [defaults["kvas"][w] for w in 1:nrw]
        transDict["connections"] = Array{Array{Int, 1}, 1}(undef, nrw)
        transDict["configuration"] = Array{String, 1}(undef, nrw)
        transDict["polarity"] = Array{Int, 1}(undef, nrw)

        for w in 1:nrw
            transDict["bus"][w] = _parse_busname(defaults["buses"][w])[1]

            conn = dyz_map[defaults["conns"][w]]
            transDict["configuration"][w] = conn

            terminals_default = conn=="wye" ? [1:nphases..., 0] : collect(1:nphases)
            terminals_w = _get_conductors_ordered_dm(defaults["buses"][w], default=terminals_default)
            transDict["connections"][w] = terminals_w
            if 0 in terminals_w
                bus = transDict["bus"][w]
                if !haskey(pmd_data["bus"][bus], "awaiting_ground")
                    pmd_data["bus"][bus]["awaiting_ground"] = []
                end
                push!(pmd_data["bus"][bus]["awaiting_ground"], transDict)
            end
            transDict["polarity"][w] = 1
            transDict["tm"][w] = fill(defaults["taps"][w], nphases)
        end

        #transDict["source_id"] = "transformer.$(defaults["name"])"
        if !isempty(defaults["bank"])
            transDict["bank"] = defaults["bank"]
        end

        # loss model (converted to SI units, referred to secondary)
        transDict["rs"] = [defaults["%rs"][w]/100 for w in 1:nrw]
        transDict["noloadloss"] = defaults["%noloadloss"]/100
        transDict["imag"] = defaults["%imag"]/100
        if nrw==2
            transDict["xsc"] = [defaults["xhl"]]/100
        elseif nrw==3
            transDict["xsc"] = [defaults[x] for x in ["xhl", "xht", "xlt"]]/100
        end

        add_virtual!(pmd_data, "transformer_nw", create_transformer_nw(;
            Dict(Symbol.(keys(transDict)).=>values(transDict))...
        ))
    end
end


"""
    _dss2pmd_reactor!(pmd_data, dss_data, import_all)

Adds PowerModels-style branch components based on DSS reactors to `pmd_data` from `dss_data`
"""
function _dss2pmd_reactor_dm!(pmd_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(pmd_data, "branch")
        pmd_data["branch"] = []
    end

    if haskey(dss_data, "reactor")
        Memento.warn(_LOGGER, "reactors as constant impedance elements is not yet supported, treating like line")
        for reactor in dss_data["reactor"]
            if haskey(reactor, "bus2")
                _apply_like!(reactor, dss_data, "reactor")
                defaults = _apply_ordered_properties(_create_reactor(reactor["bus1"], reactor["name"], reactor["bus2"]; _to_sym_keys(reactor)...), reactor)

                reactDict = Dict{String,Any}()

                nconductors = pmd_data["conductors"]

                f_bus, nodes = _parse_busname(defaults["bus1"])
                t_bus = _parse_busname(defaults["bus2"])[1]

                reactDict["name"] = defaults["name"]
                reactDict["f_bus"] = find_bus(f_bus, pmd_data)
                reactDict["t_bus"] = find_bus(t_bus, pmd_data)

                reactDict["br_r"] = _PMs.MultiConductorMatrix(_parse_matrix(diagm(0 => fill(0.2, nconductors)), nodes, nconductors))
                reactDict["br_x"] = _PMs.MultiConductorMatrix(_parse_matrix(zeros(nconductors, nconductors), nodes, nconductors))

                reactDict["g_fr"] = _PMs.MultiConductorVector(_parse_array(0.0, nodes, nconductors))
                reactDict["g_to"] = _PMs.MultiConductorVector(_parse_array(0.0, nodes, nconductors))
                reactDict["b_fr"] = _PMs.MultiConductorVector(_parse_array(0.0, nodes, nconductors))
                reactDict["b_to"] = _PMs.MultiConductorVector(_parse_array(0.0, nodes, nconductors))

                for key in ["g_fr", "g_to", "b_fr", "b_to"]
                    reactDict[key] = _PMs.MultiConductorMatrix(LinearAlgebra.diagm(0=>reactDict[key].values))
                end

                reactDict["c_rating_a"] = _PMs.MultiConductorVector(_parse_array(defaults["normamps"], nodes, nconductors))
                reactDict["c_rating_b"] = _PMs.MultiConductorVector(_parse_array(defaults["emergamps"], nodes, nconductors))
                reactDict["c_rating_c"] = _PMs.MultiConductorVector(_parse_array(defaults["emergamps"], nodes, nconductors))

                reactDict["tap"] = _PMs.MultiConductorVector(_parse_array(1.0, nodes, nconductors, NaN))
                reactDict["shift"] = _PMs.MultiConductorVector(_parse_array(0.0, nodes, nconductors))

                reactDict["br_status"] = convert(Int, defaults["enabled"])

                reactDict["angmin"] = _PMs.MultiConductorVector(_parse_array(-60.0, nodes, nconductors, -60.0))
                reactDict["angmax"] = _PMs.MultiConductorVector(_parse_array( 60.0, nodes, nconductors,  60.0))

                reactDict["transformer"] = true

                reactDict["index"] = length(pmd_data["branch"]) + 1

                nodes = .+([_parse_busname(defaults[n])[2] for n in ["bus1", "bus2"]]...)
                reactDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
                reactDict["source_id"] = "reactor.$(defaults["name"])"

                used = []
                _PMs._import_remaining!(reactDict, defaults, import_all; exclude=used)

                push!(pmd_data["branch"], reactDict)
            end
        end
    end
end


"""
    _dss2pmd_pvsystem!(pmd_data, dss_data)

Adds PowerModels-style pvsystems to `pmd_data` from `dss_data`.
"""
function _dss2pmd_pvsystem_dm!(pmd_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(pmd_data, "pvsystem")
        pmd_data["pvsystem"] = []
    end

    for pvsystem in get(dss_data, "pvsystem", [])
        _apply_like!(pvsystem, dss_data, "pvsystem")
        defaults = _apply_ordered_properties(_create_pvsystem(pvsystem["bus1"], pvsystem["name"]; _to_sym_keys(pvsystem)...), pvsystem)

        pvsystemDict = Dict{String,Any}()

        nconductors = pmd_data["conductors"]
        name, nodes = _parse_busname(defaults["bus1"])

        pvsystemDict["name"] = defaults["name"]
        pvsystemDict["pv_bus"] = find_bus(name, pmd_data)
        pvsystemDict["p"] = _PMs.MultiConductorVector(_parse_array(defaults["kw"] / 1e3, nodes, nconductors))
        pvsystemDict["q"] = _PMs.MultiConductorVector(_parse_array(defaults["kvar"] / 1e3, nodes, nconductors))
        pvsystemDict["status"] = convert(Int, defaults["enabled"])

        pvsystemDict["index"] = length(pmd_data["pvsystem"]) + 1

        pvsystemDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
        pvsystemDict["source_id"] = "pvsystem.$(defaults["name"])"

        used = ["phases", "bus1", "name"]
        _PMs._import_remaining!(pvsystemDict, defaults, import_all; exclude=used)

        push!(pmd_data["pvsystem"], pvsystemDict)
    end
end


"""
    _dss2pmd_storage!(pmd_data, dss_data, import_all)

Adds PowerModels-style storage to `pmd_data` from `dss_data`
"""
function _dss2pmd_storage_dm!(pmd_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(pmd_data, "storage")
        pmd_data["storage"] = []
    end

    for storage in get(dss_data, "storage", [])
        _apply_like!(storage, dss_data, "storage")
        defaults = _apply_ordered_properties(_create_storage(storage["bus1"], storage["name"]; _to_sym_keys(storage)...), storage)

        storageDict = Dict{String,Any}()

        nconductors = pmd_data["conductors"]
        name, nodes = _parse_busname(defaults["bus1"])

        storageDict["name"] = defaults["name"]
        storageDict["storage_bus"] = find_bus(name, pmd_data)
        storageDict["energy"] = defaults["kwhstored"] / 1e3
        storageDict["energy_rating"] = defaults["kwhrated"] / 1e3
        storageDict["charge_rating"] = defaults["%charge"] * defaults["kwrated"] / 1e3 / 100.0
        storageDict["discharge_rating"] = defaults["%discharge"] * defaults["kwrated"] / 1e3 / 100.0
        storageDict["charge_efficiency"] = defaults["%effcharge"] / 100.0
        storageDict["discharge_efficiency"] = defaults["%effdischarge"] / 100.0
        storageDict["thermal_rating"] = _PMs.MultiConductorVector(_parse_array(defaults["kva"] / 1e3 / nconductors, nodes, nconductors))
        storageDict["qmin"] = _PMs.MultiConductorVector(_parse_array(-defaults["kvar"] / 1e3 / nconductors, nodes, nconductors))
        storageDict["qmax"] = _PMs.MultiConductorVector(_parse_array( defaults["kvar"] / 1e3 / nconductors, nodes, nconductors))
        storageDict["r"] = _PMs.MultiConductorVector(_parse_array(defaults["%r"] / 100.0, nodes, nconductors))
        storageDict["x"] = _PMs.MultiConductorVector(_parse_array(defaults["%x"] / 100.0, nodes, nconductors))
        storageDict["p_loss"] = defaults["%idlingkw"] * defaults["kwrated"] / 1e3
        storageDict["q_loss"] = defaults["%idlingkvar"] * defaults["kvar"] / 1e3

        storageDict["status"] = convert(Int, defaults["enabled"])

        storageDict["ps"] = _PMs.MultiConductorVector(_parse_array(0.0, nodes, nconductors))
        storageDict["qs"] = _PMs.MultiConductorVector(_parse_array(0.0, nodes, nconductors))

        storageDict["index"] = length(pmd_data["storage"]) + 1

        storageDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
        storageDict["source_id"] = "storage.$(defaults["name"])"

        used = ["phases", "bus1", "name"]
        _PMs._import_remaining!(storageDict, defaults, import_all; exclude=used)

        push!(pmd_data["storage"], storageDict)
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
    _correct_duplicate_components!(dss_data)

Finds duplicate components in `dss_data` and merges up, meaning that older
data (lower indices) is always overwritten by newer data (higher indices).
"""
function _correct_duplicate_components!(dss_data::Dict)
    out = Dict{String,Array}()
    for (k, v) in dss_data
        if !(k in ["options"])
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
    merge!(dss_data, out)
end


"Creates a virtual branch between the `virtual_sourcebus` and `sourcebus` with the impedance given by `circuit`"
function _create_sourcebus_vbranch_dm!(pmd_data::Dict, circuit::Dict)
    #TODO convert to pu
    rs = circuit["rmatrix"]
    xs = circuit["xmatrix"]

    N = size(rs)[1]

    add_line!(pmd_data, id="_virtual_source_imp",
        f_bus="_virtual_sourcebus", t_bus="sourcebus",
        f_connections=collect(1:N), t_connections=collect(1:N),
        rs=rs, xs=xs
    )
    #vbranch = _create_vbranch!(pmd_data, sourcebus, vsourcebus; name="sourcebus_vbranch", br_r=br_r, br_x=br_x)
end


"Combines transformers with 'bank' keyword into a single transformer"
function _bank_transformers!(pmd_data::Dict)
    transformer_names = Dict(trans["name"] => n for (n, trans) in get(pmd_data, "transformer_comp", Dict()))
    bankable_transformers = [trans for trans in values(get(pmd_data, "transformer_comp", Dict())) if haskey(trans, "bank")]
    banked_transformers = Dict()
    for transformer in bankable_transformers
        bank = transformer["bank"]

        if !(bank in keys(banked_transformers))
            n = length(pmd_data["transformer_comp"])+length(banked_transformers)+1

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
        delete!(pmd_data["transformer_comp"], transformer_names[transformer["name"]])
    end

    for transformer in values(banked_transformers)
        pmd_data["transformer_comp"]["$(transformer["index"])"] = deepcopy(transformer)
    end
end


"""
    parse_options(options)

Parses options defined with the `set` command in OpenDSS.
"""
function parse_options(options)
    out = Dict{String,Any}()
    if haskey(options, "voltagebases")
        out["voltagebases"] = _parse_array(Float64, options["voltagebases"])
    end

    if !haskey(options, "defaultbasefreq")
        Memento.warn(_LOGGER, "defaultbasefreq is not defined, default for circuit set to 60 Hz")
        out["defaultbasefreq"] = 60.0
    else
        out["defaultbasefreq"] = parse(Float64, options["defaultbasefreq"])
    end

    return out
end


"Parses a Dict resulting from the parsing of a DSS file into a PowerModels usable format"
function parse_opendss_dm(dss_data::Dict; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1, bank_transformers::Bool=true)::Dict
    pmd_data = create_data_model()

    _correct_duplicate_components!(dss_data)

    parse_dss_with_dtypes!(dss_data, ["line", "linecode", "load", "generator", "capacitor",
                                      "reactor", "circuit", "transformer", "pvsystem",
                                      "storage"])

    if haskey(dss_data, "options")
        condensed_opts = [Dict{String,Any}()]
        for opt in dss_data["options"]
            merge!(condensed_opts[1], opt)
        end
        dss_data["options"] = condensed_opts
    end

    merge!(pmd_data, parse_options(get(dss_data, "options", [Dict{String,Any}()])[1]))

    #pmd_data["per_unit"] = false
    #pmd_data["source_type"] = "dss"
    #pmd_data["source_version"] = string(VersionNumber("0"))

    if haskey(dss_data, "circuit")
        circuit = dss_data["circuit"][1]
        defaults = _create_vsource(get(circuit, "bus1", "sourcebus"), circuit["name"]; _to_sym_keys(circuit)...)

        #pmd_data["name"] = defaults["name"]
        pmd_data["settings"]["v_var_scalar"] = 1E3
        pmd_data["settings"]["set_vbase_val"] = (defaults["basekv"]/sqrt(3))*1E3/pmd_data["settings"]["v_var_scalar"]
        pmd_data["settings"]["set_vbase_bus"] = "sourcebus"

        pmd_data["settings"]["set_sbase_val"] = defaults["basemva"]*1E6/pmd_data["settings"]["v_var_scalar"]
        pmd_data["settings"]["basefreq"] = pop!(pmd_data, "defaultbasefreq")
        #pmd_data["pu"] = defaults["pu"]
        #pmd_data["conductors"] = defaults["phases"]
        #pmd_data["sourcebus"] = defaults["bus1"]
    else
        Memento.error(_LOGGER, "Circuit not defined, not a valid circuit!")
    end

    _dss2pmd_bus_dm!(pmd_data, dss_data, import_all, vmin, vmax)
    _dss2pmd_line_dm!(pmd_data, dss_data, import_all)
    _dss2pmd_transformer_dm!(pmd_data, dss_data, import_all)


    _dss2pmd_load_dm!(pmd_data, dss_data, import_all)
    _dss2pmd_shunt_dm!(pmd_data, dss_data, import_all)


    #_dss2pmd_reactor!(pmd_data, dss_data, import_all)
    #_dss2pmd_gen!(pmd_data, dss_data, import_all)
    #_dss2pmd_pvsystem!(pmd_data, dss_data, import_all)
    #_dss2pmd_storage!(pmd_data, dss_data, import_all)

    #pmd_data["dcline"] = []
    #pmd_data["switch"] = []

    #InfrastructureModels.arrays_to_dicts!(pmd_data)

    if bank_transformers
        _bank_transformers!(pmd_data)
    end

    #_create_sourcebus_vbranch_dm!(pmd_data, defaults)

    _discover_terminals!(pmd_data)

    #_adjust_sourcegen_bounds!(pmd_data)

    #pmd_data["files"] = dss_data["filename"]

    return pmd_data
end

function _discover_terminals!(pmd_data)
    terminals = Dict{String, Set{Int}}([(bus["id"], Set{Int}()) for (_,bus) in pmd_data["bus"]])

    for (_,line) in pmd_data["line"]
        # ignore 0 terminal
        push!(terminals[line["f_bus"]], setdiff(line["f_connections"], [0])...)
        push!(terminals[line["t_bus"]], setdiff(line["t_connections"], [0])...)
    end

    if haskey(pmd_data, "transformer_nw")
        for (_,tr) in pmd_data["transformer_nw"]
            for w in 1:length(tr["bus"])
                # ignore 0 terminal
                push!(terminals[tr["bus"][w]], setdiff(tr["connections"][w], [0])...)
            end
        end
    end

    for (id, bus) in pmd_data["bus"]
        pmd_data["bus"][id]["terminals"] = sort(collect(terminals[id]))
    end

    # identify neutrals and propagate along cables
    bus_neutral = _find_neutrals(pmd_data)

    for (id,bus) in pmd_data["bus"]
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
                    @show comp["connections"]
                end
                #delete!(bus, "awaiting_ground")
            end
        end
        phases = haskey(bus, "neutral") ? setdiff(bus["terminals"], bus["neutral"]) : bus["terminals"]
        bus["phases"] = phases
    end
end

function _find_neutrals(pmd_data)
    vertices = [(id, t) for (id, bus) in pmd_data["bus"] for t in bus["terminals"]]
    neutrals = []
    edges = Set([((line["f_bus"], line["f_connections"][c]),(line["t_bus"], line["t_connections"][c])) for (id, line) in pmd_data["line"] for c in 1:length(line["f_connections"])])

    bus_neutrals = [(id,bus["neutral"]) for (id,bus) in pmd_data["bus"] if haskey(bus, "neutral")]
    trans_neutrals = []
    for (_, tr) in pmd_data["transformer_nw"]
        for w in 1:length(tr["connections"])
            if tr["configuration"][w] == "wye"
                push!(trans_neutrals, (tr["bus"][w], tr["connections"][w][end]))
            end
        end
    end
    load_neutrals = [(load["bus"],load["connections"][end]) for (_,load) in pmd_data["load"] if load["configuration"]=="wye"]
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
function parse_opendss_dm(io::IOStream; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1, bank_transformers::Bool=true)::Dict
    dss_data = parse_dss(io)

    return parse_opendss_dm(dss_data; import_all=import_all)
end


"""
    _get_conductors_ordered(busname; neutral=true)

Returns an ordered list of defined conductors. If ground=false, will omit any `0`
"""
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
        Memento.error("An incorrect number of nodes was specified; |$(parts[2])|!=$(length(default)).")
    end
    return ret
end
