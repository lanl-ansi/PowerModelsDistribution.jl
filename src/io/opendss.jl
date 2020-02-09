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
        Memento.error(_LOGGER, "dss_data has no branches!")
    else
        return buses
    end
end


"""
    _dss2pmd_bus!(pmd_data, dss_data)

Adds PowerModels-style buses to `pmd_data` from `dss_data`.
"""
function _dss2pmd_bus!(pmd_data::Dict, dss_data::Dict, import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1)
    if !haskey(pmd_data, "bus")
        pmd_data["bus"] = []
    end

    nconductors = pmd_data["conductors"]
    buses = _discover_buses(dss_data)
    for (n, (bus, nodes)) in enumerate(buses)
        busDict = Dict{String,Any}()

        busDict["bus_i"] = n
        busDict["index"] = n
        busDict["name"] = bus

        busDict["bus_type"] = 1

        busDict["vm"] = _parse_array(1.0, nodes, nconductors)
        busDict["va"] = _parse_array([_wrap_to_180(-rad2deg(2*pi/nconductors*(i-1))) for i in 1:nconductors], nodes, nconductors)

        busDict["vmin"] = _parse_array(vmin, nodes, nconductors, vmin)
        busDict["vmax"] = _parse_array(vmax, nodes, nconductors, vmax)

        busDict["base_kv"] = pmd_data["basekv"]

        push!(pmd_data["bus"], busDict)
    end

    # create virtual sourcebus
    circuit = _create_vsource(get(dss_data["circuit"][1], "bus1", "sourcebus"), dss_data["circuit"][1]["name"]; _to_sym_keys(dss_data["circuit"][1])...)

    busDict = Dict{String,Any}()

    nodes = Array{Bool}([1 1 1 0])
    ph1_ang = circuit["angle"]
    vm = circuit["pu"]
    vmi = circuit["pu"] - circuit["pu"] / (circuit["mvasc3"] / circuit["basemva"])
    vma = circuit["pu"] + circuit["pu"] / (circuit["mvasc3"] / circuit["basemva"])

    busDict["bus_i"] = length(pmd_data["bus"])+1
    busDict["index"] = length(pmd_data["bus"])+1
    busDict["name"] = "virtual_sourcebus"

    busDict["bus_type"] = 3

    busDict["vm"] = _parse_array(vm, nodes, nconductors)
    busDict["va"] = _parse_array([_wrap_to_180(-rad2deg(2*pi/nconductors*(i-1))+ph1_ang) for i in 1:nconductors], nodes, nconductors)

    busDict["vmin"] = _parse_array(vmi, nodes, nconductors, vmi)
    busDict["vmax"] = _parse_array(vma, nodes, nconductors, vma)

    busDict["base_kv"] = pmd_data["basekv"]

    push!(pmd_data["bus"], busDict)
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
function _dss2pmd_load!(pmd_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(pmd_data, "load")
        pmd_data["load"] = []
    end

    for load in get(dss_data, "load", [])
        _apply_like!(load, dss_data, "load")
        defaults = _apply_ordered_properties(_create_load(load["bus1"], load["name"]; _to_sym_keys(load)...), load)

        loadDict = Dict{String,Any}()

        nconductors = pmd_data["conductors"]
        name, nodes = _parse_busname(defaults["bus1"])

        kv = defaults["kv"]

        expected_kv = pmd_data["basekv"] / sqrt(pmd_data["conductors"])
        if !isapprox(kv, expected_kv; atol=expected_kv * 0.01)
            Memento.warn(_LOGGER, "Load has kv=$kv, not the expected kv=$(expected_kv). Results may not match OpenDSS")
        end

        loadDict["name"] = defaults["name"]
        loadDict["load_bus"] = find_bus(name, pmd_data)

        load_name = defaults["name"]

        cnds = [1, 2, 3, 0][nodes[1,:]]
        loadDict["conn"] = defaults["conn"]
        nph = defaults["phases"]
        delta_map = Dict([1,2]=>1, [2,3]=>2, [1,3]=>3)
        if nph==1
            # PMD convention is to specify the voltage across the load
            # this what OpenDSS does for 1-phase loads only
            loadDict["vnom_kv"] = kv
            # default is to connect betwheen L1 and N
            cnds = (isempty(cnds)) ? [1, 0] : cnds
            # if only one connection specified, implicitly connected to N
            # bus1=x.c == bus1=x.c.0
            if length(cnds)==1
                cnds = [cnds..., 0]
            end
            # if more than two, only first two are considered
            if length(cnds)>2
                # this no longer works if order is not preserved
                # throw an error instead of behaving like OpenDSS
                # cnds = cnds[1:2]
                Memento.error(_LOGGER, "A 1-phase load cannot specify more than two terminals.")
            end
            # conn property has no effect
            # delta/wye is determined by whether load connected to ground
            # or between two phases
            if cnds==[0, 0]
                pqd_premul = zeros(3)
            elseif cnds[2]==0 || cnds[1]==0
                # this is a wye load in the PMD sense
                loadDict["conn"] = "wye"
                ph = (cnds[2]==0) ? cnds[1] : cnds[2]
                pqd_premul = zeros(3)
                pqd_premul[ph] = 1
            else
                # this is a delta load in the PMD sense
                loadDict["conn"] = "delta"
                pqd_premul = zeros(3)
                pqd_premul[delta_map[cnds]] = 1
            end
        elseif nph==2
            # there are some extremely weird edge cases for this
            # the user can enter weird stuff and OpenDSS will still show some result
            # for example, take
            # nphases=3 bus1=x.1.2.0
            # this looks like a combination of a single-phase PMD delta and wye load
            # so throw an error and ask to reformulate as single and three phase loads
            Memento.error(_LOGGER, "Two-phase loads (nphases=2) are not supported, as these lead to unexpected behaviour. Reformulate this load as a combination of single-phase loads.")
        elseif nph==3
            # for 2 and 3 phase windings, kv is always in LL, also for wye
            # whilst PMD model uses actual voltage across load; so LN for wye
            if loadDict["conn"]=="wye"
                loadDict["vnom_kv"] = kv/sqrt(3)
            else
                loadDict["vnom_kv"] = kv
            end
            if cnds==[]
                pqd_premul = [1/3, 1/3, 1/3]
            else
                if (length(cnds)==3 || length(cnds)==4) && cnds==unique(cnds)
                #variations of [1, 2, 3] and [1, 2, 3, 0]
                    pqd_premul = [1/3, 1/3, 1/3]
                else
                    Memento.error(_LOGGER, "Specified connections for three-phase load $name not allowed.")
                end
            end
        else
            Memento.error(_LOGGER, "For a load, nphases should be in [1,3].")
        end
        loadDict["pd"] = pqd_premul.*defaults["kw"]./1e3
        loadDict["qd"] = pqd_premul.*defaults["kvar"]./1e3

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
        loadDict["model"] = model_int2str[model]

        loadDict["status"] = convert(Int, defaults["enabled"])

        loadDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
        loadDict["source_id"] = "load.$load_name"

        loadDict["index"] = length(pmd_data["load"]) + 1

        used = ["phases", "bus1", "name"]
        _PMs._import_remaining!(loadDict, defaults, import_all; exclude=used)

        push!(pmd_data["load"], loadDict)
    end
end


"""
    _dss2pmd_shunt!(pmd_data, dss_data, import_all)

Adds PowerModels-style shunts to `pmd_data` from `dss_data`.
"""
function _dss2pmd_shunt!(pmd_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(pmd_data, "shunt")
        pmd_data["shunt"] = []
    end

    for shunt in get(dss_data, "capacitor", [])
        _apply_like!(shunt, dss_data, "capacitor")
        defaults = _apply_ordered_properties(_create_capacitor(shunt["bus1"], shunt["name"]; _to_sym_keys(shunt)...), shunt)

        shuntDict = Dict{String,Any}()

        nconductors = pmd_data["conductors"]
        name, nodes = _parse_busname(defaults["bus1"])

        vnom_ln = defaults["kv"]
        # 'kv' is specified as phase-to-phase for phases=2/3 (unsure for 4 and more)
        if defaults["phases"] > 1
            vnom_ln = vnom_ln/sqrt(3)
        end
        # 'kvar' is specified for all phases at once; we want the per-phase one, in MVar
        qnom = (defaults["kvar"]/1E3)/defaults["phases"]
        b_cap = qnom[1]/vnom_ln^2
        #  get the base addmittance, with a LN voltage base
        Sbase = 1 # not yet pmd_data["baseMVA"] because this is done in _PMs.make_per_unit
        Ybase_ln = Sbase/(pmd_data["basekv"]/sqrt(3))^2
        # now convent b_cap to per unit
        b_cap_pu = b_cap/Ybase_ln

        b = fill(b_cap_pu, defaults["phases"])
        N = length(b)
        if defaults["conn"]=="wye"
            B = LinearAlgebra.diagm(0=>b)
        else # shunt["conn"]=="delta"
            # create delta transformation matrix Md
            Md = LinearAlgebra.diagm(0=>ones(N), 1=>-ones(N-1))
            Md[N,1] = -1
            B = Md'*LinearAlgebra.diagm(0=>b)*Md
        end

        active_phases = [n for n in 1:nconductors if nodes[n] > 0]

        if length(active_phases) < 3
            Bf = zeros(3, 3)
            Bf[active_phases, active_phases] = B
            B = Bf
        end

        shuntDict["shunt_bus"] = find_bus(name, pmd_data)
        shuntDict["name"] = defaults["name"]
        shuntDict["gs"] = fill(0.0, 3, 3)
        shuntDict["bs"] = B
        shuntDict["status"] = convert(Int, defaults["enabled"])
        shuntDict["index"] = length(pmd_data["shunt"]) + 1

        shuntDict["active_phases"] = active_phases
        shuntDict["source_id"] = "capacitor.$(defaults["name"])"

        used = ["bus1", "phases", "name"]
        _PMs._import_remaining!(shuntDict, defaults, import_all; exclude=used)

        push!(pmd_data["shunt"], shuntDict)
    end


    for shunt in get(dss_data, "reactor", [])
        if !haskey(shunt, "bus2")
            _apply_like!(shunt, dss_data, "reactor")
            defaults = _apply_ordered_properties(_create_reactor(shunt["bus1"], shunt["name"]; _to_sym_keys(shunt)...), shunt)

            shuntDict = Dict{String,Any}()

            nconductors = pmd_data["conductors"]
            name, nodes = _parse_busname(defaults["bus1"])

            Zbase = (pmd_data["basekv"] / sqrt(3.0))^2 * nconductors / pmd_data["baseMVA"]  # Use single-phase base impedance for each phase
            Gcap = Zbase * sum(defaults["kvar"]) / (nconductors * 1e3 * (pmd_data["basekv"] / sqrt(3.0))^2)

            shuntDict["shunt_bus"] = find_bus(name, pmd_data)
            shuntDict["name"] = defaults["name"]
            shuntDict["gs"] = LinearAlgebra.diagm(0=>_parse_array(0.0, nodes, nconductors))  # TODO:
            shuntDict["bs"] = LinearAlgebra.diagm(0=>_parse_array(Gcap, nodes, nconductors))
            shuntDict["status"] = convert(Int, defaults["enabled"])
            shuntDict["index"] = length(pmd_data["shunt"]) + 1

            shuntDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
            shuntDict["source_id"] = "reactor.$(defaults["name"])"

            used = ["bus1", "phases", "name"]
            _PMs._import_remaining!(shuntDict, defaults, import_all; exclude=used)

            push!(pmd_data["shunt"], shuntDict)
        end
    end
end


"""
    _dss2pmd_gen!(pmd_data, dss_data, import_all)

Adds PowerModels-style generators to `pmd_data` from `dss_data`.
"""
function _dss2pmd_gen!(pmd_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(pmd_data, "gen")
        pmd_data["gen"] = []
    end

    # sourcebus generator (created by circuit)
    circuit = dss_data["circuit"][1]
    defaults = _create_vsource(get(circuit, "bus1", "sourcebus"), "sourcegen"; _to_sym_keys(circuit)...)

    genDict = Dict{String,Any}()

    nconductors = pmd_data["conductors"]
    name, nodes = _parse_busname(defaults["bus1"])

    genDict["gen_bus"] = find_bus("virtual_sourcebus", pmd_data)
    genDict["name"] = defaults["name"]
    genDict["gen_status"] = convert(Int, defaults["enabled"])

    # TODO: populate with VSOURCE properties
    genDict["pg"] = _parse_array( 0.0, nodes, nconductors)
    genDict["qg"] = _parse_array( 0.0, nodes, nconductors)

    genDict["qmin"] = _parse_array(-NaN, nodes, nconductors)
    genDict["qmax"] = _parse_array( NaN, nodes, nconductors)

    genDict["pmin"] = _parse_array(-NaN, nodes, nconductors)
    genDict["pmax"] = _parse_array( NaN, nodes, nconductors)

    genDict["model"] = 2
    genDict["startup"] = 0.0
    genDict["shutdown"] = 0.0
    genDict["ncost"] = 3
    genDict["cost"] = [0.0, 1.0, 0.0]

    genDict["index"] = length(pmd_data["gen"]) + 1

    genDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
    genDict["source_id"] = "vsource.$(defaults["name"])"

    genDict["conn"] = "wye"

    used = ["name", "phases", "bus1"]
    _PMs._import_remaining!(genDict, defaults, import_all; exclude=used)

    push!(pmd_data["gen"], genDict)


    for gen in get(dss_data, "generator", [])
        _apply_like!(gen, dss_data, "generator")
        defaults = _apply_ordered_properties(_create_generator(gen["bus1"], gen["name"]; _to_sym_keys(gen)...), gen)

        genDict = Dict{String,Any}()

        nconductors = pmd_data["conductors"]
        name, nodes = _parse_busname(defaults["bus1"])

        genDict["gen_bus"] = find_bus(name, pmd_data)
        genDict["name"] = defaults["name"]
        genDict["gen_status"] = convert(Int, defaults["enabled"])
        genDict["pg"] = _parse_array(defaults["kw"] / (1e3 * nconductors), nodes, nconductors)
        genDict["qg"] = _parse_array(defaults["kvar"] / (1e3 * nconductors), nodes, nconductors)
        genDict["vg"] = _parse_array(defaults["kv"] / pmd_data["basekv"], nodes, nconductors)

        genDict["qmin"] = _parse_array(defaults["minkvar"] / (1e3 * nconductors), nodes, nconductors)
        genDict["qmax"] = _parse_array(defaults["maxkvar"] / (1e3 * nconductors), nodes, nconductors)

        genDict["apf"] = _parse_array(0.0, nodes, nconductors)

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

        genDict["ramp_q"] = _parse_array(max.(abs.(genDict["qmin"]), abs.(genDict["qmax"])), nodes, nconductors)
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

        genDict["conn"] = "wye"

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

        pvDict["pg"] = _parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors)
        pvDict["qg"] = _parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors)
        pvDict["vg"] = _parse_array(defaults["kv"] / pmd_data["basekv"], nodes, nconductors)

        pvDict["pmin"] = _parse_array(0.0, nodes, nconductors)
        pvDict["pmax"] = _parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors)

        pvDict["qmin"] = -_parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors)
        pvDict["qmax"] =  _parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors)

        pvDict["gen_status"] = convert(Int, defaults["enabled"])

        pvDict["model"] = 2
        pvDict["startup"] = 0.0
        pvDict["shutdown"] = 0.0
        pvDict["ncost"] = 3
        pvDict["cost"] = [0.0, 1.0, 0.0]

        pvDict["index"] = length(pmd_data["gen"]) + 1

        pvDict["active_phases"] = [nodes[n] > 0 ? 1 : 0 for n in 1:nconductors]
        pvDict["source_id"] = "pvsystem.$(defaults["name"])"

        pvDict["conn"] = "wye"

        used = ["name", "phases", "bus1"]
        _PMs._import_remaining!(pvDict, defaults, import_all; exclude=used)

        push!(pmd_data["gen"], pvDict)
    end
end


"""
    _dss2pmd_branch!(pmd_data, dss_data, import_all)

Adds PowerModels-style branches to `pmd_data` from `dss_data`.
"""
function _dss2pmd_branch!(pmd_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(pmd_data, "branch")
        pmd_data["branch"] = []
    end

    nconductors = pmd_data["conductors"]

    for line in get(dss_data, "line", [])
        _apply_like!(line, dss_data, "line")

        if haskey(line, "linecode")
            linecode = deepcopy(_get_linecode(dss_data, get(line, "linecode", "")))
            if haskey(linecode, "like")
                linecode = merge(find_component(dss_data, linecode["like"], "linecode"), linecode)
            end

            linecode["units"] = get(line, "units", "none") == "none" ? "none" : get(linecode, "units", "none")
            linecode["circuit_basefreq"] = pmd_data["basefreq"]

            linecode = _create_linecode(get(linecode, "name", ""); _to_sym_keys(linecode)...)
            delete!(linecode, "name")
        else
            linecode = Dict{String,Any}()
        end

        if haskey(line, "basefreq") && line["basefreq"] != pmd_data["basefreq"]
            Memento.warn(_LOGGER, "basefreq=$(line["basefreq"]) on line $(line["name"]) does not match circuit basefreq=$(pmd_data["basefreq"])")
            line["circuit_basefreq"] = pmd_data["basefreq"]
        end

        defaults = _apply_ordered_properties(_create_line(line["bus1"], line["bus2"], line["name"]; _to_sym_keys(line)...), line; linecode=linecode)

        bf, nodes = _parse_busname(defaults["bus1"])

        phase_order_fr = _get_conductors_ordered(defaults["bus1"]; neutral=false, nconductors=nconductors)
        phase_order_to = _get_conductors_ordered(defaults["bus2"]; neutral=false, nconductors=nconductors)

        @assert phase_order_fr == phase_order_to "Order of connections must be identical on either end of a line"

        bt = _parse_busname(defaults["bus2"])[1]

        branchDict = Dict{String,Any}()

        branchDict["name"] = defaults["name"]

        branchDict["f_bus"] = find_bus(bf, pmd_data)
        branchDict["t_bus"] = find_bus(bt, pmd_data)

        branchDict["length"] = defaults["length"]

        rmatrix = _reorder_matrix(_parse_matrix(defaults["rmatrix"], nodes, nconductors), phase_order_fr)
        xmatrix = _reorder_matrix(_parse_matrix(defaults["xmatrix"], nodes, nconductors), phase_order_fr)
        cmatrix = _reorder_matrix(_parse_matrix(defaults["cmatrix"], nodes, nconductors), phase_order_fr)

        Zbase = (pmd_data["basekv"] / sqrt(3))^2 * nconductors / (pmd_data["baseMVA"])

        Zbase = Zbase/3
        # The factor 3 here is needed to convert from a voltage base
        # in line-to-line (LL) to a voltage base in line-to-neutral (LN).
        # V_LL = √3*V_LN
        # Zbase_new = Zbase_old*(Vbase_new/Vbase_old)^2 = Zbase_old*(1/√3)^2
        # In the parser, LL voltage base is used for per unit conversion.
        # However, in the mathematical model, the voltage magnitude per phase
        # is fixed at 1. So implicitly, we later on state that the voltage base
        # is actually in LN. We compensate here for that.

        branchDict["br_r"] = rmatrix * defaults["length"] / Zbase
        branchDict["br_x"] = xmatrix * defaults["length"] / Zbase

        branchDict["g_fr"] = LinearAlgebra.diagm(0=>_parse_array(0.0, nodes, nconductors))
        branchDict["g_to"] = LinearAlgebra.diagm(0=>_parse_array(0.0, nodes, nconductors))

        branchDict["b_fr"] = Zbase * (2.0 * pi * pmd_data["basefreq"] * cmatrix * defaults["length"] / 1e9) / 2.0
        branchDict["b_to"] = Zbase * (2.0 * pi * pmd_data["basefreq"] * cmatrix * defaults["length"] / 1e9) / 2.0

        branchDict["c_rating_a"] = _parse_array(defaults["normamps"], nodes, nconductors)
        branchDict["c_rating_b"] = _parse_array(defaults["emergamps"], nodes, nconductors)
        branchDict["c_rating_c"] = _parse_array(defaults["emergamps"], nodes, nconductors)

        branchDict["tap"] = _parse_array(1.0, nodes, nconductors, 1.0)
        branchDict["shift"] = _parse_array(0.0, nodes, nconductors)

        branchDict["br_status"] = convert(Int, defaults["enabled"])

        branchDict["angmin"] = _parse_array(-60.0, nodes, nconductors, -60.0)
        branchDict["angmax"] = _parse_array( 60.0, nodes, nconductors,  60.0)

        branchDict["transformer"] = false
        branchDict["switch"] = defaults["switch"]

        branchDict["index"] = length(pmd_data["branch"]) + 1

        nodes = .+([_parse_busname(defaults[n])[2] for n in ["bus1", "bus2"]]...)
        branchDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
        branchDict["source_id"] = "line.$(defaults["name"])"

        used = ["name", "bus1", "bus2", "rmatrix", "xmatrix"]
        _PMs._import_remaining!(branchDict, defaults, import_all; exclude=used)

        push!(pmd_data["branch"], branchDict)
    end
end


"""
    _dss2pmd_transformer!(pmd_data, dss_data, import_all)

Adds ThreePhasePowerModels-style transformers to `pmd_data` from `dss_data`.
"""
function _dss2pmd_transformer!(pmd_data::Dict, dss_data::Dict, import_all::Bool)
   if !haskey(pmd_data, "transformer_comp")
        pmd_data["transformer_comp"] = Array{Any,1}()
    end

    for transformer in get(dss_data, "transformer", [])
        _apply_like!(transformer, dss_data, "transformer")
        defaults = _apply_ordered_properties(_create_transformer(transformer["name"]; _to_sym_keys(transformer)...), transformer)

        nconductors = pmd_data["conductors"]
        nrw = defaults["windings"]
        prop_suffix_w = ["", ["_$w" for w in 2:nrw]...]
        if nrw>3
            # All of the code is compatible with any number of windings,
            # except for the parsing of the loss model (the pair-wise reactance)
            Memento.error(_LOGGER, "For now parsing of xscarray is not supported. At most 3 windings are allowed, not $nrw.")
        end

        transDict = Dict{String,Any}()
        transDict["name"] = defaults["name"]
        transDict["source_id"] = "transformer.$(defaults["name"])"
        transDict["buses"] = Array{Int, 1}(undef, nrw)
        if !isempty(defaults["bank"])
            transDict["bank"] = defaults["bank"]
        end
        for i in 1:nrw
            bnstr = defaults["buses"][i]
            bus, nodes = _parse_busname(bnstr)
            active_phases = [n for n in 1:nconductors if nodes[n] > 0]
            nodes_123 = [true true true]
            if !all(nodes[1:3]) && isempty(defaults["bank"])
                Memento.warn(_LOGGER, "Only three-phase transformers are supported. The bus specification $bnstr is treated as $bus instead.")
            elseif !isempty(defaults["bank"])
                if haskey(transDict, "active_phases")
                    if transDict["active_phases"] != active_phases
                        Memento.error(_LOGGER, "Mismatched phase connections on transformer windings not supported when banking transformers")
                    end
                else
                    transDict["active_phases"] = active_phases
                end
            elseif all(nodes[1:3])
                transDict["active_phases"] = [1, 2, 3]
            else
                transDict["active_phases"] = [1, 2, 3]
                Memento.warn(_LOGGER, "Only three-phase transformers are supported. The bus specification $bnstr is treated as $bus instead.")
            end
            transDict["buses"][i] = find_bus(bus, pmd_data)
        end

        # voltage and power ratings
        #transDict["vnom_kv"] = defaults["kvs"]
        #transDict["snom_kva"] = defaults["kvas"]
        transDict["rate_a"] = [ones(nconductors)*defaults["normhkva"] for i in 1:nrw]
        transDict["rate_b"] = [ones(nconductors)*defaults["normhkva"] for i in 1:nrw]
        transDict["rate_c"] = [ones(nconductors)*defaults["emerghkva"] for i  in 1:nrw]
        # convert to 1 MVA base
        transDict["rate_a"] *= 1E-3
        transDict["rate_b"] *= 1E-3
        transDict["rate_c"] *= 1E-3
        # connection properties
        dyz_map = Dict("wye"=>"wye", "delta"=>"delta", "ll"=>"delta", "ln"=>"wye")
        dyz_primary = dyz_map[defaults["conns"][1]]
        transDict["conns"] = Array{String,1}(undef, nrw)

        transDict["config"] = Dict{Int,Any}()
        transDict["config"][1] = Dict(
            "type"=>dyz_primary,
            "polarity"=>'+',
            "cnd"=>[1, 2, 3],
            "grounded"=>true,
            "vm_nom"=>defaults["kvs"][1]
        )
        #transDict["conns"][1] = string("123+", dyz_primary)
        for w in 2:nrw
            type = dyz_map[defaults["conns"][w]]
            if dyz_primary==type
                cnd = [1,2,3]
                polarity = '+'
            else
                if defaults["leadlag"] in ["ansi", "lag"]
                    #Yd1 => (123+y,123+d)
                    #Dy1 => (123+d,231-y)
                    #pp_w = (type=="delta") ? "123+" : "231-"
                    cnd = (type=="delta") ? [1, 2, 3] : [2, 3, 1]
                    polarity = (type=="delta") ? '+' : '-'
                else # hence defaults["leadlag"] in ["euro", "lead"]
                    #Yd11 => (123+y,312-d)
                    #Dy11 => (123+d,123+y)
                    #pp_w = (type=="delta") ? "312-" : "123+"
                    cnd = (type=="delta") ? [3, 1, 2] : [1, 2, 3]
                    polarity = (type=="delta") ? '-' : '+'
                end
            end
            transDict["config"][w] = Dict(
                "type"=>type,
                "polarity"=>polarity,
                "cnd"=>cnd,
                "vm_nom"=>defaults["kvs"][w]
            )
            if type=="wye"
                transDict["config"][w]["grounded"] = true
            end
        end

        # tap properties
        transDict["tm"] = [ones(Float64,3)*defaults["taps"][i] for i in 1:nrw]
        transDict["tm_min"] = [ones(Float64,3)*defaults["mintap"] for i in 1:nrw]
        transDict["tm_max"] = [ones(Float64,3)*defaults["maxtap"] for i in 1:nrw]
        transDict["tm_step"] = [ones(Int,3)*defaults["numtaps"] for i in 1:nrw]
        transDict["fixed"] = [ones(Bool,3) for i in 1:nrw]

        # loss model (converted to SI units, referred to secondary)
        function zpn_to_abc(z, p, n; atol=1E-13)
            a = exp(im*2*pi/3)
            C = 1/sqrt(3)*[1 1 1; 1 a a^2; 1 a^2 a]
            res = inv(C)*[z 0 0; 0 p 0; 0 0 n]*C
            res = (abs.(res).>atol).*res
            return res
        end
        pos_to_abc(p) = zpn_to_abc(p, p, p)
        zbase = 1^2/(defaults["kvas"][1]/1E3)
        transDict["rs"] = Array{Matrix{Float64}, 1}(undef, nrw)
        transDict["gsh"] = Array{Matrix{Float64}, 1}(undef, nrw)
        transDict["bsh"] = Array{Matrix{Float64}, 1}(undef, nrw)
        for w in 1:nrw
            zs_w_p = defaults["%rs"][w]/100*zbase
            Zs_w = pos_to_abc(zs_w_p)

            if haskey(transformer, "rneut") || haskey(transformer, "xneut")
                #TODO handle neutral impedance
                # neutral impedance is ignored for now; all transformers are
                # grounded (that is, those with a wye and zig-zag winding).
                Memento.warn(_LOGGER, "The neutral impedance, (rneut and xneut properties), is ignored; the neutral (for wye and zig-zag windings) is connected directly to the ground.")
            end

            transDict["rs"][w] = real.(Zs_w)
            # shunt elements are added at second winding
            if w==2
                ysh_w_p = (defaults["%noloadloss"]-im*defaults["%imag"])/100/zbase
                Ysh_w = pos_to_abc(ysh_w_p)
                transDict["gsh"][w] = real.(Ysh_w)
                transDict["bsh"][w] = imag.(Ysh_w)
            else
                transDict["gsh"][w] = zeros(Float64, 3, 3)
                transDict["bsh"][w] = zeros(Float64, 3, 3)
            end
        end
        transDict["xs"] = Dict{String, Matrix{Float64}}()

        Zsc = Dict{Tuple{Int,Int}, Complex}()
        if nrw==2
            xs_map = Dict("xhl"=>(1,2))
        elseif nrw==3
            xs_map = Dict("xhl"=>(1,2), "xht"=>(1,3), "xlt"=>(2,3))
        end
        for (k,v) in xs_map
            Zsc[(v)] = im*defaults[k]/100*zbase
        end
        Zbr = _sc2br_impedance(Zsc)
        for (k,zs_ij_p) in Zbr
            Zs_ij = pos_to_abc(zs_ij_p)
            transDict["xs"]["$(k[1])-$(k[2])"] = imag.(Zs_ij)
        end

        push!(pmd_data["transformer_comp"], transDict)
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
    Y = pinv(Zb)
    Y = [-Y*ones(N-1) Y]
    Y = [-ones(1,N-1)*Y; Y]
    # extract elements
    Zbr = Dict()
    for k in keys(Zsc)
        Zbr[k] = (abs(Y[k...])==0) ? Inf : -1/Y[k...]
    end
    return Zbr
end


"""
    _dss2pmd_reactor!(pmd_data, dss_data, import_all)

Adds PowerModels-style branch components based on DSS reactors to `pmd_data` from `dss_data`
"""
function _dss2pmd_reactor!(pmd_data::Dict, dss_data::Dict, import_all::Bool)
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

                reactDict["br_r"] = _parse_matrix(diagm(0 => fill(0.2, nconductors)), nodes, nconductors)
                reactDict["br_x"] = _parse_matrix(zeros(nconductors, nconductors), nodes, nconductors)

                reactDict["g_fr"] = _parse_array(0.0, nodes, nconductors)
                reactDict["g_to"] = _parse_array(0.0, nodes, nconductors)
                reactDict["b_fr"] = _parse_array(0.0, nodes, nconductors)
                reactDict["b_to"] = _parse_array(0.0, nodes, nconductors)

                for key in ["g_fr", "g_to", "b_fr", "b_to"]
                    reactDict[key] = LinearAlgebra.diagm(0=>reactDict[key])
                end

                reactDict["c_rating_a"] = _parse_array(defaults["normamps"], nodes, nconductors)
                reactDict["c_rating_b"] = _parse_array(defaults["emergamps"], nodes, nconductors)
                reactDict["c_rating_c"] = _parse_array(defaults["emergamps"], nodes, nconductors)

                reactDict["tap"] = _parse_array(1.0, nodes, nconductors, NaN)
                reactDict["shift"] = _parse_array(0.0, nodes, nconductors)

                reactDict["br_status"] = convert(Int, defaults["enabled"])

                reactDict["angmin"] = _parse_array(-60.0, nodes, nconductors, -60.0)
                reactDict["angmax"] = _parse_array( 60.0, nodes, nconductors,  60.0)

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
function _dss2pmd_pvsystem!(pmd_data::Dict, dss_data::Dict, import_all::Bool)
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
        pvsystemDict["p"] = _parse_array(defaults["kw"] / 1e3, nodes, nconductors)
        pvsystemDict["q"] = _parse_array(defaults["kvar"] / 1e3, nodes, nconductors)
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
function _dss2pmd_storage!(pmd_data::Dict, dss_data::Dict, import_all::Bool)
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
        storageDict["thermal_rating"] = _parse_array(defaults["kva"] / 1e3 / nconductors, nodes, nconductors)
        storageDict["qmin"] = _parse_array(-defaults["kvar"] / 1e3 / nconductors, nodes, nconductors)
        storageDict["qmax"] = _parse_array( defaults["kvar"] / 1e3 / nconductors, nodes, nconductors)
        storageDict["r"] = _parse_array(defaults["%r"] / 100.0, nodes, nconductors)
        storageDict["x"] = _parse_array(defaults["%x"] / 100.0, nodes, nconductors)
        storageDict["p_loss"] = defaults["%idlingkw"] * defaults["kwrated"] / 1e3
        storageDict["q_loss"] = defaults["%idlingkvar"] * defaults["kvar"] / 1e3

        storageDict["status"] = convert(Int, defaults["enabled"])

        storageDict["ps"] = _parse_array(0.0, nodes, nconductors)
        storageDict["qs"] = _parse_array(0.0, nodes, nconductors)

        storageDict["index"] = length(pmd_data["storage"]) + 1

        storageDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
        storageDict["source_id"] = "storage.$(defaults["name"])"

        used = ["phases", "bus1", "name"]
        _PMs._import_remaining!(storageDict, defaults, import_all; exclude=used)

        push!(pmd_data["storage"], storageDict)
    end
end


"""
    _adjust_sourcegen_bounds!(pmd_data)

Changes the bounds for the sourcebus generator by checking the emergamps of all
of the branches attached to the sourcebus and taking the sum of non-infinite
values. Defaults to Inf if all emergamps connected to sourcebus are also Inf.
This method was updated to include connected transformers as well. It know
has to occur after the call to InfrastructureModels.arrays_to_dicts, so the code
was adjusted to accomodate that.
"""
function _adjust_sourcegen_bounds!(pmd_data)
    emergamps = Array{Float64,1}([0.0])
    sourcebus_n = find_bus(pmd_data["sourcebus"], pmd_data)
    for (_,line) in pmd_data["branch"]
        if (line["f_bus"] == sourcebus_n || line["t_bus"] == sourcebus_n) && !startswith(line["source_id"], "virtual")
            append!(emergamps, get(line, "c_rating_b", get(line, "rate_b", missing)))
        end
    end

    if haskey(pmd_data, "transformer")
        for (_,trans) in pmd_data["transformer"]
            if trans["f_bus"] == sourcebus_n || trans["t_bus"] == sourcebus_n
                append!(emergamps, trans["rate_b"])
            end
        end
    end

    bound = sum(emergamps)

    pmd_data["gen"]["1"]["pmin"] = fill(-bound, size(pmd_data["gen"]["1"]["pmin"]))
    pmd_data["gen"]["1"]["pmax"] = fill( bound, size(pmd_data["gen"]["1"]["pmin"]))
    pmd_data["gen"]["1"]["qmin"] = fill(-bound, size(pmd_data["gen"]["1"]["pmin"]))
    pmd_data["gen"]["1"]["qmax"] = fill( bound, size(pmd_data["gen"]["1"]["pmin"]))

    # set current rating of vbranch modelling internal impedance
    vbranch = [br for (id, br) in pmd_data["branch"] if br["name"]=="sourcebus_vbranch"][1]
    vbranch["rate_a"] = fill(bound, length(vbranch["rate_a"]))
end


"""

    function _decompose_transformers!(pmd_data)

Replaces complex transformers with a composition of ideal transformers and branches
which model losses. New buses (virtual, no physical meaning) are added.
"""
function _decompose_transformers!(pmd_data; import_all::Bool=false)
    if !haskey(pmd_data, "transformer")
        pmd_data["transformer"] = Dict{String, Any}()
    end
    ncnds = pmd_data["conductors"]
    for (tr_id, trans) in pmd_data["transformer_comp"]
        nrw = length(trans["buses"])
        endnode_id_w = Array{Int, 1}(undef, nrw)
        bus_reduce = []
        branch_reduce = []
        # sum ratings for  all windings to have internal worst-case ratings
        rate_a = sum(trans["rate_a"])
        rate_b = sum(trans["rate_b"])
        rate_c = sum(trans["rate_c"])
        for w in 1:nrw
            # 2-WINDING TRANSFORMER
            trans_dict = Dict{String, Any}()
            trans_dict["name"] = "tr$(tr_id)_w$(w)"
            trans_dict["source_id"] = "$(trans["source_id"])_$(w)"
            trans_dict["active_phases"] = [1, 2, 3]
            _push_dict_ret_key!(pmd_data["transformer"], trans_dict)
            # connection settings
            trans_dict["config_fr"] = trans["config"][w]
            trans_dict["config_to"] = Dict(
                "type"=>"wye",
                "polarity"=>'+',
                "cnd"=>[1, 2, 3],
                "grounded"=>true,
                "vm_nom"=>1.0
            )
            trans_dict["f_bus"] = trans["buses"][w]
            # make virtual bus and mark it for reduction
            vbus_tr = _create_vbus!(pmd_data, basekv=1.0, name="tr$(tr_id)_w$(w)_b1")
            trans_dict["t_bus"] = vbus_tr["index"]
            append!(bus_reduce, vbus_tr["index"])
            # convert to baseMVA, because this is not done per_unit now)
            trans_dict["rate_a"] = trans["rate_a"][w]/pmd_data["baseMVA"]
            trans_dict["rate_b"] = trans["rate_b"][w]/pmd_data["baseMVA"]
            trans_dict["rate_c"] = trans["rate_c"][w]/pmd_data["baseMVA"]
            # tap settings
            trans_dict["tm"] = trans["tm"][w]
            trans_dict["fixed"] = trans["fixed"][w]
            trans_dict["tm_max"] = trans["tm_max"][w]
            trans_dict["tm_min"] = trans["tm_min"][w]
            trans_dict["tm_step"] = trans["tm_step"][w]
            # WINDING SERIES RESISTANCE
            # make virtual bus and mark it for reduction
            vbus_br = _create_vbus!(pmd_data, basekv=1.0, name="tr$(tr_id)_w$(w)_b2")
            append!(bus_reduce, vbus_br["index"])
            # make virtual branch and mark it for reduction
            br = _create_vbranch!(
                pmd_data, vbus_tr["index"], vbus_br["index"],
                vbase=1.0,
                br_r=trans["rs"][w],
                g_fr=trans["gsh"][w],
                b_fr=trans["bsh"][w],
                rate_a=rate_a,
                rate_b=rate_b,
                rate_c=rate_c,
                name="tr$(tr_id)_w$(w)_rs"
            )
            append!(branch_reduce, br["index"])
            # save the trailing node for the reactance model
            endnode_id_w[w] = vbus_br["index"]
        end
        # now add the fully connected graph for reactances
        for w in 1:nrw
            for v in w+1:nrw
                br = _create_vbranch!(
                    pmd_data, endnode_id_w[w], endnode_id_w[v],
                    vbase=1.0,
                    br_x=trans["xs"][string(w,"-",v)],
                    rate_a=rate_a,
                    rate_b=rate_b,
                    rate_c=rate_c,
                    name="tr$(tr_id)_xs_$(w)to$(v)"
                )
                append!(branch_reduce, br["index"])
            end
        end
        _rm_redundant_pd_elements!(pmd_data, buses=string.(bus_reduce), branches=string.(branch_reduce))
    end
    # remove the transformer_comp dict unless import_all is flagged
    if !import_all
        delete!(pmd_data, "transformer_comp")
    end
end


"""
This function adds a new bus to the data model and returns its dictionary.
It is virtual in the sense that it does not correspond to a bus in the network,
but is part of the decomposition of the transformer.
"""
function _create_vbus!(pmd_data; vmin=0, vmax=Inf, basekv=pmd_data["basekv"], name="", source_id="")
    vbus = Dict{String, Any}("bus_type"=>"1", "name"=>name)
    vbus_id = _push_dict_ret_key!(pmd_data["bus"], vbus)
    vbus["bus_i"] = vbus_id
    vbus["source_id"] = source_id
    ncnds = pmd_data["conductors"]
    vbus["vm"] = ones(Float64, ncnds)
    vbus["va"] = zeros(Float64, ncnds)
    vbus["vmin"] = ones(Float64, ncnds)*vmin
    vbus["vmax"] = ones(Float64, ncnds)*vmax
    vbus["base_kv"] = basekv
    return vbus
end


"""
This function adds a new branch to the data model and returns its dictionary.
It is virtual in the sense that it does not correspond to a branch in the
network, but is part of the decomposition of the transformer.
"""
function _create_vbranch!(pmd_data, f_bus::Int, t_bus::Int;
    name="", source_id="", active_phases=[1, 2, 3],
    kwargs...)
    ncnd = pmd_data["conductors"]
    kwargs = Dict{Symbol,Any}(kwargs)
    vbase = haskey(kwargs, :vbase) ? kwargs[:vbase] : pmd_data["basekv"]
    # TODO assumes per_unit will be flagged
    sbase = haskey(kwargs, :sbase) ? kwargs[:sbase] : pmd_data["baseMVA"]
    zbase = vbase^2/sbase
    # convert to LN vbase in instead of LL vbase
    zbase *= (1/3)
    vbranch = Dict{String, Any}("f_bus"=>f_bus, "t_bus"=>t_bus, "name"=>name)
    vbranch["active_phases"] = active_phases
    vbranch["source_id"] = "virtual_branch.$name"
    for k in [:br_r, :br_x, :g_fr, :g_to, :b_fr, :b_to]
        if !haskey(kwargs, k)
            vbranch[string(k)] = zeros(ncnd, ncnd)
        else
            if k in [:br_r, :br_x]
                vbranch[string(k)] = kwargs[k]./zbase
            else
                vbranch[string(k)] = kwargs[k].*zbase
            end
        end
    end
    vbranch["angmin"] = -ones(ncnd)*60
    vbranch["angmax"] = ones(ncnd)*60
    vbranch["rate_a"] = get(kwargs, :rate_a, fill(Inf, length(active_phases)))
    vbranch["shift"] = zeros(ncnd)
    vbranch["tap"] = ones(ncnd)
    vbranch["transformer"] = false
    vbranch["switch"] = false
    vbranch["br_status"] = 1
    for k in [:rate_a, :rate_b, :rate_c, :c_rating_a, :c_rating_b, :c_rating_c]
        if haskey(kwargs, k)
            vbranch[string(k)] = kwargs[k]
        end
    end
    _push_dict_ret_key!(pmd_data["branch"], vbranch)
    return vbranch
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
This function removes zero impedance branches. Only for transformer loss model!
Branches with zero impedances are deleted, and one of the buses it connects.
For now, the implementation should only be used on the loss model of
transformers. When deleting buses, references at shunts, loads... should
be updated accordingly. In the current implementation, that is only done
for shunts. The other elements, such as loads, do not appear in the
transformer loss model.
"""
function _rm_redundant_pd_elements!(pmd_data; buses=keys(pmd_data["bus"]), branches=keys(pmd_data["branch"]))
    # temporary dictionary for pi-model shunt elements
    shunts_g = Dict{Int, Any}()
    shunts_b = Dict{Int, Any}()
    for (br_id, br) in pmd_data["branch"]
        f_bus = br["f_bus"]
        t_bus = br["t_bus"]
        # if branch is flagged
        if br_id in branches
            # flags for convenience
            is_selfloop = f_bus==t_bus # guaranteed  to be reducable because branch is flagged
            is_shorted = all(br["br_r"] .==0) && all(br["br_x"] .==0)
            is_reducable = string(f_bus) in buses || string(t_bus) in buses
            if is_shorted && is_reducable
                # choose bus to remove
                rm_bus = (f_bus in buses) ? f_bus :  t_bus
                kp_bus = (rm_bus==f_bus) ? t_bus : f_bus
            elseif is_selfloop
                kp_bus = t_bus
            else
                # nothing to do, go to next branch
                continue
            end
            # move shunts to the bus that will be left
            if !haskey(shunts_g, kp_bus)
                shunts_g[kp_bus] =  zeros(3, 3)
                shunts_b[kp_bus] =  zeros(3, 3)
            end

            # bus shunts are diagonal, but branch shunts can b e full matrices
            # ensure no data is lost by only keeping the diagonal
            # this should not be the case for the current transformer parsing
            for key in ["g_fr", "g_to", "b_fr", "b_to"]
                @assert(all(br[key]-diagm(0=>diag(br[key])).==0))
            end

            shunts_g[kp_bus] .+= br["g_fr"]
            shunts_g[kp_bus] .+= br["g_to"]
            shunts_b[kp_bus] .+= br["b_fr"]
            shunts_b[kp_bus] .+= br["b_to"]
            # remove branch from pmd_data
            delete!(pmd_data["branch"], string(br_id))
            if is_shorted && is_reducable
                # remove bus from pmd_data
                delete!(pmd_data["bus"], string(rm_bus))
                # replace bus references in branches
                for (br_id, br) in  pmd_data["branch"]
                    if br["f_bus"] == rm_bus
                        br["f_bus"] = kp_bus
                    end
                    if br["t_bus"] == rm_bus
                        br["t_bus"] = kp_bus
                    end
                end
                # replace bus references in transformers
                for (_, tr) in pmd_data["transformer"]
                    if tr["f_bus"] == rm_bus
                        tr["f_bus"] = kp_bus
                    end
                    if tr["t_bus"] == rm_bus
                        tr["t_bus"] = kp_bus
                    end
                end
                # replace bus references in gens, loads, shunts, storage
                for comp_type in ["gen", "load", "shunt", "storage"]
                    for (_, comp) in pmd_data[comp_type]
                        if comp["$(comp_type)_bus"] == rm_bus
                            comp["$(comp_type)_bus"] = kp_bus
                        end
                    end
                end
                # fix new shunt buses
                for shunts in [shunts_g, shunts_b]
                    for (bus, shunt) in shunts
                        if bus == rm_bus
                            shunts[kp_bus] .+= shunt
                            delete!(shunts, bus)
                        end
                    end
                end
                # TODO clean up other references to the removed bus
                # like for example loads, generators, ...
                # skipped  for now, not relevant for transformer loss model
                # + pvsystem
                # ...
            end
        elseif f_bus==t_bus
            # this might occur if not all buses and branches are marked for removal
            # a branch in parallel with a removed branch can turn into a self-loop
            # and if that branch is not marked for removal, we end up here
            Memento.error(_LOGGER, "Specified set of buses and branches leads to a self-loop.")
        end
    end
    # create shunts for lumped pi-model shunts
    for (bus, shunt_g) in shunts_g
        shunt_b = shunts_b[bus]
        if !all(shunt_g .==0) || !all(shunt_b  .==0)
            Memento.warn(_LOGGER, "Pi-model shunt was moved to a bus shunt. Off-diagonals will be discarded in the data model.")
            # The shunts are part of PM, and will be scaled later on by make_per_unit,
            # unlike PMD level components. The shunts here originate from PMD level
            # components which were already scaled. Therefore, we have to undo the
            # scaling here to prevent double scaling later on.
            gs = shunt_g./1*pmd_data["baseMVA"]
            bs = shunt_b./1*pmd_data["baseMVA"]
            _add_shunt!(pmd_data, bus, gs=gs,  bs=bs)
        end
    end
end


"""
Helper function to add a new shunt. The shunt element is  always inserted at the
internal bus of the second winding in OpenDSS. If one of the branches of the
loss model connected to this bus, has zero impedance (for example, if XHL==0
or XLT==0 or R[3]==0), then this bus might be removed by
_rm_redundant_pd_elements!, in which case a new shunt should be inserted at the
remaining bus of the removed branch.
"""
function _add_shunt!(pmd_data, bus; gs=zeros(3,3), bs=zeros(3,3), vbase_kv=1, sbase_mva=1)
    # TODO check whether keys are consistent with the actual data model
    shunt_dict = Dict{String, Any}("status"=>1, "shunt_bus"=>bus)
    zbase  = vbase_kv^2/sbase_mva
    shunt_dict["gs"] = gs*zbase
    shunt_dict["bs"] = bs*zbase
    _push_dict_ret_key!(pmd_data["shunt"], shunt_dict, assume_no_gaps=false)
end


"""
    function _adjust_base!(pmd_data)

Updates the voltage base at each bus, so that the ratios of the voltage bases
across a transformer are consistent with the ratios of voltage ratings of the
windings. Default behaviour is to start at the primary winding of the first
transformer, and to propagate from there. Branches are updated; the impedances
and addmittances are rescaled to be consistent with the new voltage bases.
"""
function _adjust_base!(pmd_data; start_at_first_tr_prim=false)
    # initialize arrays etc. for the recursive part
    edges_br = [(br["index"], br["f_bus"], br["t_bus"]) for (br_id_str, br) in pmd_data["branch"]]
    edges_tr = [(tr["index"], tr["f_bus"], tr["t_bus"]) for (tr_id_str, tr) in pmd_data["transformer"]]
    edges_br_visited = Dict{Int, Bool}([(edge[1], false) for edge in edges_br])
    edges_tr_visited = Dict{Int, Bool}([(edge[1], false) for edge in edges_tr])
    bus_ids = [parse(Int, x) for x in keys(pmd_data["bus"])]
    nodes_visited = Dict{Int, Bool}([(bus_id, false) for  bus_id in bus_ids])
    # retrieve old voltage bases from connected nodes before starting
    br_basekv_old = Dict([(br["index"], pmd_data["bus"][string(br["f_bus"])]["base_kv"]) for (br_id_str, br) in pmd_data["branch"]])
    # start from the primary of the first transformer
    if start_at_first_tr_prim && haskey(pmd_data, "transformer") && haskey(pmd_data["transformer"], "1")
        trans_first = pmd_data["transformer"]["1"]
        source = trans_first["f_bus"]
        base_kv_new = trans_first["config_fr"]["vm_nom"]
    else
        # start at type 3 bus if present
        buses_3 = [bus["index"] for (bus_id_str, bus) in pmd_data["bus"] if bus["bus_type"]==3]
        buses_2 = [bus["index"] for (bus_id_str, bus) in pmd_data["bus"] if bus["bus_type"]==2]
        if length(buses_3)>0
            source = buses_3[1]
        elseif length(buses_2)>0
            source = buses_2[1]
        else
            Memento.warn(_LOGGER, "No bus of type 3 found; selecting random bus instead.")
            source = parse(Int, rand(keys(pmd_data["bus"])))
        end
        base_kv_new = pmd_data["basekv"]
    end
    _adjust_base_rec!(pmd_data, source, base_kv_new, nodes_visited, edges_br, edges_br_visited, edges_tr, edges_tr_visited, br_basekv_old)
    if !all(values(nodes_visited))
        Memento.warn(_LOGGER, "The network contains buses which are not reachable from the start node for the change of voltage base.")
    end
end


"""
This is the recursive code that goes with _adjust_base!; _adjust_base!
initializes arrays and other data that is passed along in the calls to this
recursive function. For very large networks, this might have to be rewritten
to not rely on recursion.
"""
function _adjust_base_rec!(pmd_data, source::Int, base_kv_new::Float64, nodes_visited, edges_br, edges_br_visited, edges_tr, edges_tr_visited, br_basekv_old)
    source_dict = pmd_data["bus"][string(source)]
    base_kv_prev = source_dict["base_kv"]
    if !(base_kv_prev≈base_kv_new)
        # only possible when meshed; ensure consistency
        if nodes_visited[source]
            Memento.error(_LOGGER, "Transformer ratings lead to an inconsistent definition for the voltage base at bus $source.")
        end
        source_dict["base_kv"] = base_kv_new
        # update the connected shunts with the new voltage base
        source_shunts = [shunt for (sh_id_str, shunt) in pmd_data["shunt"] if shunt["shunt_bus"]==source]
        for shunt in source_shunts
            _adjust_base_shunt!(pmd_data, shunt["index"], base_kv_prev, base_kv_new)
        end
        source_name = haskey(source_dict, "name") ? source_dict["name"] : ""
        if source_dict["bus_type"]==3
            #TODO is this the desired behaviour, keep SI units for type 3 bus?
            source_dict["vm"] *= base_kv_prev/base_kv_new
            source_dict["vmax"] *= base_kv_prev/base_kv_new
            source_dict["vmin"] *= base_kv_prev/base_kv_new
            Memento.info(_LOGGER, "Rescaling vm, vmin and vmax conform with new base_kv at type 3 bus $source($source_name): $base_kv_prev => $base_kv_new")
        else
            Memento.info(_LOGGER, "Resetting base_kv at bus $source($source_name): $base_kv_prev => $base_kv_new")
        end
        # TODO rescale vmin, vmax, vm
        # what is the desired behaviour here?
        # should the p.u. set point stay the same, or the set point in SI units?
    end
    nodes_visited[source] = true
    # propagate through the connected branches
    for (br_id, f_bus, t_bus) in [edge for edge in edges_br if !edges_br_visited[edge[1]]]
        # check !edges_br_visited[edge[1]] again, might be visited by now
        if (f_bus==source || t_bus==source) && !edges_br_visited[br_id]
            # this edge will be visited
            edges_br_visited[br_id] = true
            source_new = (f_bus==source) ? t_bus : f_bus
            # assume the branch was undimensionalised with the basekv of the node
            # it is connected to; ideally this will be a property of the branch
            # itself in the future to ensure consistency
            base_kv_branch_prev = br_basekv_old[br_id]
            if base_kv_branch_prev != base_kv_new
                br = pmd_data["branch"]["$br_id"]
                br_name = haskey(br, "name") ? br["name"] : ""
                Memento.info(_LOGGER, "Rescaling impedances at branch $br_id($br_name), conform with change of voltage base: $base_kv_branch_prev => $base_kv_new")
                _adjust_base_branch!(pmd_data, br_id, base_kv_branch_prev, base_kv_new)
            end
            # follow the edge to the adjacent node and repeat
            _adjust_base_rec!(pmd_data, source_new, base_kv_new, nodes_visited, edges_br, edges_br_visited, edges_tr, edges_tr_visited, br_basekv_old)
        end
    end
    # propogate through the connected transformers
    for (tr_id, f_bus, t_bus) in [edge for edge in edges_tr if !edges_tr_visited[edge[1]]]
        if f_bus==source || t_bus==source
            # this edge is now being visited
            edges_tr_visited[tr_id] = true
            source_new = (f_bus==source) ? t_bus : f_bus
            # scale the basekv across the transformer
            trans = pmd_data["transformer"][string(tr_id)]
            base_kv_new_tr = deepcopy(base_kv_new)
            if source_new==t_bus
                base_kv_new_tr *= (trans["config_to"]["vm_nom"]/trans["config_fr"]["vm_nom"])
            else
                base_kv_new_tr *= (trans["config_fr"]["vm_nom"]/trans["config_to"]["vm_nom"])
            end
            # follow the edge to the adjacent node and repeat
            _adjust_base_rec!(pmd_data, source_new, base_kv_new_tr, nodes_visited, edges_br, edges_br_visited, edges_tr, edges_tr_visited, br_basekv_old)
        end
    end
end


"Rescales the parameters of a branch to reflect a change in voltage base"
function _adjust_base_branch!(pmd_data, br_id::Int, base_kv_old::Float64, base_kv_new::Float64)
    branch = pmd_data["branch"][string(br_id)]
    zmult = (base_kv_old/base_kv_new)^2
    branch["br_r"] *= zmult
    branch["br_x"] *= zmult
    branch["g_fr"] *= 1/zmult
    branch["b_fr"] *= 1/zmult
    branch["g_to"] *= 1/zmult
    branch["b_to"] *= 1/zmult
end


"Rescales the parameters of a shunt to reflect a change in voltage base"
function _adjust_base_shunt!(pmd_data, sh_id::Int, base_kv_old::Float64, base_kv_new::Float64)
    shunt = pmd_data["shunt"][string(sh_id)]
    zmult = (base_kv_old/base_kv_new)^2
    shunt["bs"] *= 1/zmult
    shunt["gs"] *= 1/zmult
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
function _create_sourcebus_vbranch!(pmd_data::Dict, circuit::Dict)
    sourcebus = find_bus(pmd_data["sourcebus"], pmd_data)
    vsourcebus = find_bus("virtual_sourcebus", pmd_data)

    br_r = circuit["rmatrix"]
    br_x = circuit["xmatrix"]

    vbranch = _create_vbranch!(pmd_data, sourcebus, vsourcebus; name="sourcebus_vbranch", br_r=br_r, br_x=br_x)
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
                if isa(v, Vector) && eltype(v) <: Vector
                    # most properties are arrays (indexed over the windings)
                    for w in 1:length(v)
                        banked_transformer[k][w][phase] = deepcopy(transformer[k][w][phase])
                    end
                elseif isa(v, Vector) && eltype(v) <: Matrix
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
function parse_opendss(dss_data::Dict; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1, bank_transformers::Bool=true)::Dict
    pmd_data = Dict{String,Any}()

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

    pmd_data["per_unit"] = false
    pmd_data["source_type"] = "dss"
    pmd_data["source_version"] = string(VersionNumber("0"))

    if haskey(dss_data, "circuit")
        circuit = dss_data["circuit"][1]
        defaults = _create_vsource(get(circuit, "bus1", "sourcebus"), circuit["name"]; _to_sym_keys(circuit)...)

        pmd_data["name"] = defaults["name"]
        pmd_data["basekv"] = defaults["basekv"]
        pmd_data["baseMVA"] = defaults["basemva"]
        pmd_data["basefreq"] = pop!(pmd_data, "defaultbasefreq")
        pmd_data["pu"] = defaults["pu"]
        pmd_data["conductors"] = defaults["phases"]
        pmd_data["sourcebus"] = defaults["bus1"]
    else
        Memento.error(_LOGGER, "Circuit not defined, not a valid circuit!")
    end

    _dss2pmd_bus!(pmd_data, dss_data, import_all, vmin, vmax)
    _dss2pmd_load!(pmd_data, dss_data, import_all)
    _dss2pmd_shunt!(pmd_data, dss_data, import_all)
    _dss2pmd_branch!(pmd_data, dss_data, import_all)
    _dss2pmd_transformer!(pmd_data, dss_data, import_all)
    _dss2pmd_reactor!(pmd_data, dss_data, import_all)
    _dss2pmd_gen!(pmd_data, dss_data, import_all)
    _dss2pmd_pvsystem!(pmd_data, dss_data, import_all)
    _dss2pmd_storage!(pmd_data, dss_data, import_all)

    pmd_data["dcline"] = []
    pmd_data["switch"] = []

    InfrastructureModels.arrays_to_dicts!(pmd_data)

    if bank_transformers
        _bank_transformers!(pmd_data)
    end

    for optional in ["dcline", "load", "shunt", "storage", "pvsystem", "branch"]
        if length(pmd_data[optional]) == 0
            pmd_data[optional] = Dict{String,Any}()
        end
    end

    _create_sourcebus_vbranch!(pmd_data, defaults)

    if haskey(pmd_data, "transformer_comp")
        # this has to be done before calling _adjust_sourcegen_bounds!
        _decompose_transformers!(pmd_data; import_all=import_all)
        _adjust_base!(pmd_data)
    else
        pmd_data["transformer"] = Dict{String, Any}()
    end

    _adjust_sourcegen_bounds!(pmd_data)

    pmd_data["files"] = dss_data["filename"]

    return pmd_data
end


"Parses a DSS file into a PowerModels usable format"
function parse_opendss(io::IOStream; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1, bank_transformers::Bool=true)::Dict
    dss_data = parse_dss(io)

    return parse_opendss(dss_data; import_all=import_all)
end
