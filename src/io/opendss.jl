# OpenDSS parser


"Structure representing OpenDSS `dss_source_id` giving the type of the component `dss_type`, its name `dss_name`, and the active phases `active_phases`"
struct DSSSourceId
    dss_type::AbstractString
    dss_name::AbstractString
    active_phases::Set{Int}
end


"Parses a component's OpenDSS source information into the `dss_source_id` struct"
function parse_dss_source_id(component::Dict)::DSSSourceId
    dss_type, dss_name = split(component["source_id"], '.')
    return DSSSourceId(dss_type, dss_name, Set(component["active_phases"]))
end


"returns the linecode with name `id`"
function get_linecode(dss_data::Dict, id::AbstractString)
    if haskey(dss_data, "linecode")
        for item in dss_data["linecode"]
            if item["name"] == id
                return item
            end
        end
    end
    return Dict{String,Any}()
end


"creates a starbus from a 3-winding transformer"
function create_starbus(tppm_data::Dict, transformer::Dict)::Dict
    starbus = Dict{String,Any}()

    base = convert(Int, 10^ceil(log10(abs(PMs.find_max_bus_id(tppm_data)))))
    name, nodes = parse_busname(transformer["buses"][1])
    nconductors = tppm_data["conductors"]
    starbus_id = find_bus(name, tppm_data) + base

    starbus["bus_i"] = starbus_id
    starbus["base_kv"] = 1.0
    starbus["vmin"] = PMs.MultiConductorVector(parse_array(0.9, nodes, nconductors, 0.9))
    starbus["vmax"] = PMs.MultiConductorVector(parse_array(1.1, nodes, nconductors, 1.1))
    starbus["name"] = "$(transformer["name"]) starbus"
    starbus["vm"] = PMs.MultiConductorVector(parse_array(1.0, nodes, nconductors))
    starbus["va"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))
    starbus["bus_type"] = 1
    starbus["index"] = starbus_id

    nodes = .+([parse_busname(transformer["buses"][n])[2] for n in length(transformer["buses"])]...)
    starbus["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
    starbus["source_id"] = "transformer.$(transformer["name"])"

    return starbus
end


"""
    discover_buses(dss_data)

Discovers all of the buses (not separately defined in OpenDSS), from "lines".
"""
function discover_buses(dss_data::Dict)::Array
    bus_names = []
    buses = []
    for compType in ["line", "transformer", "reactor"]
        if haskey(dss_data, compType)
            compList = dss_data[compType]
            for compObj in compList
                if compType == "transformer"
                    compObj = createTransformer(compObj["name"]; to_sym_keys(compObj)...)
                    for bus in compObj["buses"]
                        name, nodes = parse_busname(bus)
                        if !(name in bus_names)
                            push!(bus_names, name)
                            push!(buses, (name, nodes))
                        end
                    end
                elseif haskey(compObj, "bus2")
                    for key in ["bus1", "bus2"]
                        name, nodes = parse_busname(compObj[key])
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
        error(LOGGER, "dss_data has no branches!")
    else
        return buses
    end
end


"""
    dss2tppm_bus!(tppm_data, dss_data)

Adds PowerModels-style buses to `tppm_data` from `dss_data`.
"""
function dss2tppm_bus!(tppm_data::Dict, dss_data::Dict, import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1)
    if !haskey(tppm_data, "bus")
        tppm_data["bus"] = []
    end

    buses = discover_buses(dss_data)
    circuit = createVSource("", dss_data["circuit"][1]["name"]; to_sym_keys(dss_data["circuit"][1])...)

    for (n, (bus, nodes)) in enumerate(buses)
        busDict = Dict{String,Any}()

        nconductors = tppm_data["conductors"]
        ph1_ang = bus == "sourcebus" ? circuit["angle"] : 0.0
        vm = bus == "sourcebus" ? circuit["pu"] : 1.0
        vmi = bus == "sourcebus" ? circuit["pu"] - circuit["pu"] / (circuit["mvasc3"] / circuit["basemva"]) : vmin
        vma = bus == "sourcebus" ? circuit["pu"] + circuit["pu"] / (circuit["mvasc3"] / circuit["basemva"]) : vmax

        busDict["bus_i"] = n
        busDict["index"] = n
        busDict["name"] = bus

        busDict["bus_type"] = bus == "sourcebus" ? 3 : 1

        busDict["vm"] = PMs.MultiConductorVector(parse_array(vm, nodes, nconductors))
        busDict["va"] = PMs.MultiConductorVector(parse_array([rad2deg(2*pi/nconductors*(i-1))+ph1_ang for i in 1:nconductors], nodes, nconductors))

        busDict["vmin"] = PMs.MultiConductorVector(parse_array(vmi, nodes, nconductors, vmi))
        busDict["vmax"] = PMs.MultiConductorVector(parse_array(vma, nodes, nconductors, vma))

        busDict["base_kv"] = tppm_data["basekv"]

        push!(tppm_data["bus"], busDict)
    end
end


"""
    find_component(tppm_data, name, compType)

Returns the component of `compType` with `name` from `data` of type
Dict{String,Array}.
"""
function find_component(data::Dict, name::AbstractString, compType::AbstractString)::Dict
    for comp in values(data[compType])
        if comp["name"] == name
            return comp
        end
    end
    warn(LOGGER, "Could not find $compType \"$name\"")
    return Dict{String,Any}()
end


"""
    find_bus(busname, tppm_data)

Finds the index number of the bus in existing data from the given `busname`.
"""
function find_bus(busname::AbstractString, tppm_data::Dict)
    bus = find_component(tppm_data, busname, "bus")
    if haskey(bus, "bus_i")
        return bus["bus_i"]
    else
        error(LOGGER, "cannot find connected bus with id \"$busname\"")
    end
end


"""
    dss2tppm_load!(tppm_data, dss_data)

Adds PowerModels-style loads to `tppm_data` from `dss_data`.
"""
function dss2tppm_load!(tppm_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(tppm_data, "load")
        tppm_data["load"] = []
    end

    if haskey(dss_data, "load")
        for load in dss_data["load"]
            if haskey(load, "like")
                load = merge(find_component(dss_data, load["like"], "load"), load)
            end

            defaults = createLoad(load["bus1"], load["name"]; to_sym_keys(load)...)

            loadDict = Dict{String,Any}()

            nconductors = tppm_data["conductors"]
            name, nodes = parse_busname(defaults["bus1"])

            kv = defaults["kv"]
            expected_kv = tppm_data["basekv"] / sqrt(tppm_data["conductors"])
            if !isapprox(kv, expected_kv; atol=expected_kv * 0.01)
                warn(LOGGER, "Load has kv=$kv, not the expected kv=$(expected_kv). Results may not match OpenDSS")
            end

            loadDict["name"] = defaults["name"]
            loadDict["load_bus"] = find_bus(name, tppm_data)
            loadDict["pd"] = PMs.MultiConductorVector(parse_array(defaults["kw"] / 1e3, nodes, nconductors))
            loadDict["qd"] = PMs.MultiConductorVector(parse_array(defaults["kvar"] / 1e3, nodes, nconductors))
            loadDict["status"] = convert(Int, defaults["enabled"])

            loadDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
            loadDict["source_id"] = "load.$(defaults["name"])"

            loadDict["index"] = length(tppm_data["load"]) + 1

            used = ["phases", "bus1", "name"]
            PMs.import_remaining!(loadDict, defaults, import_all; exclude=used)

            push!(tppm_data["load"], loadDict)
        end
    end
end


"""
    dss2tppm_shunt!(tppm_data, dss_data)

Adds PowerModels-style shunts to `tppm_data` from `dss_data`.
"""
function dss2tppm_shunt!(tppm_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(tppm_data, "shunt")
        tppm_data["shunt"] = []
    end

    if haskey(dss_data, "capacitor")
        for shunt in dss_data["capacitor"]
            if haskey(shunt, "like")
                shunt = merge(find_component(dss_data, shunt["like"], "capacitor"), shunt)
            end

            defaults = createCapacitor(shunt["bus1"], shunt["name"]; to_sym_keys(shunt)...)

            shuntDict = Dict{String,Any}()

            nconductors = tppm_data["conductors"]
            name, nodes = parse_busname(defaults["bus1"])

            Zbase = (tppm_data["basekv"] / sqrt(3.0))^2 * nconductors / tppm_data["baseMVA"]  # Use single-phase base impedance for each phase
            Gcap = -Zbase * sum(defaults["kvar"]) / (nconductors * 1e3 * (tppm_data["basekv"] / sqrt(3.0))^2)

            shuntDict["shunt_bus"] = find_bus(name, tppm_data)
            shuntDict["name"] = defaults["name"]
            shuntDict["gs"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))  # TODO:
            shuntDict["bs"] = PMs.MultiConductorVector(parse_array(Gcap, nodes, nconductors))
            shuntDict["status"] = convert(Int, defaults["enabled"])
            shuntDict["index"] = length(tppm_data["shunt"]) + 1

            shuntDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
            shuntDict["source_id"] = "capacitor.$(defaults["name"])"

            used = ["bus1", "phases", "name"]
            PMs.import_remaining!(shuntDict, defaults, import_all; exclude=used)

            push!(tppm_data["shunt"], shuntDict)
        end
    end


    if haskey(dss_data, "reactor")
        for shunt in dss_data["reactor"]
            if !haskey(shunt, "bus2")
                if haskey(shunt, "like")
                    shunt = merge(find_component(dss_data, shunt["like"], "reactor"), shunt)
                end

                defaults = createReactor(shunt["bus1"], shunt["name"]; to_sym_keys(shunt)...)

                shuntDict = Dict{String,Any}()

                nconductors = tppm_data["conductors"]
                name, nodes = parse_busname(defaults["bus1"])

                Zbase = (tppm_data["basekv"] / sqrt(3.0))^2 * nconductors / tppm_data["baseMVA"]  # Use single-phase base impedance for each phase
                Gcap = Zbase * sum(defaults["kvar"]) / (nconductors * 1e3 * (tppm_data["basekv"] / sqrt(3.0))^2)

                shuntDict["shunt_bus"] = find_bus(name, tppm_data)
                shuntDict["name"] = defaults["name"]
                shuntDict["gs"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))  # TODO:
                shuntDict["bs"] = PMs.MultiConductorVector(parse_array(Gcap, nodes, nconductors))
                shuntDict["status"] = convert(Int, defaults["enabled"])
                shuntDict["index"] = length(tppm_data["shunt"]) + 1

                shuntDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
                shuntDict["source_id"] = "reactor.$(defaults["name"])"

                used = ["bus1", "phases", "name"]
                PMs.import_remaining!(shuntDict, defaults, import_all; exclude=used)

                push!(tppm_data["shunt"], shuntDict)
            end
        end
    end
end


"""
    dss2tppm_gen!(tppm_data, dss_data)

Adds PowerModels-style generators to `tppm_data` from `dss_data`.
"""
function dss2tppm_gen!(tppm_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(tppm_data, "gen")
        tppm_data["gen"] = []
    end

    # sourcebus generator (created by circuit)
    circuit = dss_data["circuit"][1]
    defaults = createVSource("sourcebus", "sourcebus"; to_sym_keys(circuit)...)

    genDict = Dict{String,Any}()

    nconductors = tppm_data["conductors"]
    name, nodes = parse_busname(defaults["bus1"])

    genDict["gen_bus"] = find_bus(name, tppm_data)
    genDict["name"] = defaults["name"]
    genDict["gen_status"] = convert(Int, defaults["enabled"])

    # TODO: populate with VSOURCE properties
    genDict["pg"] = PMs.MultiConductorVector(parse_array( 0.0, nodes, nconductors))
    genDict["qg"] = PMs.MultiConductorVector(parse_array( 0.0, nodes, nconductors))

    genDict["qmin"] = PMs.MultiConductorVector(parse_array(-NaN, nodes, nconductors))
    genDict["qmax"] = PMs.MultiConductorVector(parse_array( NaN, nodes, nconductors))

    genDict["pmin"] = PMs.MultiConductorVector(parse_array(-NaN, nodes, nconductors))
    genDict["pmax"] = PMs.MultiConductorVector(parse_array( NaN, nodes, nconductors))

    genDict["model"] = 2
    genDict["startup"] = 0.0
    genDict["shutdown"] = 0.0
    genDict["ncost"] = 3
    genDict["cost"] = [0.0, 1.0, 0.0]

    genDict["index"] = length(tppm_data["gen"]) + 1

    genDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
    genDict["source_id"] = "vsource.$(defaults["name"])"

    used = ["name", "phases", "bus1"]
    PMs.import_remaining!(genDict, defaults, import_all; exclude=used)

    push!(tppm_data["gen"], genDict)

    if haskey(dss_data, "generator")
        for gen in dss_data["generator"]
            if haskey(gen, "like")
                gen = merge(find_component(dss_data, gen["like"], "generator"), gen)
            end

            defaults = createGenerator(gen["bus1"], gen["name"]; to_sym_keys(gen)...)

            genDict = Dict{String,Any}()

            nconductors = tppm_data["conductors"]
            name, nodes = parse_busname(defaults["bus1"])

            genDict["gen_bus"] = find_bus(name, tppm_data)
            genDict["name"] = defaults["name"]
            genDict["gen_status"] = convert(Int, defaults["enabled"])
            genDict["pg"] = PMs.MultiConductorVector(parse_array(defaults["kw"] / (1e3 * nconductors), nodes, nconductors))
            genDict["qg"] = PMs.MultiConductorVector(parse_array(defaults["kvar"] / (1e3 * nconductors), nodes, nconductors))
            genDict["vg"] = PMs.MultiConductorVector(parse_array(defaults["kv"] / tppm_data["basekv"], nodes, nconductors))

            genDict["qmin"] = PMs.MultiConductorVector(parse_array(defaults["minkvar"] / (1e3 * nconductors), nodes, nconductors))
            genDict["qmax"] = PMs.MultiConductorVector(parse_array(defaults["maxkvar"] / (1e3 * nconductors), nodes, nconductors))

            genDict["apf"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))

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

            genDict["ramp_q"] = PMs.MultiConductorVector(parse_array(max.(abs.(genDict["qmin"].values), abs.(genDict["qmax"].values)), nodes, nconductors))
            genDict["ramp_10"] = genDict["pmax"]
            genDict["ramp_30"] = genDict["pmax"]

            genDict["control_model"] = defaults["model"]

            # if PV generator mode convert attached bus to PV bus
            if genDict["control_model"] == 3
                tppm_data["bus"][genDict["gen_bus"]]["bus_type"] = 2
            end

            genDict["model"] = 2
            genDict["startup"] = 0.0
            genDict["shutdown"] = 0.0
            genDict["ncost"] = 3
            genDict["cost"] = [0.0, 1.0, 0.0]

            genDict["index"] = length(tppm_data["gen"]) + 1

            genDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
            genDict["source_id"] = "generator.$(defaults["name"])"

            used = ["name", "phases", "bus1"]
            PMs.import_remaining!(genDict, defaults, import_all; exclude=used)

            push!(tppm_data["gen"], genDict)
        end
    end

    if haskey(dss_data, "pvsystem")
        for pv in dss_data["pvsystem"]
            warn(LOGGER, "Converting PVSystem \"$(pv["name"])\" into generator with limits determined by OpenDSS property 'kVA'")

            if haskey(pv, "like")
                pv = merge(find_component(dss_data, pv["like"], "pvsystem"), pv)
            end

            defaults = createPVSystem(pv["bus1"], pv["name"]; to_sym_keys(pv)...)

            pvDict = Dict{String,Any}()

            nconductors = tppm_data["conductors"]
            name, nodes = parse_busname(defaults["bus1"])

            pvDict["name"] = defaults["name"]
            pvDict["gen_bus"] = find_bus(name, tppm_data)

            pvDict["pg"] = PMs.MultiConductorVector(parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors))
            pvDict["qg"] = PMs.MultiConductorVector(parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors))
            pvDict["vg"] = PMs.MultiConductorVector(parse_array(defaults["kv"] / tppm_data["basekv"], nodes, nconductors))

            pvDict["pmin"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))
            pvDict["pmax"] = PMs.MultiConductorVector(parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors))

            pvDict["qmin"] = -PMs.MultiConductorVector(parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors))
            pvDict["qmax"] =  PMs.MultiConductorVector(parse_array(defaults["kva"] / (1e3 * nconductors), nodes, nconductors))

            pvDict["gen_status"] = convert(Int, defaults["enabled"])

            pvDict["model"] = 2
            pvDict["startup"] = 0.0
            pvDict["shutdown"] = 0.0
            pvDict["ncost"] = 3
            pvDict["cost"] = [0.0, 1.0, 0.0]

            pvDict["index"] = length(tppm_data["gen"]) + 1

            pvDict["active_phases"] = [nodes[n] > 0 ? 1 : 0 for n in 1:nconductors]
            pvDict["source_id"] = "pvsystem.$(defaults["name"])"

            used = ["name", "phases", "bus1"]
            PMs.import_remaining!(pvDict, defaults, import_all; exclude=used)

            push!(tppm_data["gen"], pvDict)
        end
    end
end


"""
    dss2tppm_branch!(tppm_data, dss_data)

Adds PowerModels-style branches to `tppm_data` from `dss_data`.
"""
function dss2tppm_branch!(tppm_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(tppm_data, "branch")
        tppm_data["branch"] = []
    end


    if haskey(dss_data, "line")
        for line in dss_data["line"]
            if haskey(line, "like")
                line = merge(find_component(dss_data, line["like"], "line"), line)
            end

            if haskey(line, "linecode")
                linecode = deepcopy(get_linecode(dss_data, get(line, "linecode", "")))
                if haskey(linecode, "like")
                    linecode = merge(find_component(dss_data, linecode["like"], "linecode"), linecode)
                end

                linecode["units"] = get(line, "units", "none") == "none" ? "none" : get(linecode, "units", "none")

                linecode = createLinecode(get(linecode, "name", ""); to_sym_keys(linecode)...)
                delete!(linecode, "name")
            else
                linecode = Dict{String,Any}()
            end

            if haskey(line, "basefreq") && line["basefreq"] != tppm_data["basefreq"]
                warn(LOGGER, "basefreq=$(line["basefreq"]) on line $(line["name"]) does not match circuit basefreq=$(tppm_data["basefreq"])")
                line["freq"] = deepcopy(line["basefreq"])
                line["basefreq"] = deepcopy(tppm_data["basefreq"])
            end

            defaults = createLine(line["bus1"], line["bus2"], line["name"]; to_sym_keys(line)...)
            merge!(defaults, linecode)

            bf, nodes = parse_busname(defaults["bus1"])
            bt = parse_busname(defaults["bus2"])[1]

            branchDict = Dict{String,Any}()

            nconductors = tppm_data["conductors"]

            branchDict["name"] = defaults["name"]

            branchDict["f_bus"] = find_bus(bf, tppm_data)
            branchDict["t_bus"] = find_bus(bt, tppm_data)

            branchDict["length"] = defaults["length"]

            rmatrix = parse_matrix(defaults["rmatrix"], nodes, nconductors) * 3  # why?
            xmatrix = parse_matrix(defaults["xmatrix"], nodes, nconductors) * 3  # why?
            cmatrix = parse_matrix(defaults["cmatrix"], nodes, nconductors)
            Zbase = (tppm_data["basekv"] / sqrt(3))^2 * nconductors / (tppm_data["baseMVA"])

            branchDict["br_r"] = PMs.MultiConductorMatrix(rmatrix * defaults["length"] / Zbase)
            branchDict["br_x"] = PMs.MultiConductorMatrix(xmatrix * defaults["length"] / Zbase)

            # CHECK: Do we need to reformulate to use a matrix instead of a vector for g, b?
            branchDict["g_fr"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))
            branchDict["g_to"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))

            if !isdiag(cmatrix)
                info(LOGGER, "Only diagonal elements of cmatrix are used to obtain branch values `b_fr/to`")
            end
            branchDict["b_fr"] = PMs.MultiConductorVector(diag(Zbase * (2.0 * pi * defaults["basefreq"] * cmatrix * defaults["length"] / 1e9) / 2.0))
            branchDict["b_to"] = PMs.MultiConductorVector(diag(Zbase * (2.0 * pi * defaults["basefreq"] * cmatrix * defaults["length"] / 1e9) / 2.0))

            # TODO: pick a better value for emergamps
            branchDict["rate_a"] = PMs.MultiConductorVector(parse_array(defaults["normamps"], nodes, nconductors))
            branchDict["rate_b"] = PMs.MultiConductorVector(parse_array(defaults["emergamps"], nodes, nconductors))
            branchDict["rate_c"] = PMs.MultiConductorVector(parse_array(defaults["emergamps"], nodes, nconductors))

            branchDict["tap"] = PMs.MultiConductorVector(parse_array(1.0, nodes, nconductors, 1.0))
            branchDict["shift"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))

            branchDict["br_status"] = convert(Int, defaults["enabled"])

            branchDict["angmin"] = PMs.MultiConductorVector(parse_array(-60.0, nodes, nconductors, -60.0))
            branchDict["angmax"] = PMs.MultiConductorVector(parse_array( 60.0, nodes, nconductors,  60.0))

            branchDict["transformer"] = false
            branchDict["switch"] = defaults["switch"]

            branchDict["index"] = length(tppm_data["branch"]) + 1

            nodes = .+([parse_busname(defaults[n])[2] for n in ["bus1", "bus2"]]...)
            branchDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
            branchDict["source_id"] = "line.$(defaults["name"])"

            used = ["name", "bus1", "bus2", "rmatrix", "xmatrix"]
            PMs.import_remaining!(branchDict, defaults, import_all; exclude=used)

            push!(tppm_data["branch"], branchDict)
        end
    end
end


"""
    dss2tppm_transformer!(tppm_data, dss_data, import_all)

Adds PowerModels-style transformers (branches) to `tppm_data` from `dss_data`.
"""
function dss2tppm_transformer!(tppm_data::Dict, dss_data::Dict, import_all::Bool)

    if haskey(dss_data, "transformer")
        if !haskey(tppm_data, "trans")
            tppm_data["trans"] = Array{Any,1}()
        end
        for transformer in dss_data["transformer"]
            if haskey(transformer, "like")
                transformer = merge(find_component(dss_data, transformer["like"], "transformer"), transformer)
            end

            defaults = createTransformer(transformer["name"]; to_sym_keys(transformer)...)

            nconductors = tppm_data["conductors"]
            nrw = defaults["windings"]
            if nrw>3
                # All of the code is compatible with any number of windings,
                # except for the parsing of the loss model (the pair-wise reactance)
                error(LOGGER, "For now parsing of xscarray is not supported. At most 3 windings are allowed, not $nrw.")
            end

            transDict = Dict{String,Any}()
            transDict["name"] = defaults["name"]
            transDict["buses"] = [find_bus(parse_busname(x)[1], tppm_data) for x in defaults["buses"]]

            # voltage and power ratings
            transDict["vnom_kv"] = defaults["kvs"]
            #transDict["snom_kva"] = defaults["kvas"]
            transDict["rate_a"] = defaults["normhkva"]./(tppm_data["baseMVA"]*1E3)
            transDict["rate_b"] = defaults["normhkva"]./(tppm_data["baseMVA"]*1E3)
            transDict["rate_c"] = defaults["emerghkva"]./(tppm_data["baseMVA"]*1E3)

            # connection properties
            dyz_map = Dict("wye"=>"y", "delta"=>"d", "ll"=>"d", "ln"=>"y")
            dyz_primary = dyz_map[defaults["conns"][1]]
            transDict["conns"] = Array{String,1}(undef, nrw)
            transDict["conns"][1] = string("123+", dyz_primary)
            for w in 2:nrw
                dyz_w = dyz_map[defaults["conns"][w]]
                if dyz_primary==dyz_w
                    pp_w = "123+"
                else
                    if defaults["leadlag"] in ["ansi", "lag"]
                        #Yd1 => (123+y,123+d)
                        #Dy1 => (123+d,312-y)
                        pp_w = (dyz_w=="d") ? "123+" : "312-"
                    else # hence defaults["leadlag"] in ["euro", "lead"]
                        #Yd11 => (123+y,231-d)
                        #Dy11 => (123+d,123+y)
                        pp_w = (dyz_w=="d") ? "231-" : "123+"
                    end
                end
                transDict["conns"][w] = string(pp_w, dyz_w)
            end

            # tap properties
            transDict["tapset"] = [PMs.MultiConductorVector(ones(Float64,3))*defaults["taps"][i] for i in 1:nrw]
            transDict["tapmin"] = [PMs.MultiConductorVector(ones(Float64,3))*defaults["mintap"] for i in 1:nrw]
            transDict["tapmax"] = [PMs.MultiConductorVector(ones(Float64,3))*defaults["maxtap"] for i in 1:nrw]
            transDict["tapnum"] = [PMs.MultiConductorVector(ones(Int,3))*defaults["numtaps"] for i in 1:nrw]

            # loss model (converted to SI units, referred to secondary)
            function zpn_to_abc(z, p, n; atol=1E-13)
                a = exp(im*2*pi/3)
                C = 1/sqrt(3)*[1 1 1; 1 a a^2; 1 a^2 a]
                res = inv(C)*[z 0 0; 0 p 0; 0 0 n]*C
                res = (abs.(res).>atol).*res
                return res
            end
            pos_to_abc(p) = zpn_to_abc(p, p, p)
            #TODO fix base conversion
            zbase = 1^2/(defaults["kvas"][1]/1E3)
            transDict["rs"] = Array{MultiConductorMatrix{Float64}, 1}(undef, nrw)
            transDict["gsh"] = Array{MultiConductorMatrix{Float64}, 1}(undef, nrw)
            transDict["bsh"] = Array{MultiConductorMatrix{Float64}, 1}(undef, nrw)
            for w in 1:nrw
                zs_w_p = defaults["%rs"][w]/100*zbase
                Zs_w = pos_to_abc(zs_w_p)
                #TODO handle %loadloss property
                #TODO add warning
                transDict["rs"][w] = MultiConductorMatrix(real.(Zs_w))
                # shunt elements are added at second winding
                if w==2
                    ysh_w_p = defaults["%noloadloss"]/100/zbase-im*defaults["%imag"]/100/zbase
                    Ysh_w = pos_to_abc(ysh_w_p)
                    transDict["gsh"][w] = MultiConductorMatrix(real.(Ysh_w))
                    transDict["bsh"][w] = MultiConductorMatrix(imag.(Ysh_w))
                else
                    transDict["gsh"][w] = MultiConductorMatrix(zeros(Float64, 3, 3))
                    transDict["bsh"][w] = MultiConductorMatrix(zeros(Float64, 3, 3))
                end
            end
            transDict["xs"] = Dict{String, MultiConductorMatrix{Float64}}()
            if nrw==2
                xs_map = Dict("xhl"=>"1-2")
            elseif nrw==3
                xs_map = Dict("xhl"=>"1-2", "xht"=>"1-3", "xlt"=>"2-3")
            end
            for (k,v) in xs_map
                zs_ij_p = im*defaults[k]/100*zbase
                Zs_ij = pos_to_abc(zs_ij_p)
                transDict["xs"][v] = MultiConductorMatrix(imag.(Zs_ij))
            end

            push!(tppm_data["trans"], transDict)
        end
    end
end

function dss2tppm_reactor!(tppm_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(tppm_data, "branch")
        tppm_data["branch"] = []
    end

    if haskey(dss_data, "reactor")
        warn(LOGGER, "reactors as constant impedance elements is not yet supported, treating like line")
        for reactor in dss_data["reactor"]
            if haskey(reactor, "bus2")
                if haskey(reactor, "like")
                    reactor = merge(find_component(dss_data, reactor["like"], "reactor"), reactor)
                end

                defaults = createReactor(reactor["bus1"], reactor["name"], reactor["bus2"]; to_sym_keys(reactor)...)

                reactDict = Dict{String,Any}()

                nconductors = tppm_data["conductors"]

                f_bus, nodes = parse_busname(defaults["bus1"])
                t_bus = parse_busname(defaults["bus2"])[1]

                reactDict["name"] = defaults["name"]
                reactDict["f_bus"] = find_bus(f_bus, tppm_data)
                reactDict["t_bus"] = find_bus(t_bus, tppm_data)

                reactDict["br_r"] = PMs.MultiConductorMatrix(parse_matrix(diagm(0 => fill(0.2, nconductors)), nodes, nconductors))
                reactDict["br_x"] = PMs.MultiConductorMatrix(parse_matrix(zeros(nconductors, nconductors), nodes, nconductors))

                reactDict["g_fr"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))
                reactDict["g_to"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))
                reactDict["b_fr"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))
                reactDict["b_to"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))

                reactDict["rate_a"] = PMs.MultiConductorVector(parse_array(defaults["normamps"], nodes, nconductors))
                reactDict["rate_b"] = PMs.MultiConductorVector(parse_array(defaults["emergamps"], nodes, nconductors))
                reactDict["rate_c"] = PMs.MultiConductorVector(parse_array(defaults["emergamps"], nodes, nconductors))

                reactDict["tap"] = PMs.MultiConductorVector(parse_array(1.0, nodes, nconductors, NaN))
                reactDict["shift"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))

                reactDict["br_status"] = convert(Int, defaults["enabled"])

                reactDict["angmin"] = PMs.MultiConductorVector(parse_array(-60.0, nodes, nconductors, -60.0))
                reactDict["angmax"] = PMs.MultiConductorVector(parse_array( 60.0, nodes, nconductors,  60.0))

                reactDict["transformer"] = true

                reactDict["index"] = length(tppm_data["branch"]) + 1

                nodes = .+([parse_busname(defaults[n])[2] for n in ["bus1", "bus2"]]...)
                reactDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
                reactDict["source_id"] = "reactor.$(defaults["name"])"

                used = []
                PMs.import_remaining!(reactDict, defaults, import_all; exclude=used)

                push!(tppm_data["branch"], reactDict)
            end
        end
    end
end


"""
dss2tppm_pvsystem!(tppm_data, dss_data)

Adds PowerModels-style pvsystems to `tppm_data` from `dss_data`.
"""
function dss2tppm_pvsystem!(tppm_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(tppm_data, "pvsystem")
        tppm_data["pvsystem"] = []
    end

    if haskey(dss_data, "pvsystem")
        for pvsystem in dss_data["pvsystem"]
            defaults = createPVSystem(pvsystem["bus1"], pvsystem["name"]; to_sym_keys(pvsystem)...)

            pvsystemDict = Dict{String,Any}()

            nconductors = tppm_data["conductors"]
            name, nodes = parse_busname(defaults["bus1"])

            pvsystemDict["name"] = defaults["name"]
            pvsystemDict["pv_bus"] = find_bus(name, tppm_data)
            pvsystemDict["p"] = PMs.MultiConductorVector(parse_array(defaults["kw"] / 1e3, nodes, nconductors))
            pvsystemDict["q"] = PMs.MultiConductorVector(parse_array(defaults["kvar"] / 1e3, nodes, nconductors))
            pvsystemDict["status"] = convert(Int, defaults["enabled"])

            pvsystemDict["index"] = length(tppm_data["pvsystem"]) + 1

            pvsystemDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
            pvsystemDict["source_id"] = "pvsystem.$(defaults["name"])"

            used = ["phases", "bus1", "name"]
            PMs.import_remaining!(pvsystemDict, defaults, import_all; exclude=used)

            push!(tppm_data["pvsystem"], pvsystemDict)
        end
    end
end


function dss2tppm_storage!(tppm_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(tppm_data, "storage")
        tppm_data["storage"] = []
    end

    if haskey(dss_data, "storage")
        for storage in dss_data["storage"]
            defaults = createStorage(storage["bus1"], storage["name"]; to_sym_keys(storage)...)

            storageDict = Dict{String,Any}()

            nconductors = tppm_data["conductors"]
            name, nodes = parse_busname(defaults["bus1"])

            storageDict["name"] = defaults["name"]
            storageDict["storage_bus"] = find_bus(name, tppm_data)
            storageDict["energy"] = defaults["kwhstored"] / 1e3
            storageDict["energy_rating"] = defaults["kwhrated"] / 1e3
            storageDict["charge_rating"] = defaults["%charge"] * defaults["kwrated"] / 1e3 / 100.0
            storageDict["discharge_rating"] = defaults["%discharge"] * defaults["kwrated"] / 1e3 / 100.0
            storageDict["charge_efficiency"] = defaults["%effcharge"] / 100.0
            storageDict["discharge_efficiency"] = defaults["%effdischarge"] / 100.0
            storageDict["thermal_rating"] = PMs.MultiConductorVector(parse_array(defaults["kva"] / 1e3, nodes, nconductors))
            storageDict["qmin"] = PMs.MultiConductorVector(parse_array(-defaults["kvar"] / 1e3, nodes, nconductors))
            storageDict["qmax"] = PMs.MultiConductorVector(parse_array( defaults["kvar"] / 1e3, nodes, nconductors))
            storageDict["r"] = PMs.MultiConductorVector(parse_array(defaults["%r"] / 100.0, nodes, nconductors))
            storageDict["x"] = PMs.MultiConductorVector(parse_array(defaults["%x"] / 100.0, nodes, nconductors))
            storageDict["standby_loss"] = defaults["%idlingkw"] * defaults["kwrated"] / 1e3
            storageDict["status"] = convert(Int, defaults["enabled"])

            storageDict["index"] = length(tppm_data["storage"]) + 1

            storageDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
            storageDict["source_id"] = "storage.$(defaults["name"])"

            used = ["phases", "bus1", "name"]
            PMs.import_remaining!(storageDict, defaults, import_all; exclude=used)

            push!(tppm_data["storage"], storageDict)
        end
    end
end


"""
    adjust_sourcegen_bounds!(tppm_data)

Changes the bounds for the sourcebus generator by checking the emergamps of all
of the branches attached to the sourcebus and taking the sum of non-infinite
values. Defaults to Inf if all emergamps connected to sourcebus are also Inf.
"""
function adjust_sourcegen_bounds!(tppm_data)
    emergamps = Array{Float64,1}()
    sourcebus_n = find_bus("sourcebus", tppm_data)
    for line in tppm_data["branch"]
        if line["f_bus"] == sourcebus_n || line["t_bus"] == sourcebus_n
            append!(emergamps, line["rate_b"].values)
        end
    end

    bound = sum(emergamps)

    tppm_data["gen"][1]["pmin"] = PMs.MultiConductorVector(fill(-bound, size(tppm_data["gen"][1]["pmin"])))
    tppm_data["gen"][1]["pmax"] = PMs.MultiConductorVector(fill( bound, size(tppm_data["gen"][1]["pmin"])))
    tppm_data["gen"][1]["qmin"] = PMs.MultiConductorVector(fill(-bound, size(tppm_data["gen"][1]["pmin"])))
    tppm_data["gen"][1]["qmax"] = PMs.MultiConductorVector(fill( bound, size(tppm_data["gen"][1]["pmin"])))
end

"""

    function decompose_transformers!(tppm_data)

Replaces complex transformers with a composition of ideal transformers and branches
which model losses. New buses (virtual, no physical meaning) are added.
"""
function decompose_transformers!(tppm_data)
    trans_dict = Dict{String,Any}()
    trans_dec_dict = Dict{String, Any}()
    ncnds = tppm_data["conductors"]
    for (tr_id, trans) in tppm_data["trans"]
        if !haskey(trans, "type") || trans["type"]=="complex"
            #TODO big M for flows
            #s_bigM = calc_bigM_flows(trans)
            s_bigM = [10E3 for i in 1:ncnds]
            nrw = length(trans["buses"])
            endnode_id_w = Array{Int, 1}(undef, nrw)
            trans["dec_map"] = Dict{String, Any}()
            #TODO populate these with references
            #TODO issue: ID not defined yet
            trans["dec_map"]["trans_conn"] = Array{Any,1}(undef,nrw)
            trans["dec_map"]["trans_tap"] = Array{Any,1}(undef,nrw)
            for w in 1:nrw
                last_bus_id = trans["buses"][w]
                last_bus = tppm_data["bus"][string(last_bus_id)]
                # if the voltage base is scaled correctly, these voltage bounds
                # will not further constrain the problem
                vmin = last_bus["vmin"]
                vmax = last_bus["vmax"]
                vmin = MultiConductorVector(ones(3))*0.8
                vmax = MultiConductorVector(ones(3))*1.2
                last_vnom = trans["vnom_kv"][w]
                conn = trans["conns"][w]
                # if not 123+y, we need a connection transformer
                if conn != "123+y"
                    trans_conn = Dict{String, Any}("type"=>"conn", "conn"=>conn)
                    t_bus_id = create_vbus!(tppm_data, vmin, vmax)["index"]
                    trans_conn["f_bus"] = last_bus_id
                    trans_conn["t_bus"] = t_bus_id
                    vmult = 1.0
                    if conn[5]=='d'
                        vmult = sqrt(3)
                    end
                    #TODO handle zig-zag case
                    trans_conn["vnom_kv"] = [last_vnom, last_vnom*vmult]
                    # push transformer and save reference
                    trans_conn_id = push_dict_ret_key!(trans_dict, trans_conn)
                    trans["dec_map"]["trans_conn"][w] = trans_conn_id
                    # shift the trailing bus and rated voltage
                    last_bus_id = t_bus_id
                    last_vnom = trans_conn["vnom_kv"][2]
                end
                # now and the n_w:1 transformer
                trans_tap = Dict{String, Any}("type"=>"tap")
                t_bus_id = create_vbus!(tppm_data, vmin, vmax, basekv=1.0)["index"]
                trans_tap["f_bus"] = last_bus_id
                trans_tap["t_bus"] = t_bus_id
                trans_tap["vnom_kv"] = [last_vnom, 1.0]
                trans_tap["tapset"] = trans["tapset"][w]
                trans_tap["tapfix"] = ones(Bool, 3)
                trans_tap["tapmax"] = trans["tapmax"][w]
                trans_tap["tapmin"] = trans["tapmin"][w]
                trans_tap["tapnum"] = trans["tapnum"][w]
                trans_tap_id = push_dict_ret_key!(trans_dict, trans_tap)
                trans["dec_map"]["trans_tap"][w] = trans_tap_id
                last_bus_id = t_bus_id
                # now we add the series resistance and shunt in p.u.
                t_bus_id = create_vbus!(tppm_data, vmin, vmax, basekv=1.0)["index"]
                create_vbranch!(
                    tppm_data, last_bus_id, t_bus_id,
                    vbase=1.0,
                    br_r=trans["rs"][w],
                    g_to=trans["gsh"][w],
                    b_to=trans["bsh"][w],
                    rate_a=s_bigM,
                    rate_b=s_bigM,
                    rate_c=s_bigM,
                )
                # save the last node for the reactance model
                endnode_id_w[w] = t_bus_id
            end
            # now add the fully connected graph for reactances
            # these branches are virtual; their flow limits should never be binding
            for w in 1:nrw
                for v in w+1:nrw
                    create_vbranch!(
                        tppm_data, endnode_id_w[w], endnode_id_w[v],
                        vbase=1.0,
                        br_x=trans["xs"][string(w,"-",v)],
                        rate_a=s_bigM,
                        rate_b=s_bigM,
                        rate_c=s_bigM,
                    )
                end
            end
            # save the decomposed transformer for solution building
            push_dict_ret_key!(trans_dec_dict, trans)
        else
            warn(LOGGER, string("Transformer \"", trans["name"], "\" was not decomposed!"))
            push_dict_ret_key!(trans_dict, trans)
        end
    end
    # overwrite old transformer dictionary
    tppm_data["trans"] = trans_dict
    tppm_data["trans_dec"] = trans_dec_dict
end
function create_vbus!(tppm_data, vmin, vmax; basekv=tppm_data["basekv"])
    vbus = Dict{String, Any}("bus_type"=>"1")
    vbus_id = push_dict_ret_key!(tppm_data["bus"], vbus)
    vbus["bus_i"] = vbus_id
    vbus["vmin"] = vmin
    vbus["vmax"] = vmax
    ncnds = tppm_data["conductors"]
    vbus["vm"] = MultiConductorVector(ones(Float64, ncnds))
    vbus["va"] = MultiConductorVector(zeros(Float64, ncnds))
    vbus["base_kv"] = basekv
    return vbus
end
function create_vbranch!(tppm_data, f_bus::Int, t_bus::Int; kwargs...)
    ncnd = tppm_data["conductors"]
    vbase = haskey(kwargs, :vbase) ? kwargs[:vbase] : tppm_data["basekv"]
    # TODO assumes per_unit will be flagged
    sbase = haskey(kwargs, :sbase) ? kwargs[:sbase] : tppm_data["baseMVA"]
    zbase = vbase^2/sbase
    # convert to LN vbase in instead of LL vbase
    zbase *= (1/3)
    vbranch = Dict{String, Any}("f_bus"=>f_bus, "t_bus"=>t_bus)
    for k in [:br_r, :br_x, :g_fr, :g_to, :b_fr, :b_to]
        if !haskey(kwargs, k)
            vbranch[string(k)] = MultiConductorMatrix(zeros(ncnd, ncnd))
        else
            if k in [:br_r, :br_x]
                vbranch[string(k)] = kwargs[k]./zbase
            else
                vbranch[string(k)] = kwargs[k].*zbase
            end
        end
    end
    vbranch["angmin"] = -MultiConductorVector(ones(ncnd))*60
    vbranch["angmax"] = MultiConductorVector(ones(ncnd))*60
    vbranch["shift"] = MultiConductorVector(zeros(ncnd))
    vbranch["tap"] = MultiConductorVector(ones(ncnd))
    vbranch["br_status"] = 1
    for k in [:rate_a, :rate_b, :rate_c]
        if haskey(kwargs, k)
            vbranch[string(k)] = kwargs[k]
        end
    end
    push_dict_ret_key!(tppm_data["branch"], vbranch)
    return vbranch
end
function push_dict_ret_key!(dict::Dict{String, Any}, v::Dict{String, Any})
    k = length(keys(dict))+1
    dict[string(k)] = v
    v["index"] = k
    return k
end
function calc_bigM_flows(tppm_data, trans, Shmax)
    nrw = length(trans["buses"])
    # find voltage bounds on all windings;
    # guaranteed to be the highest voltage in loss model
    Vminpu = Inf
    Vmaxpu = 0
    for bus_id in trans["buses"]
        Vmaxpu = max(Vmaxpu, maximum(tppm_data["bus"][string(bus_id)]["vmax"]))
        Vminpu = min(Vminpu, minimum(tppm_data["bus"][string(bus_id)]["vmin"]))
    end
    # loss model is referred to 1 kV
    Vbase = 1.0
    Vmax = Vmaxpu*Vbase
    Vmin = Vminpu*Vbase
    Ishmax = sum(abs.(trans["gsh"][w][:,:]+im*trans["bsh"][w][:,:]) for w in 1:nrw)*ones(3)*Vmax
    println("Ishmax:$Ishmax")
    Ihmax = Shmax./Vmin
    println("Ihmax:$Ihmax")
    Ismax = Ishmax + Ihmax
    Slossmax = (
        sum([abs.(r[:,:]) for r in trans["rs"]])+sum([abs.(x[:,:]) for (_,x) in trans["xs"]])
        )*Ismax.^2
    Smax = Shmax + Slossmax
    return Smax
end
"""

    function adjust_base!(tppm_data)

Updates the voltage base at each bus, so that the ratios of the voltage bases
across a transformer are consistent with the ratios of voltage ratings of the
windings. Default behaviour is to start at the primary winding of the first
transformer, and to propagate from there. Branches are updated; the impedances
and addmittances are rescaled to be consistent with the new voltage bases.
"""
function adjust_base!(tppm_data)
    # initialize arrays etc. for the recursive part
    edges_br = [(br["index"], br["f_bus"], br["t_bus"]) for (br_id_str, br) in tppm_data["branch"]]
    edges_tr = [(tr["index"], tr["f_bus"], tr["t_bus"]) for (tr_id_str, tr) in tppm_data["trans"]]
    edges_br_visited = zeros(Bool, length(edges_br))
    edges_tr_visited = zeros(Bool, length(edges_tr))
    nodes_visited = zeros(Bool, length(keys(tppm_data["bus"])))
    # start from the primary of the first transformer
    if haskey(tppm_data, "trans") && haskey(tppm_data["trans"], "1")
        trans_first = tppm_data["trans"]["1"]
        source = trans_first["f_bus"]
        base_kv_new = trans_first["vnom_kv"][1]
    else
        #TODO find bus of type 3?
        # Only relevant for future per-unit upgrade
        # Impossible to end up here;
        # condition checked before call to adjust_base!
    end
    adjust_base_rec!(tppm_data, source, base_kv_new, nodes_visited, edges_br, edges_br_visited, edges_tr, edges_tr_visited)
    if !all(nodes_visited)
        println(nodes_visited)
        warn(LOGGER, "The network contains buses which are not reachable from the start node for the change of voltage base.")
    end
end
function adjust_base_rec!(tppm_data, source::Int, base_kv_new::Float64, nodes_visited, edges_br, edges_br_visited, edges_tr, edges_tr_visited)
    source_dict = tppm_data["bus"][string(source)]
    base_kv_prev = source_dict["base_kv"]
    if base_kv_prev!=base_kv_new
        # only possible when meshed; ensure consistency
        if nodes_visited[source]
            error(LOGGER, "Transformer ratings lead to an inconsistent definition for the voltage base at bus $source.")
        end
        source_dict["base_kv"] = base_kv_new
        if source_dict["bus_type"]==3
            #TODO is this the desired behaviour, keep SI units for type 3 bus?
            source_dict["vm"] *= base_kv_prev/base_kv_new
            source_dict["vmax"] *= base_kv_prev/base_kv_new
            source_dict["vmin"] *= base_kv_prev/base_kv_new
            info(LOGGER, "Rescaling vm, vmin and vmax conform with new base_kv at type 3 bus $source: $base_kv_prev => $base_kv_new")
        else
            info(LOGGER, "Resetting base_kv at bus $source: $base_kv_prev => $base_kv_new")
        end
        # TODO rescale vmin, vmax, vm
        # what is the desired behaviour here?
        # should the p.u. set point stay the same, or the set point in SI units?
    end
    nodes_visited[source] = true
    # propagate through the connected branches
    for (br_id, f_bus, t_bus) in [edge for edge in edges_br if !edges_br_visited[edge[1]]]
        if f_bus==source || t_bus==source
            # this edge will be visited
            edges_br_visited[br_id] = true
            source_new = (f_bus==source) ? t_bus : f_bus
            # assume the branch was undimensiolised with the basekv of the node
            # it is connected to; ideally this will be a property of the branch
            # itself in the future to ensure consistency
            base_kv_branch_prev = tppm_data["bus"][string(source_new)]["base_kv"]
            if base_kv_branch_prev != base_kv_new
                info(LOGGER, "Rescaling impedances at branch $br_id, conform with change of voltage base: $base_kv_branch_prev => $base_kv_new")
                adjust_base_branch!(tppm_data, br_id, base_kv_prev, base_kv_new)
            end
            # follow the edge to the adjacent node and repeat
            adjust_base_rec!(tppm_data, source_new, base_kv_new, nodes_visited, edges_br, edges_br_visited, edges_tr, edges_tr_visited)
        end
    end
    # propogate through the connected transformers
    for (tr_id, f_bus, t_bus) in [edge for edge in edges_tr if !edges_tr_visited[edge[1]]]
        if f_bus==source || t_bus==source
            # this edge is now being visited
            edges_tr_visited[tr_id] = true
            source_new = (f_bus==source) ? t_bus : f_bus
            # scale the basekv across the transformer
            trans = tppm_data["trans"][string(tr_id)]
            base_kv_new_tr = deepcopy(base_kv_new)
            if source_new==t_bus
                base_kv_new_tr *= (trans["vnom_kv"][2]/trans["vnom_kv"][1])
            else
                base_kv_new_tr *= (trans["vnom_kv"][1]/trans["vnom_kv"][2])
            end
            # follow the edge to the adjacent node and repeat
            adjust_base_rec!(tppm_data, source_new, base_kv_new_tr, nodes_visited, edges_br, edges_br_visited, edges_tr, edges_tr_visited)
        end
    end
end
function adjust_base_branch!(tppm_data, br_id::Int, base_kv_old::Float64, base_kv_new::Float64)
    branch = tppm_data["branch"][string(br_id)]
    zmult = (base_kv_old/base_kv_new)^2
    branch["br_r"] *= zmult
    branch["br_x"] *= zmult
    branch["g_fr"] *= zmult
    branch["b_fr"] *= zmult
    branch["g_to"] *= zmult
    branch["b_to"] *= zmult
end

"""
    where_is_comp(data, comp_id)

Finds existing component of id `comp_id` in array of `data` and returns index.
Assumes all components in `data` are unique.
"""
function where_is_comp(data::Array, comp_id::AbstractString)::Int
    for (i, e) in enumerate(data)
        if e["name"] == comp_id
            return i
        end
    end
    return 0
end


"""
    check_duplicate_components!(dss_data)

Finds duplicate components in `dss_data` and merges up, meaning that older
data (lower indices) is always overwritten by newer data (higher indices).
"""
function check_duplicate_components!(dss_data::Dict)
    out = Dict{String,Array}()
    for (k, v) in dss_data
        if !(k in ["options"])
            out[k] = []
            for comp in v
                if isa(comp, Dict)
                    idx = where_is_comp(out[k], comp["name"])
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


"""
    parse_options(options)

Parses options defined with the `set` command in OpenDSS.
"""
function parse_options(options)
    out = Dict{String,Any}()
    if haskey(options, "voltagebases")
        out["voltagebases"] = parse_array(Float64, options["voltagebases"])
    end

    if !haskey(options, "defaultbasefreq")
        warn(LOGGER, "defaultbasefreq is not defined, default for circuit set to 60 Hz")
        out["defaultbasefreq"] = 60.0
    else
        out["defaultbasefreq"] = parse(Float64, options["defaultbasefreq"])
    end

    return out
end


"Parses a Dict resulting from the parsing of a DSS file into a PowerModels usable format."
function parse_opendss(dss_data::Dict; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1)::Dict
    tppm_data = Dict{String,Any}()

    check_duplicate_components!(dss_data)

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

    merge!(tppm_data, parse_options(get(dss_data, "options", [Dict{String,Any}()])[1]))

    tppm_data["per_unit"] = false
    tppm_data["source_type"] = "dss"
    tppm_data["source_version"] = VersionNumber("0")

    if haskey(dss_data, "circuit")
        circuit = dss_data["circuit"][1]
        defaults = createVSource("", circuit["name"]; to_sym_keys(circuit)...)

        tppm_data["name"] = defaults["name"]
        tppm_data["basekv"] = defaults["basekv"]
        tppm_data["baseMVA"] = defaults["basemva"]
        tppm_data["basefreq"] = pop!(tppm_data, "defaultbasefreq")
        tppm_data["pu"] = defaults["pu"]
        tppm_data["conductors"] = defaults["phases"]
    else
        error(LOGGER, "Circuit not defined, not a valid circuit!")
    end

    dss2tppm_bus!(tppm_data, dss_data, import_all, vmin, vmax)
    dss2tppm_load!(tppm_data, dss_data, import_all)
    dss2tppm_shunt!(tppm_data, dss_data, import_all)
    dss2tppm_branch!(tppm_data, dss_data, import_all)
    dss2tppm_transformer!(tppm_data, dss_data, import_all)
    dss2tppm_reactor!(tppm_data, dss_data, import_all)
    dss2tppm_gen!(tppm_data, dss_data, import_all)
    dss2tppm_pvsystem!(tppm_data, dss_data, import_all)
    dss2tppm_storage!(tppm_data, dss_data, import_all)

    adjust_sourcegen_bounds!(tppm_data)

    tppm_data["dcline"] = []

    InfrastructureModels.arrays_to_dicts!(tppm_data)

    for optional in ["dcline", "load", "shunt", "storage", "pvsystem"]
        if length(tppm_data[optional]) == 0
            tppm_data[optional] = Dict{String,Any}()
        end
    end

    if haskey(tppm_data, "trans")
        decompose_transformers!(tppm_data)
        adjust_base!(tppm_data)
    end

    tppm_data["files"] = dss_data["filename"]

    return tppm_data
end


"Parses a DSS file into a PowerModels usable format."
function parse_opendss(io::IOStream; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1)::Dict
    dss_data = parse_dss(io)

    return parse_opendss(dss_data; import_all=import_all)
end
