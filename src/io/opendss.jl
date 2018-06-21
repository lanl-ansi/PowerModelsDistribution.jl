# OpenDSS parser


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
    phases = tppm_data["phases"]
    starbus_id = find_bus(name, tppm_data) + base

    starbus["bus_i"] = starbus_id
    starbus["base_kv"] = 1.0
    starbus["vmin"] = PMs.MultiPhaseVector(parse_array(0.9, nodes, phases))
    starbus["vmax"] = PMs.MultiPhaseVector(parse_array(1.1, nodes, phases))
    starbus["name"] = "$(transformer["name"]) starbus"
    starbus["vm"] = PMs.MultiPhaseVector(parse_array(1.0, nodes, phases))
    starbus["va"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, phases))
    starbus["bus_type"] = 1
    starbus["index"] = starbus_id

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

        nphases = tppm_data["phases"]
        ph1_ang = bus == "sourcebus" ? circuit["angle"] : 0.0
        vm = bus == "sourcebus" ? circuit["pu"] : 1.0
        vmin = bus == "sourcebus" ? circuit["pu"] - circuit["pu"] / (circuit["mvasc3"] / circuit["basemva"]) : vmin
        vmax = bus == "sourcebus" ? circuit["pu"] + circuit["pu"] / (circuit["mvasc3"] / circuit["basemva"]) : vmax

        busDict["bus_i"] = n
        busDict["index"] = n
        busDict["name"] = bus

        busDict["bus_type"] = bus == "sourcebus" ? 3 : 1

        busDict["vm"] = PMs.MultiPhaseVector(parse_array(vm, nodes, nphases))
        busDict["va"] = PMs.MultiPhaseVector(parse_array([rad2deg(2*pi/nphases*(i-1))+ph1_ang for i in 1:nphases], nodes, nphases))

        busDict["vmin"] = PMs.MultiPhaseVector(parse_array(vmin, nodes, nphases))
        busDict["vmax"] = PMs.MultiPhaseVector(parse_array(vmax, nodes, nphases))

        busDict["base_kv"] = tppm_data["basekv"]

        push!(tppm_data["bus"], busDict)
    end
end


"""
    find_bus(busname, tppm_data)

Finds the index number of the bus in existing data from the given `busname`.
"""
function find_bus(busname, tppm_data)
    for bus in values(tppm_data["bus"])
        if bus["name"] == busname
            return bus["bus_i"]
        end
    end
    error("cannot find connected bus with id $busname")
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
            defaults = createLoad(load["bus1"], load["name"]; to_sym_keys(load)...)

            loadDict = Dict{String,Any}()

            nphases = tppm_data["phases"]
            name, nodes = parse_busname(defaults["bus1"])

            loadDict["name"] = defaults["name"]
            loadDict["load_bus"] = find_bus(name, tppm_data)
            loadDict["pd"] = PMs.MultiPhaseVector(parse_array(defaults["kw"] / 1e3, nodes, nphases))
            loadDict["qd"] = PMs.MultiPhaseVector(parse_array(defaults["kvar"] / 1e3, nodes, nphases))
            loadDict["status"] = convert(Int, defaults["enabled"])

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
            defaults = createCapacitor(shunt["bus1"], shunt["name"]; to_sym_keys(shunt)...)

            shuntDict = Dict{String,Any}()

            nphases = tppm_data["phases"]
            name, nodes = parse_busname(defaults["bus1"])

            Zbase = (tppm_data["basekv"] / sqrt(3.0))^2 * nphases / tppm_data["baseMVA"]  # Use single-phase base impedance for each phase
            Gcap = -Zbase * sum(defaults["kvar"]) / (nphases * 1e3 * (tppm_data["basekv"] / sqrt(3.0))^2)

            shuntDict["shunt_bus"] = find_bus(name, tppm_data)
            shuntDict["name"] = defaults["name"]
            shuntDict["gs"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))  # TODO:
            shuntDict["bs"] = PMs.MultiPhaseVector(parse_array(Gcap, nodes, nphases))
            shuntDict["status"] = convert(Int, defaults["enabled"])
            shuntDict["index"] = length(tppm_data["shunt"]) + 1

            used = ["bus1", "phases", "name"]
            PMs.import_remaining!(shuntDict, defaults, import_all; exclude=used)

            push!(tppm_data["shunt"], shuntDict)
        end
    end


    if haskey(dss_data, "reactor")
        for shunt in dss_data["reactor"]
            if !haskey(shunt, "bus2")
                defaults = createReactor(find_bus(shunt["bus1"]), shunt["name"]; to_sym_keys(shunt)...)

                shuntDict = Dict{String,Any}()

                nphases = tppm_data["phases"]
                name, nodes = parse_busname(defaults["bus1"])

                Zbase = (tppm_data["basekv"] / sqrt(3.0))^2 * nphases / tppm_data["baseMVA"]  # Use single-phase base impedance for each phase
                Gcap = Zbase * sum(defaults["kvar"]) / (nphases * 1e3 * (tppm_data["basekv"] / sqrt(3.0))^2)

                shuntDict["shunt_bus"] = find_bus(name, tppm_data)
                shuntDict["name"] = defaults["name"]
                shuntDict["gs"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))  # TODO:
                shuntDict["bs"] = PMs.MultiPhaseVector(parse_array(Gcap, nodes, nphases))
                shuntDict["status"] = convert(Int, defaults["enabled"])

                shuntDict["index"] = length(tppm_data["shunt"]) + 1

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
    if haskey(dss_data, "circuit")
        circuit = dss_data["circuit"][1]
        defaults = createVSource("sourcebus", "sourcebus"; to_sym_keys(circuit)...)

        genDict = Dict{String,Any}()

        nphases = tppm_data["phases"]
        name, nodes = parse_busname(defaults["bus1"])

        genDict["gen_bus"] = find_bus(name, tppm_data)
        genDict["name"] = defaults["name"]
        genDict["gen_status"] = convert(Int, defaults["enabled"])

        # TODO: populate with VSOURCE properties
        genDict["pg"] = PMs.MultiPhaseVector(parse_array( 0.0, nodes, nphases))
        genDict["qg"] = PMs.MultiPhaseVector(parse_array( 0.0, nodes, nphases))

        genDict["qmin"] = PMs.MultiPhaseVector(parse_array(-Inf, nodes, nphases))
        genDict["qmax"] = PMs.MultiPhaseVector(parse_array( Inf, nodes, nphases))

        genDict["pmax"] = PMs.MultiPhaseVector(parse_array( Inf, nodes, nphases))
        genDict["pmin"] = PMs.MultiPhaseVector(parse_array(-Inf, nodes, nphases))


        genDict["model"] = 2
        genDict["startup"] = 0.0
        genDict["shutdown"] = 0.0
        genDict["ncost"] = 3
        genDict["cost"] = [0.0, 1.0, 0.0]

        genDict["index"] = length(tppm_data["gen"]) + 1

        used = ["name", "phases", "bus1"]
        PMs.import_remaining!(genDict, defaults, import_all; exclude=used)

        push!(tppm_data["gen"], genDict)
    else
        error(LOGGER, "sourcebus, as created by circuit object, is required in opendss")
    end

    if haskey(dss_data, "generator")
        for gen in dss_data["generator"]
            defaults = createGenerator(gen["bus1"], gen["name"]; to_sym_keys(gen)...)

            genDict = Dict{String,Any}()

            nphases = tppm_data["phases"]
            name, nodes = parse_busname(defaults["bus1"])

            genDict["gen_bus"] = find_bus(name, tppm_data)
            genDict["name"] = defaults["name"]
            genDict["gen_status"] = convert(Int, defaults["enabled"])
            genDict["pg"] = PMs.MultiPhaseVector(parse_array(defaults["kw"] / (1e3 * nphases), nodes, nphases))
            genDict["qg"] = PMs.MultiPhaseVector(parse_array(defaults["kvar"] / (1e3 * nphases), nodes, nphases))
            genDict["vg"] = PMs.MultiPhaseVector(parse_array(defaults["kv"] / tppm_data["basekv"], nodes, nphases))

            genDict["qmin"] = PMs.MultiPhaseVector(parse_array(defaults["minkvar"] / (1e3 * nphases), nodes, nphases))
            genDict["qmax"] = PMs.MultiPhaseVector(parse_array(defaults["maxkvar"] / (1e3 * nphases), nodes, nphases))

            genDict["apf"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))

            genDict["pmax"] = genDict["pg"]  # Assumes generator is at rated power
            genDict["pmin"] = 0.3 * genDict["pg"]  # 30% of pmax

            genDict["pc1"] = genDict["pmax"]
            genDict["pc2"] = genDict["pmin"]
            genDict["qc1min"] = genDict["qmin"]
            genDict["qc1max"] = genDict["qmax"]
            genDict["qc2min"] = genDict["qmin"]
            genDict["qc2max"] = genDict["qmax"]

            # For distributed generation ramp rates are not usually an issue
            # and they are not supported in OpenDSS
            genDict["ramp_agc"] = genDict["pmax"]

            genDict["ramp_q"] = PMs.MultiPhaseVector(parse_array(max.(abs.(genDict["qmin"].values), abs.(genDict["qmax"].values)), nodes, nphases))
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

            used = ["name", "phases", "bus1"]
            PMs.import_remaining!(genDict, defaults, import_all; exclude=used)

            push!(tppm_data["gen"], genDict)
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
            defaults = createLine(line["bus1"], line["bus2"], line["name"]; to_sym_keys(line)...)
            linecode = createLinecode(defaults["linecode"]; to_sym_keys(get_linecode(dss_data, defaults["linecode"]))...)
            merge!(defaults, linecode)

            bf, nodes = parse_busname(defaults["bus1"])
            bt = parse_busname(defaults["bus2"])[1]

            branchDict = Dict{String,Any}()

            nphases = tppm_data["phases"]

            branchDict["name"] = defaults["name"]

            branchDict["f_bus"] = find_bus(bf, tppm_data)
            branchDict["t_bus"] = find_bus(bt, tppm_data)

            # TODO: Handle different lengths and units
            # Set length to 1.0 and units to none for now
            if defaults["length"] != 1.0 || defaults["units"] != "none"
                warn(LOGGER, "length=$(defaults["length"]) and/or units=$(defaults["units"]) not supported, setting to 1.0 and \"none\"")
                defaults["length"] = 1.0
                defaults["units"] = "none"
            end

            branchDict["length"] = defaults["length"]

            rmatrix = parse_matrix(defaults["rmatrix"], nodes, nphases)
            xmatrix = parse_matrix(defaults["xmatrix"], nodes, nphases)
            cmatrix = parse_matrix(defaults["cmatrix"], nodes, nphases)
            Zbase = (tppm_data["basekv"] * 1e3)^2 / (tppm_data["baseMVA"] * 1e6)

            branchDict["br_r"] = PMs.MultiPhaseMatrix(rmatrix * defaults["length"] / Zbase / tppm_data["baseMVA"])
            branchDict["br_x"] = PMs.MultiPhaseMatrix(xmatrix * defaults["length"] / Zbase / tppm_data["baseMVA"])

            # CHECK: Do we need to reformulate to use a matrix instead of a vector for g, b?
            branchDict["g_fr"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases) / tppm_data["baseMVA"])
            branchDict["g_to"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases) / tppm_data["baseMVA"])

            if !isdiag(cmatrix)
                info(LOGGER, "Only diagonal elements of cmatrix are used to obtain branch values `b_fr/to`")
            end
            branchDict["b_fr"] = PMs.MultiPhaseVector(diag(Zbase * (2.0 * pi * defaults["basefreq"] * cmatrix * defaults["length"] / 1e9) / 2.0) / tppm_data["baseMVA"] )
            branchDict["b_to"] = PMs.MultiPhaseVector(diag(Zbase * (2.0 * pi * defaults["basefreq"] * cmatrix * defaults["length"] / 1e9) / 2.0) / tppm_data["baseMVA"] )

            # TODO: pick a better value for emergamps
            branchDict["rate_a"] = PMs.MultiPhaseVector(parse_array(defaults["normamps"], nodes, nphases, NaN))
            branchDict["rate_b"] = PMs.MultiPhaseVector(parse_array(defaults["emergamps"], nodes, nphases, NaN))
            branchDict["rate_c"] = PMs.MultiPhaseVector(parse_array(defaults["emergamps"], nodes, nphases, NaN))

            branchDict["tap"] = PMs.MultiPhaseVector(parse_array(1.0, nodes, nphases, 1.0))
            branchDict["shift"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))

            branchDict["br_status"] = convert(Int, defaults["enabled"])

            branchDict["angmin"] = PMs.MultiPhaseVector(parse_array(-60.0, nodes, nphases, -60.0))
            branchDict["angmax"] = PMs.MultiPhaseVector(parse_array( 60.0, nodes, nphases,  60.0))

            branchDict["transformer"] = false

            branchDict["index"] = length(tppm_data["branch"]) + 1

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
    if !haskey(tppm_data, "branch")
        tppm_data["branch"] = []
    end

    if haskey(dss_data, "transformer")
        warn(LOGGER, "transformers are not yet supported, treating like non-transformer lines")
        for transformer in dss_data["transformer"]
            defaults = createTransformer(transformer["name"]; to_sym_keys(transformer)...)

            nphases = tppm_data["phases"]
            windings = defaults["windings"]

            if windings == 2
                transDict = Dict{String,Any}()
                transDict["name"] = defaults["name"]

                f_bus, nodes = parse_busname(defaults["buses"][1])
                t_bus = parse_busname(defaults["buses"][2])[1]

                transDict["f_bus"] = find_bus(f_bus, tppm_data)
                transDict["t_bus"] = find_bus(t_bus, tppm_data)

                transDict["br_r"] = PMs.MultiPhaseMatrix(parse_matrix(diagm(fill(0.2, nphases)), nodes, nphases))
                transDict["br_x"] = PMs.MultiPhaseMatrix(parse_matrix(zeros(nphases, nphases), nodes, nphases))

                transDict["g_fr"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))
                transDict["g_to"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))
                transDict["b_fr"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))
                transDict["b_to"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))

                # CHECK: unit conversion?
                transDict["rate_a"] = PMs.MultiPhaseVector(parse_array(defaults["normhkva"], nodes, nphases, NaN))
                transDict["rate_b"] = PMs.MultiPhaseVector(parse_array(defaults["emerghkva"], nodes, nphases, NaN))
                transDict["rate_c"] = PMs.MultiPhaseVector(parse_array(defaults["emerghkva"], nodes, nphases, NaN))

                transDict["tap"] = PMs.MultiPhaseVector(parse_array(/(defaults["taps"]...), nodes, nphases, 1.0))
                transDict["shift"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))

                transDict["br_status"] = convert(Int, defaults["enabled"])

                transDict["angmin"] = PMs.MultiPhaseVector(parse_array(-60.0, nodes, nphases, -60.0))
                transDict["angmax"] = PMs.MultiPhaseVector(parse_array( 60.0, nodes, nphases,  60.0))

                transDict["transformer"] = true

                transDict["index"] = length(tppm_data["branch"]) + 1

                used = []
                PMs.import_remaining!(transDict, defaults, import_all; exclude=used)

                push!(tppm_data["branch"], transDict)
            else
                warn(LOGGER, "3-winding transformers are not yet supported, treating like two non-transformer lines connected through a starbus")

                starbus = create_starbus(tppm_data, defaults)
                push!(tppm_data["bus"], starbus)

                bus1, nodes1 = parse_busname(defaults["buses"][1])
                bus2, nodes2 = parse_busname(defaults["buses"][2])
                bus3, nodes3 = parse_busname(defaults["buses"][3])

                for (m, (bus, nodes)) in enumerate(zip([bus1, bus2, bus3], [nodes1, nodes2, nodes3]))
                    transDict = Dict{String,Any}()

                    transDict["f_bus"] = find_bus(bus, tppm_data)
                    transDict["t_bus"] = starbus["bus_i"]

                    transDict["name"] = "$(defaults["name"]) winding $m"

                    transDict["br_r"] = PMs.MultiPhaseMatrix(parse_matrix(diagm(fill(0.2, nphases)), nodes, nphases))
                    transDict["br_x"] = PMs.MultiPhaseMatrix(parse_matrix(zeros(nphases, nphases), nodes, nphases))

                    transDict["g_fr"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))
                    transDict["g_to"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))
                    transDict["b_fr"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))
                    transDict["b_to"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))

                    # CHECK: unit conversion?
                    transDict["rate_a"] = PMs.MultiPhaseVector(parse_array(defaults["normhkva"], nodes, nphases, NaN))
                    transDict["rate_b"] = PMs.MultiPhaseVector(parse_array(defaults["emerghkva"], nodes, nphases, NaN))
                    transDict["rate_c"] = PMs.MultiPhaseVector(parse_array(defaults["emerghkva"], nodes, nphases, NaN))

                    transDict["tap"] = PMs.MultiPhaseVector(parse_array(defaults["taps"][m], nodes, nphases, 1.0))
                    transDict["shift"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))

                    transDict["br_status"] = convert(Int, defaults["enabled"])

                    transDict["angmin"] = PMs.MultiPhaseVector(parse_array(-60.0, nodes, nphases, -60.0))
                    transDict["angmax"] = PMs.MultiPhaseVector(parse_array( 60.0, nodes, nphases,  60.0))

                    transDict["transformer"] = true

                    transDict["index"] = length(tppm_data["branch"]) + 1

                    used = []
                    PMs.import_remaining!(transDict, defaults, import_all; exclude=used)

                    push!(tppm_data["branch"], transDict)
                end
            end
        end
    end

    if haskey(dss_data, "reactor")
        warn(LOGGER, "reactors as constant impedance elements is not yet supported, treating like line")
        for reactor in dss_data["reactor"]
            if haskey(reactor, "bus2")
                defaults = createReactor(reactor["bus1"], reactor["name"], reactor["bus2"]; to_sym_keys(reactor)...)

                reactDict = Dict{String,Any}()

                nphases = tppm_data["phases"]

                f_bus, nodes = parse_busname(defaults["bus1"])
                t_bus = parse_busname(defaults["bus2"])[1]

                reactDict["f_bus"] = find_bus(f_bus, tppm_data)
                reactDict["t_bus"] = find_bus(t_bus, tppm_data)

                reactDict["br_r"] = PMs.MultiPhaseMatrix(parse_matrix(diagm(fill(0.2, nphases)), nodes, nphases))
                reactDict["br_x"] = PMs.MultiPhaseMatrix(parse_matrix(zeros(nphases, nphases), nodes, nphases))

                reactDict["g_fr"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))
                reactDict["g_to"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))
                reactDict["b_fr"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))
                reactDict["b_to"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))

                reactDict["rate_a"] = PMs.MultiPhaseVector(parse_array(defaults["normamps"], nodes, nphases, NaN))
                reactDict["rate_b"] = PMs.MultiPhaseVector(parse_array(defaults["emergamps"], nodes, nphases, NaN))
                reactDict["rate_c"] = PMs.MultiPhaseVector(parse_array(defaults["emergamps"], nodes, nphases, NaN))

                reactDict["tap"] = PMs.MultiPhaseVector(parse_array(1.0, nodes, nphases, NaN))
                reactDict["shift"] = PMs.MultiPhaseVector(parse_array(0.0, nodes, nphases))

                reactDict["br_status"] = convert(Int, defaults["enabled"])

                reactDict["angmin"] = PMs.MultiPhaseVector(parse_array(-60.0, nodes, nphases, -60.0))
                reactDict["angmax"] = PMs.MultiPhaseVector(parse_array( 60.0, nodes, nphases,  60.0))

                reactDict["transformer"] = true

                reactDict["index"] = length(tppm_data["branch"]) + 1

                used = []
                PMs.import_remaining!(reactDict, defaults, import_all; exclude=used)

                push!(tppm_data["branch"], reactDict)
            end
        end
    end
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
    return out
end


"Parses a Dict resulting from the parsing of a DSS file into a PowerModels usable format."
function parse_opendss(dss_data::Dict; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1)::Dict
    tppm_data = Dict{String,Any}()

    check_duplicate_components!(dss_data)

    parse_dss_with_dtypes!(dss_data, ["line", "linecode", "load", "generator", "capacitor",
                                      "reactor", "circuit", "transformer"])

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
        tppm_data["pu"] = defaults["pu"]
        tppm_data["phases"] = defaults["phases"]
    else
        error(LOGGER, "Circuit not defined, not a valid circuit!")
    end

    dss2tppm_bus!(tppm_data, dss_data, import_all, vmin, vmax)
    dss2tppm_load!(tppm_data, dss_data, import_all)
    dss2tppm_shunt!(tppm_data, dss_data, import_all)
    dss2tppm_branch!(tppm_data, dss_data, import_all)
    dss2tppm_transformer!(tppm_data, dss_data, import_all)
    dss2tppm_gen!(tppm_data, dss_data, import_all)

    tppm_data["dcline"] = []

    InfrastructureModels.arrays_to_dicts!(tppm_data)

    tppm_data["files"] = dss_data["filename"]

    return tppm_data
end


"Parses a DSS file into a PowerModels usable format."
function parse_opendss(io::IOStream; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1)::Dict
    dss_data = parse_dss(io)

    return parse_opendss(dss_data; import_all=import_all)
end
