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
    if !haskey(tppm_data, "branch")
        tppm_data["branch"] = []
    end

    if haskey(dss_data, "transformer")
        warn(LOGGER, "transformers are not yet supported, treating like non-transformer lines")
        for transformer in dss_data["transformer"]
            if haskey(transformer, "like")
                transformer = merge(find_component(dss_data, transformer["like"], "transformer"), transformer)
            end

            defaults = createTransformer(transformer["name"]; to_sym_keys(transformer)...)

            nconductors = tppm_data["conductors"]
            windings = defaults["windings"]

            if windings == 2
                transDict = Dict{String,Any}()
                transDict["name"] = defaults["name"]

                f_bus, nodes = parse_busname(defaults["buses"][1])
                t_bus = parse_busname(defaults["buses"][2])[1]

                transDict["f_bus"] = find_bus(f_bus, tppm_data)
                transDict["t_bus"] = find_bus(t_bus, tppm_data)

                transDict["br_r"] = PMs.MultiConductorMatrix(parse_matrix(diagm(0 => fill(0.2, nconductors)), nodes, nconductors))
                transDict["br_x"] = PMs.MultiConductorMatrix(parse_matrix(zeros(nconductors, nconductors), nodes, nconductors))

                transDict["g_fr"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))
                transDict["g_to"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))
                transDict["b_fr"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))
                transDict["b_to"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))

                # CHECK: unit conversion?
                transDict["rate_a"] = PMs.MultiConductorVector(parse_array(defaults["normhkva"], nodes, nconductors))
                transDict["rate_b"] = PMs.MultiConductorVector(parse_array(defaults["emerghkva"], nodes, nconductors))
                transDict["rate_c"] = PMs.MultiConductorVector(parse_array(defaults["emerghkva"], nodes, nconductors))

                transDict["tap"] = PMs.MultiConductorVector(parse_array(/(defaults["taps"]...), nodes, nconductors, 1.0))
                transDict["shift"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))

                transDict["br_status"] = convert(Int, defaults["enabled"])

                transDict["angmin"] = PMs.MultiConductorVector(parse_array(-60.0, nodes, nconductors, -60.0))
                transDict["angmax"] = PMs.MultiConductorVector(parse_array( 60.0, nodes, nconductors,  60.0))

                transDict["transformer"] = true
                transDict["switch"] = false

                transDict["index"] = length(tppm_data["branch"]) + 1

                nodes = .+([parse_busname(defaults["buses"][n])[2] for n in length(defaults["buses"])]...)
                transDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
                transDict["source_id"] = "transformer.$(defaults["name"])"

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

                    transDict["br_r"] = PMs.MultiConductorMatrix(parse_matrix(diagm(0 => fill(0.2, nconductors)), nodes, nconductors))
                    transDict["br_x"] = PMs.MultiConductorMatrix(parse_matrix(zeros(nconductors, nconductors), nodes, nconductors))

                    transDict["g_fr"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))
                    transDict["g_to"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))
                    transDict["b_fr"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))
                    transDict["b_to"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))

                    # CHECK: unit conversion?
                    transDict["rate_a"] = PMs.MultiConductorVector(parse_array(defaults["normhkva"], nodes, nconductors))
                    transDict["rate_b"] = PMs.MultiConductorVector(parse_array(defaults["emerghkva"], nodes, nconductors))
                    transDict["rate_c"] = PMs.MultiConductorVector(parse_array(defaults["emerghkva"], nodes, nconductors))

                    transDict["tap"] = PMs.MultiConductorVector(parse_array(defaults["taps"][m], nodes, nconductors, 1.0))
                    transDict["shift"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))

                    transDict["br_status"] = convert(Int, defaults["enabled"])

                    transDict["angmin"] = PMs.MultiConductorVector(parse_array(-60.0, nodes, nconductors, -60.0))
                    transDict["angmax"] = PMs.MultiConductorVector(parse_array( 60.0, nodes, nconductors,  60.0))

                    transDict["transformer"] = true
                    transDict["switch"] = false

                    transDict["index"] = length(tppm_data["branch"]) + 1

                    nodes = .+([parse_busname(defaults["buses"][n])[2] for n in length(defaults["buses"])]...)
                    transDict["active_phases"] = [n for n in 1:nconductors if nodes[n] > 0]
                    transDict["source_id"] = "transformer.$(defaults["name"])"

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

    tppm_data["files"] = dss_data["filename"]

    return tppm_data
end


"Parses a DSS file into a PowerModels usable format."
function parse_opendss(io::IOStream; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1)::Dict
    dss_data = parse_dss(io)

    return parse_opendss(dss_data; import_all=import_all)
end
