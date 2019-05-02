# OpenDSS parser
import LinearAlgebra: isdiag, diag, pinv

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
        Memento.error(LOGGER, "dss_data has no branches!")
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
    Memento.warn(LOGGER, "Could not find $compType \"$name\"")
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
        Memento.error(LOGGER, "cannot find connected bus with id \"$busname\"")
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
                Memento.warn(LOGGER, "Load has kv=$kv, not the expected kv=$(expected_kv). Results may not match OpenDSS")
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
            Memento.warn(LOGGER, "Converting PVSystem \"$(pv["name"])\" into generator with limits determined by OpenDSS property 'kVA'")

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
                Memento.warn(LOGGER, "basefreq=$(line["basefreq"]) on line $(line["name"]) does not match circuit basefreq=$(tppm_data["basefreq"])")
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

            rmatrix = parse_matrix(defaults["rmatrix"], nodes, nconductors)
            xmatrix = parse_matrix(defaults["xmatrix"], nodes, nconductors)
            cmatrix = parse_matrix(defaults["cmatrix"], nodes, nconductors)
            Zbase = (tppm_data["basekv"] / sqrt(3))^2 * nconductors / (tppm_data["baseMVA"])

            Zbase = Zbase/3
            # The factor 3 here is needed to convert from a voltage base
            # in line-to-line (LL) to a voltage base in line-to-neutral (LN).
            # V_LL = √3*V_LN
            # Zbase_new = Zbase_old*(Vbase_new/Vbase_old)^2 = Zbase_old*(1/√3)^2
            # In the parser, LL voltage base is used for per unit conversion.
            # However, in the mathematical model, the voltage magnitude per phase
            # is fixed at 1. So implicitly, we later on state that the voltage base
            # is actually in LN. We compensate here for that.

            branchDict["br_r"] = PMs.MultiConductorMatrix(rmatrix * defaults["length"] / Zbase)
            branchDict["br_x"] = PMs.MultiConductorMatrix(xmatrix * defaults["length"] / Zbase)

            # CHECK: Do we need to reformulate to use a matrix instead of a vector for g, b?
            branchDict["g_fr"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))
            branchDict["g_to"] = PMs.MultiConductorVector(parse_array(0.0, nodes, nconductors))

            if !isdiag(cmatrix)
                Memento.info(LOGGER, "Only diagonal elements of cmatrix are used to obtain branch values `b_fr/to`")
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
        if !haskey(tppm_data, "trans_comp")
            tppm_data["trans_comp"] = Array{Any,1}()
        end
        for transformer in dss_data["transformer"]
            if haskey(transformer, "like")
                transformer = merge(find_component(dss_data, transformer["like"], "transformer"), transformer)
            end

            defaults = createTransformer(transformer["name"]; to_sym_keys(transformer)...)

            nconductors = tppm_data["conductors"]
            nrw = defaults["windings"]
            prop_suffix_w = ["", ["_$w" for w in 2:nrw]...]
            if nrw>3
                # All of the code is compatible with any number of windings,
                # except for the parsing of the loss model (the pair-wise reactance)
                Memento.error(LOGGER, "For now parsing of xscarray is not supported. At most 3 windings are allowed, not $nrw.")
            end

            transDict = Dict{String,Any}()
            transDict["name"] = defaults["name"]
            transDict["source_id"] = "transformer.$(defaults["name"])"
            transDict["active_phases"] = [1, 2, 3]
            transDict["buses"] = Array{Int, 1}(undef, nrw)
            for i in 1:nrw
                bnstr = defaults["buses"][i]
                bus, nodes = parse_busname(bnstr)
                nodes_0123 = [true true true true]
                nodes_123 = [true true true false]
                if !(nodes==nodes_0123 || nodes==nodes_123)
                    Memento.warn(LOGGER, "Only three-phase transformers are supported. The bus specification $bnstr is treated as $bus instead.")
                end
                transDict["buses"][i] = find_bus(bus, tppm_data)
            end

            # voltage and power ratings
            #transDict["vnom_kv"] = defaults["kvs"]
            #transDict["snom_kva"] = defaults["kvas"]
            transDict["rate_a"] = [PMs.MultiConductorVector(ones(nconductors))*defaults["normhkva"] for i in 1:nrw]
            transDict["rate_b"] = [PMs.MultiConductorVector(ones(nconductors))*defaults["normhkva"] for i in 1:nrw]
            transDict["rate_c"] = [PMs.MultiConductorVector(ones(nconductors))*defaults["emerghkva"] for i  in 1:nrw]
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
            transDict["tm"] = [PMs.MultiConductorVector(ones(Float64,3))*defaults["taps"][i] for i in 1:nrw]
            transDict["tm_min"] = [PMs.MultiConductorVector(ones(Float64,3))*defaults["mintap"] for i in 1:nrw]
            transDict["tm_max"] = [PMs.MultiConductorVector(ones(Float64,3))*defaults["maxtap"] for i in 1:nrw]
            transDict["tm_step"] = [PMs.MultiConductorVector(ones(Int,3))*defaults["numtaps"] for i in 1:nrw]
            transDict["fixed"] = [PMs.MultiConductorVector(ones(Bool,3)) for i in 1:nrw]

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
            transDict["rs"] = Array{PMs.MultiConductorMatrix{Float64}, 1}(undef, nrw)
            transDict["gsh"] = Array{PMs.MultiConductorMatrix{Float64}, 1}(undef, nrw)
            transDict["bsh"] = Array{PMs.MultiConductorMatrix{Float64}, 1}(undef, nrw)
            # deal with dual definition of rs in defaults dict
            # if non-zero, pick "%rs", else go for the per winding spec
            rs_alt = [defaults["%r$suffix"] for suffix in prop_suffix_w]
            rs = all(defaults["%rs"].==0) ? rs_alt : defaults["%rs"]
            for w in 1:nrw
                zs_w_p = rs[w]/100*zbase
                Zs_w = pos_to_abc(zs_w_p)
                #TODO handle %loadloss property
                # Problem is that for the two-winding case, both %loadloss
                # and %rs map to values for the winding series resistance.
                # Both are set to a (possibly default) value, so it is
                # impossible to know at this point which was actually specified
                # by the user. %rs is more general, so %loadloss is not
                # supported for now.
                Memento.warn(LOGGER, "The %loadloss property is ignored for now.")
                #TODO handle neutral impedance
                # neutral impedance is ignored for now; all transformers are
                # grounded (that is, those with a wye and zig-zag winding).
                Memento.warn(LOGGER, "The neutral impedance, (rg and xg properties), is ignored; the neutral (for wye and zig-zag windings) is connected directly to the ground.")
                transDict["rs"][w] = PMs.MultiConductorMatrix(real.(Zs_w))
                # shunt elements are added at second winding
                if w==2
                    ysh_w_p = (defaults["%noloadloss"]-im*defaults["%imag"])/100/zbase
                    Ysh_w = pos_to_abc(ysh_w_p)
                    transDict["gsh"][w] = PMs.MultiConductorMatrix(real.(Ysh_w))
                    transDict["bsh"][w] = PMs.MultiConductorMatrix(imag.(Ysh_w))
                else
                    transDict["gsh"][w] = PMs.MultiConductorMatrix(zeros(Float64, 3, 3))
                    transDict["bsh"][w] = PMs.MultiConductorMatrix(zeros(Float64, 3, 3))
                end
            end
            transDict["xs"] = Dict{String, PMs.MultiConductorMatrix{Float64}}()

            Zsc = Dict{Tuple{Int,Int}, Complex}()
            if nrw==2
                xs_map = Dict("xhl"=>(1,2))
            elseif nrw==3
                xs_map = Dict("xhl"=>(1,2), "xht"=>(1,3), "xlt"=>(2,3))
            end
            for (k,v) in xs_map
                Zsc[(v)] = im*defaults[k]/100*zbase
            end
            Zbr = sc2br_impedance(Zsc)
            for (k,zs_ij_p) in Zbr
                Zs_ij = pos_to_abc(zs_ij_p)
                transDict["xs"]["$(k[1])-$(k[2])"] = PMs.MultiConductorMatrix(imag.(Zs_ij))
            end

            push!(tppm_data["trans_comp"], transDict)
        end
    end
end


"""
Converts a set of short-circuit tests to an equivalent reactance network.
Reference:
R. C. Dugan, “A perspective on transformer modeling for distribution system analysis,”
in 2003 IEEE Power Engineering Society General Meeting (IEEE Cat. No.03CH37491), 2003, vol. 1, pp. 114-119 Vol. 1.
"""
function sc2br_impedance(Zsc)
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
                    Memento.error(LOGGER, "Short-circuit impedance between winding $i and $j is missing.")
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

function dss2tppm_reactor!(tppm_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(tppm_data, "branch")
        tppm_data["branch"] = []
    end

    if haskey(dss_data, "reactor")
        Memento.warn(LOGGER, "reactors as constant impedance elements is not yet supported, treating like line")
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
This method was updated to include connected transformers as well. It know
has to occur after the call to InfrastructureModels.arrays_to_dicts, so the code
was adjusted to accomodate that.
"""
function adjust_sourcegen_bounds!(tppm_data)
    emergamps = Array{Float64,1}([0.0])
    #sourcebus_n = find_bus("sourcebus", tppm_data)
    sourcebus_n = [bus["index"] for (_,bus) in tppm_data["bus"] if haskey(bus, "name") && bus["name"]=="sourcebus"][1]
    for (_,line) in tppm_data["branch"]
        if line["f_bus"] == sourcebus_n || line["t_bus"] == sourcebus_n
            append!(emergamps, line["rate_b"].values)
        end
    end
    if haskey(tppm_data, "trans")
        for (_,trans) in tppm_data["trans"]
            if trans["f_bus"] == sourcebus_n || trans["t_bus"] == sourcebus_n
                append!(emergamps, trans["rate_b"].values)
            end
        end
    end

    bound = sum(emergamps)

    tppm_data["gen"]["1"]["pmin"] = PMs.MultiConductorVector(fill(-bound, size(tppm_data["gen"]["1"]["pmin"])))
    tppm_data["gen"]["1"]["pmax"] = PMs.MultiConductorVector(fill( bound, size(tppm_data["gen"]["1"]["pmin"])))
    tppm_data["gen"]["1"]["qmin"] = PMs.MultiConductorVector(fill(-bound, size(tppm_data["gen"]["1"]["pmin"])))
    tppm_data["gen"]["1"]["qmax"] = PMs.MultiConductorVector(fill( bound, size(tppm_data["gen"]["1"]["pmin"])))
end


"""

    function decompose_transformers!(tppm_data)

Replaces complex transformers with a composition of ideal transformers and branches
which model losses. New buses (virtual, no physical meaning) are added.
"""
function decompose_transformers!(tppm_data; import_all::Bool=false)
    if !haskey(tppm_data, "trans")
        tppm_data["trans"] = Dict{String, Any}()
    end
    ncnds = tppm_data["conductors"]
    for (tr_id, trans) in tppm_data["trans_comp"]
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
            push_dict_ret_key!(tppm_data["trans"], trans_dict)
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
            vbus_tr = create_vbus!(tppm_data, basekv=1.0, name="tr$(tr_id)_w$(w)_b1")
            trans_dict["t_bus"] = vbus_tr["index"]
            append!(bus_reduce, vbus_tr["index"])
            # convert to baseMVA, because this is not done per_unit now)
            trans_dict["rate_a"] = trans["rate_a"][w]/tppm_data["baseMVA"]
            trans_dict["rate_b"] = trans["rate_b"][w]/tppm_data["baseMVA"]
            trans_dict["rate_c"] = trans["rate_c"][w]/tppm_data["baseMVA"]
            # tap settings
            trans_dict["tm"] = trans["tm"][w]
            trans_dict["fixed"] = trans["fixed"][w]
            trans_dict["tm_max"] = trans["tm_max"][w]
            trans_dict["tm_min"] = trans["tm_min"][w]
            trans_dict["tm_step"] = trans["tm_step"][w]
            # WINDING SERIES RESISTANCE
            # make virtual bus and mark it for reduction
            vbus_br = create_vbus!(tppm_data, basekv=1.0, name="tr$(tr_id)_w$(w)_b2")
            append!(bus_reduce, vbus_br["index"])
            # make virtual branch and mark it for reduction
            br = create_vbranch!(
                tppm_data, vbus_tr["index"], vbus_br["index"],
                vbase=1.0,
                br_r=trans["rs"][w],
                g_fr=diag(trans["gsh"][w]),
                b_fr=diag(trans["bsh"][w]),
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
                br = create_vbranch!(
                    tppm_data, endnode_id_w[w], endnode_id_w[v],
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
        rm_redundant_pd_elements!(tppm_data, buses=string.(bus_reduce), branches=string.(branch_reduce))
    end
    # remove the trans_comp dict unless import_all is flagged
    if !import_all
        delete!(tppm_data, "trans_comp")
    end
end


"""
This function adds a new bus to the data model and returns its dictionary.
It is virtual in the sense that it does not correspond to a bus in the network,
but is part of the decomposition of the transformer.
"""
function create_vbus!(tppm_data; vmin=0, vmax=Inf, basekv=tppm_data["basekv"], name="", source_id="")
    vbus = Dict{String, Any}("bus_type"=>"1", "name"=>name)
    vbus_id = push_dict_ret_key!(tppm_data["bus"], vbus)
    vbus["bus_i"] = vbus_id
    vbus["source_id"] = source_id
    ncnds = tppm_data["conductors"]
    vbus["vm"] = PMs.MultiConductorVector(ones(Float64, ncnds))
    vbus["va"] = PMs.MultiConductorVector(zeros(Float64, ncnds))
    vbus["vmin"] = PMs.MultiConductorVector(ones(Float64, ncnds))*vmin
    vbus["vmax"] = PMs.MultiConductorVector(ones(Float64, ncnds))*vmax
    vbus["base_kv"] = basekv
    return vbus
end

"""
This function adds a new branch to the data model and returns its dictionary.
It is virtual in the sense that it does not correspond to a branch in the
network, but is part of the decomposition of the transformer.
"""
function create_vbranch!(tppm_data, f_bus::Int, t_bus::Int; name="", source_id="", active_phases=[1, 2, 3], kwargs...)
    ncnd = tppm_data["conductors"]
    kwargs = Dict{Symbol,Any}(kwargs)
    vbase = haskey(kwargs, :vbase) ? kwargs[:vbase] : tppm_data["basekv"]
    # TODO assumes per_unit will be flagged
    sbase = haskey(kwargs, :sbase) ? kwargs[:sbase] : tppm_data["baseMVA"]
    zbase = vbase^2/sbase
    # convert to LN vbase in instead of LL vbase
    zbase *= (1/3)
    vbranch = Dict{String, Any}("f_bus"=>f_bus, "t_bus"=>t_bus, "name"=>name)
    vbranch["active_phases"] = active_phases
    vbranch["source_id"] = "virtual_branch.$name"
    for k in [:br_r, :br_x, :g_fr, :g_to, :b_fr, :b_to]
        if !haskey(kwargs, k)
            if k in [:br_r, :br_x]
                vbranch[string(k)] = PMs.MultiConductorMatrix(zeros(ncnd, ncnd))
            else
                vbranch[string(k)] = PMs.MultiConductorVector(zeros(ncnd))
            end
        else
            if k in [:br_r, :br_x]
                vbranch[string(k)] = kwargs[k]./zbase
            else
                vbranch[string(k)] = kwargs[k].*zbase
            end
        end
    end
    vbranch["angmin"] = -PMs.MultiConductorVector(ones(ncnd))*60
    vbranch["angmax"] = PMs.MultiConductorVector(ones(ncnd))*60
    vbranch["shift"] = PMs.MultiConductorVector(zeros(ncnd))
    vbranch["tap"] = PMs.MultiConductorVector(ones(ncnd))
    vbranch["transformer"] = false
    vbranch["switch"] = false
    vbranch["br_status"] = 1
    for k in [:rate_a, :rate_b, :rate_c]
        if haskey(kwargs, k)
            vbranch[string(k)] = kwargs[k]
        end
    end
    push_dict_ret_key!(tppm_data["branch"], vbranch)
    return vbranch
end


"""
This function appends a component to a component dictionary of a tppm data model.
"""
function push_dict_ret_key!(dict::Dict{String, Any}, v::Dict{String, Any}; assume_no_gaps=false)
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
function rm_redundant_pd_elements!(tppm_data; buses=keys(tppm_data["bus"]), branches=keys(tppm_data["branch"]))
    # temporary dictionary for pi-model shunt elements
    shunts_g = Dict{Int, Any}()
    shunts_b = Dict{Int, Any}()
    for (br_id, br) in tppm_data["branch"]
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
                rm_bus = (f_bus in buses)  ? f_bus :  t_bus
                kp_bus = (rm_bus==f_bus) ? t_bus : f_bus
            elseif is_selfloop
                kp_bus = t_bus
            else
                # nothing to do, go to next branch
                continue
            end
            # move shunts to the bus that will be left
            if !haskey(shunts_g, kp_bus)
                shunts_g[kp_bus] =  PMs.MultiConductorVector(zeros(3))
                shunts_b[kp_bus] =  PMs.MultiConductorVector(zeros(3))
            end
            shunts_g[kp_bus] .+= br["g_fr"]
            shunts_g[kp_bus] .+= br["g_to"]
            shunts_b[kp_bus] .+= br["b_fr"]
            shunts_b[kp_bus] .+= br["b_to"]
            # remove branch from tppm_data
            delete!(tppm_data["branch"], string(br_id))
            if is_shorted && is_reducable
                # remove bus from tppm_data
                delete!(tppm_data["bus"], string(rm_bus))
                # replace bus references in branches
                for (br_id, br) in  tppm_data["branch"]
                    if br["f_bus"] == rm_bus
                        br["f_bus"] = kp_bus
                    end
                    if br["t_bus"] == rm_bus
                        br["t_bus"] = kp_bus
                    end
                end
                # replace bus references in transformers
                for (_, tr) in tppm_data["trans"]
                    if tr["f_bus"] == rm_bus
                        tr["f_bus"] = kp_bus
                    end
                    if tr["t_bus"] == rm_bus
                        tr["t_bus"] = kp_bus
                    end
                end
                # TODO clean up other references to the removed bus
                # like for example loads, generators, ...
                # skipped  for now, not relevant for transformer loss model
                # + LOADS
                # + shunts
                # ...
            end
        elseif f_bus==t_bus
            # this might occur if not all buses and branches are marked for removal
            # a branch in parallel with a removed branch can turn into a self-loop
            # and if that branch is not marked for removal, we end up here
            Memento.error(LOGGER, "Specified set of buses and branches leads to a self-loop.")
        end
    end
    # create shunts for lumped pi-model shunts
    for (bus, shunt_g) in shunts_g
        shunt_b = shunts_b[bus]
        if !all(shunt_g .==0) || !all(shunt_b  .==0)
            Memento.warn(LOGGER, "Pi-model shunt was moved to a bus shunt. Off-diagonals will be discarded in the data model.")
            # The shunts are part of PM, and will be scaled later on by make_per_unit,
            # unlike TPPM level components. The shunts here originate from TPPM level
            # components which were already scaled. Therefore, we have to undo the
            # scaling here to prevent double scaling later on.
            gs = shunt_g/1*tppm_data["baseMVA"]
            bs = shunt_b/1*tppm_data["baseMVA"]
            add_shunt!(tppm_data, bus, gs=gs,  bs=bs)
        end
    end
end


"""
Helper function to add a new shunt. The shunt element is  always inserted at the
internal bus of the second winding in OpenDSS. If one of the branches of the
loss model connected to this bus, has zero impedance (for example, if XHL==0
or XLT==0 or R[3]==0), then this bus might be removed by
rm_redundant_pd_elements!, in which case a new shunt should be inserted at the
remaining bus of the removed branch.
"""
function add_shunt!(tppm_data, bus; gs=PMs.MultiConductorVector(zeros(3)), bs=PMs.MultiConductorVector(zeros(3)), vbase_kv=1, sbase_mva=1)
    # TODO check whether keys are consistent with the actual data model
    shunt_dict = Dict{String, Any}("status"=>1, "shunt_bus"=>bus)
    zbase  = vbase_kv^2/sbase_mva
    shunt_dict["gs"] = gs*zbase
    shunt_dict["bs"] = bs*zbase
    push_dict_ret_key!(tppm_data["shunt"], shunt_dict, assume_no_gaps=false)
end


"""

    function adjust_base!(tppm_data)

Updates the voltage base at each bus, so that the ratios of the voltage bases
across a transformer are consistent with the ratios of voltage ratings of the
windings. Default behaviour is to start at the primary winding of the first
transformer, and to propagate from there. Branches are updated; the impedances
and addmittances are rescaled to be consistent with the new voltage bases.
"""
function adjust_base!(tppm_data; start_at_first_tr_prim=true)
    # initialize arrays etc. for the recursive part
    edges_br = [(br["index"], br["f_bus"], br["t_bus"]) for (br_id_str, br) in tppm_data["branch"]]
    edges_tr = [(tr["index"], tr["f_bus"], tr["t_bus"]) for (tr_id_str, tr) in tppm_data["trans"]]
    edges_br_visited = Dict{Int, Bool}([(edge[1], false) for edge in edges_br])
    edges_tr_visited = Dict{Int, Bool}([(edge[1], false) for edge in edges_tr])
    bus_ids = [parse(Int, x) for x in keys(tppm_data["bus"])]
    nodes_visited = Dict{Int, Bool}([(bus_id, false) for  bus_id in bus_ids])
    # retrieve old voltage bases from connected nodes before starting
    br_basekv_old = Dict([(br["index"], tppm_data["bus"][string(br["f_bus"])]["base_kv"]) for (br_id_str, br) in tppm_data["branch"]])
    # start from the primary of the first transformer
    if start_at_first_tr_prim && haskey(tppm_data, "trans") && haskey(tppm_data["trans"], "1")
        trans_first = tppm_data["trans"]["1"]
        source = trans_first["f_bus"]
        base_kv_new = trans_first["config_fr"]["vm_nom"]
    else
        # start at type 3 bus if present
        buses_3 = [bus["index"] for (bus_id_str, bus) in tppm_data["bus"] if bus["bus_type"]==3]
        if length(buses_3)==0
            Memento.warn(LOGGER, "No bus of type 3 found; selecting random bus instead.")
            source = parse(Int, rand(keys(tppm_data["bus"])))
        else
            source = buses_3[1]
        end
        base_kv_new = tppm_data["basekv"]
        println(source)
        # Only relevant for future per-unit upgrade
        # Impossible to end up here;
        # condition checked before call to adjust_base!
    end
    adjust_base_rec!(tppm_data, source, base_kv_new, nodes_visited, edges_br, edges_br_visited, edges_tr, edges_tr_visited, br_basekv_old)
    if !all(values(nodes_visited))
        Memento.warn(LOGGER, "The network contains buses which are not reachable from the start node for the change of voltage base.")
    end
end


"""
This is the recursive code that goes with adjust_base!; adjust_base!
initializes arrays and other data that is passed along in the calls to this
recursive function. For very large networks, this might have to be rewritten
to not rely on recursion.
"""
function adjust_base_rec!(tppm_data, source::Int, base_kv_new::Float64, nodes_visited, edges_br, edges_br_visited, edges_tr, edges_tr_visited, br_basekv_old)
    source_dict = tppm_data["bus"][string(source)]
    base_kv_prev = source_dict["base_kv"]
    if !(base_kv_prev≈base_kv_new)
        # only possible when meshed; ensure consistency
        if nodes_visited[source]
            Memento.error(LOGGER, "Transformer ratings lead to an inconsistent definition for the voltage base at bus $source.")
        end
        source_dict["base_kv"] = base_kv_new
        # update the connected shunts with the new voltage base
        source_shunts = [shunt for (sh_id_str, shunt) in tppm_data["shunt"] if shunt["shunt_bus"]==source]
        for shunt in source_shunts
            adjust_base_shunt!(tppm_data, shunt["index"], base_kv_prev, base_kv_new)
        end
        source_name = haskey(source_dict, "name") ? source_dict["name"] : ""
        if source_dict["bus_type"]==3
            #TODO is this the desired behaviour, keep SI units for type 3 bus?
            source_dict["vm"] *= base_kv_prev/base_kv_new
            source_dict["vmax"] *= base_kv_prev/base_kv_new
            source_dict["vmin"] *= base_kv_prev/base_kv_new
            Memento.info(LOGGER, "Rescaling vm, vmin and vmax conform with new base_kv at type 3 bus $source($source_name): $base_kv_prev => $base_kv_new")
        else
            Memento.info(LOGGER, "Resetting base_kv at bus $source($source_name): $base_kv_prev => $base_kv_new")
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
                br = tppm_data["branch"]["$br_id"]
                br_name = haskey(br, "name") ? br["name"] : ""
                Memento.info(LOGGER, "Rescaling impedances at branch $br_id($br_name), conform with change of voltage base: $base_kv_branch_prev => $base_kv_new")
                adjust_base_branch!(tppm_data, br_id, base_kv_branch_prev, base_kv_new)
            end
            # follow the edge to the adjacent node and repeat
            adjust_base_rec!(tppm_data, source_new, base_kv_new, nodes_visited, edges_br, edges_br_visited, edges_tr, edges_tr_visited, br_basekv_old)
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
                base_kv_new_tr *= (trans["config_to"]["vm_nom"]/trans["config_fr"]["vm_nom"])
            else
                base_kv_new_tr *= (trans["config_fr"]["vm_nom"]/trans["config_to"]["vm_nom"])
            end
            # follow the edge to the adjacent node and repeat
            adjust_base_rec!(tppm_data, source_new, base_kv_new_tr, nodes_visited, edges_br, edges_br_visited, edges_tr, edges_tr_visited, br_basekv_old)
        end
    end
end


"""
Rescales the parameters of a branch to reflect a change in voltage base.
"""
function adjust_base_branch!(tppm_data, br_id::Int, base_kv_old::Float64, base_kv_new::Float64)
    branch = tppm_data["branch"][string(br_id)]
    zmult = (base_kv_old/base_kv_new)^2
    branch["br_r"] *= zmult
    branch["br_x"] *= zmult
    branch["g_fr"] *= 1/zmult
    branch["b_fr"] *= 1/zmult
    branch["g_to"] *= 1/zmult
    branch["b_to"] *= 1/zmult
end


"""
Rescales the parameters of a shunt to reflect a change in voltage base.
"""
function adjust_base_shunt!(tppm_data, sh_id::Int, base_kv_old::Float64, base_kv_new::Float64)
    shunt = tppm_data["shunt"][string(sh_id)]
    zmult = (base_kv_old/base_kv_new)^2
    shunt["bs"] *= 1/zmult
    shunt["gs"] *= 1/zmult
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
        Memento.warn(LOGGER, "defaultbasefreq is not defined, default for circuit set to 60 Hz")
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
        Memento.error(LOGGER, "Circuit not defined, not a valid circuit!")
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

    tppm_data["dcline"] = []

    InfrastructureModels.arrays_to_dicts!(tppm_data)

    for optional in ["dcline", "load", "shunt", "storage", "pvsystem", "branch"]
        if length(tppm_data[optional]) == 0
            tppm_data[optional] = Dict{String,Any}()
        end
    end

    if haskey(tppm_data, "trans_comp")
        # this has to be done before calling adjust_sourcegen_bounds!
        decompose_transformers!(tppm_data; import_all=import_all)
        adjust_base!(tppm_data)
    else
        tppm_data["trans"] = Dict{String, Any}()
    end

    adjust_sourcegen_bounds!(tppm_data)

    tppm_data["files"] = dss_data["filename"]

    return tppm_data
end


"Parses a DSS file into a PowerModels usable format."
function parse_opendss(io::IOStream; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1)::Dict
    dss_data = parse_dss(io)

    return parse_opendss(dss_data; import_all=import_all)
end
