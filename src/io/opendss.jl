# OpenDSS parser


""
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


"""
    split_transformer_buses!(transformer)

Normalizes the names of buses in `transformer` to "buses". OpenDSS allows for
use of 3 "bus" keywords, or a "buses" keyword, which is an array of bus names.
"""
function split_transformer_buses!(transformer::Dict)
    if haskey(transformer, "buses")
        n1, n2 = parse_array(String, transformer["buses"])
        transformer["bus"] = n1
        transformer["bus_2"] = n2
    elseif haskey(transformer, "bus")
        transformer["buses"] = [tranformer["bus"], transformer["bus_2"]]
        if haskey(transformer, "bus_3")
            push!(transformer["buses"], transformer["bus_3"])
        end
    end
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
                    split_transformer_buses!(compObj)
                end
                if haskey(compObj, "bus2")
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
    phase_on_off(data, nodes, phases=3, value=0)

Disables phase(s) that are indicated by Bool array `nodes` as being on/off by
setting the corresponding phase(s) in `data` to `value`. Currently ignores
ground.
"""
function phase_on_off(data::T, nodes::Array, phases::Int=3, value=0)::T where T
    for (i, node) in enumerate(nodes[1:phases])
        if !node
            if length(size(data)) == 2
                data[i, :] = data[:, i] = value
            else
                data[i] = value
            end
        end
    end

    return data
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

    for (n, (bus, nodes)) in enumerate(buses)
        busDict = Dict{String,Any}()

        nphases = tppm_data["phases"]

        busDict["bus_i"] = n
        busDict["index"] = n
        busDict["name"] = bus

        busDict["bus_type"] = bus == "sourcebus" ? 3 : 1

        busDict["vm"] = PMs.MultiPhaseVector(phase_on_off(ones(nphases), nodes, nphases))
        busDict["va"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))

        busDict["vmin"] = PMs.MultiPhaseVector(phase_on_off(fill(vmin, nphases), nodes, nphases))
        busDict["vmax"] = PMs.MultiPhaseVector(phase_on_off(fill(vmax, nphases), nodes, nphases))

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
            loadDict["pd"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["kw"], nphases), nodes, nphases))
            loadDict["qd"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["kvar"], nphases), nodes, nphases))
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
            shuntDict["gs"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))  # TODO:
            shuntDict["bs"] = PMs.MultiPhaseVector(phase_on_off(fill(Gcap, nphases), nodes, nphases))
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
                defaults = createReactor(shunt["bus1"], shunt["name"]; to_sym_keys(shunt)...)

                shuntDict = Dict{String,Any}()

                nphases = tppm_data["phases"]
                name, nodes = parse_busname(defaults["bus1"])

                Zbase = (tppm_data["basekv"] / sqrt(3.0))^2 * nphases / tppm_data["baseMVA"]  # Use single-phase base impedance for each phase
                Gcap = Zbase * sum(defaults["kvar"]) / (nphases * 1e3 * (tppm_data["basekv"] / sqrt(3.0))^2)

                shuntDict["shunt_bus"] = find_bus(name, tppm_data)
                shuntDict["name"] = defaults["name"]
                shuntDict["gs"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))  # TODO:
                shuntDict["bs"] = PMs.MultiPhaseVector(phase_on_off(fill(Gcap, nphases), nodes, nphases))
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


    if haskey(dss_data, "generator")
        for gen in dss_data["generator"]
            defaults = createGenerator(gen["bus1"], gen["name"]; to_sym_keys(gen)...)

            genDict = Dict{String,Any}()

            nphases = tppm_data["phases"]
            name, nodes = parse_busname(defaults["bus1"])

            genDict["gen_bus"] = find_bus(name, tppm_data)
            genDict["name"] = defaults["name"]
            genDict["gen_status"] = convert(Int, defaults["enabled"])
            genDict["pg"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["kw"] / (1e3 * nphases), nphases), nodes, nphases))
            genDict["qg"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["kvar"] / (1e3 * nphases), nphases), nodes, nphases))
            genDict["vg"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["kv"] / tppm_data["basekv"], nphases), nodes, nphases))

            genDict["qmin"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["minkvar"] / (1e3 * nphases), nphases), nodes, nphases))
            genDict["qmax"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["maxkvar"] / (1e3 * nphases), nphases), nodes, nphases))

            genDict["apf"] = PMs.MultiPhaseVector(zeros(nphases))

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

            genDict["ramp_q"] = PMs.MultiPhaseVector(phase_on_off(max.(abs.(genDict["qmin"].values), abs.(genDict["qmax"].values)), nodes, nphases))
            genDict["ramp_10"] = genDict["pmax"]
            genDict["ramp_30"] = genDict["pmax"]

            genDict["control_model"] = defaults["model"]

            # if PV generator mode convert attached bus to PV bus
            if genDict["control_model"] == 3
                tppm_data["bus"][genDict["gen_bus"]]["bus_type"] = 2
            end

            genDict["model"] = PMs.MultiPhaseVector(phase_on_off(fill(2, nphases), nodes, nphases))
            genDict["startup"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))
            genDict["shutdown"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))
            genDict["ncost"] = PMs.MultiPhaseVector(phase_on_off(fill(3, nphases), nodes, nphases))
            genDict["cost"] = PMs.MultiPhaseVector(phase_on_off([[0.0, 1.0, 0.0] for n in 1:nphases], nodes, nphases))

            genDict["index"] = length(tppm_data["gen"]) + 1

            used = ["name", "phases", "bus1"]
            PMs.import_remaining!(genDict, defaults, import_all; exclude=used)

            push!(tppm_data["gen"], genDict)
        end

        if haskey(dss_data, "vsource")
            for vsource in dss_data["vsource"]
                # TODO: Add VSOURCE AS GENERATOR
            end
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

            linecode = createLinecode(defaults["linecode"]; to_sym_keys(get_linecode(tppm_data, defaults["linecode"]))...)
            for k in ("rmatrix", "cmatrix", "xmatrix")
                defaults[k] = linecode[k]
            end

            bf, nodes = parse_busname(defaults["bus1"])
            bt = parse_busname(defaults["bus2"])[1]

            branchDict = Dict{String,Any}()

            nphases = tppm_data["phases"]

            branchDict["name"] = defaults["name"]

            branchDict["f_bus"] = find_bus(bf, tppm_data)
            branchDict["t_bus"] = find_bus(bt, tppm_data)

            defaults["length"] = 1.0
            branchDict["length"] = defaults["length"]

            c_nf = parse_matrix(Float64, defaults["cmatrix"], nphases)
            Zbase = (tppm_data["basekv"] / sqrt(3.0))^2 * nphases / tppm_data["baseMVA"]

            branchDict["br_r"] = PMs.MultiPhaseMatrix(phase_on_off(parse_matrix(Float64, defaults["rmatrix"], nphases) .* defaults["length"] ./ Zbase, nodes, nphases))
            branchDict["br_x"] = PMs.MultiPhaseMatrix(phase_on_off(parse_matrix(Float64, defaults["xmatrix"], nphases) .* defaults["length"] ./ Zbase, nodes, nphases))

            # CHECK: Do we need to reformulate to use a matrix instead of a vector for g, b?
            branchDict["g_fr"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))  # TODO: need to derive
            branchDict["g_to"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))  # TODO: need to derive

            branchDict["b_fr"] = PMs.MultiPhaseVector(phase_on_off(diag(Zbase .* (2.0 .* pi .* defaults["basefreq"] .* c_nf .* defaults["length"] ./ 1e9) ./ 2.0), nodes, nphases))
            branchDict["b_to"] = PMs.MultiPhaseVector(phase_on_off(diag(Zbase .* (2.0 .* pi .* defaults["basefreq"] .* c_nf .* defaults["length"] ./ 1e9) ./ 2.0), nodes, nphases))

            # TODO: pick a better value for emergamps
            branchDict["rate_a"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["normamps"], nphases), nodes, nphases, NaN))
            branchDict["rate_b"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["emergamps"], nphases), nodes, nphases, NaN))
            branchDict["rate_c"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["emergamps"], nphases), nodes, nphases, NaN))

            branchDict["tap"] = PMs.MultiPhaseVector(phase_on_off(ones(nphases), nodes, nphases, NaN))
            branchDict["shift"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))

            branchDict["br_status"] = convert(Int, defaults["enabled"])

            branchDict["angmin"] = PMs.MultiPhaseVector(phase_on_off(fill(-60.0, nphases), nodes, nphases, -60.0))
            branchDict["angmax"] = PMs.MultiPhaseVector(phase_on_off(fill( 60.0, nphases), nodes, nphases,  60.0))

            branchDict["transformer"] = false

            branchDict["index"] = length(tppm_data["branch"]) + 1

            used = ["name", "phases", "bus1", "bus2", "rmatrix", "xmatrix"]
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
            split_transformer_buses!(transformer)
            defaults = createTransformer(transformer["buses"], transformer["name"]; to_sym_keys(transformer)...)

            transDict = Dict{String,Any}()

            nphases = tppm_data["phases"]

            transDict["name"] = defaults["name"]


            f_bus, nodes = parse_busname(defaults["buses"][1])
            t_bus = parse_busname(defaults["buses"][2])[1]

            transDict["f_bus"] = find_bus(f_bus, tppm_data)
            transDict["t_bus"] = find_bus(t_bus, tppm_data)

            transDict["br_r"] = PMs.MultiPhaseMatrix(phase_on_off(diagm(fill(0.2, nphases)), nodes, nphases))
            transDict["br_x"] = PMs.MultiPhaseMatrix(phase_on_off(zeros(nphases, nphases), nodes, nphases))

            transDict["g_fr"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))
            transDict["g_to"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))
            transDict["b_fr"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))
            transDict["b_to"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))

            # CHECK: unit conversion?
            transDict["rate_a"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["normhkva"], nphases), nodes, nphases, NaN))
            transDict["rate_b"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["emerghkva"], nphases), nodes, nphases, NaN))
            transDict["rate_c"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["emerghkva"], nphases), nodes, nphases, NaN))

            transDict["tap"] = PMs.MultiPhaseVector(phase_on_off(ones(nphases), nodes, nphases, NaN))
            transDict["shift"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))

            transDict["br_status"] = convert(Int, defaults["enabled"])

            transDict["angmin"] = PMs.MultiPhaseVector(phase_on_off(fill(-60.0, nphases), nodes, nphases, -60.0))
            transDict["angmax"] = PMs.MultiPhaseVector(phase_on_off(fill( 60.0, nphases), nodes, nphases,  60.0))

            transDict["transformer"] = true

            transDict["index"] = length(tppm_data["branch"]) + 1

            used = []
            PMs.import_remaining!(transDict, defaults, import_all; exclude=used)

            push!(tppm_data["branch"], transDict)
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

                reactDict["br_r"] = PMs.MultiPhaseMatrix(phase_on_off(diagm(fill(0.2, nphases)), nodes, nphases))
                reactDict["br_x"] = PMs.MultiPhaseMatrix(phase_on_off(zeros(nphases, nphases), nodes, nphases))

                reactDict["g_fr"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))
                reactDict["g_to"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))
                reactDict["b_fr"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))
                reactDict["b_to"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))

                reactDict["rate_a"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["normamps"], nphases), nodes, nphases, NaN))
                reactDict["rate_b"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["emergamps"], nphases), nodes, nphases, NaN))
                reactDict["rate_c"] = PMs.MultiPhaseVector(phase_on_off(fill(defaults["emergamps"], nphases), nodes, nphases, NaN))

                reactDict["tap"] = PMs.MultiPhaseVector(phase_on_off(ones(nphases), nodes, nphases, NaN))
                reactDict["shift"] = PMs.MultiPhaseVector(phase_on_off(zeros(nphases), nodes, nphases))

                reactDict["br_status"] = convert(Int, defaults["enabled"])

                reactDict["angmin"] = PMs.MultiPhaseVector(phase_on_off(fill(-60.0, nphases), nodes, nphases, -60.0))
                reactDict["angmax"] = PMs.MultiPhaseVector(phase_on_off(fill( 60.0, nphases), nodes, nphases,  60.0))

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
function parse_opendss(filename::String; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1)::Dict
    dss_data = parse_dss(filename)

    return parse_opendss(dss_data; import_all=import_all)
end
