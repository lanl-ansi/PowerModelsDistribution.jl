# OpenDSS parser
import LinearAlgebra: diagm


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


"Adds nodes as buses to `data_eng` from `data_dss`"
function _dss2eng_bus!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool=false)
    buses = _discover_buses(data_dss)

    for bus in buses
        @assert !startswith(bus, "_virtual") "Bus $bus: identifiers should not start with _virtual"

        eng_obj = create_bus(status=1, bus_type=1)

        if !haskey(data_eng, "bus")
            data_eng["bus"] = Dict{String,Any}()
        end

        data_eng["bus"][bus] = eng_obj
    end
end


"Adds loadshapes to `data_eng` from `data_dss`"
function _dss2eng_loadshape!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool=false)
    for (name, dss_obj) in get(data_dss, "loadshape", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "loadshape")
        defaults = _apply_ordered_properties(_create_loadshape(name; _to_kwargs(dss_obj)...), dss_obj)

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


"""
Adds loads to `data_eng` from `data_dss`

Constant can still be scaled by other settings, fixed cannot
Note that in the current feature set, fixed therefore equals constant

1: Constant P and Q, default
2: Constant Z
3: Constant P and quadratic Q
4: Exponential
5: Constant I
6: Constant P and fixed Q
# 7: Constant P and quadratic Q (i.e., fixed reactance)
# 8: ZIP
"""
function _dss2eng_load!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool, ground_terminal::Int=4)
    for (name, dss_obj) in get(data_dss, "load", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "load")
        defaults = _apply_ordered_properties(_create_load(dss_obj["bus1"], dss_obj["name"]; _to_kwargs(dss_obj)...), dss_obj)

        # parse the model
        model = defaults["model"]

        if model == 3
            Memento.warn(_LOGGER, "$load_name: load model 3 not supported. Treating as model 1.")
            model = 1
        elseif model == 4
            Memento.warn(_LOGGER, "$load_name: load model 4 not supported. Treating as model 1.")
            model = 1
        elseif model == 6
            Memento.warn(_LOGGER, "$load_name: load model 6 identical to model 1 in current feature set. Treating as model 1.")
            model = 1
        elseif model == 7
            Memento.warn(_LOGGER, "$load_name: load model 7 not supported. Treating as model 1.")
            model = 1
        elseif model == 8
            Memento.warn(_LOGGER, "$load_name: load model 8 not supported. Treating as model 1.")
            model = 1
        end

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
        eng_obj = Dict{String,Any}(
            "name" => name,
            "bus" => bus,
            "model" => model,
            "configuration" => conf,
            "connections" => connections,
            "source_id" => "load.$name",
            "status" => convert(Int, defaults["enabled"])
        )

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

        eng_obj["vnom"] = kv

        eng_obj["pd"] = fill(defaults["kw"]/nphases, nphases)
        eng_obj["qd"] = fill(defaults["kvar"]/nphases, nphases)

        if import_all
            _import_all!(eng_obj, defaults, dss_obj["prop_order"])
        end

        if !haskey(data_eng, "load")
            data_eng["load"] = Dict{String,Any}()
        end

        data_eng["load"][name] = eng_obj
    end
end


"Adds capacitors to `data_eng` from `data_dss`"
function _dss2eng_capacitor!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "capacitor", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "capacitor")
        defaults = _apply_ordered_properties(_create_capacitor(dss_obj["bus1"], dss_obj["name"]; _to_kwargs(dss_obj)...), dss_obj)

        nphases = defaults["phases"]

        conn = defaults["conn"]

        bus_name = _parse_busname(defaults["bus1"])[1]
        bus2_name = _parse_busname(defaults["bus2"])[1]

        if bus_name!=bus2_name
            Memento.error(_LOGGER, "Capacitor $(name): bus1 and bus2 should connect to the same bus.")
        end

        f_terminals = _get_conductors_ordered_dm(defaults["bus1"], default=collect(1:nphases))
        if conn=="wye"
            t_terminals = _get_conductors_ordered_dm(defaults["bus2"], default=fill(0,nphases))
        else
            # if delta connected, ignore bus2 and generate t_terminals such that
            # it corresponds to a delta winding
            t_terminals = [f_terminals[2:end]..., f_terminals[1]]
        end

        eng_obj = Dict{String,Any}()

        eng_obj["bus"] = bus_name

        eng_obj["configuration"] = conn

        eng_obj["f_terminals"] = f_terminals
        eng_obj["t_terminals"] = t_terminals

        eng_obj["kv"] = defaults["kv"]
        eng_obj["kvar"] = defaults["kvar"]

        eng_obj["phases"] = defaults["phases"]

        eng_obj["status"] = convert(Int, defaults["enabled"])

        eng_obj["source_id"] = "capacitor.$name"

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
    for (name, dss_obj) in get(data_dss, "reactor", Dict{String,Any}())
        if !haskey(dss_obj, "bus2")
            _apply_like!(dss_obj, data_dss, "reactor")
            defaults = _apply_ordered_properties(_create_reactor(dss_obj["bus1"], dss_obj["name"]; _to_kwargs(dss_obj)...), dss_obj)

            eng_obj = Dict{String,Any}()

            eng_obj["phases"] = defaults["phases"]

            eng_obj["configuration"] = defaults["conn"]
            connections_default = eng_obj["configuration"] == "wye" ? [collect(1:eng_obj["phases"])..., 0] : collect(1:eng_obj["phases"])
            eng_obj["connections"] = _get_conductors_ordered_dm(defaults["bus1"], default=connections_default, check_length=false)

            eng_obj["bus"] = _parse_busname(defaults["bus1"])[1]
            eng_obj["kvar"] = defaults["kvar"]
            eng_obj["status"] = convert(Int, defaults["enabled"])
            eng_obj["source_id"] = "reactor.$name"

            if import_all
                _import_all!(eng_obj, defaults, dss_obj["prop_order"])
            end

            if !haskey(data_eng, "shunt_reactor")
                data_eng["shunt_reactor"] = Dict{String,Any}()
            end

            data_eng["shunt_reactor"][name] = eng_obj
        end
    end
end


"Adds generators to `data_eng` from `data_dss`"
function _dss2eng_generator!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "generator", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "generator")
        defaults = _apply_ordered_properties(_create_generator(dss_obj["bus1"], dss_obj["name"]; _to_kwargs(dss_obj)...), dss_obj)

        eng_obj = Dict{String,Any}()

        eng_obj["phases"] = defaults["phases"]
        eng_obj["connections"] = _get_conductors_ordered_dm(defaults["bus1"], check_length=false)

        eng_obj["bus"] = _parse_busname(defaults["bus1"])[1]

        eng_obj["kw"] = defaults["kw"]
        eng_obj["kvar"] = defaults["kvar"]
        eng_obj["kv"] = defaults["kv"]

        eng_obj["kvar_min"] = defaults["minkvar"]
        eng_obj["kvar_max"] = defaults["maxkvar"]

        eng_obj["control_model"] = defaults["model"]

        eng_obj["model"] = 2
        eng_obj["startup"] = 0.0
        eng_obj["shutdown"] = 0.0
        eng_obj["ncost"] = 3
        eng_obj["cost"] = [0.0, 1.0, 0.0]

        eng_obj["configuration"] = "wye"

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


"Adds vsources to `data_eng` from `data_dss`"
function _dss2eng_vsource!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "vsource", Dict{<:Any,<:Any}())
        eng_obj = Dict{String,Any}()

        if !haskey(data_eng, "vsource")
            data_eng["vsource"] = Dict{String,Dict{String,Any}}()
        end

        data_eng["vsource"][name] = eng_obj
    end
end


"Adds lines to `data_eng` from `data_dss`"
function _dss2eng_linecode!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "linecode", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "linecode")

        dss_obj["units"] = get(dss_obj, "units", "none")
        dss_obj["circuit_basefreq"] = data_eng["settings"]["basefreq"]

        defaults = _apply_ordered_properties(_create_linecode(name; _to_kwargs(dss_obj)...), dss_obj)

        eng_obj = Dict{String,Any}()

        nphases = defaults["nphases"]

        eng_obj["rs"] = reshape(defaults["rmatrix"], nphases, nphases)
        eng_obj["xs"] = reshape(defaults["xmatrix"], nphases, nphases)

        eng_obj["b_fr"] = reshape(defaults["cmatrix"], nphases, nphases) ./ 2.0
        eng_obj["b_to"] = reshape(defaults["cmatrix"], nphases, nphases) ./ 2.0
        eng_obj["g_fr"] = fill(0.0, nphases, nphases)
        eng_obj["g_to"] = fill(0.0, nphases, nphases)

        if !haskey(data_eng, "linecode")
            data_eng["linecode"] = Dict{String,Any}()
        end

        data_eng["linecode"][name] = eng_obj
    end
end


"Adds lines to `data_eng` from `data_dss`"
function _dss2eng_line!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "line", Dict())
        _apply_like!(dss_obj, data_dss, "line")

        if haskey(dss_obj, "basefreq") && dss_obj["basefreq"] != data_eng["settings"]["basefreq"]
            Memento.warn(_LOGGER, "basefreq=$(dss_obj["basefreq"]) on line $(dss_obj["name"]) does not match circuit basefreq=$(data_eng["settings"]["basefreq"])")
            dss_obj["circuit_basefreq"] = data_eng["settings"]["basefreq"]
        end

        defaults = _apply_ordered_properties(_create_line(dss_obj["bus1"], dss_obj["bus2"], dss_obj["name"]; _to_kwargs(dss_obj)...), dss_obj)

        eng_obj = Dict{String,Any}()

        eng_obj["f_bus"] = _parse_busname(defaults["bus1"])[1]
        eng_obj["t_bus"] = _parse_busname(defaults["bus2"])[1]

        eng_obj["length"] = defaults["length"]

        nphases = defaults["phases"]
        eng_obj["n_conductors"] = nphases

        if haskey(dss_obj, "linecode")
            eng_obj["linecode"] = dss_obj["linecode"]
        end

        if any(haskey(dss_obj, key) for key in ["r0", "r1", "rg", "rmatrix"]) || !haskey(dss_obj, "linecode")
            eng_obj["rs"] = reshape(defaults["rmatrix"], nphases, nphases)
        end

        if any(haskey(dss_obj, key) for key in ["x0", "x1", "xg", "xmatrix"]) || !haskey(dss_obj, "linecode")
            eng_obj["xs"] = reshape(defaults["xmatrix"], nphases, nphases)
        end

        if any(haskey(dss_obj, key) for key in ["b0", "b1", "c0", "c1", "cmatrix"]) || !haskey(dss_obj, "linecode")
            eng_obj["b_fr"] = reshape(defaults["cmatrix"], nphases, nphases) ./ 2.0
            eng_obj["b_to"] = reshape(defaults["cmatrix"], nphases, nphases) ./ 2.0
            eng_obj["g_fr"] = fill(0.0, nphases, nphases)
            eng_obj["g_to"] = fill(0.0, nphases, nphases)
        end

        eng_obj["f_connections"] = _get_conductors_ordered_dm(defaults["bus1"], default=collect(1:nphases))
        eng_obj["t_connections"] = _get_conductors_ordered_dm(defaults["bus2"], default=collect(1:nphases))

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


"Adds transformers to `data_eng` from `data_dss`"
function _dss2eng_transformer!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
   if !haskey(data_eng, "transformer")
        data_eng["transformer"] = Dict{String,Any}()
    end

    for (name, dss_obj) in get(data_dss, "transformer", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "transformer")
        defaults = _apply_ordered_properties(_create_transformer(dss_obj["name"]; _to_kwargs(dss_obj)...), dss_obj)

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

        eng_obj["nphases"] = defaults["phases"]

        for w in 1:nrw
            eng_obj["bus"][w] = _parse_busname(defaults["buses"][w])[1]

            conf = dyz_map[defaults["conns"][w]]
            eng_obj["configuration"][w] = conf

            terminals_default = conf=="wye" ? [1:nphases..., 0] : collect(1:nphases)

            # append ground if connections one too short
            terminals_w = _get_conductors_ordered_dm(defaults["buses"][w], default=terminals_default, pad_ground=(conf=="wye"))
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

            if w>1
                prim_conf = eng_obj["configuration"][1]
                if defaults["leadlag"] in ["ansi", "lag"]
                    if prim_conf=="delta" && conf=="wye"
                        eng_obj["polarity"][w] = -1
                        eng_obj["connections"][w] = [_barrel_roll(eng_obj["connections"][w][1:end-1], 1)..., eng_obj["connections"][w][end]]
                    end
                else # hence defaults["leadlag"] in ["euro", "lead"]
                    if prim_conf=="wye" && conf=="delta"
                        eng_obj["polarity"][w] = -1
                        eng_obj["connections"][w] = _barrel_roll(eng_obj["connections"][w], -1)
                    end

                end
            end
        end

        eng_obj["source_id"] = "transformer.$(name)"

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

        eng_obj = create_transformer_nw(; Dict(Symbol.(keys(eng_obj)).=>values(eng_obj))...)

        if !haskey(data_eng, "transformer")
            data_eng["transformer"] = Dict{String,Any}()
        end

        if import_all
            _import_all!(data_eng, defaults, dss_obj["prop_order"])
        end

        data_eng["transformer"][name] = eng_obj
        # add_virtual!(data_eng, "transformer", eng_obj)
    end
end


"Adds line reactors to `data_eng` from `data_dss`"
function _dss2eng_line_reactor!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "reactor", Dict{String,Any}())
        if haskey(dss_obj, "bus2")
            Memento.warn(_LOGGER, "reactors as constant impedance elements is not yet supported, treating reactor.$name like line")
            _apply_like!(dss_obj, data_dss, "reactor")
            defaults = _apply_ordered_properties(_create_reactor(dss_obj["bus1"], name, dss_obj["bus2"]; _to_kwargs(dss_obj)...), dss_obj)

            eng_obj = Dict{String,Any}()

            nphases = defaults["phases"]
            eng_obj["phases"] = nphases

            eng_obj["f_bus"] = _parse_busname(defaults["bus1"])[1]
            eng_obj["t_bus"] = _parse_busname(defaults["bus2"])[1]

            eng_obj["br_r"] = diagm(0 => fill(0.2, nphases))
            eng_obj["br_x"] = zeros(nphases, nphases)

            eng_obj["g_fr"] = fill(0.0, nphases)
            eng_obj["g_to"] = fill(0.0, nphases)
            eng_obj["b_fr"] = fill(0.0, nphases)
            eng_obj["b_to"] = fill(0.0, nphases)

            for key in ["g_fr", "g_to", "b_fr", "b_to"]
                eng_obj[key] = diagm(0=>eng_obj[key])
            end

            eng_obj["c_rating_a"] = defaults["normamps"]
            eng_obj["c_rating_b"] = defaults["emergamps"]
            eng_obj["c_rating_c"] = defaults["emergamps"]

            eng_obj["tap"] = fill(1.0, nphases)
            eng_obj["shift"] = fill(0.0, nphases)

            eng_obj["br_status"] = convert(Int, defaults["enabled"])

            eng_obj["angmin"] = fill(-60.0, nphases)
            eng_obj["angmax"] = fill( 60.0, nphases)

            eng_obj["transformer"] = true

            eng_obj["source_id"] = "reactor.$(name)"

            if import_all
                _import_all!(eng_obj, defaults, dss_obj["prop_order"])
            end

            if !haskey(data_eng, "line_reactor")
                data_eng["line_reactor"] = Dict{String,Any}()
            end

            data_eng["line_reactor"][name] = eng_obj
        end
    end
end


"Adds pvsystems to `data_eng` from `data_dss`"
function _dss2eng_pvsystem!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "pvsystem", Dict{String,Any}())
        Memento.warn(_LOGGER, "Converting PVSystem \"$(dss_obj["name"])\" into generator with limits determined by OpenDSS property 'kVA'")

        _apply_like!(dss_obj, data_dss, "pvsystem")
        defaults = _apply_ordered_properties(_create_pvsystem(dss_obj["bus1"], dss_obj["name"]; _to_kwargs(dss_obj)...), dss_obj)

        eng_obj = Dict{String,Any}()

        eng_obj["name"] = name
        eng_obj["bus"] = _parse_busname(defaults["bus1"])[1]

        eng_obj["phases"] = defaults["phases"]
        eng_obj["connections"] = _get_conductors_ordered_dm(defaults["bus1"], check_length=false)

        eng_obj["kva"] = defaults["kva"]
        eng_obj["kv"] = defaults["kv"]

        eng_obj["model"] = 2
        eng_obj["startup"] = 0.0
        eng_obj["shutdown"] = 0.0
        eng_obj["ncost"] = 3
        eng_obj["cost"] = [0.0, 1.0, 0.0]

        eng_obj["configuration"] = "wye"

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


"Adds storage to `data_eng` from `data_dss`"
function _dss2eng_storage!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "storage", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "storage")
        defaults = _apply_ordered_properties(_create_storage(dss_obj["bus1"], dss_obj["name"]; _to_kwargs(dss_obj)...), dss_obj)

        eng_obj = Dict{String,Any}()

        eng_obj["phases"] = defaults["phases"]
        eng_obj["connections"] = _get_conductors_ordered_dm(defaults["bus1"], check_length=false)

        eng_obj["name"] = name
        eng_obj["bus"] = _parse_busname(defaults["bus1"])[1]

        eng_obj["kwhstored"] = defaults["kwhstored"]
        eng_obj["kwhrated"] = defaults["kwhrated"]
        eng_obj["kwrated"] = defaults["kwrated"]

        eng_obj["%charge"] = defaults["%charge"]
        eng_obj["%discharge"] = defaults["%discharge"]
        eng_obj["%effcharge"] = defaults["%effcharge"]
        eng_obj["%effdischarge"] = defaults["%effdischarge"]
        eng_obj["kva"] = defaults["kva"]
        eng_obj["kvar"] = defaults["kvar"]
        eng_obj["%r"] = defaults["%r"]
        eng_obj["%x"] = defaults["%x"]
        eng_obj["%idlingkw"] = defaults["%idlingkw"]
        eng_obj["%idlingkvar"] = defaults["%idlingkvar"]

        eng_obj["status"] = convert(Int, defaults["enabled"])


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


"Adds sourcebus as a voltage source to `data_eng` from `data_dss`"
function _dss2eng_sourcebus!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool=false)
    circuit = _create_vsource(get(data_dss["circuit"], "bus1", "sourcebus"), data_dss["circuit"]["name"]; _to_kwargs(data_dss["circuit"])...)

    nodes = Array{Bool}([1 1 1 0])
    ph1_ang = circuit["angle"]
    vm_pu = circuit["pu"]

    phases = circuit["phases"]
    vnom = data_eng["settings"]["set_vbase_val"]

    vm = fill(vm_pu, phases)*vnom
    va = rad2deg.(_wrap_to_pi.([-2*pi/phases*(i-1)+deg2rad(ph1_ang) for i in 1:phases]))

    eng_obj = Dict{String,Any}(
        "bus" => circuit["bus1"],
        "connections" => collect(1:phases),
        "vm" => vm,
        "va" => va,
        "rs" => circuit["rmatrix"],
        "xs" => circuit["xmatrix"],
        "status" => 1
    )

    if import_all
        _import_all!(eng_obj, circuit, data_dss["circuit"]["prop_order"])
    end

    if !haskey(data_eng, "voltage_source")
        data_eng["voltage_source"] = Dict{String,Any}()
    end

    data_eng["voltage_source"][circuit["name"]] = eng_obj
end


"Parses a DSS file into a PowerModels usable format"
function parse_opendss_dm(io::IOStream; import_all::Bool=false, bank_transformers::Bool=true)::Dict
    data_dss = parse_dss(io)

    return parse_opendss_dm(data_dss; import_all=import_all)
end


"Parses a Dict resulting from the parsing of a DSS file into a PowerModels usable format"
function parse_opendss_dm(data_dss::Dict{String,<:Any}; import_all::Bool=false, bank_transformers::Bool=true)::Dict{String,Any}
    data_eng = Dict{String,Any}(
        "source_type" => data_dss["source_type"],
        "settings" => Dict{String,Any}(),
    )

    if import_all
        data_eng["dss_options"] = data_dss["options"]
    end

    if haskey(data_dss, "circuit")
        circuit = data_dss["circuit"]
        defaults = _create_vsource(get(circuit, "bus1", "sourcebus"), circuit["name"]; _to_kwargs(circuit)...)

        data_eng["name"] = circuit["name"]
        data_eng["sourcebus"] = defaults["bus1"]
        data_eng["data_model"] = "engineering"

        # TODO rename fields
        # TODO fix scale factors
        data_eng["settings"]["v_var_scalar"] = 1e3
        # data_eng["settings"]["set_vbase_val"] = defaults["basekv"] / sqrt(3) * 1e3
        data_eng["settings"]["set_vbase_val"] = defaults["basekv"]
        data_eng["settings"]["set_vbase_bus"] = data_eng["sourcebus"]
        data_eng["settings"]["set_sbase_val"] = defaults["basemva"]
        data_eng["settings"]["basefreq"] = get(data_dss["options"], "defaultbasefreq", 60.0)

        data_eng["files"] = data_dss["filename"]
    else
        Memento.error(_LOGGER, "Circuit not defined, not a valid circuit!")
    end

    _dss2eng_bus!(data_eng, data_dss, import_all)
    _dss2eng_buscoords!(data_eng, data_dss)

    _dss2eng_linecode!(data_eng, data_dss, import_all)
    _dss2eng_line!(data_eng, data_dss, import_all)

    # _dss2eng_xfrmcode!(data_eng, data_dss, import_all)  # TODO
    _dss2eng_transformer!(data_eng, data_dss, import_all)

    _dss2eng_line_reactor!(data_eng, data_dss, import_all)

    _dss2eng_loadshape!(data_eng, data_dss, import_all)
    _dss2eng_load!(data_eng, data_dss, import_all)

    _dss2eng_capacitor!(data_eng, data_dss, import_all)
    _dss2eng_shunt_reactor!(data_eng, data_dss, import_all)

    _dss2eng_sourcebus!(data_eng, data_dss, import_all)
    _dss2eng_generator!(data_eng, data_dss, import_all)
    _dss2eng_pvsystem!(data_eng, data_dss, import_all)
    _dss2eng_storage!(data_eng, data_dss, import_all)

    _dss2eng_vsource!(data_eng, data_dss, import_all)

    if bank_transformers
        _bank_transformers!(data_eng)
    end

    _discover_terminals!(data_eng)

    return data_eng
end
