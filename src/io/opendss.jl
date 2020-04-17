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

        eng_obj = Dict{String,Any}(
            "hour" => defaults["hour"],
            "pmult" => defaults["pmult"],
            "qmult" => defaults["qmult"],
            "use_actual" => defaults["useactual"],
            "source_id" => "loadshape.$name",
        )

        if import_all
            _import_all!(eng_obj, dss_obj, dss_obj["prop_order"])
        end

        _add_eng_obj!(data_eng, "loadshape", name, eng_obj)
    end
end


"Adds xycurve objects to `data_eng` from `data_dss`"
function _dss2eng_xycurve!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool=false)
    for (name, dss_obj) in get(data_dss, "xycurve", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "xycurve")
        defaults = _apply_ordered_properties(_create_xycurve(name; _to_kwargs(dss_obj)...), dss_obj)

        xarray = defaults["xarray"] .* defaults["xscale"] .+ defaults["xshift"]
        yarray = defaults["yarray"] .* defaults["yscale"] .+ defaults["yshift"]

        @assert length(xarray) >= 2 && length(yarray) >= 2 "XYCurve data must have two or more points"

        k = length(xarray) - 1

        eng_obj = Dict{String,Any}(
            "interpolated_curve" => Spline1D(xarray, yarray; k=min(k, 5)),
            "source_id" => "xycurve.$name",
        )

        if import_all
            _import_all!(eng_obj, dss_obj, dss_obj["prop_order"])
        end

        _add_eng_obj!(data_eng, "xycurve", name, eng_obj)
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
        defaults = _apply_ordered_properties(_create_load(name; _to_kwargs(dss_obj)...), dss_obj)

        # parse the model
        model = defaults["model"]

        if model == 3
            Memento.warn(_LOGGER, "$name: load model 3 not supported. Treating as model 1.")
            model = 1
        elseif model == 4
            Memento.warn(_LOGGER, "$name: load model 4 not supported. Treating as model 1.")
            model = 1
        elseif model == 6
            Memento.warn(_LOGGER, "$name: load model 6 identical to model 1 in current feature set. Treating as model 1.")
            model = 1
        elseif model == 7
            Memento.warn(_LOGGER, "$name: load model 7 not supported. Treating as model 1.")
            model = 1
        elseif model == 8
            Memento.warn(_LOGGER, "$name: load model 8 not supported. Treating as model 1.")
            model = 1
        end

        model_int2str = Dict(1=>"constant_power", 2=>"constant_impedance", 5=>"constant_current")
        model = model_int2str[model]

        nphases = defaults["phases"]
        conf = defaults["conn"]

        if conf=="delta"
            @assert(nphases in [1, 3], "$name: only 1 and 3-phase delta loads are supported!")
        end

        # connections
        bus = _parse_busname(defaults["bus1"])[1]
        connections_default = conf=="wye" ? [collect(1:nphases)..., 0] : nphases==1 ? [1,2] : [1,2,3]
        connections = _get_conductors_ordered(defaults["bus1"], default=connections_default, pad_ground=(conf=="wye"))

        @assert(length(unique(connections))==length(connections), "$name: connections cannot be made to a terminal more than once.")

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
            _register_awaiting_ground!(data_eng["bus"][bus], connections)
        end

        kv = defaults["kv"]
        if conf=="wye" && nphases in [2, 3]
            kv = kv/sqrt(3)
        end

        eng_obj["vnom"] = kv

        eng_obj["pd_nom"] = fill(defaults["kw"]/nphases, nphases)
        eng_obj["qd_nom"] = fill(defaults["kvar"]/nphases, nphases)

        if import_all
            _import_all!(eng_obj, dss_obj, dss_obj["prop_order"])
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
        defaults = _apply_ordered_properties(_create_capacitor(name; _to_kwargs(dss_obj)...), dss_obj)

        nphases = defaults["phases"]
        conn = defaults["conn"]

        f_terminals = _get_conductors_ordered(defaults["bus1"], default=collect(1:nphases))
        if conn=="wye"
            t_terminals = _get_conductors_ordered(defaults["bus2"], default=fill(0,nphases))
        else
            # if delta connected, ignore bus2 and generate t_terminals such that
            # it corresponds to a delta winding
            t_terminals = [f_terminals[2:end]..., f_terminals[1]]
        end

        eng_obj = Dict{String,Any}(
            "configuration" => conn,
            "status" => convert(Int, defaults["enabled"]),
            "source_id" => "capacitor.$name",
            "phases" => nphases,
            "f_connections" => f_terminals,
            "t_connections" => t_terminals,
        )

        if import_all
            _import_all!(eng_obj, dss_obj, dss_obj["prop_order"])
        end

        bus_name = _parse_busname(defaults["bus1"])[1]
        bus2_name = _parse_busname(defaults["bus2"])[1]
        if bus_name == bus2_name
            eng_obj["bus"] = bus_name

            # 'kv' is specified as phase-to-phase for phases=2/3 (unsure for 4 and more)
            #TODO figure out for more than 3 phases
            vnom_ln = defaults["kv"]
            if defaults["phases"] in [2,3]
                vnom_ln = vnom_ln/sqrt(3)
            end
            defaults["kv"] = fill(vnom_ln, nphases)

            # 'kvar' is specified for all phases at once; we want the per-phase one
            defaults["kvar"] = fill(defaults["kvar"] / nphases, nphases)

            # TODO check unit conversion on qnom/b
            vnom_ln = defaults["kv"]
            qnom = defaults["kvar"] ./ 1e3
            b = qnom ./ vnom_ln.^2

            # convert to a shunt matrix
            terminals, B = _calc_shunt(f_terminals, t_terminals, b)

            # if one terminal is ground (0), reduce shunt addmittance matrix
            terminals, B = _calc_ground_shunt_admittance_matrix(terminals, B, 0)

            eng_obj["bs"] = B

            _add_eng_obj!(data_eng, "shunt_capacitor", name, eng_obj)
        else
            eng_obj["f_bus"] = bus_name
            eng_obj["t_bus"] = bus2_name

            eng_obj["kv"] = fill(defaults["kv"] / nphases, nphases)
            eng_obj["kvar"] = fill(defaults["kvar"] / nphases, nphases)

            _add_eng_obj!(data_eng, "series_capacitor", name, eng_obj)
        end
    end
end



"Adds shunt reactors to `data_eng` from `data_dss`"
function _dss2eng_reactor!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "reactor", Dict{String,Any}())
        if !haskey(dss_obj, "bus2")
            _apply_like!(dss_obj, data_dss, "reactor")
            defaults = _apply_ordered_properties(_create_reactor(name; _to_kwargs(dss_obj)...), dss_obj)

            nphases = defaults["phases"]

            eng_obj = Dict{String,Any}(
                "configuration" => defaults["conn"],
                "bus" => _parse_busname(defaults["bus1"])[1],
                "status" => convert(Int, defaults["enabled"]),
                "source_id" => "reactor.$name",
            )

            connections_default = eng_obj["configuration"] == "wye" ? [collect(1:nphases)..., 0] : collect(1:nphases)
            eng_obj["connections"] = _get_conductors_ordered(defaults["bus1"], default=connections_default, check_length=false)

            # if the ground is used directly, register
            if 0 in eng_obj["connections"]
                _register_awaiting_ground!(data_eng["bus"][eng_obj["bus"]], eng_obj["connections"])
            end

            # TODO Check unit conversion on Gcap
            Gcap = sum(defaults["kvar"]) / (nphases * 1e3 * (data_eng["settings"]["vbase"])^2)

            eng_obj["bs"] = diagm(0=>fill(Gcap, nphases))

            if import_all
                _import_all!(eng_obj, dss_obj, dss_obj["prop_order"])
            end

            _add_eng_obj!(data_eng, "shunt_reactor", name, eng_obj)
        else
            Memento.warn(_LOGGER, "reactors as constant impedance elements is not yet supported, treating reactor.$name like line")
            _apply_like!(dss_obj, data_dss, "reactor")
            defaults = _apply_ordered_properties(_create_reactor(name; _to_kwargs(dss_obj)...), dss_obj)

            eng_obj = Dict{String,Any}()

            nphases = defaults["phases"]

            eng_obj["f_bus"] = _parse_busname(defaults["bus1"])[1]
            eng_obj["t_bus"] = _parse_busname(defaults["bus2"])[1]

            eng_obj["length"] = 1.0

            eng_obj["rs"] = diagm(0 => fill(0.2, nphases))

            for key in ["xs", "g_fr", "g_to", "b_fr", "b_to"]
                eng_obj[key] = zeros(nphases, nphases)
            end

            eng_obj["f_connections"] = _get_conductors_ordered(defaults["bus1"], default=collect(1:nphases))
            eng_obj["t_connections"] = _get_conductors_ordered(defaults["bus2"], default=collect(1:nphases))

            # if the ground is used directly, register
            if 0 in eng_obj["f_connections"]
                _register_awaiting_ground!(data_eng["bus"][eng_obj["f_bus"]], eng_obj["f_connections"])
            end
            if 0 in eng_obj["t_connections"]
                _register_awaiting_ground!(data_eng["bus"][eng_obj["t_bus"]], eng_obj["t_connections"])
            end

            eng_obj["status"] = convert(Int, defaults["enabled"])

            eng_obj["source_id"] = "reactor.$(name)"

            if import_all
                _import_all!(eng_obj, dss_obj, dss_obj["prop_order"])
            end

            _add_eng_obj!(data_eng, "line_reactor", name, eng_obj)
        end
    end
end


"Adds generators to `data_eng` from `data_dss`"
function _dss2eng_generator!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "generator", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "generator")
        defaults = _apply_ordered_properties(_create_generator(name; _to_kwargs(dss_obj)...), dss_obj)

        nphases = defaults["phases"]

        eng_obj = Dict{String,Any}(
            "phases" => nphases,
            "connections" => _get_conductors_ordered(defaults["bus1"], pad_ground=true, default=collect(1:defaults["phases"]+1)),
            "bus" => _parse_busname(defaults["bus1"])[1],
            "pg" => fill(defaults["kw"] / nphases, nphases),
            "qg" => fill(defaults["kvar"] / nphases, nphases),
            "vg" => fill(defaults["kv"], nphases),
            "qg_lb" => fill(defaults["minkvar"] / nphases, nphases),
            "qg_ub" => fill(defaults["maxkvar"] / nphases, nphases),
            "control_mode" => defaults["model"],
            "configuration" => "wye",
            "status" => convert(Int, defaults["enabled"]),
            "source_id" => "generator.$(name)"
        )

        # if the ground is used directly, register load
        if 0 in eng_obj["connections"]
            _register_awaiting_ground!(data_eng["bus"][eng_obj["bus"]], eng_obj["connections"])
        end

        if import_all
            _import_all!(eng_obj, dss_obj, dss_obj["prop_order"])
        end

        _add_eng_obj!(data_eng, "generator", name, eng_obj)
    end
end


"Adds vsources to `data_eng` from `data_dss`"
function _dss2eng_vsource!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "vsource", Dict{<:Any,<:Any}())
        _apply_like!(dss_obj, data_dss, "vsource")
        defaults = _create_vsource(name; _to_kwargs(dss_obj)...)


        ph1_ang = defaults["angle"]
        vm_pu = defaults["pu"]

        phases = defaults["phases"]
        vnom = data_eng["settings"]["vbase"]

        vm = fill(vm_pu, phases)*vnom
        va = rad2deg.(_wrap_to_pi.([-2*pi/phases*(i-1)+deg2rad(ph1_ang) for i in 1:phases]))

        eng_obj = Dict{String,Any}(
            "bus" => _parse_busname(defaults["bus1"])[1],
            "connections" => collect(1:phases),
            "vm" => vm,
            "va" => va,
            "rs" => defaults["rmatrix"],
            "xs" => defaults["xmatrix"],
            "source_id" => "vsource.$name",
            "status" => convert(Int, defaults["enabled"])
        )

        # if the ground is used directly, register load
        if 0 in eng_obj["connections"]
            _register_awaiting_ground!(data_eng["bus"][eng_obj["bus"]], eng_obj["connections"])
        end

        if import_all
            _import_all!(eng_obj, dss_obj, dss_obj["prop_order"])
        end

        _add_eng_obj!(data_eng, "voltage_source", name, eng_obj)
    end
end


"Adds lines to `data_eng` from `data_dss`"
function _dss2eng_linecode!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "linecode", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "linecode")

        dss_obj["circuit_basefreq"] = data_eng["settings"]["base_frequency"]

        defaults = _apply_ordered_properties(_create_linecode(name; _to_kwargs(dss_obj)...), dss_obj)

        eng_obj = Dict{String,Any}()

        nphases = defaults["nphases"]

        eng_obj["rs"] = reshape(defaults["rmatrix"], nphases, nphases)
        eng_obj["xs"] = reshape(defaults["xmatrix"], nphases, nphases)

        eng_obj["b_fr"] = reshape(defaults["cmatrix"], nphases, nphases) ./ 2.0
        eng_obj["b_to"] = reshape(defaults["cmatrix"], nphases, nphases) ./ 2.0
        eng_obj["g_fr"] = fill(0.0, nphases, nphases)
        eng_obj["g_to"] = fill(0.0, nphases, nphases)

        if import_all
            _import_all!(eng_obj, dss_obj, dss_obj["prop_order"])
        end

        _add_eng_obj!(data_eng, "linecode", name, eng_obj)
    end
end


"Adds lines to `data_eng` from `data_dss`"
function _dss2eng_line!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "line", Dict())
        _apply_like!(dss_obj, data_dss, "line")

        if haskey(dss_obj, "basefreq") && dss_obj["basefreq"] != data_eng["settings"]["base_frequency"]
            Memento.warn(_LOGGER, "basefreq=$(dss_obj["basefreq"]) on line $(dss_obj["name"]) does not match circuit basefreq=$(data_eng["settings"]["base_frequency"])")
        end

        defaults = _apply_ordered_properties(_create_line(name; _to_kwargs(dss_obj)...), dss_obj)

        _like_is_switch = haskey(dss_obj, "like") && get(get(data_dss["line"], dss_obj["like"], Dict{String,Any}()), "switch", false)
        nphases = defaults["phases"]

        eng_obj = Dict{String,Any}(
            "f_bus" => _parse_busname(defaults["bus1"])[1],
            "t_bus" => _parse_busname(defaults["bus2"])[1],
            "length" => defaults["switch"] || _like_is_switch ? 0.001 : defaults["length"],
            "f_connections" => _get_conductors_ordered(defaults["bus1"], default=collect(1:nphases)),
            "t_connections" => _get_conductors_ordered(defaults["bus2"], default=collect(1:nphases)),
            "status" => convert(Int, defaults["enabled"]),
            "source_id" => "line.$name"
        )

        if haskey(dss_obj, "linecode")
            eng_obj["linecode"] = dss_obj["linecode"]
        end

        if any(haskey(dss_obj, key) && _is_after_linecode(dss_obj["prop_order"], key) for key in ["r0", "r1", "rg", "rmatrix"]) || !haskey(dss_obj, "linecode")
            eng_obj["rs"] = reshape(defaults["rmatrix"], nphases, nphases)
        end

        if any(haskey(dss_obj, key) && _is_after_linecode(dss_obj["prop_order"], key) for key in ["x0", "x1", "xg", "xmatrix"]) || !haskey(dss_obj, "linecode")
            eng_obj["xs"] = reshape(defaults["xmatrix"], nphases, nphases)
        end

        if any(haskey(dss_obj, key) && _is_after_linecode(dss_obj["prop_order"], key) for key in ["b0", "b1", "c0", "c1", "cmatrix"]) || !haskey(dss_obj, "linecode")
            eng_obj["b_fr"] = reshape(defaults["cmatrix"], nphases, nphases) ./ 2.0
            eng_obj["b_to"] = reshape(defaults["cmatrix"], nphases, nphases) ./ 2.0
            eng_obj["g_fr"] = fill(0.0, nphases, nphases)
            eng_obj["g_to"] = fill(0.0, nphases, nphases)
        end

        # if the ground is used directly, register
        for prefix in ["f_", "t_"]
            if 0 in eng_obj["$(prefix)connections"]
                _register_awaiting_ground!(data_eng["bus"][eng_obj["$(prefix)bus"]], eng_obj["$(prefix)connections"])
            end
        end

        if import_all
            _import_all!(eng_obj, dss_obj, dss_obj["prop_order"])
        end

        if defaults["switch"]
            _add_eng_obj!(data_eng, "switch", name, eng_obj)
        else
            _add_eng_obj!(data_eng, "line", name, eng_obj)
        end
    end
end


"Adds transformers to `data_eng` from `data_dss`"
function _dss2eng_xfmrcode!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "xfmrcode", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "xfmrcode")
        defaults = _apply_ordered_properties(_create_xfmrcode(name; _to_kwargs(dss_obj)...), dss_obj)

        nphases = defaults["phases"]
        nrw = defaults["windings"]

        eng_obj = Dict{String,Any}(
            "tm" => Vector{Vector{Float64}}([fill(tap, nphases) for tap in defaults["taps"]]),
            "tm_min" => Vector{Vector{Float64}}(fill(fill(defaults["mintap"], nphases), nrw)),
            "tm_max" => Vector{Vector{Float64}}(fill(fill(defaults["maxtap"], nphases), nrw)),
            "tm_fix" => Vector{Vector{Int}}([fill(1, nphases) for w in 1:nrw]),
            "tm_step" => Vector{Vector{Float64}}(fill(fill(1/32, nphases), nrw)),
            "fixed" => Vector{Int}(fill(1, nrw)),
            "vnom" => Vector{Float64}(defaults["kvs"]),
            "snom" => Vector{Float64}(defaults["kvas"]),
            "configuration" => Vector{String}(defaults["conns"]),
            "rs" => Vector{Float64}(defaults["%rs"] ./ 100),
            "noloadloss" => defaults["%noloadloss"] / 100,
            "imag" => defaults["%imag"] / 100,
            "xsc" => nrw == 2 ? [defaults["xhl"] / 100] : [defaults["xhl"], defaults["xht"], defaults["xlt"]] ./ 100,
            "source_id" => "xfmrcode.$name",
        )

        if import_all
            _import_all!(eng_obj, dss_obj, dss_obj["prop_order"])
        end

        _add_eng_obj!(data_eng, "xfmrcode", name, eng_obj)
    end
end


"Adds transformers to `data_eng` from `data_dss`"
function _dss2eng_transformer!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "transformer", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "transformer")
        defaults = _apply_ordered_properties(_create_transformer(name; _to_kwargs(dss_obj)...), dss_obj)

        eng_obj = Dict{String, Any}(
            "source_id" => "transformer.$name"
        )

        # no way around checking xfmrcode for these properties
        _shared = Dict{String,Any}(
            "leadlag" => defaults["leadlag"],
            "conns" => defaults["conns"],
            "phases" => defaults["phases"],
            "windings" => defaults["windings"]
        )
        if haskey(dss_obj, "xfmrcode")
            xfmrcode_dss_obj = deepcopy(data_dss["xfmrcode"][dss_obj["xfmrcode"]])
            _apply_like!(xfmrcode_dss_obj, data_dss, "xfmrcode")
            xfmrcode = _apply_ordered_properties(_create_xfmrcode(string(dss_obj["xfmrcode"]); _to_kwargs(xfmrcode_dss_obj)...), xfmrcode_dss_obj)

            for key in ["leadlag", "conns", "phases", "windings"]
                if haskey(dss_obj, key) && _is_after_xfmrcode(dss_obj["prop_order"], key)
                    _shared[key] = defaults[key]
                else
                    _shared[key] = xfmrcode[key]
                end
            end
        end

        leadlag = _shared["leadlag"]
        confs = _shared["conns"]
        nphases = _shared["phases"]
        nrw = _shared["windings"]

        # taps
        if isempty(defaults["xfmrcode"]) || (haskey(dss_obj, "taps") && _is_after_xfmrcode(dss_obj["prop_order"], "taps")) || all(haskey(dss_obj, k) && _is_after_xfmrcode(dss_obj["prop_order"], k) for k in ["tap", "tap_2", "tap_3"])
            eng_obj["tm"] = [fill(defaults["taps"][w], nphases) for w in 1:nrw]
        else
            for (w, key_suffix) in enumerate(["", "_2", "_3"])
                if haskey(dss_obj, "tap$(key_suffix)") && _is_after_xfmrcode(dss_obj["prop_order"], "tap$(key_suffix)")
                    if !haskey(eng_obj, "tm")
                        eng_obj["tm"] = Vector{Any}(missing, nrw)
                    end
                    eng_obj["tm"][w] = fill(defaults["taps"][defaults["wdg$(key_suffix)"]], nphases)
                end
            end
        end

        # kvs, kvas
        for (fr_key, to_key) in zip(["kv", "kva"], ["vnom", "snom"])
            if isempty(defaults["xfmrcode"]) || (haskey(dss_obj, "$(fr_key)s") && _is_after_xfmrcode(dss_obj["prop_order"], "$(fr_key)s")) || all(haskey(dss_obj, "$(fr_key)$(key_suffix)") && _is_after_xfmrcode(dss_obj["prop_order"], "$(fr_key)$(key_suffix)") for key_suffix in ["", "_2", "_3"])
                eng_obj[to_key] = defaults["$(fr_key)s"]
            else
                for (w, key_suffix) in enumerate(["", "_2", "_3"])
                    if haskey(dss_obj, "$(fr_key)$(key_suffix)") && _is_after_xfmrcode(dss_obj["prop_order"], "$(fr_key)$(key_suffix)")
                        if !haskey(eng_obj, to_key)
                            eng_obj[to_key] = Vector{Any}(missing, nrw)
                        end
                        eng_obj[to_key][w] = defaults["$(fr_key)s"][defaults["wdg$(key_suffix)"]]
                    end
                end
            end
        end

        # mintap, maxtap
        for (fr_key, to_key) in zip(["mintap", "maxtap"], ["tm_min", "tm_max"])
            if isempty(defaults["xfmrcode"]) || (haskey(dss_obj, fr_key) && _is_after_xfmrcode(dss_obj["prop_order"], fr_key))
                eng_obj[to_key] = fill(fill(defaults[fr_key], nphases), nrw)
            end
        end

        # %noloadloss, %imag
        for (fr_key, to_key) in zip(["%noloadloss", "%imag"], ["noloadloss", "imag"])
            if isempty(defaults["xfmrcode"]) || (haskey(dss_obj, fr_key) && _is_after_xfmrcode(dss_obj["prop_order"], fr_key))
                eng_obj[to_key] = defaults[fr_key] / 100
            end
        end

        # %rs
        if isempty(defaults["xfmrcode"]) || (haskey(dss_obj, "%rs") && _is_after_xfmrcode(dss_obj["prop_order"], "%rs")) || all(haskey(dss_obj, k) && _is_after_xfmrcode(dss_obj["prop_order"], k) for k in ["%r", "%r_2", "%r_3"])
            eng_obj["rs"] = defaults["%rs"] / 100
        else
            for (w, key_suffix) in enumerate(["", "_2", "_3"])
                if haskey(dss_obj, "%r$(key_suffix)") && _is_after_xfmrcode(dss_obj["prop_order"], "%r$(key_suffix)")
                    if !haskey(eng_obj, "rs")
                        eng_obj["rs"] = Vector{Any}(missing, nrw)
                    end
                    eng_obj["rs"][w] = defaults["%rs"][defaults["wdg$(key_suffix)"]] / 100
                end
            end
        end

        # loss model (converted to SI units, referred to secondary)
        if nrw == 2
            if isempty(defaults["xfmrcode"]) || (haskey(dss_obj, "xhl") && _is_after_xfmrcode(dss_obj["prop_order"], "xhl"))
                eng_obj["xsc"] = [defaults["xhl"]] / 100
            end
        elseif nrw == 3
            for (w, key) in enumerate(["xhl", "xht", "xlt"])
                if isempty(defaults["xfmrcode"]) || (haskey(dss_obj, key) && _is_after_xfmrcode(dss_obj["prop_order"], key))
                    if !haskey(eng_obj, "xsc")
                        eng_obj["xsc"] = Vector{Any}(missing, 3)
                    end
                    eng_obj["xsc"][w] = defaults[key] / 100
                end
            end
        end

        # tm_fix, tm_step don't appear in opendss
        if isempty(defaults["xfmrcode"])
            eng_obj["tm_fix"] = [fill(1, nphases) for w in 1:nrw]
            eng_obj["tm_step"] = fill(fill(1/32, nphases), nrw)
            eng_obj["fixed"] = fill(1, nrw)
        end

        # always required
        eng_obj["bus"] = Array{String, 1}(undef, nrw)
        eng_obj["connections"] = Array{Array{Int, 1}, 1}(undef, nrw)
        eng_obj["polarity"] = fill(1, nrw)


        if !haskey(dss_obj, "xfmrcode")
            eng_obj["configuration"] = confs
        end

        # test if this transformer conforms with limitations
        if nphases<3 && "delta" in confs
            Memento.error(_LOGGER, "Transformers with delta windings should have at least 3 phases to be well-defined: $name.")
        end
        if nrw>3
            # All of the code is compatible with any number of windings,
            # except for the parsing of the loss model (the pair-wise reactance)
            Memento.error(_LOGGER, "For now parsing of xscarray is not supported. At most 3 windings are allowed, not $nrw.")
        end

        for w in 1:nrw
            eng_obj["bus"][w] = _parse_busname(defaults["buses"][w])[1]

            conf = confs[w]
            terminals_default = conf=="wye" ? [1:nphases..., 0] : collect(1:nphases)

            # append ground if connections one too short
            terminals_w = _get_conductors_ordered(defaults["buses"][w], default=terminals_default, pad_ground=(conf=="wye"))
            eng_obj["connections"][w] = terminals_w

            if 0 in terminals_w
                _register_awaiting_ground!(data_eng["bus"][eng_obj["bus"][w]], eng_obj["connections"][w])
            end

            if w>1
                prim_conf = confs[1]
                if leadlag in ["ansi", "lag"]
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

        for key in ["bank", "xfmrcode"]
            if !isempty(defaults[key])
                eng_obj[key] = defaults[key]
            end
        end

        if !haskey(data_eng, "transformer")
            data_eng["transformer"] = Dict{String,Any}()
        end

        if import_all
            _import_all!(eng_obj, dss_obj, dss_obj["prop_order"])
        end

        _add_eng_obj!(data_eng, "transformer", name, eng_obj)
    end
end


"Adds pvsystems to `data_eng` from `data_dss`"
function _dss2eng_pvsystem!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "pvsystem", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "pvsystem")
        defaults = _apply_ordered_properties(_create_pvsystem(name; _to_kwargs(dss_obj)...), dss_obj)

        # TODO pick parameters for solar objects

        nphases = defaults["phases"]

        eng_obj = Dict{String,Any}(
            "bus" => _parse_busname(defaults["bus1"])[1],
            "configuration" => defaults["conn"],
            "connections" => _get_conductors_ordered(defaults["bus1"], pad_ground=true, default=collect(1:defaults["phases"]+1)),
            "kva" => fill(defaults["kva"] / nphases, nphases),
            "kvar" => fill(defaults["kvar"] / nphases, nphases),
            "kv" => fill(defaults["kv"] / nphases, nphases),
            "max_rated_power" => fill(defaults["pmpp"] / nphases, nphases),
            "irradiance" => defaults["irradiance"],
            "temperature" => defaults["temperature"],
            "p-t_curve" => defaults["p-tcurve"],
            "efficiency_curve" => defaults["effcurve"],
            "rs" => diagm(0 => fill(defaults["%r"] / 100., nphases)),
            "xs" => diagm(0 => fill(defaults["%x"] / 100., nphases)),
            "status" => convert(Int, defaults["enabled"]),
            "source_id" => "pvsystem.$name",
        )

        # if the ground is used directly, register load
        if 0 in eng_obj["connections"]
            _register_awaiting_ground!(data_eng["bus"][eng_obj["bus"]], eng_obj["connections"])
        end

        if import_all
            _import_all!(eng_obj, dss_obj, dss_obj["prop_order"])
        end

        _add_eng_obj!(data_eng, "solar", name, eng_obj)
    end
end


"Adds storage to `data_eng` from `data_dss`"
function _dss2eng_storage!(data_eng::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, import_all::Bool)
    for (name, dss_obj) in get(data_dss, "storage", Dict{String,Any}())
        _apply_like!(dss_obj, data_dss, "storage")
        defaults = _apply_ordered_properties(_create_storage(name; _to_kwargs(dss_obj)...), dss_obj)

        nphases = defaults["phases"]

        eng_obj = Dict{String,Any}(
            "bus" => _parse_busname(defaults["bus1"])[1],
            "connections" => _get_conductors_ordered(defaults["bus1"], check_length=false),
            "configuration" => "wye",
            "energy" => defaults["kwhstored"],
            "energy_ub" => defaults["kwrated"],
            "charge_ub" => defaults["%charge"] / 100.0 * defaults["kwrated"],
            "discharge_ub" => defaults["%discharge"] / 100.0 * defaults["kwrated"],
            "cm_ub" => fill(defaults["kva"] / nphases, nphases),
            "charge_efficiency" => defaults["%effcharge"],
            "discharge_efficiency" => defaults["%effdischarge"],
            "qs_lb" => -fill(defaults["kvar"] / nphases, nphases),
            "qs_ub" =>  fill(defaults["kvar"] / nphases, nphases),
            "rs" => fill(defaults["%r"] / nphases / 100.0, nphases),
            "xs" => fill(defaults["%x"] / nphases / 100.0, nphases),
            "pex" => defaults["%idlingkw"] .* defaults["kwrated"],
            "qex" => defaults["%idlingkvar"] .* defaults["kvar"],
            "status" => convert(Int, defaults["enabled"]),
            "source_id" => "storage.$(name)",
        )

        # if the ground is used directly, register load
        if 0 in eng_obj["connections"]
            _register_awaiting_ground!(data_eng["bus"][eng_obj["bus"]], eng_obj["connections"])
        end

        if import_all
            _import_all(eng_obj, defaults, dss_obj["prop_order"])
        end

        _add_eng_obj!(data_eng, "storage", name, eng_obj)
    end
end


"Parses a DSS file into a PowerModels usable format"
function parse_opendss(io::IOStream; import_all::Bool=false, bank_transformers::Bool=true)::Dict
    data_dss = parse_dss(io)

    return parse_opendss(data_dss; import_all=import_all, bank_transformers=bank_transformers)
end


"Parses a Dict resulting from the parsing of a DSS file into a PowerModels usable format"
function parse_opendss(data_dss::Dict{String,<:Any}; import_all::Bool=false, bank_transformers::Bool=true)::Dict{String,Any}
    data_eng = Dict{String,Any}(
        "data_model" => "engineering",
        "settings" => Dict{String,Any}(),
    )

    if import_all
        data_eng["dss_options"] = data_dss["options"]
    end

    if haskey(data_dss, "vsource") && haskey(data_dss["vsource"], "source") && haskey(data_dss, "circuit")
        defaults = _create_vsource("source"; _to_kwargs(data_dss["vsource"]["source"])...)
        source_bus = _parse_busname(defaults["bus1"])[1]

        data_eng["name"] = data_dss["circuit"]
        data_eng["data_model"] = "engineering"

        # TODO rename fields
        data_eng["settings"]["v_var_scalar"] = 1e3
        data_eng["settings"]["vbase"] = defaults["basekv"] / sqrt(3)
        data_eng["settings"]["base_bus"] = source_bus
        data_eng["settings"]["sbase"] = defaults["basemva"] * 1e3
        data_eng["settings"]["base_frequency"] = get(get(data_dss, "options", Dict{String,Any}()), "defaultbasefreq", 60.0)

        data_eng["files"] = data_dss["filename"]
    else
        Memento.error(_LOGGER, "Circuit not defined, not a valid circuit!")
    end

    _dss2eng_bus!(data_eng, data_dss, import_all)
    _dss2eng_buscoords!(data_eng, data_dss)

    _dss2eng_linecode!(data_eng, data_dss, import_all)
    _dss2eng_line!(data_eng, data_dss, import_all)

    _dss2eng_xfmrcode!(data_eng, data_dss, import_all)
    _dss2eng_transformer!(data_eng, data_dss, import_all)

    _dss2eng_capacitor!(data_eng, data_dss, import_all)
    _dss2eng_reactor!(data_eng, data_dss, import_all)

    _dss2eng_loadshape!(data_eng, data_dss, import_all)
    _dss2eng_load!(data_eng, data_dss, import_all)

    _dss2eng_xycurve!(data_eng, data_dss, import_all)

    _dss2eng_vsource!(data_eng, data_dss, import_all)
    _dss2eng_generator!(data_eng, data_dss, import_all)
    _dss2eng_pvsystem!(data_eng, data_dss, import_all)
    _dss2eng_storage!(data_eng, data_dss, import_all)

    _discover_terminals!(data_eng)

    if bank_transformers
        _bank_transformers!(data_eng)
    end

    return data_eng
end
