# OpenDSS parser
"Parses buscoords lon,lat (if present) into their respective buses"
function _dss2eng_buscoords!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel)
    for (id, coords) in data_dss.buscoordinates
        if haskey(data_eng["bus"], id)
            bus = data_eng["bus"][id]
            bus["lon"] = coords["x"]
            bus["lat"] = coords["y"]
        end
    end
end


"Adds nodes as buses to `data_eng` from `data_dss`"
function _dss2eng_bus!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool=false)
    buses = _discover_buses(data_dss)

    for id in buses
        @assert !startswith(id, "_virtual") "Bus $id: identifiers should not start with _virtual"

        eng_obj = create_bus(status=ENABLED)

        _add_eng_obj!(data_eng, "bus", id, eng_obj)
    end
end


"Adds loadshapes to `data_eng` from `data_dss`"
function _dss2eng_loadshape!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool=false)
    for (id, dss_obj) in get(data_dss, "loadshape", Dict{String,Any}())
        eng_obj = Dict{String,Any}(
            "time" => dss_obj["hour"],
            "offset" => 0.0,
            "replace" => dss_obj["useactual"],
            "values" => dss_obj["pmult"],
            "source_id" => "loadshape.$id",
        )

        if import_all
            _import_all!(eng_obj, dss_obj)
        end

        if _is_loadshape_split(dss_obj)
            @info "Loadshape '$id' contains mismatched pmult and qmult, splitting into `time_series` ids '$(id)_p' and '$(id)_q'"
            _add_eng_obj!(data_eng, "time_series", "$(id)_p", eng_obj)

            eng_obj_qmult = deepcopy(eng_obj)
            eng_obj_qmult["values"] = dss_obj["qmult"]

            _add_eng_obj!(data_eng, "time_series", "$(id)_q", eng_obj_qmult)
        else
            _add_eng_obj!(data_eng, "time_series", id, eng_obj)
        end
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
function _dss2eng_load!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool, time_series::String="daily")
    for (id, dss_obj) in get(data_dss, "load", Dict{String,Any}())
        nphases = dss_obj["phases"]
        bus = _parse_bus_id(dss_obj["bus1"])[1]
        conf = nphases==1 && dss_obj["kv"]==0.24 ? DELTA : dss_obj["conn"] # check if load is connected between split-phase terminals of triplex node (nominal line-line voltage=240V), TODO: better generalization

        if conf==DELTA
            @assert(nphases in [1, 3], "$id: only 1 and 3-phase delta loads are supported!")
        end

        # connections
        connections_default = conf==WYE ? [collect(1:nphases)..., 0] : nphases==1 ? [1,2] : [1,2,3]
        connections = _get_conductors_ordered(dss_obj["bus1"], default=connections_default, pad_ground=(conf==WYE))

        @assert(length(unique(connections))==length(connections), "$id: connections cannot be made to a terminal more than once.")

        # now we can create the load; if you do not have the correct model,
        # pd/qd fields will be populated by default (should not happen for constant current/impedance)
        eng_obj = Dict{String,Any}(
            "bus" => bus,
            "model" => dss_obj["model"],
            "configuration" => conf,
            "connections" => connections,
            "dispatchable" => NO,
            "source_id" => "load.$id",
            "status" => dss_obj["enabled"]
        )

        # if the ground is used directly, register load
        if 0 in connections
            _register_awaiting_ground!(data_eng["bus"][bus], connections)
        end

        kv = dss_obj["kv"]
        if conf==WYE && nphases in [2, 3]
            kv = kv/sqrt(3)
        end

        eng_obj["vm_nom"] = kv

        eng_obj["pd_nom"] = fill(dss_obj["kw"]/nphases, nphases)
        eng_obj["qd_nom"] = fill(dss_obj["kvar"]/nphases, nphases)

        # if ZIP load, include weighting factors and cut-off voltage
        if eng_obj["model"]==ZIP
            eng_obj["zipv"] = collect(dss_obj["zipv"])
        end

        _build_time_series_reference!(eng_obj, dss_obj, data_dss, time_series, "pd_nom", "qd_nom")

        if import_all
            _import_all!(eng_obj, dss_obj)
        end

        _add_eng_obj!(data_eng, "load", id, eng_obj)
    end
end


"Adds capacitors to `data_eng` from `data_dss`"
function _dss2eng_capacitor!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool)
    for (id, dss_obj) in get(data_dss, "capacitor", Dict{String,Any}())
        nphases = dss_obj["phases"]
        conn = dss_obj["conn"]

        f_terminals = _get_conductors_ordered(dss_obj["bus1"], default=collect(1:nphases))
        if conn==WYE
            t_terminals = _get_conductors_ordered(dss_obj["bus2"], default=fill(0,nphases))
        else
            # if delta connected, ignore bus2 and generate t_terminals such that
            # it corresponds to a delta winding
            t_terminals = [f_terminals[2:end]..., f_terminals[1]]
        end

        eng_obj = Dict{String,Any}(
            "configuration" => conn,
            "model" => CAPACITOR,
            "dispatchable" => NO,
            "status" => dss_obj["enabled"],
            "source_id" => "capacitor.$id",
        )

        if import_all
            _import_all!(eng_obj, dss_obj)
        end

        bus_name = _parse_bus_id(dss_obj["bus1"])[1]
        bus2_name = _parse_bus_id(dss_obj["bus2"])[1]
        if bus_name == bus2_name
            eng_obj["bus"] = bus_name

            # 'kv' is specified as phase-to-phase for phases=2/3 (unsure for 4 and more)
            #TODO figure out for more than 3 phases
            vnom_ln = dss_obj["kv"]
            if dss_obj["phases"] in [2,3]
                vnom_ln = vnom_ln/sqrt(3)
            end
            kv = fill(vnom_ln, nphases)

            # 'kvar' is specified for all phases at once; we want the per-phase one
            kvar = fill(dss_obj["kvar"] / nphases, nphases)

            vnom_ln = kv
            qnom = kvar ./ 1e3
            b = qnom ./ vnom_ln.^2

            # convert to a shunt matrix
            terminals, B = _calc_shunt(f_terminals, t_terminals, b)

            # if one terminal is ground (0), reduce shunt admittance matrix
            terminals, B = _calc_ground_shunt_admittance_matrix(terminals, B, 0)

            eng_obj["gs"] = zeros(size(B))
            eng_obj["bs"] = B
            eng_obj["connections"] = terminals

            _add_eng_obj!(data_eng, "shunt", id, eng_obj)
        else
            @warn "capacitors as constant impedance elements is not yet supported, treating reactor.$id like line"
            _eng_obj = Dict{String,Any}(
                "f_bus" => _parse_bus_id(dss_obj["bus1"])[1],
                "t_bus" => _parse_bus_id(dss_obj["bus2"])[1],
                "f_connections" => _get_conductors_ordered(dss_obj["bus1"], default=collect(1:nphases)),
                "t_connections" => _get_conductors_ordered(dss_obj["bus2"], default=collect(1:nphases)),
                "length" => 1.0,
                "rs" => LinearAlgebra.diagm(0 => fill(0.2, nphases)),
                "xs" => zeros(nphases, nphases),
                "g_fr" => zeros(nphases, nphases),
                "b_fr" => zeros(nphases, nphases),
                "g_to" => zeros(nphases, nphases),
                "b_to" => zeros(nphases, nphases),
                "status" => dss_obj["enabled"],
                "source_id" => "capacitor.$id",
            )

            merge!(eng_obj, _eng_obj)

            # if the ground is used directly, register
            if 0 in eng_obj["f_connections"]
                _register_awaiting_ground!(data_eng["bus"][eng_obj["f_bus"]], eng_obj["f_connections"])
            end
            if 0 in eng_obj["t_connections"]
                _register_awaiting_ground!(data_eng["bus"][eng_obj["t_bus"]], eng_obj["t_connections"])
            end

            if import_all
                _import_all!(eng_obj, dss_obj)
            end

            _add_eng_obj!(data_eng, "line", id, eng_obj)
        end
    end
end


"Adds shunt reactors to `data_eng` from `data_dss`"
function _dss2eng_reactor!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool)
    for (id, dss_obj) in get(data_dss, "reactor", Dict{String,Any}())
        if isempty(dss_obj["bus2"])
            nphases = dss_obj["phases"]

            eng_obj = Dict{String,Any}(
                "bus" => _parse_bus_id(dss_obj["bus1"])[1],
                "configuration" => dss_obj["conn"],
                "model" => REACTOR,
                "dispatchable" => NO,
                "status" => dss_obj["enabled"],
                "source_id" => "reactor.$id",
            )

            connections_default = eng_obj["configuration"] == WYE ? [collect(1:nphases)..., 0] : collect(1:nphases)
            eng_obj["connections"] = _get_conductors_ordered(dss_obj["bus1"], default=connections_default, check_length=false)

            # if the ground is used directly, register
            if 0 in eng_obj["connections"]
                _register_awaiting_ground!(data_eng["bus"][eng_obj["bus"]], eng_obj["connections"])
            end

            # TODO Check unit conversion on Gcap
            Gcap = sum(dss_obj["kvar"]) / (nphases * 1e3 * (dss_obj["kv"] / sqrt(nphases))^2)

            eng_obj["bs"] = LinearAlgebra.diagm(0=>fill(Gcap, nphases))
            eng_obj["gs"] = zeros(size(eng_obj["bs"]))

            if import_all
                _import_all!(eng_obj, dss_obj)
            end

            _add_eng_obj!(data_eng, "shunt", id, eng_obj)
        else
            @warn "reactors as constant impedance elements is not yet supported, treating reactor.$id like line"

            nphases = dss_obj["phases"]

            eng_obj = Dict{String,Any}(
                "f_bus" => _parse_bus_id(dss_obj["bus1"])[1],
                "t_bus" => _parse_bus_id(dss_obj["bus2"])[1],
                "f_connections" => _get_conductors_ordered(dss_obj["bus1"], default=collect(1:nphases)),
                "t_connections" => _get_conductors_ordered(dss_obj["bus2"], default=collect(1:nphases)),
                "length" => 1.0,
                "rs" => dss_obj["rmatrix"],
                "xs" => dss_obj["xmatrix"],
                "g_fr" => zeros(nphases, nphases),
                "b_fr" => zeros(nphases, nphases),
                "g_to" => zeros(nphases, nphases),
                "b_to" => zeros(nphases, nphases),
                "status" => dss_obj["enabled"],
                "source_id" => "reactor.$id",
            )

            # if the ground is used directly, register
            if 0 in eng_obj["f_connections"]
                _register_awaiting_ground!(data_eng["bus"][eng_obj["f_bus"]], eng_obj["f_connections"])
            end
            if 0 in eng_obj["t_connections"]
                _register_awaiting_ground!(data_eng["bus"][eng_obj["t_bus"]], eng_obj["t_connections"])
            end

            if import_all
                _import_all!(eng_obj, dss_obj)
            end

            _add_eng_obj!(data_eng, "line", id, eng_obj)
        end
    end
end


"Adds generators to `data_eng` from `data_dss`"
function _dss2eng_generator!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool, time_series::String="daily")
    for (id, dss_obj) in get(data_dss, "generator", Dict{String,Any}())
        nphases = dss_obj["phases"]
        conf = nphases==1 && dss_obj["kv"]==0.24 ? DELTA : WYE # check if generator is connected between split-phase terminals of triplex node (nominal line-line voltage=240V), TODO: better generalization

        eng_obj = Dict{String,Any}(
            "phases" => nphases,
            "connections" => _get_conductors_ordered(dss_obj["bus1"], pad_ground=dss_obj["conn"] == WYE, default=dss_obj["conn"] == WYE ? [collect(1:dss_obj["phases"])..., 0] : nphases == 1 ? [1,0] : collect(1:nphases)),
            "bus" => _parse_bus_id(dss_obj["bus1"])[1],
            "pg" => fill(dss_obj["kw"] / nphases, nphases),
            "qg" => fill(dss_obj["kvar"] / nphases, nphases),
            "vg" => fill(dss_obj["kv"] / sqrt(nphases), nphases),
            "qg_lb" => fill(dss_obj["minkvar"] / nphases, nphases),
            "qg_ub" => fill(dss_obj["maxkvar"] / nphases, nphases),
            "pg_lb" => fill(0.0, nphases),
            "pg_ub" => fill(dss_obj["kw"] / nphases, nphases),
            "control_mode" => FREQUENCYDROOP,
            "configuration" => conf,
            "status" => dss_obj["enabled"],
            "source_id" => "generator.$id"
        )

        # if the ground is used directly, register load
        if 0 in eng_obj["connections"]
            _register_awaiting_ground!(data_eng["bus"][eng_obj["bus"]], eng_obj["connections"])
        end

        # _build_time_series_reference!(eng_obj, dss_obj, data_dss, dss_obj, time_series, "pg", "qg")

        if import_all
            _import_all!(eng_obj, dss_obj)
        end

        _add_eng_obj!(data_eng, "generator", id, eng_obj)
    end
end


"Adds vsources to `data_eng` from `data_dss`"
function _dss2eng_vsource!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool, time_series::String="daily")
    for (id, dss_obj) in get(data_dss, "vsource", Dict{String,Any}())
        ph1_ang = dss_obj["angle"]
        vm_pu = dss_obj["pu"]

        phases = dss_obj["phases"]
        vnom = dss_obj["basekv"] / sqrt(phases)

        data_eng["settings"]["vbases_default"][_parse_bus_id(dss_obj["bus1"])[1]] = vnom

        vm = fill(vm_pu, phases)*vnom
        va = rad2deg.(_wrap_to_pi.([-2*pi/phases*(i-1)+deg2rad(ph1_ang) for i in 1:phases]))

        eng_obj = Dict{String,Any}(
            "bus" => _parse_bus_id(dss_obj["bus1"])[1],
            "connections" => _get_conductors_ordered(dss_obj["bus1"]; default=[collect(1:phases)..., 0], pad_ground=true),
            "configuration" => WYE,
            "source_id" => "vsource.$id",
            "status" => dss_obj["enabled"]
        )

        # some values require addition of neutral by default
        n_conductors = length(eng_obj["connections"])
        eng_obj["rs"] = fill(dss_obj["r_mutual"], n_conductors, n_conductors)
        eng_obj["rs"][LinearAlgebra.diagind(eng_obj["rs"])] .= dss_obj["r_self"]

        eng_obj["xs"] = fill(dss_obj["x_mutual"], n_conductors, n_conductors)
        eng_obj["xs"][LinearAlgebra.diagind(eng_obj["xs"])] .= dss_obj["x_self"]

        eng_obj["vm"] = zeros(n_conductors)
        eng_obj["vm"][1:phases] = vm

        eng_obj["va"] = zeros(n_conductors)
        eng_obj["va"][1:phases] = va

        # if the ground is used directly, register load
        if 0 in eng_obj["connections"]
            _register_awaiting_ground!(data_eng["bus"][eng_obj["bus"]], eng_obj["connections"])
        end

        # _build_time_series_reference!(eng_obj, dss_obj, data_dss, time_series, "pg", "qg")

        if import_all
            _import_all!(eng_obj, dss_obj)
        end

        _add_eng_obj!(data_eng, "voltage_source", id, eng_obj)
    end
end


"Adds lines to `data_eng` from `data_dss`"
function _dss2eng_linecode!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool)
    for (id, dss_obj) in get(data_dss, "linecode", Dict{String,Any}())
        eng_obj = Dict{String,Any}()

        nphases = dss_obj["nphases"]

        eng_obj["rs"] = reshape(dss_obj["rmatrix"], nphases, nphases)
        eng_obj["xs"] = reshape(dss_obj["xmatrix"], nphases, nphases)

        eng_obj["b_fr"] = reshape(dss_obj["cmatrix"], nphases, nphases) ./ 2.0
        eng_obj["b_to"] = reshape(dss_obj["cmatrix"], nphases, nphases) ./ 2.0
        eng_obj["g_fr"] = fill(0.0, nphases, nphases)
        eng_obj["g_to"] = fill(0.0, nphases, nphases)

        eng_obj["cm_ub"] = fill(dss_obj["emergamps"], nphases)

        if import_all
            _import_all!(eng_obj, dss_obj)
        end

        _add_eng_obj!(data_eng, "linecode", id, eng_obj)
    end
end


"Adds lines to `data_eng` from `data_dss`"
function _dss2eng_line!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool)
    for (id, dss_obj) in get(data_dss, "line", Dict())
        if haskey(dss_obj, "basefreq") && dss_obj["basefreq"] != data_eng["settings"]["base_frequency"]
            @warn "basefreq=$(dss_obj["basefreq"]) on line.$id does not match circuit basefreq=$(data_eng["settings"]["base_frequency"])"
        end

        nphases = dss_obj["phases"]

        f_connections = _get_conductors_ordered(dss_obj["bus1"], default=collect(1:nphases))
        t_connections = _get_conductors_ordered(dss_obj["bus2"], default=collect(1:nphases))

        ncond = length(f_connections)

        eng_obj = Dict{String,Any}(
            "f_bus" => _parse_bus_id(dss_obj["bus1"])[1],
            "t_bus" => _parse_bus_id(dss_obj["bus2"])[1],
            "length" => dss_obj["switch"] ? 0.001 : dss_obj["length"],
            "f_connections" => f_connections,
            "t_connections" => t_connections,
            "status" => dss_obj["enabled"],
            "source_id" => "line.$id"
        )

        if haskey(dss_obj, "linecode")
            eng_obj["linecode"] = dss_obj["linecode"]
        end

        if (haskey(dss_obj, "emergamps") && _is_after(dss_obj.raw_dss, "emergamps", "linecode")) || isempty(dss_obj["linecode"])
            eng_obj["cm_ub"] = fill(dss_obj["emergamps"], ncond)
        end

        from_geometry = haskey(dss_obj, "geometry") && _is_after(dss_obj.raw_dss, "geometry", "linecode")
        from_spacing = haskey(dss_obj, "wires") && _is_after(dss_obj.raw_dss, "wires", "linecode")
        from_cncables = haskey(dss_obj, "cndata") && _is_after(dss_obj.raw_dss, "cncables", "linecode")
        from_tscables = haskey(dss_obj, "tsdata") && _is_after(dss_obj.raw_dss, "tscables", "linecode")
        if from_geometry || from_spacing || from_cncables || from_tscables
            z, y = calculate_line_constants(data_dss, dss_obj)

            rs, xs = real(z), imag(z)
            g, b = real(y), imag(y)

            eng_obj["rs"] = rs
            eng_obj["xs"] = xs
            eng_obj["b_fr"] = b ./ 2.0
            eng_obj["b_to"] = b ./ 2.0
            eng_obj["g_fr"] = g ./ 2.0
            eng_obj["g_to"] = g ./ 2.0
        end

        if any(haskey(dss_obj, key) && (_is_after(dss_obj.raw_dss, key, "linecode") && _is_after(dss_obj.raw_dss, key, "geometry")) for key in ["r0", "r1", "rg", "rmatrix"]) || (isempty(dss_obj["linecode"]) && isempty(dss_obj["geometry"]) && isempty(dss_obj["wires"]))
            eng_obj["rs"] = reshape(dss_obj["rmatrix"], nphases, nphases)
        end

        if any(haskey(dss_obj, key) && _is_after(dss_obj.raw_dss, key, "linecode") && _is_after(dss_obj.raw_dss, key, "geometry") for key in ["x0", "x1", "xg", "xmatrix"]) || (isempty(dss_obj["linecode"]) && isempty(dss_obj["geometry"]) && isempty(dss_obj["wires"]))
            eng_obj["xs"] = reshape(dss_obj["xmatrix"], nphases, nphases)
        end

        if any(haskey(dss_obj, key) && _is_after(dss_obj.raw_dss, key, "linecode") && _is_after(dss_obj.raw_dss, key, "geometry") for key in ["b0", "b1", "c0", "c1", "cmatrix"]) || (isempty(dss_obj["linecode"]) && isempty(dss_obj["geometry"]) && isempty(dss_obj["wires"]))
            eng_obj["b_fr"] = reshape(dss_obj["cmatrix"], nphases, nphases) ./ 2.0
            eng_obj["b_to"] = reshape(dss_obj["cmatrix"], nphases, nphases) ./ 2.0
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
            _import_all!(eng_obj, dss_obj)
        end

        if dss_obj["switch"]
            eng_obj["state"] = eng_obj["status"] == ENABLED ? CLOSED : OPEN
            eng_obj["dispatchable"] = YES # default
            eng_obj["status"] = ENABLED

            if haskey(eng_obj, "linecode")
                _apply_linecode!(eng_obj, data_eng)
            end
            delete!(eng_obj, "linecode")

            # ENGINEERING model switches are zero-length objects
            for k in ["b_fr", "b_to", "g_fr", "g_to", "rs", "xs"]
                if haskey(eng_obj, k)
                    eng_obj[k] .*= get(eng_obj, "length", 1.0)
                end
            end
            delete!(eng_obj, "length")

            _add_eng_obj!(data_eng, "switch", id, eng_obj)
        else
            _add_eng_obj!(data_eng, "line", id, eng_obj)
        end
    end
end


"Adds transformers to `data_eng` from `data_dss`"
function _dss2eng_xfmrcode!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool, sm_ub::String="emergency")
    @assert sm_ub in ["emergency", "normal"] "Unrecognized sm_ub '$sm_ub'. Must be either 'emergency' or 'normal'"
    sm_ub = sm_ub == "emergency" ? "emerghkva" : "normhkva"

    for (id, dss_obj) in get(data_dss, "xfmrcode", Dict{String,Any}())
        nphases = dss_obj["phases"]
        nrw = dss_obj["windings"]

        eng_obj = Dict{String,Any}(
            "tm_set" => Vector{Vector{Float64}}([fill(tap, nphases) for tap in dss_obj["taps"]]),
            "tm_lb" => Vector{Vector{Float64}}(fill(fill(dss_obj["mintap"], nphases), nrw)),
            "tm_ub" => Vector{Vector{Float64}}(fill(fill(dss_obj["maxtap"], nphases), nrw)),
            "tm_fix" => Vector{Vector{Bool}}(fill(ones(Bool, nphases), nrw)),
            "tm_step" => Vector{Vector{Float64}}(fill(fill(1/32, nphases), nrw)),
            "vm_nom" => Vector{Float64}(dss_obj["kvs"]),
            "sm_nom" => Vector{Float64}(dss_obj["kvas"]),
            "sm_ub" => dss_obj[sm_ub],
            "configuration" => Vector{ConnConfig}(dss_obj["conns"]),
            "rw" => Vector{Float64}(dss_obj["%rs"] ./ 100),
            "noloadloss" => dss_obj["%noloadloss"] / 100,
            "cmag" => dss_obj["%imag"] / 100,
            "xsc" => nrw == 2 ? [dss_obj["xhl"] / 100] : [dss_obj["xhl"], dss_obj["xht"], dss_obj["xlt"]] ./ 100,
            "source_id" => "xfmrcode.$id",
        )

        if import_all
            _import_all!(eng_obj, dss_obj)
        end

        _add_eng_obj!(data_eng, "xfmrcode", id, eng_obj)
    end
end


"Adds transformers to `data_eng` from `data_dss`"
function _dss2eng_transformer!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool, sm_ub::String)
    @assert sm_ub in ["emergency", "normal"] "Unrecognized sm_ub '$sm_ub'. Must be either 'emergency' or 'normal'"
    sm_ub = sm_ub == "emergency" ? "emerghkva" : "normhkva"

    for (id, dss_obj) in get(data_dss, "transformer", Dict{String,Any}())
        eng_obj = Dict{String, Any}(
            "source_id" => "transformer.$id",
            "status" => dss_obj["enabled"],
            "name" => id,
        )

        shared = Dict{String,Any}(
            "phases" => dss_obj.phases,
            "windings" => dss_obj.windings,
            "leadlag" => dss_obj.leadlag,
            "conns" => dss_obj.conns,
        )

        if !isempty(dss_obj.xfmrcode) && !ismissing(data_dss)
            xfmrcode = data_dss.xfmrcode[dss_obj.xfmrcode]

            for key in ["leadlag", "conns", "phases", "windings"]
                if !_is_after(dss_obj.raw_dss, key, "xfmrcode")
                    shared[key] = xfmrcode[key]
                end
            end
        end

        # two-phase delta transformers have single coil
        if all(conf==DELTA for conf in shared["conns"]) && shared["phases"]==2
            ncoils = 1
        else
            ncoils = shared["phases"]
        end

        xfmrcode_idx = findfirst(x->x.first=="xfmrcode", dss_obj.raw_dss)
        xfmrcode_idx = isnothing(xfmrcode_idx) ? Inf : xfmrcode_idx

        wdgs_after_xfmrcode = Int[parse(Int, property_pair.second) for (i,property_pair) in enumerate(dss_obj.raw_dss) if property_pair.first == "wdg" && i > xfmrcode_idx]

        # taps
        if isempty(dss_obj["xfmrcode"]) || _is_after(dss_obj.raw_dss, "taps", "xfmrcode") || all(_is_after(dss_obj.raw_dss, "tap", "xfmrcode", w) for w in 1:shared["windings"])
            eng_obj["tm_set"] = [fill(dss_obj["taps"][w], ncoils) for w in 1:shared["windings"]]
        else
            if !haskey(eng_obj, "tm_set")
                eng_obj["tm_set"] = Vector{Union{Vector{Real},Missing}}(missing, shared["windings"])
            end

            for w in wdgs_after_xfmrcode
                eng_obj["tm_set"][w] = fill(dss_obj.taps[w], ncoils)
            end
        end

        # kvs, kvas
        for (fr_key, to_key) in zip(["kv", "kva"], ["vm_nom", "sm_nom"])
            if isempty(dss_obj.xfmrcode) || _is_after(dss_obj.raw_dss, "$(fr_key)s", "xfmrcode") || all(_is_after(dss_obj.raw_dss, fr_key, "xfmrcode", w) for w in 1:shared["windings"])
                eng_obj[to_key] = dss_obj["$(fr_key)s"]
            else
                for w in wdgs_after_xfmrcode
                    if !haskey(eng_obj, to_key)
                        eng_obj[to_key] = Vector{Union{Real,Missing}}(missing, shared["windings"])
                    end
                    eng_obj[to_key][w] = dss_obj["$(fr_key)s"][w]
                end
            end
        end

        # mintap, maxtap
        for (fr_key, to_key) in zip(["mintap", "maxtap"], ["tm_lb", "tm_ub"])
            if isempty(dss_obj.xfmrcode) || _is_after(dss_obj.raw_dss, "$(fr_key)", "xfmrcode")
                eng_obj[to_key] = fill(fill(dss_obj[fr_key], ncoils), shared["windings"])
            end
        end

        # %noloadloss, %imag
        for (fr_key, to_key) in zip(["%noloadloss", "%imag"], ["noloadloss", "cmag"])
            if isempty(dss_obj.xfmrcode) || _is_after(dss_obj.raw_dss, "$(fr_key)", "xfmrcode")
                eng_obj[to_key] = dss_obj[fr_key] / 100
            end
        end

        # %rs
        if isempty(dss_obj.xfmrcode) || _is_after(dss_obj.raw_dss, "%rs", "xfmrcode") || all(_is_after(dss_obj.raw_dss, "%r", "xfmrcode", w) for w in 1:shared["windings"])
            eng_obj["rw"] = dss_obj["%rs"] / 100
        else
            for w in wdgs_after_xfmrcode
                if !haskey(eng_obj, "rw")
                    eng_obj["rw"] = Vector{Union{Real,Missing}}(missing, shared["windings"])
                end
                eng_obj["rw"][w] = dss_obj["%rs"][w] / 100
            end
        end

        # emerghkva
        if isempty(dss_obj["xfmrcode"]) || (_is_after(dss_obj.raw_dss, sm_ub, "xfmrcode"))
            eng_obj["sm_ub"] = dss_obj[sm_ub]
        end

        # loss model (converted to SI units, referred to secondary)
        if shared["windings"] == 2
            if isempty(dss_obj.xfmrcode) || _is_after(dss_obj.raw_dss, "xhl", "xfmrcode")
                eng_obj["xsc"] = [dss_obj.xhl] / 100
            end
        elseif shared["windings"] == 3
            for (w, key) in enumerate(["xhl", "xht", "xlt"])
                if isempty(dss_obj.xfmrcode) || _is_after(dss_obj.raw_dss, key, "xfmrcode")
                    if !haskey(eng_obj, "xsc")
                        eng_obj["xsc"] = Vector{Union{Float64,Missing}}(missing, 3)
                    end
                    eng_obj["xsc"][w] = dss_obj[key] / 100
                end
            end
        end

        # tm_fix, tm_step don't appear in opendss
        if isempty(dss_obj.xfmrcode)
            eng_obj["tm_fix"] = fill(ones(Bool, ncoils), shared["windings"])
            eng_obj["tm_step"] = fill(fill(1/32, ncoils), shared["windings"])
        end

        # always required
        eng_obj["bus"] = Array{String, 1}(undef, shared["windings"])
        eng_obj["connections"] = Array{Array{Int, 1}, 1}(undef, shared["windings"])
        eng_obj["polarity"] = fill(1, shared["windings"])

        if isempty(dss_obj.xfmrcode)
            eng_obj["configuration"] = shared["conns"]
        end

        # test if this transformer conforms with limitations
        if shared["phases"]<3 && DELTA in shared["conns"]
            # error("Transformers with delta windings should have at least 3 phases to be well-defined: $id.")
        end
        if shared["windings"]>3
            # All of the code is compatible with any number of windings,
            # except for the parsing of the loss model (the pair-wise reactance)
            error("For now parsing of xscarray is not supported. At most 3 windings are allowed, not $nrw.")
        end

        for w in 1:shared["windings"]
            eng_obj["bus"][w] = _parse_bus_id(dss_obj.buses[w])[1]

            conf = shared["conns"][w]
            terminals_default = conf==WYE ? [1:shared["phases"]..., 0] : collect(1:shared["phases"])

            # append ground if connections one too short
            eng_obj["connections"][w] = _get_conductors_ordered(dss_obj.buses[w], default=terminals_default, pad_ground=(conf==WYE))

            if w>1
                prim_conf = shared["conns"][1]
                if shared["leadlag"] in ["ansi", "lag"]
                    if prim_conf==DELTA && conf==WYE
                        eng_obj["polarity"][w] = -1
                        eng_obj["connections"][w] = [_barrel_roll(eng_obj["connections"][w][1:end-1], 1)..., eng_obj["connections"][w][end]]
                    end
                else # hence dss_obj.leadlag in ["euro", "lead"]
                    if prim_conf==WYE && conf==DELTA
                        eng_obj["polarity"][w] = -1
                        eng_obj["connections"][w] = _barrel_roll(eng_obj["connections"][w], -1)
                    end
                end
                if w==3 && eng_obj["connections"][2][2]==0 && eng_obj["connections"][3][1]==0 # center-tap transformers
                    eng_obj["polarity"][w] = -1
                end
            end

            if 0 in eng_obj["connections"][w]
                _register_awaiting_ground!(data_eng["bus"][eng_obj["bus"][w]], eng_obj["connections"][w])
            end
        end

        for key in ["bank", "xfmrcode"]
            if !isempty(dss_obj[key])
                eng_obj[key] = dss_obj[key]
            end
        end

        if import_all
            _import_all!(eng_obj, dss_obj)
        end

        _add_eng_obj!(data_eng, "transformer", id, eng_obj)
    end
end


"Adds pvsystems to `data_eng` from `data_dss`"
function _dss2eng_pvsystem!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool, time_series::String="daily")
    for (id, dss_obj) in get(data_dss, "pvsystem", Dict{String,Any}())
        # TODO pick parameters for solar objects

        nphases = dss_obj["phases"]
        bus = _parse_bus_id(dss_obj["bus1"])[1]
        conf = nphases==1 && dss_obj["kv"]==0.24 ? DELTA : dss_obj["conn"] # check if solar is connected between split-phase terminals of triplex node (nominal line-line voltage=240V), TODO: better generalization

        eng_obj = Dict{String,Any}(
            "bus" => bus,
            "configuration" => conf,
            "connections" => _get_conductors_ordered(dss_obj["bus1"], pad_ground=conf == WYE, default=conf == WYE ? [collect(1:dss_obj["phases"])..., 0] : nphases == 1 ? [1,0] : collect(1:nphases)),
            "pg" => fill(min(dss_obj["pmpp"] * dss_obj["irradiance"], dss_obj["kva"]) / nphases, nphases),
            "qg" => fill(dss_obj["kvar"] / nphases, nphases),
            "vg" => fill(dss_obj["kv"] / sqrt(nphases), nphases),
            "pg_lb" => fill(0.0, nphases),
            "pg_ub" => fill( dss_obj["kva"]  / nphases, nphases),
            "qg_lb" => fill(-dss_obj["kvar"] / nphases, nphases),
            "qg_ub" => fill( dss_obj["kvar"] / nphases, nphases),
            # "sm_ub" => fill(dss_obj["pmpp"] / nphases, nphases), # TODO add irradiance model
            # "irradiance" => dss_obj["irradiance"],
            # "temperature" => dss_obj["temperature"],
            # "p-t_curve" => dss_obj["p-tcurve"],
            # "efficiency_curve" => dss_obj["effcurve"],
            # "rs" => LinearAlgebra.diagm(0 => fill(dss_obj["%r"] / 100., nphases)),
            # "xs" => LinearAlgebra.diagm(0 => fill(dss_obj["%x"] / 100., nphases)),
            "status" => dss_obj["enabled"],
            "source_id" => "pvsystem.$id",
        )

        # if the ground is used directly, register load
        if 0 in eng_obj["connections"]
            _register_awaiting_ground!(data_eng["bus"][eng_obj["bus"]], eng_obj["connections"])
        end

        _build_time_series_reference!(eng_obj, dss_obj, data_dss, time_series, "pg", "qg")
        _build_time_series_reference!(eng_obj, dss_obj, data_dss, time_series, "pg_ub", "qg_ub")

        if import_all
            _import_all!(eng_obj, dss_obj)
        end

        _add_eng_obj!(data_eng, "solar", id, eng_obj)
    end
end


"Adds storage to `data_eng` from `data_dss`"
function _dss2eng_storage!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool, time_series::String="daily")
    for (id, dss_obj) in get(data_dss, "storage", Dict{String,Any}())
        nphases = dss_obj["phases"]

        eng_obj = Dict{String,Any}(
            "bus" => _parse_bus_id(dss_obj["bus1"])[1],
            "connections" => _get_conductors_ordered(dss_obj["bus1"], pad_ground=true, default=[collect(1:dss_obj["phases"])..., 0]),
            "configuration" => WYE,
            "energy" => dss_obj["kwhstored"],
            "energy_ub" => dss_obj["kwhrated"],
            "charge_ub" => dss_obj["%charge"] / 100.0 * dss_obj["kwrated"],
            "discharge_ub" => dss_obj["%discharge"] / 100.0 * dss_obj["kwrated"],
            "sm_ub" => dss_obj["kva"],
            "charge_efficiency" => dss_obj["%effcharge"],
            "discharge_efficiency" => dss_obj["%effdischarge"],
            "qs_lb" => -dss_obj["kva"],  # The storage element can also produce or absorb reactive power (vars) within the kVA rating of the inverter
            "qs_ub" =>  dss_obj["kva"],  # The storage element can also produce or absorb reactive power (vars) within the kVA rating of the inverter
            "rs" => dss_obj["%r"] / 100.0,
            "xs" => dss_obj["%x"] / 100.0,
            "pex" => dss_obj["%idlingkw"] ./ 100.0 .* dss_obj["kwrated"],
            "qex" => dss_obj["%idlingkvar"] ./ 100.0 .* dss_obj["kwrated"], # percent of kwrated consumed as kvar
            "ps" => -dss_obj["kw"],  # PMD convention is opposite from dss
            "qs" => -dss_obj["kvar"],  # PMD convention is opposite from dss
            "status" => dss_obj["enabled"],
            "source_id" => "storage.$id",
        )

        # if the ground is used directly, register load
        if 0 in eng_obj["connections"]
            _register_awaiting_ground!(data_eng["bus"][eng_obj["bus"]], eng_obj["connections"])
        end

        if import_all
            _import_all!(eng_obj, dss_obj)
        end

        _add_eng_obj!(data_eng, "storage", id, eng_obj)
    end
end


"Adds regcontrol to `data_eng` from `data_dss`"
function _dss2eng_regcontrol!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool)
    for (id, dss_obj) in get(data_dss, "regcontrol", Dict{String,Any}())

        nrw = get(data_dss["transformer"]["$(dss_obj["transformer"])"],"windings",2)
        nphases = get(data_dss["transformer"]["$(dss_obj["transformer"])"], "phases", 3)

        eng_obj = Dict{String,Any}(
            "vreg" => [[w == dss_obj["winding"] && p == dss_obj["ptphase"] ? dss_obj["vreg"] : 0.0 for p in 1:nphases] for w in 1:nrw],
            "band" => [[w == dss_obj["winding"] && p == dss_obj["ptphase"] ? dss_obj["band"] : 0.0 for p in 1:nphases] for w in 1:nrw],
            "ptratio" => [[w == dss_obj["winding"] && p == dss_obj["ptphase"] ? dss_obj["ptratio"] : 0.0 for p in 1:nphases] for w in 1:nrw],
            "ctprim" => [[w == dss_obj["winding"] && p == dss_obj["ptphase"] ? dss_obj["ctprim"] : 0.0 for p in 1:nphases] for w in 1:nrw],
            "r" => [[w == dss_obj["winding"] && p == dss_obj["ptphase"] ? dss_obj["r"] : 0.0 for p in 1:nphases] for w in 1:nrw],
            "x" => [[w == dss_obj["winding"] && p == dss_obj["ptphase"] ? dss_obj["x"] : 0.0 for p in 1:nphases] for w in 1:nrw]
            )

        if import_all
            _import_all!(eng_obj, dss_obj)
        end

        # add regcontrol items to transformer if present
        data_eng["transformer"]["$(dss_obj["transformer"])"]["controls"] = eng_obj
        if haskey(data_eng["transformer"]["$(dss_obj["transformer"])"],"tm_fix")
            data_eng["transformer"]["$(dss_obj["transformer"])"]["tm_fix"] = [[w == dss_obj["winding"] && p == dss_obj["ptphase"] ? false : true for p in 1:nphases] for w in 1:nrw]
        end
    end
end


"Adds capcontrol to `data_eng` from `data_dss`"
function _dss2eng_capcontrol!(data_eng::Dict{String,<:Any}, data_dss::OpenDssDataModel, import_all::Bool)
    for (id, dss_obj) in get(data_dss, "capcontrol", Dict{String,Any}())
        type = dss_obj["type"]
        nphases = data_dss["capacitor"]["$(dss_obj["capacitor"])"]["phases"]

        if !haskey(data_eng["shunt"]["$(dss_obj["capacitor"])"],"controls")
            eng_obj = Dict{String,Any}("element" => dss_obj["element"],
                    "voltoverride" => type == CAP_REACTIVE_POWER ? dss_obj["voltoverride"] : [p == dss_obj["ptphase"] ? dss_obj["voltoverride"] : false for p in 1:nphases]
                    )

            if type == CAP_REACTIVE_POWER
                eng_obj["type"] = type
                eng_obj["terminal"] =  dss_obj["terminal"]
                eng_obj["onsetting"] =  dss_obj["onsetting"]
                eng_obj["offsetting"] =  dss_obj["offsetting"]
            elseif type == CAP_VOLTAGE
                eng_obj["type"] = [p == dss_obj["ptphase"] ? type : _dss2pmd_capcontrol_type[""] for p in 1:nphases]
                eng_obj["terminal"] = [p == dss_obj["ptphase"] ? dss_obj["terminal"] : 0 for p in 1:nphases]
                eng_obj["onsetting"] = [p == dss_obj["ptphase"] ? dss_obj["onsetting"] : 0.0 for p in 1:nphases]
                eng_obj["offsetting"] = [p == dss_obj["ptphase"] ? dss_obj["offsetting"] : 0.0 for p in 1:nphases]
                eng_obj["ptratio"] = [p == dss_obj["ptphase"] ? dss_obj["ptratio"] : 0.0 for p in 1:nphases]
            elseif type == CAP_CURRENT
                eng_obj["type"] = [p == dss_obj["ctphase"] ? type : _dss2pmd_capcontrol_type[""] for p in 1:nphases]
                eng_obj["terminal"] = [p == dss_obj["ctphase"] ? dss_obj["terminal"] : 0 for p in 1:nphases]
                eng_obj["onsetting"] = [p == dss_obj["ctphase"] ? dss_obj["onsetting"] : 0.0 for p in 1:nphases]
                eng_obj["offsetting"] = [p == dss_obj["ctphase"] ? dss_obj["offsetting"] : 0.0 for p in 1:nphases]
                eng_obj["ctratio"] = [p == dss_obj["ctphase"] ? dss_obj["ctratio"] : 0.0 for p in 1:nphases]
            elseif type == CAP_TIME
                eng_obj["type"] = type
                eng_obj["terminal"] =  dss_obj["terminal"]
            end

            if dss_obj["voltoverride"]
                eng_obj["vmin"] = type == CAP_REACTIVE_POWER ? dss_obj["vmin"] : [p == dss_obj["ptphase"] ? dss_obj["vmin"] : 0.0 for p in 1:nphases]
                eng_obj["vmax"] = type == CAP_REACTIVE_POWER ? dss_obj["vmax"] : [p == dss_obj["ptphase"] ? dss_obj["vmax"] : 0.0 for p in 1:nphases]
                eng_obj["ptratio"] = type == CAP_REACTIVE_POWER ? dss_obj["ptratio"] : [p == dss_obj["ptphase"] ? dss_obj["ptratio"] : 0.0 for p in 1:nphases]
            end

            if import_all
                _import_all!(eng_obj, dss_obj)
            end

            # add capcontrol items to capacitor if present
            data_eng["shunt"]["$(dss_obj["capacitor"])"]["controls"] = eng_obj
        else
            eng_obj = data_eng["shunt"]["$(dss_obj["capacitor"])"]["controls"]
            if type == CAP_VOLTAGE
                eng_obj["type"][dss_obj["ptphase"]] = type
                eng_obj["terminal"][dss_obj["ptphase"]] = dss_obj["terminal"]
                eng_obj["onsetting"][dss_obj["ptphase"]] = dss_obj["onsetting"]
                eng_obj["offsetting"][dss_obj["ptphase"]] = dss_obj["offsetting"]
                eng_obj["ptratio"][dss_obj["ptphase"]] = dss_obj["ptratio"]
            end
            if type == CAP_CURRENT
                eng_obj["type"][dss_obj["ctphase"]] = type
                eng_obj["terminal"][dss_obj["ctphase"]] = dss_obj["terminal"]
                eng_obj["onsetting"][dss_obj["ctphase"]] = dss_obj["onsetting"]
                eng_obj["offsetting"][dss_obj["ctphase"]] = dss_obj["offsetting"]
                eng_obj["ctratio"][dss_obj["ptphase"]] = dss_obj["ctratio"]
            end
            if type == CAP_TIME
                eng_obj["type"] = type
                eng_obj["terminal"] = dss_obj["terminal"]
            end

            if dss_obj["voltoverride"]
                eng_obj["voltoverride"][dss_obj["ptphase"]] = dss_obj["voltoverride"]
                eng_obj["ptratio"][dss_obj["ptphase"]] = dss_obj["ptratio"]
                eng_obj["vmin"][dss_obj["ptphase"]] = dss_obj["vmin"]
                eng_obj["vmax"][dss_obj["ptphase"]] = dss_obj["vmax"]
            end

            if import_all
                _import_all!(eng_obj, dss_obj)
            end
        end
    end
end


"""
    parse_opendss(
        io::IO;
        import_all::Bool=false,
        bank_transformers::Bool=true,
        time_series::String="daily",
        dss2eng_extensions::Vector{<:Function}=Function[],
    )::Dict{String,Any}

Parses an IO, into raw dss dictionary via [`parse_dss`](@ref parse_dss), into the `ENGINEERING` [`DataModel`](@ref DataModel)

See [`parse_opendss`](@ref parse_opendss)
"""
function parse_opendss(
    io::IO;
    import_all::Bool=false,
    bank_transformers::Bool=true,
    time_series::String="daily",
    dss2eng_extensions::Vector{<:Function}=Function[],
    )::Dict{String,Any}

    data_dss = parse_dss(io)

    return parse_opendss(data_dss;
        import_all=import_all,
        bank_transformers=bank_transformers,
        time_series=time_series,
        dss2eng_extensions=dss2eng_extensions,
    )
end


"""
    parse_opendss(
        data_dss::OpenDssDataModel;
        import_all::Bool=false,
        bank_transformers::Bool=true,
        time_series::String="daily",
        dss2eng_extensions::Vector{<:Function}=Function[]
    )::Dict{String,Any}

Parses a raw dss data structure (dictionary), resulting from the parsing of a `DSS` file, into the `ENGINEERING` [`DataModel`](@ref DataModel)

If `import_all` is true, all raw dss properties will be included in the final dictionary under `"dss"`.

If `bank_transformers` is true (default), transformers that are indicated to be part of a bank in dss will be combined into a single multiphase
transformer.

`time_series` defines which property the time series will be taken from, `"daily"` or "yearly". More complex parsing of time series data
should be performed with `dss2eng_extensions`.

# `dss2eng_extensions`

If a user wishes to parse additional components that are not yet natively supported by PowerModelsDistribution, `dss2eng_extensions` can
be utilized. Custom user functions provided under `dss2eng_extensions` will be excuted __after__ all built-in dss2eng transformations
have been performed and transformers have been banked together (if `bank_transformers==true`). dss2eng_extension functions should have the
following function signature:

    dss2eng_func!(data_eng, data_dss)

where `data_eng` is a non-multinetwork ENGINEERING data model (_i.e._, time series data has not yet been expanded into a multinetwork
structure), and `data_dss` is the raw dss data parsed by [`parse_dss`](@ref parse_dss).
"""
function parse_opendss(
    data_dss::OpenDssDataModel;
    import_all::Bool=false,
    bank_transformers::Bool=true,
    time_series::String="daily",
    dss2eng_extensions::Vector{<:Function}=Function[],
    )::Dict{String,Any}

    data_eng = Dict{String,Any}(
        "data_model" => ENGINEERING,
        "settings" => Dict{String,Any}(),
    )

    if import_all
        data_eng["dss_options"] = data_dss["options"]
    end

    if haskey(data_dss, "vsource") && haskey(data_dss["vsource"], "source") && !isempty(data_dss["circuit"])
        dss_obj = data_dss["vsource"]["source"]
        source_bus = _parse_bus_id(dss_obj["bus1"])[1]

        data_eng["name"] = first(data_dss["circuit"]).first

        data_eng["settings"]["voltage_scale_factor"] = 1e3
        data_eng["settings"]["power_scale_factor"] = 1e3
        data_eng["settings"]["vbases_default"] = Dict{String,Real}()
        dss_obj["basemva"] == 100.0 && @info "basemva=100 is the default value, you may want to adjust sbase_default for better convergence"
        data_eng["settings"]["sbase_default"] = dss_obj["basemva"] * 1e3
        data_eng["settings"]["base_frequency"] = data_dss.options.defaultbasefrequency

        # collect turns the Set into Array, making it serializable
        data_eng["files"] = collect(data_dss["filename"])
    else
        error("Circuit not defined, not a valid circuit!")
    end

    _dss2eng_bus!(data_eng, data_dss, import_all)
    _dss2eng_buscoords!(data_eng, data_dss)

    _dss2eng_linecode!(data_eng, data_dss, import_all)
    _dss2eng_line!(data_eng, data_dss, import_all)

    _dss2eng_xfmrcode!(data_eng, data_dss, import_all, "emergency")
    _dss2eng_transformer!(data_eng, data_dss, import_all, "emergency")

    _dss2eng_capacitor!(data_eng, data_dss, import_all)
    _dss2eng_reactor!(data_eng, data_dss, import_all)

    _dss2eng_regcontrol!(data_eng, data_dss, import_all)
    _dss2eng_capcontrol!(data_eng, data_dss, import_all)

    _dss2eng_loadshape!(data_eng, data_dss, import_all)
    _dss2eng_load!(data_eng, data_dss, import_all, time_series)

    _dss2eng_vsource!(data_eng, data_dss, import_all, time_series)
    _dss2eng_generator!(data_eng, data_dss, import_all, time_series)
    _dss2eng_pvsystem!(data_eng, data_dss, import_all, time_series)
    _dss2eng_storage!(data_eng, data_dss, import_all, time_series)

    _discover_terminals!(data_eng)

    if bank_transformers
        _bank_transformers!(data_eng)
    end

    for dss2eng_func! in dss2eng_extensions
        dss2eng_func!(data_eng, data_dss)
    end

    find_conductor_ids!(data_eng)

    return data_eng
end


"""
    add_voltage_starts!(eng::Dict{String,<:Any}, voltages::Dict{String,<:Any})

Function to add vm_start and va_start properties to buses from a voltages dictionary with the formats

```julia
Dict{String,Any}(
    "bus" => Dict{String,Any}(
        "terminals" => Int[],
        "vm" => Real[],
        "va" => Real[],
        "vbase" => Real,
    )
)
```

`"vm_start"`, `"va_start"`, and `"vbase"` are expected to be in SI units. `"vbase"` is optional.
"""
function add_voltage_starts!(eng::Dict{String,<:Any}, voltages::Dict{String,<:Any})
    for (bus_id, obj) in voltages
        eng_obj = eng["bus"][bus_id]

        eng_obj["vm_start"] = Real[t in obj["terminals"] ? obj["vm"][findfirst(isequal(t), obj["terminals"])] : 0.0 for t in eng_obj["terminals"]] ./ eng["settings"]["voltage_scale_factor"]
        eng_obj["va_start"] = Real[t in obj["terminals"] ? obj["va"][findfirst(isequal(t), obj["terminals"])] : 0.0 for t in eng_obj["terminals"]]
        if haskey(obj, "vbase")
            eng_obj["vbase"] = obj["vbase"] ./ eng["settings"]["voltage_scale_factor"]
        end
    end
end


"""
    add_voltage_starts!(eng::Dict{String,<:Any}, voltages_file::String)

Function to add `vm_start` and `va_start` properties to buses from a voltages csv file exported from OpenDSS,
using [`parse_dss_voltages_export`](@ref parse_dss_voltages_export)
"""
function add_voltage_starts!(eng::Dict{String,<:Any}, voltages_file::String)
    add_voltage_starts!(eng, parse_dss_voltages_export(voltages_file))
end
