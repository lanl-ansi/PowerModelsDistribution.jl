"""
"""
function create_eng_object(::Type{T}, dss_obj::DssCapacitor; import_all::Bool=false)::T where T <: EngLine
    @info "treating capacitor.$(dss_obj) like line"

    nphases = dss_obj.phases
    conn = dss_obj.conn

    f_terminals = _get_conductors_ordered(dss_obj.bus1, default=collect(1:nphases))
    if conn==WYE
        t_terminals = _get_conductors_ordered(dss_obj.bus2, default=fill(0,nphases))
    else
        # if delta connected, ignore bus2 and generate t_terminals such that
        # it corresponds to a delta winding
        t_terminals = [f_terminals[2:end]..., f_terminals[1]]
    end

    bus1, _ = _parse_bus_id(dss_obj.bus1)
    bus2, _ = _parse_bus_id(dss_obj.bus2)

    T(;
        f_bus = bus1,
        t_bus = bus2,
        f_connections=_get_conductors_ordered(dss_obj.bus1, default=collect(1:nphases)),
        t_connections=_get_conductors_ordered(dss_obj.bus2, default=collect(1:nphases)),
        length=1.0,
        rs=diagm(0 => dss_obj.r),
        xs=diagm(0 => dss_obj.xl),
        g_fr=zeros(nphases,nphases),
        b_fr=zeros(nphases,nphases),
        g_to=zeros(nphases,nphases),
        b_to=zeros(nphases,nphases),
        status=dss_obj.enabled,
        source_id="capacitor.$(dss_obj.name)",
        dss = import_all ? dss_obj : missing,
    )
end


"""
"""
function create_eng_object(::Type{T}, dss_obj::DssReactor; import_all::Bool=false)::T where T <: EngLine
    @warn "treating reactor.$(dss_obj.name) like line"

    nphases = dss_obj.phases

    T(;
        name = dss_obj.name,
        f_bus = _parse_bus_id(dss_obj.bus1)[1],
        t_bus = _parse_bus_id(dss_obj.bus2)[1],
        f_connections = _get_conductors_ordered(dss_obj.bus1, default=collect(1:nphases)),
        t_connections = _get_conductors_ordered(dss_obj.bus2, default=collect(1:nphases)),
        length = 1.0,
        rs = dss_obj.rmatrix,
        xs = dss_obj.xmatrix,
        g_fr = zeros(nphases, nphases),
        b_fr = zeros(nphases, nphases),
        g_to = zeros(nphases, nphases),
        b_to = zeros(nphases, nphases),
        status = dss_obj.enabled,
        source_id = "reactor.$(dss_obj.name)",
        dss = import_all ? dss_obj : missing,
    )
end


"""
"""
function create_eng_object(::Type{T}, dss_obj::DssLine; import_all::Bool=false, cm_ub::String="emergency", dss::Union{Missing,OpenDssDataModel}=missing)::T where T <: EngLine
    @assert cm_ub in ["emergency", "normal"] "Unrecognized cm_ub '$cm_ub'. Must be either 'emergency' or 'normal'"

    cm_ub = cm_ub == "emergency" ? "emergamps" : "normamps"

    nphases = dss_obj.phases

    f_connections = _get_conductors_ordered(dss_obj.bus1, default=collect(1:nphases))
    t_connections = _get_conductors_ordered(dss_obj.bus2, default=collect(1:nphases))

    ncond = length(f_connections)

    eng_obj = T(;
        name = dss_obj.name,
        f_bus = _parse_bus_id(dss_obj.bus1)[1],
        t_bus = _parse_bus_id(dss_obj.bus2)[1],
        length = dss_obj.length,
        f_connections = f_connections,
        t_connections = t_connections,
        status = dss_obj.enabled,
        linecode = dss_obj.linecode,
        source_id = "line.$(dss_obj.name)",
        dss = import_all ? dss_obj : missing,
    )

    if isempty(eng_obj.linecode) || _is_after(dss_obj.raw_dss, "cm_ub", "linecode")
        eng_obj.cm_ub = fill(dss_obj[cm_ub], ncond)
    end

    from_geometry = !isempty(dss_obj.geometry) && _is_after(dss_obj.raw_dss, "geometry", "linecode")
    from_spacing = !isempty(dss_obj.wires) && _is_after(dss_obj.raw_dss, "wires", "linecode")
    from_cncables = !isempty(dss_obj.cncables) && _is_after(dss_obj.raw_dss, "cncables", "linecode")
    from_tscables = !isempty(dss_obj.tscables) && _is_after(dss_obj.raw_dss, "tscables", "linecode")
    if from_geometry || from_spacing || from_cncables || from_tscables
        z, y = calculate_line_constants(dss, dss_obj)

        rs, xs = real(z), imag(z)
        g, b = real(y), imag(y)

        eng_obj.rs = rs
        eng_obj.xs = xs
        eng_obj.b_fr = b ./ 2.0
        eng_obj.b_to = b ./ 2.0
        eng_obj.g_fr = g ./ 2.0
        eng_obj.g_to = g ./ 2.0
    end


    if any(!isempty(getproperty(dss_obj, Symbol(key))) && (_is_after(dss_obj.raw_dss, key, "linecode") && _is_after(dss_obj.raw_dss, key, "geometry")) for key in ["r0", "r1", "rg", "rmatrix"]) || (isempty(dss_obj.linecode) && isempty(dss_obj.geometry) && isempty(dss_obj.wires))
        eng_obj.rs = dss_obj.rmatrix
    end

    if any(!isempty(getproperty(dss_obj, Symbol(key))) && (_is_after(dss_obj.raw_dss, key, "linecode") && _is_after(dss_obj.raw_dss, key, "geometry")) for key in ["x0", "x1", "xg", "xmatrix"]) || (isempty(dss_obj.linecode) && isempty(dss_obj.geometry) && isempty(dss_obj.wires))
        eng_obj.xs = dss_obj.xmatrix
    end

    if any(!isempty(getproperty(dss_obj, Symbol(key))) && (_is_after(dss_obj.raw_dss, key, "linecode") && _is_after(dss_obj.raw_dss, key, "geometry")) for key in ["b0", "b1", "c0", "c1", "cmatrix"]) || (isempty(dss_obj.linecode) && isempty(dss_obj.geometry) && isempty(dss_obj.wires))
        eng_obj.b_fr = dss_obj.cmatrix ./ 2.0
        eng_obj.b_to = dss_obj.cmatrix ./ 2.0
        eng_obj.g_fr = zeros(Float64, ncond, ncond)
        eng_obj.g_to = zeros(Float64, ncond, ncond)
    end

    return eng_obj
end


"""
"""
function create_eng_object(::Type{T}, dss_obj::DssLine; import_all::Bool=false, cm_ub::String="emergency")::T where T <: EngSwitch
    @assert cm_ub in ["emergency", "normal"] "Unrecognized cm_ub '$cm_ub'. Must be either 'emergency' or 'normal'"

    cm_ub = cm_ub == "emergency" ? "emergamps" : "normamps"

    nphases = dss_obj.phases

    f_connections = _get_conductors_ordered(dss_obj.bus1, default=collect(1:nphases))
    t_connections = _get_conductors_ordered(dss_obj.bus2, default=collect(1:nphases))

    ncond = length(f_connections)

    eng_obj = T(;
        name = dss_obj.name,
        f_bus = _parse_bus_id(dss_obj.bus1)[1],
        t_bus = _parse_bus_id(dss_obj.bus2)[1],
        length = 0.001,
        f_connections = f_connections,
        t_connections = t_connections,
        status = dss_obj.enabled,
        linecode = dss_obj.linecode,
        state = CLOSED,
        dispatchable = YES,
        source_id = "line.$(dss_obj.name)",
        dss = import_all ? dss_obj : missing,
    )

    if isempty(eng_obj.linecode) || _is_after(dss_obj.raw_dss, cm_ub, "linecode")
        eng_obj.cm_ub = fill(dss_obj[cm_ub], ncond)
    end

    if isempty(eng_obj.linecode) || any(_is_after(dss_obj.raw_dss, prop, "linecode") for prop in ["r0", "r1", "rg", "rmatrix"])
        eng_obj.rs = dss_obj.rmatrix
    end

    if isempty(eng_obj.linecode) || any(_is_after(dss_obj.raw_dss, prop, "linecode") for prop in ["x0", "x1", "xg", "xmatrix"])
        eng_obj.xs = dss_obj.xmatrix
    end

    return eng_obj
end


"""
"""
function create_eng_object(::Type{T}, dss_obj::DssTransformer; import_all::Bool=false, sm_ub::String="emergency", dss::Union{Missing,OpenDssDataModel}=missing)::T where T <: EngTransformer
    shared = Dict{String,Any}(
        "phases" => dss_obj.phases,
        "windings" => dss_obj.windings,
        "leadlag" => dss_obj.leadlag,
        "conns" => dss_obj.conns,
    )

    @assert sm_ub in ["emergency", "normal"] "Unrecognized cm_ub '$cm_ub'. Must be either 'emergency' or 'normal'"

    sm_ub = sm_ub == "emergency" ? "emerghkva" : "normhkva"


    if !isempty(dss_obj.xfmrcode) && !ismissing(dss)
        xfmrcode = dss.xfmrcode[dss_obj.xfmrcode]

        for key in ["leadlag", "conns", "phases", "windings"]
            if !_is_after(dss_obj.raw_dss, key, "xfmrcode")
                shared[key] = xfmrcode[key]
            end
        end
    end

    eng_obj = Dict{String,Any}()

    # two-phase delta transformers have single coil
    if all(conf==DELTA for conf in shared["conns"]) && shared["phases"]==2
        ncoils = 1
    else
        ncoils = shared["phases"]
    end

    xfmrcode_idx = findfirst(x->x.first=="xfmrcode", dss_obj.raw_dss)
    xfmrcode_idx = isnothing(xfmrcode_idx) ? Inf : xfmrcode_idx

    wdgs_after_xfmrcode = Int[property_pair.second for (i,property_pair) in enumerate(dss_obj.raw_dss) if property_pair.first == "wdg" && i > xfmrcode_idx]

    # taps
    if isempty(dss_obj["xfmrcode"]) || _is_after(dss_obj.raw_dss, "taps", "xfmrcode") || all(_is_after(dss_obj.raw_dss, "tap", "xfmrcode", w) for w in 1:shared["windings"])
        eng_obj["tm_set"] = [fill(dss_obj["taps"][w], ncoils) for w in 1:shared["windings"]]
    else
        if !haskey(eng_obj, "tm_set")
            eng_obj["tm_set"] = fill(fill(NaN, ncoils), shared["windings"])
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
                    eng_obj[to_key] = fill(NaN, shared["windings"])
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
                eng_obj["rw"] = fill(NaN, shared["windings"])
            end
            eng_obj["rw"][w] = dss_obj["%rs"][w] / 100
        end
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
                    eng_obj["xsc"] = Vector{Float64}(NaN, 3)
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
        eng_obj["configurations"] = shared["conns"]
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
    end

    eng_obj = T(;
        name = dss_obj.name,
        bus = eng_obj["bus"],
        connections = eng_obj["connections"],
        configurations = get(eng_obj, "configurations", missing),
        xfmrcode = dss_obj.xfmrcode,
        xsc = get(eng_obj, "xsc", missing),
        rw = get(eng_obj, "rw", missing),
        cmag = get(eng_obj, "cmag", missing),
        noloadloss = get(eng_obj, "noloadloss", missing),
        tm_nom = get(eng_obj, "tm_nom", missing),
        tm_ub = get(eng_obj, "tm_ub", missing),
        tm_lb = get(eng_obj, "tm_lb", missing),
        tm_set = get(eng_obj, "tm_set", missing),
        tm_fix = get(eng_obj, "tm_fix", missing),
        tm_step = get(eng_obj, "tm_step", missing),
        vm_nom = get(eng_obj, "vm_nom", missing),
        sm_nom = get(eng_obj, "sm_nom", missing),
        polarity = eng_obj["polarity"],
        bank = dss_obj.bank,
        status = dss_obj.enabled,
        source_id = "transformer.$(dss_obj.name)",
        dss = import_all ? dss_obj : missing,
    )

    if isempty(eng_obj.xfmrcode) || _is_after(dss_obj.raw_dss, sm_ub, "xfmrcode")
        eng_obj.sm_ub = dss_obj[sm_ub]
    end

    return eng_obj
end
