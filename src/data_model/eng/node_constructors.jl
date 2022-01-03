"""
Creates EngLoad objects

Constant can still be scaled by other settings, fixed cannot
Note that in the current feature set, fixed therefore equals constant

1: Constant P and Q (default)
2: Constant Z
3: Constant P and quadratic Q
4: Exponential
5: Constant I
6: Constant P and fixed Q
7: Constant P and quadratic Q (i.e., fixed reactance) (Not yet supported)
8: ZIP (Not yet supported)
"""
function create_eng_object(::Type{T}, dss_obj::DssLoad; import_all::Bool=false, time_series::String="daily")::T where T <: EngLoad
    nphases = dss_obj["phases"]
    conf = dss_obj["conn"]

    if conf == DELTA
        @assert nphases in [1, 3] "$(dss_obj.name): only 1 and 3-phase delta loads are supported!"
    end

    # connections
    bus, _ = _parse_bus_id(dss_obj.bus1)
    connections_default = conf == WYE ? [collect(1:nphases)..., 0] : nphases==1 ? [1,2] : [1,2,3]
    connections = _get_conductors_ordered(dss_obj.bus1, default=connections_default, pad_ground=(conf==WYE))

    @assert length(unique(connections))==length(connections) "$(dss_obj.name): connections cannot be made to a terminal more than once."

    kv = dss_obj.kv
    if conf==WYE && nphases in [2, 3]
        kv = kv/sqrt(3)
    end

    T(;
        name = dss_obj.name,
        bus=bus,
        model=dss_obj.model,
        configuration=conf,
        connections=connections,
        dispatchable=NO,
        source_id="load.$(dss_obj.name)",
        status=dss_obj.enabled,
        vm_nom = kv,
        pd_nom = fill(dss_obj.kw/nphases, nphases),
        qd_nom = fill(dss_obj.kvar/nphases, nphases),
        time_series = Dict{String,String}(
            "pd_nom" => getproperty(dss_obj, Symbol(time_series)),
            "qd_nom" => getproperty(dss_obj, Symbol(time_series)),
        ),
        dss = import_all ? dss_obj : missing,
    )
end


"""
"""
function create_eng_object(::Type{T}, dss_obj::DssCapacitor; import_all::Bool)::T where T <: EngShunt
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

    bus_name, _ = _parse_bus_id(dss_obj.bus1)
    # 'kv' is specified as phase-to-phase for phases=2/3 (unsure for 4 and more)
    #TODO figure out for more than 3 phases
    vnom_ln = dss_obj.kv
    if dss_obj.phases in [2,3]
        vnom_ln = vnom_ln/sqrt(3)
    end

    kv = fill(vnom_ln, nphases)

    # 'kvar' is specified for all phases at once; we want the per-phase one
    kvar = fill(dss_obj.kvar / nphases, nphases)

    vnom_ln = kv
    qnom = kvar ./ 1e3
    b = qnom ./ vnom_ln.^2

    # convert to a shunt matrix
    terminals, B = _calc_shunt(f_terminals, t_terminals, b)

    # if one terminal is ground (0), reduce shunt admittance matrix
    terminals, B = _calc_ground_shunt_admittance_matrix(terminals, B, 0)

    T(;
        name = dss_obj.name,
        bus=bus_name,
        model=CAPACITOR,
        dispatchable=NO,
        status=dss_obj.enabled,
        source_id="capacitor.$(dss_obj.name)",
        gs=zeros(size(B)),
        bs=B,
        connections=terminals,
        dss = import_all ? dss_obj : missing,
    )
end


"""
"""
function create_eng_object(::Type{T}, dss_obj::DssReactor; import_all::Bool=false)::T where T <: EngShunt
    nphases = dss_obj.phases
    bus, _ = _parse_bus_id(dss_obj.bus1)

    connections_default = dss_obj.conn == WYE ? [collect(1:nphases)..., 0] : collect(1:nphases)
    connections = _get_conductors_ordered(dss_obj.bus1, default=connections_default, check_length=false)

    # TODO Check unit conversion on Gcap
    Gcap = sum(dss_obj.kvar) / (nphases * 1e3 * (dss_obj.kv / sqrt(nphases))^2)

    T(;
        name = dss_obj.name,
        bus = bus,
        configuration = dss_obj.conn,
        connections = connections,
        model = REACTOR,
        dispatchable = NO,
        status = dss_obj.enabled,
        source_id = "reactor.$(dss_obj.name)",
        bs = diagm(0 => fill(Gcap, nphases)),
        gs = zeros(nphases, nphases),
        dss = import_all ? dss_obj : missing,
    )
end


"""
"""
function create_eng_object(::Type{T}, dss_obj::DssGenerator; import_all::Bool=false, time_series::String="daily")::T where T <: EngGenerator
    nphases = dss_obj.phases

    T(;
        name = dss_obj.name,
        connections = _get_conductors_ordered(dss_obj.bus1, pad_ground=dss_obj.conn == WYE, default=dss_obj.conn == WYE ? [collect(1:dss_obj.phases)..., 0] : nphases == 1 ? [1,0] : collect(1:nphases)),
        bus = _parse_bus_id(dss_obj.bus1)[1],
        pg = fill(dss_obj.kw / nphases, nphases),
        qg = fill(dss_obj.kvar / nphases, nphases),
        vg = fill(dss_obj.kv / sqrt(nphases), nphases),
        qg_lb = fill(dss_obj.minkvar / nphases, nphases),
        qg_ub = fill(dss_obj.maxkvar / nphases, nphases),
        pg_lb = fill(0.0, nphases),
        pg_ub = fill(dss_obj.kw / nphases, nphases),
        control_mode = FREQUENCYDROOP,
        configuration = WYE,
        status = dss_obj.enabled,
        cost_pg_model = 2,
        cost_pg_parameters = [0.0, 1.0, 0.0],
        source_id = "generator.$(dss_obj.name)",
        time_series = Dict{String,String}(
            "pg" => getproperty(dss_obj, Symbol(time_series)),
            "qg" => getproperty(dss_obj, Symbol(time_series))
        ),
        dss = import_all ? dss_obj : missing,
    )
end


"""
"""
function create_eng_object(::Type{T}, dss_obj::DssVsource; import_all::Bool=false)::T where T <: EngVoltageSource
    ph1_ang = dss_obj.angle
    vm_pu = dss_obj.pu

    phases = dss_obj.phases
    vnom = dss_obj.basekv / sqrt(phases)

    vm = fill(vm_pu, phases)*vnom
    va = rad2deg.(_wrap_to_pi.([-2*pi/phases*(i-1)+deg2rad(ph1_ang) for i in 1:phases]))

    connections = _get_conductors_ordered(dss_obj.bus1; default=[collect(1:phases)..., 0], pad_ground=true)
    n_conductors = length(connections)

    if dss_obj.name == "source" && isempty(dss_obj.bus1)
        dss_obj.bus1 = "sourcebus"
    end

    eng_obj = EngVoltageSource(;
        name = dss_obj.name,
        bus = _parse_bus_id(dss_obj.bus1)[1],
        connections = connections,
        source_id = "vsource.$(dss_obj.name)",
        status = dss_obj.enabled,
        rs = zeros(n_conductors, n_conductors),
        xs = zeros(n_conductors, n_conductors),
        vm = zeros(n_conductors),
        va = zeros(n_conductors),
        time_series = Dict{String,String}(),
        dss = import_all ? dss_obj : missing,
    )

    # some values require addition of neutral by default
    eng_obj.rs[1:phases, 1:phases] .= dss_obj.rmatrix
    eng_obj.xs[1:phases, 1:phases] .= dss_obj.xmatrix
    eng_obj.vm[1:phases] .= vm
    eng_obj.va[1:phases] .= va

    return eng_obj
end


"""
"""
function create_eng_object(::Type{T}, dss_obj::DssPvsystem; import_all::Bool=false, time_series="daily")::T where T <: EngSolar
    # TODO pick parameters for solar objects

    nphases = dss_obj.phases

    T(;
        name = dss_obj.name,
        bus = _parse_bus_id(dss_obj.bus1)[1],
        configuration = dss_obj.conn,
        connections = _get_conductors_ordered(dss_obj.bus1, pad_ground=dss_obj.conn == WYE, default=dss_obj.conn == WYE ? [collect(1:dss_obj.phases)..., 0] : nphases == 1 ? [1,0] : collect(1:nphases)),
        pg = fill(min(dss_obj.pmpp * dss_obj.irradiance, dss_obj.kva) / nphases, nphases),
        qg = fill(dss_obj.kvar / nphases, nphases),
        vg = fill(dss_obj.kv / sqrt(nphases), nphases),
        pg_lb = fill(0.0, nphases),
        pg_ub = fill( dss_obj.kva  / nphases, nphases),
        qg_lb = fill(-dss_obj.kvar / nphases, nphases),
        qg_ub = fill( dss_obj.kvar / nphases, nphases),
        # sm_ub = fill(dss_obj.pmpp / nphases, nphases), # TODO add irradiance model
        # irradiance = dss_obj.irradiance,
        # temperature = dss_obj.temperature,
        # p-t_curve = dss_obj.p-tcurve,
        # efficiency_curve = dss_obj.effcurve,
        # r => LinearAlgebra.diagm(0 => fill(dss_obj[Symbol("%r")] / 100., nphases)),
        # x => LinearAlgebra.diagm(0 => fill(dss_obj[Symbol("%x")] / 100., nphases)),
        status = dss_obj.enabled,
        source_id = "pvsystem.$(dss_obj.name)",
        time_series = Dict{String,String}(
            "pg" => getproperty(dss_obj, Symbol(time_series)),
            "qg" => getproperty(dss_obj, Symbol(time_series)),
            "pg_ub" => getproperty(dss_obj, Symbol(time_series)),
            "qg_ub" => getproperty(dss_obj, Symbol(time_series)),
        ),
        dss = import_all ? dss_obj : missing,
    )
end


"""
"""
function create_eng_object(::Type{T}, dss_obj::DssStorage; import_all::Bool=false, time_series::String="daily")::T where T <: EngStorage
    nphases = dss_obj.phases

    T(;
        name = dss_obj.name,
        bus = _parse_bus_id(dss_obj.bus1)[1],
        connections = _get_conductors_ordered(dss_obj.bus1, pad_ground=true, default=[collect(1:nphases)..., 0]),
        configuration = WYE,
        energy = dss_obj.kwhstored,
        energy_ub = dss_obj.kwhrated,
        charge_ub = dss_obj["%charge"] / 100.0 * dss_obj.kwrated,
        discharge_ub = dss_obj["%discharge"] / 100.0 * dss_obj.kwrated,
        sm_ub = dss_obj.kva,
        charge_efficiency = dss_obj["%effcharge"],
        discharge_efficiency = dss_obj["%effdischarge"],
        qs_lb = -dss_obj.kva,
        qs_ub =  dss_obj.kva,
        rs = dss_obj["%r"] / 100.0,
        xs = dss_obj["%x"] / 100.0,
        pex = dss_obj["%idlingkw"] ./ 100.0 .* dss_obj.kwrated,
        qex = dss_obj["%idlingkvar"] ./ 100.0 .* dss_obj.kwrated,
        ps = -dss_obj.kw,
        qs = -dss_obj.kvar,
        status = dss_obj.enabled,
        source_id = "storage.$(dss_obj.name)",
        time_series = Dict{String,String}(
            "ps_ub" => getproperty(dss_obj, Symbol(time_series)),
            "qs_ub" => getproperty(dss_obj, Symbol(time_series)),
        ),
        dss = import_all ? dss_obj : missing,
    )
end
