"adds kwargs that were specified but unused by the required defaults to the component"
function _add_unused_kwargs!(object::Dict{String,<:Any}, kwargs)
    for (property, value) in kwargs
        if !haskey(object, string(property))
            object[string(property)] = value
        end
    end
end


"""
    add_object!(data_eng::Dict{String,<:Any}, obj_type::String, obj_id::String, object::Dict{String,<:Any})

Generic add function to add components to an engineering data model
"""
function add_object!(data_eng::Dict{String,<:Any}, obj_type::String, obj_id::String, object::Dict{String,<:Any})
    if !haskey(data_eng, obj_type)
        data_eng[obj_type] = Dict{String,Any}()
    end

    if !haskey(object, "source_id")
        object["source_id"] = "$obj_type.$obj_id"
    end

    if obj_type == "voltage_source"
        if !haskey(data_eng["settings"], "base_bus")
            data_eng["settings"]["base_bus"] = object["bus"]
        end
    end

    for bus_key in ["f_", "t_", ""]
        if haskey(object, "$(bus_key)bus")
            if !haskey(data_eng, "bus")
                data_eng["bus"] = Dict{String,Any}()
            end

            if obj_type == "transformer"
                if haskey(object, "f_bus") && haskey(object, "t_bus")
                    if !haskey(data_eng["bus"], object["f_bus"])
                        data_eng["bus"][object["f_bus"]] = create_bus(; terminals=object["f_connections"])
                    end

                    if !haskey(data_eng["bus"], object["t_bus"])
                        data_eng["bus"][object["t_bus"]] = create_bus(; terminals=object["t_connections"])
                    end
                else
                    for (wdg, bus_id) in enumerate(object["bus"])
                        if !haskey(data_eng["bus"], bus_id)
                            data_eng["bus"][bus_id] = create_bus(; terminals=object["connections"][wdg])
                        end
                    end
                end
            else
                if !haskey(data_eng["bus"], object["$(bus_key)bus"])
                    data_eng["bus"][object["$(bus_key)bus"]] = create_bus(; terminals=object["$(bus_key)connections"])
                end
            end
        end
    end

    data_eng[obj_type][obj_id] = object
end


"""
    Model(model_type::DataModel)

Instantiates a PowerModelsDistribution data model
"""
function Model(model_type::DataModel=ENGINEERING; kwargs...)::Dict{String,Any}
    if model_type == ENGINEERING
        data_model = Dict{String,Any}(
            "data_model" => model_type,
            "per_unit" => false,
            "settings" => Dict{String,Any}(
                "voltage_scale_factor" => get(kwargs, :voltage_scale_factor, 1e3),
                "power_scale_factor" => get(kwargs, :power_scale_factor, 1e3),
                "vbases_default" => get(kwargs, :vbases_default, Dict{Any,Real}()),
                "sbase_default" => get(kwargs, :sbase_default, 1.0),
                "base_frequency" => get(kwargs, :basefreq, 60.0),
            )
        )

        _add_unused_kwargs!(data_model["settings"], kwargs)
    elseif model_type == MATHEMATICAL
        @warn "There are not currently any helper functions to help build a mathematical model, this will only instantiate required fields."
        data_model = Dict{String,Any}(
            "bus" => Dict{String,Any}(),
            "load" => Dict{String,Any}(),
            "shunt" => Dict{String,Any}(),
            "gen" => Dict{String,Any}(),
            "storage" => Dict{String,Any}(),
            "branch" => Dict{String,Any}(),
            "switch" => Dict{String,Any}(),
            "per_unit" => false,
            "data_model" => model_type
        )

        _add_unused_kwargs!(data_model, kwargs)
    else
        error("Model type '$model_type' not recognized")
    end

    return data_model
end


"""
    create_linecode(
        rs::Matrix{<:Real},
        xs::Matrix{<:Real};
        g_fr::Union{Matrix{<:Real},Missing}=missing,
        b_fr::Union{Matrix{<:Real},Missing}=missing,
        g_to::Union{Matrix{<:Real},Missing}=missing,
        b_to::Union{Matrix{<:Real},Missing}=missing,
        cm_ub::Union{Vector{<:Real},Missing}=missing,
        kwargs...
    )::Dict{String,Any}

creates a linecode with some defaults
"""
function create_linecode(rs::Matrix{<:Real}, xs::Matrix{<:Real};
    g_fr::Union{Matrix{<:Real},Missing}=missing,
    b_fr::Union{Matrix{<:Real},Missing}=missing,
    g_to::Union{Matrix{<:Real},Missing}=missing,
    b_to::Union{Matrix{<:Real},Missing}=missing,
    cm_ub::Union{Vector{<:Real},Missing}=missing,
    kwargs...
        )::Dict{String,Any}

    shape = size(rs)

    for v in [rs, xs, g_fr, g_to, b_fr, b_to]
        if !ismissing(v)
            @assert size(v) == shape "not all of the properties are the same size, aborting linecode creation"
        end
    end

    linecode = Dict{String,Any}(
        "rs" => rs,
        "xs" => xs,
        "g_fr" => !ismissing(g_fr) ? g_fr : fill(0.0, shape...),
        "b_fr" => !ismissing(b_fr) ? b_fr : fill(0.0, shape...),
        "g_to" => !ismissing(g_to) ? g_to : fill(0.0, shape...),
        "b_to" => !ismissing(b_to) ? b_to : fill(0.0, shape...),
    )

    if !ismissing(cm_ub)
        linecode["cm_ub"] = cm_ub
    end

    _add_unused_kwargs!(linecode, kwargs)

    return linecode
end


"""
    create_line(
        f_bus::String,
        t_bus::String,
        f_connections::Vector{Int},
        t_connections::Vector{Int};
        linecode::Union{String,Missing}=missing,
        rs::Union{Matrix{<:Real},Missing}=missing,
        xs::Union{Matrix{<:Real},Missing}=missing,
        g_fr::Union{Matrix{<:Real},Missing}=missing,
        b_fr::Union{Matrix{<:Real},Missing}=missing,
        g_to::Union{Matrix{<:Real},Missing}=missing,
        b_to::Union{Matrix{<:Real},Missing}=missing,
        length::Real=1.0,
        cm_ub::Union{Vector{<:Real},Missing}=missing,
        sm_ub::Union{Vector{<:Real},Missing}=missing,
        vad_lb::Union{Vector{<:Real},Missing}=missing,
        vad_ub::Union{Vector{<:Real},Missing}=missing,
        status::Status=ENABLED,
        kwargs...
    )::Dict{String,Any}

Create a line with some default values
"""
function create_line(f_bus::String, t_bus::String, f_connections::Vector{Int}, t_connections::Vector{Int};
    linecode::Union{String,Missing}=missing,
    rs::Union{Matrix{<:Real},Missing}=missing,
    xs::Union{Matrix{<:Real},Missing}=missing,
    g_fr::Union{Matrix{<:Real},Missing}=missing,
    b_fr::Union{Matrix{<:Real},Missing}=missing,
    g_to::Union{Matrix{<:Real},Missing}=missing,
    b_to::Union{Matrix{<:Real},Missing}=missing,
    length::Real=1.0,
    cm_ub::Union{Vector{<:Real},Missing}=missing,
    sm_ub::Union{Vector{<:Real},Missing}=missing,
    vad_lb::Union{Vector{<:Real},Missing}=missing,
    vad_ub::Union{Vector{<:Real},Missing}=missing,
    status::Status=ENABLED,
    kwargs...
        )::Dict{String,Any}

    n_conductors = size(f_connections)[1]
    shape = (n_conductors, n_conductors)

    for v in [rs, xs, g_fr, b_fr, g_to, b_to, cm_ub, sm_ub, vad_lb, vad_ub]
        if !ismissing(v)
            if isa(v, Matrix)
                @assert size(v) == shape
            else
                @assert size(v)[1] == n_conductors
            end
        end
    end

    line = Dict{String,Any}(
        "f_bus" => f_bus,
        "t_bus" => t_bus,
        "status" => status,
        "f_connections" => f_connections,
        "t_connections" => t_connections,
        "vad_lb" => !ismissing(vad_lb) ? vad_lb : fill(-5.0, n_conductors),
        "vad_ub" => !ismissing(vad_lb) ? vad_lb : fill( 5.0, n_conductors),
        "length" => length,
    )

    if ismissing(linecode)
        if !ismissing(rs) && !ismissing(xs)
            line["rs"] = rs
            line["xs"] = xs

        else
            error("A linecode or rs & xs must be specified to create a valid line object")
        end

        line["g_fr"] = !ismissing(g_fr) ? g_fr : fill(0.0, shape...)
        line["b_fr"] = !ismissing(b_fr) ? b_fr : fill(0.0, shape...)
        line["g_to"] = !ismissing(g_to) ? g_to : fill(0.0, shape...)
        line["b_to"] = !ismissing(b_to) ? b_to : fill(0.0, shape...)
    else
        line["linecode"] = linecode
        for (k,v) in [("rs", rs), ("xs", xs), ("g_fr", g_fr), ("b_fr", b_fr), ("g_to", g_to), ("b_to", b_to)]
            if !ismissing(v)
                line[k] = v
            end
        end
    end

    for (k,v) in [("cm_ub", cm_ub),  ("sm_ub", sm_ub)]
        if !ismissing(v)
            line[k] = v
        end
    end

    _add_unused_kwargs!(line, kwargs)

    return line
end


"""
    create_switch(
        f_bus::String,
        t_bus::String,
        f_connections::Vector{Int},
        t_connections::Vector{Int};
        cm_ub::Union{Vector{<:Real},Missing}=missing,
        sm_ub::Union{Vector{<:Real},Missing}=missing,
        linecode::Union{String,Missing}=missing,
        rs::Union{Matrix{<:Real},Missing}=missing,
        xs::Union{Matrix{<:Real},Missing}=missing,
        dispatchable::Dispatchable=NO,
        state::SwitchState=CLOSED,
        status::Status=ENABLED,
        kwargs...
    )::Dict{String,Any}

creates a switch object with some defaults
"""
function create_switch(f_bus::String, t_bus::String, f_connections::Vector{Int}, t_connections::Vector{Int};
    cm_ub::Union{Vector{<:Real},Missing}=missing,
    sm_ub::Union{Vector{<:Real},Missing}=missing,
    linecode::Union{String,Missing}=missing,
    rs::Union{Matrix{<:Real},Missing}=missing,
    xs::Union{Matrix{<:Real},Missing}=missing,
    dispatchable::Dispatchable=NO,
    state::SwitchState=CLOSED,
    status::Status=ENABLED,
    kwargs...
        )::Dict{String,Any}

    eng_obj = Dict{String,Any}(
        "f_bus" => f_bus,
        "t_bus" => t_bus,
        "f_connections" => f_connections,
        "t_connections" => t_connections,
        "dispatchable" => dispatchable,
        "state" => state,
        "status" => status,
    )

    for (k,v) in [("cm_ub", cm_ub), ("sm_ub", sm_ub), ("linecode", linecode), ("rs", rs), ("xs", xs)]
        if !ismissing(v)
            eng_obj[k] = v
        end
    end

    _add_unused_kwargs!(eng_obj, kwargs)

    return eng_obj
end


"""
    create_bus(;
        status::Status=ENABLED,
        terminals::Vector{Int}=Int[],
        grounded::Vector{Int}=Int[],
        rg::Vector{<:Real}=Float64[],
        xg::Vector{<:Real}=Float64[],
        kwargs...
    )::Dict{String,Any}

creates a bus object with some defaults
"""
function create_bus(;
    status::Status=ENABLED,
    terminals::Vector{Int}=Int[],
    grounded::Vector{Int}=Int[],
    rg::Vector{<:Real}=Float64[],
    xg::Vector{<:Real}=Float64[],
    kwargs...
        )::Dict{String,Any}

    # grounded = Vector{Bool}([terminal in grounded for terminal in terminals])

    bus = Dict{String,Any}(
        "status" => status,
        "terminals" => terminals,
        "grounded" => grounded,
        "rg" => isempty(rg) ? fill(0.0, length(grounded)) : rg,
        "xg" => isempty(xg) ? fill(0.0, length(grounded)) : xg,
    )

    _add_unused_kwargs!(bus, kwargs)

    return bus
end


"""
    create_load(
        bus::String,
        connections::Vector{Int};
        configuration::ConnConfig=WYE,
        model::LoadModel=POWER,
        pd_nom::Union{Vector{<:Real},Missing}=missing,
        qd_nom::Union{Vector{<:Real},Missing}=missing,
        vm_nom::Real=1.0,
        dispatchable::Dispatchable=NO,
        status::Status=ENABLED,
        kwargs...
    )::Dict{String,Any}

creates a load object with some defaults
"""
function create_load(bus::String, connections::Vector{Int};
    configuration::ConnConfig=WYE,
    model::LoadModel=POWER,
    pd_nom::Union{Vector{<:Real},Missing}=missing,
    qd_nom::Union{Vector{<:Real},Missing}=missing,
    vm_nom::Real=1.0,
    dispatchable::Dispatchable=NO,
    status::Status=ENABLED,
    kwargs...
        )::Dict{String,Any}

    n_conductors = configuration == WYE ? length(connections)-1 : length(connections)

    for v in [pd_nom, qd_nom]
        if !ismissing(v)
            @assert length(v) == n_conductors
        end
    end

    load = Dict{String,Any}(
        "bus" => bus,
        "connections" => connections,
        "configuration" => configuration,
        "model" => model,
        "pd_nom" => !ismissing(pd_nom) ? pd_nom : fill(0.0, n_conductors),
        "qd_nom" => !ismissing(qd_nom) ? qd_nom : fill(0.0, n_conductors),
        "vm_nom" => vm_nom,
        "dispatchable" => dispatchable,
        "status" => status,
    )

    _add_unused_kwargs!(load, kwargs)

    return load
end


"""
    create_generator(
        bus::String,
        connections::Vector{Int};
        configuration::ConnConfig=WYE,
        pg::Union{Vector{<:Real},Missing}=missing,
        qg::Union{Vector{<:Real},Missing}=missing,
        vg::Union{Vector{<:Real},Missing}=missing,
        pg_lb::Union{Vector{<:Real},Missing}=missing,
        pg_ub::Union{Vector{<:Real},Missing}=missing,
        qg_lb::Union{Vector{<:Real},Missing}=missing,
        qg_ub::Union{Vector{<:Real},Missing}=missing,
        control_mode::ControlMode=FREQUENCYDROOP,
        status::Status=ENABLED,
        kwargs...
    )::Dict{String,Any}

creates a generator object with some defaults
"""
function create_generator(bus::String, connections::Vector{Int};
    configuration::ConnConfig=WYE,
    pg::Union{Vector{<:Real},Missing}=missing,
    qg::Union{Vector{<:Real},Missing}=missing,
    vg::Union{Vector{<:Real},Missing}=missing,
    pg_lb::Union{Vector{<:Real},Missing}=missing,
    pg_ub::Union{Vector{<:Real},Missing}=missing,
    qg_lb::Union{Vector{<:Real},Missing}=missing,
    qg_ub::Union{Vector{<:Real},Missing}=missing,
    control_mode::ControlMode=FREQUENCYDROOP,
    status::Status=ENABLED,
    kwargs...
        )::Dict{String,Any}

    n_conductors = length(connections)

    generator = Dict{String,Any}(
        "bus" => bus,
        "connections" => connections,
        "configuration" => configuration,
        "control_mode" => control_mode,
        "status" => status,
    )

    for (k,v) in [("pg", pg), ("qg", qg), ("vg", vg), ("pg_lb", pg_lb), ("pg_ub", pg_ub), ("qg_lb", qg_lb), ("qg_ub", qg_ub)]
        if !ismissing(v)
            @assert length(v) == n_conductors
            generator[k] = v
        end
    end

    _add_unused_kwargs!(generator, kwargs)

    return generator
end


"""
    create_xfmrcode(;
        configurations::Union{Vector{ConnConfig},Missing}=missing,
        xsc::Union{Vector{<:Real},Missing}=missing,
        rw::Union{Vector{<:Real},Missing}=missing,
        tm_nom::Union{Vector{<:Real},Missing}=missing,
        tm_lb::Union{Vector{Vector{<:Real}},Missing}=missing,
        tm_ub::Union{Vector{Vector{<:Real}},Missing}=missing,
        tm_set::Union{Vector{Vector{<:Real}},Missing}=missing,
        tm_fix::Union{Vector{Vector{<:Real}},Missing}=missing,
        kwargs...
    )::Dict{String,Any}

creates transformer code with some defaults
"""
function create_xfmrcode(;
    configurations::Union{Vector{ConnConfig},Missing}=missing,
    xsc::Union{Vector{<:Real},Missing}=missing,
    rw::Union{Vector{<:Real},Missing}=missing,
    tm_nom::Union{Vector{<:Real},Missing}=missing,
    tm_lb::Union{Vector{Vector{<:Real}},Missing}=missing,
    tm_ub::Union{Vector{Vector{<:Real}},Missing}=missing,
    tm_set::Union{Vector{Vector{<:Real}},Missing}=missing,
    tm_fix::Union{Vector{Vector{<:Real}},Missing}=missing,
    kwargs...
)::Dict{String,Any}

    n_windings = 0
    for v in [configurations, rw, tm_nom, tm_lb, tm_ub, tm_set, tm_fix]
        if !ismissing(v)
            n_windings = length(v)
            break
        end
    end

    @assert n_windings >= 2 "Cannot determine valid number of windings"
    @assert all(length(v) == n_windings for v in [configurations, rw, tm_nom, tm_lb, tm_ub, tm_set, tm_fix]) "Number of windings inconsistent between parameters"

    n_phases = 0
    for v in [tm_lb, tm_ub, tm_set, tm_fix]
        if !ismissing(v)
            n_windings = length(v[1])
            break
        end
    end

    @assert n_phases >= 1 "Cannot determine valid number of phases"

    eng_obj = Dict{String,Any}(
        "configurations" => !ismissing(configurations) ? configurations : fill(WYE, n_windings),
        "xsc" => !ismissing(xsc) ? xsc : zeros(Int(n_windings * (n_windings-1)//2)),
        "rw" => !ismissing(rw) ? rw : zeros(n_windings),
        "tm_nom" => !ismissing(tm_nom) ? tm_nom : ones(n_windings),
        "tm_set" => !ismissing(tm_set) ? tm_set : fill(fill(1.0, ))
    )

    return eng_obj
end


"""
    create_transformer(
        buses::Vector{String},
        connections::Vector{Vector{Int}};
        configurations::Union{Vector{ConnConfig},Missing}=missing,
        xfmrcode::Union{String,Missing}=missing,
        xsc::Union{Vector{<:Real},Missing}=missing,
        rw::Union{Vector{<:Real},Missing}=missing,
        imag::Real=0.0,
        noloadloss::Real=0.0,
        tm_nom::Union{Vector{<:Real},Missing}=missing,
        tm_lb::Union{Vector{Vector{<:Real}},Missing}=missing,
        tm_ub::Union{Vector{Vector{<:Real}},Missing}=missing,
        tm_set::Union{Vector{Vector{<:Real}},Missing}=missing,
        tm_fix::Union{Vector{Vector{Bool}},Missing}=missing,
        polarity::Union{Vector{Int},Missing}=missing,
        vm_nom::Union{Vector{<:Real},Missing}=missing,
        sm_nom::Union{Vector{<:Real},Missing}=missing,
        status::Status=ENABLED,
        kwargs...
    )::Dict{String,Any}

creates a n-winding transformer object with some defaults
"""
function create_transformer(buses::Vector{String}, connections::Vector{Vector{Int}};
    configurations::Union{Vector{ConnConfig},Missing}=missing,
    xfmrcode::Union{String,Missing}=missing,
    xsc::Union{Vector{<:Real},Missing}=missing,
    rw::Union{Vector{<:Real},Missing}=missing,
    imag::Real=0.0,
    noloadloss::Real=0.0,
    tm_nom::Union{Vector{<:Real},Missing}=missing,
    tm_lb::Union{Vector{Vector{<:Real}},Missing}=missing,
    tm_ub::Union{Vector{Vector{<:Real}},Missing}=missing,
    tm_set::Union{Vector{Vector{<:Real}},Missing}=missing,
    tm_fix::Union{Vector{Vector{Bool}},Missing}=missing,
    polarity::Union{Vector{Int},Missing}=missing,
    vm_nom::Union{Vector{<:Real},Missing}=missing,
    sm_nom::Union{Vector{<:Real},Missing}=missing,
    status::Status=ENABLED,
    kwargs...
        )::Dict{String,Any}

    n_windings = length(buses)
    n_conductors = length(connections[1])

    transformer = Dict{String,Any}(
        "buses" => buses,
        "connections" => connections,
        "configurations" => !ismissing(configurations) ? configurations : fill(WYE, n_windings),
        "xsc" => !ismissing(xsc) ? xsc : zeros(Int(n_windings * (n_windings-1)//2)),
        "rw" => !ismissing(rw) ? rw : zeros(n_windings),
        "cmag" => imag,
        "noloadloss" => noloadloss,
        "tm_nom" => !ismissing(tm_nom) ? tm_nom : ones(n_windings),
        "tm_set" => !ismissing(tm_set) ? tm_set : fill(fill(1.0, n_conductors), n_windings),
        "tm_fix" => !ismissing(tm_fix) ? tm_fix : fill(fill(true, n_conductors), n_windings),
        "polarity" => !ismissing(polarity) ? polarity : fill(1, n_windings),
        "status" => status,
    )

    for (k,v) in [("tm_lb", tm_lb), ("tm_ub", tm_ub), ("vm_nom", vm_nom), ("sm_nom", sm_nom)]
        if !ismissing(v)
            transformer[k] = v
        end
    end

    _add_unused_kwargs!(transformer, kwargs)

    return transformer
end


"""
    create_al2w_transformer(
        f_bus::String,
        t_bus::String,
        f_connections::Vector{Int},
        t_connections::Vector{Int};
        configuration::ConnConfig=WYE,
        tm_nom::Real=1.0,
        tm_lb::Union{Vector{<:Real},Missing}=missing,
        tm_ub::Union{Vector{<:Real},Missing}=missing,
        tm_set::Union{Vector{<:Real},Missing}=missing,
        tm_fix::Union{Vector{Bool},Missing}=missing,
        status::Status=ENABLED,
        kwargs...
    )::Dict{String,Any}

creates a aysmmetric lossless 2-winding transformer object with some defaults
"""
function create_al2w_transformer(f_bus::String, t_bus::String, f_connections::Vector{Int}, t_connections::Vector{Int};
    configuration::ConnConfig=WYE,
    tm_nom::Real=1.0,
    tm_lb::Union{Vector{<:Real},Missing}=missing,
    tm_ub::Union{Vector{<:Real},Missing}=missing,
    tm_set::Union{Vector{<:Real},Missing}=missing,
    tm_fix::Union{Vector{Bool},Missing}=missing,
    status::Status=ENABLED,
    kwargs...
        )::Dict{String,Any}

    n_conductors = length(f_connections)
    n_conductors = configuration == WYE ? n_conductors-1 : n_conductors

    transformer = Dict{String,Any}(
        "f_bus" => f_bus,
        "t_bus" => t_bus,
        "f_connections" => f_connections,
        "t_connections" => t_connections,
        "configuration" => configuration,
        "tm_nom" => tm_nom,
        "tm_set" => !ismissing(tm_set) ? tm_set : fill(1.0, n_conductors),
        "tm_fix" => !ismissing(tm_fix) ? tm_fix : fill(true, n_conductors),
        "status" => status,
    )

    for (k,v) in [("tm_lb", tm_lb), ("tm_ub", tm_ub)]
        if !ismissing(v)
            transformer[k] = v
        end
    end

    _add_unused_kwargs!(transformer, kwargs)

    return transformer
end


"""
    create_shunt(
        bus::String,
        connections::Vector{Int};
        gs::Union{Matrix{<:Real},Missing}=missing,
        bs::Union{Matrix{<:Real},Missing}=missing,
        model::ShuntModel=GENERIC,
        dispatchable::Dispatchable=NO,
        status::Status=ENABLED,
        kwargs...
    )::Dict{String,Any}

creates a generic shunt with some defaults
"""
function create_shunt(bus::String, connections::Vector{Int};
    gs::Union{Matrix{<:Real},Missing}=missing,
    bs::Union{Matrix{<:Real},Missing}=missing,
    model::ShuntModel=GENERIC,
    dispatchable::Dispatchable=NO,
    status::Status=ENABLED,
    kwargs...
        )::Dict{String,Any}

    n_conductors = length(connections)

    shunt = Dict{String,Any}(
        "bus" => bus,
        "connections" => connections,
        "gs" => !ismissing(gs) ? gs : fill(0.0, n_conductors, n_conductors),
        "bs" => !ismissing(bs) ? bs : fill(0.0, n_conductors, n_conductors),
        "model" => model,
        "dispatchable" => dispatchable,
        "status" => status,
    )

    _add_unused_kwargs!(shunt, kwargs)

    return shunt
end


"""
    create_solar(
        bus::String,
        connections::Vector{Int};
        configuration::ConnConfig=WYE,
        pg_lb::Union{Vector{<:Real},Missing}=missing,
        pg_ub::Union{Vector{<:Real},Missing}=missing,
        qg_lb::Union{Vector{<:Real},Missing}=missing,
        qg_ub::Union{Vector{<:Real},Missing}=missing,
        pg::Union{Vector{<:Real},Missing}=missing,
        qg::Union{Vector{<:Real},Missing}=missing,
        status::Status=ENABLED,
        kwargs...
    )::Dict{String,Any}

creates a solar generator with some defaults
"""
function create_solar(bus::String, connections::Vector{Int};
    configuration::ConnConfig=WYE,
    pg_lb::Union{Vector{<:Real},Missing}=missing,
    pg_ub::Union{Vector{<:Real},Missing}=missing,
    qg_lb::Union{Vector{<:Real},Missing}=missing,
    qg_ub::Union{Vector{<:Real},Missing}=missing,
    pg::Union{Vector{<:Real},Missing}=missing,
    qg::Union{Vector{<:Real},Missing}=missing,
    status::Status=ENABLED,
    kwargs...
        )::Dict{String,Any}

    eng_obj = Dict{String,Any}(
        "bus" => bus,
        "connections" => connections,
        "configuration" => configuration,
        "status" => status
    )

    for (k,v) in [("pg_lb", pg_lb), ("pg_ub", pg_ub), ("qg_lb", qg_lb), ("qg_ub", qg_ub), ("pg", pg), ("qg", qg)]
        if !ismissing(v)
            eng_obj[k] = v
        end
    end

    _add_unused_kwargs!(eng_obj, kwargs)

    return eng_obj
end


"""
    create_storage(
        configuration::ConnConfig=WYE,
        energy::Real=0.0,
        energy_ub::Real=0.0,
        charge_ub::Real=0.0,
        discharge_ub::Real=0.0,
        sm_ub::Union{Real,Missing}=missing,
        cm_ub::Union{Real,Missing}=missing,
        charge_efficiency::Real=1.0,
        discharge_efficiency::Real=1.0,
        qs_lb::Union{Real,Missing}=missing,
        qs_ub::Union{Real,Missing}=missing,
        rs::Real=0.0,
        xs::Real=0.0,
        pex::Real=0.0,
        qex::Real=0.0,
        ps::Union{Real,Vector{<:Real},Missing}=missing,
        qs::Union{Real,Vector{<:Real},Missing}=missing,
        status::Status=ENABLED,
        kwargs...
        )::Dict{String,Any}

creates energy storage object with some defaults
"""
function create_storage(bus::String, connections::Vector{Int};
    configuration::ConnConfig=WYE,
    energy::Real=0.0,
    energy_ub::Real=0.0,
    charge_ub::Real=0.0,
    discharge_ub::Real=0.0,
    sm_ub::Union{Real,Missing}=missing,
    cm_ub::Union{Real,Missing}=missing,
    charge_efficiency::Real=1.0,
    discharge_efficiency::Real=1.0,
    qs_lb::Union{Real,Missing}=missing,
    qs_ub::Union{Real,Missing}=missing,
    rs::Real=0.0,
    xs::Real=0.0,
    pex::Real=0.0,
    qex::Real=0.0,
    ps::Union{Real,Vector{<:Real},Missing}=missing,
    qs::Union{Real,Vector{<:Real},Missing}=missing,
    status::Status=ENABLED,
    kwargs...
        )::Dict{String,Any}

    n_conductors = length(connections)

    eng_obj = Dict{String,Any}(
        "bus" => bus,
        "connections" => connections,
        "configuration" => configuration,
        "energy" => energy,
        "energy_ub" => energy_ub,
        "charge_ub" => charge_ub,
        "discharge_ub" => discharge_ub,
        "charge_efficiency" => charge_efficiency,
        "discharge_efficiency" => discharge_efficiency,
        "rs" => rs,
        "xs" => xs,
        "pex" => pex,
        "qex" => qex,
        "status" => status
    )

    for (k,v) in [("sm_ub", sm_ub), ("cm_ub", cm_ub), ("qs_lb", qs_lb), ("qs_ub", qs_ub), ("ps", ps), ("qs", qs)]
        if !ismissing(v)
            eng_obj[k] = v
        end
    end

    _add_unused_kwargs!(eng_obj, kwargs)

    return eng_obj
end


"""
    create_voltage_source(
        bus::String,
        connections::Vector{Int};
        configuration::ConnConfig=WYE,
        vm::Union{Vector{<:Real},Missing}=missing,
        va::Union{Vector{<:Real},Missing}=missing,
        vm_lb::Union{Vector{<:Real},Missing}=missing,
        vm_ub::Union{Vector{<:Real},Missing}=missing,
        rs::Union{Vector{<:Real},Missing}=missing,
        xs::Union{Vector{<:Real},Missing}=missing,
        status::Status=ENABLED,
        kwargs...
    )::Dict{String,Any}

creates a voltage source with some defaults
"""
function create_voltage_source(bus::String, connections::Vector{Int};
    configuration::ConnConfig=WYE,
    vm::Union{Vector{<:Real},Missing}=missing,
    va::Union{Vector{<:Real},Missing}=missing,
    vm_lb::Union{Vector{<:Real},Missing}=missing,
    vm_ub::Union{Vector{<:Real},Missing}=missing,
    rs::Union{Vector{<:Real},Missing}=missing,
    xs::Union{Vector{<:Real},Missing}=missing,
    status::Status=ENABLED,
    kwargs...
        )::Dict{String,Any}

    n_conductors = length(connections)
    nphases = configuration == WYE ? n_conductors - 1 : n_conductors

    voltage_source = Dict{String,Any}(
        "bus" => bus,
        "connections" => connections,
        "configuration" => configuration,
        "vm" => !ismissing(vm) ? vm : ones(n_conductors),
        "va" => !ismissing(va) ? va : [[0., -120., 120.]..., zeros(n_conductors)...][connections],
        "status" => status,
    )

    for (k,v) in [("rs", rs), ("xs", xs), ("vm_lb", vm_lb), ("vm_ub", vm_ub)]
        if !ismissing(v)
            if isa(v, Vector)
                @assert length(v) == n_conductors "$k is the wrong length, expected $n_conductors but got $(length(v))"
            elseif isa(v, Matrix)
                @assert size(v) == (n_conductors, n_conductors) "$k is the wrong size, expected ($n_conductors, $n_conductors) but got $(size(v))"
            end
            voltage_source[k] = v
        end
    end

    _add_unused_kwargs!(voltage_source, kwargs)

    return voltage_source
end


"""
    delete_component!(data_eng::Dict{String,<:Any}, component_type::String, component_id::String)

deletes a component from the engineering data model
"""
function delete_component!(data_eng::Dict{String,<:Any}, component_type::String, component_id::String)
    delete!(data_eng[component_type], component_id)
    if isempty(data_eng[component_type])
        delete!(data_eng, component_type)
    end
end


"""
    add_vbase_default!(data_eng::Dict{String,<:Any}, bus::String, vbase::Real)

Function to add default vbase for a bus
"""
function add_vbase_default!(data_eng::Dict{String,<:Any}, bus::String, vbase::Real)
    if !haskey(data_eng, "settings")
        data_eng["settings"] = Dict{String,Any}()
    end

    if !haskey(data_eng["settings"], "vbases_default")
        data_eng["settings"]["vbases_default"] = Dict{Any,Real}()
    end

    data_eng["settings"]["vbases_default"][bus] = vbase
end


# Data objects
add_bus!(data_eng::Dict{String,<:Any}, id::String; kwargs...) = add_object!(data_eng, "bus", id, create_bus(; kwargs...))
add_linecode!(data_eng::Dict{String,<:Any}, id::String, rs::Matrix{<:Real}, xs::Matrix{<:Real}; kwargs...) = add_object!(data_eng, "linecode", id, create_linecode(rs, xs; kwargs...))
add_xfmrcode!(data_eng::Dict{String,<:Any}, id::String; kwargs...) = add_object!(data_eng, "xfmrcode", id, create_xfmrcode(; kwargs...))
# add_time_series!(data_eng::Dict{String,<:Any}, id::String; kwargs...) = add_object!(data_eng, "time_series", id, create_timeseries(; kwargs...))

@doc "adds a bus to provided ENGINEERING model, see [`create_bus`](@ref create_bus)" add_bus!
@doc "adds a linecode to provided ENGINEERING model, see [`create_linecode`](@ref create_linecode)" add_linecode!
@doc "adds a transformer code (xmfrcode) to provided ENGINEERING model, see [`create_xfmrcode`](@ref create_xfmrcode)" add_xfmrcode!
# @doc "adds a time series to provided ENGINEERING model, see [`create_timeseries`](@ref create_timeseries)" add_time_series!

# Edge objects
add_line!(data_eng::Dict{String,<:Any}, id::String, f_bus::String, t_bus::String, f_connections::Vector{Int}, t_connections::Vector{Int}; kwargs...) = add_object!(data_eng, "line", id, create_line(f_bus, t_bus, f_connections, t_connections; kwargs...))
add_transformer!(data_eng::Dict{String,<:Any}, id::String, buses::Vector{<:String}, connections::Vector{Vector{Int}}; kwargs...) = add_object!(data_eng, "transformer", id, create_transformer(buses, connections; kwargs...))
add_transformer!(data_eng::Dict{String,<:Any}, id::String, f_bus::String, t_bus::String, f_connections::Vector{Int}, t_connections::Vector{Int}; kwargs...) = add_object!(data_eng, "transformer", id, create_al2w_transformer(f_bus, t_bus, f_connections, t_connections; kwargs...))
add_switch!(data_eng::Dict{String,<:Any}, id::String, f_bus::String, t_bus::String, f_connections::Vector{Int}, t_connections::Vector{Int}; kwargs...) = add_object!(data_eng, "switch", id, create_switch(f_bus, t_bus, f_connections, t_connections; kwargs...))

@doc "adds a line to provided ENGINEERING model, see [`create_line`](@ref create_line)" add_line!
@doc "adds a transformer to provided ENGINEERING model, see [`create_transformer`](@ref create_transformer) and [`create_al2w_transformer`](@ref create_al2w_transformer)" add_transformer!
@doc "adds a switch to provided ENGINEERING model, see [`create_switch`](@ref create_switch)" add_switch!

# Node objects
add_load!(data_eng::Dict{String,<:Any}, id::String, bus::String, connections::Vector{Int}; kwargs...) = add_object!(data_eng, "load", id, create_load(bus, connections; kwargs...))
add_shunt!(data_eng::Dict{String,<:Any}, id::String, bus::String, connections::Vector{Int}; kwargs...) = add_object!(data_eng, "shunt", id, create_shunt(bus, connections; kwargs...))
add_voltage_source!(data_eng::Dict{String,<:Any}, id::String, bus::String, connections::Vector{Int}; kwargs...) = add_object!(data_eng, "voltage_source", id, create_voltage_source(bus, connections; kwargs...))
add_generator!(data_eng::Dict{String,<:Any}, id::String, bus::String, connections::Vector{Int}; kwargs...) = add_object!(data_eng, "generator", id, create_generator(bus, connections; kwargs...))
add_storage!(data_eng::Dict{String,<:Any}, id::String, bus::String, connections::Vector{Int}; kwargs...) = add_object!(data_eng, "storage", id, create_storage(bus, connections; kwargs...))
add_solar!(data_eng::Dict{String,<:Any}, id::String, bus::String, connections::Vector{Int}; kwargs...) = add_object!(data_eng, "solar", id, create_solar(bus, connections; kwargs...))

@doc "adds a load to provided ENGINEERING model, see [`create_load`](@ref create_load)" add_load!
@doc "adds a shunt to provided ENGINEERING model, see [`create_shunt`](@ref create_shunt)" add_shunt!
@doc "adds a voltage source to provided ENGINEERING model, see [`create_voltage_source`](@ref create_voltage_source)" add_voltage_source!
@doc "adds a generator to provided ENGINEERING model, see [`create_generator`](@ref create_generator)" add_generator!
@doc "adds a storage to provided ENGINEERING model, see [`create_storage`](@ref create_storage)" add_storage!
@doc "adds a PV to provided ENGINEERING model, see [`create_solar`](@ref create_solar)" add_solar!
