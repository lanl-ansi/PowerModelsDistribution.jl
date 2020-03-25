#TODO
# Can buses in a voltage zone have different terminals?
# Add current/power bounds to data model


"adds kwargs that were specified but unused by the required defaults to the component"
function _add_unused_kwargs!(object::Dict{String,<:Any}, kwargs::Dict{Symbol,<:Any})
    for (property, value) in kwargs
        if !haskey(object, string(property))
            object[string(property)] = value
        end
    end
end


"Generic add function to add components to an engineering data model"
function add_object!(data_eng::Dict{String,<:Any}, obj_type::String, obj_id::Any, object::Dict{String,<:Any})
    if !haskey(data_eng, obj_type)
        data_eng[obj_type] = Dict{Any,Any}()
    end

    if !haskey(object, "source_id")
        object["source_id"] = "$obj_type.$obj_id"
    end

    if obj_type == "voltage_source"
        if !haskey(data_eng["settings"], "set_vbase_bus")
            data_eng["settings"]["set_vbase_bus"] = object["bus"]
        end
    end

    for bus_key in ["f_", "t_", ""]
        if haskey(object, "$(bus_key)bus")
            if !haskey(data_eng, "bus")
                data_eng["bus"] = Dict{Any,Any}()
            end

            if obj_type == "transformer"
                for (wdg, bus_id) in enumerate(data_eng["bus"])
                    if !haskey(data_eng["bus"], bus_id)
                        data_eng["bus"][bus_id] = create_bus(; terminals=object["connections"][wdg])
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


"Instantiates a PowerModelsDistribution data model"
function Model(model_type::String="engineering"; kwargs...)::Dict{String,Any}
    kwargs = Dict{Symbol,Any}(kwargs)

    if model_type == "engineering"
        data_model = Dict{String,Any}(
            "data_model" => "engineering",
            "per_unit" => false,
            "settings" => Dict{String,Any}(
                "v_var_scalar" => get(kwargs, :v_var_scalar, 1e3),
                "set_vbase_val" => get(kwargs, :basekv, 1.0),
                "set_sbase_val" => get(kwargs, :baseMVA, 1.0),
                "basefreq" => get(kwargs, :basefreq, 60.0),
            )
        )

        _add_unused_kwargs!(data_model["settings"], kwargs)
    elseif model_type == "mathematical"
        Memento.warn(_LOGGER, "There are not currently any helper functions to help build a mathematical model, this will only instantiate required fields.")
        data_model = Dict{String,Any}(
            "bus" => Dict{String,Any}(),
            "load" => Dict{String,Any}(),
            "shunt" => Dict{String,Any}(),
            "gen" => Dict{String,Any}(),
            "storage" => Dict{String,Any}(),
            "branch" => Dict{String,Any}(),
            "switch" => Dict{String,Any}(),
            "dcline" => Dict{String,Any}(),
            "per_unit" => false,
            "baseMVA" => 100.0,
            "basekv" => 1.0,
            "data_model" => "mathematical"
        )

        _add_unused_kwargs!(data_model, kwargs)
    else
        Memento.error(_LOGGER, "Model type '$model_type' not recognized")
    end

    return data_model
end


""
function create_linecode(; kwargs...)
    n_conductors = 0
    for key in [:rs, :xs, :g_fr, :g_to, :b_fr, :b_to]
        if haskey(kwargs, key)
            n_conductors = size(kwargs[key])[1]
        end
    end

    linecode = Dict{String,Any}(
        "rs" => get(kwargs, :rs, fill(0.0, n_conductors, n_conductors)),
        "xs" => get(kwargs, :xs, fill(0.0, n_conductors, n_conductors)),
        "g_fr" => get(kwargs, :g_fr, fill(0.0, n_conductors, n_conductors)),
        "b_fr" => get(kwargs, :b_fr, fill(0.0, n_conductors, n_conductors)),
        "g_to" => get(kwargs, :g_to, fill(0.0, n_conductors, n_conductors)),
        "b_to" => get(kwargs, :b_to, fill(0.0, n_conductors, n_conductors)),
    )

    _add_unused_kwargs!(linecode, kwargs)

    return linecode
end


""
function create_line(; kwargs...)
    kwargs = Dict{Symbol,Any}(kwargs)

    @assert haskey(kwargs, :f_bus) && haskey(kwargs, :t_bus) "Line must at least have f_bus and t_bus specified"

    N = length(get(kwargs, :f_connections, collect(1:4)))

    # if no linecode, then populate loss parameters with zero
    if !haskey(kwargs, :linecode)
        n_conductors = N
        for key in [:rs, :xs, :g_fr, :g_to, :b_fr, :b_to]
            if haskey(kwargs, key)
                n_conductors = size(kwargs[key])[1]
            end
        end
    end

    line = Dict{String,Any}(
        "f_bus" => kwargs[:f_bus],
        "t_bus" => kwargs[:t_bus],
        "status" => get(kwargs, :status, 1),
        "f_connections" => get(kwargs, :f_connections, collect(1:4)),
        "t_connections" => get(kwargs, :t_connections, collect(1:4)),
        "angmin" => get(kwargs, :angmin, fill(-60/180*pi, N)),
        "angmax" => get(kwargs, :angmax, fill( 60/180*pi, N)),
        "length" => get(kwargs, :length, 1.0),
        "rs" => get(kwargs, :rs, diagm(0 => fill(0.01, n_conductors))),
        "xs" => get(kwargs, :xs, diagm(0 => fill(0.01, n_conductors))),
        "g_fr" => get(kwargs, :g_fr, diagm(0 => fill(0.0, n_conductors))),
        "b_fr" => get(kwargs, :b_fr, diagm(0 => fill(0.0, n_conductors))),
        "g_to" => get(kwargs, :g_to, diagm(0 => fill(0.0, n_conductors))),
        "b_to" => get(kwargs, :b_to, diagm(0 => fill(0.0, n_conductors))),
    )

    _add_unused_kwargs!(line, kwargs)

    return line
end


"creates a bus object with some defaults"
function create_bus(; kwargs...)
    kwargs = Dict{Symbol,Any}(kwargs)

    bus = Dict{String,Any}(
        "status" => get(kwargs, :status, 1),
        "terminals" => get(kwargs, :terminals, collect(1:4)),
        "grounded" => get(kwargs, :grounded, []),
        "bus_type" => get(kwargs, :bus_type, 1),
        "rg" => get(kwargs, :rg, Array{Float64, 1}()),
        "xg" => get(kwargs, :xg, Array{Float64, 1}()),
    )

    _add_unused_kwargs!(bus, kwargs)

    return bus
end


"creates a load object with some defaults"
function create_load(; kwargs...)
    kwargs = Dict{Symbol,Any}(kwargs)

    load = Dict{String,Any}(
        "status" => get(kwargs, :status, 1),
        "configuration" => get(kwargs, :configuration, "wye"),
        "model" => get(kwargs, :model, "constant_power"),
        "connections" => get(kwargs, :connections, get(kwargs, :configuration, "wye")=="wye" ? [1, 2, 3, 4] : [1, 2, 3]),
        "vnom" => get(kwargs, :vnom, 1.0)
    )

    if load["model"]=="constant_power"
        load["pd"] = get(kwargs, :pd, fill(0.0, 3))
        load["qd"] = get(kwargs, :qd, fill(0.0, 3))
    else
        load["pd_ref"] = get(kwargs, :pd_ref, fill(0.0, 3))
        load["qd_ref"] = get(kwargs, :qd_ref, fill(0.0, 3))
    end

    _add_unused_kwargs!(load, kwargs)

    return load
end


"creates a generator object with some defaults"
function create_generator(; kwargs...)
    kwargs = Dict{Symbol,Any}(kwargs)

    generator = Dict{String,Any}(
        "status" => get(kwargs, :status, 1),
        "configuration" => get(kwargs, :configuration, "wye"),
        "cost" => get(kwargs, :cost, [1.0, 0.0]*1E-3),
    )

    generator["connections"] = get(kwargs, :connections, generator["configuration"]=="wye" ? [1, 2, 3, 4] : [1, 2, 3])

    _add_unused_kwargs!(generator, kwargs)

    return generator
end


"creates a n-winding transformer object with some defaults"
function create_transformer(; kwargs...)
    kwargs = Dict{Symbol,Any}(kwargs)

    @assert haskey(kwargs, :bus) "bus must be defined at the very least"
    n_windings = length(kwargs[:bus])

    transformer = Dict{String,Any}(
        "status" => get(kwargs, :status, 1),
        "configuration" => get(kwargs, :configuration, fill("wye", n_windings)),
        "polarity" => get(kwargs, :polarity, fill(true, n_windings)),
        "rs" => get(kwargs, :rs, zeros(n_windings)),
        "xsc" => get(kwargs, :xsc, zeros(n_windings^2-n_windings)),
        "noloadloss" => get(kwargs, :noloadloss, 0.0),
        "imag" => get(kwargs, :imag, 0.0),
        "tm" => get(kwargs, :tm, fill(fill(1.0, 3), n_windings)),
        "tm_min" => get(kwargs, :tm_min, fill(fill(0.9, 3), n_windings)),
        "tm_max" => get(kwargs, :tm_max, fill(fill(1.1, 3), n_windings)),
        "tm_step" => get(kwargs, :tm_step, fill(fill(1/32, 3), n_windings)),
        "tm_fix" => get(kwargs, :tm_fix, fill(fill(true, 3), n_windings)),
        "connections" => get(kwargs, :connections, fill(collect(1:4), n_windings)),
    )

    _add_unused_kwargs!(transformer, kwargs)

    return transformer
end


"creates a shunt capacitor object with some defaults"
function create_shunt_capacitor(; kwargs...)
    shunt_capacitor = Dict{String,Any}(
        "status" => get(kwargs, :status, 1),
        "configuration" => get(kwargs, :configuration, "wye"),
        "connections" => get(kwargs, :connections, collect(1:4)),
        "qd_ref" => get(kwargs, :qd_ref, fill(0.0, 3)),
    )

    _add_unused_kwargs!(shunt_capacitor, kwargs)

    return shunt_capacitor
end


"creates a generic shunt with some defaults"
function create_shunt(; kwargs...)
    kwargs = Dict{Symbol,Any}(kwargs)

    N = length(get(kwargs, :connections, collect(1:4)))

    shunt = Dict{String,Any}(
        "status" => get(kwargs, :status, 1),
        "g_sh" => get(kwargs, :g_sh, fill(0.0, N, N)),
        "b_sh" => get(kwargs, :b_sh, fill(0.0, N, N)),
    )

    _add_unused_kwargs!(shunt, kwargs)

    return shunt
end


"creates a voltage source with some defaults"
function create_voltage_source(; kwargs...)
    kwargs = Dict{Symbol,Any}(kwargs)

    voltage_source = Dict{String,Any}(
        "status" => get(kwargs, :status, 1),
        "connections" => get(kwargs, :connections, collect(1:3)),
    )

    nphases = length(voltage_source["connections"])
    voltage_source["vm"] = get(kwargs, :vm, fill(1.0, nphases))
    voltage_source["va"] = deg2rad.(get(kwargs, :va, rad2deg.(_wrap_to_pi.([-2*pi/nphases*(i-1) for i in 1:nphases]))))
    voltage_source["rs"] = fill(0.1, nphases, nphases)
    voltage_source["xs"] = fill(0.0, nphases, nphases)

    _add_unused_kwargs!(voltage_source, kwargs)

    return voltage_source
end


"deletes a component from the engineering data model"
function delete_component!(data_eng::Dict{String,<:Any}, component_type::String, component_id::Any)
    delete!(data_eng[component_type], component_id)
    if isempty(data_eng[component_type])
        delete!(data_eng, component_type)
    end
end

# Data objects
add_bus!(data_eng::Dict{String,<:Any}, id::Any; kwargs...) = add_object!(data_eng, "bus", id, create_bus(; kwargs...))
add_linecode!(data_eng::Dict{String,<:Any}, id::Any; kwargs...) = add_object!(data_eng, "linecode", id, create_linecode(; kwargs...))
# add_xfmrcode!(data_eng::Dict{String,<:Any}, id::Any; kwargs...) = add_object!(data_eng, "xfmrcode", id, create_xfmrcode(; kwargs...))
# add_timeseries!(data_eng::Dict{String,<:Any}, id::Any; kwargs...) = add_object!(data_eng, "timeseries", id, create_timeseries(; kwargs...))

# Edge objects
add_line!(data_eng::Dict{String,<:Any}, id::Any, f_bus::Any, t_bus::Any; kwargs...) = add_object!(data_eng, "line", id, create_line(; f_bus=f_bus, t_bus=t_bus, kwargs...))
add_transformer!(data_eng::Dict{String,<:Any}, id::Any, bus::Vector{Any}; kwargs...) = add_object!(data_eng, "transformer", id, create_transformer(; bus=bus, kwargs...))
# add_switch!(data_eng::Dict{String,<:Any}, id::Any, f_bus::Any, t_bus::Any; kwargs...) = add_object!(data_eng, "switch", id, create_switch(; f_bus=f_bus, t_bus=t_bus, kwargs...))
# add_series_capacitor!(data_eng::Dict{String,<:Any}, id::Any, f_bus::Any, t_bus::Any; kwargs...) = add_object!(data_eng, "series_capacitor", id, create_series_capacitor(; f_bus=f_bus, t_bus=t_bus, kwargs...))
# add_line_reactor!(data_eng::Dict{String,<:Any}, id::Any, f_bus::Any, t_bus::Any; kwargs...) = add_object!(data_eng, "line_reactor", id, create_line_reactor(; f_bus=f_bus, t_bus=t_bus, kwargs...))

# Node objects
add_load!(data_eng::Dict{String,<:Any}, id::Any, bus::Any; kwargs...) = add_object!(data_eng, "load", id, create_load(; bus=bus, kwargs...))
add_shunt_capacitor!(data_eng::Dict{String,<:Any}, id::Any, bus::Any; kwargs...) = add_object!(data_eng, "shunt_capacitor", id, create_shunt_capacitor(; bus=bus, kwargs...))
# add_shunt_reactor!(data_eng::Dict{String,<:Any}, id::Any, bus::Any; kwargs...) = add_object!(data_eng, "shunt_reactor", id, create_shunt_reactor(; bus=bus, kwargs...))
add_shunt!(data_eng::Dict{String,<:Any}, id::Any, bus::Any; kwargs...) = add_object!(data_eng, "shunt", id, create_shunt(; bus=bus, kwargs...))
add_voltage_source!(data_eng::Dict{String,<:Any}, id::Any, bus::Any; kwargs...) = add_object!(data_eng, "voltage_source", id, create_voltage_source(; bus=bus, kwargs...))
add_generator!(data_eng::Dict{String,<:Any}, id::Any, bus::Any; kwargs...) = add_object!(data_eng, "generator", id, create_generator(; bus=bus, kwargs...))
# add_storage!(data_eng::Dict{String,<:Any}, id::Any, bus::Any; kwargs...) = add_object!(data_eng, "storage", id, create_storage(; bus=bus, kwargs...))
# add_solar!(data_eng::Dict{String,<:Any}, id::Any, bus::Any; kwargs...) = add_object!(data_eng, "solar", id, create_solar(; bus=bus, kwargs...))
# add_wind!(data_eng::Dict{String,<:Any}, id::Any, bus::Any; kwargs...) = add_object!(data_eng, "wind", id, create_wind(; bus=bus, kwargs...))