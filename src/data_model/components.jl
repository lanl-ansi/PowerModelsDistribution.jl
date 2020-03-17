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

    data_eng[obj_type][obj_id] = object
end


""
function create_linecode(; kwargs...)
    linecode = Dict{String,Any}(
        "rs" => get(kwargs, :rs, fill(0.0, n_conductors, n_conductors)),
        "xs" => get(kwargs, :xs, fill(0.0, n_conductors, n_conductors)),
        "g_fr" => get(kwargs, :g_fr, fill(0.0, n_conductors, n_conductors)),
        "b_fr" => get(kwargs, :b_fr, fill(0.0, n_conductors, n_conductors)),
        "g_to" => get(kwargs, :g_to, fill(0.0, n_conductors, n_conductors)),
        "b_to" => get(kwargs, :b_to, fill(0.0, n_conductors, n_conductors)),
)

    n_conductors = 0
    for key in [:rs, :xs, :g_fr, :g_to, :b_fr, :b_to]
        if haskey(kwargs, key)
            n_conductors = size(kwargs[key])[1]
        end
    end

    _add_unused_kwargs!(linecode, kwargs)

    return linecode
end


"creates an engineering model"
function create_eng_model(name; kwargs...)::Dict{String,Any}
    kwargs = Dict{Symbol,Any}(kwargs)

    data_model = Dict{String,Any}(
        "name" => name,
        "data_model" => "engineering",
        "settings" => Dict{String,Any}(
            "v_var_scalar" => get(kwargs, :v_var_scalar, 1e3)
        )
    )

    _add_unused_kwargs!(data_model["settings"], kwargs)

    return data_model
end


""
function create_line(; kwargs...)
    kwargs = Dict{Symbol,Any}(kwargs)

    N = length(get(kwargs, :f_connections, collect(1:4)))

    # if no linecode, then populate loss parameters with zero
    if !haskey(kwargs, :linecode)
        n_conductors = 0
        for key in [:rs, :xs, :g_fr, :g_to, :b_fr, :b_to]
            if haskey(kwargs, key)
                n_conductors = size(kwargs[key])[1]
            end
        end
    end

    line = Dict{String,Any}(
        "status" => get(kwargs, :status, 1),
        "f_connections" => get(kwargs, :f_connections, collect(1:4)),
        "t_connections" => get(kwargs, :t_connections, collect(1:4)),
        "angmin" => get(kwargs, :angmin, fill(-60/180*pi, N)),
        "angmax" => get(kwargs, :angmax, fill( 60/180*pi, N)),
        "rs" => get(kwargs, :rs, fill(0.0, n_conductors, n_conductors)),
        "xs" => get(kwargs, :xs, fill(0.0, n_conductors, n_conductors)),
        "g_fr" => get(kwargs, :g_fr, fill(0.0, n_conductors, n_conductors)),
        "b_fr" => get(kwargs, :b_fr, fill(0.0, n_conductors, n_conductors)),
        "g_to" => get(kwargs, :g_to, fill(0.0, n_conductors, n_conductors)),
        "b_to" => get(kwargs, :b_to, fill(0.0, n_conductors, n_conductors)),
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
    )

    load["connections"] = get(kwargs, :connections, load["configuration"]=="wye" ? [1, 2, 3, 4] : [1, 2, 3]),
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


# create add_{component_type}! methods that will create a component and add it to the engineering data model
for comp in keys(_eng_model_dtypes)
    eval(Meta.parse("add_$(comp)!(data_model, name; kwargs...) = add_object!(data_model, \"$comp\", name, create_$comp(; kwargs...))"))
end
