#TODO
# Can buses in a voltage zone have different terminals?
# Add current/power bounds to data model


function copy_kwargs_to_dict_if_present!(dict, kwargs, args)
    for arg in args
        if haskey(kwargs, arg)
            dict[string(arg)] = kwargs[arg]
        end
    end
end

function add_kwarg!(dict, kwargs, name, default)
    if haskey(kwargs, name)
        dict[string(name)] = kwargs[name]
    else
        dict[string(name)] = default
    end
end

function component_dict_from_list!(list)
    dict = Dict{String, Any}()
    for object in list
        dict[object["obj_name"]] = object
    end
    return dict
end

const _eng_model_dtypes = Dict{Symbol,Dict{Symbol,Type}}(
    :linecode => Dict{Symbol,Type}(
        :obj_name => Any,
        :rs => Array{<:Real, 2},
        :xs => Array{<:Real, 2},
        :g_fr => Array{<:Real, 2},
        :g_to => Array{<:Real, 2},
        :b_fr => Array{<:Real, 2},
        :b_to => Array{<:Real, 2}
    ),
    :line => Dict{Symbol,Type}(
        :obj_name => Any,
        :status => Int,
        :f_bus => AbstractString,
        :t_bus => AbstractString,
        :f_connections => Vector{<:Int},
        :t_connections => Vector{<:Int},
        :linecode => AbstractString,
        :length => Real,
        :c_rating =>Vector{<:Real},
        :s_rating =>Vector{<:Real},
        :angmin=>Vector{<:Real},
        :angmax=>Vector{<:Real},
        :rs => Array{<:Real, 2},
        :xs => Array{<:Real, 2},
        :g_fr => Array{<:Real, 2},
        :g_to => Array{<:Real, 2},
        :b_fr => Array{<:Real, 2},
        :b_to => Array{<:Real, 2},
    ),
    :bus => Dict{Symbol,Type}(
        :obj_name => Any,
        :status => Int,
        :bus_type => Int,
        :terminals => Array{<:Any},
        :phases => Array{<:Int},
        :neutral => Union{Int, Missing},
        :grounded => Array{<:Any},
        :rg => Array{<:Real},
        :xg => Array{<:Real},
        :vm_pn_min => Real,
        :vm_pn_max => Real,
        :vm_pp_min => Real,
        :vm_pp_max => Real,
        :vm_min => Array{<:Real, 1},
        :vm_max => Array{<:Real, 1},
        :vm_fix => Array{<:Real, 1},
        :va_fix => Array{<:Real, 1},
    ),
    :load => Dict{Symbol,Type}(
        :obj_name => Any,
        :status => Int,
        :bus => Any,
        :connections => Vector,
        :configuration => String,
        :model => String,
        :pd => Array{<:Real, 1},
        :qd => Array{<:Real, 1},
        :pd_ref => Array{<:Real, 1},
        :qd_ref => Array{<:Real, 1},
        :vnom => Array{<:Real, 1},
        :alpha => Array{<:Real, 1},
        :beta => Array{<:Real, 1},
    ),
    :generator => Dict{Symbol,Type}(
        :obj_name => Any,
        :status => Int,
        :bus => Any,
        :connections => Vector,
        :configuration => String,
        :cost => Vector{<:Real},
        :pg => Array{<:Real, 1},
        :qg => Array{<:Real, 1},
        :pg_min => Array{<:Real, 1},
        :pg_max => Array{<:Real, 1},
        :qg_min => Array{<:Real, 1},
        :qg_max => Array{<:Real, 1},
    ),
    :transformer_nw => Dict{Symbol,Type}(
        :obj_name => Any,
        :status => Int,
        :bus => Array{<:AbstractString, 1},
        :connections => Vector,
        :vnom => Array{<:Real, 1},
        :snom => Array{<:Real, 1},
        :configuration => Array{String, 1},
        :polarity => Array{Bool, 1},
        :xsc => Array{<:Real, 1},
        :rs => Array{<:Real, 1},
        :noloadloss => Real,
        :imag => Real,
        :tm_fix => Array{Array{Bool, 1}, 1},
        :tm => Array{<:Array{<:Real, 1}, 1},
        :tm_min => Array{<:Array{<:Real, 1}, 1},
        :tm_max => Array{<:Array{<:Real, 1}, 1},
        :tm_step => Array{<:Array{<:Real, 1}, 1},
    ),
    :capacitor => Dict{Symbol,Type}(
        :obj_name => Any,
        :status => Int,
        :bus => Any,
        :connections => Vector,
        :configuration => String,
        :qd_ref => Array{<:Real, 1},
        :vnom => Real,
    ),
    :shunt => Dict{Symbol,Type}(
        :obj_name => Any,
        :status => Int,
        :bus => Any,
        :connections => Vector,
        :g_sh => Array{<:Real, 2},
        :b_sh => Array{<:Real, 2},
    ),
    :voltage_source => Dict{Symbol,Type}(
        :obj_name => Any,
        :status => Int,
        :bus => Any,
        :connections => Vector,
        :vm =>Array{<:Real},
        :va =>Array{<:Real},
        :pg_max =>Array{<:Real},
        :pg_min =>Array{<:Real},
        :qg_max =>Array{<:Real},
        :qg_min =>Array{<:Real},
    ),
    :ev => Dict{Symbol,Type}(),
    :storage => Dict{Symbol,Type}(),
    :pv => Dict{Symbol,Type}(),
    :wind => Dict{Symbol,Type}(),
    :switch => Dict{Symbol,Type}(),
    :autotransformer => Dict{Symbol,Type}(),
    :synchronous_generator => Dict{Symbol,Type}(),
    :zip_load => Dict{Symbol,Type}(),
    :grounding => Dict{Symbol,Type}(
        :bus => Any,
        :rg => Real,
        :xg => Real,
    ),
    :boundary => Dict{Symbol,Type}(),
    :meter => Dict{Symbol,Type}()
)

const _eng_model_req_fields= Dict{Symbol,Array{Symbol,1}}(
    :linecode => Array{Symbol,1}([:obj_name, :rs, :xs, :g_fr, :g_to, :b_fr, :b_to]),
    :line => Array{Symbol,1}([:obj_name, :status, :f_bus, :f_connections, :t_bus, :t_connections, :linecode, :length]),
    :bus => Array{Symbol,1}([:obj_name, :status, :terminals, :grounded, :rg, :xg]),
    :load => Array{Symbol,1}([:obj_name, :status, :bus, :connections, :configuration]),
    :generator => Array{Symbol,1}([:obj_name, :status, :bus, :connections]),
    :transformer_nw => Array{Symbol,1}([:obj_name, :status, :bus, :connections, :vnom, :snom, :configuration, :polarity, :xsc, :rs, :noloadloss, :imag, :tm_fix, :tm, :tm_min, :tm_max, :tm_step]),
    :capacitor => Array{Symbol,1}([:obj_name, :status, :bus, :connections, :configuration, :qd_ref, :vnom]),
    :shunt => Array{Symbol,1}([:obj_name, :status, :bus, :connections, :g_sh, :b_sh]),
    :voltage_source => Array{Symbol,1}([:obj_name, :status, :bus, :connections, :vm, :va]),
    :ev => Array{Symbol,1}([]),
    :storage => Array{Symbol,1}([]),
    :pv => Array{Symbol,1}([]),
    :wind => Array{Symbol,1}([]),
    :switch => Array{Symbol,1}([]),
    :autotransformer => Array{Symbol,1}([]),
    :synchronous_generator => Array{Symbol,1}([]),
    :zip_load => Array{Symbol,1}([]),
    :grounding => Array{Symbol,1}([]),
    :boundary => Array{Symbol,1}([]),
    :meter => Array{Symbol,1}([])
)

_eng_model_checks = Dict()


""
function check__eng_model_dtypes(dict, _eng_model_dtypes, comp_type, obj_name)
    for key in keys(dict)
        symb = Symbol(key)
        if haskey(_eng_model_dtypes, symb)
            @assert(isa(dict[key], _eng_model_dtypes[symb]), "$comp_type $obj_name: the property $key should be a $(_eng_model_dtypes[symb]), not a $(typeof(dict[key])).")
        else
            #@assert(false, "$comp_type $obj_name: the property $key is unknown.")
        end
    end
end


""
function add!(data_eng::Dict{String,<:Any}, obj_type::AbstractString, obj_name::AbstractString, object::Dict{String,<:Any})
    # @assert(haskey(object, "obj_name"), "The component does not have an obj_name defined.")
    if !haskey(data_eng, obj_type)
        data_eng[obj_type] = Dict{String,Any}()
    else
        @assert(!haskey(data_eng[obj_type], obj_name), "There is already a $obj_type with name $obj_name.")
    end

    data_eng[obj_type][obj_name] = object
end


""
function _add_unused_kwargs!(object, kwargs)
    for (prop, val) in kwargs
        if !haskey(object, "$prop")
            object["$prop"] = val
        end
    end
end


""
function check_data_model(data)
    for component in keys(_eng_model_dtypes)
        if haskey(data, string(component))
            for (obj_name, object) in data[string(component)]
                if haskey(_eng_model_req_fields, component)
                    for field in _eng_model_req_fields[component]
                        @assert(haskey(object, string(field)), "The property \'$field\' is missing for $component $obj_name.")
                    end
                end
                if haskey(_eng_model_dtypes, component)
                    check__eng_model_dtypes(object, _eng_model_dtypes[component], component, obj_name)
                end
                if haskey(_eng_model_checks, component)
                    _eng_model_checks[component](data, object)
                end
            end
        end
    end
end


""
function create_data_model(; kwargs...)
    data_model = Dict{String, Any}("settings"=>Dict{String, Any}())

    add_kwarg!(data_model["settings"], kwargs, :kv_kvar_scalar, 1e3)

    _add_unused_kwargs!(data_model["settings"], kwargs)

    return data_model
end


""
function _check_same_size(data, keys; context=missing)
    size_comp = size(data[string(keys[1])])
    for key in keys[2:end]
        @assert(all(size(data[string(key)]).==size_comp), "$context: the property $key should have the same size as $(keys[1]).")
    end
end


""
function _check_has_size(data, keys, size_comp; context=missing, allow_missing=true)
    for key in keys
        if haskey(data, key) || !allow_missing
            @assert(all(size(data[string(key)]).==size_comp), "$context: the property $key should have as size $size_comp.")
        end
    end
end


""
function _check_connectivity(data, object; context=missing)
    if haskey(object, "f_bus")
        # two-port element
        _check_bus_and_terminals(data, object["f_bus"], object["f_connections"], context)
        _check_bus_and_terminals(data, object["t_bus"], object["t_connections"], context)
    elseif haskey(object, "bus")
        if isa(object["bus"], Vector)
            for i in 1:length(object["bus"])
                _check_bus_and_terminals(data, object["bus"][i], object["connections"][i], context)
            end
        else
            _check_bus_and_terminals(data, object["bus"], object["connections"], context)
        end
    end
end


""
function _check_bus_and_terminals(data, bus_obj_name, terminals, context=missing)
    @assert(haskey(data, "bus") && haskey(data["bus"], bus_obj_name), "$context: the bus $bus_obj_name is not defined.")
    bus = data["bus"][bus_obj_name]
    for t in terminals
        @assert(t in bus["terminals"], "$context: bus $(bus["obj_name"]) does not have terminal \'$t\'.")
    end
end


""
function _check_has_keys(object, keys; context=missing)
    for key in keys
        @assert(haskey(object, key), "$context: the property $key is missing.")
    end
end


""
function _check_configuration_infer_dim(object; context=missing)
    conf = object["configuration"]
    @assert(conf in ["delta", "wye"], "$context: the configuration should be \'delta\' or \'wye\', not \'$conf\'.")
    return conf=="wye" ? length(object["connections"])-1 : length(object["connections"])
end


# linecode
_eng_model_checks[:linecode] = function check_linecode(data, linecode)
    _check_same_size(linecode, [:rs, :xs, :g_fr, :g_to, :b_fr, :b_to])
end

function create_linecode(; kwargs...)
    linecode = Dict{String,Any}()

    n_conductors = 0
    for key in [:rs, :xs, :g_fr, :g_to, :b_fr, :b_to]
        if haskey(kwargs, key)
            n_conductors = size(kwargs[key])[1]
        end
    end
    add_kwarg!(linecode, kwargs, :rs, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(linecode, kwargs, :xs, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(linecode, kwargs, :g_fr, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(linecode, kwargs, :b_fr, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(linecode, kwargs, :g_to, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(linecode, kwargs, :b_to, fill(0.0, n_conductors, n_conductors))

    _add_unused_kwargs!(linecode, kwargs)

    return linecode
end

# line
_eng_model_checks[:line] = function check_line(data, line)
    i = line["obj_name"]

    # for now, always require a line code
    if haskey(line, "linecode")
        # line is defined with a linecode
        @assert(haskey(line, "length"), "line $i: a line defined through a linecode, should have a length property.")

        linecode_obj_name = line["linecode"]
        @assert(haskey(data, "linecode") && haskey(data["linecode"], "$linecode_obj_name"), "line $i: the linecode $linecode_obj_name is not defined.")
        linecode = data["linecode"]["$linecode_obj_name"]

        for key in ["n_conductors", "rs", "xs", "g_fr", "g_to", "b_fr", "b_to"]
            @assert(!haskey(line, key), "line $i: a line with a linecode, should not specify $key; this is already done by the linecode.")
        end

        N = size(linecode["rs"])[1]
        @assert(length(line["f_connections"])==N, "line $i: the number of terminals should match the number of conductors in the linecode.")
        @assert(length(line["t_connections"])==N, "line $i: the number of terminals should match the number of conductors in the linecode.")
    else
        # normal line
        @assert(!haskey(line, "length"), "line $i: length only makes sense for linees defined through linecodes.")
        for key in ["n_conductors", "rs", "xs", "g_fr", "g_to", "b_fr", "b_to"]
            @assert(haskey(line, key), "line $i: a line without linecode, should specify $key.")
        end
    end

    _check_connectivity(data, line, context="line $(line["obj_name"])")
end


function create_line(; kwargs...)
    line = Dict{String,Any}()

    add_kwarg!(line, kwargs, :status, 1)
    add_kwarg!(line, kwargs, :f_connections, collect(1:4))
    add_kwarg!(line, kwargs, :t_connections, collect(1:4))

    N = length(line["f_connections"])
    add_kwarg!(line, kwargs, :angmin, fill(-60/180*pi, N))
    add_kwarg!(line, kwargs, :angmax, fill( 60/180*pi, N))

    # if no linecode, then populate loss parameters with zero
    if !haskey(kwargs, :linecode)
        n_conductors = 0
        for key in [:rs, :xs, :g_fr, :g_to, :b_fr, :b_to]
            if haskey(kwargs, key)
                n_conductors = size(kwargs[key])[1]
            end
        end
        add_kwarg!(line, kwargs, :rs, fill(0.0, n_conductors, n_conductors))
        add_kwarg!(line, kwargs, :xs, fill(0.0, n_conductors, n_conductors))
        add_kwarg!(line, kwargs, :g_fr, fill(0.0, n_conductors, n_conductors))
        add_kwarg!(line, kwargs, :b_fr, fill(0.0, n_conductors, n_conductors))
        add_kwarg!(line, kwargs, :g_to, fill(0.0, n_conductors, n_conductors))
        add_kwarg!(line, kwargs, :b_to, fill(0.0, n_conductors, n_conductors))
    end

    _add_unused_kwargs!(line, kwargs)

    return line
end

# Bus
_eng_model_checks[:bus] = function check_bus(data, bus)
    obj_name = bus["obj_name"]

    _check_same_size(bus, ["grounded", "rg", "xg"], context="bus $obj_name")

    N = length(bus["terminals"])
    _check_has_size(bus, ["vm_max", "vm_min", "vm", "va"], N, context="bus $obj_name")

    if haskey(bus, "neutral")
        assert(haskey(bus, "phases"), "bus $obj_name: has a neutral, but no phases.")
    end
end

function create_bus(; kwargs...)
    bus = Dict{String,Any}()

    add_kwarg!(bus, kwargs, :status, 1)
    add_kwarg!(bus, kwargs, :terminals, collect(1:4))
    add_kwarg!(bus, kwargs, :grounded, [])
    add_kwarg!(bus, kwargs, :bus_type, 1)
    add_kwarg!(bus, kwargs, :rg, Array{Float64, 1}())
    add_kwarg!(bus, kwargs, :xg, Array{Float64, 1}())

    _add_unused_kwargs!(bus, kwargs)

    return bus
end

# Load
_eng_model_checks[:load] = function check_load(data, load)
    obj_name = load["obj_name"]

    N = _check_configuration_infer_dim(load; context="load $obj_name")

    model = load["model"]
    @assert(model in ["constant_power", "constant_impedance", "constant_current", "exponential"])
    if model=="constant_power"
        _check_has_keys(load, ["pd", "qd"], context="load $obj_name, $model:")
        _check_has_size(load, ["pd", "qd"], N, context="load $obj_name, $model:")
    elseif model=="exponential"
        _check_has_keys(load, ["pd_ref", "qd_ref", "vnom", "alpha", "beta"], context="load $obj_name, $model")
        _check_has_size(load, ["pd_ref", "qd_ref", "vnom", "alpha", "beta"], N, context="load $obj_name, $model:")
    else
        _check_has_keys(load, ["pd_ref", "qd_ref", "vnom"], context="load $obj_name, $model")
        _check_has_size(load, ["pd_ref", "qd_ref", "vnom"], N, context="load $obj_name, $model:")
    end

    _check_connectivity(data, load; context="load $obj_name")
end


function create_load(; kwargs...)
    load = Dict{String,Any}()

    add_kwarg!(load, kwargs, :status, 1)
    add_kwarg!(load, kwargs, :configuration, "wye")
    add_kwarg!(load, kwargs, :connections, load["configuration"]=="wye" ? [1, 2, 3, 4] : [1, 2, 3])
    add_kwarg!(load, kwargs, :model, "constant_power")
    if load["model"]=="constant_power"
        add_kwarg!(load, kwargs, :pd, fill(0.0, 3))
        add_kwarg!(load, kwargs, :qd, fill(0.0, 3))
    else
        add_kwarg!(load, kwargs, :pd_ref, fill(0.0, 3))
        add_kwarg!(load, kwargs, :qd_ref, fill(0.0, 3))
    end

    _add_unused_kwargs!(load, kwargs)

    return load
end

# generator
_eng_model_checks[:generator] = function check_generator(data, generator)
    obj_name = generator["obj_name"]

    N = _check_configuration_infer_dim(generator; context="generator $obj_name")
    _check_has_size(generator, ["pd", "qd", "pd_min", "pd_max", "qd_min", "qd_max"], N, context="generator $obj_name")

    _check_connectivity(data, generator; context="generator $obj_name")
end

function create_generator(; kwargs...)
    generator = Dict{String,Any}()

    add_kwarg!(generator, kwargs, :status, 1)
    add_kwarg!(generator, kwargs, :configuration, "wye")
    add_kwarg!(generator, kwargs, :cost, [1.0, 0.0]*1E-3)
    add_kwarg!(generator, kwargs, :connections, generator["configuration"]=="wye" ? [1, 2, 3, 4] : [1, 2, 3])

    _add_unused_kwargs!(generator, kwargs)

    return generator
end


# Transformer, n-windings three-phase lossy
_eng_model_checks[:transformer_nw] = function check_transformer_nw(data, trans)
    obj_name = trans["obj_name"]
    nrw = length(trans["bus"])
    _check_has_size(trans, ["bus", "connections", "vnom", "snom", "configuration", "polarity", "rs", "tm_fix", "tm_set", "tm_min", "tm_max", "tm_step"], nrw, context="trans $obj_name")
    @assert(length(trans["xsc"])==(nrw^2-nrw)/2)

    nphs = []
    for w in 1:nrw
        @assert(trans["configuration"][w] in ["wye", "delta"])
        conf = trans["configuration"][w]
        conns = trans["connections"][w]
        nph = conf=="wye" ? length(conns)-1 : length(conns)
        @assert(all(nph.==nphs), "transformer $obj_name: winding $w has a different number of phases than the previous ones.")
        push!(nphs, nph)
        #TODO check length other properties
    end

    _check_connectivity(data, trans; context="transformer_nw $obj_name")
end


function create_transformer_nw(; kwargs...)
    trans = Dict{String,Any}()

    @assert(haskey(kwargs, :bus), "You have to specify at least the buses.")
    n_windings = length(kwargs[:bus])
    add_kwarg!(trans, kwargs, :status, 1)
    add_kwarg!(trans, kwargs, :configuration, fill("wye", n_windings))
    add_kwarg!(trans, kwargs, :polarity, fill(true, n_windings))
    add_kwarg!(trans, kwargs, :rs, fill(0.0, n_windings))
    add_kwarg!(trans, kwargs, :xsc, fill(0.0, n_windings^2-n_windings))
    add_kwarg!(trans, kwargs, :noloadloss, 0.0)
    add_kwarg!(trans, kwargs, :imag, 0.0)
    add_kwarg!(trans, kwargs, :tm, fill(fill(1.0, 3), n_windings))
    add_kwarg!(trans, kwargs, :tm_min, fill(fill(0.9, 3), n_windings))
    add_kwarg!(trans, kwargs, :tm_max, fill(fill(1.1, 3), n_windings))
    add_kwarg!(trans, kwargs, :tm_step, fill(fill(1/32, 3), n_windings))
    add_kwarg!(trans, kwargs, :tm_fix, fill(fill(true, 3), n_windings))

    _add_unused_kwargs!(trans, kwargs)

    return trans
end

#
# # Transformer, two-winding three-phase
#
# _eng_model_dtypes[:transformer_2w_obj_nameeal] = Dict(
#     :obj_name => Any,
#     :f_bus => String,
#     :t_bus => String,
#     :configuration => String,
#     :f_terminals => Array{Int, 1},
#     :t_terminals => Array{Int, 1},
#     :tm_nom => Real,
#     :tm_set => Real,
#     :tm_min => Real,
#     :tm_max => Real,
#     :tm_step => Real,
#     :tm_fix => Real,
# )
#
#
# _eng_model_checks[:transformer_2w_obj_nameeal] = function check_transformer_2w_obj_nameeal(data, trans)
# end
#
#
# function create_transformer_2w_obj_nameeal(obj_name, f_bus, t_bus, tm_nom; kwargs...)
#     trans = Dict{String,Any}()
#     trans["obj_name"] = obj_name
#     trans["f_bus"] = f_bus
#     trans["t_bus"] = t_bus
#     trans["tm_nom"] = tm_nom
#     add_kwarg!(trans, kwargs, :configuration, "wye")
#     add_kwarg!(trans, kwargs, :f_terminals, trans["configuration"]=="wye" ? [1, 2, 3, 4] : [1, 2, 3])
#     add_kwarg!(trans, kwargs, :t_terminals, [1, 2, 3, 4])
#     add_kwarg!(trans, kwargs, :tm_set, 1.0)
#     add_kwarg!(trans, kwargs, :tm_min, 0.9)
#     add_kwarg!(trans, kwargs, :tm_max, 1.1)
#     add_kwarg!(trans, kwargs, :tm_step, 1/32)
#     add_kwarg!(trans, kwargs, :tm_fix, true)
#     return trans
# end


# Capacitor
_eng_model_checks[:capacitor] = function check_capacitor(data, cap)
    obj_name = cap["obj_name"]
    N = length(cap["connections"])
    config = cap["configuration"]
    if config=="wye"
        @assert(length(cap["qd_ref"])==N-1, "capacitor $obj_name: qd_ref should have $(N-1) elements.")
    else
        @assert(length(cap["qd_ref"])==N, "capacitor $obj_name: qd_ref should have $N elements.")
    end
    @assert(config in ["delta", "wye", "wye-grounded", "wye-floating"])
    if config=="delta"
        @assert(N>=3, "Capacitor $obj_name: delta-connected capacitors should have at least 3 elements.")
    end

    _check_connectivity(data, cap; context="capacitor $obj_name")
end


function create_capacitor(; kwargs...)
    cap = Dict{String,Any}()

    add_kwarg!(cap, kwargs, :status, 1)
    add_kwarg!(cap, kwargs, :configuration, "wye")
    add_kwarg!(cap, kwargs, :connections, collect(1:4))
    add_kwarg!(cap, kwargs, :qd_ref, fill(0.0, 3))

    _add_unused_kwargs!(cap, kwargs)

    return cap
end


# Shunt
_eng_model_checks[:shunt] = function check_shunt(data, shunt)
    _check_connectivity(data, shunt; context="shunt $obj_name")

end


function create_shunt(; kwargs...)
    shunt = Dict{String,Any}()

    N = length(kwargs[:connections])

    add_kwarg!(shunt, kwargs, :status, 1)
    add_kwarg!(shunt, kwargs, :g_sh, fill(0.0, N, N))
    add_kwarg!(shunt, kwargs, :b_sh, fill(0.0, N, N))

    _add_unused_kwargs!(shunt, kwargs)

    return shunt
end


# voltage source
_eng_model_checks[:voltage_source] = function check_voltage_source(data, vs)
    obj_name = vs["obj_name"]
    _check_connectivity(data, vs; context="voltage source $obj_name")
    N = length(vs["connections"])
    _check_has_size(vs, ["vm", "va", "pg_max", "pg_min", "qg_max", "qg_min"], N, context="voltage source $obj_name")

end


function create_voltage_source(; kwargs...)
    vs = Dict{String,Any}()

    add_kwarg!(vs, kwargs, :status, 1)
    add_kwarg!(vs, kwargs, :connections, collect(1:3))

    _add_unused_kwargs!(vs, kwargs)

    return vs
end


# create add_comp! methods
for comp in keys(_eng_model_dtypes)
    eval(Meta.parse("add_$(comp)!(data_model, name; kwargs...) = add!(data_model, \"$comp\", name, create_$comp(; kwargs...))"))
end
