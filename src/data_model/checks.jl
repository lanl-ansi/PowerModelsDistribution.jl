"data check functions for the engineering data model"
const _eng_model_checks = Dict{Symbol,Symbol}(
    :bus => :_check_bus,
    :linecode => :_check_linecode,
    # :xfmrcode => :_check_xfmrcode,
    :line => :_check_line,
    :transformer => :_check_transformer,
    # :switch => :_check_switch,
    # :line_reactor => :_check_line_reactor,
    # :series_capacitor => :_check_series_capacitor,
    :load => :_check_load,
    :shunt_capacitor => :_check_shunt_capacitor,
    # :shunt_reactor => :_check_shunt_reactor,
    :shunt => :_check_shunt,
    :generator => :_check_generator,
    :voltage_source => :_check_voltage_source,
    # :solar => :_check_solar,
    # :storage => :_check_storage,
    # :grounding => :_check_grounding,
)

"Data types of accepted fields in the engineering data model"
const _eng_model_dtypes = Dict{Symbol,Dict{Symbol,Type}}(
    :linecode => Dict{Symbol,Type}(
        :rs => Array{<:Real, 2},
        :xs => Array{<:Real, 2},
        :g_fr => Array{<:Real, 2},
        :g_to => Array{<:Real, 2},
        :b_fr => Array{<:Real, 2},
        :b_to => Array{<:Real, 2}
    ),
    :line => Dict{Symbol,Type}(
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
    :transformer => Dict{Symbol,Type}(
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
    :shunt_capacitor => Dict{Symbol,Type}(
        :status => Int,
        :bus => Any,
        :connections => Vector,
        :configuration => String,
        :qd_ref => Array{<:Real, 1},
        :vnom => Real,
    ),
    :shunt => Dict{Symbol,Type}(
        :status => Int,
        :bus => Any,
        :connections => Vector,
        :g_sh => Array{<:Real, 2},
        :b_sh => Array{<:Real, 2},
    ),
    :voltage_source => Dict{Symbol,Type}(
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
    :solar => Dict{Symbol,Type}(),
    :storage => Dict{Symbol,Type}(),
    :switch => Dict{Symbol,Type}(),
    :grounding => Dict{Symbol,Type}(
        :bus => Any,
        :rg => Real,
        :xg => Real,
    ),
    # :ev => Dict{Symbol,Type}(),
    # :wind => Dict{Symbol,Type}(),
    # :autotransformer => Dict{Symbol,Type}(),
    # :synchronous_generator => Dict{Symbol,Type}(),
    # :zip_load => Dict{Symbol,Type}(),
    # :boundary => Dict{Symbol,Type}(),
    # :meter => Dict{Symbol,Type}()
)

"required fields in the engineering data model"
const _eng_model_req_fields= Dict{Symbol,Vector{Symbol}}(
    :bus => Vector{Symbol}([
        :status, :terminals, :grounded, :rg, :xg
    ]),
    :linecode => Vector{Symbol}([
        :rs, :xs, :g_fr, :g_to, :b_fr, :b_to
     ]),
    :xfmrcode => Vector{Symbol}([
        :vnom, :snom, :configuration, :polarity, :xsc, :rs, :noloadloss, :imag,
        :tm_fix, :tm, :tm_min, :tm_max, :tm_step
    ]),
    :line => Vector{Symbol}([
        :status, :f_bus, :f_connections, :t_bus, :t_connections, :linecode,
        :length
    ]),
    :transformer => Vector{Symbol}([
        :status, :bus, :connections, :vnom, :snom, :configuration, :polarity,
        :xsc, :rs, :noloadloss, :imag, :tm_fix, :tm, :tm_min, :tm_max,
        :tm_step
    ]),
    :switch => Vector{Symbol}([]),
    :line_reactor => Vector{Symbol}([]),
    :series_capacitor => Vector{Symbol}([]),
    :load => Vector{Symbol}([
        :status, :bus, :connections, :configuration
    ]),
    :shunt_capacitor => Vector{Symbol}([
        :status, :bus, :connections, :configuration, :qd_ref, :vnom
    ]),
    :shunt_reactor => Vector{Symbol}([]),
    :shunt => Vector{Symbol}([
        :status, :bus, :connections, :g_sh, :b_sh
    ]),
    :generator => Vector{Symbol}([
        :status, :bus, :connections
    ]),
    :voltage_source => Vector{Symbol}([
        :status, :bus, :connections, :vm, :va
        ]),
    :solar => Vector{Symbol}([]),
    :storage => Vector{Symbol}([]),
    :grounding => Vector{Symbol}([]),
    # Future Components
    # :ev => Vector{Symbol}([]),
    # :wind => Vector{Symbol}([]),
    # :autotransformer => Vector{Symbol}([]),
    # :zip_load => Vector{Symbol}([]),
    # :synchronous_generator => Vector{Symbol}([]),
    # :boundary => Vector{Symbol}([]),
    # :meter => Vector{Symbol}([])
)


"checks that an engineering model component has the correct data types"
function _check_eng_component_dtypes(data_eng::Dict{String,<:Any}, component_type::String, component_name::String; additional_dtypes=Dict{Symbol,Type}())
    if haskey(_eng_model_dtypes, Symbol(component_type))
        dtypes = merge(_eng_model_dtypes[Symbol(component_type)], additional_dtypes)
    else
        dtypes = additional_dtypes
    end

    if haskey(data_eng, component_type) && haskey(data_eng[component_type], component_name)
        component = data_eng[component_type][component_name]

        for (field, dtype) in dtypes
            if haskey(component, string(field))
                @assert isa(component[field], dtype) "$component_type $component_name: the property $field should be a $dtype, not a $(typeof(component[field]))"
            end
        end
    else
        Memento.warn(_LOGGER, "$component_type $component_name does not exist")
    end
end



"check that all data in `fields` have the same size"
function _check_same_size(component::Dict{String,<:Any}, fields::Vector{String}; context::Union{String,Missing}=missing)
    @assert length(unique([size(component[string(field)]) for field in fields])) == 1 "$context: not all properties are the same size"
end


"check that `fields` has size `data_size`"
function _check_has_size(component::Dict{String,<:Any}, fields::Vector{String}, data_size::Tuple; context::Union{String,Missing}=missing, allow_missing::Bool=true)
    for key in keys
        if haskey(component, key) || !allow_missing
            @assert all(size(component[string(key)]).==data_size) "$context: the property $key should have as size $data_size"
        end
    end
end


"checks connectivity of object"
function _check_connectivity(data_eng::Dict{String,<:Any}, object::Dict{String,<:Any}; context::Union{String,Missing}=missing)
    if haskey(object, "f_bus")
        # two-port element
        _check_bus_and_terminals(data_eng, object["f_bus"], object["f_connections"], context)
        _check_bus_and_terminals(data_eng, object["t_bus"], object["t_connections"], context)
    elseif haskey(object, "bus")
        if isa(object["bus"], Vector)
            for i in 1:length(object["bus"])
                _check_bus_and_terminals(data_eng, object["bus"][i], object["connections"][i], context)
            end
        else
            _check_bus_and_terminals(data_eng, object["bus"], object["connections"], context)
        end
    end
end


"checks `bus_name` exists and has `terminals`"
function _check_bus_and_terminals(data_eng::Dict{String,<:Any}, bus_name::Any, terminals::Vector{Int}, context::Union{String,Missing}=missing)
    @assert haskey(data_eng, "bus") && haskey(data_eng["bus"], bus_name) "$context: the bus $bus_name is not defined."

    bus = data_eng["bus"][bus_name]
    for t in terminals
        @assert t in bus["terminals"] "$context: bus $(bus["obj_name"]) does not have terminal \'$t\'."
    end
end


"checks that a component has `fields`"
function _check_has_keys(object::Dict{String,<:Any}, fields::Vector{String}; context::Union{String,Missing}=missing)
    for key in fields
        @assert haskey(object, key) "$context: the property $key is missing."
    end
end


"checks the connection configuration and infers the dimensions of the connection (number of connected terminals)"
function _check_configuration_infer_dim(object::Dict{String,<:Any}; context::Union{String,Missing}=missing)::Int
    conf = object["configuration"]
    @assert conf in ["delta", "wye"] "$context: the configuration should be \'delta\' or \'wye\', not \'$conf\'."

    return conf=="wye" ? length(object["connections"])-1 : length(object["connections"])
end


"checks the engineering data model for correct data types, required fields and applies default checks"
function check_eng_data_model(data_eng::Dict{String,<:Any})
    for (component_type, components) in data_eng
        if isa(components, Dict)
            for (name, component) in keys(components)
                _check_eng_component_dtypes(data_eng, component_type, name)

                for field in get(_eng_model_req_fields, Symbol(component_type), Vector{Symbol}([]))
                    @assert haskey(component, string(field)) "The property \'$field\' is missing on $component_type $name"
                end

                for check in get(_eng_model_checks, Symbol(component_type), missing)
                    if !ismissing(check)
                        @eval $(check)(data_eng, name)
                    end
                end
            end
        end
    end
end


"bus data checks"
function _check_bus(data_eng::Dict{String,<:Any}, name::Any)
    bus = data_eng["bus"][name]

    _check_same_size(bus, ["grounded", "rg", "xg"], context="bus $name")

    N = length(bus["terminals"])
    _check_has_size(bus, ["vm_max", "vm_min", "vm", "va"], N, context="bus $name")

    if haskey(bus, "neutral")
        @assert haskey(bus, "phases") "bus $name: has a neutral, but no phases."
    end
end


"load data checks"
function _check_load(data_eng::Dict{String,<:Any}, name::Any)
    load = data_eng["load"][name]

    N = _check_configuration_infer_dim(load; context="load $name")

    model = load["model"]
    @assert model in ["constant_power", "constant_impedance", "constant_current", "exponential"]

    if model=="constant_power"
        _check_has_keys(load, ["pd", "qd"], context="load $name, $model:")
        _check_has_size(load, ["pd", "qd"], N, context="load $name, $model:")
    elseif model=="exponential"
        _check_has_keys(load, ["pd_ref", "qd_ref", "vnom", "alpha", "beta"], context="load $name, $model")
        _check_has_size(load, ["pd_ref", "qd_ref", "vnom", "alpha", "beta"], N, context="load $name, $model:")
    else
        _check_has_keys(load, ["pd_ref", "qd_ref", "vnom"], context="load $name, $model")
        _check_has_size(load, ["pd_ref", "qd_ref", "vnom"], N, context="load $name, $model:")
    end

    _check_connectivity(data_eng, load; context="load $name")
end


"linecode data checks"
function _check_linecode(data_eng::Dict{String,<:Any}, name::Any)
    _check_same_size(data_eng["linecode"][name], [:rs, :xs, :g_fr, :g_to, :b_fr, :b_to])
end


"line data checks"
function _check_line(data_eng::Dict{String,<:Any}, name::Any)
    line = data_eng["line"][name]

    # for now, always require a line code
    if haskey(line, "linecode")
        # line is defined with a linecode
        @assert haskey(line, "length") "line $name: a line defined through a linecode, should have a length property."

        linecode_obj_name = line["linecode"]
        @assert haskey(data_eng, "linecode") && haskey(data_eng["linecode"], "$linecode_obj_name")  "line $name: the linecode $linecode_obj_name is not defined."
        linecode = data_eng["linecode"]["$linecode_obj_name"]

        for key in ["n_conductors", "rs", "xs", "g_fr", "g_to", "b_fr", "b_to"]
            @assert !haskey(line, key) "line $name: a line with a linecode, should not specify $key; this is already done by the linecode."
        end

        N = size(linecode["rs"])[1]
        @assert length(line["f_connections"])==N "line $name: the number of terminals should match the number of conductors in the linecode."
        @assert length(line["t_connections"])==N "line $name: the number of terminals should match the number of conductors in the linecode."
    else
        # normal line
        @assert !haskey(line, "length") "line $name: length only makes sense for linees defined through linecodes."
        for key in ["n_conductors", "rs", "xs", "g_fr", "g_to", "b_fr", "b_to"]
            @assert haskey(line, key) "line $name: a line without linecode, should specify $key."
        end
    end

    _check_connectivity(data_eng, line, context="line $(name)")
end


"generator data checks"
function _check_generator(data_eng::Dict{String,<:Any}, name::Any)
    generator = data_eng["generator"][name]

    N = _check_configuration_infer_dim(generator; context="generator $name")
    _check_has_size(generator, ["pd", "qd", "pd_min", "pd_max", "qd_min", "qd_max"], N, context="generator $name")

    _check_connectivity(data_eng, generator; context="generator $name")
end


"Transformer, n-windings three-phase lossy data checks"
function _check_transformer(data_eng::Dict{String,<:Any}, name::Any)
    transformer = data_eng["transformer"][name]

    nrw = length(transformer["bus"])
    _check_has_size(transformer, ["bus", "connections", "vnom", "snom", "configuration", "polarity", "rs", "tm_fix", "tm_set", "tm_min", "tm_max", "tm_step"], nrw, context="trans $name")

    @assert length(transformer["xsc"])==(nrw^2-nrw)/2

    nphs = []
    for w in 1:nrw
        @assert transformer["configuration"][w] in ["wye", "delta"]

        conf = transformer["configuration"][w]
        conns = transformer["connections"][w]
        nph = conf=="wye" ? length(conns)-1 : length(conns)
        @assert all(nph.==nphs) "transformer $name: winding $w has a different number of phases than the previous ones."

        push!(nphs, nph)
        #TODO check length other properties
    end

    _check_connectivity(data_eng, transformer; context="transformer_nw $name")
end


"shunt capacitor data checks"
function _check_shunt_capacitor(data_eng::Dict{String,<:Any}, name::Any)
    shunt_capacitor = data_eng["shunt_capacitor"][name]

    N = length(shunt_capacitor["connections"])
    config = shunt_capacitor["configuration"]
    if config=="wye"
        @assert length(shunt_capacitor["qd_ref"])==N-1 "capacitor $name: qd_ref should have $(N-1) elements."
    else
        @assert length(shunt_capacitor["qd_ref"])==N "capacitor $name: qd_ref should have $N elements."
    end

    @assert config in ["delta", "wye", "wye-grounded", "wye-floating"]

    if config=="delta"
        @assert N>=3 "Capacitor $name: delta-connected capacitors should have at least 3 elements."
    end

    _check_connectivity(data_eng, shunt_capacitor; context="capacitor $name")
end


"shunt data checks"
function _check_shunt(data_eng::Dict{String,<:Any}, name::Any)
    shunt = data_eng["shunt"][name]

    _check_connectivity(data_eng, shunt; context="shunt $name")
end


"voltage source data checks"
function _check_voltage_source(data_eng::Dict{String,<:Any}, name::Any)
    voltage_source = data_eng["voltage_source"][name]

    _check_connectivity(data_eng, voltage_source; context="voltage source $name")

    N = length(voltage_source["connections"])
    _check_has_size(voltage_source, ["vm", "va", "pg_max", "pg_min", "qg_max", "qg_min"], N, context="voltage source $name")
end
