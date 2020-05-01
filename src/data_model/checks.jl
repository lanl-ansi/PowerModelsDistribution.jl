"data check functions for the engineering data model"
const _eng_model_checks = Dict{Symbol,Symbol}(
    :bus => :_check_bus,
    :linecode => :_check_linecode,
    # :xfmrcode => :_check_xfmrcode,
    :line => :_check_line,
    :transformer => :_check_transformer,
    # :switch => :_check_switch,
    :load => :_check_load,
    :shunt => :_check_shunt,
    :generator => :_check_generator,
    :voltage_source => :_check_voltage_source,
    # :solar => :_check_solar,
    # :storage => :_check_storage,
)

"Data types of accepted fields in the engineering data model"
const _eng_model_dtypes = Dict{Symbol,Dict{Symbol,Type}}(
    :bus => Dict{Symbol,Type}(
        :status => Int,
        :terminals => Vector{Any},
        :phases => Vector{Any},
        :neutral => Any,
        :grounded => Vector{Any},
        :rg => Vector{<:Real},
        :xg => Vector{<:Real},
        :vm_pn_lb => Real,
        :vm_pn_ub => Real,
        :vm_pp_lb => Real,
        :vm_pp_ub => Real,
        :vm_lb => Vector{<:Real},
        :vm_ub => Vector{<:Real},
        :vm => Vector{<:Real},
        :va => Vector{<:Real},
    ),
    :line => Dict{Symbol,Type}(
        :status => Int,
        :f_bus => Any,
        :t_bus => Any,
        :f_connections => Vector{Any},
        :t_connections => Vector{Any},
        :linecode => String,
        :length => Real,
        :cm_ub =>Vector{<:Real},
        :sm_ub =>Vector{<:Real},
        :vad_lb=>Vector{<:Real},
        :vad_ub=>Vector{<:Real},
        :rs => Matrix{<:Real},
        :xs => Matrix{<:Real},
        :g_fr => Matrix{<:Real},
        :g_to => Matrix{<:Real},
        :b_fr => Matrix{<:Real},
        :b_to => Matrix{<:Real},
    ),
    :transformer => Dict{Symbol,Type}(
        :status => Int,
        :bus => Vector{Any},
        :connections => Vector{Any},
        :vnom => Vector{<:Real},
        :snom => Vector{<:Real},
        :configuration => Vector{String},
        :polarity => Vector{Int},
        :xsc => Vector{<:Real},
        :rs => Vector{<:Real},
        :noloadloss => Real,
        :imag => Real,
        :tm_fix => Vector{Union{Vector{Int},Int}},
        :tm => Vector{Union{Vector{<:Real},<:Real}},
        :tm_min => Vector{Union{Vector{<:Real},<:Real}},
        :tm_max => Vector{Union{Vector{<:Real},<:Real}},
        :tm_step => Vector{Union{Vector{<:Real},<:Real}},
        :tm_nom => Union{Vector{<:Real}, Real},
        :f_bus => Any,
        :t_bus => Any,
        :f_connections => Vector{Any},
        :t_connections => Vector{Any},
        :configuration => String,
        :xfmrcode => String,
    ),
    :switch => Dict{Symbol,Type}(
        :status => Int,
        :f_bus => Any,
        :t_bus => Any,
        :f_connections => Vector{Any},
        :t_connections => Vector{Any},
        :cm_ub => Vector{<:Real},
        :sm_ub => Vector{<:Real},
        :rs => Matrix{<:Real},
        :xs => Matrix{<:Real},
        :g_fr => Matrix{<:Real},
        :b_fr => Matrix{<:Real},
        :g_to => Matrix{<:Real},
        :b_to => Matrix{<:Real},
        :state => Int,
    ),
    :fuse => Dict{Symbol,Type}(
        :status => Int,
        :f_bus => Any,
        :t_bus => Any,
        :f_connections => Vector{Any},
        :t_connections => Vector{Any},
        :cm_ub => Vector{<:Real},
        :sm_ub => Vector{<:Real},
        :rs => Matrix{<:Real},
        :xs => Matrix{<:Real},
        :g_fr => Matrix{<:Real},
        :b_fr => Matrix{<:Real},
        :g_to => Matrix{<:Real},
        :b_to => Matrix{<:Real},
        :state => Int,
        :fuse_curve => String,
        :minimum_melting_curve => String,
    ),
    :line_reactor => Dict{Symbol,Type}(),
    :series_capacitor => Dict{Symbol,Type}(),
    :shunt => Dict{Symbol,Type}(
        :status => Int,
        :bus => Any,
        :connections => Vector{Any},
        :configuration => String,
        :gs => Matrix{<:Real},
        :bs => Matrix{<:Real},
        :vnom => Real,
    ),
    :shunt_capacitor => Dict{Symbol,Type}(
        :status => Int,
        :bus => Any,
        :connections => Vector{Any},
        :configuration => String,
        :bs => Matrix{<:Real},
        :vnom => Real,
    ),
    :shunt_reactor => Dict{Symbol,Type}(
        :status => Int,
        :bus => Any,
        :connections => Vector{Any},
        :configuration => String,
        :bs => Matrix{<:Real},
        :vnom => Real,
    ),
    :load => Dict{Symbol,Type}(
        :status => Int,
        :bus => Any,
        :connections => Vector{Any},
        :configuration => String,
        :model => String,
        :pd_nom => Vector{<:Real},
        :qd_nom => Vector{<:Real},
        :vnom => Real,
        :pd_exp => Real,
        :qd_exp => Real,
        :pd_nom_z => Real,
        :pd_nom_i => Real,
        :pd_nom_p => Real,
        :qd_nom_z => Real,
        :qd_nom_i => Real,
        :qd_nom_p => Real,
    ),
    :generator => Dict{Symbol,Type}(
        :status => Int,
        :bus => Any,
        :connections => Vector{Any},
        :configuration => String,
        :model => Int,
        :pg => Vector{<:Real},
        :qg => Vector{<:Real},
        :pg_lb => Vector{<:Real},
        :pg_ub => Vector{<:Real},
        :qg_lb => Vector{<:Real},
        :qg_ub => Vector{<:Real},
        :cost_pg_parameters => Vector{<:Real},
        :cost_pg_model => Int,
    ),
    :solar => Dict{Symbol,Type}(
        :status => Int,
        :bus => Any,
        :connections => Vector{Any},
        :configuration => String,
        :pg => Vector{<:Real},
        :qg => Vector{<:Real},
        :pg_lb => Vector{<:Real},
        :pg_ub => Vector{<:Real},
        :qg_lb => Vector{<:Real},
        :qg_ub => Vector{<:Real},
        :cost_pg_parameters => Vector{<:Real},
        :cost_pg_model => Int,
    ),
    :storage => Dict{Symbol,Type}(
        :status => Int,
        :bus => Any,
        :connections => Vector{Any},
        :configuration => String,
        :energy => Real,
        :energy_ub => Real,
        :charge_ub => Real,
        :sm_ub => Vector{<:Real},
        :cm_ub => Vector{<:Real},
        :charge_efficiency => Real,
        :discharge_efficiency => Real,
        :qs_lb => Vector{<:Real},
        :qs_ub => Vector{<:Real},
        :rs => Vector{<:Real},
        :xs => Vector{<:Real},
        :pex => Real,
        :qex => Real,
    ),
    :voltage_source => Dict{Symbol,Type}(
        :status => Int,
        :bus => Any,
        :connections => Vector{Any},
        :configuration => String,
        :vm => Vector{<:Real},
        :va => Real,
        :rs => Matrix{<:Real},
        :xs => Matrix{<:Real},
    ),
    :linecode => Dict{Symbol,Type}(
        :rs => Matrix{<:Real},
        :xs => Matrix{<:Real},
        :g_fr => Matrix{<:Real},
        :g_to => Matrix{<:Real},
        :b_fr => Matrix{<:Real},
        :b_to => Matrix{<:Real}
    ),
    :xfmrcode => Dict{Symbol,Type}(
        :configurations => Vector{String},
        :xsc => Vector{<:Real},
        :rs => Vector{<:Real},
        :tm_nom => Vector{<:Real},
        :tm_ub => Vector{<:Real},
        :tm_lb => Vector{<:Real},
        :tm_step => Vector{<:Real},
        :tm_set => Vector{<:Real},
        :tm_fix => Vector{<:Real},
    ),
    :curve => Dict{Symbol,Type}(
        :curve => Function,
    ),
    :time_series => Dict{Symbol,Type}(
        :time => Vector{<:Real},
        :values => Vector{<:Real},
        :replace => Bool,
    ),
    # Future Components
    # :ev => Dict{Symbol,Type}(),
    # :wind => Dict{Symbol,Type}(),
    # :autotransformer => Dict{Symbol,Type}(),
    # :meter => Dict{Symbol,Type}()
)

"required fields in the engineering data model"
const _eng_model_req_fields= Dict{Symbol,Vector{Symbol}}(
    :bus => Vector{Symbol}([
        :status, :terminals, :grounded, :rg, :xg,
    ]),
    :line => Vector{Symbol}([
        :status, :f_bus, :t_bus, :f_connections, :t_connections, :length,
    ]),
    :transformer => Vector{Symbol}([
        :status, :configurations, :vnom, :snom, :polarity, :xsc, :rs,
        :noloadloss, :imag, :tm_fix, :tm_set, :tm_step,
    ]),
    :switch => Vector{Symbol}([
        :status, :f_bus, :t_bus, :f_connections, :t_connections,
    ]),
    :fuse => Vector{Symbol}([
        :status, :f_bus, :t_bus, :f_connections, :t_connections,
    ]),
    :line_reactor => Vector{Symbol}([
        :status, :f_bus, :t_bus, :f_connections, :t_connections,
    ]),
    :series_capacitor => Vector{Symbol}([
        :status, :f_bus, :t_bus, :f_connections, :t_connections,
    ]),
    :shunt => Vector{Symbol}([
        :status, :bus, :connections, :configuration, :gs, :gs, :vnom,
    ]),
    :shunt_capacitor => Vector{Symbol}([
        :status, :bus, :connections, :configuration, :bs, :vnom,
    ]),
    :shunt_reactor => Vector{Symbol}([
        :status, :bus, :connections, :configuration, :bs, :vnom,
    ]),
    :load => Vector{Symbol}([
        :status, :bus, :connections, :configuration, :model, :pd_nom, :qd_nom,
        :vnom,
    ]),
    :generator => Vector{Symbol}([
        :status, :bus, :connections, :configuration, :model,
    ]),
    :solar => Vector{Symbol}([
        :status, :bus, :connections, :configuration,
    ]),
    :storage => Vector{Symbol}([
        :status, :bus, :connections, :configuration, :energy,
        :charge_efficiency, :discharge_efficiency, :rs, :xs, :pex, :qex,
    ]),
    :voltage_source => Vector{Symbol}([
        :status, :bus, :connections, :configuration, :vm, :va,
    ]),
    :linecode => Vector{Symbol}([
        :rs, :xs, :g_fr, :g_to, :b_fr, :b_to, :cm_ub,
     ]),
    :xfmrcode => Vector{Symbol}([
        :status, :vnom, :snom, :xsc, :rs, :noloadloss, :imag, :tm_fix, :tm,
        :tm_min, :tm_max, :tm_step,
    ]),
    :grounding => Vector{Symbol}([]),

    # Future Components
    # :ev => Vector{Symbol}([]),
    # :wind => Vector{Symbol}([]),
    # :autotransformer => Vector{Symbol}([]),
    # :meter => Vector{Symbol}([])
)


"checks the engineering data model for correct data types, required fields and applies default checks"
function check_eng_data_model(data_eng::Dict{String,<:Any})
    for (component_type, components) in data_eng
        if isa(components, Dict)
            for (name, component) in components
                _check_eng_component_dtypes(data_eng, component_type, name)

                for field in get(_eng_model_req_fields, Symbol(component_type), Vector{Symbol}([]))
                    @assert haskey(component, string(field)) "The property \'$field\' is missing on $component_type $name"
                end

                check = get(_eng_model_checks, Symbol(component_type), missing)
                if !ismissing(check)
                    getfield(PowerModelsDistribution, check)(data_eng, name)
                end
            end
        end
    end
end


"checks that an engineering model component has the correct data types"
function _check_eng_component_dtypes(data_eng::Dict{String,<:Any}, component_type::String, component_name::Any; additional_dtypes=Dict{Symbol,Type}())
    if haskey(_eng_model_dtypes, Symbol(component_type))
        dtypes = merge(_eng_model_dtypes[Symbol(component_type)], additional_dtypes)
    else
        dtypes = additional_dtypes
    end

    if haskey(data_eng, component_type) && haskey(data_eng[component_type], component_name)
        component = data_eng[component_type][component_name]

        for (field, dtype) in dtypes
            if haskey(component, string(field))
                @assert isa(component[string(field)], dtype) "$component_type $component_name: the property $field should be a $dtype, not a $(typeof(component[string(field)]))"
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
function _check_has_size(component::Dict{String,<:Any}, fields::Vector{String}, data_size::Union{Int, Tuple}; context::Union{String,Missing}=missing, allow_missing::Bool=true)
    for key in fields
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
    @assert conf in [DELTA, WYE] "$context: the configuration should be \'delta\' or \'wye\', not \'$conf\'."

    return conf==WYE ? length(object["connections"])-1 : length(object["connections"])
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
    @assert model in [POWER, IMPEDANCE, CURRENT, EXPONENTIAL]

    if model==POWER
        _check_has_keys(load, ["pd", "qd"], context="load $name, $model:")
        _check_has_size(load, ["pd", "qd"], N, context="load $name, $model:")
    elseif model==EXPONENTIAL
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
    _check_same_size(data_eng["linecode"][name], string.([:rs, :xs, :g_fr, :g_to, :b_fr, :b_to]))
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
    _check_has_size(transformer, ["bus", "connections", "vnom", "snom", "configuration", "polarity", "rs", "tm_fix", "tm_set", "tm_lb", "tm_ub", "tm_step"], nrw, context="trans $name")

    @assert length(transformer["xsc"])==(nrw^2-nrw)/2

    nphs = []
    for w in 1:nrw
        @assert transformer["configuration"][w] in [WYE, DELTA]

        conf = transformer["configuration"][w]
        conns = transformer["connections"][w]
        nph = conf==WYE ? length(conns)-1 : length(conns)
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
    if config==WYE
        @assert length(shunt_capacitor["qd_ref"])==N-1 "capacitor $name: qd_ref should have $(N-1) elements."
    else
        @assert length(shunt_capacitor["qd_ref"])==N "capacitor $name: qd_ref should have $N elements."
    end

    @assert config in [DELTA, WYE, "wye-grounded", "wye-floating"]

    if config==DELTA
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
