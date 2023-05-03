"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssEdgeObject, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    add_bus!(eng, create_eng_object(UnsupportedEngEdgeObject, dss_obj; import_all=import_all))
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssNodeObject, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    add_bus!(eng, create_eng_object(UnsupportedEngNodeObject, dss_obj; import_all=import_all))
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssOptions, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    eng.settings.base_frequency = dss_obj.defaultbasefrequency

    if import_all
        eng.settings.dss = dss_obj
    end
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssBuscoords, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    if !haskey(eng.extra, "bus")
        eng.extra["bus"] = Dict{String,Any}()
    end

    if !haskey(eng.extra["bus"], dss_obj.bus)
        eng.extra["bus"][dss_obj.bus] = Dict{String,Any}()
    end

    eng.extra["bus"][dss_obj.bus]["x"] = dss_obj.x
    eng.extra["bus"][dss_obj.bus]["y"] = dss_obj.y
end


"""
"""
function add_bus!(eng::EngineeringDataModel, bus_id_string::String)::String
    bus_name, terminals = _parse_bus_id(bus_id_string)
    if !isempty(bus_name)
        phases = length(terminals)
        terminals = _get_conductors_ordered(bus_id_string; default=[collect(1:phases)..., 0], pad_ground=false)

        if !haskey(eng.bus, bus_name)
            eng.bus[bus_name] = EngBusObj(;
                name = bus_name,
                terminals = terminals,
            )
        else
            eng.bus[bus_name].terminals = sort(union(eng.bus[bus_name].terminals, terminals))
        end

        if 0 in terminals
            _register_awaiting_ground!(eng, eng.bus[bus_name], terminals)
        end
    end

    return bus_id_string
end


"""
"""
function add_bus!(eng::EngineeringDataModel, eng_obj::T)::T where T <: EngNodeObject
    if !haskey(eng.bus, eng_obj.bus)
        eng.bus[eng_obj.bus] = EngBusObj(;
            name = eng_obj.bus,
            terminals = eng_obj.connections,
        )
    else
        eng.bus[eng_obj.bus].terminals = sort(union(eng.bus[eng_obj.bus].terminals, eng_obj.connections))
    end

    if 0 in eng_obj.connections
        _register_awaiting_ground!(eng, eng.bus[eng_obj.bus], eng_obj.connections)
    end

    return eng_obj
end


"""
"""
function add_bus!(eng::EngineeringDataModel, eng_obj::T)::T where T <: EngEdgeObject
    if !haskey(eng.bus, eng_obj.f_bus)
        eng.bus[eng_obj.f_bus] = EngBusObj(;
            name = eng_obj.f_bus,
            terminals = eng_obj.f_connections,
        )
    else
        eng.bus[eng_obj.f_bus].terminals = sort(union(eng.bus[eng_obj.f_bus].terminals, eng_obj.f_connections))
    end

    if 0 in eng_obj.f_connections
        _register_awaiting_ground!(eng, eng.bus[eng_obj.f_bus], eng_obj.f_connections)
    end

    if !haskey(eng.bus, eng_obj.t_bus)
        eng.bus[eng_obj.t_bus] = EngBusObj(;
            name = eng_obj.t_bus,
            terminals = eng_obj.t_connections,
        )
    else
        eng.bus[eng_obj.t_bus].terminals = sort(union(eng.bus[eng_obj.t_bus].terminals, eng_obj.t_connections))
    end

    if 0 in eng_obj.t_connections
        _register_awaiting_ground!(eng, eng.bus[eng_obj.t_bus], eng_obj.t_connections)
    end

    return eng_obj
end


"""
"""
function add_bus!(eng::EngineeringDataModel, eng_obj::T)::T where T <: EngTransformer
    for (w, bus_id) in enumerate(eng_obj.bus)
        if !haskey(eng.bus, bus_id)
            eng.bus[bus_id] = EngBusObj(;
                name = bus_id,
                terminals = eng_obj.connections[w],
            )
        else
            eng.bus[bus_id].terminals = sort(union(eng.bus[bus_id].terminals, eng_obj.connections[w]))
        end

        if 0 in eng_obj.connections[w]
            _register_awaiting_ground!(eng, eng.bus[bus_id], eng_obj.connections[w])
        end
    end

    return eng_obj
end


"""
"""
function add_voltage_base_defaults!(eng::EngineeringDataModel)
    for (_,eng_obj) in eng.voltage_source
        eng.settings.vbases_default[eng_obj.bus] = eng_obj.vm[1]
    end
end


"""
"""
function add_bus_vbases!(eng::EngineeringDataModel)
    bus_vbases = calc_voltage_bases(eng)

    for (id,vbase) in bus_vbases
        eng.bus[id].vbase = vbase
    end
end


"""
"""
function update_time_series_references!(eng::EngineeringDataModel)
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssLoadshape, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    eng_obj_p, eng_obj_q = create_eng_object(EngTimeSeriesObj, dss_obj; import_all=import_all)

    if eng_obj_p.values == eng_obj_q.values
        eng_obj_p.name = "$(dss_obj.name)"

        eng.time_series[eng_obj_p.name] = eng_obj_p
    else
        eng.time_series[eng_obj_p.name] = eng_obj_p
        eng.time_series[eng_obj_q.name] = eng_obj_q
    end
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssLoad, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    eng.load[dss_obj.name] = add_bus!(eng, create_eng_object(EngLoadObj, dss_obj; import_all=import_all, time_series=time_series))
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssCapacitor, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    bus1_name = _parse_bus_id(dss_obj.bus1)[1]
    bus2_name = _parse_bus_id(dss_obj.bus2)[1]

    if bus1_name == bus2_name
        eng.shunt[dss_obj.name] = add_bus!(eng, create_eng_object(EngShuntObj, dss_obj; import_all=import_all))
    else
        @info "capacitors as constant impedance elements is not supported, treating capacitor.$(dss_obj.name) like line"
        eng.line[dss_obj.name] = add_bus!(eng, create_eng_object(EngLineObj, dss_obj; import_all=import_all))
    end
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssReactor, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    if isempty(dss_obj.bus2)
        eng.shunt[dss_obj.name] = add_bus!(eng, create_eng_object(EngShuntObj, dss_obj; import_all=import_all))
    else
        @info "reactors as constant impedance elements is not explicitly supported, treating reactor.$(dss_obj.name) like line"
        eng.line[dss_obj.name] = add_bus!(eng, create_eng_object(EngLineObj, dss_obj; import_all=import_all))
    end
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssGenerator, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    eng.generator[dss_obj.name] = add_bus!(eng, create_eng_object(EngGeneratorObj, dss_obj; import_all=import_all, time_series=time_series))
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssVsource, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    eng.voltage_source[dss_obj.name] = add_bus!(eng, create_eng_object(EngVoltageSourceObj, dss_obj; import_all=import_all))

    if dss_obj.name == "source"
        eng.settings.sbase_default = dss_obj.basemva
    end

    if !isempty(dss_obj.bus2) && _parse_bus_id(dss_obj.bus2)[1] != _parse_bus_id(dss_obj.bus1)[1]
        add_bus!(eng, dss_obj.bus2)
    end

    return eng.voltage_source[dss_obj.name]
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssLinecode, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    eng.linecode[dss_obj.name] = create_eng_object(EngLinecodeObj, dss_obj; import_all=import_all)
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssLine, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    eng_obj = create_eng_object(dss_obj.switch ? EngSwitchObj : EngLineObj, dss_obj; import_all=import_all, dss=dss)

    if isa(eng_obj, EngSwitch)
        eng.switch[eng_obj.name] = add_bus!(eng, eng_obj)
    else
        eng.line[eng_obj.name] = add_bus!(eng, eng_obj)
    end
end

# const _dss_to_eng_types::Dict{Type,Type} = Dict{Type,Type}(
#     DssTransformer => EngTransformer
# )

# function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::T; import_all::Bool=false) where T <: DssObject
#     obj_type = string(T)
#     eng_type = _dss_to_eng_types[T]
#     eng_pn = Symbol(lowercase(string(eng_type)[4:end]))

#     eng[eng_pn][dss_obj.name] = add_bus!(eng, create_eng_object())
# end

"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssXfmrcode, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    eng.xfmrcode[dss_obj.name] = create_eng_object(EngXfmrcodeObj, dss_obj; import_all=import_all)
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssTransformer, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    eng.transformer[dss_obj.name] = add_bus!(eng, create_eng_object(EngTransformerObj, dss_obj; import_all=import_all, dss=dss))
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssPvsystem, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    eng.solar[dss_obj.name] = add_bus!(eng, create_eng_object(EngSolarObj, dss_obj; import_all=import_all, time_series=time_series))
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssStorage, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    eng.storage[dss_obj.name] = add_bus!(eng, create_eng_object(EngStorageObj, dss_obj; import_all=import_all, time_series=time_series))
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssRegcontrol, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    add_controls!(eng, dss_obj.transformer, create_eng_object(EngTransformerControlsObj, dss_obj; import_all=import_all))
end


"""
"""
function add_controls!(eng::EngineeringDataModel, transformer_id::String, eng_obj::EngTransformerControls)
    if haskey(eng.transformer, transformer_id)
        _correct_ptphase!(eng_obj, eng.transformer[transformer_id])
        if ismissing(eng.transformer[transformer_id].controls)
            eng.transformer[transformer_id].controls = eng_obj
        else
            merge!(eng.transformer[transformer_id].controls, eng_obj)
        end
    else
        @warn "there is no transformer object '$(transformer_id)', cannot add controls"
    end
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssCapcontrol, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily")
    add_controls!(eng, dss_obj.capacitor, create_eng_object(EngShuntControlsObj, dss_obj; import_all=import_all))
end


"""
"""
function add_controls!(eng::EngineeringDataModel, capacitor_id::String, eng_obj::EngShuntControls)
    if haskey(eng.shunt, capacitor_id)
        if ismissing(eng.shunt[capacitor_id].controls)
            eng.shunt[capacitor_id].controls = eng_obj
        else
            merge!(eng.shunt[capacitor_id].controls, eng_obj)
        end
    else
        @warn "there is no shunt object '$(capacitor_id)', cannot add controls"
    end
end

convert_dss2eng!(::EngineeringModel, ::DssDataObject, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily") = nothing
convert_dss2eng!(::EngineeringModel, ::DssControlObject, dss::Union{Missing,DssModel}=missing, import_all::Bool=false, time_series::String="daily") = nothing


"""
"""
function transform_data_model(
    ::Type{EngineeringModel{NetworkModel}},
    dss::DssModel;
    import_all::Bool=false,
    time_series::String="daily",
    bank_transformers::Bool=true,
    dss2eng_extensions::Vector{<:Function}=Function[],
    transformations::Vector{<:Any}=[],
    disable_dss_voltage_bounds::Bool=false,
    kwargs...
    )::EngineeringDataModel

    eng = EngineeringDataModel()

    for (_, dss_objects) in dss
        if !isa(dss_objects, DssObject)
            for (_, dss_obj) in dss_objects
                convert_dss2eng!(eng, dss_obj, dss, import_all, time_series)
            end
        else
            convert_dss2eng!(eng, dss_objects, dss, import_all, time_series)
        end
    end

    update_time_series_references!(eng)

    add_voltage_base_defaults!(eng)

    discover_terminals!(eng)

    bank_transformers && bank_transformers!(eng)

    add_bus_vbases!(eng)

    !disable_dss_voltage_bounds && add_bus_bounds_from_dss!(eng, dss)

    # user extensions
    for dss2eng_func! in dss2eng_extensions
        dss2eng_func!(eng.extra, eng, dss)
    end

    for transform! in transformations
        @assert isa(transform!, Function) || isa(transform!, Tuple{<:Function,Vararg{Pair{String,<:Any}}})

        if isa(transform!, Tuple)
            transform![1](eng; [Symbol(k)=>v for (k,v) in transform![2:end]]...)
        else
            transform!(eng)
        end
    end

    return eng
end


transform_data_model(::Type{EngineeringModel{U}}, dss::DssRawModel; kwargs...) where U = transform_data_model(EngineeringModel{U}, transform_data_model(DssModel, dss); kwargs...)
transform_data_model(dss::DssModel; kwargs...) = transform_data_model(EngineeringModel{NetworkModel}, dss; kwargs...)


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


"""
"""
function add_bus_bounds_from_dss!(eng::EngineeringModel, dss::OpenDssDataModel)
    for dss_type in ["load", "generator", "storage", "pvsystem"]
        for (id,dss_obj) in getproperty(dss, Symbol(dss_type))
            bus_id = getproperty(eng, Symbol(dss_type == "pvsystem" ? "solar" : dss_type))[id].bus
            bus = eng.bus[bus_id]
            bus_terminals = bus.terminals
            bus_grounded = bus.grounded
            bus_vbase = bus.vbase
            obj_connections = getproperty(eng, Symbol(dss_type == "pvsystem" ? "solar" : dss_type))[id].connections

            if !ismissing(bus_vbase)
                _vm_lb = dss_obj.vminpu * bus_vbase
                _vm_ub = dss_obj.vmaxpu * bus_vbase

                vm_lb = ismissing(bus.vm_lb) ? fill(0.0, length(bus_terminals)) : bus.vm_lb
                vm_ub = ismissing(bus.vm_ub) ? fill(Inf, length(bus_terminals)) : bus.vm_ub

                for (idx,t) in enumerate(bus_terminals)
                    if t ∈ obj_connections && t ∉ bus_grounded
                        if vm_lb[idx] < _vm_lb
                            vm_lb[idx] = _vm_lb
                        end

                        if vm_ub[idx] > _vm_ub
                            vm_ub[idx] = _vm_ub
                        end
                    end
                end

                eng.bus[bus_id].vm_lb = vm_lb
                eng.bus[bus_id].vm_ub = vm_ub
            end
        end
    end
end
