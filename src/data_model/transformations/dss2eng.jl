"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssOptions; import_all::Bool=false)
    eng.settings.base_frequency = dss_obj.defaultbasefrequency

    if import_all
        eng.settings.dss = dss_obj
    end
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssBuscoords; import_all::Bool=false)
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
function add_bus!(eng::EngineeringDataModel, eng_obj::T)::T where T <: EngNodeObject
    if !haskey(eng.bus, eng_obj.bus)
        eng.bus[eng_obj.bus] = EngBus(;
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
        eng.bus[eng_obj.f_bus] = EngBus(;
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
        eng.bus[eng_obj.t_bus] = EngBus(;
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
            eng.bus[bus_id] = EngBus(;
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
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssLoadshape; import_all::Bool=false)
    eng_obj_p, eng_obj_q = create_eng_object(EngTimeSeries, dss_obj; import_all=import_all)

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
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssLoad; import_all::Bool=false, time_series::String="daily")
    eng.load[dss_obj.name] = add_bus!(eng, create_eng_object(EngLoad, dss_obj; import_all=import_all, time_series=time_series))
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssCapacitor; import_all::Bool=false)
    bus1_name = _parse_bus_id(dss_obj.bus1)[1]
    bus2_name = _parse_bus_id(dss_obj.bus2)[1]

    if bus1_name == bus2_name
        eng.shunt[dss_obj.name] = add_bus!(eng, create_eng_object(EngShunt, dss_obj; import_all=import_all))
    else
        @info "capacitors as constant impedance elements is not supported, treating capacitor.$(dss_obj.name) like line"
        eng.line[dss_obj.name] = add_bus!(eng, create_eng_object(EngLine, dss_obj; import_all=import_all))
    end
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssReactor; import_all::Bool=false)
    if isempty(dss_obj.bus2)
        eng.shunt[dss_obj.name] = add_bus!(eng, create_eng_object(EngShunt, dss_obj; import_all=import_all))
    else
        @info "reactors as constant impedance elements is not explicitly supported, treating reactor.$(dss_obj.name) like line"
        eng.line[dss_obj.name] = add_bus!(eng, create_eng_object(EngLine, dss_obj; import_all=import_all))
    end
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssGenerator; import_all::Bool=false, time_series::String="daily")
    eng.generator[dss_obj.name] = add_bus!(eng, create_eng_object(EngGenerator, dss_obj; import_all=import_all, time_series=time_series))
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssVsource; import_all::Bool=false)
    eng.voltage_source[dss_obj.name] = add_bus!(eng, create_eng_object(EngVoltageSource, dss_obj; import_all=import_all))

    if dss_obj.name == "source"
        eng.settings.sbase_default = dss_obj.basemva
    end

    return eng.voltage_source[dss_obj.name]
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssLinecode; import_all::Bool=false)
    eng.linecode[dss_obj.name] = create_eng_object(EngLinecode, dss_obj; import_all=import_all)
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssLine; import_all::Bool=false)
    eng_obj = create_eng_object(dss_obj.switch ? EngSwitch : EngLine, dss_obj; import_all=import_all)

    if isa(eng_obj, EngSwitch)
        eng.switch[eng_obj.name] = add_bus!(eng, eng_obj)
    else
        eng.line[eng_obj.name] = add_bus!(eng, eng_obj)
    end
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssXfmrcode; import_all::Bool=false)
    eng.xfmrcode[dss_obj.name] = create_eng_object(EngXfmrcode, dss_obj; import_all=import_all)
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssTransformer; import_all::Bool=false, dss::Union{Missing,OpenDssDataModel}=missing)
    eng.transformer[dss_obj.name] = add_bus!(eng, create_eng_object(EngTransformer, dss_obj; import_all=import_all, dss=dss))
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssPvsystem; import_all::Bool=false, time_series::String="daily")
    eng.solar[dss_obj.name] = add_bus!(eng, create_eng_object(EngSolar, dss_obj; import_all=import_all, time_series=time_series))
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssStorage; import_all::Bool=false, time_series::String="daily")
    eng.storage[dss_obj.name] = add_bus!(eng, create_eng_object(EngStorage, dss_obj; import_all=import_all, time_series=time_series))
end


"""
"""
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssRegcontrol; import_all::Bool=false)
    add_controls!(eng, dss_obj.transformer, create_eng_object(EngTransformerControls, dss_obj; import_all=import_all))
end


"""
"""
function add_controls!(eng::EngineeringDataModel, transformer_id::String, eng_obj::EngTransformerControls)
    if haskey(eng.transformer, transformer_id)
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
function convert_dss2eng!(eng::EngineeringDataModel, dss_obj::DssCapcontrol; import_all::Bool=false)
    add_controls!(eng, dss_obj.capacitor, create_eng_object(EngShuntControls, dss_obj; import_all=import_all))
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

convert_dss2eng!(::EngineeringModel, ::DssDataObject; import_all::Bool=false) = nothing
convert_dss2eng!(::EngineeringModel, ::DssControlObject; import_all::Bool=false) = nothing


"""
"""
function transform_data_model(
    ::Type{EngineeringModel},
    dss::OpenDssDataModel;
    import_all::Bool=false,
    time_series::String="daily",
    bank_transformers::Bool=true,
    dss2eng_extensions::Vector{Function}=Function[],
    )::EngineeringDataModel

    eng = EngineeringDataModel()

    for (_, dss_objects) in dss
        if !isa(dss_objects, OpenDssObject)
            for (_, dss_obj) in dss_objects
                if isa(dss_obj, DssTimeSeriesObjects)
                    convert_dss2eng!(eng, dss_obj; import_all=import_all, time_series=time_series)
                elseif isa(dss_obj, DssTransformer)
                    convert_dss2eng!(eng, dss_obj; import_all=import_all, dss=dss)
                else
                    convert_dss2eng!(eng, dss_obj; import_all=import_all)
                end
            end
        else
            convert_dss2eng!(eng, dss_objects; import_all=import_all)
        end
    end

    update_time_series_references!(eng)

    add_voltage_base_defaults!(eng)

    discover_terminals!(eng)

    bank_transformers && bank_transformers!(eng)

    add_bus_vbases!(eng)

    # user extensions
    for func! in dss2eng_extensions
        func!(eng.extra, eng, dss)
    end

    return eng
end


transform_data_model(::EngineeringModel, dss::DssRawModel; kwargs...) = transform_data_model(EngineeringModel, transform_data_model(DssModel, dss); kwargs...)

parse_file(file::String; data_model=EngineeringModel) = transform_data_model(EngineeringModel, parse_dss(file))
parse_file(io::IO; data_model=EngineeringModel) = transform_data_model(EngineeringModel, parse_dss(io))


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
