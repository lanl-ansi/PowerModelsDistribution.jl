"cim-ravens to math object mapping"
const _math_to_ravens = Dict{String,String}(
    "bus" => "connectivity_node",
    "transformer" => "power_transformer",
    "switch" => "switch",
    # "shunt" => "shunt_compensator",
    "load" => "energy_consumer",
    # "generator" => "rotating_machine",
    # "solar" => "photovoltaic_unit",
    # "storage" => "battery_unit",
    "voltage_source" => "energy_source",
)


"list of nodal type elements in the ravens model"
const _ravens_node_elements = String[
    "energy_consumer", "energy_source"
]

"list of edge type elements in the ravens model"
const _ravens_edge_elements = String[
    "conductor"
]

"list of all ravens asset types"
const pmd_ravens_asset_types = String[
    "connectivity_node", _ravens_edge_elements..., _ravens_node_elements...
]


function transform_data_model_ravens(
    data::Dict{String,<:Any};
    kron_reduce::Bool=true,
    phase_project::Bool=false,
    multinetwork::Bool=false,
    global_keys::Set{String}=Set{String}(),
    ravens2math_passthrough::Dict{String,<:Vector{<:String}}=Dict{String,Vector{String}}(),
    ravens2math_extensions::Vector{<:Function}=Function[],
    make_pu::Bool=true,
    make_pu_extensions::Vector{<:Function}=Function[],
    correct_network_data::Bool=true,
    )::Dict{String,Any}

    current_data_model = get(data, "data_model", MATHEMATICAL)

    ## TODO
    # if multinetwork && !ismultinetwork(data)
    #     data = make_multinetwork(data; global_keys=global_keys)
    # end

    data_math = _map_ravens2math(
        data;
        kron_reduce=kron_reduce,
        phase_project=phase_project,
        ravens2math_extensions=ravens2math_extensions,
        ravens2math_passthrough=ravens2math_passthrough,
        global_keys=global_keys,
    )

    # TODO: Correct network data transforms a lot of the values of lines/branches (other values maybe too)
    correct_network_data && correct_network_data!(data_math; make_pu=make_pu, make_pu_extensions=make_pu_extensions)

    return data_math

end


"base function for converting ravens model to mathematical model"
function _map_ravens2math(
    data_ravens::Dict{String,<:Any};
    kron_reduce::Bool=true,
    phase_project::Bool=false,
    ravens2math_extensions::Vector{<:Function}=Function[],
    ravens2math_passthrough::Dict{String,Vector{String}}=Dict{String,Vector{String}}(),
    global_keys::Set{String}=Set{String}(),
    )::Dict{String,Any}

    _data_ravens = deepcopy(data_ravens)

    ## TODO
    # if kron_reduce
    #     apply_kron_reduction!(_data_ravens)
    #     apply_phase_projection_delta!(_data_ravens)
    # end

    # if phase_project
    #     apply_phase_projection!(_data_ravens)
    # end


    # # TODO: Add vbases to settings
    # _vbases_default = Dict("vbases_default" => Dict{String, Real}())
    # for (name, ravens_obj) in get(data_ravens, "BaseVoltage", Dict{String,Any}())
    #     _vbases_default = Dict(name => ravens_obj["BaseVoltage.nominalVoltage"])
    # end

    # TODO: Add settings (defaults)
    basemva = 100
    _settings = Dict("sbase_default" => basemva * 1e3,
                    "voltage_scale_factor" => 1e3,
                    "power_scale_factor" => 1e3,
                    "base_frequency" => get(_data_ravens, "BaseFrequency", 60.0),
                    "vbases_default" => Dict{String,Real}()
    )

    ## TODO: Multinetwork
    if ismultinetwork(data_ravens)
        data_math = Dict{String,Any}(
            "name" => get(_data_ravens, "name", ""),
            "data_model" => MATHEMATICAL,
            "nw" => Dict{String,Any}(
                n => Dict{String,Any}(
                    "per_unit" => get(_data_ravens, "per_unit", false),
                    "is_projected" => get(nw, "is_projected", false),
                    "is_kron_reduced" => get(nw, "is_kron_reduced", false),
                    "settings" => deepcopy(_settings),
                    "time_elapsed" => get(nw, "time_elapsed", 1.0),
                ) for (n,nw) in _data_ravens["nw"]
            ),
            "multinetwork" => ismultinetwork(data_ravens),
            [k => data_ravens[k] for k in global_keys if haskey(data_ravens, k)]...
        )
    else
        data_math = Dict{String,Any}(
            "name" => get(_data_ravens, "name", ""),
            "per_unit" => get(_data_ravens, "per_unit", false),
            "data_model" => MATHEMATICAL,
            "is_projected" => get(_data_ravens, "is_projected", false),
            "is_kron_reduced" => get(_data_ravens, "is_kron_reduced", false),
            "settings" => deepcopy(_settings),
            "time_elapsed" => get(_data_ravens, "time_elapsed", 1.0),
        )
    end

    apply_pmd!(_map_ravens2math_nw!, data_math, _data_ravens; ravens2math_passthrough=ravens2math_passthrough, ravens2math_extensions=ravens2math_extensions)

    if ismultinetwork(data_ravens)
        _collect_nw_maps!(data_math)
        _collect_nw_bus_lookups!(data_math)
    end

    return data_math
end


"""
"""
function _map_ravens2math_nw!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; ravens2math_passthrough::Dict{String,Vector{String}}=Dict{String,Vector{String}}(), ravens2math_extensions::Vector{<:Function}=Function[])
    data_math["map"] = Vector{Dict{String,Any}}([
        Dict{String,Any}("unmap_function" => "_map_math2ravens_root!")
    ])

    _init_base_components!(data_math)

    ## TODO
    # for property in get(ravens2math_passthrough, "root", String[])
    #     if haskey(data_ravens, property)
    #         data_math[property] = deepcopy(data_ravens[property])
    #     end
    # end

    for type in pmd_ravens_asset_types
        getfield(PowerModelsDistribution, Symbol("_map_ravens2math_$(type)!"))(data_math, data_ravens; pass_props=get(ravens2math_passthrough, type, String[]))
    end

    # Custom ravens2math transformation functions
    for ravens2math_func! in ravens2math_extensions
        ravens2math_func!(data_math, data_ravens)
    end

    ## TODO: See if this is still neccesary for RAVENS
    # # post fix
    # if !get(data_math, "is_kron_reduced", false)
    #     #TODO fix this in place / throw error instead? IEEE8500 leads to switches
    #     # with 3x3 R matrices but only 1 phase
    #     #NOTE: Don't do this when kron-reducing, it will undo the padding
    #     _slice_branches!(data_math)
    # end

    find_conductor_ids!(data_math)
    _map_conductor_ids!(data_math)

    # TODO: See if this is still neccesary for RAVENS
    _map_settings_vbases_default!(data_math)

end


"converts ravens connectivity_node components into mathematical bus components"
function _map_ravens2math_connectivity_node!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[])

    _voltage_scale_factor_sqrt3 = data_math["settings"]["voltage_scale_factor"]*sqrt(3)

    for (name, ravens_obj) in get(data_ravens, "ConnectivityNode", Dict{String,Any}())

        math_obj = _init_math_obj_ravens("bus", name, ravens_obj, length(data_math["bus"])+1; pass_props=pass_props)
        math_obj["bus_i"] = math_obj["index"]
        math_obj["source_id"] = "bus.$name"

        # default bus_type, change as elements are added (e.g., load, generator, source).
        math_obj["bus_type"] = 1

        # TODO: needed? - default grounded
        math_obj["grounded"] = Bool[0, 0, 0]

        ## TODO
        # # take care of grounding; convert to shunt if lossy
        # grounded_perfect, shunts = _convert_grounding(ravens_obj["terminals"], ravens_obj["grounded"], ravens_obj["rg"], ravens_obj["xg"])

        # math_obj["grounded"] = grounded_perfect
        # to_sh = []
        # for (sh_connections, sh_y) in shunts
        #     sh_index = length(data_math["shunt"]) + 1
        #     data_math["shunt"]["$sh_index"] = Dict(
        #         "index" => sh_index,
        #         "shunt_bus" => math_obj["bus_i"],
        #         "connections" => sh_connections,
        #         "gs" => real.(sh_y),
        #         "bs" => imag.(sh_y),
        #     )
        #     push!(to_sh, "shunt.$sh_index")
        # end

        if haskey(ravens_obj, "SvVoltage.v")
            math_obj["vm"] = (ravens_obj["SvVoltage.v"]/_voltage_scale_factor_sqrt3)
        end

        if haskey(ravens_obj, "SvVoltage.angle")
            math_obj["va"] = ravens_obj["SvVoltage.angle"]
        end

        # TODO: add vm_pair_lb/ub (may not be needed!)
        math_obj["vm_pair_lb"] = Tuple{Any, Any, Real}[]
        math_obj["vm_pair_ub"] = Tuple{Any, Any, Real}[]

        data_math["bus"]["$(math_obj["index"])"] = math_obj

        if !haskey(data_math, "bus_lookup")
            data_math["bus_lookup"] = Dict{Any,Int}()
        end

        data_math["bus_lookup"][name] = math_obj["index"]

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "bus.$(math_obj["index"])",
            "unmap_function" => "_map_math2ravens_bus!",
        ))
    end
end


"converts ravensineering conductors (e.g., ACLineSegments) into mathematical branches"
function _map_ravens2math_conductor!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[])

    _data_ravens_conductor = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]["Conductor"]

    for (name, ravens_obj) in get(_data_ravens_conductor, "ACLineSegment", Dict{Any,Dict{String,Any}}())

        math_obj = _init_math_obj_ravens("ac_line_segment", name, ravens_obj, length(data_math["branch"])+1; pass_props=pass_props)

        nphases = length(ravens_obj["ACLineSegment.ACLineSegmentPhase"])

        f_connectivity_node = replace(split(ravens_obj["ConductingEquipment.Terminals"][1]["Terminal.ConnectivityNode"], "::")[2], "'" => "")
        t_connectivity_node = replace(split(ravens_obj["ConductingEquipment.Terminals"][2]["Terminal.ConnectivityNode"], "::")[2], "'" => "")
        math_obj["f_bus"] = data_math["bus_lookup"][f_connectivity_node]
        math_obj["t_bus"] = data_math["bus_lookup"][t_connectivity_node]

        # define the number of terminals (i.e., phases) that the bus/es have, based on the lines
        if nphases >= 3
            data_math["bus"][string(math_obj["f_bus"])]["terminals"]= [i for i = 1:nphases]
            data_math["bus"][string(math_obj["t_bus"])]["terminals"]= [i for i = 1:nphases]
        else
            _terminals = Vector{Int64}()
            for phase_info in ravens_obj["ACLineSegment.ACLineSegmentPhase"]
                phase = phase_info["ACLineSegmentPhase.phase"]
                if phase == "SinglePhaseKind.A"
                    push!(_terminals, 1)
                elseif phase == "SinglePhaseKind.B"
                    push!(_terminals, 2)
                elseif phase == "SinglePhaseKind.C"
                    push!(_terminals, 3)
                else
                    @error("Terminals/Phases for buses '$(f_connectivity_node)' and '$(t_connectivity_node)' not recognized. Check your model!")
                end
            end
            data_math["bus"][string(math_obj["f_bus"])]["terminals"]= _terminals
            data_math["bus"][string(math_obj["t_bus"])]["terminals"]= _terminals
        end

        # Initilize vmin and vmax in buses
        data_math["bus"][string(math_obj["f_bus"])]["vmin"] = fill(0.0, nphases)
        data_math["bus"][string(math_obj["f_bus"])]["vmax"] = fill(Inf, nphases)

        data_math["bus"][string(math_obj["t_bus"])]["vmin"] = fill(0.0, nphases)
        data_math["bus"][string(math_obj["t_bus"])]["vmax"] = fill(Inf, nphases)

        #-------------------------------------------------------------------------------

        # f_ and t_connections for the lines/branches
        math_obj["f_connections"] = data_math["bus"][string(math_obj["f_bus"])]["terminals"]
        math_obj["t_connections"] = data_math["bus"][string(math_obj["t_bus"])]["terminals"]


        _perlengthimpedance_name = replace(split(ravens_obj["ACLineSegment.PerLengthImpedance"], "::")[2], "'" => "")
        _perlengthimpedance_data = data_ravens["PerLengthLineParameter"]["PerLengthImpedance"]["PerLengthPhaseImpedance"][_perlengthimpedance_name]

        math_obj["br_r"] = _impedance_conversion_ravens(_perlengthimpedance_data, ravens_obj, "PhaseImpedanceData.r")
        math_obj["br_x"] = _impedance_conversion_ravens(_perlengthimpedance_data, ravens_obj, "PhaseImpedanceData.x")

        # b is given in mhos in CIM-RAVENS schema
        math_obj["b_fr"] = _admittance_conversion_ravens(_perlengthimpedance_data, ravens_obj, "PhaseImpedanceData.b"; freq=data_math["settings"]["base_frequency"])
        math_obj["b_to"] = _admittance_conversion_ravens(_perlengthimpedance_data, ravens_obj, "PhaseImpedanceData.b"; freq=data_math["settings"]["base_frequency"])

        # g is given in mhos in CIM-RAVENS schema
        math_obj["g_fr"] = _admittance_conversion_ravens(_perlengthimpedance_data, ravens_obj, "PhaseImpedanceData.g"; freq=data_math["settings"]["base_frequency"])
        math_obj["g_to"] = _admittance_conversion_ravens(_perlengthimpedance_data, ravens_obj, "PhaseImpedanceData.g"; freq=data_math["settings"]["base_frequency"])

        ## TODO:: these do not appear in JSON schema, are they needed? if yes, we need to represent them on CIM-RAVENS schema
        math_obj["angmin"] = get(ravens_obj, "vad_lb", fill(-60.0, nphases))
        math_obj["angmax"] = get(ravens_obj, "vad_ub", fill( 60.0, nphases))

        # Operational Limits
        _oplimitset_id = replace(split(ravens_obj["ConductingEquipment.Terminals"][1]["ACDCTerminal.OperationalLimitSet"], "::")[2], "'" => "")
        _oplimitset = data_ravens["OperationalLimitSet"][_oplimitset_id]["OperationalLimitSet.OperationalLimitValue"][2] # [1] is Normal amps, [2] is Emergency amps - by default use emerg amps

        for (f_key, t_key) in [("CurrentLimit.value", "c_rating_a"), ("CurrentLimit.value", "c_rating_b"), ("CurrentLimit.value", "c_rating_c"),
            ("ApparentPowerLimit.value", "rate_a"), ("ApparentPowerLimit.value", "rate_b"), ("ApparentPowerLimit.value", "rate_c")]
            math_obj[t_key] = haskey(_oplimitset, f_key) ? fill(_oplimitset[f_key], nphases) : fill(Inf, nphases)
        end

        math_obj["br_status"] = get(ravens_obj, "ConductingEquipment.SvStatus", 1)

        data_math["branch"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "branch.$(math_obj["index"])",
            "unmap_function" => "_map_math2ravens_line!",
        ))
    end

end

# TODO: Transformers need a lot of changes/refactors!!!
"converts ravensineering n-winding transformers into mathematical ideal 2-winding lossless transformer branches and impedance branches to represent the loss model"
function _map_ravens2math_power_transformer!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[])

    _data_ravens_transformer = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]

    for (name, ravens_obj) in get(_data_ravens_transformer, "PowerTransformer", Dict{Any,Dict{String,Any}}())

        # Build map first, so we can update it as we decompose the transformer
        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => String[],
            "unmap_function" => "_map_math2ravens_transformer!",
        ))

        to_map = data_math["map"][end]["to"]

        # _apply_xfmrcode!(ravens_obj, data_ravens)

        # PowerTransformerTank
        if haskey(ravens_obj, "PowerTransformer.TransformerTank")

            @info "PowerTransformerTank ENTER"
            @assert all(haskey(ravens_obj, k) for k in ["f_bus", "t_bus", "f_connections", "t_connections"]) "Incomplete definition of AL2W tranformer $name, aborting ravens2math conversion"

            nphases = length(ravens_obj["f_connections"])

            math_obj = Dict{String,Any}(
                "name" => name,
                "source_id" => ravens_obj["source_id"],
                "f_bus" => data_math["bus_lookup"][ravens_obj["f_bus"]],
                "t_bus" => data_math["bus_lookup"][ravens_obj["t_bus"]],
                "f_connections" => ravens_obj["f_connections"],
                "t_connections" => ravens_obj["t_connections"],
                "configuration" => get(ravens_obj, "configuration", WYE),
                "tm_nom" => get(ravens_obj, "tm_nom", 1.0),
                "tm_set" => get(ravens_obj, "tm_set", fill(1.0, nphases)),
                "tm_fix" => get(ravens_obj, "tm_fix", fill(true, nphases)),
                "polarity" => get(ravens_obj, "polarity", -1),
                "sm_ub" => get(ravens_obj, "sm_ub", Inf),
                "cm_ub" => get(ravens_obj, "cm_ub", Inf),
                "status" => Int(get(ravens_obj, "status", ENABLED)),
                "index" => length(data_math["transformer"])+1
            )

            for k in [["tm_lb", "tm_ub"]; pass_props]
                if haskey(ravens_obj, k)
                    math_obj[k] = ravens_obj[k]
                end
            end

            data_math["transformer"]["$(math_obj["index"])"] = math_obj

            push!(to_map, "transformer.$(math_obj["index"])")


        # PowerTransformerEnd
        else

            nrw = length(ravens_obj["bus"])

            for w in 1:length(ravens_obj["PowerTransformer.PowerTransformerEnd"])

                vnom = ravens_obj["PowerTransformer.PowerTransformerEnd"][w]["PowerTransformerEnd.ratedU"]
                snom = ravens_obj["PowerTransformer.PowerTransformerEnd"][w]["PowerTransformerEnd.ratedS"]

                # calculate zbase in which the data is specified, and convert to SI
                zbase = (vnom.^2) ./ snom

                # x_sc is specified with respect to first winding
                x_sc = ravens_obj["xsc"] .* zbase[1]

                # rs is specified with respect to each winding
                r_s = ravens_obj["rw"] .* zbase

                g_sh =  (ravens_obj["noloadloss"]*snom[1])/vnom[1]^2
                b_sh = -(ravens_obj["cmag"]*snom[1])/vnom[1]^2

                # data is measured externally, but we now refer it to the internal side
                ratios = vnom/data_ravens["settings"]["voltage_scale_factor"]
                x_sc = x_sc./ratios[1]^2
                r_s = r_s./ratios.^2
                g_sh = g_sh*ratios[1]^2
                b_sh = b_sh*ratios[1]^2

                # convert x_sc from list of upper triangle elements to an explicit dict
                y_sh = g_sh + im*b_sh
                z_sc = Dict([(key, im*x_sc[i]) for (i,key) in enumerate([(i,j) for i in 1:nrw for j in i+1:nrw])])

                dims = length(ravens_obj["tm_set"][1])
                transformer_t_bus_w = _build_loss_model!(data_math, name, to_map, r_s, z_sc, y_sh,ravens_obj["connections"][1]; nphases=dims, status=Int(ravens_obj["status"] == ENABLED))

                # 2-WINDING TRANSFORMER
                # make virtual bus and mark it for reduction
                tm_nom = ravens_obj["configuration"][w]==DELTA ? ravens_obj["vm_nom"][w]*sqrt(3) : ravens_obj["vm_nom"][w]
                transformer_2wa_obj = Dict{String,Any}(
                    "name"          => "_virtual_transformer.$name.$w",
                    "source_id"     => "_virtual_transformer.$(ravens_obj["source_id"]).$w",
                    "f_bus"         => data_math["bus_lookup"][ravens_obj["bus"][w]],
                    "t_bus"         => transformer_t_bus_w[w],
                    "tm_nom"        => tm_nom,
                    "f_connections" => ravens_obj["connections"][w],
                    "t_connections" => get(data_math, "is_kron_reduced", false) ? ravens_obj["connections"][1] : collect(1:dims+1),
                    "configuration" => ravens_obj["configuration"][w],
                    "polarity"      => ravens_obj["polarity"][w],
                    "tm_set"        => ravens_obj["tm_set"][w],
                    "tm_fix"        => ravens_obj["tm_fix"][w],
                    "sm_ub"         => get(ravens_obj, "sm_ub", Inf),
                    "cm_ub"         => get(ravens_obj, "cm_ub", Inf),
                    "status"        => ravens_obj["status"] == DISABLED ? 0 : 1,
                    "index"         => length(data_math["transformer"])+1
                )

                for prop in [["tm_lb", "tm_ub", "tm_step"]; pass_props]
                    if haskey(ravens_obj, prop)
                        transformer_2wa_obj[prop] = ravens_obj[prop][w]
                    end
                end

                data_math["transformer"]["$(transformer_2wa_obj["index"])"] = transformer_2wa_obj

                # add regcontrol items to math model
                if haskey(ravens_obj,"controls") && !all(data_math["transformer"]["$(transformer_2wa_obj["index"])"]["tm_fix"])
                    reg_obj = Dict{String,Any}(
                        "vreg" => ravens_obj["controls"]["vreg"][w],
                        "band" => ravens_obj["controls"]["band"][w],
                        "ptratio" => ravens_obj["controls"]["ptratio"][w],
                        "ctprim" => ravens_obj["controls"]["ctprim"][w],
                        "r" => ravens_obj["controls"]["r"][w],
                        "x" => ravens_obj["controls"]["x"][w],
                    )
                    data_math["transformer"]["$(transformer_2wa_obj["index"])"]["controls"] = reg_obj
                end

                if w==3 && ravens_obj["polarity"][w]==-1 # identify center-tapped transformer and mark all secondary-side nodes as triplex by adding va_start
                    default_va = [0, -120, 120][ravens_obj["connections"][1][1]]
                    data_math["bus"]["$(transformer_2wa_obj["f_bus"])"]["va_start"] = haskey(data_ravens["bus"][ravens_obj["bus"][w]],"va_start") ? data_ravens["bus"][ravens_obj["bus"][w]]["va_start"] : [default_va, (default_va+180)]
                    idx = 0
                    bus_ids = []
                    t_bus = haskey(data_ravens, "line") ? [data["t_bus"] for (_,data) in data_ravens["line"] if data["f_bus"] == ravens_obj["bus"][w]] : []
                    while length(t_bus)>0 || idx<length(bus_ids)
                        for bus_idx in t_bus
                            bus_id = data_math["bus_lookup"]["$bus_idx"]
                            push!(bus_ids, bus_id)
                            default_va = [0, -120, 120][ravens_obj["connections"][1][1]]
                            data_math["bus"]["$bus_id"]["va_start"] = haskey(data_ravens["bus"]["$bus_idx"],"va_start") ? data_ravens["bus"]["$bus_idx"]["va_start"] : [default_va, (default_va+180)]
                        end
                        idx += 1
                        t_bus = [data["t_bus"] for (_,data) in data_ravens["line"] if data["f_bus"] == data_math["bus"]["$(bus_ids[idx])"]["name"]]
                    end
                end

                push!(to_map, "transformer.$(transformer_2wa_obj["index"])")

            end
        end
    end
end



"converts ravensineering load components into mathematical load components"
function _map_ravens2math_energy_consumer!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[])

    _data_ravens_energyconnection = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]["EnergyConnection"]
    _power_scale_factor = data_math["settings"]["power_scale_factor"]
    _voltage_scale_factor = data_math["settings"]["voltage_scale_factor"]

    for (name, ravens_obj) in get(_data_ravens_energyconnection, "EnergyConsumer", Dict{Any,Dict{String,Any}}())

        # Get load response characteristic
        load_response_characts = replace(split(ravens_obj["EnergyConsumer.LoadResponseCharacteristic"], "::")[2], "'" => "")

        # Initialize energy consumer object
        math_obj = _init_math_obj_ravens("energy_consumer", name, ravens_obj, length(data_math["load"])+1; pass_props=pass_props)

        # Get Bus/ConnectivityNode
        connectivity_node = replace(split(ravens_obj["ConductingEquipment.Terminals"][1]["Terminal.ConnectivityNode"], "::")[2], "'" => "")
        math_obj["load_bus"] = data_math["bus_lookup"][connectivity_node]

        # Model of the Load/EnergyConsumer
        if load_response_characts == "Constant Z"

        elseif load_response_characts == "Motor"

        elseif load_response_characts == "Mix Motor/Res"

        elseif load_response_characts == "Constant I"

        elseif load_response_characts == "Variable P, Fixed Q"

        elseif load_response_characts == "Variable P, Fixed X"

        # "Constant kVA" : default
        else

            if load_response_characts != "Constant kVA"
                @warn("Load model (response characteristic) for $(name) not supported! Defaulting to 'Constant kVA'")
            end

            # Add the model of the load based on loadresponsecharacteristic
            math_obj["model"] = POWER

            # P and Q consumption
            math_obj["pd"] = [ravens_obj["EnergyConsumer.p"]/_power_scale_factor]
            math_obj["qd"] = [ravens_obj["EnergyConsumer.q"]/_power_scale_factor]

        end

        # Vnom
        base_voltage_ref = replace(split(ravens_obj["ConductingEquipment.BaseVoltage"], "::")[2], "'" => "")
        math_obj["vnom_kv"] = (data_ravens["BaseVoltage"][base_voltage_ref]["BaseVoltage.nominalVoltage"]/_voltage_scale_factor)/(sqrt(3)/2)

        # Get voltage bounds for specific bus connected (TODO: see if it can be coverted to standalone function to avoid repetition)
        bus_info = string(math_obj["load_bus"])
        bus_conn = data_math["bus"][bus_info]

        base_voltage_id = replace(split(ravens_obj["ConductingEquipment.BaseVoltage"], "::")[2], "'" => "")
        base_voltage = data_ravens["BaseVoltage"][base_voltage_id]["BaseVoltage.nominalVoltage"]
        op_limit_id = replace(split(data_ravens["ConnectivityNode"][connectivity_node]["ConnectivityNode.OperationalLimitSet"], "::")[2], "'" => "")
        # op_limit_max = data_ravens["OperationalLimitSet"][op_limit_id]["OperationalLimitSet.OperationalLimitValue"][1]["VoltageLimit.value"]./base_voltage
        # op_limit_min = data_ravens["OperationalLimitSet"][op_limit_id]["OperationalLimitSet.OperationalLimitValue"][2]["VoltageLimit.value"]./base_voltage
        op_limit_max = data_ravens["OperationalLimitSet"][op_limit_id]["OperationalLimitSet.OperationalLimitValue"][1]["VoltageLimit.value"]
        op_limit_min = data_ravens["OperationalLimitSet"][op_limit_id]["OperationalLimitSet.OperationalLimitValue"][2]["VoltageLimit.value"]

        _connections = Vector{Int64}()
        if haskey(ravens_obj, "EnergyConsumer.EnergyConsumerPhase")
            for phase_info in ravens_obj["EnergyConsumer.EnergyConsumerPhase"]
                phase = phase_info["EnergyConsumerPhase.phase"]
                if phase == "SinglePhaseKind.A"
                    phase_index = findfirst(==(1), bus_conn["terminals"])
                    bus_conn["vmax"][phase_index] = op_limit_max
                    bus_conn["vmin"][phase_index] = op_limit_min
                    push!(_connections, 1)
                elseif phase == "SinglePhaseKind.B"
                    phase_index = findfirst(==(2), bus_conn["terminals"])
                    bus_conn["vmax"][phase_index] = op_limit_max
                    bus_conn["vmin"][phase_index] = op_limit_min
                    push!(_connections, 2)
                elseif phase == "SinglePhaseKind.C"
                    phase_index = findfirst(==(3), bus_conn["terminals"])
                    bus_conn["vmax"][phase_index] = op_limit_max
                    bus_conn["vmin"][phase_index] = op_limit_min
                    push!(_connections, 3)
                else
                    @error("Terminals/Phases for buses '$(f_connectivity_node)' and '$(t_connectivity_node)' not recognized. Check your model!")
                end
            end
            math_obj["connections"] = _connections
        else
            # assumes it is a three-phase connection
            N = length(bus_conn["terminals"])
            bus_conn["vmax"] = fill(op_limit_max, N)
            bus_conn["vmin"] = fill(op_limit_min, N)
            math_obj["connections"] = bus_conn["terminals"]
        end

        # TODO: Configuration
        _config = ravens_obj["EnergyConsumer.phaseConnection"]
        if _config == "PhaseShuntConnectionKind.Y"
            math_obj["configuration"] = WYE
        elseif _config == "PhaseShuntConnectionKind.D"
            math_obj["configuration"] = DELTA
        elseif _config == "PhaseShuntConnectionKind.Yn"
        elseif _config == "PhaseShuntConnectionKind.I"
        elseif _config == "PhaseShuntConnectionKind.G"
        else
            @error("Configuration of load $(name) is not supported.")
        end

        # Status
        math_obj["status"] = ravens_obj["ConductingEquipment.SvStatus"]

        # TODO: Dispatchable? default=0
        math_obj["dispatchable"] = 0

        # Index
        data_math["load"]["$(math_obj["index"])"] = math_obj

        # TODO: Assign 'grounding' to corresponding connectivity node
        if ravens_obj["EnergyConsumer.grounded"] == "true"
            bus_conn["grounded"] = Bool[0, 0, 0]
        end

        # Revise bus_type of connectivity node to PQ bus
        bus_conn["bus_type"] = 1

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "load.$(math_obj["index"])",
            "unmap_function" => "_map_math2ravens_load!",
        ))

   end
end


#TODO:
"converts ravensineering voltage sources into mathematical generators and (if needed) impedance branches to represent the loss model"
function _map_ravens2math_energy_source!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[])

    _data_ravens_energyconnection = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]["EnergyConnection"]
    _voltage_scale_factor = data_math["settings"]["voltage_scale_factor"]
    _voltage_scale_factor_sqrt3 = data_math["settings"]["voltage_scale_factor"]*sqrt(3)

    for (name, ravens_obj) in get(_data_ravens_energyconnection, "EnergySource", Dict{String,Any}())

        math_obj = _init_math_obj_ravens("energy_source", name, ravens_obj, length(data_math["gen"])+1; pass_props=pass_props)
        math_obj["name"] = "_virtual_gen.energy_source.$name"

        # get connectivity node info. (bus info)
        connectivity_node = replace(split(ravens_obj["ConductingEquipment.Terminals"][1]["Terminal.ConnectivityNode"], "::")[2], "'" => "")
        math_obj["gen_bus"] = gen_bus = data_math["bus_lookup"][connectivity_node]
        bus_conn = data_math["bus"][string(gen_bus)]

        # Revise bus_type of connectivity node to PV bus
        bus_conn["bus_type"] = 3

        nconductors = 0
        _connections = Vector{Int64}()
        if haskey(ravens_obj, "EnergySource.EnergySourcePhase")

            nconductors = length(ravens_obj["EnergySource.EnergySourcePhase"])
            # TODO: how to redefine the kron reduction?
            # nphases = get(ravens_obj, "configuration", WYE) == WYE && !get(data_ravens, "is_kron_reduced", false) ? nconductors - 1 : nconductors
            nphases = nconductors

            for phase_info in ravens_obj["EnergySource.EnergySourcePhase"]
                phase = phase_info["EnergySourcePhase.phase"]
                if phase == "SinglePhaseKind.A"
                    push!(_connections, 1)
                elseif phase == "SinglePhaseKind.B"
                    push!(_connections, 2)
                elseif phase == "SinglePhaseKind.C"
                    push!(_connections, 3)
                else
                    @error("Connections for energy source '$(name)' not recognized. Check your model!")
                end
            end

            math_obj["connections"] = _connections
        else
            # assumes it is a three-phase connection
            nconductors = length(bus_conn["terminals"])
            # TODO: how to redefine the kron reduction? we need
            # nphases = get(ravens_obj, "configuration", WYE) == WYE && !get(data_ravens, "is_kron_reduced", false) ? nconductors - 1 : nconductors
            nphases = nconductors
            math_obj["connections"] = bus_conn["terminals"]
        end

        math_obj["gen_status"] = status = Int(ravens_obj["ConductingEquipment.SvStatus"])

        # Vnom and add vbases_default from energy source
        base_voltage_ref = replace(split(ravens_obj["ConductingEquipment.BaseVoltage"], "::")[2], "'" => "")
        vnom = data_ravens["BaseVoltage"][base_voltage_ref]["BaseVoltage.nominalVoltage"] / sqrt(nphases)
        data_math["settings"]["vbases_default"][connectivity_node] = vnom/_voltage_scale_factor

        # P, Q, Vg, etc.
        math_obj["pg"] = get(ravens_obj, "EnergySource.activePower", fill(0.0, nphases)).*fill(1.0, nphases)
        math_obj["qg"] = get(ravens_obj, "EnergySource.reactivePower", fill(0.0, nphases)).*fill(1.0, nphases)
        math_obj["vg"] = fill(get(ravens_obj, "EnergySource.voltageMagnitude", _voltage_scale_factor_sqrt3)/_voltage_scale_factor_sqrt3, nphases)
        math_obj["pmin"] = get(ravens_obj, "EnergySource.pMin", fill(-Inf, nphases)).*fill(1.0, nphases)
        math_obj["pmax"] = get(ravens_obj, "EnergySource.pMax", fill( Inf, nphases)).*fill(1.0, nphases)
        math_obj["qmin"] = get(ravens_obj, "EnergySource.qMin", fill(-Inf, nphases)).*fill(1.0, nphases)
        math_obj["qmax"] = get(ravens_obj, "EnergySource.qMax", fill( Inf, nphases)).*fill(1.0, nphases)

        # configuration
        math_obj["configuration"] = get(ravens_obj, "EnergySource.connectionKind", WYE)

        # TODO: Do we need a control_mode parameter for this? add to RAVENS schema
        math_obj["control_mode"] = control_mode = Int(get(ravens_obj, "EnergySource.connectionKind", ISOCHRONOUS))
        math_obj["source_id"] = "energy_source.$name"

        # TODO: inside this function, there are elements that do not exist in the RAVENS schema, so the default is used.
        _add_gen_cost_model!(math_obj, ravens_obj)

        map_to = "gen.$(math_obj["index"])"

        if !all(isapprox.(get(ravens_obj, "EnergySource.r", zeros(1, 1)), 0)) && !all(isapprox.(get(ravens_obj, "EnergySource.x", zeros(1, 1)), 0))

            # Revise bus_type of connectivity node to PV bus (virtual bus becomes the new slack bus)
            bus_conn["bus_type"] = 1

            f_bus = deepcopy(data_math["bus"]["$(math_obj["gen_bus"])"])

            bus_obj = Dict{String,Any}(
                "bus_i" => length(data_math["bus"])+1,
                "index" => length(data_math["bus"])+1,
                "terminals" => math_obj["connections"],
                # TODO: grounded energyosurce default?
                "grounded" => Bool[0, 0, 0],
                "name" => "_virtual_bus.energy_source.$name",
                "bus_type" => status == 0 ? 4 : control_mode == Int(ISOCHRONOUS) ? 3 : 2,
                "vm" => fill(ravens_obj["EnergySource.voltageMagnitude"]/_voltage_scale_factor_sqrt3, nphases),
                "va" => rad2deg.(_wrap_to_pi.([-2*pi/nphases*(i-1)+deg2rad(ravens_obj["EnergySource.voltageAngle"]) for i in 1:nphases])),
                "vmin" => fill(ravens_obj["EnergySource.voltageMagnitude"]/_voltage_scale_factor_sqrt3, nphases),
                "vmax" => fill(ravens_obj["EnergySource.voltageMagnitude"]/_voltage_scale_factor_sqrt3, nphases),
                "vm_pair_lb" => deepcopy(get(ravens_obj, "EnergySource.vpairMin", Tuple{Any,Any,Real}[])),
                "vm_pair_ub" => deepcopy(get(ravens_obj, "EnergySource.vpairMax", Tuple{Any,Any,Real}[])),
                "source_id" => "energy_source.$name",
            )

            # TODO:
            # for (i,t) in enumerate(math_obj["connections"])
            #     if data_math["bus"]["$(data_math["bus_lookup"][ravens_obj["bus"]])"]["grounded"][i]
            #         bus_obj["vm"][i] = 0
            #         bus_obj["vmin"][i] = 0
            #         bus_obj["vmax"][i] = Inf
            #     end
            # end

            math_obj["gen_bus"] = gen_bus = bus_obj["bus_i"]

            data_math["bus"]["$(bus_obj["index"])"] = bus_obj

            branch_obj = Dict{String,Any}(
                "name" => "_virtual_branch.energy_source.$name",
                "source_id" => "energy_source.$name",
                "f_bus" => bus_obj["bus_i"],
                "t_bus" => data_math["bus_lookup"][connectivity_node],
                "f_connections" => math_obj["connections"],
                "t_connections" => math_obj["connections"],
                "angmin" => fill(-10.0, nconductors),
                "angmax" => fill( 10.0, nconductors),
                "c_rating_a" => fill(Inf, nconductors),
                "br_status" => status,
                "br_r" => _impedance_conversion_ravens_energy_source(data_ravens, ravens_obj, "EnergySource.r", "EnergySource.r0"),
                "br_x" => _impedance_conversion_ravens_energy_source(data_ravens, ravens_obj, "EnergySource.x", "EnergySource.x0"),
                "g_fr" => zeros(nconductors, nconductors),
                "g_to" => zeros(nconductors, nconductors),
                "b_fr" => zeros(nconductors, nconductors),
                "b_to" => zeros(nconductors, nconductors),
                "index" => length(data_math["branch"])+1
            )

            data_math["branch"]["$(branch_obj["index"])"] = branch_obj

            map_to = [map_to, "bus.$(bus_obj["index"])", "branch.$(branch_obj["index"])"]
        else
            vm_lb = control_mode == Int(ISOCHRONOUS) ? ravens_obj["EnergySource.voltageMagnitude"] : get(ravens_obj, "EnergySource.vMin", fill(0.0, nphases))
            vm_ub = control_mode == Int(ISOCHRONOUS) ? ravens_obj["EnergySource.voltageMagnitude"] : get(ravens_obj, "EnergySource.vMax", fill(Inf, nphases))

            data_math["bus"]["$gen_bus"]["vmin"] = [vm_lb..., [0.0 for n in 1:(nconductors-nphases)]...]
            data_math["bus"]["$gen_bus"]["vmax"] = [vm_ub..., [Inf for n in 1:(nconductors-nphases)]...]
            data_math["bus"]["$gen_bus"]["vm"] = [ravens_obj["EnergySource.voltageMagnitude"]/_voltage_scale_factor_sqrt3..., [0.0 for n in 1:(nconductors-nphases)]...]
            data_math["bus"]["$gen_bus"]["va"] = [ravens_obj["EnergySource.voltageAngle"]..., [0.0 for n in 1:(nconductors-nphases)]...]

            bus_type = data_math["bus"]["$gen_bus"]["bus_type"]
            data_math["bus"]["$gen_bus"]["bus_type"] = _compute_bus_type(bus_type, status, control_mode)
        end

        data_math["gen"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => map_to,
            "unmap_function" => "_map_math2ravens_energy_source!",
        ))
    end
end


"converts ravensineering switches into mathematical switches and (if neeed) impedance branches to represent loss model"
function _map_ravens2math_switch!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    # TODO enable real switches (right now only using vitual lines)
    for (name, ravens_obj) in get(data_ravens, "switch", Dict{Any,Dict{String,Any}}())
      # TODO
    end
end
