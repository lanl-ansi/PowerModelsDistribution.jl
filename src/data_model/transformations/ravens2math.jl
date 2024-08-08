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

    find_conductor_ids!(data_math)
    _map_conductor_ids!(data_math)
    _map_settings_vbases_default!(data_math)

end


"""
Converts ravens connectivity_node components into mathematical bus components.
"""
function _map_ravens2math_connectivity_node!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    voltage_scale_factor_sqrt3 = data_math["settings"]["voltage_scale_factor"] * sqrt(3)
    connectivity_nodes = get(data_ravens, "ConnectivityNode", Dict{String,Any}())

    for (name, ravens_obj) in connectivity_nodes
        index = length(data_math["bus"]) + 1
        math_obj = _init_math_obj_ravens("bus", name, ravens_obj, index; pass_props=pass_props)

        # Set basic bus properties
        math_obj["bus_i"] = index
        math_obj["source_id"] = "bus.$name"
        math_obj["bus_type"] = 1  # Default bus_type, will be modified as needed
        math_obj["grounded"] = Bool[0, 0, 0]
        math_obj["vm_pair_lb"] = Tuple{Any, Any, Real}[]
        math_obj["vm_pair_ub"] = Tuple{Any, Any, Real}[]

        # Set voltage magnitude and angle
        if haskey(ravens_obj, "SvVoltage.v")
            math_obj["vm"] = (ravens_obj["SvVoltage.v"] / voltage_scale_factor_sqrt3)
        end

        if haskey(ravens_obj, "SvVoltage.angle")
            math_obj["va"] = ravens_obj["SvVoltage.angle"]
        end

        # Store the mathematical bus object
        data_math["bus"]["$(index)"] = math_obj

        # Update bus lookup if necessary
        data_math["bus_lookup"] = get(data_math, "bus_lookup", Dict{Any,Int}())
        data_math["bus_lookup"][name] = index

        # Map the ravens object to the math object
        push!(data_math["map"], Dict(
            "from" => name,
            "to" => "bus.$index",
            "unmap_function" => "_map_math2ravens_bus!"
        ))
    end
end


"""
Converts ravens conductors (e.g., ACLineSegments) into mathematical branches.
"""
function _map_ravens2math_conductor!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    conductors = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]["Conductor"]

    for (name, ravens_obj) in get(conductors, "ACLineSegment", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj_ravens("ac_line_segment", name, ravens_obj, length(data_math["branch"]) + 1; pass_props=pass_props)

        nphases = length(ravens_obj["ACLineSegment.ACLineSegmentPhase"])
        terminals = ravens_obj["ConductingEquipment.Terminals"]

        f_node = replace(split(terminals[1]["Terminal.ConnectivityNode"], "::")[2], "'" => "")
        t_node = replace(split(terminals[2]["Terminal.ConnectivityNode"], "::")[2], "'" => "")

        math_obj["f_bus"] = data_math["bus_lookup"][f_node]
        math_obj["t_bus"] = data_math["bus_lookup"][t_node]

        phase_map = Dict("SinglePhaseKind.A" => 1, "SinglePhaseKind.B" => 2, "SinglePhaseKind.C" => 3)
        bus_terminals = nphases >= 3 ? collect(1:nphases) : [phase_map[phase["ACLineSegmentPhase.phase"]] for phase in ravens_obj["ACLineSegment.ACLineSegmentPhase"]]

        for bus in [math_obj["f_bus"], math_obj["t_bus"]]
            data_math["bus"][string(bus)]["terminals"] = bus_terminals
            data_math["bus"][string(bus)]["vmin"] = fill(0.0, nphases)
            data_math["bus"][string(bus)]["vmax"] = fill(Inf, nphases)
        end

        math_obj["f_connections"] = bus_terminals
        math_obj["t_connections"] = bus_terminals

        impedance_name = replace(split(ravens_obj["ACLineSegment.PerLengthImpedance"], "::")[2], "'" => "")
        impedance_data = data_ravens["PerLengthLineParameter"]["PerLengthImpedance"]["PerLengthPhaseImpedance"][impedance_name]

        math_obj["br_r"] = _impedance_conversion_ravens(impedance_data, ravens_obj, "PhaseImpedanceData.r")
        math_obj["br_x"] = _impedance_conversion_ravens(impedance_data, ravens_obj, "PhaseImpedanceData.x")

        base_freq = data_math["settings"]["base_frequency"]
        for (key, param) in [("b_fr", "PhaseImpedanceData.b"), ("b_to", "PhaseImpedanceData.b"), ("g_fr", "PhaseImpedanceData.g"), ("g_to", "PhaseImpedanceData.g")]
            math_obj[key] = _admittance_conversion_ravens(impedance_data, ravens_obj, param; freq=base_freq)
        end

        math_obj["angmin"] = get(ravens_obj, "vad_lb", fill(-60.0, nphases))
        math_obj["angmax"] = get(ravens_obj, "vad_ub", fill(60.0, nphases))

        oplimitset_id = replace(split(terminals[1]["ACDCTerminal.OperationalLimitSet"], "::")[2], "'" => "")
        oplimitset = data_ravens["OperationalLimitSet"][oplimitset_id]["OperationalLimitSet.OperationalLimitValue"][2]

        limit_keys = [("CurrentLimit.value", "c_rating_a"), ("CurrentLimit.value", "c_rating_b"), ("CurrentLimit.value", "c_rating_c"),
                      ("ApparentPowerLimit.value", "rate_a"), ("ApparentPowerLimit.value", "rate_b"), ("ApparentPowerLimit.value", "rate_c")]

        for (f_key, t_key) in limit_keys
            math_obj[t_key] = haskey(oplimitset, f_key) ? fill(oplimitset[f_key], nphases) : fill(Inf, nphases)
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
    end
end


"""
Converts ravens load components into mathematical load components.
"""
function _map_ravens2math_energy_consumer!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    conducting_equipment = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]["EnergyConnection"]
    power_scale_factor = data_math["settings"]["power_scale_factor"]
    voltage_scale_factor = data_math["settings"]["voltage_scale_factor"]
    voltage_scale_factor_sqrt3 = voltage_scale_factor * sqrt(3)

    for (name, ravens_obj) in get(conducting_equipment, "EnergyConsumer", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj_ravens("energy_consumer", name, ravens_obj, length(data_math["load"]) + 1; pass_props=pass_props)

        # Set the load bus based on connectivity node
        connectivity_node = replace(split(ravens_obj["ConductingEquipment.Terminals"][1]["Terminal.ConnectivityNode"], "::")[2], "'" => "")
        math_obj["load_bus"] = data_math["bus_lookup"][connectivity_node]

        # Handle Load Response Characteristics
        load_response_characts = replace(split(ravens_obj["EnergyConsumer.LoadResponseCharacteristic"], "::")[2], "'" => "")
        if load_response_characts == "Constant Z"
            math_obj["model"] = IMPEDANCE
        elseif load_response_characts == "Motor"
            @error("Load model not supported yet!")
        elseif load_response_characts == "Mix Motor/Res"
            @error("Load model not supported yet!")
        elseif load_response_characts == "Constant I"
            math_obj["model"] = CURRENT
        elseif load_response_characts == "Variable P, Fixed Q"
            @error("Load model not supported yet!")
        elseif load_response_characts == "Variable P, Fixed X"
            @error("Load model not supported yet!")
        else
            if load_response_characts != "Constant kVA"
                @warn("Load model (response characteristic) for $(name) not supported! Defaulting to 'Constant kVA'")
            end
            # Set default model and consumption values
            math_obj["model"] = POWER
        end

        # Set p and q
        math_obj["pd"] = [ravens_obj["EnergyConsumer.p"] / power_scale_factor]
        math_obj["qd"] = [ravens_obj["EnergyConsumer.q"] / power_scale_factor]

        # Set the nominal voltage
        base_voltage_ref = replace(split(ravens_obj["ConductingEquipment.BaseVoltage"], "::")[2], "'" => "")
        base_voltage = data_ravens["BaseVoltage"][base_voltage_ref]["BaseVoltage.nominalVoltage"]
        math_obj["vnom_kv"] = (base_voltage / voltage_scale_factor) / (sqrt(3) / 2)

        # Set voltage bounds for the bus connected
        bus_info = string(math_obj["load_bus"])
        bus_conn = data_math["bus"][bus_info]
        op_limit_id = replace(split(data_ravens["ConnectivityNode"][connectivity_node]["ConnectivityNode.OperationalLimitSet"], "::")[2], "'" => "")
        op_limits = data_ravens["OperationalLimitSet"][op_limit_id]["OperationalLimitSet.OperationalLimitValue"]
        op_limit_max = op_limits[1]["VoltageLimit.value"] / voltage_scale_factor_sqrt3
        op_limit_min = op_limits[2]["VoltageLimit.value"] / voltage_scale_factor_sqrt3

        # Handle phase-specific or three-phase connection
        phase_map = Dict("SinglePhaseKind.A" => 1, "SinglePhaseKind.B" => 2, "SinglePhaseKind.C" => 3)
        if haskey(ravens_obj, "EnergyConsumer.EnergyConsumerPhase")
            connections = Vector{Int64}()
            for phase_info in ravens_obj["EnergyConsumer.EnergyConsumerPhase"]
                phase = phase_map[phase_info["EnergyConsumerPhase.phase"]]
                phase_index = findfirst(==(phase), bus_conn["terminals"])
                bus_conn["vmax"][phase_index] = op_limit_max
                bus_conn["vmin"][phase_index] = op_limit_min
                push!(connections, phase)
            end
            math_obj["connections"] = connections
        else
            N = length(bus_conn["terminals"])
            bus_conn["vmax"] = fill(op_limit_max, N)
            bus_conn["vmin"] = fill(op_limit_min, N)
            math_obj["connections"] = bus_conn["terminals"]
        end

        # Set the configuration
        # TODO: ADD: "PhaseShuntConnectionKind.Yn", "PhaseShuntConnectionKind.I", "PhaseShuntConnectionKind.G"
        config_map = Dict("PhaseShuntConnectionKind.Y" => WYE, "PhaseShuntConnectionKind.D" => DELTA)
        config = get(config_map, ravens_obj["EnergyConsumer.phaseConnection"], nothing)
        if config !== nothing
            math_obj["configuration"] = config
        else
            @error("Configuration of load $(name) is not supported.")
        end

        # Set status, dispatchable flag, and index
        math_obj["status"] = ravens_obj["ConductingEquipment.SvStatus"]
        math_obj["dispatchable"] = 0
        data_math["load"]["$(math_obj["index"])"] = math_obj

        # Handle grounding
        if ravens_obj["EnergyConsumer.grounded"] == "true"
            bus_conn["grounded"] = Bool[0, 0, 0]
        end

        # Set bus type to PQ bus
        bus_conn["bus_type"] = 1

        # Map the object
        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "load.$(math_obj["index"])",
            "unmap_function" => "_map_math2ravens_load!",
        ))
    end
end


"""
Converts ravens voltage sources into mathematical generators and (if needed) impedance branches to represent the loss model.
"""
function _map_ravens2math_energy_source!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    conducting_equipment = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]["EnergyConnection"]
    voltage_scale_factor = data_math["settings"]["voltage_scale_factor"]
    voltage_scale_factor_sqrt3 = voltage_scale_factor * sqrt(3)

    for (name, ravens_obj) in get(conducting_equipment, "EnergySource", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj_ravens("energy_source", name, ravens_obj, length(data_math["gen"]) + 1; pass_props=pass_props)
        math_obj["name"] = "_virtual_gen.energy_source.$name"

        # Get connectivity node info (bus info)
        connectivity_node = replace(split(ravens_obj["ConductingEquipment.Terminals"][1]["Terminal.ConnectivityNode"], "::")[2], "'" => "")
        gen_bus = data_math["bus_lookup"][connectivity_node]
        math_obj["gen_bus"] = gen_bus
        bus_conn = data_math["bus"][string(gen_bus)]
        bus_conn["bus_type"] = 3  # Set bus type to PV bus

        # Handle phase-specific or three-phase connection
        connections = Vector{Int64}()
        nconductors = length(get(ravens_obj, "EnergySource.EnergySourcePhase", bus_conn["terminals"]))

        if haskey(ravens_obj, "EnergySource.EnergySourcePhase")
            phase_map = Dict("SinglePhaseKind.A" => 1, "SinglePhaseKind.B" => 2, "SinglePhaseKind.C" => 3)
            for phase_info in ravens_obj["EnergySource.EnergySourcePhase"]
                phase = phase_map[phase_info["EnergySourcePhase.phase"]]
                push!(connections, phase)
            end
            math_obj["connections"] = connections
        else
            math_obj["connections"] = bus_conn["terminals"]
        end

        # Generator status and configuration
        math_obj["gen_status"] = Int(ravens_obj["ConductingEquipment.SvStatus"])
        math_obj["configuration"] = get(ravens_obj, "EnergySource.connectionKind", WYE)

        # Vnom and vbases_default
        base_voltage_ref = replace(split(ravens_obj["ConductingEquipment.BaseVoltage"], "::")[2], "'" => "")
        vnom = data_ravens["BaseVoltage"][base_voltage_ref]["BaseVoltage.nominalVoltage"] / sqrt(nconductors)
        data_math["settings"]["vbases_default"][connectivity_node] = vnom / voltage_scale_factor

        # Power, voltage, and limits
        nphases = nconductors  # You can adjust nphases based on your specific kron reduction logic if needed
        fill_values = (v) -> fill(v, nphases)
        math_obj["pg"] = get(ravens_obj, "EnergySource.activePower", fill_values(0.0))
        math_obj["qg"] = get(ravens_obj, "EnergySource.reactivePower", fill_values(0.0))
        math_obj["vg"] = fill(get(ravens_obj, "EnergySource.voltageMagnitude", voltage_scale_factor_sqrt3) / voltage_scale_factor_sqrt3, nphases)
        math_obj["pmin"] = get(ravens_obj, "EnergySource.pMin", fill_values(-Inf))
        math_obj["pmax"] = get(ravens_obj, "EnergySource.pMax", fill_values(Inf))
        math_obj["qmin"] = get(ravens_obj, "EnergySource.qMin", fill_values(-Inf))
        math_obj["qmax"] = get(ravens_obj, "EnergySource.qMax", fill_values(Inf))

        # Control mode and source ID
        math_obj["control_mode"] = Int(get(ravens_obj, "EnergySource.connectionKind", ISOCHRONOUS))
        math_obj["source_id"] = "energy_source.$name"

        # Add generator cost model
        _add_gen_cost_model!(math_obj, ravens_obj)

        # Check for impedance and adjust bus type if necessary
        map_to = "gen.$(math_obj["index"])"
        if !all(isapprox.(get(ravens_obj, "EnergySource.r", zeros(1, 1)), 0)) && !all(isapprox.(get(ravens_obj, "EnergySource.x", zeros(1, 1)), 0))
            bus_conn["bus_type"] = 1  # Virtual bus becomes the new slack bus

            bus_obj = Dict(
                "bus_i" => length(data_math["bus"]) + 1,
                "index" => length(data_math["bus"]) + 1,
                "terminals" => math_obj["connections"],
                "grounded" => Bool[0, 0, 0],
                "name" => "_virtual_bus.energy_source.$name",
                "bus_type" => math_obj["gen_status"] == 0 ? 4 : math_obj["control_mode"] == Int(ISOCHRONOUS) ? 3 : 2,
                "vm" => fill(ravens_obj["EnergySource.voltageMagnitude"] / voltage_scale_factor_sqrt3, nphases),
                "va" => rad2deg.(_wrap_to_pi.([-2 * π / nphases * (i - 1) + deg2rad(ravens_obj["EnergySource.voltageAngle"]) for i in 1:nphases])),
                "vmin" => fill(ravens_obj["EnergySource.voltageMagnitude"] / voltage_scale_factor_sqrt3, nphases),
                "vmax" => fill(ravens_obj["EnergySource.voltageMagnitude"] / voltage_scale_factor_sqrt3, nphases),
                "vm_pair_lb" => deepcopy(get(ravens_obj, "EnergySource.vpairMin", Tuple{Any,Any,Real}[])),
                "vm_pair_ub" => deepcopy(get(ravens_obj, "EnergySource.vpairMax", Tuple{Any,Any,Real}[])),
                "source_id" => "energy_source.$name",
            )

            math_obj["gen_bus"] = bus_obj["bus_i"]
            data_math["bus"]["$(bus_obj["index"])"] = bus_obj

            branch_obj = Dict(
                "name" => "_virtual_branch.energy_source.$name",
                "source_id" => "energy_source.$name",
                "f_bus" => bus_obj["bus_i"],
                "t_bus" => gen_bus,
                "f_connections" => math_obj["connections"],
                "t_connections" => math_obj["connections"],
                "angmin" => fill(-10.0, nconductors),
                "angmax" => fill(10.0, nconductors),
                "c_rating_a" => fill(Inf, nconductors),
                "br_status" => math_obj["gen_status"],
                "br_r" => _impedance_conversion_ravens_energy_source(data_ravens, ravens_obj, "EnergySource.r", "EnergySource.r0"),
                "br_x" => _impedance_conversion_ravens_energy_source(data_ravens, ravens_obj, "EnergySource.x", "EnergySource.x0"),
                "g_fr" => zeros(nconductors, nconductors),
                "g_to" => zeros(nconductors, nconductors),
                "b_fr" => zeros(nconductors, nconductors),
                "b_to" => zeros(nconductors, nconductors),
                "index" => length(data_math["branch"]) + 1
            )

            data_math["branch"]["$(branch_obj["index"])"] = branch_obj
            map_to = [map_to, "bus.$(bus_obj["index"])", "branch.$(branch_obj["index"])"]
        else
            # Handle bus voltage limits if no impedance is present
            vm_lb = math_obj["control_mode"] == Int(ISOCHRONOUS) ? ravens_obj["EnergySource.voltageMagnitude"] : get(ravens_obj, "EnergySource.vMin", fill(0.0, nphases))
            vm_ub = math_obj["control_mode"] == Int(ISOCHRONOUS) ? ravens_obj["EnergySource.voltageMagnitude"] : get(ravens_obj, "EnergySource.vMax", fill(Inf, nphases))

            data_math["bus"]["$gen_bus"]["vmin"] = [vm_lb..., fill(0.0, nconductors - nphases)...]
            data_math["bus"]["$gen_bus"]["vmax"] = [vm_ub..., fill(Inf, nconductors - nphases)...]
            data_math["bus"]["$gen_bus"]["vm"] = fill(ravens_obj["EnergySource.voltageMagnitude"] / voltage_scale_factor_sqrt3, nphases)
            data_math["bus"]["$gen_bus"]["va"] = fill(ravens_obj["EnergySource.voltageAngle"], nphases)

            data_math["bus"]["$gen_bus"]["bus_type"] = _compute_bus_type(bus_conn["bus_type"], math_obj["gen_status"], math_obj["control_mode"])
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
