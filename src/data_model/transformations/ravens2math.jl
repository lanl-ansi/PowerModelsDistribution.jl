"cim-ravens to math object mapping"
const _math_to_ravens = Dict{String,String}(
    "bus" => "connectivity_node",
    "transformer" => "power_transformer",
    "switch" => "switch",
    "shunt" => "shunt_compensator",
    "load" => "energy_consumer",
    "generator" => "rotating_machine",
    "solar" => "photovoltaic_unit",
    "storage" => "battery_unit",
    "voltage_source" => "energy_source",
)


"list of nodal type elements in the ravens model"
const _ravens_node_elements = String[
    "energy_consumer", "shunt_compensator", "rotating_machine", "power_electronics", "energy_source"
]

"list of edge type elements in the ravens model"
const _ravens_edge_elements = String[
    "conductor", "switch", "power_transformer"
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

    data_math = _map_ravens2math(
        data;
        multinetwork=multinetwork,
        kron_reduce=kron_reduce,
        phase_project=phase_project,
        ravens2math_extensions=ravens2math_extensions,
        ravens2math_passthrough=ravens2math_passthrough,
        global_keys=global_keys,
    )

    correct_network_data && correct_network_data!(data_math; make_pu=make_pu, make_pu_extensions=make_pu_extensions)

    return data_math

end


"base function for converting ravens model to mathematical model"
function _map_ravens2math(
    data_ravens::Dict{String,<:Any};
    multinetwork::Bool=false,
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
                    "vbases_default" => Dict{String,Real}(),
                    "vbases_network" => Dict{String,Real}(),
                    "vbases_buses" => Dict{String,Real}()

    )

    # Multinetwork
    if multinetwork
        data_math = Dict{String,Any}(
            "name" => get(_data_ravens, "name", ""),
            "data_model" => MATHEMATICAL,
            "multinetwork" => multinetwork,
            "nw" => Dict{String,Any}()
        )

        if haskey(data_ravens, "BasicIntervalSchedule")
            schdls = get(data_ravens, "BasicIntervalSchedule", Dict{String,Any}()) # Get schedules/timeseries

            # Check for shortest timeseries. Use that length for mn
            min_length = Inf
            for (name, ravens_obj) in schdls
                if haskey(ravens_obj, "EnergyConsumerSchedule.RegularTimePoints")
                    points = ravens_obj["EnergyConsumerSchedule.RegularTimePoints"]
                    length_points = length(points)

                    if length_points < min_length
                        min_length = length_points
                    end
                end
            end

            # Throw error if no schedule/timeseries is found.
            if min_length == Inf
                throw("Minimum length is equal to Inf. Cannot do multinetwork with Inf timepoints.")
            end

            # Vector to store nws dictionaries
            nws_vect = Vector{Dict{String,Any}}(undef, min_length)

            # Multithreaded loop to create each nw dictionary. Store them in vector (to allow multithreading)
            Threads.@threads for n=1:1:min_length
                nw_dict = Dict{String,Any}(
                    string(n) => Dict{String,Any}(
                        "per_unit" => get(_data_ravens, "per_unit", false),
                        "is_projected" => get(_data_ravens, "is_projected", false),
                        "is_kron_reduced" => get(_data_ravens, "is_kron_reduced", true), # TODO: Kron reduction?
                        "settings" => deepcopy(_settings),
                        "time_elapsed" => get(_data_ravens, "time_elapsed", 1.0),
                    )
                )

                # Store nw dict in vector
                nws_vect[n] = nw_dict

                # Perform conversion ravens2math
                apply_pmd!(_map_ravens2math_nw!, nws_vect[n], _data_ravens; ravens2math_passthrough=ravens2math_passthrough, ravens2math_extensions=ravens2math_extensions, nw=n)

            end

            # Merge dict in vector into data_math dictionary (other for loop to allow multithreading)
            for nw_dict in nws_vect
                merge!(data_math["nw"], nw_dict)
            end

        else
            @error("No timeseries and/or multinetwork information detected.")
        end

    else # No multinetwork
        data_math = Dict{String,Any}(
            "name" => get(_data_ravens, "name", ""),
            "per_unit" => get(_data_ravens, "per_unit", false),
            "data_model" => MATHEMATICAL,
            "is_projected" => get(_data_ravens, "is_projected", false),
            "is_kron_reduced" => get(_data_ravens, "is_kron_reduced", true), # TODO: Kron reduction?
            "settings" => deepcopy(_settings),
            "time_elapsed" => get(_data_ravens, "time_elapsed", 1.0),
        )
        apply_pmd!(_map_ravens2math_nw!, data_math, _data_ravens; ravens2math_passthrough=ravens2math_passthrough, ravens2math_extensions=ravens2math_extensions)
    end

    # multinetwork collections
    if multinetwork
        _collect_nw_maps!(data_math)
        _collect_nw_bus_lookups!(data_math)
    end

    return data_math
end


"""
"""
function _map_ravens2math_nw!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; ravens2math_passthrough::Dict{String,Vector{String}}=Dict{String,Vector{String}}(), ravens2math_extensions::Vector{<:Function}=Function[], nw::Int=nw_id_default)

    if nw==0

        data_math["map"] = Vector{Dict{String,Any}}([
            Dict{String,Any}("unmap_function" => "_map_math2ravens_root!")
        ])

        _init_base_components!(data_math)

        for property in get(ravens2math_passthrough, "root", String[])
            if haskey(data_ravens, property)
                data_math[property] = deepcopy(data_ravens[property])
            end
        end

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

    else

        data_math[string(nw)]["map"] = Vector{Dict{String,Any}}([
            Dict{String,Any}("unmap_function" => "_map_math2ravens_root!")
        ])

        _init_base_components!(data_math[string(nw)])

        for property in get(ravens2math_passthrough, "root", String[])
            if haskey(data_ravens, property)
                data_math[property] = deepcopy(data_ravens[property])
            end
        end

        for type in pmd_ravens_asset_types
            getfield(PowerModelsDistribution, Symbol("_map_ravens2math_$(type)!"))(data_math[string(nw)], data_ravens; pass_props=get(ravens2math_passthrough, type, String[]), nw=nw)
        end

        # Custom ravens2math transformation functions
        for ravens2math_func! in ravens2math_extensions
            ravens2math_func!(data_math[string(nw)], data_ravens)
        end

        find_conductor_ids!(data_math[string(nw)])
        _map_conductor_ids!(data_math[string(nw)])
        _map_settings_vbases_default!(data_math[string(nw)])

    end

end


"""
Converts ravens connectivity_node components into mathematical bus components.
"""
function _map_ravens2math_connectivity_node!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[], nw::Int=nw_id_default)
    voltage_scale_factor_sqrt3 = data_math["settings"]["voltage_scale_factor"] * sqrt(3)
    connectivity_nodes = get(data_ravens, "ConnectivityNode", Dict{String,Any}())

    for (name, ravens_obj) in connectivity_nodes
        index = length(data_math["bus"]) + 1
        math_obj = _init_math_obj_ravens("bus", name, ravens_obj, index; pass_props=pass_props)

        # Set basic bus properties
        math_obj["bus_i"] = index
        math_obj["source_id"] = "bus.$name"
        math_obj["bus_type"] = 1  # Default bus_type, will be modified as needed
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
function _map_ravens2math_conductor!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[], nw::Int=nw_id_default)

    if haskey(data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"], "Conductor")
        conductors = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]["Conductor"]

        for (name, ravens_obj) in get(conductors, "ACLineSegment", Dict{Any,Dict{String,Any}}())
            math_obj = _init_math_obj_ravens("conductor", name, ravens_obj, length(data_math["branch"]) + 1; pass_props=pass_props)
            nconds = length(ravens_obj["ACLineSegment.ACLineSegmentPhase"]) # number of conductors/wires
            nphases = 0 # init number of phases
            terminals = ravens_obj["ConductingEquipment.Terminals"]

            f_node = _extract_name(terminals[1]["Terminal.ConnectivityNode"])
            t_node = _extract_name(terminals[2]["Terminal.ConnectivityNode"])

            math_obj["f_bus"] = data_math["bus_lookup"][f_node]
            math_obj["t_bus"] = data_math["bus_lookup"][t_node]

            bus_terminals = nconds >= 3 ? collect(1:nconds) : [_phase_map[phase["ACLineSegmentPhase.phase"]] for phase in ravens_obj["ACLineSegment.ACLineSegmentPhase"]]

            # TODO: Kron reduce bus terminals by removing conn 4
            reduce = false # flag for Kron reduction
            if 4 in bus_terminals
                reduce = true
                bus_terminals = filter!(x -> x != 4, bus_terminals)
                nphases = nconds - 1
            else
                nphases = nconds
            end

            # Add vmin/vmax/terminals info to fbus and tbus if missing
            for bus in [math_obj["f_bus"], math_obj["t_bus"]]
                bus_data = data_math["bus"][string(bus)]
                if !(haskey(bus_data, "terminals")) || (length(bus_data["terminals"]) < length(bus_terminals))
                    bus_data["terminals"] = bus_terminals
                    bus_data["vmin"] = fill(0.0, nphases)
                    bus_data["vmax"] = fill(Inf, nphases)
                    bus_data["grounded"] = zeros(Bool, nphases)
                end
            end

            math_obj["f_connections"] = bus_terminals
            math_obj["t_connections"] = bus_terminals

            # System frequency
            base_freq = data_math["settings"]["base_frequency"]

            if (haskey(ravens_obj, "ACLineSegment.PerLengthImpedance"))

                impedance_name = _extract_name(ravens_obj["ACLineSegment.PerLengthImpedance"])
                impedance_data = data_ravens["PerLengthLineParameter"]["PerLengthImpedance"]["PerLengthPhaseImpedance"][impedance_name]

                math_obj["br_r"] = _impedance_conversion_ravens(impedance_data, ravens_obj, "PhaseImpedanceData.r")
                math_obj["br_x"] = _impedance_conversion_ravens(impedance_data, ravens_obj, "PhaseImpedanceData.x")

                for (key, param) in [("b_fr", "PhaseImpedanceData.b"), ("b_to", "PhaseImpedanceData.b"), ("g_fr", "PhaseImpedanceData.g"), ("g_to", "PhaseImpedanceData.g")]
                    math_obj[key] = _admittance_conversion_ravens(impedance_data, ravens_obj, param)
                end

            elseif (haskey(ravens_obj, "ACLineSegment.WireSpacingInfo"))

                # Get WireSpacingInfo
                spacinginfo_name = _extract_name(ravens_obj["ACLineSegment.WireSpacingInfo"])
                spacinginfo_data = data_ravens["AssetInfo"]["WireSpacingInfo"][spacinginfo_name]
                wire_positions = spacinginfo_data["WireSpacingInfo.WirePositions"]
                num_of_wires = length(wire_positions)

                # Coordinates
                x_coords = Vector{Float64}(undef, num_of_wires)
                y_coords = Vector{Float64}(undef, num_of_wires)

                for i in 1:1:num_of_wires
                    x_coords[i] = get(wire_positions[i], "WirePosition.xCoord", 0.0)
                    y_coords[i] = get(wire_positions[i], "WirePosition.yCoord", 0.0)
                end

                # angular frequency
                ω = 2π * base_freq
                ω₀ = 2π * base_freq

                # Get data for each specific ACLineSegmentPhase
                segmentphase_data = ravens_obj["ACLineSegment.ACLineSegmentPhase"]

                # Wire Info.
                gmr = Vector{Float64}(undef, nconds)   # gmr of Wire, default: radius of wire * 0.7788
                radius = Vector{Float64}(undef, nconds)    # radius of Wire
                rac = Vector{Float64}(undef, nconds)   # AC resistance
                rdc = Vector{Float64}(undef, nconds)   # DC resistance, default: AC resistance / 1.02
                dcable = Vector{Float64}(undef, nconds)   # diameter of Wire: radius of wire * 2

                # Concentric Neutrals Info.
                rstrand = Vector{Float64}(undef, nconds)  # resistance of CN cable
                nstrand = Vector{Float64}(undef, nconds)  # number of CN conductors
                dstrand = Vector{Float64}(undef, nconds)  # diameter of CN conductor
                gmrstrand = Vector{Float64}(undef, nconds)    # gmr of CN conductor, default: radius of CN * 0.7788

                # insulation info.
                dins = Vector{Float64}(undef, nconds) # diameter over insulation (over jacket)
                tins = Vector{Float64}(undef, nconds) # thickness of insulation

                # tape shield info.
                diashield = Vector{Float64}(undef, nconds)
                tapelayer = Vector{Float64}(undef, nconds)
                tapelap = Vector{Float64}(undef, nconds)

                for i in 1:1:nconds

                    wireinfo_name = _extract_name(segmentphase_data[i]["PowerSystemResource.AssetDatasheet"])
                    wireinfo_data = data_ravens["AssetInfo"]["WireInfo"][wireinfo_name]

                    radius[i] = get(wireinfo_data, "WireInfo.radius", NaN)
                    @assert  radius[i] != NaN "WireInfo radius not found! using NaN. Revise data."

                    # Note: gets rewritten as missing if not needed
                    dcable[i] = radius[i] * 2.0

                    gmr[i] = get(wireinfo_data, "WireInfo.gmr", radius[i] * 0.778)

                    if wireinfo_data["Ravens.cimObjectType"] == "OverheadWireInfo"
                        rac[i] = get(wireinfo_data, "WireInfo.rAC25", NaN)
                        @assert rac[i] != NaN "WireInfo AC25 resistance is not found! using NaN. Revise input data."
                        rdc[i] = rac[i] / 1.02
                    elseif wireinfo_data["Ravens.cimObjectType"] == "ConcentricNeutralCableInfo"
                        rdc[i] = get(wireinfo_data, "WireInfo.rDC20", NaN)
                        @assert rdc[i] != NaN "WireInfo rDC20 resistance is not found! using NaN. Revise input data."
                        rac[i] = rdc[i] * 1.02
                    else
                        @error("Cable type not supported. Resistances (AC or DC) not found!")
                    end

                    # Concentric Neutrals Information.
                    rstrand[i] = get(wireinfo_data, "ConcentricNeutralCableInfo.neutralStrandRDC20", NaN)
                    nstrand[i] = get(wireinfo_data, "ConcentricNeutralCableInfo.neutralStrandCount", NaN)
                    dstrand[i] = get(wireinfo_data, "ConcentricNeutralCableInfo.neutralStrandRadius", NaN) * 2.0
                    gmrstrand[i] = get(wireinfo_data, "ConcentricNeutralCableInfo.neutralStrandGmr", (dstrand[i]/2.0) * 0.778)

                    # insulation information
                    dins[i] = get(wireinfo_data, "CableInfo.diameterOverJacket", NaN)
                    tins[i] = get(wireinfo_data, "WireInfo.insulationThickness", NaN)

                    # TODO: tape shielded cables information
                    diashield[i] = NaN
                    tapelayer[i] = NaN
                    tapelap[i] = NaN

                end

                # Check for NaNs and replace with missing.
                rstrand = findfirst(isnan, rstrand) !== nothing ? missing : rstrand
                nstrand = findfirst(isnan, nstrand) !== nothing ? missing : nstrand
                dcable = findfirst(isnan, dstrand) !== nothing ? missing : dcable   # use dstrand as signal for dcable to be missing
                dstrand = findfirst(isnan, dstrand) !== nothing ? missing : dstrand
                gmrstrand = findfirst(isnan, gmrstrand) !== nothing ? missing : gmrstrand
                epsr = findfirst(isnan, dins) !== nothing ? missing : ones(nconds).*2.3 # use dins as signal for epsr to be missing
                dins = findfirst(isnan, dins) !== nothing ? missing : dins
                tins = findfirst(isnan, tins) !== nothing ? missing : tins
                diashield = findfirst(isnan, diashield) !== nothing ? missing : diashield
                tapelayer = findfirst(isnan, tapelayer) !== nothing ? missing : tapelayer
                tapelap = findfirst(isnan, tapelap) !== nothing ? missing : tapelap

                # TODO: earth model (using default)
                earth_model = "deri"

                # rho (default) - ρ = earth resistivity = 100 Ω-m
                rho = 100

                # Calculate line constants
                z, y =  calculate_line_constants(
                    x_coords,
                    y_coords,
                    ω,
                    gmr,
                    radius,
                    nconds,
                    earth_model,
                    rac,
                    ω₀,
                    rdc,
                    rho,
                    nphases,
                    rstrand,
                    nstrand,
                    dcable,
                    dstrand,
                    gmrstrand,
                    epsr,
                    dins,
                    tins,
                    diashield,
                    tapelayer,
                    tapelap
                )

                # Kron reduction
                if reduce
                    z, y = _kron(z, y, nphases)
                end

                rs, xs = real(z), imag(z)
                g, b = real(y), imag(y)

                b_fr = (b ./ 2.0) .* base_freq
                b_to = (b ./ 2.0) .* base_freq
                g_fr = (g ./ 2.0) .* base_freq
                g_to = (g ./ 2.0) .* base_freq

                math_obj["br_r"] = _impedance_conversion_ravens(ravens_obj, rs)
                math_obj["br_x"] = _impedance_conversion_ravens(ravens_obj, xs)

                math_obj["b_fr"] = _admittance_conversion_ravens(ravens_obj, b_fr)
                math_obj["b_to"] = _admittance_conversion_ravens(ravens_obj, b_to)

                math_obj["g_fr"] = _admittance_conversion_ravens(ravens_obj, g_fr)
                math_obj["g_to"] = _admittance_conversion_ravens(ravens_obj, g_to)

            end

            math_obj["angmin"] = get(ravens_obj, "vad_lb", fill(-60.0, nphases))
            math_obj["angmax"] = get(ravens_obj, "vad_ub", fill(60.0, nphases))

            if (haskey(terminals[1], "ACDCTerminal.OperationalLimitSet"))
                oplimitset_id = _extract_name(terminals[1]["ACDCTerminal.OperationalLimitSet"])
                oplimitset = data_ravens["OperationalLimitSet"][oplimitset_id]["OperationalLimitSet.OperationalLimitValue"][1]  # [1] emerg, [2] normal
            else
                oplimitset = Dict()
            end

            limit_keys = [("CurrentLimit.value", "c_rating_a"), ("CurrentLimit.value", "c_rating_b"), ("CurrentLimit.value", "c_rating_c"),
                        ("ApparentPowerLimit.value", "rate_a"), ("ApparentPowerLimit.value", "rate_b"), ("ApparentPowerLimit.value", "rate_c")]

            for (f_key, t_key) in limit_keys
                math_obj[t_key] = haskey(oplimitset, f_key) ? fill(oplimitset[f_key], nphases) : fill(Inf, nphases)
            end

            math_obj["br_status"] = get(ravens_obj, "Equipment.inService", true) == true ? 1 : 0
            data_math["branch"]["$(math_obj["index"])"] = math_obj

            push!(data_math["map"], Dict{String,Any}(
                "from" => name,
                "to" => "branch.$(math_obj["index"])",
                "unmap_function" => "_map_math2ravens_line!",
            ))

        end

    end
end


# TODO: Transformers need a lot of changes/refactors!!!
"converts ravens n-winding transformers into mathematical ideal 2-winding lossless transformer branches and impedance branches to represent the loss model"
function _map_ravens2math_power_transformer!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[], nw::Int=nw_id_default)

    conducting_equipment = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]
    power_scale_factor = data_math["settings"]["power_scale_factor"]
    voltage_scale_factor = data_math["settings"]["voltage_scale_factor"]
    voltage_scale_factor_sqrt3 = voltage_scale_factor * sqrt(3)

    for (name, ravens_obj) in get(conducting_equipment, "PowerTransformer", Dict{Any,Dict{String,Any}}())

        # Build map first, so we can update it as we decompose the transformer
        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => String[],
            "unmap_function" => "_map_math2ravens_transformer!",
        ))

        to_map = data_math["map"][end]["to"]

        if haskey(ravens_obj, "PowerTransformer.PowerTransformerEnd")

            # Get nrw: number of windings
            wdgs = ravens_obj["PowerTransformer.PowerTransformerEnd"]
            nrw = length(wdgs)

            # connections
            connections = Vector{Vector{Int64}}(undef, nrw)

            # configurations
            wdgs_confs = Vector{ConnConfig}(undef, nrw)

            # RegulatorControls flag
            reg_controls = [false for _ in 1:nrw]
            reg_obj = [Dict() for _ in 1:nrw]

            # Transformer data for each winding
            vnom = Vector{Float64}(undef, nrw)
            snom = Vector{Float64}(undef, nrw)
            zbase = Vector{Float64}(undef, nrw)
            x_sc = Vector{Float64}(undef, nrw)
            r_s = Vector{Float64}(undef, nrw)
            g_sh = zeros(Float64, nrw)
            b_sh = zeros(Float64, nrw)

            # Init RatioTapChanger data (default)
            tm_set = Vector{Vector{Float64}}(undef, nrw)
            tm_lb = Vector{Vector{Float64}}(undef, nrw)
            tm_ub = Vector{Vector{Float64}}(undef, nrw)
            tm_fix = Vector{Vector{Bool}}(undef, nrw)
            tm_step = Vector{Vector{Float64}}(undef, nrw)

            for wdg_id in 1:nrw

                # wdg phasecode & terminals
                wdg_terminals = ravens_obj["ConductingEquipment.Terminals"][wdg_id]
                wdg_phasecode = wdg_terminals["Terminal.phases"]

                # wdg endNumber
                wdg_endNumber = wdgs[wdg_id]["TransformerEnd.endNumber"]

                # Connections (based on _phasecode_map)
                if haskey(_phasecode_map, wdg_phasecode)
                    connections[wdg_endNumber] = _phasecode_map[wdg_phasecode]
                else
                    @error("PhaseCode not supported yet!")
                end

                # nphases
                nphases = length(connections[wdg_endNumber])

                # Add terminals and voltage limits info. if missing
                node = _extract_name(wdg_terminals["Terminal.ConnectivityNode"])
                bus = data_math["bus_lookup"][node]
                bus_data = data_math["bus"][string(bus)]
                if !(haskey(bus_data, "terminals")) || (length(bus_data["terminals"]) < length(connections[wdg_endNumber]))
                    bus_data["terminals"] = connections[wdg_endNumber]
                    bus_data["vmin"] = fill(0.0, nphases)
                    bus_data["vmax"] = fill(Inf, nphases)
                    bus_data["grounded"] = zeros(Bool, nphases)
                end

                # wdgs configurations
                if wdgs[wdg_endNumber]["PowerTransformerEnd.connectionKind"] == "WindingConnection.Y" || wdgs[wdg_endNumber]["PowerTransformerEnd.connectionKind"] == "WindingConnection.Yn"
                    wdgs_confs[wdg_endNumber] = WYE
                elseif wdgs[wdg_endNumber]["PowerTransformerEnd.connectionKind"] == "WindingConnection.D"
                    wdgs_confs[wdg_endNumber] = DELTA
                else
                    @error("PowerTransformer ConnectionKind not supported yet!")
                end

                # Transformer data for each winding
                vnom[wdg_endNumber] = wdgs[wdg_endNumber]["PowerTransformerEnd.ratedU"]
                snom[wdg_endNumber] = wdgs[wdg_endNumber]["PowerTransformerEnd.ratedS"]

                # Add vnom info to vbases
                data_math["settings"]["vbases_network"][string(bus)] = deepcopy(vnom[wdg_endNumber]/voltage_scale_factor)

                # resistance
                transf_star_impedance = get(wdgs[wdg_endNumber], "TransformerEnd.StarImpedance", Dict())
                r_s[wdg_endNumber] = get(wdgs[wdg_endNumber], "PowerTransformerEnd.r",
                                        get(transf_star_impedance, "TransformerStarImpedance.r", 0.0))

                # reactance
                transf_mesh_impedance = get(wdgs[wdg_endNumber], "TransformerEnd.MeshImpedance", Dict())
                x_sc[wdg_endNumber] = get(transf_mesh_impedance, "TransformerMeshImpedance.x",
                                        get(transf_star_impedance, "TransformerStarImpedance.x", 0.0))

                # admittance
                transf_core_impedance = get(wdgs[wdg_endNumber], "TransformerEnd.CoreAdmittance", Dict())
                g_sh[wdg_id] =  get(transf_core_impedance, "TransformerCoreAdmittance.g", 0.0)
                b_sh[wdg_id] = - get(transf_core_impedance, "TransformerCoreAdmittance.b", 0.0)


                # Set RatioTapChanger in specific wdg
                if haskey(wdgs[wdg_endNumber], "TransformerEnd.RatioTapChanger")

                    rtc_name = _extract_name(wdgs[wdg_endNumber]["TransformerEnd.RatioTapChanger"])
                    rtc_data = data_ravens["PowerSystemResource"]["TapChanger"]["RatioTapChanger"][rtc_name]

                    # tm_step
                    hstep = get(rtc_data, "TapChanger.highStep", 16)
                    lstep = get(rtc_data, "TapChanger.lowStep", -16)
                    step_dist = abs(hstep) + abs(lstep)
                    step_tap = 1/step_dist
                    tm_step[wdg_endNumber] = fill(step_tap, nphases)

                    # tm_set
                    step = get(rtc_data, "TapChanger.step", 1.0)    # Starting Tap changer position/step
                    tm_set[wdg_endNumber] = fill(step, nphases)

                    # tm_fix
                    ltcFlag = get(rtc_data, "TapChanger.ltcFlag", false)
                    if (ltcFlag == true)
                        tm_fix[wdg_endNumber] = zeros(Bool, nphases)
                    else
                        tm_fix[wdg_endNumber] = ones(Bool, nphases)
                    end

                    # tm_ub/tm_lb
                    neutralVoltPu = get(rtc_data, "TapChanger.neutralU", vnom[wdg_endNumber])/vnom[wdg_endNumber]
                    step_volt_increment = get(rtc_data, "RatioTapChanger.stepVoltageIncrement", 100.0)
                    volt_lb = neutralVoltPu + step_tap * (step_volt_increment/100.0) * lstep
                    volt_ub = neutralVoltPu + step_tap * (step_volt_increment/100.0) * hstep
                    tm_lb[wdg_endNumber] = fill(volt_lb, nphases)
                    tm_ub[wdg_endNumber] = fill(volt_ub, nphases)

                    # Regulator Control
                    if haskey(rtc_data, "TapChanger.TapChangerControl") && !all(tm_fix[wdg_endNumber])
                        reg_controls[wdg_endNumber] = true
                        reg_obj[wdg_endNumber] = Dict{String,Any}(
                                "vreg" => fill(rtc_data["TapChanger.TapChangerControl"]["RegulatingControl.targetValue"], nphases),
                                "band" =>  fill(rtc_data["TapChanger.TapChangerControl"]["RegulatingControl.targetDeadband"], nphases),
                                "ptratio" => fill(rtc_data["TapChanger.ptRatio"], nphases),
                                "ctprim" => fill(rtc_data["TapChanger.ctRating"], nphases),
                                "r" => fill(rtc_data["TapChanger.TapChangerControl"]["TapChangerControl.lineDropR"], nphases),
                                "x" => fill(rtc_data["TapChanger.TapChangerControl"]["TapChangerControl.lineDropX"], nphases)
                            )
                    end

                else # default
                    tm_set[wdg_id] = fill(1.0, nphases)
                    tm_lb[wdg_id] = fill(0.9, nphases)
                    tm_ub[wdg_id] = fill(1.1, nphases)
                    tm_fix[wdg_id] = ones(Bool, nphases)
                    tm_step[wdg_id] = fill(1/32, nphases)
                end

            end

            # data is measured externally, but we now refer it to the internal side - some values are referred to wdg 1
            ratios = vnom/voltage_scale_factor
            x_sc = (x_sc[1]./ratios[1]^2)./100.0
            r_s = r_s./ratios.^2
            g_sh = g_sh[1]*ratios[1]^2
            b_sh = b_sh[1]*ratios[1]^2

            # convert x_sc from list of upper triangle elements to an explicit dict
            y_sh = g_sh + im*b_sh
            z_sc = Dict([(key, im*x_sc[i]) for (i,key) in enumerate([(i,j) for i in 1:nrw for j in i+1:nrw])])

            # dimesions
            dims = length(tm_set[1])

            # init polarity
            polarity = fill(1, nrw)

            # Status
            status = haskey(ravens_obj, "Equipment.inService") ? ravens_obj["Equipment.inService"] : true
            status = status == true ? 1 : 0

            # Build loss model
            transformer_t_bus_w = _build_loss_model!(data_math, name, to_map, r_s, z_sc, y_sh, connections[1]; nphases=dims, status=status)

            # Mathematical model for transformer
            for wdg_id in 1:nrw

                # 2-WINDING TRANSFORMER
                # correct polarity and connections
                if wdg_id > 1
                    if wdgs_confs[1] == DELTA && wdgs_confs[wdg_id] == WYE
                        polarity[wdg_id] = -1
                        connections[wdg_id] = _barrel_roll(connections[wdg_id][1:end], 1)
                    end
                    if wdgs_confs[1] == WYE && wdgs_confs[wdg_id] == DELTA
                        polarity[wdg_id] = -1
                        connections[wdg_id] = _barrel_roll(connections[wdg_id], -1)
                    end
                end

                # make virtual bus and mark it for reduction
                tm_nom = wdgs_confs[wdg_id]==DELTA ? wdgs[wdg_id]["PowerTransformerEnd.ratedU"]*sqrt(3)/voltage_scale_factor : wdgs[wdg_id]["PowerTransformerEnd.ratedU"]/voltage_scale_factor

                # Get correct f_node for winding
                wdg_term = ravens_obj["ConductingEquipment.Terminals"][wdg_id]
                f_node_wdgterm = _extract_name(wdg_term["Terminal.ConnectivityNode"])

                # Transformer Object
                transformer_2wa_obj = Dict{String,Any}(
                    "name"          => "_virtual_transformer.$name.$wdg_id",
                    "source_id"     => "_virtual_transformer.transformer.$name.$wdg_id",
                    "f_bus"         => data_math["bus_lookup"][f_node_wdgterm],
                    "t_bus"         => transformer_t_bus_w[wdg_id],
                    "tm_nom"        => tm_nom,
                    "f_connections" => connections[wdg_id],
                    "t_connections" => connections[1],
                    "configuration" => wdgs_confs[wdg_id],
                    "polarity"      => polarity[wdg_id],
                    "tm_set"        => tm_set[wdg_id],
                    "tm_fix"        => tm_fix[wdg_id],
                    "sm_ub"         => get(wdgs[wdg_id], "PowerTransformerEnd.ratedS", Inf)/power_scale_factor,
                    "cm_ub"         => get(wdgs[wdg_id], "PowerTransformerEnd.ratedI", Inf),
                    "status"        => status,
                    "index"         => length(data_math["transformer"])+1
                )


                transformer_2wa_obj["tm_lb"] = tm_lb[wdg_id]
                transformer_2wa_obj["tm_ub"] = tm_ub[wdg_id]
                transformer_2wa_obj["tm_step"] = tm_step[wdg_id]

                data_math["transformer"]["$(transformer_2wa_obj["index"])"] = transformer_2wa_obj

                # Add Regulator Controls (only if flag is true)
                if (reg_controls[wdg_id] == true)
                    data_math["transformer"]["$(transformer_2wa_obj["index"])"]["controls"] = reg_obj[wdg_id]
                end

                # TODO: Center-Tapped Transformers (3 Windings)
                # if w==3 && eng_obj["polarity"][w]==-1 # identify center-tapped transformer and mark all secondary-side nodes as triplex by adding va_start
                # end

                push!(to_map, "transformer.$(transformer_2wa_obj["index"])")

            end

        elseif haskey(ravens_obj, "PowerTransformer.TransformerTank")

            # Get tanks data
            tanks = ravens_obj["PowerTransformer.TransformerTank"]

            # TODO: flag for debugging TransformerTanks models
            combine_tanks = true

            # Combine tanks into a single transformer model (default)
            if (combine_tanks == true)

                # number of tanks
                ntanks = length(tanks)

                # TODO: IMPORTANT ASSUMPTIONS
                # 1) assume there is at least 1 tank and that all tanks have the same number of windings (i.e., TransformerTankEnds)
                # 2) assume the number of phases is equal to the number of tanks
                nphases = length(tanks) # assume nphases == ntanks
                nrw = length(tanks[1]["TransformerTank.TransformerTankEnd"])

                # init connections vector for combined transformer windings
                connections = [zeros(Int64, nphases) for _ in 1:nrw]

                # init nodes vector for combined transformer windings
                nodes = ["" for _ in 1:nrw]

                # init rs, x_sc, g_sh, and b_sh data per wdg/tank(phase)
                r_s = [zeros(Float64, nphases) for _ in 1:nrw]
                x_sc = [zeros(Float64, nphases) for _ in 1:nrw]
                g_sh = zeros(Float64, nphases)
                b_sh = zeros(Float64, nphases)

                # init sm_ub and cm_ub
                sm_ub = zeros(Float64, nrw)
                cm_ub = zeros(Float64, nrw)

                # init configuration - default WYE-WYE
                configuration = [WYE for _ in 1:nrw]

                # RegulatorControls flag
                reg_controls = [false for _ in 1:nrw]
                reg_obj = [Dict() for _ in 1:nrw]

                # init vnom for all windings
                vnom = zeros(Float64, nrw)

                # temp store previous for checking
                nodes_prev = []
                configuration_prev = []
                vnom_prev = []

                # Init RatioTapChanger data (default)
                tm_set = Vector{Vector{Float64}}(fill(fill(1.0, nphases), nrw))
                tm_lb = Vector{Vector{Float64}}(fill(fill(0.9, nphases), nrw))
                tm_ub = Vector{Vector{Float64}}(fill(fill(1.1, nphases), nrw))
                tm_fix = Vector{Vector{Bool}}(fill(ones(Bool, nphases), nrw))
                tm_step = Vector{Vector{Float64}}(fill(fill(1/32, nphases), nrw))

                for tank_id in 1:ntanks

                    # Get wdg data
                    wdgs = tanks[tank_id]["TransformerTank.TransformerTankEnd"]

                    # Tank Asset
                    tank_asset_name = _extract_name(tanks[tank_id]["PowerSystemResource.AssetDatasheet"])
                    tank_asset_data = data_ravens["AssetInfo"]["PowerTransformerInfo"][tank_asset_name]

                    for wdg_id in 1:nrw

                        # wdg terminals & phasecode
                        wdg_terminals = ravens_obj["ConductingEquipment.Terminals"][wdg_id]
                        wdg_phasecode = wdg_terminals["Terminal.phases"]

                        # wdg endNumber
                        wdg_endNumber = wdgs[wdg_id]["TransformerEnd.endNumber"]

                        # from-and-to-nodes for wdgs
                        nodes[wdg_endNumber] = _extract_name(wdg_terminals["Terminal.ConnectivityNode"])

                        # Connections (based on _phasecode_map)
                        if haskey(_phasecode_map, wdg_phasecode)
                            phasecode_conns = _phasecode_map[wdg_phasecode]
                            if !(length(phasecode_conns)>1)
                                connections[wdg_endNumber][tank_id] = phasecode_conns[1]
                            else
                                connections[wdg_endNumber] = phasecode_conns
                            end
                        else
                            @error("PhaseCode not supported yet!")
                        end

                        # transformer tank end info.
                        transf_end_info = tank_asset_data["PowerTransformerInfo.TransformerTankInfos"][tank_asset_name]["TransformerTankInfo.TransformerEndInfos"]
                        vnom_wdg = transf_end_info[wdg_endNumber]["TransformerEndInfo.ratedU"]
                        snom_wdg = transf_end_info[wdg_endNumber]["TransformerEndInfo.ratedS"]
                        zbase = (vnom_wdg^2) / snom_wdg
                        ratios = vnom_wdg/voltage_scale_factor

                        # assign vnom_wdg to vnom for transformer
                        vnom[wdg_endNumber] = vnom_wdg

                        # resistance computation
                        transf_star_impedance = get(transf_end_info[wdg_endNumber], "TransformerEndInfo.TransformerStarImpedance", Dict())
                        r_s[wdg_endNumber][tank_id] = get(transf_end_info[wdg_endNumber], "TransformerEndInfo.r",
                                            get(transf_star_impedance, "TransformerStarImpedance.r", 0.0))
                        r_s[wdg_endNumber][tank_id] = r_s[wdg_endNumber][tank_id]/ratios^2

                        # reactance computation
                        x_sc[wdg_endNumber][tank_id] = get(transf_end_info[wdg_endNumber], "TransformerEndInfo.x",
                                            get(transf_star_impedance, "TransformerStarImpedance.x", 0.0))
                        # -- alternative computation of xsc using sc tests
                        if haskey(transf_end_info[wdg_endNumber], "TransformerEndInfo.EnergisedEndShortCircuitTests")
                            leak_impedance_wdg = transf_end_info[wdg_endNumber]["TransformerEndInfo.EnergisedEndShortCircuitTests"][1]["ShortCircuitTest.leakageImpedance"]
                            x_sc[wdg_endNumber][tank_id] = (sqrt((leak_impedance_wdg / zbase)^2 - ((r_s[wdg_endNumber][tank_id]*100.0 / zbase)*2)^2)/100.0)*zbase
                        end
                        x_sc[wdg_endNumber][tank_id] = (x_sc[wdg_endNumber][tank_id]/ratios^2)

                        # g_sh always with respect to wdg #1 always
                        if wdg_endNumber == 1
                            transf_end_noloadtest = get(transf_end_info[wdg_endNumber], "TransformerEndInfo.EnergisedEndNoLoadTests", [Dict()])
                            loss = get(transf_end_noloadtest[1], "NoLoadTest.loss", 0.0)
                            g_sh_tank =  (loss*snom_wdg)/zbase
                            exct_current = get(transf_end_noloadtest[1], "NoLoadTest.excitingCurrent", 0.0)
                            b_sh_tank = -((sqrt(abs((exct_current)^2 - (loss/(0.01*snom_wdg))^2)))/(100.0*zbase))
                            # data is measured externally, but we now refer it to the internal side
                            g_sh[tank_id] = g_sh_tank*ratios^2
                            b_sh[tank_id] = b_sh_tank*ratios^2
                        end

                        # configuration
                        conf = transf_end_info[wdg_endNumber]["TransformerEndInfo.connectionKind"]
                        if conf == "WindingConnection.Y" || conf == "WindingConnection.I" ||  conf == "WindingConnection.Yn"
                            configuration[wdg_endNumber] =  WYE
                        elseif conf == "WindingConnection.D"
                            configuration[wdg_endNumber] = DELTA
                        else
                            @error("TransformerTank ConnectionKind not supported yet!")
                        end

                        # add sm_ub if greater than existing (assumes the greatest value as the ratings for all phases in wdg)
                        semerg_wdg = get(transf_end_info[wdg_endNumber], "TransformerEndInfo.emergencyS", Inf)
                        if semerg_wdg > sm_ub[wdg_endNumber]
                            sm_ub[wdg_endNumber] = semerg_wdg
                        end

                        # add cm_ub if greater than existing for winding (assumes the greatest value as the ratings for all phases in wdg)
                        cm_wdg = get(transf_end_info[wdg_endNumber], "TransformerEndInfo.ratedI", Inf)
                        if cm_wdg > cm_ub[wdg_endNumber]
                            cm_ub[wdg_endNumber] = cm_wdg
                        end

                        # Set RatioTapChanger in specific wdg
                        if haskey(wdgs[wdg_endNumber], "TransformerEnd.RatioTapChanger")

                            rtc_name = _extract_name(wdgs[wdg_endNumber]["TransformerEnd.RatioTapChanger"])
                            rtc_data = data_ravens["PowerSystemResource"]["TapChanger"]["RatioTapChanger"][rtc_name]

                            # tm_step
                            hstep = get(rtc_data, "TapChanger.highStep", 16)
                            lstep = get(rtc_data, "TapChanger.lowStep", -16)
                            step_dist = abs(hstep) + abs(lstep)
                            step_tap = 1/step_dist
                            tm_step[wdg_endNumber] = fill(step_tap, nphases)

                            # tm_set
                            step = get(rtc_data, "TapChanger.step", 1.0)    # Starting Tap changer position/step
                            tm_set[wdg_endNumber] = fill(step, nphases)

                            # tm_fix
                            ltcFlag = get(rtc_data, "TapChanger.ltcFlag", false)
                            if (ltcFlag == true)
                                tm_fix[wdg_endNumber] = zeros(Bool, nphases)
                            else
                                tm_fix[wdg_endNumber] = ones(Bool, nphases)
                            end

                            # tm_ub/tm_lb
                            neutralVoltPu = get(rtc_data, "TapChanger.neutralU", vnom[wdg_endNumber])/vnom[wdg_endNumber]
                            step_volt_increment = get(rtc_data, "RatioTapChanger.stepVoltageIncrement", 100.0)
                            volt_lb = neutralVoltPu + step_tap * (step_volt_increment/100.0) * lstep
                            volt_ub = neutralVoltPu + step_tap * (step_volt_increment/100.0) * hstep
                            tm_lb[wdg_endNumber] = fill(volt_lb, nphases)
                            tm_ub[wdg_endNumber] = fill(volt_ub, nphases)

                            # Regulator Control
                            if haskey(rtc_data, "TapChanger.TapChangerControl") && !all(tm_fix[wdg_endNumber])
                                reg_controls[wdg_endNumber] = true
                                reg_obj[wdg_endNumber] = Dict{String,Any}(
                                        "vreg" => fill(rtc_data["TapChanger.TapChangerControl"]["RegulatingControl.targetValue"], nphases),
                                        "band" =>  fill(rtc_data["TapChanger.TapChangerControl"]["RegulatingControl.targetDeadband"], nphases),
                                        "ptratio" => fill(rtc_data["TapChanger.ptRatio"], nphases),
                                        "ctprim" => fill(rtc_data["TapChanger.ctRating"], nphases),
                                        "r" => fill(rtc_data["TapChanger.TapChangerControl"]["TapChangerControl.lineDropR"], nphases),
                                        "x" => fill(rtc_data["TapChanger.TapChangerControl"]["TapChangerControl.lineDropX"], nphases)
                                    )
                            end

                        end

                    end

                    ### --- Consistency checks across tanks ---
                    # check that nodes are the same after first tank iter
                    if tank_id != 1
                        @assert nodes == nodes_prev "nodes are not the same for all tanks! check ConnectivityNodes."  # check if node names are the same as expected
                    else
                        nodes_prev = deepcopy(nodes)
                    end

                    # check that configurations across tanks are consistent
                    if tank_id != 1
                        @assert configuration == configuration_prev "Configurations (e.g., WYE, DELTA) are not the same for all tanks and windings! check Configurations."  # check if node names are the same as expected
                    else
                        configuration_prev = deepcopy(configuration)
                    end

                    # check that vnoms across tanks are consistent for wdgs
                    if tank_id != 1
                        @assert vnom == vnom_prev "rated Voltages are not consistent for all tanks and windings! check TransformerEndInfo.ratedU values."  # check if node names are the same as expected
                    else
                        vnom_prev = deepcopy(vnom)
                    end

                end

                 # Add information about bus/node if missing
                for i in 1:length(nodes)
                    node = nodes[i]
                    bus = data_math["bus_lookup"][node]
                    bus_data = data_math["bus"][string(bus)]
                    if !(haskey(bus_data, "terminals")) || (length(bus_data["terminals"]) < length(connections[i]))
                        bus_data["terminals"] = connections[i]
                        bus_data["vmin"] = fill(0.0, nphases)
                        bus_data["vmax"] = fill(Inf, nphases)
                        bus_data["grounded"] = zeros(Bool, nphases)
                    end
                    # Add vnom info to bus
                    data_math["settings"]["vbases_network"][string(bus)] = deepcopy(vnom[i]/voltage_scale_factor)
                end

                # wdg i, tank 1  - assumes tank 1 always exists
                r_s = [r_s[i][1] for i in 1:nrw]
                x_sc = [x_sc[1][1]] # wrt to wdg 1
                g_sh = g_sh[1]      # wrt to wdg 1
                b_sh = b_sh[1]      # wrt to wdg 1

                # convert x_sc from list of upper triangle elements to an explicit dict
                y_sh = g_sh + im*b_sh
                z_sc = Dict([(key, im*x_sc[i]) for (i,key) in enumerate([(i,j) for i in 1:nrw for j in i+1:nrw])])

                # init Polarity
                polarity = fill(1, nrw)

                # Status
                status = haskey(ravens_obj, "Equipment.inService") ? ravens_obj["Equipment.inService"] : true
                status = status == true ? 1 : 0

                # Build loss model
                transformer_t_bus_w = _build_loss_model!(data_math, name, to_map, r_s, z_sc, y_sh, connections[1]; nphases=nphases, status=status)

                # Mathematical model for transformer
                for wdg_id in 1:nrw
                    # 2-WINDING TRANSFORMER

                    # correct polarity and connections
                    if wdg_id>1
                        if configuration[1] == DELTA && configuration[wdg_id] == WYE
                            polarity[wdg_id] = -1
                            connections[wdg_id] = _barrel_roll(connections[wdg_id][1:end], 1)
                        end
                        if configuration[1] == WYE && configuration[wdg_id] == DELTA
                            polarity[wdg_id] = -1
                            connections[wdg_id] = _barrel_roll(connections[wdg_id], -1)
                        end
                    end

                    # tm_nom depending on wdg configuration
                    tm_nom = configuration[wdg_id]==DELTA ? vnom[wdg_id]*sqrt(3)/voltage_scale_factor : vnom[wdg_id]/voltage_scale_factor

                    # Transformer Object
                    transformer_2wa_obj = Dict{String,Any}(
                        "name"          => "_virtual_transformer.$name.$wdg_id",
                        "source_id"     => "_virtual_transformer.transformer.$name.$wdg_id",
                        "f_bus"         => data_math["bus_lookup"][nodes[wdg_id]],
                        "t_bus"         => transformer_t_bus_w[wdg_id],
                        "tm_nom"        => tm_nom,
                        "f_connections" => connections[wdg_id],
                        "t_connections" => connections[1],
                        "configuration" => configuration[wdg_id],
                        "polarity"      => polarity[wdg_id],
                        "tm_set"        => tm_set[wdg_id],
                        "tm_fix"        => tm_fix[wdg_id],
                        "sm_ub"         => sm_ub[wdg_id]/power_scale_factor,
                        "cm_ub"         => cm_ub[wdg_id], # TODO: this may need scaling
                        "status"        => status,
                        "index"         => length(data_math["transformer"])+1
                    )

                    # RatioTapChanger
                    transformer_2wa_obj["tm_lb"] = tm_lb[wdg_id]
                    transformer_2wa_obj["tm_ub"] = tm_ub[wdg_id]
                    transformer_2wa_obj["tm_step"] = tm_step[wdg_id]

                    data_math["transformer"]["$(transformer_2wa_obj["index"])"] = transformer_2wa_obj

                    # Add Regulator Controls (only if flag is true)
                    if (reg_controls[wdg_id] == true)
                        data_math["transformer"]["$(transformer_2wa_obj["index"])"]["controls"] = reg_obj[wdg_id]
                    end

                    # TODO: Center-Tapped Transformers (3 Windings)
                    # if w==3 && eng_obj["polarity"][w]==-1 # identify center-tapped transformer and mark all secondary-side nodes as triplex by adding va_start
                    # end

                    push!(to_map, "transformer.$(transformer_2wa_obj["index"])")

                end

            else # Create a transformer for each tank
                for tank_id in 1:length(tanks)

                    tank_data = tanks[tank_id]   # tank data
                    wdgs_data = tank_data["TransformerTank.TransformerTankEnd"]     # wdgs data
                    tank_asset_name = _extract_name(tanks[tank_id]["PowerSystemResource.AssetDatasheet"])   # tank asset name
                    tank_asset_data = data_ravens["AssetInfo"]["PowerTransformerInfo"][tank_asset_name]     # tank asset data
                    nrw = length(tank_data["TransformerTank.TransformerTankEnd"])   # number of windings
                    nphases = 0 # init nphases var

                    # per tank windings connections
                    connections = Vector{Vector{Int64}}(undef, nrw)

                    # wdgs data vectors
                    vnom = zeros(Float64, nrw)
                    # init rs, x_sc, g_sh, and b_sh data per wdg
                    r_s = zeros(Float64, nrw)
                    x_sc = zeros(Float64, nrw)
                    g_sh = 0.0
                    b_sh = 0.0

                    # configurations
                    wdgs_confs = Vector{ConnConfig}(undef, nrw)

                    # RegulatorControls flag
                    reg_controls = [false for _ in 1:nrw]
                    reg_obj = [Dict() for _ in 1:nrw] # init

                    # Init RatioTapChanger data (default)
                    tm_set = Vector{Vector{Float64}}(undef, nrw)
                    tm_lb = Vector{Vector{Float64}}(undef, nrw)
                    tm_ub = Vector{Vector{Float64}}(undef, nrw)
                    tm_fix = Vector{Vector{Bool}}(undef, nrw)
                    tm_step = Vector{Vector{Float64}}(undef, nrw)

                    for wdg_id in 1:nrw

                        wdg_terminals = ravens_obj["ConductingEquipment.Terminals"][wdg_id]
                        wdg_phasecode = wdg_terminals["Terminal.phases"]
                        wdg_endNumber = wdgs_data[wdg_id]["TransformerEnd.endNumber"]

                        # connections (based on _phasecode_map)
                        if haskey(_phasecode_map, wdg_phasecode)
                            wdg_connections = _phasecode_map[wdg_phasecode]
                            connections[wdg_endNumber] = wdg_connections
                        else
                            @error("PhaseCode not supported yet!")
                        end

                        # from-and-to-nodes for wdg
                        nphases = length(wdg_connections)
                        node = _extract_name(wdg_terminals["Terminal.ConnectivityNode"])
                        bus = data_math["bus_lookup"][node]
                        bus_data = data_math["bus"][string(bus)]

                        # Add vmin/vmax/terminals info. if missing
                        if !(haskey(bus_data, "terminals")) || (length(bus_data["terminals"]) < length(wdg_connections))
                            bus_data["terminals"] = wdg_connections
                            bus_data["vmin"] = fill(0.0, nphases)
                            bus_data["vmax"] = fill(Inf, nphases)
                            bus_data["grounded"] = zeros(Bool, nphases)
                        end

                        # transformer tank end info.
                        transf_end_info = tank_asset_data["PowerTransformerInfo.TransformerTankInfos"][tank_asset_name]["TransformerTankInfo.TransformerEndInfos"]
                        vnom_wdg = transf_end_info[wdg_endNumber]["TransformerEndInfo.ratedU"]
                        snom_wdg = transf_end_info[wdg_endNumber]["TransformerEndInfo.ratedS"]

                        zbase = (vnom_wdg^2) / snom_wdg
                        ratios = vnom_wdg/voltage_scale_factor

                        # assign vnom_wdg to vnom for transformer
                        vnom[wdg_endNumber] = vnom_wdg

                        # Add vnom info to bus
                        data_math["settings"]["vbases_network"][string(bus)] = deepcopy(vnom_wdg/voltage_scale_factor)

                        # resistance computation
                        transf_star_impedance = get(transf_end_info[wdg_endNumber], "TransformerEndInfo.TransformerStarImpedance", Dict())
                        r_s[wdg_endNumber] = get(transf_end_info[wdg_endNumber], "TransformerEndInfo.r",
                                            get(transf_star_impedance, "TransformerStarImpedance.r", 0.0))
                        r_s[wdg_endNumber] = r_s[wdg_endNumber]/ratios^2

                        # reactance computation
                        x_sc[wdg_endNumber] = get(transf_end_info[wdg_endNumber], "TransformerEndInfo.x",
                        get(transf_star_impedance, "TransformerStarImpedance.x", 0.0))
                        # -- alternative computation of xsc using sc tests
                        if haskey(transf_end_info[wdg_endNumber], "TransformerEndInfo.EnergisedEndShortCircuitTests")
                        leak_impedance_wdg = transf_end_info[wdg_endNumber]["TransformerEndInfo.EnergisedEndShortCircuitTests"][1]["ShortCircuitTest.leakageImpedance"]
                        x_sc[wdg_endNumber] = (sqrt((leak_impedance_wdg / zbase)^2 - ((r_s[wdg_endNumber]*100.0 / zbase)*2)^2)/100.0)*zbase
                        end
                        x_sc[wdg_endNumber] = (x_sc[wdg_endNumber]/ratios^2)

                        # g_sh always with respect to wdg #1 always
                        if wdg_endNumber == 1
                            transf_end_noloadtest = get(transf_end_info[wdg_endNumber], "TransformerEndInfo.EnergisedEndNoLoadTests", [Dict()])
                            loss = get(transf_end_noloadtest[1], "NoLoadTest.loss", 0.0)
                            g_sh_tank =  (loss*snom_wdg)/zbase
                            exct_current = get(transf_end_noloadtest[1], "NoLoadTest.excitingCurrent", 0.0)
                            b_sh_tank = -((sqrt(abs((exct_current)^2 - (loss/(0.01*snom_wdg))^2)))/(100.0*zbase))
                            # data is measured externally, but we now refer it to the internal side
                            g_sh = g_sh_tank*ratios^2
                            b_sh = b_sh_tank*ratios^2
                        end

                        # wdgs configurations
                        wdg_conf = transf_end_info[wdg_endNumber]["TransformerEndInfo.connectionKind"] # extract wdg conf
                        if wdg_conf == "WindingConnection.Y" || wdg_conf == "WindingConnection.I" ||  wdg_conf == "WindingConnection.Yn"
                            wdgs_confs[wdg_endNumber] = WYE
                        elseif wdg_conf == "WindingConnection.D"
                            wdgs_confs[wdg_endNumber] = DELTA
                        else
                            @error("PowerTransformer ConnectionKind not supported yet!")
                        end

                        # Set RatioTapChanger in specific wdg
                        if haskey(wdgs_data[wdg_endNumber], "TransformerEnd.RatioTapChanger")

                            rtc_name = _extract_name(wdgs_data[wdg_endNumber]["TransformerEnd.RatioTapChanger"])
                            rtc_data = data_ravens["PowerSystemResource"]["TapChanger"]["RatioTapChanger"][rtc_name]

                            # tm_step
                            hstep = get(rtc_data, "TapChanger.highStep", 16)
                            lstep = get(rtc_data, "TapChanger.lowStep", -16)
                            step_dist = abs(hstep) + abs(lstep)
                            step_tap = 1/step_dist
                            tm_step[wdg_endNumber] = fill(step_tap, nphases)

                            # tm_set
                            step = get(rtc_data, "TapChanger.step", 1.0)    # Starting Tap changer position/step
                            tm_set[wdg_endNumber] = fill(step, nphases)

                            # tm_fix
                            ltcFlag = get(rtc_data, "TapChanger.ltcFlag", false)
                            if (ltcFlag == true)
                                tm_fix[wdg_endNumber] = zeros(Bool, nphases)
                            else
                                tm_fix[wdg_endNumber] = ones(Bool, nphases)
                            end

                            # tm_ub/tm_lb
                            neutralVoltPu = get(rtc_data, "TapChanger.neutralU", vnom[wdg_endNumber])/vnom[wdg_endNumber]
                            step_volt_increment = get(rtc_data, "RatioTapChanger.stepVoltageIncrement", 100.0)
                            volt_lb = neutralVoltPu + step_tap * (step_volt_increment/100.0) * lstep
                            volt_ub = neutralVoltPu + step_tap * (step_volt_increment/100.0) * hstep
                            tm_lb[wdg_endNumber] = fill(volt_lb, nphases)
                            tm_ub[wdg_endNumber] = fill(volt_ub, nphases)

                            # Regulator Control
                            if haskey(rtc_data, "TapChanger.TapChangerControl") && !all(tm_fix[wdg_endNumber])
                                reg_controls[wdg_endNumber] = true
                                reg_obj[wdg_endNumber] = Dict{String,Any}(
                                        "vreg" => fill(rtc_data["TapChanger.TapChangerControl"]["RegulatingControl.targetValue"], nphases),
                                        "band" =>  fill(rtc_data["TapChanger.TapChangerControl"]["RegulatingControl.targetDeadband"], nphases),
                                        "ptratio" => fill(rtc_data["TapChanger.ptRatio"], nphases),
                                        "ctprim" => fill(rtc_data["TapChanger.ctRating"], nphases),
                                        "r" => fill(rtc_data["TapChanger.TapChangerControl"]["TapChangerControl.lineDropR"], nphases),
                                        "x" => fill(rtc_data["TapChanger.TapChangerControl"]["TapChangerControl.lineDropX"], nphases)
                                    )
                            end

                        else # default
                            tm_set[wdg_endNumber] = fill(1.0, nphases)
                            tm_lb[wdg_endNumber] = fill(0.9, nphases)
                            tm_ub[wdg_endNumber] = fill(1.1, nphases)
                            tm_fix[wdg_endNumber] = ones(Bool, nphases)
                            tm_step[wdg_endNumber] = fill(1/32, nphases)
                        end

                    end

                    # convert x_sc from list of upper triangle elements to an explicit dict
                    y_sh = g_sh + im*b_sh
                    z_sc = Dict([(key, im*x_sc[i]) for (i,key) in enumerate([(i,j) for i in 1:nrw for j in i+1:nrw])])

                    # dimensions
                    dims = length(tm_set[1])

                    # init polarity
                    polarity = fill(1, nrw)

                    # Status
                    status = haskey(tanks[tank_id], "Equipment.inService") ? tanks[tank_id]["Equipment.inService"] : true
                    status = status == true ? 1 : 0

                    # Build loss model
                    transformer_t_bus_w = _build_loss_model!(data_math, name, to_map, r_s, z_sc, y_sh, connections[1]; nphases=dims, status=status)

                    # Mathematical model for transformer
                    for wdg_id in 1:1:nrw
                        # 2-WINDING TRANSFORMER

                        # correct polarity and connections
                        if wdg_id>1
                            if wdgs_confs[1] == DELTA && wdgs_confs[wdg_id] == WYE
                                polarity[wdg_id] = -1
                                connections[wdg_id] = _barrel_roll(connections[wdg_id][1:end], 1)
                            end
                            if wdgs_confs[1] == WYE && wdgs_confs[wdg_id] == DELTA
                                polarity[wdg_id] = -1
                                connections[wdg_id] = _barrel_roll(connections[wdg_id], -1)
                            end
                        end

                        # tank asset data
                        wdg_info = tank_asset_data["PowerTransformerInfo.TransformerTankInfos"][tank_asset_name]["TransformerTankInfo.TransformerEndInfos"]

                        # make virtual bus and mark it for reduction
                        tm_nom = wdgs_confs[wdg_id]==DELTA ? vnom[wdg_id]*sqrt(3)/voltage_scale_factor : vnom[wdg_id]/voltage_scale_factor

                        # Get correct f_node for winding
                        wdg_term = ravens_obj["ConductingEquipment.Terminals"][wdg_id]
                        f_node_wdgterm = _extract_name(wdg_term["Terminal.ConnectivityNode"])

                        # Transformer Object
                        transformer_2wa_obj = Dict{String,Any}(
                            "name"          => "_virtual_transformer.$name.$wdg_id.$(connections[wdg_id])",
                            "source_id"     => "_virtual_transformer.transformer.$name.$wdg_id.$(connections[wdg_id])",
                            "f_bus"         => data_math["bus_lookup"][f_node_wdgterm],
                            "t_bus"         => transformer_t_bus_w[wdg_id],
                            "tm_nom"        => tm_nom,
                            "f_connections" => connections[wdg_id],
                            "t_connections" => connections[1],
                            "configuration" => wdgs_confs[wdg_id],
                            "polarity"      => polarity[wdg_id],
                            "tm_set"        => tm_set[wdg_id],
                            "tm_fix"        => tm_fix[wdg_id],
                            "sm_ub"         => get(wdg_info[wdg_id], "TransformerEndInfo.emergencyS", Inf)/power_scale_factor,
                            "cm_ub"         => get(wdg_info[wdg_id], "TransformerEndInfo.ratedI", Inf),
                            "status"        => status,
                            "index"         => length(data_math["transformer"])+1
                        )

                        # RatioTapChanger
                        transformer_2wa_obj["tm_lb"] = tm_lb[wdg_id]
                        transformer_2wa_obj["tm_ub"] = tm_ub[wdg_id]
                        transformer_2wa_obj["tm_step"] = tm_step[wdg_id]

                        data_math["transformer"]["$(transformer_2wa_obj["index"])"] = transformer_2wa_obj

                        # Add Regulator Controls (only if flag is true)
                        if (reg_controls[wdg_id] == true)
                            data_math["transformer"]["$(transformer_2wa_obj["index"])"]["controls"] = reg_obj[wdg_id]
                        end

                        # TODO: Center-Tapped Transformers (3 Windings)
                        # if w==3 && eng_obj["polarity"][w]==-1 # identify center-tapped transformer and mark all secondary-side nodes as triplex by adding va_start
                        # end

                        push!(to_map, "transformer.$(transformer_2wa_obj["index"])")

                    end
                end
            end
        end
    end

    # Checks and calculates voltage bases for elements that do not have Voltage Bases.
    # TODO: Revise if this is the best way to be calculate vbases for missing elements
    data_math["settings"]["vbases_buses"] = calc_math_voltage_bases(data_math, data_math["settings"]["vbases_network"])[1] # [1] bus_vbase, [2] edge_vbase

end


"""
Converts ravens load components into mathematical load components.
"""
function _map_ravens2math_energy_consumer!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[], nw::Int=nw_id_default)
    energy_connections = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]["EnergyConnection"]
    power_scale_factor = data_math["settings"]["power_scale_factor"]
    voltage_scale_factor = data_math["settings"]["voltage_scale_factor"]
    voltage_scale_factor_sqrt3 = voltage_scale_factor * sqrt(3)

    for (name, ravens_obj) in get(energy_connections, "EnergyConsumer", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj_ravens("energy_consumer", name, ravens_obj, length(data_math["load"]) + 1; pass_props=pass_props)

        # Set the load bus based on connectivity node
        connectivity_node = _extract_name(ravens_obj["ConductingEquipment.Terminals"][1]["Terminal.ConnectivityNode"])
        math_obj["load_bus"] = data_math["bus_lookup"][connectivity_node]

        # Handle Load Response Characteristics
        load_response_characts = _extract_name(ravens_obj["EnergyConsumer.LoadResponse"])
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

        # Set voltage bounds for the bus connected
        bus_info = string(math_obj["load_bus"])
        bus_conn = data_math["bus"][bus_info]

        # Set the nominal voltage
        if haskey(ravens_obj, "ConductingEquipment.BaseVoltage")
            base_voltage_ref = _extract_name(ravens_obj["ConductingEquipment.BaseVoltage"])
            base_voltage = data_ravens["BaseVoltage"][base_voltage_ref]["BaseVoltage.nominalVoltage"]
            math_obj["vnom_kv"] = (base_voltage / voltage_scale_factor_sqrt3)
        else
            math_obj["vnom_kv"] = data_math["settings"]["vbases_buses"][bus_info]
        end

        if haskey(data_ravens["ConnectivityNode"][connectivity_node], "ConnectivityNode.OperationalLimitSet")
            op_limit_id = _extract_name(data_ravens["ConnectivityNode"][connectivity_node]["ConnectivityNode.OperationalLimitSet"])
            op_limits = data_ravens["OperationalLimitSet"][op_limit_id]["OperationalLimitSet.OperationalLimitValue"]

            # Loop through op. limits
            for lim in op_limits

                if haskey(data_ravens, "OperationalLimitType")
                    lim_type_name = _extract_name(lim["OperationalLimit.OperationalLimitType"])
                    lim_type = data_ravens["OperationalLimitType"][lim_type_name]["OperationalLimitType.direction"]
                else
                    lim_type = lim["OperationalLimit.OperationalLimitType"]["OperationalLimitType.direction"]
                end

                if lim_type == "OperationalLimitDirectionKind.high"
                    op_limit_max = lim["VoltageLimit.value"] / voltage_scale_factor_sqrt3
                elseif lim_type == "OperationalLimitDirectionKind.low"
                    op_limit_min = lim["VoltageLimit.value"] / voltage_scale_factor_sqrt3
                end
            end

        else
            op_limit_max = Inf
            op_limit_min = 0.0
        end

        # Handle phase-specific or three-phase connection
        nphases = 0
        if haskey(ravens_obj, "EnergyConsumer.EnergyConsumerPhase")
            connections = Vector{Int64}()
            for phase_info in ravens_obj["EnergyConsumer.EnergyConsumerPhase"]
                phase = _phase_map[phase_info["EnergyConsumerPhase.phase"]]
                phase_index = findfirst(==(phase), bus_conn["terminals"])
                bus_conn["vmax"][phase_index] = op_limit_max
                bus_conn["vmin"][phase_index] = op_limit_min
                push!(connections, phase)
            end
            math_obj["connections"] = connections
            nphases = length(math_obj["connections"])
        else
            nphases = length(bus_conn["terminals"])
            bus_conn["vmax"] = fill(op_limit_max, nphases)
            bus_conn["vmin"] = fill(op_limit_min, nphases)
            math_obj["connections"] = bus_conn["terminals"]
        end

        # Set p and q (w/ multinetwork support)
        if nw==0
            math_obj["pd"] = fill(get(ravens_obj, "EnergyConsumer.p", 0.0) / (power_scale_factor*nphases), nphases)
            math_obj["qd"] = fill(get(ravens_obj, "EnergyConsumer.q", 0.0) / (power_scale_factor*nphases), nphases)
        else
            # Get timeseries schedule
            if haskey(ravens_obj, "EnergyConsumer.LoadProfile")

                # get active (P) and reactive power (Q) of the load
                active_power = zeros(Float64, nphases)
                reactive_power = zeros(Float64, nphases)

                if haskey(ravens_obj, "EnergyConsumer.EnergyConsumerPhase")
                    for id in 1:nphases
                        phase_info = ravens_obj["EnergyConsumer.EnergyConsumerPhase"][id]
                        active_power[id] = phase_info["EnergyConsumerPhase.p"]
                        reactive_power[id] = phase_info["EnergyConsumerPhase.q"]
                    end
                else
                    active_power = fill(get(ravens_obj, "EnergyConsumer.p", 0.0) / (power_scale_factor*nphases), nphases)
                    reactive_power = fill(get(ravens_obj, "EnergyConsumer.p", 0.0) / (power_scale_factor*nphases), nphases)
                end

                schdl_name = _extract_name(ravens_obj["EnergyConsumer.LoadProfile"])
                schdl = data_ravens["BasicIntervalSchedule"][schdl_name]

                # units and multiplier modifiers
                if haskey(schdl, "BasicIntervalSchedule.value1Multiplier")
                    value1_multiplier = _multipliers_map[schdl["BasicIntervalSchedule.value1Multiplier"]]
                else
                    value1_multiplier = 1.0
                end

                if haskey(schdl, "BasicIntervalSchedule.value2Multiplier")
                    value2_multiplier = _multipliers_map[schdl["BasicIntervalSchedule.value2Multiplier"]]
                else
                    value2_multiplier = 1.0
                end

                if haskey(schdl, "BasicIntervalSchedule.value1Unit")
                    unit_symbol = schdl["BasicIntervalSchedule.value1Unit"]
                    value1_unit = lowercase(unit_symbol[findfirst(isequal('.'), unit_symbol) + 1:end])
                    if value1_unit == "w"
                        math_obj["pd"] = fill(get(schdl["EnergyConsumerSchedule.RegularTimePoints"][nw], "RegularTimePoint.value1", 0.0) * value1_multiplier / power_scale_factor, nphases)
                    end
                    if value1_unit == "var"
                        math_obj["qd"] = fill(get(schdl["EnergyConsumerSchedule.RegularTimePoints"][nw], "RegularTimePoint.value1", 0.0) * value1_multiplier / power_scale_factor, nphases)
                    end
                end

                if haskey(schdl, "BasicIntervalSchedule.value2Unit")
                    unit_symbol = schdl["BasicIntervalSchedule.value2Unit"]
                    value2_unit = lowercase(unit_symbol[findfirst(isequal('.'), unit_symbol) + 1:end])
                    if value2_unit == "w"
                        math_obj["pd"] = fill(get(schdl["EnergyConsumerSchedule.RegularTimePoints"][nw], "RegularTimePoint.value2", 0.0) * value2_multiplier / power_scale_factor, nphases)
                    end
                    if value2_unit == "var"
                        math_obj["qd"] = fill(get(schdl["EnergyConsumerSchedule.RegularTimePoints"][nw], "RegularTimePoint.value2", 0.0) * value2_multiplier / power_scale_factor, nphases)
                    end
                end

                # Multipliers instead of actual values
                if !haskey(schdl, "BasicIntervalSchedule.value1Unit")
                    math_obj["pd"] = get(schdl["EnergyConsumerSchedule.RegularTimePoints"][nw], "RegularTimePoint.value1", 1.0) .* active_power ./ power_scale_factor
                    math_obj["qd"] = get(schdl["EnergyConsumerSchedule.RegularTimePoints"][nw], "RegularTimePoint.value1", 1.0) .* reactive_power ./ power_scale_factor
                end

            else
                @error("No timeseries, load forecast or multinetwork information found!")
            end
        end


        # Set the configuration
        # TODO: ADD: "PhaseShuntConnectionKind.I", "PhaseShuntConnectionKind.G"
        config_map = Dict("PhaseShuntConnectionKind.Y" => WYE, "PhaseShuntConnectionKind.D" => DELTA, "PhaseShuntConnectionKind.Yn" => WYE)
        config = get(config_map, ravens_obj["EnergyConsumer.phaseConnection"], nothing)
        if config !== nothing
            math_obj["configuration"] = config
        else
            @error("Configuration of load $(name) is not supported.")
        end

        # Correct (if needed) single-phase DELTA connections
        if (math_obj["configuration"] == DELTA) && (nphases == 1)
            math_obj["configuration"] = WYE
            @warn "EnergyConsumer (load): $name has DELTA configuration but only 1 connection (phase). DELTA configurations must have at least 2 or 3 connections!. EnergyConsumer converted to WYE connection."
        end

        # Set status, dispatchable flag, and index
        math_obj["status"] = haskey(ravens_obj, "Equipment.inService") ? ravens_obj["Equipment.inService"] : true
        math_obj["status"] = math_obj["status"] == true ? 1 : 0
        math_obj["dispatchable"] = 0
        data_math["load"]["$(math_obj["index"])"] = math_obj

        # Handle grounding
        if !(haskey(bus_conn, "grounded"))
            if ravens_obj["EnergyConsumer.grounded"] == true
                bus_conn["grounded"] = ones(Bool, length(math_obj["connections"]))
            else
                bus_conn["grounded"] = zeros(Bool, length(math_obj["connections"]))
            end
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
function _map_ravens2math_energy_source!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[], nw::Int=nw_id_default)
    energy_connections = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]["EnergyConnection"]
    voltage_scale_factor = data_math["settings"]["voltage_scale_factor"]
    voltage_scale_factor_sqrt3 = voltage_scale_factor * sqrt(3)

    for (name, ravens_obj) in get(energy_connections, "EnergySource", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj_ravens("energy_source", name, ravens_obj, length(data_math["gen"]) + 1; pass_props=pass_props)
        math_obj["name"] = "_virtual_gen.energy_source.$name"

        # Get connectivity node info (bus info)
        connectivity_node = _extract_name(ravens_obj["ConductingEquipment.Terminals"][1]["Terminal.ConnectivityNode"])
        gen_bus = data_math["bus_lookup"][connectivity_node]
        math_obj["gen_bus"] = gen_bus
        bus_conn = data_math["bus"][string(gen_bus)]
        bus_conn["bus_type"] = 3  # Set bus type to PV bus

        # Handle phase-specific or three-phase connection
        connections = Vector{Int64}()
        nconductors = length(get(ravens_obj, "EnergySource.EnergySourcePhase", bus_conn["terminals"]))

        if haskey(ravens_obj, "EnergySource.EnergySourcePhase")
            for phase_info in ravens_obj["EnergySource.EnergySourcePhase"]
                phase = _phase_map[phase_info["EnergySourcePhase.phase"]]
                push!(connections, phase)
            end
            math_obj["connections"] = connections
        else
            # Terminal Phases
            if haskey(ravens_obj["ConductingEquipment.Terminals"][1], "Terminal.phases")
                phasecode = ravens_obj["ConductingEquipment.Terminals"][1]["Terminal.phases"]
                math_obj["connections"] = _phasecode_map[phasecode]
            else
                math_obj["connections"] = bus_conn["terminals"]
            end
        end

        # Generator status and configuration
        math_obj["gen_status"] = haskey(ravens_obj, "Equipment.inService") ? ravens_obj["Equipment.inService"] : true
        math_obj["gen_status"] = math_obj["gen_status"] == true ? 1 : 0

        math_obj["configuration"] = get(ravens_obj, "EnergySource.connectionKind", WYE)

        # Set the nominal voltage
        if haskey(ravens_obj, "ConductingEquipment.BaseVoltage")
            base_voltage_ref = _extract_name(ravens_obj["ConductingEquipment.BaseVoltage"])
            vnom = data_ravens["BaseVoltage"][base_voltage_ref]["BaseVoltage.nominalVoltage"] / sqrt(nconductors)
            data_math["settings"]["vbases_default"][connectivity_node] = vnom / voltage_scale_factor
        else
            vnom = ravens_obj["EnergySource.nominalVoltage"] / sqrt(nconductors)
            data_math["settings"]["vbases_default"][connectivity_node] = vnom / voltage_scale_factor
        end

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


        rs = fill(get(ravens_obj, "EnergySource.r", zeros(1, 1)), nconductors, nconductors)
        xs = fill(get(ravens_obj, "EnergySource.x", zeros(1, 1)), nconductors, nconductors)

        # Check for impedance and adjust bus type if necessary
        map_to = "gen.$(math_obj["index"])"
        if !all(isapprox.(rs, 0)) && !all(isapprox.(xs, 0))
            bus_conn["bus_type"] = 1  # Virtual bus becomes the new slack bus

            bus_obj = Dict(
                "bus_i" => length(data_math["bus"]) + 1,
                "index" => length(data_math["bus"]) + 1,
                "terminals" => math_obj["connections"],
                "grounded" => zeros(Bool, nphases),
                "name" => "_virtual_bus.energy_source.$name",
                "bus_type" => math_obj["gen_status"] == 0 ? 4 : math_obj["control_mode"] == Int(ISOCHRONOUS) ? 3 : 2,
                "vm" => fill(ravens_obj["EnergySource.voltageMagnitude"] / voltage_scale_factor_sqrt3, nphases),
                "va" => rad2deg.(_wrap_to_pi.([-2 * π / nphases * (i - 1) + get(ravens_obj, "EnergySource.voltageAngle", 0.0) for i in 1:nphases])),
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
            vm_lb = math_obj["control_mode"] == Int(ISOCHRONOUS) ? fill(ravens_obj["EnergySource.voltageMagnitude"] / voltage_scale_factor_sqrt3, nphases) : get(ravens_obj, "EnergySource.vMin", fill(1.0, nphases))
            vm_ub = math_obj["control_mode"] == Int(ISOCHRONOUS) ? fill(ravens_obj["EnergySource.voltageMagnitude"] / voltage_scale_factor_sqrt3, nphases) : get(ravens_obj, "EnergySource.vMax", fill(1.0, nphases))

            data_math["bus"]["$gen_bus"]["vmin"] = [vm_lb..., fill(0.0, nconductors - nphases)...]
            data_math["bus"]["$gen_bus"]["vmax"] = [vm_ub..., fill(Inf, nconductors - nphases)...]
            data_math["bus"]["$gen_bus"]["vm"] = fill(ravens_obj["EnergySource.voltageMagnitude"] / voltage_scale_factor_sqrt3, nphases)
            data_math["bus"]["$gen_bus"]["va"] = rad2deg.(_wrap_to_pi.([-2 * π / nphases * (i - 1) + get(ravens_obj, "EnergySource.voltageAngle", 0.0) for i in 1:nphases]))
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



"converts engineering generators into mathematical generators"
function _map_ravens2math_rotating_machine!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[], nw::Int=nw_id_default)
    energy_connections = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]["EnergyConnection"]

    if haskey(energy_connections, "RegulatingCondEq")

        regulating_cond_eq = energy_connections["RegulatingCondEq"]
        power_scale_factor = data_math["settings"]["power_scale_factor"]
        voltage_scale_factor = data_math["settings"]["voltage_scale_factor"]

        for (name, ravens_obj) in get(regulating_cond_eq, "RotatingMachine", Dict{Any,Dict{String,Any}}())

            math_obj = _init_math_obj_ravens("rotating_machine", name, ravens_obj, length(data_math["gen"])+1; pass_props=pass_props)

            # Connections/phases obtained from Terminals
            connections = _phasecode_map[get(ravens_obj["ConductingEquipment.Terminals"][1], "Terminal.phases", "PhaseCode.ABC")]

            nconductors = length(connections)
            math_obj["connections"] = connections

            connectivity_node = _extract_name(ravens_obj["ConductingEquipment.Terminals"][1]["Terminal.ConnectivityNode"])
            math_obj["gen_bus"] = data_math["bus_lookup"][connectivity_node]
            math_obj["gen_status"] = get(ravens_obj, "Equipment.inService", true)
            math_obj["gen_status"] = status = math_obj["gen_status"] == true ? 1 : 0

            # TODO: control mode do not exist in the RAVENS-CIM (Need to be added)
            math_obj["control_mode"] = control_mode = Int(get(ravens_obj, "control_mode", FREQUENCYDROOP))

            # Set Pmax and Pmin for generator
            if haskey(ravens_obj, "RotatingMachine.GeneratingUnit")
                math_obj["pmin"] = ((get(ravens_obj["RotatingMachine.GeneratingUnit"], "GeneratingUnit.minOperatingP", 0) * ones(nconductors)) ./ nconductors)./(power_scale_factor)
                math_obj["pmax"] = ((get(ravens_obj["RotatingMachine.GeneratingUnit"], "GeneratingUnit.maxOperatingP", Inf) * ones(nconductors)) ./ nconductors)./(power_scale_factor)
            else
                math_obj["pmin"] = (zeros(nconductors) ./ nconductors)./(power_scale_factor)
                math_obj["pmax"] = ((get(ravens_obj, "RotatingMachine.ratedS", Inf) * ones(nconductors)) ./ nconductors)./(power_scale_factor)
            end

            # Set bus type
            bus_type = data_math["bus"]["$(math_obj["gen_bus"])"]["bus_type"]
            data_math["bus"]["$(math_obj["gen_bus"])"]["bus_type"] = _compute_bus_type(bus_type, status, control_mode)

           # Set the nominal voltage
           bus_conn =  data_math["bus"]["$(math_obj["gen_bus"])"]
           if haskey(ravens_obj, "ConductingEquipment.BaseVoltage")
               base_voltage_ref = _extract_name(ravens_obj["ConductingEquipment.BaseVoltage"])
               nominal_voltage = data_ravens["BaseVoltage"][base_voltage_ref]["BaseVoltage.nominalVoltage"]
               base_voltage =  nominal_voltage / sqrt(nconductors)
               math_obj["vbase"] =  base_voltage / voltage_scale_factor
           else
                math_obj["vbase"] = data_math["settings"]["vbases_buses"][string(math_obj["gen_bus"])]
           end

            if control_mode == Int(ISOCHRONOUS) && status == 1
                data_math["bus"]["$(math_obj["gen_bus"])"]["vm"] = ((get(ravens_obj, "RotatingMachine.ratedU", nominal_voltage))/nominal_voltage)* ones(nconductors)
                data_math["bus"]["$(math_obj["gen_bus"])"]["vmax"] = ((get(ravens_obj, "RotatingMachine.ratedU", nominal_voltage))/nominal_voltage)* ones(nconductors)
                data_math["bus"]["$(math_obj["gen_bus"])"]["vmin"] = ((get(ravens_obj, "RotatingMachine.ratedU", nominal_voltage))/nominal_voltage)* ones(nconductors)
                data_math["bus"]["$(math_obj["gen_bus"])"]["va"] = [0.0, -120, 120, zeros(length(data_math["bus"]["$(math_obj["gen_bus"])"]) - 3)...][data_math["bus"]["$(math_obj["gen_bus"])"]["terminals"]]
            end

            # Set min and max Q
            if haskey(ravens_obj, "RotatingMachine.minQ")
                math_obj["qmin"] = ((ravens_obj["RotatingMachine.minQ"] * ones(nconductors)) ./ nconductors)./(power_scale_factor)
            elseif haskey(ravens_obj, "SynchronousMachine.minQ")
                math_obj["qmin"] = ((ravens_obj["SynchronousMachine.minQ"] * ones(nconductors)) ./ nconductors)./(power_scale_factor)
            else
                math_obj["qmin"] = fill(-Inf, nconductors)
            end

            if haskey(ravens_obj, "RotatingMachine.maxQ")
                math_obj["qmax"] = ((ravens_obj["RotatingMachine.maxQ"] * ones(nconductors)) ./ nconductors)./(power_scale_factor)
            elseif haskey(ravens_obj, "SynchronousMachine.maxQ")
                math_obj["qmax"] = ((ravens_obj["SynchronousMachine.maxQ"] * ones(nconductors)) ./ nconductors)./(power_scale_factor)
            else
                math_obj["qmax"] = fill(Inf, nconductors)
            end

            # Set pg and qg
            math_obj["pg"] = (get(ravens_obj, "RotatingMachine.p", 0.0) * ones(nconductors) ./ nconductors)./(power_scale_factor)
            math_obj["qg"] = (get(ravens_obj, "RotatingMachine.q", 0.0) * ones(nconductors) ./ nconductors)./(power_scale_factor)

            # TODO: add a polynomial parameters to be added to gen cost
            _add_gen_cost_model!(math_obj, ravens_obj)

            # TODO: configuration for generators is not available on CIM (yet)
            math_obj["configuration"] = get(ravens_obj, "configuration", WYE)

            # Set index
            data_math["gen"]["$(math_obj["index"])"] = math_obj

            push!(data_math["map"], Dict{String,Any}(
                "from" => name,
                "to" => "gen.$(math_obj["index"])",
                "unmap_function" => "_map_math2ravens_rotating_machine!",
            ))
        end

    end

end


"converts ravens power_electronics units such as PVs and Batteries into mathematical components"
function _map_ravens2math_power_electronics!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[], nw::Int=nw_id_default)
    energy_connections = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]["EnergyConnection"]

    if haskey(energy_connections, "RegulatingCondEq")

        regulating_cond_eq = energy_connections["RegulatingCondEq"]
        power_scale_factor = data_math["settings"]["power_scale_factor"]
        voltage_scale_factor = data_math["settings"]["voltage_scale_factor"]

        for (name, ravens_obj) in get(regulating_cond_eq, "PowerElectronicsConnection", Dict{Any,Dict{String,Any}}())

            # Get type of PowerElectronicsUnit
            pec_type = get(ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"], "Ravens.cimObjectType", "")

            if (pec_type == "PhotoVoltaicUnit")

                math_obj = _init_math_obj_ravens("photovoltaic_unit", name, ravens_obj, length(data_math["gen"])+1; pass_props=pass_props)

                # Connections/phases
                connections = _phasecode_map[get(ravens_obj["ConductingEquipment.Terminals"][1], "Terminal.phases", "PhaseCode.ABC")]

                nconductors = length(connections)
                math_obj["connections"] = connections

                connectivity_node = _extract_name(ravens_obj["ConductingEquipment.Terminals"][1]["Terminal.ConnectivityNode"])
                math_obj["gen_bus"] = data_math["bus_lookup"][connectivity_node]
                math_obj["gen_status"] = get(ravens_obj, "Equipment.inService", true)
                math_obj["gen_status"] = status = math_obj["gen_status"] == true ? 1 : 0

                # TODO: control mode do not exist in the RAVENS-CIM (Need to be added)
                math_obj["control_mode"] = control_mode = Int(get(ravens_obj, "control_mode", FREQUENCYDROOP))

                # Set bus type
                bus_type = data_math["bus"]["$(math_obj["gen_bus"])"]["bus_type"]
                data_math["bus"]["$(math_obj["gen_bus"])"]["bus_type"] = _compute_bus_type(bus_type, status, control_mode)

                # Set the nominal voltage
                base_voltage_ref = _extract_name(ravens_obj["ConductingEquipment.BaseVoltage"])
                nominal_voltage = data_ravens["BaseVoltage"][base_voltage_ref]["BaseVoltage.nominalVoltage"]
                base_voltage =  nominal_voltage / sqrt(nconductors)
                math_obj["vbase"] =  base_voltage / voltage_scale_factor

                if control_mode == Int(ISOCHRONOUS) && status == 1
                    data_math["bus"]["$(math_obj["gen_bus"])"]["vm"] = ((get(ravens_obj, "PowerElectronicsConnection.ratedU", nominal_voltage))/nominal_voltage)* ones(nconductors)
                    data_math["bus"]["$(math_obj["gen_bus"])"]["vmax"] = ((get(ravens_obj, "PowerElectronicsConnection.ratedU", nominal_voltage))/nominal_voltage)* ones(nconductors)
                    data_math["bus"]["$(math_obj["gen_bus"])"]["vmin"] = ((get(ravens_obj, "PowerElectronicsConnection.ratedU", nominal_voltage))/nominal_voltage)* ones(nconductors)
                    data_math["bus"]["$(math_obj["gen_bus"])"]["va"] = [0.0, -120, 120, zeros(length(data_math["bus"]["$(math_obj["gen_bus"])"]) - 3)...][data_math["bus"]["$(math_obj["gen_bus"])"]["terminals"]]
                    data_math["bus"]["$(math_obj["gen_bus"])"]["bus_type"] = 3
                end

                # Set vg
                for (fr_k, to_k) in [("PowerElectronicsConnection.ratedU", "vg")]
                    if haskey(ravens_obj, fr_k)
                        math_obj[to_k] = (ravens_obj[fr_k]/nominal_voltage)*ones(nconductors)/voltage_scale_factor
                    end
                end

                # TODO: configuration for generators is not available on CIM (yet)
                math_obj["configuration"] = get(ravens_obj, "configuration", WYE)


                # TODO: refactor the calculation of N when connections and configuration issues are solved.
                N = math_obj["configuration"]==DELTA && length(connections)==1 ? 1 : _infer_int_dim(connections,  math_obj["configuration"], false) # if solar is delta-connected to triplex node, N can be equal to 1

                # Set pmax
                if !haskey(ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"], "PowerElectronicsUnit.maxP")
                    math_obj["pmax"] = ((get(ravens_obj, "PowerElectronicsConnection.ratedS", Inf) * ones(nconductors)) ./ nconductors)./(power_scale_factor)
                else
                    math_obj["pmax"] = ((get(ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"], "PowerElectronicsUnit.maxP", Inf) * ones(nconductors)) ./ nconductors)./(power_scale_factor)
                end
                # Set pmin
                math_obj["pmin"] = ((get(ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"], "PowerElectronicsUnit.minP", 0) * ones(nconductors)) ./ nconductors)./(power_scale_factor)
                # Set qmin
                math_obj["qmin"] = ((get(ravens_obj, "PowerElectronicsConnection.minQ", -0) * ones(nconductors)) ./ nconductors)./(power_scale_factor)
                # Set qmax
                math_obj["qmax"] = ((get(ravens_obj, "PowerElectronicsConnection.maxQ", 0) * ones(nconductors)) ./ nconductors)./(power_scale_factor)


                # Set pg and qg
                math_obj["pg"] = (get(ravens_obj, "PowerElectronicsConnection.p", 0.0) * ones(nconductors) ./ nconductors)./(power_scale_factor)
                math_obj["qg"] = (get(ravens_obj, "PowerElectronicsConnection.q", 0.0) * ones(nconductors) ./ nconductors)./(power_scale_factor)

                # TODO: add a polynomial parameters to be added to gen cost
                _add_gen_cost_model!(math_obj, ravens_obj)

                # Set index
                data_math["gen"]["$(math_obj["index"])"] = math_obj

                push!(data_math["map"], Dict{String,Any}(
                    "from" => name,
                    "to" => "gen.$(math_obj["index"])",
                    "unmap_function" => "_map_math2ravens_photovoltaic_unit!",
                ))

            elseif (pec_type == "BatteryUnit")

                math_obj = _init_math_obj_ravens("storage", name, ravens_obj, length(data_math["storage"])+1; pass_props=pass_props)

                # TODO: connections/phases do not exist in the RAVENS-CIM (Need to be added) - should come from terminals
                connections = [1, 2, 3] # TODO
                nconductors = length(connections)
                math_obj["connections"] = connections

                # Set the bus
                connectivity_node = _extract_name(ravens_obj["ConductingEquipment.Terminals"][1]["Terminal.ConnectivityNode"])
                math_obj["storage_bus"] = data_math["bus_lookup"][connectivity_node]
                math_obj["status"] = get(ravens_obj, "Equipment.inService", true)
                math_obj["status"] = status = math_obj["status"] == true ? 1 : 0

                # TODO: configuration for generators is not available on CIM (yet)
                math_obj["configuration"] = get(ravens_obj, "configuration", WYE)

                # Set battery parameters
                math_obj["energy"] = ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"]["BatteryUnit.storedE"]/power_scale_factor

                if !haskey(ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"], "BatteryUnit.BatteryUnitEfficiency")
                    math_obj["energy_rating"] = ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"]["BatteryUnit.ratedE"]/power_scale_factor
                    math_obj["charge_efficiency"] = 100.0
                    math_obj["discharge_efficiency"] = 100.0
                    math_obj["p_loss"] = 0
                    math_obj["q_loss"] = 0
                else
                    math_obj["energy_rating"] = ((get(ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"]["BatteryUnit.BatteryUnitEfficiency"], "BatteryUnitEfficiency.limitEnergy", 100.0)/100.0)*ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"]["BatteryUnit.ratedE"])/power_scale_factor
                    math_obj["charge_efficiency"] = get(ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"], "BatteryUnitEfficiency.efficiencyCharge", 100.0) / 100.0
                    math_obj["discharge_efficiency"] = get(ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"], "BatteryUnitEfficiency.efficiencyDischarge", 100.0) / 100.0
                    # TODO: These are still missing from the RAVENS Schema
                    math_obj["p_loss"] = get(ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"], "BatteryUnitEfficiency.idlingActivePower", 0)./(power_scale_factor)
                    math_obj["q_loss"] = get(ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"], "BatteryUnitEfficiency.idlingReactivePower", 0)./(power_scale_factor)
                end

                if !haskey(ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"], "PowerElectronicsUnit.maxP")
                    math_obj["charge_rating"] = (get(ravens_obj, "PowerElectronicsConnection.ratedS", Inf))./(power_scale_factor)
                    math_obj["discharge_rating"] = math_obj["charge_rating"]
                else
                    math_obj["charge_rating"] = -(get(ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"], "PowerElectronicsUnit.minP", Inf))./(power_scale_factor)
                    math_obj["discharge_rating"] = (get(ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"], "PowerElectronicsUnit.maxP", Inf))./(power_scale_factor)
                end

                math_obj["thermal_rating"] = get(ravens_obj, "PowerElectronicsConnection.ratedS", Inf)/power_scale_factor

                math_obj["qmin"] = (get(ravens_obj, "PowerElectronicsConnection.minQ", -math_obj["discharge_rating"]*power_scale_factor))./(power_scale_factor)
                math_obj["qmax"] = (get(ravens_obj, "PowerElectronicsConnection.maxQ", math_obj["charge_rating"]*power_scale_factor))./(power_scale_factor)

                # TODO: verify that these CIM terms are equivalent to the needed values.
                math_obj["r"] = get(ravens_obj, "PowerElectronicsConnection.r", 0)
                math_obj["x"] = get(ravens_obj, "PowerElectronicsConnection.x", 0)

                # TODO: control mode do not exist in the RAVENS-CIM (Need to be added)
                math_obj["control_mode"] = control_mode = Int(get(ravens_obj, "control_mode", FREQUENCYDROOP))

                # Set the ps and qs
                math_obj["ps"] = (get(ravens_obj, "PowerElectronicsConnection.p", 0.0))./(power_scale_factor)
                math_obj["qs"] = (get(ravens_obj, "PowerElectronicsConnection.q", 0.0))./(power_scale_factor)

                # Set bus type
                bus_type = data_math["bus"]["$(math_obj["storage_bus"])"]["bus_type"]
                data_math["bus"]["$(math_obj["storage_bus"])"]["bus_type"] = _compute_bus_type(bus_type, status, control_mode)

                if control_mode == Int(ISOCHRONOUS) && math_obj["status"] == 1
                    data_math["bus"]["$(math_obj["storage_bus"])"]["va"] = [0.0, -120, 120, zeros(length(data_math["bus"]["$(math_obj["storage_bus"])"]) - 3)...][data_math["bus"]["$(math_obj["storage_bus"])"]["terminals"]]
                end

                data_math["storage"]["$(math_obj["index"])"] = math_obj

                push!(data_math["map"], Dict{String,Any}(
                    "from" => name,
                    "to" => "storage.$(math_obj["index"])",
                    "unmap_function" => "_map_math2ravens_battery_unit!",
                ))
            end
        end
    end
end


"converts ravens switches into mathematical switches and (if neeed) impedance branches to represent loss model"
function _map_ravens2math_switch!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[], nw::Int=nw_id_default)
    conducting_equipment = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]

    for (name, ravens_obj) in get(conducting_equipment, "Switch", Dict{Any,Dict{String,Any}}())

        math_obj = _init_math_obj_ravens("switch", name, ravens_obj, length(data_math["switch"])+1; pass_props=pass_props)

        # Terminals and phases
        terminals = ravens_obj["ConductingEquipment.Terminals"]

        # Loop through terminals
        f_conns = []
        t_conns = []
        for term in terminals
            if haskey(term, "Terminal.phases")
                phasecode = term["Terminal.phases"]
                if term["ACDCTerminal.sequenceNumber"] == 1
                    f_conns = _phasecode_map[phasecode]
                else
                    t_conns = _phasecode_map[phasecode]
                end
            else
                f_conns = [1, 2, 3]
                t_conns = [1, 2, 3]
            end
        end

        # Verify connections are correct.
        @assert f_conns == t_conns "f_conns are not equal to t_conns!. Revise connections/phases in Switch terminals"

        math_obj["f_connections"] = f_conns
        math_obj["t_connections"] = t_conns

        # Phases
        nphases = length(f_conns)

        # Connectivity Nodes
        f_node = _extract_name(terminals[1]["Terminal.ConnectivityNode"])
        t_node = _extract_name(terminals[2]["Terminal.ConnectivityNode"])
        math_obj["f_bus"] = data_math["bus_lookup"][f_node]
        math_obj["t_bus"] = data_math["bus_lookup"][t_node]

        # Add vmin/vmax/terminals info to fbus and tbus if missing
        for bus in [math_obj["f_bus"], math_obj["t_bus"]]
            bus_data = data_math["bus"][string(bus)]
            if !(haskey(bus_data, "terminals")) || (length(bus_data["terminals"]) < length(f_conns))
                bus_data["terminals"] = f_conns
                bus_data["vmin"] = fill(0.0, nphases)
                bus_data["vmax"] = fill(Inf, nphases)
                bus_data["grounded"] = zeros(Bool, nphases)
            end
        end

        # Status
        math_obj["status"] = get(ravens_obj, "Equipment.inService", true)
        math_obj["status"] = status = math_obj["status"] == true ? 1 : 0

        # State
        sw_state = get(ravens_obj, "Switch.open", false)
        sw_state = sw_state == false ? CLOSED : OPEN
        math_obj["state"] = Int(sw_state)

        # TODO: Dispatchable
        math_obj["dispatchable"] = Int(get(ravens_obj, "dispatchable", YES))

        # Current and Power Limits
        if haskey(ravens_obj, "PowerSystemResource.AssetDatasheet")
            swinfo_name = _extract_name(ravens_obj["PowerSystemResource.AssetDatasheet"])
            swinfo_data = data_ravens["AssetInfo"]["SwitchInfo"][swinfo_name]
            math_obj["current_rating"] = fill(get(swinfo_data, "SwitchInfo.ratedCurrent", Inf), nphases)
            math_obj["sm_ub"] = math_obj["current_rating"] .* get(swinfo_data, "SwitchInfo.ratedVoltage", Inf)
        else
            math_obj["current_rating"] = fill(get(ravens_obj, "Switch.ratedCurrent", Inf), nphases)
            math_obj["sm_ub"] = math_obj["current_rating"] .* get(ravens_obj, "Switch.ratedVoltage", Inf)
        end

        # TODO: not found on CIM - kron reductions
        for (f_key, t_key) in [("cm_ub_b", "c_rating_b"), ("cm_ub_c", "c_rating_c"), ("sm_ub_b", "rate_b"), ("sm_ub_c", "rate_c")]
            math_obj[t_key] = haskey(ravens_obj, f_key) ? fill(ravens_obj[f_key], nphases) : fill(Inf, nphases)
        end

        # Map index
        map_to = "switch.$(math_obj["index"])"

        # Push Mapping
        data_math["switch"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => map_to,
            "unmap_function" => "_map_math2ravens_switch!",
        ))

    end

end


"converts ravens generic shunt components into mathematical shunt components"
function _map_ravens2math_shunt_compensator!(data_math::Dict{String,<:Any}, data_ravens::Dict{String,<:Any}; pass_props::Vector{String}=String[], nw::Int=nw_id_default)
    energy_connections = data_ravens["PowerSystemResource"]["Equipment"]["ConductingEquipment"]["EnergyConnection"]

    if haskey(energy_connections, "RegulatingCondEq")
        regulating_cond_eq = energy_connections["RegulatingCondEq"]
        power_scale_factor = data_math["settings"]["power_scale_factor"]
        voltage_scale_factor = data_math["settings"]["voltage_scale_factor"]

        for (name, ravens_obj) in get(regulating_cond_eq, "ShuntCompensator", Dict{Any,Dict{String,Any}}())

            math_obj = _init_math_obj("shunt", name, ravens_obj, length(data_math["shunt"])+1; pass_props=pass_props)

            # Get connectivity node info (bus info)
            connectivity_node = _extract_name(ravens_obj["ConductingEquipment.Terminals"][1]["Terminal.ConnectivityNode"])
            math_obj["shunt_bus"] = data_math["bus_lookup"][connectivity_node]

            # Status
            status = haskey(ravens_obj, "Equipment.inService") ? ravens_obj["Equipment.inService"] : true
            math_obj["status"] = status == true ? 1 : 0

            # Connections/phases obtained from Terminals
            connections = _phasecode_map[get(ravens_obj["ConductingEquipment.Terminals"][1], "Terminal.phases", "PhaseCode.ABC")]

            math_obj["connections"] = connections
            terminals = connections

            # TODO: dispatchable
            math_obj["dispatchable"] = 0

            # bs - TODO: make sure b matrix is being calculated correctly
            b = ravens_obj["LinearShuntCompensator.bPerSection"]
            B = _calc_shunt_admittance_matrix(terminals, b)
            math_obj["bs"] = B

            # gs
            if haskey(ravens_obj, "LinearShuntCompensator.gPerSection")
                g = ravens_obj["LinearShuntCompensator.gPerSection"]
                G = _calc_shunt_admittance_matrix(terminals, g)
                math_obj["gs"] = G
            else
                math_obj["gs"] = zeros(size(math_obj["bs"]))
            end

            # Index
            data_math["shunt"]["$(math_obj["index"])"] = math_obj

            # TODO: Add CapControl
            # .....

            push!(data_math["map"], Dict{String,Any}(
                "from" => name,
                "to" => "shunt.$(math_obj["index"])",
                "unmap_function" => "_map_math2ravens_shunt!",
            ))
        end
    end
end
