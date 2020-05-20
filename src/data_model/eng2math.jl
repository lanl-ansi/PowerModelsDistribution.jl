import LinearAlgebra: diagm

"items that are mapped one-to-one from engineering to math models"
const _1to1_maps = Dict{String,Vector{String}}(
    "bus" => ["vm", "va", "terminals", "phases", "neutral", "vm_pn_lb", "vm_pn_ub", "vm_pp_lb", "vm_pp_ub", "vm_ng_ub", "dss"],
    "line" => ["f_connections", "t_connections", "source_id", "dss"],
    "transformer" => ["f_connections", "t_connections", "source_id", "dss"],
    "switch" => ["status", "f_connections", "t_connections", "source_id", "dss"],
    "shunt" => ["status", "dispatchable", "gs", "bs", "connections", "source_id", "dss"],
    "load" => ["model", "configuration", "connections", "dispatchable", "status", "source_id", "dss"],
    "generator" => ["pg", "qg", "vg", "configuration", "connections", "source_id", "dss"],
    "solar" => ["pg", "qg", "configuration", "connections", "source_id", "dss"],
    "storage" => ["status", "energy", "ps", "qs", "connections", "source_id", "dss"],
    "voltage_source" => ["source_id", "dss"],
)

"list of nodal type elements in the engineering model"
const _eng_node_elements = Vector{String}([
    "load", "shunt", "generator", "solar", "storage", "voltage_source"
])

"list of edge type elements in the engineering model"
const _eng_edge_elements = Vector{String}([
    "line", "switch", "transformer"
])

"list of nodal type elements in the engineering model"
const _math_node_elements = Vector{String}([
    "load", "shunt", "gen", "storage"
])

"list of edge type elements in the engineering model"
const _math_edge_elements = Vector{String}([
    "branch", "switch", "transformer", "dcline"
])

"list of multinetwork keys that belong at the root level"
const _pmd_math_global_keys = Set{String}([
    "data_model", "per_unit", "name", "settings", "map", "bus_lookup"
])


"converts a engineering multinetwork to a math multinetwork"
function _map_eng2math_multinetwork(data_eng_mn::Dict{String,Any}; kron_reduced::Bool=kron_reduced)::Dict{String,Any}
    data_math_mn = Dict{String,Any}(
        "nw" => Dict{String,Any}(),
        "multinetwork" => true
    )
    for (n, nw) in data_eng_mn["nw"]
        for k in _pmd_eng_global_keys
            nw[k] = data_eng_mn[k]
        end

        data_math_mn["nw"][n] = _map_eng2math(nw; kron_reduced=kron_reduced)

        for k in _pmd_math_global_keys
            data_math_mn[k] = data_math_mn["nw"][n][k]
            delete!(data_math_mn["nw"][n], k)
        end
    end

    return data_math_mn
end


"base function for converting engineering model to mathematical model"
function _map_eng2math(data_eng; kron_reduced::Bool=true)
    @assert get(data_eng, "data_model", MATHEMATICAL) == ENGINEERING

    data_math = Dict{String,Any}(
        "name" => get(data_eng, "name", ""),
        "per_unit" => get(data_eng, "per_unit", false),
        "data_model" => MATHEMATICAL,
        "settings" => deepcopy(data_eng["settings"]),
    )

    #TODO the PM tests break for branches which are not of the size indicated by conductors;
    # for now, set to 1 to prevent this from breaking when not kron-reduced
    data_math["conductors"] = kron_reduced ? 3 : 1

    data_math["map"] = Vector{Dict{String,Any}}([
        Dict{String,Any}("unmap_function" => "_map_math2eng_root!")
    ])

    _init_base_components!(data_math)

    # convert buses
    _map_eng2math_bus!(data_math, data_eng; kron_reduced=kron_reduced)

    # convert edges
    _map_eng2math_line!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_switch!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_transformer!(data_math, data_eng; kron_reduced=kron_reduced)

    # convert nodes
    _map_eng2math_load!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_shunt!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_generator!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_solar!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_storage!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_voltage_source!(data_math, data_eng; kron_reduced=kron_reduced)

    # post fix
    if kron_reduced
        #TODO move this out when kron-reducing becomes a transformation
        _kron_reduce_buses!(data_math)
    else
        #TODO fix this in place / throw error instead? IEEE8500 leads to switches
        # with 3x3 R matrices but only 1 phase
        #NOTE: Don't do this when kron-reducing, it will undo the padding
        _slice_branches!(data_math)
    end

    return data_math
end


"converts engineering bus components into mathematical bus components"
function _map_eng2math_bus!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "bus", Dict{String,Any}())
        terminals = eng_obj["terminals"]
        nconductors = data_math["conductors"]

        math_obj = _init_math_obj("bus", name, eng_obj, length(data_math["bus"])+1)

        math_obj["bus_i"] = math_obj["index"]
        math_obj["bus_type"] = _bus_type_conversion(data_eng, eng_obj, "status")

        # take care of grounding; convert to shunt if lossy
        grounded_perfect, shunts = _convert_grounding(eng_obj["terminals"], eng_obj["grounded"], eng_obj["rg"], eng_obj["xg"])

        math_obj["grounded"] = grounded_perfect
        to_sh = []
        for (sh_connections, sh_y) in shunts
            sh_index = length(data_math["shunt"]) + 1
            data_math["shunt"]["$sh_index"] = Dict(
                "index" => sh_index,
                "shunt_bus" => math_obj["bus_i"],
                "connections" => sh_connections,
                "gs" => real.(sh_y),
                "bs" => real.(sh_y),
            )
            push!(to_sh, "shunt.$sh_index")
        end

        if haskey(eng_obj, "vm")
            math_obj["vm"] = eng_obj["vm"]
        end
        if haskey(eng_obj, "va")
            math_obj["va"] = eng_obj["va"]
        end

        math_obj["vmin"] = get(eng_obj, "vm_lb", fill(0.0, length(terminals)))
        math_obj["vmax"] = get(eng_obj, "vm_ub", fill(Inf, length(terminals)))

        if kron_reduced
            filter = terminals.!=kr_neutral
            terminals_kr = terminals[filter]
            @assert all(t in kr_phases for t in terminals_kr) "bus $name has terminals $(terminals), outside of $kr_phases, cannot be kron reduced"

            _apply_filter!(math_obj, ["vm", "va", "vmin", "vmax"], filter)
            _pad_properties!(math_obj, ["vm", "va", "vmin", "vmax"], terminals_kr, kr_phases)
        end

        data_math["bus"]["$(math_obj["index"])"] = math_obj

        if !haskey(data_math, "bus_lookup")
            data_math["bus_lookup"] = Dict{Any,Int}()
        end

        data_math["bus_lookup"][name] = math_obj["index"]

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "bus.$(math_obj["index"])",
            "unmap_function" => "_map_math2eng_bus!",
        ))
    end
end


"converts engineering lines into mathematical branches"
function _map_eng2math_line!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "line", Dict{Any,Dict{String,Any}}())
        _apply_linecode!(eng_obj, data_eng)

        math_obj = _init_math_obj("line", name, eng_obj, length(data_math["branch"])+1)

        nphases = length(eng_obj["f_connections"])
        nconductors = data_math["conductors"]

        math_obj["f_bus"] = data_math["bus_lookup"][eng_obj["f_bus"]]
        math_obj["t_bus"] = data_math["bus_lookup"][eng_obj["t_bus"]]

        math_obj["br_r"] = _impedance_conversion(data_eng, eng_obj, "rs")
        math_obj["br_x"] = _impedance_conversion(data_eng, eng_obj, "xs")

        math_obj["g_fr"] = _admittance_conversion(data_eng, eng_obj, "g_fr")
        math_obj["g_to"] = _admittance_conversion(data_eng, eng_obj, "g_to")

        math_obj["b_fr"] = _admittance_conversion(data_eng, eng_obj, "b_fr")
        math_obj["b_to"] = _admittance_conversion(data_eng, eng_obj, "b_to")

        math_obj["angmin"] = get(eng_obj, "vad_lb", fill(-60.0, nphases))
        math_obj["angmax"] = get(eng_obj, "vad_ub", fill( 60.0, nphases))

        for (f_key, t_key) in [("cm_ub", "c_rating_a"), ("cm_ub_b", "c_rating_b"), ("cm_ub_c", "c_rating_c"),
            ("sm_ub", "rate_a"), ("sm_ub_b", "rate_b"), ("sm_ub_c", "rate_c")]
            if haskey(eng_obj, f_key)
                math_obj[t_key] = eng_obj[f_key]
            end
        end

        math_obj["transformer"] = false
        math_obj["shift"] = zeros(nphases)
        math_obj["tap"] = ones(nphases)

        if kron_reduced
            @assert all(eng_obj["f_connections"].==eng_obj["t_connections"]) "Kron reduction is only supported if f_connections == t_connections"
            filter = _kron_reduce_branch!(math_obj,
                ["br_r", "br_x"], ["g_fr", "b_fr", "g_to", "b_to"],
                eng_obj["f_connections"], kr_neutral
            )
            _apply_filter!(math_obj, ["angmin", "angmax", "tap", "shift"], filter)
            connections = eng_obj["f_connections"][filter]
            _pad_properties!(math_obj, ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to", "shift"], connections, kr_phases)
            _pad_properties!(math_obj, ["angmin"], connections, kr_phases; pad_value=-60.0)
            _pad_properties!(math_obj, ["angmax"], connections, kr_phases; pad_value=60.0)
            _pad_properties!(math_obj, ["tap"], connections, kr_phases; pad_value=1.0)

            for key in ["c_rating_a", "c_rating_b", "c_rating_c", "rate_a", "rate_b", "rate_c"]
                if haskey(math_obj, key)
                    _apply_filter!(math_obj, [key], filter)
                    _pad_properties!(math_obj, [key], connections, kr_phases)
                end
            end
        else
            math_obj["f_connections"] = eng_obj["f_connections"]
            math_obj["t_connections"] = eng_obj["t_connections"]
        end

        math_obj["switch"] = false

        math_obj["br_status"] = Int(eng_obj["status"])

        data_math["branch"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "branch.$(math_obj["index"])",
            "unmap_function" => "_map_math2eng_line!",
        ))
    end
end


"converts engineering n-winding transformers into mathematical ideal 2-winding lossless transformer branches and impedance branches to represent the loss model"
function _map_eng2math_transformer!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "transformer", Dict{Any,Dict{String,Any}}())
        # Build map first, so we can update it as we decompose the transformer
        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => Vector{String}([]),
            "unmap_function" => "_map_math2eng_transformer!",
        ))

        to_map = data_math["map"][end]["to"]

        _apply_xfmrcode!(eng_obj, data_eng)

        if haskey(eng_obj, "f_bus") && haskey(eng_obj, "t_bus")
            @assert all(haskey(eng_obj, k) for k in ["f_bus", "t_bus", "f_connections", "t_connections"]) "Incomplete definition of AL2W tranformer $name, aborting eng2math conversion"

            nphases = length(eng_obj["f_connections"])

            math_obj = Dict{String,Any}(
                "name" => name,
                "source_id" => eng_obj["source_id"],
                "f_bus" => data_math["bus_lookup"][eng_obj["f_bus"]],
                "t_bus" => data_math["bus_lookup"][eng_obj["t_bus"]],
                "f_connections" => eng_obj["f_connections"],
                "t_connections" => eng_obj["t_connections"],
                "configuration" => get(eng_obj, "configuration", WYE),
                "tm_nom" => get(eng_obj, "tm_nom", 1.0),
                "tm_set" => get(eng_obj, "tm_set", fill(1.0, nphases)),
                "tm_fix" => get(eng_obj, "tm_fix", fill(true, nphases)),
                "polarity" => fill(1, nphases),
                "status" => Int(get(eng_obj, "status", ENABLED)),
                "index" => length(data_math["transformer"])+1
            )

            for k in ["tm_lb", "tm_ub"]
                if haskey(eng_obj, k)
                    math_obj[k] = eng_obj[k]
                end
            end

            if kron_reduced
                # TODO fix how padding works, this is a workaround to get bank working
                if eng_obj["configuration"] == WYE
                    f_connections = math_obj["f_connections"]
                    _pad_properties!(math_obj, ["tm_lb", "tm_ub", "tm_set"], f_connections[f_connections.!=kr_neutral], kr_phases; pad_value=1.0)
                    _pad_properties!(math_obj, ["tm_fix"], f_connections[f_connections.!=kr_neutral], kr_phases; pad_value=false)

                    math_obj["f_connections"] = math_obj["t_connections"]
                end
            end

            data_math["transformer"]["$(math_obj["index"])"] = math_obj

            push!(to_map, "transformer.$(math_obj["index"])")
        else
            vnom = eng_obj["vm_nom"] * data_eng["settings"]["voltage_scale_factor"]
            snom = eng_obj["sm_nom"] * data_eng["settings"]["power_scale_factor"]

            nrw = length(eng_obj["bus"])

            # calculate zbase in which the data is specified, and convert to SI
            zbase = (vnom.^2) ./ snom

            # x_sc is specified with respect to first winding
            x_sc = eng_obj["xsc"] .* zbase[1]

            # rs is specified with respect to each winding
            r_s = eng_obj["rw"] .* zbase

            g_sh =  (eng_obj["noloadloss"]*snom[1])/vnom[1]^2
            b_sh = -(eng_obj["cmag"]*snom[1])/vnom[1]^2

            # data is measured externally, but we now refer it to the internal side
            ratios = vnom/data_eng["settings"]["voltage_scale_factor"]
            x_sc = x_sc./ratios[1]^2
            r_s = r_s./ratios.^2
            g_sh = g_sh*ratios[1]^2
            b_sh = b_sh*ratios[1]^2

            # convert x_sc from list of upper triangle elements to an explicit dict
            y_sh = g_sh + im*b_sh
            z_sc = Dict([(key, im*x_sc[i]) for (i,key) in enumerate([(i,j) for i in 1:nrw for j in i+1:nrw])])

            #TODO remove once moving out kron-reduction
            dims = kron_reduced ? 3 : length(eng_obj["tm_set"][1])
            transformer_t_bus_w = _build_loss_model!(data_math, name, to_map, r_s, z_sc, y_sh, nphases=dims)

            for w in 1:nrw
                # 2-WINDING TRANSFORMER
                # make virtual bus and mark it for reduction
                tm_nom = eng_obj["configuration"][w]==DELTA ? eng_obj["vm_nom"][w]*sqrt(3) : eng_obj["vm_nom"][w]
                transformer_2wa_obj = Dict{String,Any}(
                    "name"          => "_virtual_transformer.$name.$w",
                    "source_id"     => "_virtual_transformer.$(eng_obj["source_id"]).$w",
                    "f_bus"         => data_math["bus_lookup"][eng_obj["bus"][w]],
                    "t_bus"         => transformer_t_bus_w[w],
                    "tm_nom"        => tm_nom,
                    "f_connections" => eng_obj["connections"][w],
                    "t_connections" => collect(1:dims+1),
                    "configuration" => eng_obj["configuration"][w],
                    "polarity"      => eng_obj["polarity"][w],
                    "tm_set"        => eng_obj["tm_set"][w],
                    "tm_fix"        => eng_obj["tm_fix"][w],
                    "status"        => Int(get(eng_obj, "status", ENABLED)),
                    "index"         => length(data_math["transformer"])+1
                )

                for prop in ["tm_lb", "tm_ub", "tm_step"]
                    if haskey(eng_obj, prop)
                        transformer_2wa_obj[prop] = eng_obj[prop][w]
                    end
                end

                if kron_reduced
                    # TODO fix how padding works, this is a workaround to get bank working
                    if all(conf==WYE for conf in eng_obj["configuration"])
                        f_connections = transformer_2wa_obj["f_connections"]
                        _pad_properties!(transformer_2wa_obj, ["tm_lb", "tm_ub", "tm_set"], f_connections[f_connections.!=kr_neutral], kr_phases; pad_value=1.0)
                        _pad_properties!(transformer_2wa_obj, ["tm_fix"], f_connections[f_connections.!=kr_neutral], kr_phases; pad_value=false)

                        transformer_2wa_obj["f_connections"] = transformer_2wa_obj["t_connections"]
                    end
                end

                data_math["transformer"]["$(transformer_2wa_obj["index"])"] = transformer_2wa_obj

                push!(to_map, "transformer.$(transformer_2wa_obj["index"])")
            end
        end
    end
end


"converts engineering switches into mathematical switches and (if neeed) impedance branches to represent loss model"
function _map_eng2math_switch!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    # TODO enable real switches (right now only using vitual lines)
    for (name, eng_obj) in get(data_eng, "switch", Dict{Any,Dict{String,Any}}())
        nphases = length(eng_obj["f_connections"])
        nconductors = data_math["conductors"]

        math_obj = _init_math_obj("switch", name, eng_obj, length(data_math["switch"])+1)

        math_obj["f_bus"] = data_math["bus_lookup"][eng_obj["f_bus"]]
        math_obj["t_bus"] = data_math["bus_lookup"][eng_obj["t_bus"]]

        math_obj["state"] = get(eng_obj, "state", CLOSED)

        # OPF bounds
        for (fr_key, to_key) in [("cm_ub", "c_rating")]
            if haskey(eng_obj, fr_key)
                math_obj[to_key] = eng_obj[fr_key]
            end
        end

        map_to = "switch.$(math_obj["index"])"

        if haskey(eng_obj, "linecode")
            _apply_linecode!(eng_obj, data_eng)
        end

        if true
        # if !all(isapprox.(get(eng_obj, "rs", zeros(1, 1)), 0)) && !all(isapprox.(get(eng_obj, "xs", zeros(1, 1)), 0)) # TODO enable real switches
            # build virtual bus

            f_bus = deepcopy(data_math["bus"]["$(math_obj["f_bus"])"])
            t_bus = deepcopy(data_math["bus"]["$(math_obj["t_bus"])"])

            bus_obj = Dict{String,Any}(
                "name" => "_virtual_bus.switch.$name",
                "bus_i" => length(data_math["bus"])+1,
                "bus_type" => get(eng_obj, "state", CLOSED) == OPEN ? 4 : 1,
                "terminals" => t_bus["terminals"],  # connected to the switch on the to-side
                "grounded" => t_bus["grounded"],  # connected to the switch on the to-side
                "vmin" => f_bus["vmin"],
                "vmax" => f_bus["vmax"],
                "index" => length(data_math["bus"])+1,
            )

            #= TODO enable real switches
            # math_obj["t_bus"] = bus_obj["bus_i"]
            # data_math["bus"]["$(bus_obj["index"])"] = bus_obj
            =#

            # TODO remove after enabling real switches
            if all(isapprox.(get(eng_obj, "rs", zeros(nphases, nphases)), 0))
                eng_obj["rs"] = fill(1e-4, nphases, nphases)
            end

            if all(isapprox.(get(eng_obj, "xs", zeros(nphases, nphases)), 0))
                eng_obj["xs"] = fill(1e-3, nphases, nphases)
            end

            branch_obj = _init_math_obj("line", name, eng_obj, length(data_math["branch"])+1)

            _branch_obj = Dict{String,Any}(
                "name" => "_virtual_branch.switch.$name",
                "source_id" => "_virtual_branch.switch.$name",
                # "f_bus" => bus_obj["bus_i"],  # TODO enable real switches
                "f_bus" => data_math["bus_lookup"][eng_obj["f_bus"]],
                "t_bus" => data_math["bus_lookup"][eng_obj["t_bus"]],
                "f_connections" => eng_obj["t_connections"],  # the virtual branch connects to the switch on the to-side
                "t_connections" => eng_obj["t_connections"],  # should be identical to the switch's to-side connections
                "br_r" => _impedance_conversion(data_eng, eng_obj, "rs"),
                "br_x" => _impedance_conversion(data_eng, eng_obj, "xs"),
                "g_fr" => zeros(nphases, nphases),
                "g_to" => zeros(nphases, nphases),
                "b_fr" => zeros(nphases, nphases),
                "b_to" => zeros(nphases, nphases),
                "angmin" => fill(-60.0, nphases),
                "angmax" => fill( 60.0, nphases),
                "transformer" => false,
                "shift" => zeros(nphases),
                "tap" => ones(nphases),
                "switch" => false,
                "br_status" => get(eng_obj, "state", CLOSED) == OPEN ? 0 : 1,
            )

            merge!(branch_obj, _branch_obj)

            if kron_reduced
                @assert all(eng_obj["f_connections"].==eng_obj["t_connections"]) "Kron reduction is only supported if f_connections == t_connections"
                filter = _kron_reduce_branch!(branch_obj,
                    ["br_r", "br_x"], ["g_fr", "b_fr", "g_to", "b_to"],
                    eng_obj["f_connections"], kr_neutral
                )
                _apply_filter!(branch_obj, ["angmin", "angmax", "tap", "shift"], filter)
                connections = eng_obj["f_connections"][filter]
                _pad_properties!(branch_obj, ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to", "shift"], connections, kr_phases)
                _pad_properties!(branch_obj, ["angmin"], connections, kr_phases; pad_value=-60.0)
                _pad_properties!(branch_obj, ["angmax"], connections, kr_phases; pad_value=60.0)
                _pad_properties!(branch_obj, ["tap"], connections, kr_phases; pad_value=1.0)
            else
                branch_obj["f_connections"] = eng_obj["f_connections"]
                branch_obj["f_connections"] = eng_obj["t_connections"]
            end

            data_math["branch"]["$(branch_obj["index"])"] = branch_obj

            map_to = ["branch.$(branch_obj["index"])"]
            # map_to = [map_to, "bus.$(bus_obj["index"])", "branch.$(branch_obj["index"])"]  # TODO enable real switches
        end

        # data_math["switch"]["$(math_obj["index"])"] = math_obj  # TODO enable real switches

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => map_to,
            "unmap_function" => "_map_math2eng_switch!",
        ))
    end
end


"converts engineering generic shunt components into mathematical shunt components"
function _map_eng2math_shunt!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "shunt", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("shunt", name, eng_obj, length(data_math["shunt"])+1)

        # TODO change to new capacitor shunt calc logic
        math_obj["shunt_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        math_obj["gs"] = get(eng_obj, "gs", zeros(size(eng_obj["bs"])))

        if kron_reduced
            filter = _kron_reduce_branch!(math_obj,
                Vector{String}([]), ["gs", "bs"],
                eng_obj["connections"], kr_neutral
            )
            connections = eng_obj["connections"][filter]
            _pad_properties!(math_obj, ["gs", "bs"], connections, kr_phases)
        else
            math_obj["connections"] = eng_obj["connections"]
        end

        data_math["shunt"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "shunt.$(math_obj["index"])",
            "unmap_function" => "_map_math2eng_shunt!",
        ))
    end
end


"converts engineering load components into mathematical load components"
function _map_eng2math_load!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "load", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("load", name, eng_obj, length(data_math["load"])+1)

        connections = eng_obj["connections"]

        math_obj["load_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        math_obj["pd"] = eng_obj["pd_nom"]
        math_obj["qd"] = eng_obj["qd_nom"]

        if kron_reduced
            if math_obj["configuration"]==WYE
                @assert connections[end]==kr_neutral "for wye-connected loads, if kron_reduced the connections list should end with a neutral"
                _pad_properties!(math_obj, ["pd", "qd"], connections[connections.!=kr_neutral], kr_phases)
            else
                _pad_properties_delta!(math_obj, ["pd", "qd"], connections, kr_phases)
            end
        else
            math_obj["connections"] = connections
        end

        math_obj["vnom_kv"] = eng_obj["vm_nom"]

        data_math["load"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "load.$(math_obj["index"])",
            "unmap_function" => "_map_math2eng_load!",
        ))
   end
end


"converts engineering generators into mathematical generators"
function _map_eng2math_generator!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "generator", Dict{String,Any}())
        math_obj = _init_math_obj("generator", name, eng_obj, length(data_math["gen"])+1)

        connections = eng_obj["connections"]
        nconductors = data_math["conductors"]

        math_obj["gen_bus"] = data_math["bus_lookup"][eng_obj["bus"]]
        math_obj["gen_status"] = Int(eng_obj["status"])
        math_obj["control_mode"] = get(eng_obj, "control_mode", FREQUENCYDROOP)
        math_obj["pmax"] = get(eng_obj, "pg_ub", fill(Inf, nconductors))

        for (f_key, t_key) in [("qg_lb", "qmin"), ("qg_ub", "qmax"), ("pg_lb", "pmin")]
            if haskey(eng_obj, f_key)
                math_obj[t_key] = eng_obj[f_key]
            end
        end

        _add_gen_cost_model!(math_obj, eng_obj)

        math_obj["configuration"] = get(eng_obj, "configuration", WYE)

        if kron_reduced
            if math_obj["configuration"]==WYE
                @assert connections[end]==kr_neutral "For WYE connected generators, if kron_reduced the connections list should end with a neutral conductor"
                _pad_properties!(math_obj, ["pg", "qg", "vg", "pmin", "pmax", "qmin", "qmax"], connections[1:end-1], kr_phases)
            else
                _pad_properties_delta!(math_obj, ["pg", "qg", "vg", "pmin", "pmax", "qmin", "qmax"], connections, kr_phases)
            end
        else
            math_obj["connections"] = connections
        end

        # if PV generator mode convert attached bus to PV bus
        if math_obj["control_mode"] == ISOCHRONOUS
            data_math["bus"]["$(data_math["bus_lookup"][eng_obj["bus"]])"]["bus_type"] = 2
        end

        data_math["gen"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "gen.$(math_obj["index"])",
            "unmap_function" => "_map_math2eng_generator!",
        ))
    end
end


"converts engineering solar components into mathematical generators"
function _map_eng2math_solar!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "solar", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("solar", name, eng_obj, length(data_math["gen"])+1)

        connections = eng_obj["connections"]
        nconductors = data_math["conductors"]

        math_obj["gen_bus"] = data_math["bus_lookup"][eng_obj["bus"]]
        math_obj["gen_status"] = Int(eng_obj["status"])

        for (fr_k, to_k) in [("vg", "vg"), ("pg_lb", "pmin"), ("pg_ub", "pmax"), ("qg_lb", "qmin"), ("qg_ub", "qmax")]
            if haskey(eng_obj, fr_k)
                math_obj[to_k] = eng_obj[fr_k]
            end
        end

        _add_gen_cost_model!(math_obj, eng_obj)

        if kron_reduced
            if math_obj["configuration"]==WYE
                @assert(connections[end]==kr_neutral)
                _pad_properties!(math_obj, ["pg", "qg", "vg", "pmin", "pmax", "qmin", "qmax"], connections[1:end-1], kr_phases)
            else
                _pad_properties_delta!(math_obj, ["pg", "qg", "vg", "pmin", "pmax", "qmin", "qmax"], connections, kr_phases)
            end
        else
            math_obj["connections"] = connections
        end

        data_math["gen"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "gen.$(math_obj["index"])",
            "unmap_function" => "_map_math2eng_solar!",
        ))
    end
end


"converts engineering storage into mathematical storage"
function _map_eng2math_storage!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "storage", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("storage", name, eng_obj, length(data_math["storage"])+1)

        connections = eng_obj["connections"]
        nconductors = data_math["conductors"]

        math_obj["storage_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        # needs to be in units MW
        math_obj["energy"] = eng_obj["energy"] * data_eng["settings"]["power_scale_factor"] / 1e6
        #TODO is scale factor correct?
        math_obj["energy_rating"] = eng_obj["energy_ub"] * data_eng["settings"]["power_scale_factor"] / 1e6
        math_obj["charge_rating"] = eng_obj["charge_ub"] * data_eng["settings"]["power_scale_factor"] / 1e6
        math_obj["discharge_rating"] = eng_obj["discharge_ub"] * data_eng["settings"]["power_scale_factor"] / 1e6
        math_obj["charge_efficiency"] = eng_obj["charge_efficiency"] / 100.0
        math_obj["discharge_efficiency"] = eng_obj["discharge_efficiency"] / 100.0
        #TODO is scale factor correct? what  should be the  unit?
        math_obj["thermal_rating"] = eng_obj["cm_ub"] .* data_eng["settings"]["power_scale_factor"] ./ 1e6
        math_obj["qmin"] = eng_obj["qs_lb"] .* data_eng["settings"]["power_scale_factor"] ./ 1e6
        math_obj["qmax"] =  eng_obj["qs_ub"] .* data_eng["settings"]["power_scale_factor"] ./ 1e6
        math_obj["r"] = eng_obj["rs"]
        math_obj["x"] = eng_obj["xs"]
        math_obj["p_loss"] = eng_obj["pex"] .* data_eng["settings"]["power_scale_factor"] ./ 1e6
        math_obj["q_loss"] = eng_obj["qex"] .* data_eng["settings"]["power_scale_factor"] ./ 1e6

        math_obj["ps"] = get(eng_obj, "ps", zeros(size(eng_obj["cm_ub"])))
        math_obj["qs"] = get(eng_obj, "qs", zeros(size(eng_obj["cm_ub"])))

        if kron_reduced
            _pad_properties!(math_obj, ["thermal_rating", "qmin", "qmax", "r", "x", "ps", "qs"], connections, kr_phases)
        end

        data_math["storage"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "storage.$(math_obj["index"])",
            "unmap_function" => "_map_math2eng_storage!",
        ))
    end
end


"converts engineering voltage sources into mathematical generators and (if needed) impedance branches to represent the loss model"
function _map_eng2math_voltage_source!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1,2,3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "voltage_source", Dict{Any,Any}())
        nconductors = length(eng_obj["vm"])

        math_obj = _init_math_obj("voltage_source", name, eng_obj, length(data_math["gen"])+1)

        math_obj["name"] = "_virtual_gen.voltage_source.$name"
        math_obj["gen_bus"] = gen_bus = data_math["bus_lookup"][eng_obj["bus"]]
        math_obj["gen_status"] = Int(eng_obj["status"])
        math_obj["pg"] = fill(0.0, nconductors)
        math_obj["qg"] = fill(0.0, nconductors)
        math_obj["pmin"] = get(eng_obj, "pg_lb", fill(-Inf, nconductors))
        math_obj["pmax"] = get(eng_obj, "pg_ub", fill( Inf, nconductors))
        math_obj["qmin"] = get(eng_obj, "qg_lb", fill(-Inf, nconductors))
        math_obj["qmax"] = get(eng_obj, "qg_ub", fill( Inf, nconductors))
        math_obj["configuration"] = WYE
        math_obj["source_id"] = "_virtual_gen.$(eng_obj["source_id"])"

        _add_gen_cost_model!(math_obj, eng_obj)

        map_to = "gen.$(math_obj["index"])"

        if !all(isapprox.(get(eng_obj, "rs", zeros(1, 1)), 0)) && !all(isapprox.(get(eng_obj, "xs", zeros(1, 1)), 0))
            bus_obj = Dict{String,Any}(
                "bus_i" => length(data_math["bus"])+1,
                "index" => length(data_math["bus"])+1,
                "terminals" => collect(1:nconductors+1),
                "grounded" => [fill(false, nconductors)..., true],
                "name" => "_virtual_bus.voltage_source.$name",
                "bus_type" => 3,
                "vm" => [eng_obj["vm"]..., 0.0],
                "va" => [eng_obj["va"]..., 0.0],
                "vmin" => [eng_obj["vm"]..., 0.0],
                "vmax" => [eng_obj["vm"]..., 0.0]
            )

            math_obj["gen_bus"] = gen_bus = bus_obj["bus_i"]

            data_math["bus"]["$(bus_obj["index"])"] = bus_obj

            branch_obj = Dict{String,Any}(
                "name" => "_virtual_branch.voltage_source.$name",
                "source_id" => "_virtual_branch.$(eng_obj["source_id"])",
                "f_bus" => bus_obj["bus_i"],
                "t_bus" => data_math["bus_lookup"][eng_obj["bus"]],
                "f_connections" => eng_obj["connections"],
                "t_connections" => eng_obj["connections"],
                "angmin" => fill(-60.0, nconductors),
                "angmax" => fill( 60.0, nconductors),
                "shift" => fill(0.0, nconductors),
                "tap" => fill(1.0, nconductors),
                "tranformer" => false,
                "switch" => false,
                "br_status" => 1,
                "br_r" => _impedance_conversion(data_eng, eng_obj, "rs"),
                "br_x" => _impedance_conversion(data_eng, eng_obj, "xs"),
                "g_fr" => zeros(nconductors, nconductors),
                "g_to" => zeros(nconductors, nconductors),
                "b_fr" => zeros(nconductors, nconductors),
                "b_to" => zeros(nconductors, nconductors),
                "index" => length(data_math["branch"])+1
            )

            # finally, we have to set a neutral for the virtual generator
            neutral = _get_ground_math!(data_math["bus"]["$gen_bus"], exclude_terminals=[1:nconductors...])
            math_obj["connections"] = [collect(1:nconductors)..., neutral]

            data_math["branch"]["$(branch_obj["index"])"] = branch_obj

            map_to = [map_to, "bus.$(bus_obj["index"])", "branch.$(branch_obj["index"])"]
        else
            data_math["bus"]["$gen_bus"]["vmin"] = [eng_obj["vm"]..., 0.0]
            data_math["bus"]["$gen_bus"]["vmax"] = [eng_obj["vm"]..., 0.0]
            data_math["bus"]["$gen_bus"]["vm"] = [eng_obj["vm"]..., 0.0]
            data_math["bus"]["$gen_bus"]["va"] = [eng_obj["va"]..., 0.0]
            data_math["bus"]["$gen_bus"]["bus_type"] = 3
        end

        data_math["gen"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => map_to,
            "unmap_function" => "_map_math2eng_voltage_source!",
        ))
    end
end
