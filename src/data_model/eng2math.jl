import LinearAlgebra: diagm

"items that are mapped one-to-one from engineering to math models"
const _1to1_maps = Dict{String,Vector{String}}(
    "bus" => ["vm", "va", "terminals", "phases", "neutral", "dss"],
    "line" => ["f_connections", "t_connections", "source_id", "dss"],
    "transformer" => ["f_connections", "t_connections", "source_id", "dss"],
    "switch" => ["status", "state", "f_connections", "t_connections", "source_id", "dss"],
    "line_reactor" => ["f_connections", "t_connections", "source_id", "dss"],
    "series_capacitor" => ["f_connections", "t_connections", "source_id", "dss"],
    "shunt" => ["status", "gs", "bs", "connections", "source_id", "dss"],
    "shunt_capacitor" => ["status", "bs", "connections", "source_id", "dss"],
    "shunt_reactor" => ["status", "connections", "source_id", "dss"],
    "load" => ["model", "configuration", "connections", "status", "source_id", "dss"],
    "generator" => ["pg", "qg", "configuration", "connections", "source_id", "dss"],
    "solar" => ["configuration", "connections", "source_id", "dss"],
    "storage" => ["status", "energy", "ps", "qs", "connections", "source_id", "dss"],
    "voltage_source" => ["source_id", "dss"],
)

"list of nodal type elements in the engineering model"
const _node_elements = Vector{String}([
    "load", "capacitor", "shunt_reactor", "generator", "solar", "storage", "vsource"
])

"list of edge type elements in the engineering model"
const _edge_elements = Vector{String}([
    "line", "switch", "transformer", "line_reactor", "series_capacitor"
])

"list of time-series supported parameters that map one-to-one"
const _time_series_parameters = Dict{String,Dict{String,Tuple{Function, String}}}(
    "switch" => Dict{String,Tuple{Function, String}}(
        "status" => (_no_conversion, "status"),
        "state" => (_no_conversion, "state"),
        "c_rating" => (_no_conversion, "cm_ub"),
        "s_rating" => (_no_conversion, "sm_ub"),
        "br_r" => (_impedance_conversion, "rs"),
        "br_x" => (_impedance_conversion, "xs"),
        "g_fr" => (_admittance_conversion, "g_fr"),
        "g_to" => (_admittance_conversion, "g_to"),
        "b_fr" => (_admittance_conversion, "b_fr"),
        "b_to" => (_admittance_conversion, "b_to")
    ),
    "fuse" => Dict{String,Tuple{Function, String}}(
        "status" => (_no_conversion, "status"),
        "state" => (_no_conversion, "state"),
        "c_rating" => (_no_conversion, "cm_ub"),
        "s_rating" => (_no_conversion, "sm_ub"),
        "br_r" => (_impedance_conversion, "rs"),
        "br_x" => (_impedance_conversion, "xs"),
        "g_fr" => (_admittance_conversion, "g_fr"),
        "g_to" => (_admittance_conversion, "g_to"),
        "b_fr" => (_admittance_conversion, "b_fr"),
        "b_to" => (_admittance_conversion, "b_to")
    ),
    "line" => Dict{String,Tuple{Function, String}}(
        "br_status" => (_no_conversion, "status"),
        "c_rating" => (_no_conversion, "cm_ub"),
        "s_rating" => (_no_conversion, "sm_ub"),
        "br_r" => (_impedance_conversion, "rs"),
        "br_x" => (_impedance_conversion, "xs"),
        "g_fr" => (_admittance_conversion, "g_fr"),
        "g_to" => (_admittance_conversion, "g_to"),
        "b_fr" => (_admittance_conversion, "b_fr"),
        "b_to" => (_admittance_conversion, "b_to")
    ),
    "transformer" => Dict{String,Tuple{Function, String}}(
        "status" => (_no_conversion, "status"),
        # TODO need to figure out how to convert values for time series for decomposed transformers
    ),
    "bus" => Dict{String,Tuple{Function, String}}(
        "bus_type" => (_bus_type_conversion, "status"),
        "vmin" => (_no_conversion, "vm_lb"),
        "vmax" => (_no_conversion, "vm_ub"),
        "vm" => (_no_conversion, "vm"),
        "va" => (_no_conversion, "va")
    ),
    "shunt" => Dict{String,Tuple{Function, String}}(
        "status" => (_no_conversion, "status"),
        "gs" => (_no_conversion, "gs"),
        "bs" => (_no_conversion, "bs"),
    ),
    "shunt_capacitor" => Dict{String,Tuple{Function, String}}(
        "status" => (_no_conversion, "status"),
        "gs" => (_no_conversion, "gs"),
        "bs" => (_no_conversion, "bs"),
    ),
    "shunt_reactor" => Dict{String,Tuple{Function, String}}(
        "status" => (_no_conversion, "status"),
        "gs" => (_no_conversion, "gs"),
        "bs" => (_no_conversion, "bs"),
    ),
    "load" => Dict{String,Tuple{Function, String}}(
        "status" => (_no_conversion, "status"),
        "pd_ref" => (_no_conversion, "pd_nom"),
        "qd_ref" => (_no_conversion, "qd_nom"),
    ),
    "generator" => Dict{String,Tuple{Function, String}}(
        "gen_status" => (_no_conversion, "status"),
        "pg" => (_no_conversion, "pg"),
        "qg" => (_no_conversion, "qg"),
        "vg" => (_vnom_conversion, "vg"),
        "pmin" => (_no_conversion, "pg_lb"),
        "pmax" => (_no_conversion, "pg_ub"),
        "qmin" => (_no_conversion, "qg_lb"),
        "qmax" => (_no_conversion, "qg_ub"),
    ),
    "solar" => Dict{String,Tuple{Function, String}}(
        "gen_status" => (_no_conversion, "status"),
        "pg" => (_no_conversion, "pg"),
        "qg" => (_no_conversion, "qg"),
        "vg" => (_vnom_conversion, "vg"),
        "pmin" => (_no_conversion, "pg_lb"),
        "pmax" => (_no_conversion, "pg_ub"),
        "qmin" => (_no_conversion, "qg_lb"),
        "qmax" => (_no_conversion, "qg_ub"),
    ),
    "storage" => Dict{String,Tuple{Function, String}}(
        "status" => (_no_conversion, "status"),
        "energy" => (_no_conversion, "energy"),
        "energy_rating" => (_no_conversion, "energy_ub"),
        "charge_rating" => (_no_conversion, "charge_ub"),
        "discharge_rating" => (_no_conversion, "discharge_ub"),
        "charge_efficiency" => (_no_conversion, "charge_efficiency"),
        "discharge_efficiency" => (_no_conversion, "discharge_efficiency"),
        "thermal_rating" => (_no_conversion, "cm_ub"),
        "qmin" => (_no_conversion, "qs_lb"),
        "qmax" => (_no_conversion, "qs_ub"),
        "r" => (_no_conversion, "rs"),
        "x" => (_no_conversion, "xs"),
        "p_loss" => (_no_conversion, "pex"),
        "q_loss" => (_no_conversion, "qex"),
        "ps" => (_no_conversion, "ps"),
        "qs" => (_no_conversion, "qs"),
    ),
    "voltage_source" => Dict{String,Tuple{Function, String}}(
        "gen_status" => (_no_conversion, "status"),
        "vm" => (_no_conversion, "vm"),
        "va" => (_angle_shift_conversion, "va"),
    ),

)

"base function for converting engineering model to mathematical model"
function _map_eng2math(data_eng; kron_reduced::Bool=true)
    @assert get(data_eng, "data_model", "mathematical") == "engineering"

    data_math = Dict{String,Any}(
        "name" => get(data_eng, "name", ""),
        "per_unit" => get(data_eng, "per_unit", false),
        "data_model" => "mathematical",
        "settings" => data_eng["settings"],
    )

    #TODO the PM tests break for branches which are not of the size indicated by conductors;
    # for now, set to 1 to prevent this from breaking when not kron-reduced
    data_math["conductors"] = kron_reduced ? 3 : 1
    data_math["basekv"] = data_eng["settings"]["vbase"]
    data_math["baseMVA"] = data_eng["settings"]["sbase"]*data_eng["settings"]["v_var_scalar"]/1E6

    data_math["map"] = Dict{Int,Dict{Symbol,Any}}(
        1 => Dict{Symbol,Any}(
            :unmap_function => :_map_math2eng_root!,
        )
    )

    _init_base_components!(data_math)

    # convert buses
    _map_eng2math_bus!(data_math, data_eng; kron_reduced=kron_reduced)

    # convert edges
    _map_eng2math_line!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_switch!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_transformer!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_line_reactor!(data_math, data_eng; kron_reduced=kron_reduced)
    # _map_eng2math_series_capacitor(data_math, data_eng; kron_reduced=kron_reduced)  # TODO build conversion for series capacitors

    # convert nodes
    _map_eng2math_load!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_shunt!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_shunt_capacitor!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_shunt_reactor!(data_math, data_eng; kron_reduced=kron_reduced)

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

        math_obj["base_kv"] = data_eng["settings"]["vbase"]

        if kron_reduced
            filter = terminals.!=kr_neutral
            terminals_kr = terminals[filter]
            @assert(all(t in kr_phases for t in terminals_kr))

            _apply_filter!(math_obj, ["vm", "va", "vmin", "vmax"], filter)
            _pad_properties!(math_obj, ["vm", "va", "vmin", "vmax"], terminals_kr, kr_phases)
        end

        data_math["bus"]["$(math_obj["index"])"] = math_obj

        if !haskey(data_math, "bus_lookup")
            data_math["bus_lookup"] = Dict{Any,Int}()
        end

        data_math["bus_lookup"][name] = math_obj["index"]

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => name,
            :to => "bus.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_bus!,
        )

        # time series
        for (fr, (f, to)) in _time_series_parameters["bus"]
            if haskey(eng_obj, "$(fr)_time_series")
                time_series = data_eng["time_series"][eng_obj["$(fr)_time_series"]]
                _parse_time_series_parameter!(data_math, time_series, eng_obj, "bus", "$(math_obj["index"])", fr, to, f)
            end
        end
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

        math_obj["angmin"] = fill(-60.0, nphases)
        math_obj["angmax"] = fill( 60.0, nphases)

        math_obj["transformer"] = false
        math_obj["shift"] = zeros(nphases)
        math_obj["tap"] = ones(nphases)

        f_bus = data_eng["bus"][eng_obj["f_bus"]]
        t_bus = data_eng["bus"][eng_obj["t_bus"]]

        if kron_reduced
            @assert(all(eng_obj["f_connections"].==eng_obj["t_connections"]), "Kron reduction is only supported if f_connections is the same as t_connections.")
            filter = _kron_reduce_branch!(math_obj,
                ["br_r", "br_x"], ["g_fr", "b_fr", "g_to", "b_to"],
                eng_obj["f_connections"], kr_neutral
            )
            _apply_filter!(math_obj, ["angmin", "angmax", "tap", "shift"], filter)
            connections = eng_obj["f_connections"][filter]
            _pad_properties!(math_obj, ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to", "angmin", "angmax", "tap", "shift"], connections, kr_phases)
        else
            math_obj["f_connections"] = eng_obj["f_connections"]
            math_obj["t_connections"] = eng_obj["t_connections"]
        end

        math_obj["switch"] = false

        math_obj["br_status"] = eng_obj["status"]

        data_math["branch"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => name,
            :to => "branch.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_line!,
        )

        # time series
        # TODO
        for (fr, to) in zip(["status"], ["status"])
            if haskey(eng_obj, "$(fr)_time_series")
                time_series = data_eng["time_series"][eng_obj["$(fr)_time_series"]]
                _parse_time_series_parameter!(data_math, time_series, eng_obj[fr], "branch", "$(math_obj["index"])", to)
            end
        end
    end
end


"converts engineering n-winding transformers into mathematical ideal 2-winding lossless transformer branches and impedance branches to represent the loss model"
function _map_eng2math_transformer!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "transformer", Dict{Any,Dict{String,Any}}())
        # Build map first, so we can update it as we decompose the transformer
        map_idx = length(data_math["map"])+1
        data_math["map"][map_idx] = Dict{Symbol,Any}(
            :from => name,
            :to => Vector{String}([]),
            :unmap_function => :_map_math2eng_transformer!,
        )

        to_map = data_math["map"][map_idx][:to]

        _apply_xfmrcode!(eng_obj, data_eng)

        vnom = eng_obj["vnom"] * data_eng["settings"]["v_var_scalar"]
        snom = eng_obj["snom"] * data_eng["settings"]["v_var_scalar"]

        nrw = length(eng_obj["bus"])

        # calculate zbase in which the data is specified, and convert to SI
        zbase = (vnom.^2) ./ snom

        # x_sc is specified with respect to first winding
        x_sc = eng_obj["xsc"] .* zbase[1]

        # rs is specified with respect to each winding
        r_s = eng_obj["rs"] .* zbase

        g_sh =  (eng_obj["noloadloss"]*snom[1])/vnom[1]^2
        b_sh = -(eng_obj["imag"]*snom[1])/vnom[1]^2

        # data is measured externally, but we now refer it to the internal side
        ratios = vnom/data_eng["settings"]["v_var_scalar"]
        x_sc = x_sc./ratios[1]^2
        r_s = r_s./ratios.^2
        g_sh = g_sh*ratios[1]^2
        b_sh = b_sh*ratios[1]^2

        # convert x_sc from list of upper triangle elements to an explicit dict
        y_sh = g_sh + im*b_sh
        z_sc = Dict([(key, im*x_sc[i]) for (i,key) in enumerate([(i,j) for i in 1:nrw for j in i+1:nrw])])

        dims = length(eng_obj["tm"][1])
        transformer_t_bus_w = _build_loss_model!(data_math, name, to_map, r_s, z_sc, y_sh, nphases=dims)

        for w in 1:nrw
            # 2-WINDING TRANSFORMER
            # make virtual bus and mark it for reduction
            tm_nom = eng_obj["configuration"][w]=="delta" ? eng_obj["vnom"][w]*sqrt(3) : eng_obj["vnom"][w]
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
                "tm"            => eng_obj["tm"][w],
                "fixed"         => eng_obj["fixed"][w],
                "index"         => length(data_math["transformer"])+1
            )

            for prop in ["tm_min", "tm_max", "tm_step"]
                if haskey(eng_obj, prop)
                    transformer_2wa_obj[prop] = eng_obj[prop][w]
                end
            end

            if kron_reduced
                # TODO fix how padding works, this is a workaround to get bank working
                if all(eng_obj["configuration"] .== "wye")
                    f_connections = transformer_2wa_obj["f_connections"]
                    _pad_properties!(transformer_2wa_obj, ["tm_min", "tm_max", "tm", "fixed"], f_connections[f_connections.!=kr_neutral], kr_phases; pad_value=1.0)

                    transformer_2wa_obj["f_connections"] = transformer_2wa_obj["t_connections"]
                end
            end

            data_math["transformer"]["$(transformer_2wa_obj["index"])"] = transformer_2wa_obj

            push!(to_map, "transformer.$(transformer_2wa_obj["index"])")
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
        math_obj["t_bus"] = data_math["bus_lookup"][eng_obj["f_bus"]]

        # OPF bounds
        for (fr_key, to_key) in zip(["cm_ub"], ["c_rating"])
            if haskey(eng_obj, fr_key)
                math_obj[to_key] = eng_obj[fr_key]
            end
        end

        map_to = "switch.$(math_obj["index"])"

        # time series
        # TODO switch time series
        for (fr, to) in zip(["status", "state"], ["status", "state"])
            if haskey(eng_obj, "$(fr)_time_series")
                time_series = data_eng["time_series"][eng_obj["$(fr)_time_series"]]
                _parse_time_series_parameter!(data_math, time_series, eng_obj[fr], "switch", "$(math_obj["index"])", to)
            end
        end

        if any(haskey(eng_obj, k) for k in ["rs", "xs", "linecode"])
            # build virtual bus

            f_bus = data_math["bus"]["$(math_obj["f_bus"])"]

            bus_obj = Dict{String,Any}(
                "name" => "_virtual_bus.switch.$name",
                "bus_i" => length(data_math["bus"])+1,
                "bus_type" => 1,
                "vmin" => f_bus["vmin"],
                "vmax" => f_bus["vmax"],
                "base_kv" => f_bus["base_kv"],
                "status" => 1,
                "index" => length(data_math["bus"])+1,
            )

            #= TODO enable real switches
            # math_obj["t_bus"] = bus_obj["bus_i"]
            # data_math["bus"]["$(bus_obj["index"])"] = bus_obj
            =#

            # build virtual branch
            _apply_linecode!(eng_obj, data_eng)

            branch_obj = _init_math_obj("line", name, eng_obj, length(data_math["branch"])+1)

            _branch_obj = Dict{String,Any}(
                "name" => "_virtual_branch.switch.$name",
                "source_id" => "_virtual_branch.switch.$name",
                # "f_bus" => bus_obj["bus_i"],  # TODO enable real switches
                "f_bus" => data_math["bus_lookup"][eng_obj["f_bus"]],
                "t_bus" => data_math["bus_lookup"][eng_obj["t_bus"]],
                "br_r" => _impedance_conversion(data_eng, eng_obj, "rs"),
                "br_x" => _impedance_conversion(data_eng, eng_obj, "xs"),
                "g_fr" => _admittance_conversion(data_eng, eng_obj, "g_fr"),
                "g_to" => _admittance_conversion(data_eng, eng_obj, "g_to"),
                "b_fr" => _admittance_conversion(data_eng, eng_obj, "b_fr"),
                "b_to" => _admittance_conversion(data_eng, eng_obj, "b_to"),
                "angmin" => fill(-60.0, nphases),
                "angmax" => fill( 60.0, nphases),
                "transformer" => false,
                "shift" => zeros(nphases),
                "tap" => ones(nphases),
                "switch" => false,
                "br_status" => 1,
            )

            merge!(branch_obj, _branch_obj)

            if kron_reduced
                @assert(all(eng_obj["f_connections"].==eng_obj["t_connections"]), "Kron reduction is only supported if f_connections is the same as t_connections.")
                filter = _kron_reduce_branch!(branch_obj,
                    ["br_r", "br_x"], ["g_fr", "b_fr", "g_to", "b_to"],
                    eng_obj["f_connections"], kr_neutral
                )
                _apply_filter!(branch_obj, ["angmin", "angmax", "tap", "shift"], filter)
                connections = eng_obj["f_connections"][filter]
                _pad_properties!(branch_obj, ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to", "angmin", "angmax", "tap", "shift"], connections, kr_phases)
            else
                branch_obj["f_connections"] = eng_obj["f_connections"]
                branch_obj["f_connections"] = eng_obj["t_connections"]
            end

            data_math["branch"]["$(branch_obj["index"])"] = branch_obj

            # build switch
            switch_obj = Dict{String,Any}(
                "name" => name,
                "source_id" => eng_obj["source_id"],
                "f_bus" => data_math["bus_lookup"][eng_obj["f_bus"]],
                "t_bus" => bus_obj["bus_i"],
                "status" => eng_obj["status"],
                "index" => length(data_math["switch"])+1
            )

            # data_math["switch"]["$(switch_obj["index"])"] = switch_obj

            map_to = ["branch.$(branch_obj["index"])"]
            # map_to = [map_to, "bus.$(bus_obj["index"])", "branch.$(branch_obj["index"])"]  # TODO enable real switches
        end

        # data_math["switch"]["$(math_obj["index"])"] = math_obj  # TODO enable real switches

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => name,
            :to => map_to,
            :unmap_function => :_map_math2eng_switch!,
        )

    end
end


"converts engineering line reactors into mathematical branches"
function _map_eng2math_line_reactor!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    # TODO support line reactors natively, currently treated like branches
    for (name, eng_obj) in get(data_eng, "line_reactor", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("line", name, eng_obj, length(data_math["branch"])+1)
        math_obj["name"] = "_virtual_branch.$(eng_obj["source_id"])"

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

        math_obj["angmin"] = fill(-60.0, nphases)
        math_obj["angmax"] = fill( 60.0, nphases)

        math_obj["transformer"] = false
        math_obj["shift"] = zeros(nphases)
        math_obj["tap"] = ones(nphases)

        f_bus = data_eng["bus"][eng_obj["f_bus"]]
        t_bus = data_eng["bus"][eng_obj["t_bus"]]

        if kron_reduced
            @assert(all(eng_obj["f_connections"].==eng_obj["t_connections"]), "Kron reduction is only supported if f_connections is the same as t_connections.")
            filter = _kron_reduce_branch!(math_obj,
                ["br_r", "br_x"], ["g_fr", "b_fr", "g_to", "b_to"],
                eng_obj["f_connections"], kr_neutral
            )
            _apply_filter!(math_obj, ["angmin", "angmax", "tap", "shift"], filter)
            connections = eng_obj["f_connections"][filter]
            _pad_properties!(math_obj, ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to", "angmin", "angmax", "tap", "shift"], connections, kr_phases)
        else
            math_obj["f_connections"] = eng_obj["f_connections"]
            math_obj["t_connections"] = eng_obj["t_connections"]
        end

        math_obj["switch"] = false

        math_obj["br_status"] = eng_obj["status"]

        data_math["branch"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => name,
            :to => "branch.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_line_reactor!,
        )

        # time series
        # TODO
        for (fr, to) in zip(["status"], ["status"])
            if haskey(eng_obj, "$(fr)_time_series")
                time_series = data_eng["time_series"][eng_obj["$(fr)_time_series"]]
                _parse_time_series_parameter!(data_math, time_series, eng_obj[fr], "branch", "$(math_obj["index"])", to)
            end
        end
    end
end


"converts engineering generic shunt components into mathematical shunt components"
function _map_eng2math_shunt!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "shunt", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("shunt", name, eng_obj, length(data_math["shunt"])+1)

        # TODO change to new capacitor shunt calc logic
        math_obj["shunt_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        math_obj["gs"] = eng_obj["gs"]
        math_obj["bs"] = eng_obj["bs"]

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

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => name,
            :to => "shunt.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_capacitor!,
        )

        # time series
        for (fr, (f, to)) in _time_series_parameters["shunt"]
            if haskey(eng_obj, "$(fr)_time_series")
                time_series = data_eng["time_series"][eng_obj["$(fr)_time_series"]]
                _parse_time_series_parameter!(data_math, time_series, eng_obj[fr], "shunt", "$(math_obj["index"])", to)
            end
        end
    end
end


"converts engineering shunt capacitors into mathematical shunts"
function _map_eng2math_shunt_capacitor!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "shunt_capacitor", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("shunt_capacitor", name, eng_obj, length(data_math["shunt"])+1)

        math_obj["shunt_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        math_obj["gs"] = zeros(size(eng_obj["bs"]))

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

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => name,
            :to => "shunt.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_shunt_capacitor!,
        )

        # time series
        # TODO
        for (fr, (f, to)) in _time_series_parameters["shunt_capacitor"]
            if haskey(eng_obj, "$(fr)_time_series")
                time_series = data_eng["time_series"][eng_obj["$(fr)_time_series"]]
                _parse_time_series_parameter!(data_math, time_series, eng_obj[fr], "shunt", "$(math_obj["index"])", to)
            end
        end
    end
end


"converts engineering shunt reactors into mathematical shunts"
function _map_eng2math_shunt_reactor!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "shunt_reactor", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("shunt_reactor", name, eng_obj, length(data_math["shunt"])+1)

        connections = eng_obj["connections"]
        nconductors = data_math["conductors"]

        math_obj["shunt_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        math_obj["gs"] = fill(0.0, size(eng_obj["bs"])...)

        if kron_reduced
            if eng_obj["configuration"] == "wye"
                _pad_properties!(math_obj, ["gs", "bs"], connections[1:end-1], kr_phases)
            else
                _pad_properties!(math_obj, ["gs", "bs"], connections, kr_phases)
            end
        end

        data_math["shunt"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => name,
            :to => "shunt.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_shunt_reactor!,
        )

        # time series
        # TODO
        for (fr, (f, to)) in _time_series_parameters["shunt_reactor"]
            if haskey(eng_obj, "$(fr)_time_series")
                time_series = data_eng["time_series"][eng_obj["$(fr)_time_series"]]
                _parse_time_series_parameter!(data_math, time_series, eng_obj[fr], "shunt", "$(math_obj["index"])", to)
            end
        end
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
            if math_obj["configuration"]=="wye"
                @assert(connections[end]==kr_neutral)
                _pad_properties!(math_obj, ["pd", "qd"], connections[connections.!=kr_neutral], kr_phases)
            else
                _pad_properties_delta!(math_obj, ["pd", "qd"], connections, kr_phases)
            end
        else
            math_obj["connections"] = connections
        end

        math_obj["vnom_kv"] = eng_obj["vnom"]

        data_math["load"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => name,
            :to => "load.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_load!,
        )

        # time series
        for (fr, (f, to)) in _time_series_parameters["load"]
            if haskey(eng_obj, "$(fr)_time_series")
                time_series = data_eng["time_series"][eng_obj["$(fr)_time_series"]]
                _parse_time_series_parameter!(data_math, time_series, eng_obj[fr], "load", "$(math_obj["index"])", to)
            end
        end
    end
end


"converts engineering generators into mathematical generators"
function _map_eng2math_generator!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "generator", Dict{String,Any}())
        math_obj = _init_math_obj("generator", name, eng_obj, length(data_math["gen"])+1)

        phases = eng_obj["phases"]
        connections = eng_obj["connections"]
        nconductors = data_math["conductors"]

        math_obj["gen_bus"] = data_math["bus_lookup"][eng_obj["bus"]]
        math_obj["gen_status"] = eng_obj["status"]

        math_obj["vg"] = eng_obj["vg"] ./ data_math["basekv"]

        math_obj["qmin"] = eng_obj["qg_lb"]
        math_obj["qmax"] = eng_obj["qg_ub"]

        math_obj["pmax"] = eng_obj["pg"]
        math_obj["pmin"] = zeros(phases)

        _add_gen_cost_model!(math_obj, eng_obj)

        math_obj["configuration"] = eng_obj["configuration"]

        if kron_reduced
            if math_obj["configuration"]=="wye"
                @assert(connections[end]==kr_neutral)
                _pad_properties!(math_obj, ["pg", "qg", "vg", "pmin", "pmax", "qmin", "qmax"], connections[1:end-1], kr_phases)
            else
                _pad_properties_delta!(math_obj, ["pg", "qg", "vg", "pmin", "pmax", "qmin", "qmax"], connections, kr_phases)
            end
        else
            math_obj["connections"] = connections
        end

        # if PV generator mode convert attached bus to PV bus
        if eng_obj["control_mode"] == 3
            data_math["bus"]["$(data_math["bus_lookup"][eng_obj["bus"]])"]["bus_type"] = 2
        end

        data_math["gen"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => name,
            :to => "gen.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_generator!,
        )

        # time series
        # TODO
        for (fr, (f, to)) in _time_series_parameters["generator"]
            if haskey(eng_obj, "$(fr)_time_series")
                time_series = data_eng["time_series"][eng_obj["$(fr)_time_series"]]
                _parse_time_series_parameter!(data_math, time_series, eng_obj[fr], "gen", "$(math_obj["index"])", to)
            end
        end
    end
end


"converts engineering solar components into mathematical generators"
function _map_eng2math_solar!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "solar", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("solar", name, eng_obj, length(data_math["gen"])+1)

        connections = eng_obj["connections"]
        nconductors = data_math["conductors"]

        math_obj["gen_bus"] = data_math["bus_lookup"][eng_obj["bus"]]
        math_obj["gen_status"] = eng_obj["status"]

        math_obj["pg"] = eng_obj["kva"]
        math_obj["qg"] = eng_obj["kvar"]
        math_obj["vg"] = eng_obj["kv"]

        math_obj["pmin"] = get(eng_obj, "minkva", zeros(size(eng_obj["kva"])))
        math_obj["pmax"] = get(eng_obj, "maxkva", eng_obj["kva"])

        math_obj["qmin"] =  get(eng_obj, "minkvar", -eng_obj["kvar"])
        math_obj["qmax"] =  get(eng_obj, "maxkvar",  eng_obj["kvar"])

        _add_gen_cost_model!(math_obj, eng_obj)

        if kron_reduced
            if math_obj["configuration"]=="wye"
                @assert(connections[end]==kr_neutral)
                _pad_properties!(math_obj, ["pg", "qg", "vg", "pmin", "pmax", "qmin", "qmax"], connections[1:end-1], kr_phases)
            else
                _pad_properties_delta!(math_obj, ["pg", "qg", "vg", "pmin", "pmax", "qmin", "qmax"], connections, kr_phases)
            end
        else
            math_obj["connections"] = connections
        end

        data_math["gen"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => "solar.$name",
            :to => "gen.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_solar!,
        )

        # time series
        # TODO
        for (fr, (f, to)) in _time_series_parameters["solar"]
            if haskey(eng_obj, "$(fr)_time_series")
                time_series = data_eng["time_series"][eng_obj["$(fr)_time_series"]]
                _parse_time_series_parameter!(data_math, time_series, eng_obj, "gen", "$(math_obj["index"])", fr, to, f)
            end
        end
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
        math_obj["energy"] = eng_obj["energy"] * data_eng["settings"]["v_var_scalar"] / 1e6
        math_obj["energy_rating"] = eng_obj["energy_ub"] * data_eng["settings"]["v_var_scalar"] / 1e6
        math_obj["charge_rating"] = eng_obj["charge_ub"] * data_eng["settings"]["v_var_scalar"] / 1e6
        math_obj["discharge_rating"] = eng_obj["discharge_ub"] * data_eng["settings"]["v_var_scalar"] / 1e6
        math_obj["charge_efficiency"] = eng_obj["charge_efficiency"] / 100.0
        math_obj["discharge_efficiency"] = eng_obj["discharge_efficiency"] / 100.0
        math_obj["thermal_rating"] = eng_obj["cm_ub"] .* data_eng["settings"]["v_var_scalar"] ./ 1e6
        math_obj["qmin"] = eng_obj["qs_lb"] .* data_eng["settings"]["v_var_scalar"] ./ 1e6
        math_obj["qmax"] =  eng_obj["qs_ub"] .* data_eng["settings"]["v_var_scalar"] ./ 1e6
        math_obj["r"] = eng_obj["rs"]
        math_obj["x"] = eng_obj["xs"]
        math_obj["p_loss"] = eng_obj["pex"] .* data_eng["settings"]["v_var_scalar"] ./ 1e6
        math_obj["q_loss"] = eng_obj["qex"] .* data_eng["settings"]["v_var_scalar"] ./ 1e6

        math_obj["ps"] = get(eng_obj, "ps", zeros(size(eng_obj["cm_ub"])))
        math_obj["qs"] = get(eng_obj, "qs", zeros(size(eng_obj["cm_ub"])))

        if kron_reduced
            _pad_properties!(math_obj, ["thermal_rating", "qmin", "qmax", "r", "x", "ps", "qs"], connections, kr_phases)
        end

        data_math["storage"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => name,
            :to => "storage.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_storage!,
        )

        # time series
        # TODO
        for (fr, (f, to)) in _time_series_parameters["storage"]
            if haskey(eng_obj, "$(fr)_time_series")
                time_series = data_eng["time_series"][eng_obj["$(fr)_time_series"]]
                _parse_time_series_parameter!(data_math, time_series, eng_obj[fr], "storage", "$(math_obj["index"])", to)
            end
        end
    end
end


"converts engineering voltage sources into mathematical generators and (if needed) impedance branches to represent the loss model"
function _map_eng2math_voltage_source!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1,2,3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "voltage_source", Dict{Any,Any}())
        nconductors = length(eng_obj["vm"])

        math_obj = _init_math_obj("voltage_source", name, eng_obj, length(data_math["gen"])+1)

        math_obj["name"] = "_virtual_gen.voltage_source.$name"
        math_obj["gen_bus"] = gen_bus = data_math["bus_lookup"][eng_obj["bus"]]
        math_obj["gen_status"] = eng_obj["status"]
        math_obj["pg"] = fill(0.0, nconductors)
        math_obj["qg"] = fill(0.0, nconductors)
        math_obj["configuration"] = "wye"
        math_obj["source_id"] = "_virtual_gen.$(eng_obj["source_id"])"

        _add_gen_cost_model!(math_obj, eng_obj)

        map_to = "gen.$(math_obj["index"])"

        if haskey(eng_obj, "rs") && haskey(eng_obj, "xs")
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
                "vmax" => [eng_obj["vm"]..., 0.0],
                "basekv" => data_math["basekv"]
            )

            if kron_reduced
                bus_obj["terminals"] = bus_obj["terminals"][1:end-1]
                bus_obj["grounded"] = bus_obj["grounded"][1:end-1]
                bus_obj["vm"] = bus_obj["vm"][1:end-1]
                bus_obj["va"] = bus_obj["va"][1:end-1]
                bus_obj["vmin"] = bus_obj["vmin"][1:end-1]
                bus_obj["vmax"] = bus_obj["vmax"][1:end-1]
            end

            math_obj["gen_bus"] = gen_bus = bus_obj["bus_i"]

            data_math["bus"]["$(bus_obj["index"])"] = bus_obj

            branch_obj = Dict{String,Any}(
                "name" => "_virtual_branch.voltage_source.$name",
                "source_id" => "_virtual_branch.$(eng_obj["source_id"])",
                "f_bus" => bus_obj["bus_i"],
                "t_bus" => data_math["bus_lookup"][eng_obj["bus"]],
                "f_connections" => collect(1:nconductors),
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
        end

        data_math["gen"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => name,
            :to => map_to,
            :unmap_function => :_map_math2eng_voltage_source!,
        )
    end
end
