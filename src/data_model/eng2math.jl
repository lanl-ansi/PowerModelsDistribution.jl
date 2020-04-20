import LinearAlgebra: diagm


const _1to1_maps = Dict{String,Vector{String}}(
    "bus" => ["vm", "va", "vmin", "vmax", "terminals", "phases", "neutral"],
    "load" => ["model", "configuration", "status", "source_id", "connections"],
    "shunt_capacitor" => ["status", "source_id"],
    "series_capacitor" => [],
    "shunt" => ["status", "source_id"],
    "shunt_reactor" => ["status", "source_id"],
    "generator" => ["source_id", "configuration"],
    "solar" => ["source_id", "configuration"],
    "storage" => ["status", "source_id"],
    "line" => ["source_id"],
    "line_reactor" => ["source_id"],
    "switch" => ["source_id", "state", "status"],
    "line_reactor" => ["source_id"],
    "transformer" => ["source_id"],
    "voltage_source" => ["source_id"],
)

const _node_elements = ["load", "capacitor", "shunt_reactor", "generator", "solar", "storage", "vsource"]

const _edge_elements = ["line", "switch", "transformer", "line_reactor", "series_capacitor"]


""
function _map_eng2math(data_eng; kron_reduced::Bool=true, lossless::Bool=false, adjust_bounds::Bool=true)
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

    _map_eng2math_bus!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_load!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_shunt_capacitor!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_shunt_reactor!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_shunt!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_generator!(data_math, data_eng; kron_reduced=kron_reduced, lossless=lossless)
    _map_eng2math_solar!(data_math, data_eng; kron_reduced=kron_reduced, lossless=lossless)
    _map_eng2math_storage!(data_math, data_eng; kron_reduced=kron_reduced, lossless=lossless)
    _map_eng2math_voltage_source!(data_math, data_eng; kron_reduced=kron_reduced, lossless=lossless)

    _map_eng2math_line!(data_math, data_eng; kron_reduced=kron_reduced)
    # _map_eng2math_series_capacitor(data_math, data_eng; kron_reduced=kron_reduced)  # TODO
    _map_eng2math_line_reactor!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_switch!(data_math, data_eng; kron_reduced=kron_reduced, lossless=lossless)

    _map_eng2math_transformer!(data_math, data_eng; kron_reduced=kron_reduced)

    if !adjust_bounds
        # _find_new_bounds(data_math) # TODO
    end

    return data_math
end


""
function _init_math_obj(obj_type::String, eng_obj::Dict{String,<:Any}, index::Int)::Dict{String,Any}
    math_obj = Dict{String,Any}()

    for key in _1to1_maps[obj_type]
        if haskey(eng_obj, key)
            math_obj[key] = eng_obj[key]
        end
    end

    math_obj["index"] = index

    return math_obj
end


""
function _init_base_components!(data_math::Dict{String,<:Any})
    for key in ["bus", "load", "shunt", "gen", "branch", "switch", "transformer", "storage", "dcline"]
        if !haskey(data_math, key)
            data_math[key] = Dict{String,Any}()
        end
    end
end


""
function _init_lookup!(data_math::Dict{String,<:Any})


    for key in keys(_1to1_maps)
        if !haskey(data_math["lookup"], key)
            data_math["lookup"][key] = Dict{Any,Int}()
        end
    end
end


""
function _map_eng2math_bus!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "bus", Dict{String,Any}())
        # TODO fix vnom
        terminals = eng_obj["terminals"]
        nconductors = data_math["conductors"]

        math_obj = _init_math_obj("bus", eng_obj, length(data_math["bus"])+1)
        math_obj["name"] = name

        math_obj["bus_i"] = math_obj["index"]
        math_obj["bus_type"] = eng_obj["status"] == 1 ? 1 : 4

        # take care of grounding; convert to shunt if lossy
        grounded_perfect, shunts = _convert_grounding(eng_obj["terminals"], eng_obj["grounded"], eng_obj["rg"], eng_obj["xg"])
        math_obj["grounded"] = grounded_perfect
        to_sh = []
        for (sh_connections, sh_y) in shunts
            sh_index = length(data_math["shunt"])+1
            data_math["shunt"]["$sh_index"] = Dict(
                "index" => sh_index,
                "shunt_bus" => bus_index,
                "connections" => sh_connections,
                "gs" => real.(sh_y),
                "bs" => real.(sh_y),
            )
            push!(to_sh, "shunt.$sh_index")
        end

        math_obj["vmin"] = get(eng_obj, "vmin", fill(0.0, length(terminals)))
        math_obj["vmax"] = get(eng_obj, "vmax", fill(Inf, length(terminals)))

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
    end
end


""
function _map_eng2math_load!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "load", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("load", eng_obj, length(data_math["load"])+1)
        math_obj["name"] = name

        connections = eng_obj["connections"]
        nconductors = data_math["conductors"]

        math_obj["load_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        math_obj["status"] = eng_obj["status"]

        math_obj["pd"] = eng_obj["pd_nom"]
        math_obj["qd"] = eng_obj["qd_nom"]

        math_obj["configuration"] = eng_obj["configuration"]

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
    end
end


""
function _map_eng2math_shunt_capacitor!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "shunt_capacitor", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("shunt_capacitor", eng_obj, length(data_math["shunt"])+1)

        # TODO change to new capacitor shunt calc logic
        math_obj["name"] = name
        math_obj["shunt_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        f_terminals = eng_obj["f_connections"]
        t_terminals = eng_obj["t_connections"]

        # TODO change in OpenDS parser
        b = diag(eng_obj["bs"])

        # convert to a shunt matrix
        terminals, B = _calc_shunt(f_terminals, t_terminals, b)

        # if one terminal is ground (0), reduce shunt addmittance matrix
        terminals, B = _calc_ground_shunt_admittance_matrix(terminals, B, 0)

        math_obj["bs"] = B
        math_obj["gs"] = zeros(size(B))

        if kron_reduced
            filter = _kron_reduce_branch!(math_obj,
                Vector{String}([]), ["gs", "bs"],
                eng_obj["f_connections"], kr_neutral
            )
            connections = eng_obj["f_connections"][filter]
            _pad_properties!(math_obj, ["gs", "bs"], connections, kr_phases)
        else
            math_obj["connections"] = terminals
        end

        data_math["shunt"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => name,
            :to => "shunt.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_shunt_capacitor!,
        )
    end
end


""
function _map_eng2math_shunt!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "shunt", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("shunt", eng_obj, length(data_math["shunt"])+1)

        # TODO change to new capacitor shunt calc logic
        math_obj["name"] = name
        math_obj["shunt_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        math_obj["gs"] = eng_obj["g_sh"]
        math_obj["bs"] = eng_obj["b_sh"]

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
    end
end


""
function _map_eng2math_shunt_reactor!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "shunt_reactor", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("shunt_reactor", eng_obj, length(data_math["shunt"])+1)

        nphases = eng_obj["phases"]
        connections = eng_obj["connections"]
        nconductors = data_math["conductors"]

        Gcap = sum(eng_obj["kvar"]) / (nphases * 1e3 * (data_math["basekv"])^2)

        math_obj["name"] = name
        math_obj["shunt_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        math_obj["gs"] = fill(0.0, nphases, nphases)
        math_obj["bs"] = diagm(0=>fill(Gcap, nphases))

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
    end
end


""
function _map_eng2math_generator!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4, lossless::Bool=false)
    for (name, eng_obj) in get(data_eng, "generator", Dict{String,Any}())
        math_obj = _init_math_obj("generator", eng_obj, length(data_math["gen"])+1)

        phases = eng_obj["phases"]
        connections = eng_obj["connections"]
        nconductors = data_math["conductors"]

        math_obj["gen_bus"] = data_math["bus_lookup"][eng_obj["bus"]]
        math_obj["name"] = name
        math_obj["gen_status"] = eng_obj["status"]

        math_obj["pg"] = eng_obj["kw"]
        math_obj["qg"] = eng_obj["kvar"]
        math_obj["vg"] = eng_obj["kv"] ./ data_math["basekv"]

        math_obj["qmin"] = eng_obj["kvar_min"]
        math_obj["qmax"] = eng_obj["kvar_max"]

        math_obj["pmax"] = eng_obj["kw"]
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
        if eng_obj["control_model"] == 3
            data_math["bus"]["$(data_math["bus_lookup"][eng_obj["bus"]])"]["bus_type"] = 2
        end

        data_math["gen"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => name,
            :to => "gen.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_generator!,
        )
    end
end


""
function _map_eng2math_solar!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4, lossless::Bool=false)
    for (name, eng_obj) in get(data_eng, "solar", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("solar", eng_obj, length(data_math["gen"])+1)

        connections = eng_obj["connections"]
        nconductors = data_math["conductors"]

        math_obj["name"] = name
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
    end
end


""
function _map_eng2math_storage!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4, lossless::Bool=false)
    for (name, eng_obj) in get(data_eng, "storage", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("storage", eng_obj, length(data_math["storage"])+1)

        connections = eng_obj["connections"]
        nconductors = data_math["conductors"]

        math_obj["name"] = name
        math_obj["storage_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        math_obj["energy"] = eng_obj["kwhstored"] / 1e3
        math_obj["energy_rating"] = eng_obj["kwhrated"] / 1e3
        math_obj["charge_rating"] = eng_obj["%charge"] * eng_obj["kwrated"] / 1e3 / 100.0
        math_obj["discharge_rating"] = eng_obj["%discharge"] * eng_obj["kwrated"] / 1e3 / 100.0
        math_obj["charge_efficiency"] = eng_obj["%effcharge"] / 100.0
        math_obj["discharge_efficiency"] = eng_obj["%effdischarge"] / 100.0
        math_obj["thermal_rating"] = eng_obj["kva"] ./ 1e3
        math_obj["qmin"] = -eng_obj["kvar"] ./ 1e3
        math_obj["qmax"] =  eng_obj["kvar"] ./ 1e3
        math_obj["r"] = eng_obj["%r"] ./ 100.0
        math_obj["x"] = eng_obj["%x"] ./ 100.0
        math_obj["p_loss"] = eng_obj["%idlingkw"] .* eng_obj["kwrated"] ./ 1e3
        math_obj["q_loss"] = eng_obj["%idlingkvar"] * sum(eng_obj["kvar"]) / 1e3

        math_obj["ps"] = get(eng_obj, "ps", zeros(size(eng_obj["kva"])))
        math_obj["qs"] = get(eng_obj, "qs", zeros(size(eng_obj["kva"])))

        if kron_reduced
            _pad_properties!(math_obj, ["thermal_rating", "qmin", "qmax", "r", "x", "ps", "qs"], connections, kr_phases)
        end

        data_math["storage"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => name,
            :to => "storage.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_storage!,
        )
    end
end


""
function _map_eng2math_line!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    for (name, eng_obj) in get(data_eng, "line", Dict{Any,Dict{String,Any}}())
        _apply_linecode!(eng_obj, data_eng)

        math_obj = _init_math_obj("line", eng_obj, length(data_math["branch"])+1)
        math_obj["name"] = name

        nphases = length(eng_obj["f_connections"])
        nconductors = data_math["conductors"]

        math_obj["f_bus"] = data_math["bus_lookup"][eng_obj["f_bus"]]
        math_obj["t_bus"] = data_math["bus_lookup"][eng_obj["t_bus"]]

        math_obj["br_r"] = eng_obj["rs"] * eng_obj["length"]
        math_obj["br_x"] = eng_obj["xs"] * eng_obj["length"]

        math_obj["g_fr"] = (2.0 * pi * data_eng["settings"]["base_frequency"] * eng_obj["g_fr"] * eng_obj["length"] / 1e9)
        math_obj["g_to"] = (2.0 * pi * data_eng["settings"]["base_frequency"] * eng_obj["g_to"] * eng_obj["length"] / 1e9)

        math_obj["b_fr"] = (2.0 * pi * data_eng["settings"]["base_frequency"] * eng_obj["b_fr"] * eng_obj["length"] / 1e9)
        math_obj["b_to"] = (2.0 * pi * data_eng["settings"]["base_frequency"] * eng_obj["b_to"] * eng_obj["length"] / 1e9)

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
    end
end


""
function _map_eng2math_line_reactor!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4)
    # TODO support line reactors natively, currently treated like branches
    for (name, eng_obj) in get(data_eng, "line_reactor", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("line", eng_obj, length(data_math["branch"])+1)
        math_obj["name"] = "_virtual_branch.$(eng_obj["source_id"])"

        nphases = length(eng_obj["f_connections"])
        nconductors = data_math["conductors"]

        math_obj["f_bus"] = data_math["bus_lookup"][eng_obj["f_bus"]]
        math_obj["t_bus"] = data_math["bus_lookup"][eng_obj["t_bus"]]

        math_obj["br_r"] = eng_obj["rs"] * eng_obj["length"]
        math_obj["br_x"] = eng_obj["xs"] * eng_obj["length"]

        math_obj["g_fr"] = (2.0 * pi * data_eng["settings"]["base_frequency"] * eng_obj["g_fr"] * eng_obj["length"] / 1e9)
        math_obj["g_to"] = (2.0 * pi * data_eng["settings"]["base_frequency"] * eng_obj["g_to"] * eng_obj["length"] / 1e9)

        math_obj["b_fr"] = (2.0 * pi * data_eng["settings"]["base_frequency"] * eng_obj["b_fr"] * eng_obj["length"] / 1e9)
        math_obj["b_to"] = (2.0 * pi * data_eng["settings"]["base_frequency"] * eng_obj["b_to"] * eng_obj["length"] / 1e9)

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
    end
end


""
function _map_eng2math_switch!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1, 2, 3], kr_neutral::Int=4, lossless::Bool=false)
    # TODO enable real switches (right now only using vitual lines)
    for (name, eng_obj) in get(data_eng, "switch", Dict{Any,Dict{String,Any}}())
        nphases = length(eng_obj["f_connections"])
        nconductors = data_math["conductors"]

        if lossless
            math_obj = _init_math_obj("switch", eng_obj, length(data_math["switch"])+1)
            math_obj["name"] = name

            math_obj["f_bus"] = data_math["bus_lookup"][eng_obj["f_bus"]]
            math_obj["t_bus"] = data_math["bus_lookup"][eng_obj["f_bus"]]

            data_math["switch"]["$(math_obj["index"])"] = math_obj

            data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
                :from => name,
                :to => "switch.$(math_obj["index"])",
                :unmap_function => :_map_math2eng_switch!,
            )
        else
            # build virtual bus
            f_bus = data_math["bus"]["$(data_math["bus_lookup"][eng_obj["f_bus"]])"]

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

            # data_math["bus"]["$(bus_obj["index"])"] = bus_obj

            # build virtual branch
            if haskey(eng_obj, "linecode")
                linecode = data_eng["linecode"][eng_obj["linecode"]]

                for property in ["rs", "xs", "g_fr", "g_to", "b_fr", "b_to"]
                    if !haskey(eng_obj, property) && haskey(linecode, property)
                        eng_obj[property] = linecode[property]
                    end
                end
            end

            branch_obj = _init_math_obj("line", eng_obj, length(data_math["branch"])+1)

            Zbase = (data_math["basekv"])^2 / (data_math["baseMVA"])
            Zbase = 1

            _branch_obj = Dict{String,Any}(
                "name" => "_virtual_branch.switch.$name",
                "source_id" => "_virtual_branch.switch.$name",
                # "f_bus" => bus_obj["bus_i"],  # TODO enable real switches
                "f_bus" => data_math["bus_lookup"][eng_obj["f_bus"]],
                "t_bus" => data_math["bus_lookup"][eng_obj["t_bus"]],
                "br_r" => eng_obj["rs"] * eng_obj["length"] / Zbase,
                "br_x" => eng_obj["xs"] * eng_obj["length"] / Zbase,
                "g_fr" => Zbase * (2.0 * pi * data_eng["settings"]["base_frequency"] * eng_obj["g_fr"] * eng_obj["length"] / 1e9),
                "g_to" => Zbase * (2.0 * pi * data_eng["settings"]["base_frequency"] * eng_obj["g_to"] * eng_obj["length"] / 1e9),
                "b_fr" => Zbase * (2.0 * pi * data_eng["settings"]["base_frequency"] * eng_obj["b_fr"] * eng_obj["length"] / 1e9),
                "b_to" => Zbase * (2.0 * pi * data_eng["settings"]["base_frequency"] * eng_obj["b_to"] * eng_obj["length"] / 1e9),
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

            data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
                :from => name,
                # :to_id => ["switch.$(switch_obj["index"])", "bus.$(bus_obj["index"])", "branch.$(branch_obj["index"])"],  # TODO enable real switches
                :to => ["branch.$(branch_obj["index"])"],
                :unmap_function => :_map_math2eng_switch!,
            )
        end
    end
end


""
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

        nphases = length(eng_obj["fixed"][1])
        transformer_t_bus_w = _build_loss_model!(data_math, name, to_map, r_s, z_sc, y_sh, nphases=nphases, kron_reduced=kron_reduced)

        for w in 1:nrw
            # 2-WINDING TRANSFORMER
            # make virtual bus and mark it for reduction
            tm_nom = eng_obj["vnom"][w]
            # three-phase delta transformers need a correction of sqrt(3)
            if eng_obj["configuration"][w]=="delta" && length(eng_obj["connections"][w])==3
                tm_nom = eng_obj["vnom"][w]*sqrt(3)
            end
            # tm_nom = eng_obj["configuration"][w]=="delta" ? eng_obj["vnom"][w]*sqrt(3) : eng_obj["vnom"][w]


            transformer_2wa_obj = Dict{String,Any}(
                "name"          => "_virtual_transformer.$name.$w",
                "source_id"     => "_virtual_transformer.$(eng_obj["source_id"]).$w",
                "f_bus"         => data_math["bus_lookup"][eng_obj["bus"][w]],
                "t_bus"         => transformer_t_bus_w[w],
                "tm_nom"        => tm_nom,
                "f_connections" => eng_obj["connections"][w],
                "t_connections" => collect(1:nphases+1),
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


""
function _map_eng2math_voltage_source!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true, kr_phases::Vector{Int}=[1,2,3], kr_neutral::Int=4, lossless::Bool=false)
    for (name, eng_obj) in get(data_eng, "voltage_source", Dict{Any,Any}())
        nconductors = length(eng_obj["vm"])

        if lossless
            #TODO add unreduced logic here
            math_obj = _init_math_obj("voltage_source", eng_obj, length(data_math["gen"])+1)

            math_obj["name"] = name
            math_obj["gen_bus"] = data_math["bus_lookup"][eng_obj["bus"]]
            math_obj["gen_status"] = eng_obj["status"]
            math_obj["pg"] = fill(0.0, nconductors)
            math_obj["qg"] = fill(0.0, nconductors)
            math_obj["configuration"] = "wye"

            _add_gen_cost_model!(math_obj, eng_obj)

            data_math["gen"]["$(math_obj["index"])"] = math_obj

            data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
                :from => name,
                :to => "gen.$(math_obj["index"])",
                :unmap_function => :_map_math2eng_voltage_source!,
            )
        else
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

            data_math["bus"]["$(bus_obj["index"])"] = bus_obj

            gen_obj = Dict{String,Any}(
                "gen_bus" => bus_obj["bus_i"],
                "connections" => collect(1:nconductors+1),
                "name" => "_virtual_gen.voltage_source.$name",
                "gen_status" => eng_obj["status"],
                "pg" => fill(0.0, nconductors),
                "qg" => fill(0.0, nconductors),
                "model" => 2,
                "startup" => 0.0,
                "shutdown" => 0.0,
                "ncost" => 3,
                "cost" => [0.0, 1.0, 0.0],
                "configuration" => "wye",
                "index" => length(data_math["gen"]) + 1,
                "source_id" => "_virtual_gen.$(eng_obj["source_id"])"
            )

            data_math["gen"]["$(gen_obj["index"])"] = gen_obj

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
                "br_r" => eng_obj["rs"],
                "br_x" => eng_obj["xs"],
                "g_fr" => zeros(nconductors, nconductors),
                "g_to" => zeros(nconductors, nconductors),
                "b_fr" => zeros(nconductors, nconductors),
                "b_to" => zeros(nconductors, nconductors),
                "index" => length(data_math["branch"])+1
            )

            data_math["branch"]["$(branch_obj["index"])"] = branch_obj

            data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
                :from => name,
                :to => ["gen.$(gen_obj["index"])", "bus.$(bus_obj["index"])", "branch.$(branch_obj["index"])"],
                :unmap_function => :_map_math2eng_voltage_source!,
            )
        end
    end
end
