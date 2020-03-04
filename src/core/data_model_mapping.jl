import LinearAlgebra: diagm


const _1to1_maps = Dict{String,Vector{String}}(
    "bus" => ["vm", "va", "vmin", "vmax"],
    "load" => ["model", "configuration", "status", "source_id"],
    "capacitor" => ["status"],
    "shunt_reactor" => ["status"],
    "generator" => ["model", "startup", "shutdown", "ncost", "cost", "source_id"],
    "pvsystem" => ["model", "startup", "shutdown", "ncost", "cost", "source_id"],
    "storage" => ["status", "source_id"],
    "line" => ["source_id"],
    "switch" => ["source_id"],
    "line_reactor" => ["source_id"],
    "transformer" => ["source_id"],
    "voltage_source" => ["source_id"],
)

const _extra_eng_data = Dict{String,Vector{String}}(
    "root" => ["sourcebus", "files", "dss_options", "settings"],
    "bus" => ["grounded", "neutral", "awaiting_ground", "xg", "phases", "rg", "terminals"],
    "load" => [],
    "capacitor" => [],
    "shunt_reactor" => [],
    "generator" => ["control_model"],
    "pvsystem" => [],
    "storage" => [],
    "line" => ["f_connections", "t_connections", "linecode"],
    "switch" => [],
    "line_reactor" => [],
    "transformer" => [],
    "voltage_source" => [],
)

const _node_elements = ["load", "capacitor", "shunt_reactor", "generator", "pvsystem", "storage", "vsource"]

const _edge_elements = ["line", "switch", "transformer"]


""
function _map_eng2math(data_eng; kron_reduced::Bool=true)
    @assert get(data_eng, "data_model", "mathematical") == "engineering"

    data_math = Dict{String,Any}(
        "name" => data_eng["name"],
        "per_unit" => get(data_eng, "per_unit", false)
    )

    data_math["conductors"] = kron_reduced ? 3 : 4
    data_math["basekv"] = data_eng["settings"]["set_vbase_val"]
    data_math["baseMVA"] = data_eng["settings"]["set_sbase_val"]

    data_math["data_model"] = "mathematical"

    data_math["map"] = Dict{Int,Dict{Symbol,Any}}(
        1 => Dict{Symbol,Any}(
            :component_type => "root",
            :unmap_function => :_map_math2eng_root!,
            :extra => Dict{String,Any}((k,v) for (k,v) in data_eng if k in _extra_eng_data["root"])
        )
    )

    _init_base_components!(data_math)

    _map_eng2math_bus!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_load!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_capacitor!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_shunt_reactor!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_generator!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_pvsystem!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_storage!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_vsource!(data_math, data_eng; kron_reduced=kron_reduced)  # TODO

    _map_eng2math_line!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_switch!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_transformer!(data_math, data_eng; kron_reduced=kron_reduced)  # TODO

    _map_eng2math_sourcebus!(data_math, data_eng; kron_reduced=kron_reduced)

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
function _map_eng2math_bus!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    for (name, eng_obj) in get(data_eng, "bus", Dict{String,Any}())
        # TODO fix vnom
        phases = get(eng_obj, "phases", [1, 2, 3])
        neutral = get(eng_obj, "neutral", 4)
        terminals = eng_obj["terminals"]

        @assert all(t in [phases..., neutral] for t in terminals)

        math_obj = _init_math_obj("bus", eng_obj, length(data_math["bus"])+1)

        math_obj["bus_i"] = math_obj["index"]
        math_obj["bus_type"] = eng_obj["status"] == 1 ? 1 : 4

        math_obj["vm"] = get(eng_obj, "vm", fill(1.0, length(terminals)))
        math_obj["va"] = get(eng_obj, "va", [_wrap_to_180(-rad2deg(2*pi/length(phases)*(i-1))) for i in terminals])

        math_obj["vmin"] = fill(0.0, length(terminals))
        math_obj["vmax"] = fill(Inf, length(terminals))

        math_obj["base_kv"] = data_eng["settings"]["set_vbase_val"]

        if kron_reduced
            for key in ["vm", "va", "vmin", "vmax"]
                if haskey(math_obj, key)
                    math_obj[key] = math_obj[key][terminals.!=neutral]
                end
            end
        end

        data_math["bus"]["$(math_obj["index"])"] = math_obj

        if !haskey(data_math, "bus_lookup")
            data_math["bus_lookup"] = Dict{Any,Int}()
        end

        data_math["bus_lookup"][name] = math_obj["index"]

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => "bus.$name",
            :to => ["bus.$(math_obj["index"])"],
            :unmap_function => :_map_math2eng_bus!,
            :kron_reduced => kron_reduced,
            :extra => Dict{String,Any}((k,v) for (k,v) in eng_obj if k in _extra_eng_data["bus"])
        )
    end
end


""
function _map_eng2math_load!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    # TODO add delta loads
    for (name, eng_obj) in get(data_eng, "load", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("load", eng_obj, length(data_math["load"])+1)

        phases = get(eng_obj, "phases", [1, 2, 3])
        conductors = collect(1:data_math["conductors"])

        if eng_obj["configuration"] == "wye"
            bus = data_eng["bus"][eng_obj["bus"]]
            # TODO add message for failure
            @assert length(bus["grounded"]) == 1 && bus["grounded"][1] == eng_obj["connections"][end]
        else
            # TODO add message for failure
            @assert all(load["connections"] .== phases)
        end

        math_obj["load_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        math_obj["pd"] = eng_obj["pd"] / 1e3
        math_obj["qd"] = eng_obj["qd"] / 1e3

        math_obj["vnom_kv"] = get(eng_obj, "vnom", data_eng["settings"]["set_vbase_val"]/sqrt(3)*1e3)

        _pad_properties!(math_obj, ["pd", "qd"], eng_obj["connections"], phases; neutral=data_eng["bus"][eng_obj["bus"]]["neutral"])

        data_math["load"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => "load.$name",
            :to => "load.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_load!,
            :kron_reduced => kron_reduced,
            :extra => Dict{String,Any}((k,v) for (k,v) in eng_obj if k in _extra_eng_data["load"])
        )
    end
end


""
function _map_eng2math_capacitor!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    for (name, eng_obj) in get(data_eng, "capacitor", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("capacitor", eng_obj, length(data_math["shunt"])+1)

        # TODO change to new capacitor shunt calc logic
        math_obj["shunt_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        nphases = eng_obj["phases"]

        vnom_ln = eng_obj["kv"]
        # 'kv' is specified as phase-to-phase for phases=2/3 (unsure for 4 and more)
        if eng_obj["phases"] > 1
            vnom_ln = vnom_ln/sqrt(3)
        end
        # 'kvar' is specified for all phases at once; we want the per-phase one, in MVar
        qnom = (eng_obj["kvar"]/1E3)/nphases
        # indexing qnom[1] is a dirty fix to support both kvar=[x] and kvar=x
        # TODO fix this in a clear way, in dss_structs.jl
        b_cap = qnom[1]/vnom_ln^2
        #  get the base addmittance, with a LN voltage base
        Sbase = 1 # not yet pmd_data["baseMVA"] because this is done in _PMs.make_per_unit
        Ybase_ln = Sbase/(data_math["basekv"]/sqrt(3))^2
        # now convent b_cap to per unit
        b_cap_pu = b_cap/Ybase_ln

        b = fill(b_cap_pu, nphases)
        N = length(b)

        if eng_obj["configuration"] == "wye"
            B = diagm(0=>b)
        else
            # create delta transformation matrix Md
            Md = diagm(0=>ones(N), 1=>-ones(N-1))
            Md[N,1] = -1
            B = Md'*diagm(0=>b)*Md
        end

        math_obj["gs"] = fill(0.0, nphases, nphases)
        math_obj["bs"] = B

        neutral = get(data_eng["bus"][eng_obj["bus"]], "neutral", 4)

        _pad_properties!(math_obj, ["gs", "bs"], eng_obj["f_terminals"], collect(1:data_math["conductors"]); neutral=neutral)

        data_math["shunt"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => "capacitor.$name",
            :to => "shunt.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_capacitor!,
            :kron_reduced => kron_reduced,
            :extra => Dict{String,Any}((k,v) for (k,v) in eng_obj if k in _extra_eng_data["capacitor"])
        )
    end
end


""
function _map_eng2math_shunt_reactor!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
end


""
function _map_eng2math_generator!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    for (name, eng_obj) in get(data_eng, "generator", Dict{String,Any}())
        math_obj = _init_math_obj("generator", eng_obj, length(data_math["gen"])+1)

        phases = eng_obj["phases"]

        math_obj["gen_bus"] = data_math["bus_lookup"][eng_obj["bus"]]
        math_obj["name"] = name
        math_obj["gen_status"] = eng_obj["status"]

        math_obj["pg"] = fill(eng_obj["kw"] / 1e3 / phases, phases)
        math_obj["qg"] = fill(eng_obj["kvar"] / 1e3 / phases, phases)
        math_obj["vg"] = fill(eng_obj["kv"] / data_math["basekv"])

        math_obj["qmin"] = fill(eng_obj["kvar_min"] / 1e3 / phases, phases)
        math_obj["qmax"] = fill(eng_obj["kvar_max"] / 1e3 / phases, phases)

        math_obj["pmax"] = fill(eng_obj["kw"] / 1e3 / phases, phases)
        math_obj["pmin"] = fill(0.0, phases)

        math_obj["conn"] = eng_obj["configuration"]

        # if PV generator mode convert attached bus to PV bus
        if eng_obj["control_model"] == 3
            data_math["bus"][data_math["bus_lookup"][eng_obj["bus"]]]["bus_type"] = 2
        end

        data_math["gen"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :component_type => "generator",
            :from => "generator.$name",
            :to => "gen.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_generator!,
            :kron_reduced => kron_reduced,
            :extra => Dict{String,Any}((k,v) for (k,v) in eng_obj if k in _extra_eng_data["generator"])
        )
    end
end


""
function _map_eng2math_pvsystem!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    for (name, eng_obj) in get(data_eng, "pvsystem", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("pvsystem", eng_obj, length(data_math["gen"])+1)

        phases = eng_obj["phases"]

        math_obj["name"] = name
        math_obj["gen_bus"] = data_math["bus_lookup"][eng_obj["bus"]]
        math_obj["gen_status"] = eng_obj["status"]

        math_obj["pg"] = fill(eng_obj["kva"] / 1e3 / phases, phases)
        math_obj["qg"] = fill(eng_obj["kva"] / 1e3 / phases, phases)
        math_obj["vg"] = fill(eng_obj["kv"] / data_math["basekv"], phases)

        math_obj["pmin"] = fill(0.0, phases)
        math_obj["pmax"] = fill(eng_obj["kva"] / 1e3 / phases, phases)

        math_obj["qmin"] =  fill(-eng_obj["kva"] / 1e3 / phases, phases)
        math_obj["qmax"] =  fill( eng_obj["kva"] / 1e3 / phases, phases)

        math_obj["conn"] = eng_obj["configuration"]

        data_math["gen"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => "pvsystem.$name",
            :to => "gen.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_pvsystem!,
            :kron_reduced => kron_reduced,
            :extra => Dict{String,Any}((k,v) for (k,v) in eng_obj if k in _extra_eng_data["pvsystem"])
        )
    end
end


""
function _map_eng2math_storage!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    for (name, eng_obj) in get(data_eng, "storage", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("storage", eng_obj, length(data_math["storage"])+1)

        phases = eng_obj["phases"]

        math_obj["name"] = name
        math_obj["storage_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        math_obj["energy"] = eng_obj["kwhstored"] / 1e3
        math_obj["energy_rating"] = eng_obj["kwhrated"] / 1e3
        math_obj["charge_rating"] = eng_obj["%charge"] * eng_obj["kwrated"] / 1e3 / 100.0
        math_obj["discharge_rating"] = eng_obj["%discharge"] * eng_obj["kwrated"] / 1e3 / 100.0
        math_obj["charge_efficiency"] = eng_obj["%effcharge"] / 100.0
        math_obj["discharge_efficiency"] = eng_obj["%effdischarge"] / 100.0
        math_obj["thermal_rating"] = fill(eng_obj["kva"] / 1e3 / phases, phases)
        math_obj["qmin"] = fill(-eng_obj["kvar"] / 1e3 / phases, phases)
        math_obj["qmax"] = fill( eng_obj["kvar"] / 1e3 / phases, phases)
        math_obj["r"] = fill(eng_obj["%r"] / 100.0, phases)
        math_obj["x"] = fill(eng_obj["%x"] / 100.0, phases)
        math_obj["p_loss"] = eng_obj["%idlingkw"] * eng_obj["kwrated"] / 1e3
        math_obj["q_loss"] = eng_obj["%idlingkvar"] * eng_obj["kvar"] / 1e3

        math_obj["ps"] = fill(0.0, phases)
        math_obj["qs"] = fill(0.0, phases)

        data_math["storage"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => "storage.$name",
            :to => "storage.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_storage!,
            :kron_reduced => kron_reduced,
            :extra => Dict{String,Any}((k,v) for (k,v) in eng_obj if k in _extra_eng_data["storage"])
        )
    end
end


""
function _map_eng2math_vsource!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    for (name, eng_obj) in get(data_eng, "voltage_source", Dict{Any,Dict{String,Any}}())
        if eng_obj["bus"] != data_eng["sourcebus"]
        end
    end
end


""
function _map_eng2math_line!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    for (name, eng_obj) in get(data_eng, "line", Dict{Any,Dict{String,Any}}())
        if haskey(eng_obj, "linecode")
            linecode = data_eng["linecode"][eng_obj["linecode"]]

            for property in ["rs", "xs", "g_fr", "g_to", "b_fr", "b_to"]
                if !haskey(eng_obj, property) && haskey(linecode, property)
                    eng_obj[property] = linecode[property]
                end
            end
        end

        math_obj = _init_math_obj("line", eng_obj, length(data_math["branch"])+1)

        nphases = length(eng_obj["f_connections"])
        nconductors = data_math["conductors"]

        Zbase = (data_math["basekv"] / sqrt(3))^2 * nconductors / (data_math["baseMVA"])
        Zbase = Zbase / 3

        math_obj["f_bus"] = data_math["bus_lookup"][eng_obj["f_bus"]]
        math_obj["t_bus"] = data_math["bus_lookup"][eng_obj["t_bus"]]

        math_obj["br_r"] = eng_obj["rs"] * eng_obj["length"] / Zbase
        math_obj["br_x"] = eng_obj["xs"] * eng_obj["length"] / Zbase

        math_obj["g_fr"] = Zbase * (2.0 * pi * data_eng["settings"]["basefreq"] * eng_obj["g_fr"] * eng_obj["length"] / 1e9)
        math_obj["g_to"] = Zbase * (2.0 * pi * data_eng["settings"]["basefreq"] * eng_obj["g_to"] * eng_obj["length"] / 1e9)

        math_obj["b_fr"] = Zbase * (2.0 * pi * data_eng["settings"]["basefreq"] * eng_obj["b_fr"] * eng_obj["length"] / 1e9)
        math_obj["b_to"] = Zbase * (2.0 * pi * data_eng["settings"]["basefreq"] * eng_obj["b_to"] * eng_obj["length"] / 1e9)

        math_obj["angmin"] = fill(-60.0, nphases)
        math_obj["angmax"] = fill( 60.0, nphases)

        math_obj["transformer"] = false
        math_obj["shift"] = zeros(nphases)
        math_obj["tap"] = ones(nphases)

        f_bus = data_eng["bus"][eng_obj["f_bus"]]
        t_bus = data_eng["bus"][eng_obj["t_bus"]]
        neutral = haskey(f_bus, "neutral") ? f_bus["neutral"] : haskey(t_bus, "neutral") ? t_bus["neutral"] : nphases+1

        _pad_properties!(math_obj, ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to", "angmin", "angmax", "tap", "shift"], eng_obj["f_connections"], collect(1:nconductors); neutral=neutral)

        math_obj["switch"] = false

        math_obj["br_status"] = eng_obj["status"]

        data_math["branch"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :from => "line.$name",
            :to => "branch.$(math_obj["index"])",
            :unmap_function => :_map_math2eng_load!,
            :kron_reduced => kron_reduced,
            :extra => Dict{String,Any}((k,v) for (k,v) in eng_obj if k in _extra_eng_data["line"])
        )
    end
end


""
function _map_eng2math_switch!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    # TODO enable real switches (right now only using vitual lines)
    for (name, eng_obj) in get(data_eng, "switch", Dict{Any,Dict{String,Any}}())
        nphases = length(eng_obj["f_connections"])
        nconductors = data_math["conductors"]

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

        branch_obj = _init_math_obj("switch", eng_obj, length(data_math["branch"])+1)

        Zbase = (data_math["basekv"] / sqrt(3))^2 * nconductors / (data_math["baseMVA"])
        Zbase = Zbase / 3

        _branch_obj = Dict{String,Any}(
            "name" => "_virtual_branch.switch.$name",
            "source_id" => "_virtual_branch.switch.$name",
            # "f_bus" => bus_obj["bus_i"],  # TODO enable real switches
            "f_bus" => data_math["bus_lookup"][eng_obj["f_bus"]],
            "t_bus" => data_math["bus_lookup"][eng_obj["t_bus"]],
            "br_r" => eng_obj["rs"] * eng_obj["length"] / Zbase,
            "br_x" => eng_obj["xs"] * eng_obj["length"] / Zbase,
            "g_fr" => Zbase * (2.0 * pi * data_eng["settings"]["basefreq"] * eng_obj["g_fr"] * eng_obj["length"] / 1e9),
            "g_to" => Zbase * (2.0 * pi * data_eng["settings"]["basefreq"] * eng_obj["g_to"] * eng_obj["length"] / 1e9),
            "b_fr" => Zbase * (2.0 * pi * data_eng["settings"]["basefreq"] * eng_obj["b_fr"] * eng_obj["length"] / 1e9),
            "b_to" => Zbase * (2.0 * pi * data_eng["settings"]["basefreq"] * eng_obj["b_to"] * eng_obj["length"] / 1e9),
            "angmin" => fill(-60.0, nphases),
            "angmax" => fill( 60.0, nphases),
            "transformer" => false,
            "shift" => zeros(nphases),
            "tap" => ones(nphases),
            "switch" => false,
            "br_status" => 1,
        )

        merge!(branch_obj, _branch_obj)

        f_bus = data_eng["bus"][eng_obj["f_bus"]]
        t_bus = data_eng["bus"][eng_obj["t_bus"]]
        neutral = haskey(f_bus, "neutral") ? f_bus["neutral"] : haskey(t_bus, "neutral") ? t_bus["neutral"] : nphases+1

        _pad_properties!(branch_obj, ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to", "angmin", "angmax", "tap", "shift"], eng_obj["f_connections"], collect(1:nconductors); neutral=neutral)

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
            :from_id => "switch.$name",
            # :to_id => ["switch.$(switch_obj["index"])", "bus.$(bus_obj["index"])", "branch.$(branch_obj["index"])"],  # TODO enable real switches
            :to_id => ["branch.$(branch_obj["index"])"],
            :unmap_function => :_map_math2eng_switch!,
            :kron_reduced => kron_reduced,
            :extra => Dict{String,Any}((k,v) for (k,v) in eng_obj if k in _extra_eng_data["switch"])
        )
    end
end


""
function _map_eng2math_transformer!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    for (name, eng_obj) in get(data_eng, "transformer", Dict{Any,Dict{String,Any}}())
        # Build map first, so we can update it as we decompose the transformer
        map_idx = length(data_math["map"])+1
        data_math["map"][map_idx] = Dict{Symbol,Any}(
            :from_id => "transformer.$name",
            :to_id => Vector{String}([]),
            :unmap_function => :_map_math2eng_transformer!,
            :kron_reduced => kron_reduced,
            :extra => Dict{String,Any}((k,v) for (k,v) in eng_obj if k in _extra_eng_data["transformer"])
        )

        to_map = data_math["map"][map_idx][:to_id]

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

        transformer_t_bus_w = _build_loss_model!(data_math, name, to_map, r_s, z_sc, y_sh)

        for w in 1:nrw
            # 2-WINDING TRANSFORMER
            # make virtual bus and mark it for reduction
            tm_nom = eng_obj["configuration"][w]=="delta" ? eng_obj["vnom"][w]*sqrt(3) : eng_obj["vnom"][w]

            transformer_2wa_obj = Dict{String,Any}(
                "name"          => "_virtual_transformer.$name.$w",
                "f_bus"         => data_math["bus_lookup"][eng_obj["bus"][w]],
                "t_bus"         => transformer_t_bus_w[w],
                "tm_nom"        => tm_nom,
                "f_connections" => eng_obj["connections"][w],
                "t_connections" => collect(1:4),
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

            data_math["transformer"]["$(transformer_2wa_obj["index"])"] = transformer_2wa_obj

            push!(to_map, "transformer.$(transformer_2wa_obj["index"])")
        end
    end
end


""
function _map_eng2math_sourcebus!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    # TODO create option for lossy vs lossless sourcebus connection
    for (name, eng_obj) in data_eng["voltage_source"]
        if eng_obj["bus"] == data_eng["sourcebus"]
            nconductors = data_math["conductors"]

            # TODO fix per unit problem
            bus_obj = Dict{String,Any}(
                "bus_i" => length(data_math["bus"])+1,
                "index" => length(data_math["bus"])+1,
                "name" => "_virtual_sourcebus",
                "bus_type" => 3,
                "vm" => eng_obj["vm"] ./ data_eng["settings"]["set_vbase_val"],
                "va" => eng_obj["va"],
                "vmin" => eng_obj["vm"] ./ data_eng["settings"]["set_vbase_val"],
                "vmax" => eng_obj["vm"] ./ data_eng["settings"]["set_vbase_val"],
                "basekv" => data_math["basekv"]
            )

            data_math["bus"]["$(bus_obj["index"])"] = bus_obj

            gen_obj = Dict{String,Any}(
                "gen_bus" => bus_obj["bus_i"],
                "name" => "_virtual_sourcebus",
                "gen_status" => 1,
                "pg" => fill(0.0, nconductors),
                "qg" => fill(0.0, nconductors),
                "model" => 2,
                "startup" => 0.0,
                "shutdown" => 0.0,
                "ncost" => 3,
                "cost" => [0.0, 1.0, 0.0],
                "conn" => "wye",  # TODO change name to configuration
                "index" => length(data_math["gen"]) + 1,
                "source_id" => "vsource._virtual_sourcebus"
            )

            data_math["gen"]["$(gen_obj["index"])"] = gen_obj

            vbase = data_math["basekv"]
            sbase = data_math["baseMVA"]
            zbase = vbase^2 / sbase / 3

            branch_obj = Dict{String,Any}(
                "name" => "_virtual_sourcebus",
                "source_id" => "vsource._virtual_sourcebus",
                "f_bus" => bus_obj["bus_i"],
                "t_bus" => data_math["bus_lookup"][data_eng["sourcebus"]],
                "angmin" => fill(-60.0, nconductors),
                "angmax" => fill( 60.0, nconductors),
                "shift" => fill(0.0, nconductors),
                "tap" => fill(1.0, nconductors),
                "tranformer" => false,
                "switch" => false,
                "br_status" => 1,
                "br_r" => eng_obj["rs"]./zbase,
                "br_x" => eng_obj["xs"]./zbase,
                "g_fr" => zeros(nconductors, nconductors),
                "g_to" => zeros(nconductors, nconductors),
                "b_fr" => zeros(nconductors, nconductors),
                "b_to" => zeros(nconductors, nconductors),
                "index" => length(data_math["branch"])+1
            )

            data_math["branch"]["$(branch_obj["index"])"] = branch_obj

            data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
                :from_id => "voltage_source.$name",
                :to_id => ["gen.$(gen_obj["index"])", "bus.$(bus_obj["index"])", "branch.$(branch_obj["index"])"],
                :unmap_function => :_map_math2eng_sourcebus!,
                :kron_reduced => kron_reduced,
                :extra => Dict{String,Any}((k,v) for (k,v) in eng_obj if k in _extra_eng_data["voltage_source"])
            )
        end
    end
end


""
function _map_math2eng!(data_math)
    @assert get(data_math, "data_model", "mathematical") == "mathematical" "Cannot map data to engineering model: provided data is not a mathematical model"
    @assert haskey(data_math, "map") "Cannot map data to engineering model: no mapping from mathematical to engineering data model is provided"

    data_eng = Dict{<:Any,<:Any}()

    map_keys = sort(keys(data_math["map"]); reverse=true)
    for map in map_keys
        # TODO
    end

end


""
function _expand_linecode!(data_model)
    # expand line codes
    for (id, line) in data_model["line"]
        if haskey(line, "linecode")
            linecode = data_model["linecode"][line["linecode"]]
            for key in ["rs", "xs", "g_fr", "g_to", "b_fr", "b_to"]
                line[key] = line["length"]*linecode[key]
            end
            delete!(line, "linecode")
            delete!(line, "length")
        end
    end
    delete!(data_model, "linecode")
end


""
function _lossy_ground_to_shunt!(data_model)
    mappings = []

    if haskey(data_model, "bus")
        for (id, bus) in data_model["bus"]
            grounding_lossy_inds = [i for (i,t) in enumerate(bus["grounded"]) if bus["rg"][i]!=0 || bus["xg"][i]!=0]
            grounding_lossy = bus["grounded"][grounding_lossy_inds]
            grounding_perfect = bus["grounded"][setdiff(1:length(bus["grounded"]), grounding_lossy_inds)]

            if !isempty(grounding_lossy)
                zg = bus["rg"][grounding_lossy_inds].+im*bus["xg"][grounding_lossy_inds]
                Y_sh = diagm(0=>inv.(zg)) # diagonal matrix, so matrix inverse is element-wise inverse
                add_virtual!(data_model, "shunt", create_shunt(bus=bus["id"], connections=grounding_lossy,
                    g_sh=real.(Y_sh), b_sh=imag.(Y_sh)
                ))
            end
        end
    end
    return mappings
end


""
function _load_to_shunt!(data_model)
    mappings = []
    if haskey(data_model, "load")
        for (id, load) in data_model["load"]
            if load["model"]=="constant_impedance"
                b = load["qd_ref"]./load["vnom"].^2*1E3
                g = load["pd_ref"]./load["vnom"].^2*1E3
                y = b.+im*g
                N = length(b)

                if load["configuration"]=="delta"
                    # create delta transformation matrix Md
                    Md = diagm(0=>ones(N), 1=>-ones(N-1))
                    Md[N,1] = -1
                    Y = Md'*diagm(0=>y)*Md

                else # load["configuration"]=="wye"
                    Y_fr = diagm(0=>y)
                    # B = [[b]; -1'*[b]]*[I -1]
                    Y = vcat(Y_fr, -ones(N)'*Y_fr)*hcat(diagm(0=>ones(N)),  -ones(N))
                end

                shunt = add_virtual!(data_model, "shunt", create_shunt(bus=load["bus"], connections=load["connections"], b_sh=imag.(Y), g_sh=real.(Y)))

                delete_component!(data_model, "load", load)

                push!(mappings, Dict(
                    "load" => load,
                    "shunt_id" => shunt["id"],
                ))
            end
        end
    end

    return mappings
end


""
function _capacitor_to_shunt!(data_model)
    mappings = []

    if haskey(data_model, "capacitor")
        for (id, cap) in data_model["capacitor"]
            b = cap["qd_ref"]./cap["vnom"]^2*1E3
            N = length(b)

            if cap["configuration"]=="delta"
                # create delta transformation matrix Md
                Md = diagm(0=>ones(N), 1=>-ones(N-1))
                Md[N,1] = -1
                B = Md'*diagm(0=>b)*Md

            elseif cap["configuration"]=="wye-grounded"
                B = diagm(0=>b)

            elseif cap["configuration"]=="wye-floating"
                # this is a floating wye-segment
                # B = [b]*(I-1/(b'*1)*[b';...;b'])
                B = diagm(0=>b)*(diagm(0=>ones(N)) - 1/sum(b)*repeat(b',N,1))

            else # cap["configuration"]=="wye"
                B_fr = diagm(0=>b)
                # B = [[b]; -1'*[b]]*[I -1]
                B = vcat(B_fr, -ones(N)'*B_fr)*hcat(diagm(0=>ones(N)),  -ones(N))
            end

            shunt = create_shunt(NaN, cap["bus"], cap["connections"], b_sh=B)
            add_virtual!(data_model, "shunt", shunt)
            delete_component!(data_model, "capacitor", cap)

            push!(mappings, Dict(
                "capacitor" => cap,
                "shunt_id" => shunt["id"],
            ))
        end
    end

    return mappings
end



"""
Converts a set of short-circuit tests to an equivalent reactance network.
Reference:
R. C. Dugan, “A perspective on transformer modeling for distribution system analysis,”
in 2003 IEEE Power Engineering Society General Meeting (IEEE Cat. No.03CH37491), 2003, vol. 1, pp. 114-119 Vol. 1.
"""
function _sc2br_impedance(Zsc)
    N = maximum([maximum(k) for k in keys(Zsc)])
    # check whether no keys are missing
    # Zsc should contain tupples for upper triangle of NxN
    for i in 1:N
        for j in i+1:N
            if !haskey(Zsc, (i,j))
                if haskey(Zsc, (j,i))
                    # Zsc is symmetric; use value of lower triangle if defined
                    Zsc[(i,j)] =  Zsc[(j,i)]
                else
                    Memento.error(_LOGGER, "Short-circuit impedance between winding $i and $j is missing.")
                end
            end
        end
    end
    # make Zb
    Zb = zeros(Complex{Float64}, N-1,N-1)
    for i in 1:N-1
        Zb[i,i] = Zsc[(1,i+1)]
    end
    for i in 1:N-1
        for j in 1:i-1
            Zb[i,j] = (Zb[i,i]+Zb[j,j]-Zsc[(j+1,i+1)])/2
            Zb[j,i] = Zb[i,j]
        end
    end
    # get Ybus
    Y = LinearAlgebra.pinv(Zb)
    Y = [-Y*ones(N-1) Y]
    Y = [-ones(1,N-1)*Y; Y]
    # extract elements
    Zbr = Dict()
    for k in keys(Zsc)
        Zbr[k] = (abs(Y[k...])==0) ? Inf : -1/Y[k...]
    end
    return Zbr
end


""
function _build_loss_model!(data_math::Dict{String,<:Any}, transformer_name::String, to_map::Vector{String}, r_s::Vector{Float64}, zsc::Dict{Tuple{Int,Int},Complex{Float64}}, ysh::Complex{Float64}; nphases::Int=3)
    # precompute the minimal set of buses and lines
    N = length(r_s)
    tr_t_bus = collect(1:N)
    buses = Set(1:2*N)

    edges = [[[i,i+N] for i in 1:N]..., [[i+N,j+N] for (i,j) in keys(zsc)]...]
    lines = Dict(enumerate(edges))

    z = Dict(enumerate([r_s..., values(zsc)...]))

    shunts = Dict(2=>ysh)

    # remove Inf lines
    for (l,edge) in lines
        if real(z[l])==Inf || imag(z[l])==Inf
            delete!(lines, l)
            delete!(z, l)
        end
    end

    # merge short circuits
    stack = Set(keys(lines))

    while !isempty(stack)
        l = pop!(stack)
        if z[l] == 0
            (i,j) = lines[l]

            # remove line
            delete!(lines, l)

            # remove  bus j
            delete!(buses, j)

            # update lines
            for (k,(edge)) in lines
                if edge[1] == j
                    edge[1] = i
                end
                if edge[2] == j
                    edge[2] = i
                end
                if edge[1]==edge[2]
                    delete!(lines, k)
                    delete!(stack, k)
                end
            end

            # move shunts
            if haskey(shunts, j)
                if haskey(shunts, i)
                    shunts[i] += shunts[j]
                else
                    shunts[i] = shunts[j]
                end
            end

            # update transformer buses
            for w in 1:N
                if tr_t_bus[w] == j
                    tr_t_bus[w] = i
                end
            end
        end
    end

    bus_ids = Dict()
    for bus in buses
        bus_obj = Dict{String,Any}(
            "name" => "_virtual_bus.transformer.$(transformer_name)_$(bus)",
            "bus_i" => length(data_math["bus"])+1,
            "vm" => fill(1.0, nphases),
            "va" => fill(0.0, nphases),
            "vmin" => fill(0.0, nphases),
            "vmax" => fill(Inf, nphases),
            "base_kv" => 1.0,
            "bus_type" => 1,
            "status" => 1,
            "index" => length(data_math["bus"])+1,
        )

        data_math["bus"]["$(bus_obj["index"])"] = bus_obj

        bus_ids[bus] = bus_obj["bus_i"]

        push!(to_map, "bus.$(bus_obj["index"])")
    end


    for (l,(i,j)) in lines
        # merge the shunts into the shunts of the pi model of the line
        g_fr = b_fr = g_to = b_to = 0

        if haskey(shunts, i)
            g_fr = real(shunts[i])
            b_fr = imag(shunts[i])
            delete!(shunts, i)
        end

        if haskey(shunts, j)
            g_fr = real(shunts[j])
            b_fr = imag(shunts[j])
            delete!(shunts, j)
        end

        branch_obj = Dict{String,Any}(
            "name" => "_virtual_branch.transformer.$(transformer_name)_$(l)",
            "source_id" => "_virtual_branch.transformer.$(transformer_name)_$(l)",
            "index" => length(data_math["branch"])+1,
            "br_status"=>1,
            "f_bus"=>bus_ids[i],
            "t_bus"=>bus_ids[j],
            "f_connections"=>collect(1:nphases),
            "t_connections"=>collect(1:nphases),
            "br_r" => diagm(0=>fill(real(z[l]), nphases)),
            "br_x" => diagm(0=>fill(imag(z[l]), nphases)),
            "g_fr" => diagm(0=>fill(g_fr, nphases)),
            "b_fr" => diagm(0=>fill(b_fr, nphases)),
            "g_to" => diagm(0=>fill(g_to, nphases)),
            "b_to" => diagm(0=>fill(b_to, nphases)),
            "angmin" => fill(-60.0, nphases),
            "angmax" => fill( 60.0, nphases),
            "shift" => zeros(nphases),
            "tap" => ones(nphases)
        )

        data_math["branch"]["$(branch_obj["index"])"] = branch_obj

        push!(to_map, "branch.$(branch_obj["index"])")
    end

    return [bus_ids[bus] for bus in tr_t_bus]
end


""
function _pad_properties!(object::Dict{<:Any,<:Any}, properties::Vector{String}, connections::Vector{Int}, phases::Vector{Int}; neutral::Int=4, kron_reduced::Bool=true)
    if kron_reduced
        pos = Dict((x,i) for (i,x) in enumerate(phases))
        inds = [pos[x] for x in connections[connections.!=neutral]]
    else
        # TODO
    end

    for property in properties
        if haskey(object, property)
            if isa(object[property], Vector)
                tmp = zeros(length(phases))
                tmp[inds] = object[property]
                object[property] = tmp
            elseif isa(object[property], Matrix)
                tmp = zeros(length(phases), length(phases))
                tmp[inds, inds] = object[property]
                object[property] = tmp
            end
        end
    end
end


# MAP SOLUTION UP
""
function solution_unmap!(solution::Dict, data_model::Dict)
    for (name, data) in reverse(data_model["mappings"])
        if name=="decompose_transformer_nw"
            for bus_id in values(data["vbuses"])
                delete!(solution["bus"], bus_id)
            end

            for line_id in values(data["vlines"])
                delete!(solution["branch"], line_id)
            end

            pt = [solution["transformer"][tr_id]["pf"] for tr_id in data["trans_2wa"]]
            qt = [solution["transformer"][tr_id]["qf"] for tr_id in data["trans_2wa"]]
            for tr_id in data["trans_2wa"]
                delete!(solution["transformer"], tr_id)
            end

            add_solution!(solution, "transformer_nw", data["trans"]["id"], Dict("pt"=>pt, "qt"=>qt))
        elseif name=="capacitor_to_shunt"
            # shunt has no solutions defined
            delete_solution!(solution, "shunt", data["shunt_id"])
            add_solution!(solution, "capacitor", data["capacitor"]["id"], Dict())
        elseif name=="load_to_shunt"
            # shunt has no solutions, but a load should have!
            delete!(solution, "shunt", data["shunt_id"])
            add_solution!(solution, "load", data["load"]["id"], Dict())
        elseif name=="decompose_voltage_source"
            gen = solution["gen"][data["gen_id"]]
            delete_solution!(solution, "gen", data["gen_id"])
            add_solution!(solution, "voltage_source", data["voltage_source"]["id"], Dict("pg"=>gen["pg"], "qg"=>gen["qg"]))

        end
    end

    # remove component dicts if empty
    for (comp_type, comp_dict) in solution
        if isa(comp_dict, Dict) && isempty(comp_dict)
            delete!(solution, comp_type)
        end
    end
end



"""
Given a set of addmittances 'y' connected from the conductors 'f_cnds' to the
conductors 't_cnds', this method will return a list of conductors 'cnd' and a
matrix 'Y', which will satisfy I[cnds] = Y*V[cnds].
"""
function calc_shunt(f_cnds::Vector{Int}, t_cnds::Vector{Int}, y::Vector{Float64})
    cnds = unique([f_cnds..., t_cnds...])
    e(f,t) = reshape([c==f ? 1 : c==t ? -1 : 0 for c in cnds], length(cnds), 1)
    Y = sum([e(f_cnds[i], t_cnds[i])*y[i]*e(f_cnds[i], t_cnds[i])' for i in 1:length(y)])
    return (cnds, Y)
end



"""
Given a set of terminals 'cnds' with associated shunt addmittance 'Y', this
method will calculate the reduced addmittance matrix if terminal 'ground' is
grounded.
"""
function _calc_ground_shunt_admittance_matrix(cnds::Vector{Int}, Y::Matrix{Float64}, ground::Int)
    # TODO add types
    if ground in cnds
        cndsr = setdiff(cnds, ground)
        cndsr_inds = _get_idxs(cnds, cndsr)
        Yr = Y[cndsr_inds, cndsr_inds]
        return (cndsr, Yr)
    else
        return cnds, Y
    end
end


""
function _rm_floating_cnd(cnds::Vector{Int}, Y::Matrix{Float64}, f::Int)
    P = setdiff(cnds, f)
    f_inds = _get_idxs(cnds, [f])
    P_inds = _get_idxs(cnds, P)
    Yrm = Y[P_inds,P_inds]-(1/Y[f_inds,f_inds][1])*Y[P_inds,f_inds]*Y[f_inds,P_inds]
    return (P,Yrm)
end
