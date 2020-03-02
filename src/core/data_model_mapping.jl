import LinearAlgebra


const _1to1_maps = Dict{String,Array{String,1}}(
    "bus" => ["vm", "va", "vmin", "vmax"],
    "load" => ["pd", "qd", "model", "configuration", "status"],
    "capacitor" => ["status"],
    "shunt_reactor" => ["status"],
    "generator" => ["configuration", "status"],
    "pvsystem" => ["status"],
    "storage" => ["status"],
    "line" => [],
    "switch" => [],
    "transformer" => [],
    "vsource" => [],
)

const _extra_eng_data = Dict{String,Array{String,1}}(
    "root" => ["sourcebus", "files", "dss_options", "settings"],
    "bus" => ["grounded", "neutral", "awaiting_ground", "xg", "phases", "rg", "terminals"],
    "load" => [],
    "line" => ["f_connections", "t_connections", "linecode"],
)


const _node_elements = ["load", "capacitor", "shunt_reactor", "generator", "pvsystem", "storage", "vsource"]

const _edge_elements = ["line", "switch", "transformer"]

# MAP DATA MODEL DOWN

function _map_eng2math(data_eng; kron_reduced::Bool=true)
    @assert get(data_eng, "data_model", "mathematical") == "engineering"

    data_model_make_pu!(data_eng)

    data_math = Dict{String,Any}(
        "name" => data_eng["name"],
        "per_unit" => get(data_eng, "per_unit", false)
    )

    data_math["map"] = Dict{Int,Dict{Symbol,Any}}(
        1 => Dict{Symbol,Any}(
            :component_type => "root",
            :unmap_function => :_map_math2eng_root!,
            :extra => Dict{String,Any}((k,v) for (k,v) in data_eng if k in _extra_eng_data["root"])
        )
    )

    data_math["settings"] = deepcopy(data_eng["settings"])

    data_math["lookup"] = Dict{String,Dict{Any,Int}}()

    data_math["conductors"] = kron_reduced ? 3 : 4
    data_math["basekv"] = data_eng["settings"]["set_vbase_val"]
    data_math["baseMVA"] = data_eng["settings"]["set_sbase_val"]

    _map_eng2math_bus!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_load!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_capacitor!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_shunt_reactor!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_generator!(data_math, data_eng; kron_reduced=kron_reduced)
    # _map_eng2math_pvsystem!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_storage!(data_math, data_eng; kron_reduced=kron_reduced)
    # _map_eng2math_vsource!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_line!(data_math, data_eng; kron_reduced=kron_reduced)
    _map_eng2math_switch!(data_math, data_eng; kron_reduced=kron_reduced)

    # _map_eng2math_transformer!(data_math, data_eng; kron_reduced=kron_reduced)

    _map_eng2math_sourcebus!(data_math, data_eng; kron_reduced=kron_reduced)

    # # needs to happen before _expand_linecode, as it might contain a linecode for the internal impedance
    # add_mappings!(data_eng, "decompose_voltage_source", _decompose_voltage_source!(data_eng))
    # _expand_linecode!(data_eng)
    # # creates shunt of 4x4; disabled for now (incompatible 3-wire kron-reduced)
    # #add_mappings!(data_eng, "load_to_shunt", _load_to_shunt!(data_eng))
    # add_mappings!(data_eng, "capacitor_to_shunt", _capacitor_to_shunt!(data_eng))
    # add_mappings!(data_eng, "decompose_transformer_nw", _decompose_transformer_nw!(data_eng))
    # add_mappings!(data_eng, "_lossy_ground_to_shunt", _lossy_ground_to_shunt!(data_eng))

    # # add low level component types if not present yet
    # for comp_type in ["load", "generator", "bus", "line", "shunt", "transformer", "storage", "switch"]
    #     if !haskey(data_eng, comp_type)
    #         data_eng[comp_type] = Dict{String, Any}()
    #     end
    # end

    data_math["dcline"] = Dict{String,Any}()

    # data_math["per_unit"] = false

    # data_model_make_pu!(data_math)

    delete!(data_math, "lookup")


    data_math["data_model"] = "mathematical"

    return data_math
end


""
function _map_defaults(eng_obj::Dict{String,Any}, component_type::String, component_name::Any, kron_reduced::Bool=true; phases::Vector{Int}=[1, 2, 3], neutral::Int=4)::Dict{String,Any}
    math_obj = Dict{String,Any}()

    math_obj["name"] = component_name

    if component_type in _node_elements
        math_obj["source_id"] = eng_obj["source_id"]
        terminals = eng_obj["connections"]
    elseif component_type in _edge_elements
        f_terminals = eng_obj["f_connections"]
        t_terminals = eng_obj["t_connections"]
    elseif component_type == "bus"
        terminals = eng_obj["terminals"]
    end

    # TODO clean this up
    for key in _1to1_maps[component_type]
        if haskey(eng_obj, key)
            if kron_reduced
                if component_type == "bus"
                    terminals = eng_obj["terminals"]
                    math_obj[key] = eng_obj[key][terminals.!=neutral]
                elseif component_type in _node_elements
                    math_obj[key] = eng_obj[key]
                    _pad_properties!(math_obj, [key], eng_obj["connections"], phases)
                elseif component_type in _edge_elements
                    math_obj[key] = eng_obj[key]
                    # TODO
                end
            else
                math_obj[key] = eng_obj[key]
            end
        end
    end

    return math_obj
end


""
function _map_eng2math_bus!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    if !haskey(data_math, "bus")
        data_math["bus"] = Dict{String,Any}()
    end

    for (name, eng_obj) in get(data_eng, "bus", Dict{String,Any}())

        phases = get(eng_obj, "phases", [1, 2, 3])
        neutral = get(eng_obj, "neutral", 4)
        terminals = eng_obj["terminals"]

        @assert all(t in [phases..., neutral] for t in terminals)

        math_obj = _map_defaults(eng_obj, "bus", name, kron_reduced; phases=phases, neutral=neutral)

        math_obj["vm"] = get(eng_obj, "vm", fill(1.0, length(phases)))
        math_obj["va"] = get(eng_obj, "va", [_wrap_to_180(-rad2deg(2*pi/length(phases)*(i-1))) for i in phases])

        math_obj["vmin"] = fill(NaN, length(phases))
        math_obj["vmax"] = fill(NaN, length(phases))

        math_obj["base_kv"] = eng_obj["vbase"]

        math_obj["bus_type"] = eng_obj["status"] == 1 ? 1 : 4

        math_obj["index"] = length(data_math["bus"]) + 1
        math_obj["bus_i"] = math_obj["index"]

        data_math["bus"][string(math_obj["index"])] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :component_type => "bus",
            :from => name,
            :to => "$(math_obj["index"])",
            :unmap_function => :_map_math2eng_bus!,
            :kron_reduced => kron_reduced,
            :extra => Dict{String,Any}((k,v) for (k,v) in eng_obj if k in _extra_eng_data["bus"])
        )

        if !haskey(data_math["lookup"], "bus")
            data_math["lookup"]["bus"] = Dict{Any,Int}()
        end

        data_math["lookup"]["bus"][name] = math_obj["index"]
    end
end


""
function _map_eng2math_load!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    if !haskey(data_math, "load")
        data_math["load"] = Dict{String,Any}()
    end

    for (name, eng_obj) in get(data_eng, "load", Dict{Any,Dict{String,Any}}())
        math_obj = _map_defaults(eng_obj, "load", name, kron_reduced)

        if eng_obj["configuration"] == "wye"
            bus = data_eng["bus"][eng_obj["bus"]]
            # TODO add message for failure
            @assert length(bus["grounded"]) == 1 && bus["grounded"][1] == eng_obj["connections"][end]
        else
            # TODO add message for failure
            @assert all(load["connections"] .== phases)
        end

        math_obj["load_bus"] = data_math["lookup"]["bus"][eng_obj["bus"]]

        math_obj["index"] = length(data_math["load"]) + 1

        data_math["load"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :component_type => "load",
            :from => name,
            :to => "$(math_obj["index"])",
            :unmap_function => :_map_math2eng_load!,
            :kron_reduced => kron_reduced,
            :extra => Dict{String,Any}((k,v) for (k,v) in eng_obj if k in _extra_eng_data["load"])
        )

        if !haskey(data_math["lookup"], "load")
            data_math["lookup"]["load"] = Dict{Any,Int}()
        end

        data_math["lookup"]["load"][name] = math_obj["index"]
    end
end


""
function _map_eng2math_capacitor!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    if !haskey(data_math, "shunt")
        data_math["shunt"] = Dict{String,Any}()
    end

    for (name, eng_obj) in get(data_eng, "capacitor", Dict{Any,Dict{String,Any}}())
        math_obj = _map_defaults(eng_obj, "load", name, kron_reduced)

        math_obj["shunt_bus"] = data_math["lookup"]["bus"][eng_obj["bus"]]

        math_obj["index"] = length(data_math["shunt"]) + 1

        data_math["shunt"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :component_type => "capacitor",
            :from => name,
            :to => "$(math_obj["index"])",
            :unmap_function => :_map_math2eng_load!,
            :kron_reduced => kron_reduced,
            :extra => Dict{String,Any}((k,v) for (k,v) in eng_obj if k in _extra_eng_data["load"])
        )

        if !haskey(data_math["lookup"], "capacitor")
            data_math["lookup"]["capacitor"] = Dict{Any,Int}()
        end

        data_math["lookup"]["capacitor"][name] = math_obj["index"]
    end
end


""
function _map_eng2math_shunt_reactor!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    if !haskey(data_math, "shunt")
        data_math["shunt"] = Dict{String,Any}()
    end

end


""
function _map_eng2math_generator!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    if !haskey(data_math, "gen")
        data_math["gen"] = Dict{String,Any}()
    end
end


""
function _map_eng2math_pvsystem!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    if !haskey(data_math, "gen")
        data_math["gen"] = Dict{String,Any}()
    end

end


""
function _map_eng2math_storage!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    if !haskey(data_math, "storage")
        data_math["storage"] = Dict{String,Any}()
    end
end


""
function _map_eng2math_vsource!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    if !haskey(data_math, "gen")
        data_math["gen"] = Dict{String,Any}()
    end


end


""
function _map_eng2math_line!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    if !haskey(data_math, "branch")
        data_math["branch"] = Dict{String,Any}()
    end

    for (name, eng_obj) in get(data_eng, "line", Dict{Any,Dict{String,Any}}())
        if haskey(eng_obj, "linecode")
            linecode = data_eng["linecode"][eng_obj["linecode"]]

            for property in ["rmatrix", "xmatrix", "cmatrix"]
                if !haskey(eng_obj, property) && haskey(linecode, property)
                    eng_obj[property] = linecode[property]
                end
            end
        end

        nphases = length(eng_obj["f_connections"])

        math_obj = _map_defaults(eng_obj, "line", name, kron_reduced)

        math_obj["f_bus"] = data_math["lookup"]["bus"][eng_obj["f_bus"]]
        math_obj["t_bus"] = data_math["lookup"]["bus"][eng_obj["t_bus"]]

        math_obj["br_r"] = eng_obj["rmatrix"] * eng_obj["length"]
        math_obj["br_x"] = eng_obj["xmatrix"] * eng_obj["length"]

        math_obj["g_fr"] = fill(0.0, nphases, nphases)
        math_obj["g_to"] = fill(0.0, nphases, nphases)

        math_obj["b_fr"] = (2.0 * pi * data_eng["settings"]["basefreq"] * eng_obj["cmatrix"] * eng_obj["length"] / 1e9) / 2.0
        math_obj["b_to"] = (2.0 * pi * data_eng["settings"]["basefreq"] * eng_obj["cmatrix"] * eng_obj["length"] / 1e9) / 2.0

        math_obj["angmin"] = fill(-60.0, nphases)
        math_obj["angmax"] = fill( 60.0, nphases)

        math_obj["transformer"] = false
        math_obj["shift"] = zeros(nphases)
        math_obj["tap"] = ones(nphases)

        math_obj["switch"] = false

        math_obj["br_status"] = eng_obj["status"]

        math_obj["index"] = length(data_math["branch"])+1

        data_math["branch"]["$(math_obj["index"])"] = math_obj

        data_math["map"][length(data_math["map"])+1] = Dict{Symbol,Any}(
            :component_type => "line",
            :from => name,
            :to => "$(math_obj["index"])",
            :unmap_function => :_map_math2eng_load!,
            :kron_reduced => kron_reduced,
            :extra => Dict{String,Any}((k,v) for (k,v) in eng_obj if k in _extra_eng_data["line"])
        )

        if !haskey(data_math["lookup"], "line")
            data_math["lookup"]["line"] = Dict{Any,Int}()
        end

        data_math["lookup"]["line"][name] = math_obj["index"]

    end
end


""
function _map_eng2math_switch!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    if !haskey(data_math, "switch")
        data_math["switch"] = Dict{String,Any}()
    end
end


""
function _map_eng2math_transformer!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)

end


""
function _map_eng2math_sourcebus!(data_math::Dict{String,<:Any}, data_eng::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    if !haskey(data_math, "bus")
        data_math["bus"] = Dict{String,Any}()
    end

    if !haskey(data_math, "gen")
        data_math["gen"] = Dict{String,Any}()
    end

    if !haskey(data_math, "branch")
        data_math["branch"] = Dict{String,Any}()
    end

    sourcebus = data_eng["sourcebus"]
    sourcebus_vsource = data_eng["voltage_source"][sourcebus]

    nconductors = data_math["conductors"]

    # TODO fix per unit problem
    bus_obj = Dict{String,Any}(
        "bus_i" => length(data_math["bus"])+1,
        "index" => length(data_math["bus"])+1,
        "name" => "_virtual_sourcebus",
        "bus_type" => 3,
        "vm" => sourcebus_vsource["vm"],
        "va" => sourcebus_vsource["va"],
        "vmin" => sourcebus_vsource["vm"],
        "vmax" => sourcebus_vsource["vm"],
        "basekv" => data_math["basekv"]
    )

    data_math["bus"]["$(bus_obj["index"])"] = bus_obj

    gen_obj = Dict{String,Any}(
        "gen_bus" => bus_obj["bus_i"],
        "name" => "_virtual_sourcebus",
        "gen_status" => sourcebus_vsource["status"],
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
        "t_bus" => data_math["lookup"]["bus"][sourcebus],
        "angmin" => fill(-60.0, nconductors),
        "angmax" => fill( 60.0, nconductors),
        "shift" => fill(0.0, nconductors),
        "tap" => fill(1.0, nconductors),
        "tranformer" => false,
        "switch" => false,
        "br_status" => 1,
        "br_r" => sourcebus_vsource["rmatrix"]./zbase,
        "br_x" => sourcebus_vsource["xmatrix"]./zbase,
        "g_fr" => zeros(nconductors, nconductors),
        "g_to" => zeros(nconductors, nconductors),
        "b_fr" => zeros(nconductors, nconductors),
        "b_to" => zeros(nconductors, nconductors),
        "index" => length(data_math["branch"])+1
    )
    # branch_obj = _create_vbranch(data_math, data_math["lookup"]["bus"][sourcebus], bus_obj["bus_i"]; name="_virtual_sourcebus", br_r=sourcebus_vsource["rs"]/1e3, br_x=sourcebus_vsource["xs"]/1e3)

    data_math["branch"]["$(branch_obj["index"])"] = branch_obj
end


"""
This function adds a new branch to the data model and returns its dictionary.
It is virtual in the sense that it does not correspond to a branch in the
network, but is part of the decomposition of the transformer.
"""
function _create_vbranch(data_math::Dict{<:Any,<:Any}, f_bus::Int, t_bus::Int; name::String="", source_id::String="", active_phases::Vector{Int}=[1, 2, 3], kwargs...)
    ncnd = data_math["conductors"]

    kwargs = Dict{Symbol,Any}(kwargs)

    vbase = haskey(kwargs, :vbase) ? kwargs[:vbase] : data_math["basekv"]
    # TODO assumes per_unit will be flagged
    sbase = haskey(kwargs, :sbase) ? kwargs[:sbase] : data_math["baseMVA"]
    zbase = vbase^2/sbase
    # convert to LN vbase in instead of LL vbase
    zbase *= (1/3)

    vbranch = Dict{String, Any}("f_bus"=>f_bus, "t_bus"=>t_bus, "name"=>name)

    vbranch["active_phases"] = active_phases
    vbranch["source_id"] = "virtual_branch.$name"

    for k in [:br_r, :br_x, :g_fr, :g_to, :b_fr, :b_to]
        if !haskey(kwargs, k)
            vbranch[string(k)] = zeros(ncnd, ncnd)
        else
            if k in [:br_r, :br_x]
                vbranch[string(k)] = kwargs[k]./zbase
            else
                vbranch[string(k)] = kwargs[k].*zbase
            end
        end
    end

    vbranch["angmin"] = -ones(ncnd)*60
    vbranch["angmax"] = ones(ncnd)*60

    vbranch["rate_a"] = get(kwargs, :rate_a, fill(Inf, length(active_phases)))

    vbranch["shift"] = zeros(ncnd)
    vbranch["tap"] = ones(ncnd)

    vbranch["transformer"] = false
    vbranch["switch"] = false
    vbranch["br_status"] = 1

    for k in [:rate_a, :rate_b, :rate_c, :c_rating_a, :c_rating_b, :c_rating_c]
        if haskey(kwargs, k)
            vbranch[string(k)] = kwargs[k]
        end
    end

    vbranch["index"] = length(data_math["branch"])+1

    return vbranch
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
                Y_sh = LinearAlgebra.diagm(0=>inv.(zg)) # diagonal matrix, so matrix inverse is element-wise inverse
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
                    Md = LinearAlgebra.diagm(0=>ones(N), 1=>-ones(N-1))
                    Md[N,1] = -1
                    Y = Md'*LinearAlgebra.diagm(0=>y)*Md

                else # load["configuration"]=="wye"
                    Y_fr = LinearAlgebra.diagm(0=>y)
                    # B = [[b]; -1'*[b]]*[I -1]
                    Y = vcat(Y_fr, -ones(N)'*Y_fr)*hcat(LinearAlgebra.diagm(0=>ones(N)),  -ones(N))
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
                Md = LinearAlgebra.diagm(0=>ones(N), 1=>-ones(N-1))
                Md[N,1] = -1
                B = Md'*LinearAlgebra.diagm(0=>b)*Md

            elseif cap["configuration"]=="wye-grounded"
                B = LinearAlgebra.diagm(0=>b)

            elseif cap["configuration"]=="wye-floating"
                # this is a floating wye-segment
                # B = [b]*(I-1/(b'*1)*[b';...;b'])
                B = LinearAlgebra.diagm(0=>b)*(LinearAlgebra.diagm(0=>ones(N)) - 1/sum(b)*repeat(b',N,1))

            else # cap["configuration"]=="wye"
                B_fr = LinearAlgebra.diagm(0=>b)
                # B = [[b]; -1'*[b]]*[I -1]
                B = vcat(B_fr, -ones(N)'*B_fr)*hcat(LinearAlgebra.diagm(0=>ones(N)),  -ones(N))
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


""
function _decompose_voltage_source!(data_model)
    mappings = []

    if haskey(data_model, "voltage_source")
        for (id, vs) in data_model["voltage_source"]

            bus = data_model["bus"][vs["bus"]]

            line_kwargs = Dict(Symbol(prop)=>vs[prop] for prop in ["rs", "xs", "g_fr", "b_fr", "g_to", "b_to", "linecode", "length"] if haskey(vs, prop))
            lossy = !isempty(line_kwargs)

            # if any loss parameters (or linecode) were supplied, then create a line and internal bus
            if lossy
                sourcebus = add_virtual!(data_model, "bus", create_bus(terminals=deepcopy(vs["connections"])))

                line = add_virtual!(data_model, "line", create_line(;
                    f_bus=sourcebus["id"], f_connections=vs["connections"], t_bus=bus["id"], t_connections=vs["connections"],
                    line_kwargs...
                ))
            else
                sourcebus = bus
            end

            ground = _get_ground!(sourcebus)
            gen = create_generator(bus=sourcebus["id"], connections=[vs["connections"]..., ground])

            for prop in ["pg_max", "pg_min", "qg_max", "qg_min"]
                if haskey(vs, prop)
                    gen[prop] = vs[prop]
                end
            end

            add_virtual!(data_model, "generator", gen)

            conns = vs["connections"]
            terminals = bus["terminals"]

            tmp = Dict(enumerate(conns))
            sourcebus["vm"] = sourcebus["vmax"] = sourcebus["vmin"] = [haskey(tmp, t) ? vs["vm"][tmp[t]] : NaN for t in terminals]
            sourcebus["va"] = [haskey(tmp, t) ? vs["va"][tmp[t]] : NaN for t in terminals]
            sourcebus["bus_type"] = 3

            delete_component!(data_model, "voltage_source", vs["id"])
            push!(mappings, Dict("voltage_source"=>vs, "gen_id"=>gen["id"],
                "vbus_id"  => lossy ? sourcebus["id"] : nothing,
                "vline_id" => lossy ? line["id"]      : nothing,
            ))
        end
    end

    return mappings
end


"""
Replaces complex transformers with a composition of ideal transformers and lines
which model losses. New buses (virtual, no physical meaning) are added.
"""
function _decompose_transformer_nw!(data_model)
    mappings = []

    if haskey(data_model, "transformer_nw")
        for (tr_id, trans) in data_model["transformer_nw"]

            vnom = trans["vnom"]*data_model["settings"]["kv_kvar_scalar"]
            snom = trans["snom"]*data_model["settings"]["kv_kvar_scalar"]

            nrw = length(trans["bus"])

            # calculate zbase in which the data is specified, and convert to SI
            zbase = (vnom.^2)./snom
            # x_sc is specified with respect to first winding
            x_sc = trans["xsc"].*zbase[1]
            # rs is specified with respect to each winding
            r_s = trans["rs"].*zbase

            g_sh =  (trans["noloadloss"]*snom[1])/vnom[1]^2
            b_sh = -(trans["imag"]*snom[1])/vnom[1]^2

            # data is measured externally, but we now refer it to the internal side
            ratios = vnom/data_model["settings"]["kv_kvar_scalar"]
            x_sc = x_sc./ratios[1]^2
            r_s = r_s./ratios.^2
            g_sh = g_sh*ratios[1]^2
            b_sh = b_sh*ratios[1]^2

            # convert x_sc from list of upper triangle elements to an explicit dict
            y_sh = g_sh + im*b_sh
            z_sc = Dict([(key, im*x_sc[i]) for (i,key) in enumerate([(i,j) for i in 1:nrw for j in i+1:nrw])])

            vbuses, vlines, trans_t_bus_w = _build_loss_model!(data_model, r_s, z_sc, y_sh)

            trans_ids_w = Array{String, 1}(undef, nrw)
            for w in 1:nrw
                # 2-WINDING TRANSFORMER
                # make virtual bus and mark it for reduction
                tm_nom = trans["configuration"][w]=="delta" ? trans["vnom"][w]*sqrt(3) : trans["vnom"][w]
                trans_2wa = add_virtual!(data_model, "transformer", Dict(
                    "f_bus"         => trans["bus"][w],
                    "t_bus"         => trans_t_bus_w[w],
                    "tm_nom"        => tm_nom,
                    "f_connections" => trans["connections"][w],
                    "t_connections" => collect(1:4),
                    "configuration" => trans["configuration"][w],
                    "polarity"      => trans["polarity"][w],
                    "tm"            => trans["tm"][w],
                    "fixed"         => trans["fixed"][w],
                ))

                for prop in ["tm_min", "tm_max", "tm_step"]
                    if haskey(trans, prop)
                        trans_2wa[prop] = trans[prop][w]
                    end
                end

                trans_ids_w[w] = trans_2wa["id"]
            end

            delete_component!(data_model, "transformer_nw", trans)

            push!(mappings, Dict(
                "trans"=>trans,
                "trans_2wa"=>trans_ids_w,
                "vlines"=>vlines,
                "vbuses"=>vbuses,
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
function _build_loss_model!(data_model, r_s, zsc, ysh; n_phases=3)
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
                if edge[1]==j
                    edge[1] = i
                end
                if edge[2]==j
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
                if tr_t_bus[w]==j
                    tr_t_bus[w] = i
                end
            end
        end
    end

    bus_ids = Dict()
    for bus in buses
        bus_ids[bus] = add_virtual_get_id!(data_model, "bus", create_bus(id=""))
    end
    line_ids = Dict()
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
        line_ids[l] = add_virtual_get_id!(data_model, "line", Dict(
            "status"=>1,
            "f_bus"=>bus_ids[i], "t_bus"=>bus_ids[j],
            "f_connections"=>collect(1:n_phases),
            "t_connections"=>collect(1:n_phases),
            "rs"=>LinearAlgebra.diagm(0=>fill(real(z[l]), n_phases)),
            "xs"=>LinearAlgebra.diagm(0=>fill(imag(z[l]), n_phases)),
            "g_fr"=>LinearAlgebra.diagm(0=>fill(g_fr, n_phases)),
            "b_fr"=>LinearAlgebra.diagm(0=>fill(b_fr, n_phases)),
            "g_to"=>LinearAlgebra.diagm(0=>fill(g_to, n_phases)),
            "b_to"=>LinearAlgebra.diagm(0=>fill(b_to, n_phases)),
        ))
    end

    return bus_ids, line_ids, [bus_ids[bus] for bus in tr_t_bus]
end


""
function _alias!(dict, fr, to)
    if haskey(dict, fr)
        dict[to] = dict[fr]
    end
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


""
function data_model_make_compatible_v8!(data_model; phases=[1, 2, 3], neutral=4)
    data_model["conductors"] = 3
    data_model["buspairs"] = nothing
    for (_, bus) in data_model["bus"]
        bus["bus_i"] = bus["index"]
        terminals = bus["terminals"]
        @assert(all(t in [phases..., neutral] for t in terminals))
        for prop in ["vm", "va", "vmin", "vmax"]
            if haskey(bus, prop)
                if length(bus[prop])==4
                    val = bus[prop]
                    bus[prop] = val[terminals.!=neutral]
                end
            end
        end
    end

    for (_, load) in data_model["load"]
        # remove neutral
        if load["configuration"]=="wye"
            bus = data_model["bus"][string(load["bus"])]
            @assert(length(bus["grounded"])==1 && bus["grounded"][1]==load["connections"][end])
            load["connections"] = load["connections"][1:end-1]
            _pad_properties!(load, ["pd", "qd"], load["connections"], phases)
        else
            # three-phase loads can only be delta-connected
            #@assert(all(load["connections"].==phases))
        end
        _alias!(load, "bus", "load_bus")
    end

    data_model["gen"] = data_model["generator"]

    # has to be three-phase
    for (_, gen) in data_model["gen"]
        if gen["configuration"]=="wye"
            @assert(all(gen["connections"].==[phases..., neutral]))
        else
            @assert(all(gen["connections"].==phases))
        end

        _alias!(gen, "status", "gen_status")
        _alias!(gen, "bus", "gen_bus")
        _alias!(gen, "pg_min", "pmin")
        _alias!(gen, "qg_min", "qmin")
        _alias!(gen, "pg_max", "pmax")
        _alias!(gen, "qg_max", "qmax")
        _alias!(gen, "configuration", "conn")

        gen["model"] = 2
    end

    data_model["branch"] = data_model["line"]
    for (_, br) in data_model["branch"]
        @assert(all(x in phases for x in br["f_connections"]))
        @assert(all(x in phases for x in br["t_connections"]))
        @assert(all(br["f_connections"].==br["t_connections"]))

        _pad_properties!(br, ["rs", "xs", "g_fr", "g_to", "b_fr", "b_to", "s_rating", "c_rating"], br["f_connections"], phases)

        # rename
        _alias!(br, "status", "br_status")
        _alias!(br, "rs", "br_r")
        _alias!(br, "xs", "br_x")

        br["tap"] = 1.0
        br["shift"] = 0

        if !haskey(br, "angmin")
            N = size(br["br_r"])[1]
            br["angmin"] = fill(-pi/2, N)
            br["angmax"] = fill(pi/2, N)
        end
    end

    for (_, shunt) in data_model["shunt"]
        @assert(all(x in phases for x in shunt["connections"]))
        _pad_properties!(shunt, ["g_sh", "b_sh"], shunt["connections"], phases)
        _alias!(shunt, "bus", "shunt_bus")
        _alias!(shunt, "g_sh", "gs")
        _alias!(shunt, "b_sh", "bs")
    end

    data_model["dcline"] = Dict()
    data_model["transformer"] = data_model["transformer"]

    data_model["per_unit"] = true
    data_model["baseMVA"] = data_model["settings"]["sbase"]*data_model["settings"]["kv_kvar_scalar"]/1E6
    data_model["name"] = "IDC"


    return data_model
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


""
function transform_solution!(solution, data_model)
    solution_make_si!(solution, data_model)
    solution_identify!(solution, data_model)
    solution_unmap!(solution, data_model)
end
