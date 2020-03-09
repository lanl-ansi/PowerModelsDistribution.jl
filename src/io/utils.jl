"all edge types that can help define buses"
const _dss_edge_components = ["line", "transformer", "reactor"]

"components currently supported for automatic data type parsing"
const _dss_supported_components = ["line", "linecode", "load", "generator", "capacitor", "reactor", "circuit", "transformer", "pvsystem", "storage", "loadshape"]

"two number operators for reverse polish notation"
_double_operators = Dict(
    "+" => +,
    "-" => -,
    "*" => *,
    "/" => /,
    "^" => ^,
    "atan2" => (x, y) -> rad2deg(atan(y, x))
)

"single number operators in reverse polish notation"
_single_operators = Dict(
    "sqr" => x -> x * x,
    "sqrt" => sqrt,
    "inv" => inv,
    "ln" => log,
    "exp" => exp,
    "log10" => log10,
    "sin" => sind,
    "cos" => cosd,
    "tan" => tand,
    "asin" => asind,
    "acos" => acosd,
    "atan" => atand
)

"different acceptable delimiters for arrays"
const _array_delimiters = ['\"', '\'', '[', '{', '(', ']', '}', ')']

"properties that should be excluded from being overwritten during the application of `like`"
const _like_exclusions = Dict{String,Array}(
    "all" => ["name", "bus1", "bus2", "phases", "nphases", "enabled"],
    "line" => ["switch"],
    "transformer" => ["bank", "bus", "bus_2", "bus_3", "buses", "windings", "wdg", "wdg_2", "wdg_3"],
    "linegeometry" => ["nconds"]
)

"data types of various dss option inputs"
const _dss_option_dtypes = Dict{String,Type}(
    "defaultbasefreq" => Float64,
    "voltagebases" => Float64,
    "tolerance" => Float64
)


"detects if `expr` is Reverse Polish Notation expression"
function _isa_rpn(expr::AbstractString)::Bool
    expr = split(strip(expr, _array_delimiters))
    op_keys = keys(merge(_double_operators, _single_operators))
    for item in expr
        if item in op_keys
            return true
        end
    end
    return false
end


"parses Reverse Polish Notation `expr`"
function _parse_rpn(expr::AbstractString, dtype::Type=Float64)
    clean_expr = strip(expr, _array_delimiters)

    if occursin("rollup", clean_expr) || occursin("rolldn", clean_expr) || occursin("swap", clean_expr)
        Memento.warn(_LOGGER, "_parse_rpn does not support \"rollup\", \"rolldn\", or \"swap\", leaving as String")
        return expr
    end

    stack = []
    split_expr = occursin(",", clean_expr) ? split(clean_expr, ',') : split(clean_expr)

    for item in split_expr
        try
            if haskey(_double_operators, item)
                b = pop!(stack)
                a = pop!(stack)
                push!(stack, _double_operators[item](a, b))
            elseif haskey(_single_operators, item)
                push!(stack, _single_operators[item](pop!(stack)))
            else
                if item == "pi"
                    push!(stack, pi)
                else
                    push!(stack, parse(dtype, item))
                end
            end
        catch error
            if isa(error, ArgumentError)
                Memento.warn(_LOGGER, "\"$expr\" is not valid Reverse Polish Notation, leaving as String")
                return expr
            else
                throw(error)
            end
        end
    end
    if length(stack) > 1
        Memento.warn(_LOGGER, "\"$expr\" is not valid Reverse Polish Notation, leaving as String")
        return expr
    else
        return stack[1]
    end
end


"checks is a string is a connection by checking the values"
function _isa_conn(expr::AbstractString)::Bool
    if expr in ["wye", "y", "ln", "delta", "ll"]
        return true
    else
        return false
    end
end


"parses connection \"conn\" specification reducing to wye or delta"
function _parse_conn(conn::String)::String
    if conn in ["wye", "y", "ln"]
        return "wye"
    elseif conn in ["delta", "ll"]
        return "delta"
    else
        Memento.warn(_LOGGER, "Unsupported connection $conn, defaulting to \"wye\"")
        return "wye"
    end
end


"checks if `data` is an opendss-style matrix string"
function _isa_matrix(data::AbstractString)::Bool
    if occursin("|", data)
        return true
    else
        return false
    end
end


"""
Parses a OpenDSS style triangular matrix string `data` into a two dimensional
array of type `dtype`. Matrix strings are capped by either parenthesis or
brackets, rows are separated by "|", and columns are separated by spaces.
"""
function _parse_matrix(dtype::Type, data::AbstractString)::Array
    rows = []
    for line in split(strip(data, _array_delimiters), '|')
        cols = []
        for item in split(line)
            push!(cols, parse(dtype, item))
        end
        push!(rows, cols)
    end

    nphases = maximum([length(row) for row in rows])

    if dtype == AbstractString || dtype == String
        matrix = fill("", nphases, nphases)
    elseif dtype == Char
        matrix = fill(' ', nphases, nphases)
    else
        matrix = zeros(dtype, nphases, nphases)
    end

    if length(rows) == 1
        for i in 1:nphases
            matrix[i, i] = rows[1][1]
        end
    elseif all([length(row) for row in rows] .== [i for i in 1:nphases])
        for (i, row) in enumerate(rows)
            for (j, col) in enumerate(row)
                matrix[i, j] = matrix[j, i] = col
            end
        end
    elseif all([length(row) for row in rows] .== nphases)
        for (i, row) in enumerate(rows)
            for (j, col) in enumerate(row)
                matrix[i, j] = col
            end
        end
    end

    return matrix
end


"checks if `data` is an opendss-style array string"
function _isa_array(data::AbstractString)::Bool
    clean_data = strip(data)
    if !occursin("|", clean_data)
        if occursin(",", clean_data)
            return true
        elseif startswith(clean_data, "[") && endswith(clean_data, "]")
            return true
        elseif startswith(clean_data, "\"") && endswith(clean_data, "\"")
            return true
        elseif startswith(clean_data, "\'") && endswith(clean_data, "\'")
            return true
        elseif startswith(clean_data, "(") && endswith(clean_data, ")")
            return true
        elseif startswith(clean_data, "{") && endswith(clean_data, "}")
            return true
        else
            return false
        end
    else
        return false
    end
end


"""
Parses a OpenDSS style array string `data` into a one dimensional array of type
`dtype`. Array strings are capped by either brackets, single quotes, or double
quotes, and elements are separated by spaces.
"""
function _parse_array(dtype::Type, data::AbstractString)::Vector
    if occursin(",", data)
        split_char = ','
    else
        split_char = ' '
    end

    if _isa_rpn(data)
        matches = collect((m.match for m = eachmatch(Regex(string("[",join(_array_delimiters, '\\'),"]")), data, overlap=false)))
        if length(matches) == 2
            if dtype == String
                return data
            else
                return _parse_rpn(data, dtype)
            end

        else
            elements = _parse_properties(data[2:end-1])
        end
    else
        elements = split(strip(data, _array_delimiters), split_char)
        elements = [strip(el) for el in elements if strip(el) != ""]
    end

    if dtype == String || dtype == AbstractString || dtype == Char
        array = Vector{String}([])
        for el in elements
            push!(array, el)
        end
    else
        array = zeros(dtype, length(elements))
        for (i, el) in enumerate(elements)
            if _isa_rpn(data)
                array[i] = _parse_rpn(el, dtype)
            else
                array[i] = parse(dtype, el)
            end
        end
    end

    return array
end


"Combines transformers with 'bank' keyword into a single transformer"
function _bank_transformers!(data_eng::Dict{String,<:Any})
    banks = Dict{String,Array{String,1}}()
    for (name, transformer) in get(data_eng, "transformer", Dict{String,Any}())
        if haskey(transformer, "bank")
            if !haskey(banks, transformer["bank"])
                banks[transformer["bank"]] = Array{String,1}([name])
            else
                push!(banks[transformer["bank"]], name)
            end
        end
    end

    banked = Dict{String,Any}()
    for (bank, names) in banks
        transformers = [data_eng["transformer"][name] for name in names]

        total_phases = sum(Int[transformer["nphases"] for transformer in transformers])

        if !haskey(banked, bank)
            eng_obj = Dict{String,Any}(
                "source_id" => "transformer.$bank"
            )

            for transformer in transformers
                for key in ["polarity", "bus", "configuration", "status", "bank", "rs", "noloadloss", "tm_step", "xsc", "imag", "snom", "vnom"]
                    if !haskey(eng_obj, key)
                        eng_obj[key] = transformer[key]
                    else
                        @assert transformer[key] == eng_obj[key] "$key property of transformers does not match, cannot bank"
                    end
                end

                if !haskey(eng_obj, "connections")
                    eng_obj["connections"] = transformer["connections"]
                else
                    @assert length(transformer["connections"]) == length(eng_obj["connections"]) "inconsistent number of windings, cannot bank"

                    for (i, wdg) in enumerate(transformer["connections"])
                        for conn in wdg
                            if !(conn in eng_obj["connections"][i])
                                push!(eng_obj["connections"][i], conn)
                            end
                        end
                        eng_obj["connections"][i] = _roll(sort(eng_obj["connections"][i]), count(j->j==0, eng_obj["connections"][i]); right=false)
                    end
                end

                if !haskey(eng_obj, "nphases")
                    eng_obj["nphases"] = transformer["nphases"]
                else
                    eng_obj["nphases"] += transformer["nphases"]
                end

                for key in ["tm", "tm_max", "tm_min", "fixed"]
                    if !haskey(eng_obj, key)
                        eng_obj[key] = [zeros(3), zeros(3)]
                    end

                    for (i, wdg) in enumerate(transformer[key])
                        for (j, v) in enumerate(wdg)
                            eng_obj[key][i][transformer["connections"][i][j]] = v
                        end
                    end
                end

            end

            banked[bank] = eng_obj
        end
    end

    for (bank, names) in banks
        if haskey(banked, bank)
            for name in names
                delete!(data_eng["transformer"], name)
            end
        end

        data_eng["transformer"][bank] = banked[bank]
    end
end


""
function _discover_terminals!(data_eng::Dict{String,<:Any})
    terminals = Dict{String, Set{Int}}([(name, Set{Int}()) for (name,bus) in data_eng["bus"]])

    for (_,eng_obj) in data_eng["line"]
        # ignore 0 terminal
        push!(terminals[eng_obj["f_bus"]], setdiff(eng_obj["f_connections"], [0])...)
        push!(terminals[eng_obj["t_bus"]], setdiff(eng_obj["t_connections"], [0])...)
    end

    if haskey(data_eng, "transformer")
        for (_,tr) in data_eng["transformer"]
            for w in 1:length(tr["bus"])
                # ignore 0 terminal
                push!(terminals[tr["bus"][w]], setdiff(tr["connections"][w], [0])...)
            end
        end
    end

    for comp_type in [x for x in ["voltage_source", "load", "shunt", "gen"] if haskey(data_eng, x)]
        for comp in values(data_eng[comp_type])
            push!(terminals[comp["bus"]], setdiff(comp["connections"], [0])...)
        end
    end

    for (id, bus) in data_eng["bus"]
        data_eng["bus"][id]["terminals"] = sort(collect(terminals[id]))
    end

    for (id,bus) in data_eng["bus"]
        if haskey(bus, "awaiting_ground")
            neutral = !(4 in bus["terminals"]) ? 4 : maximum(bus["terminals"])+1
            push!(bus["terminals"], neutral)

            bus["grounded"] = [neutral]
            bus["rg"] = [0.0]
            bus["xg"] = [0.0]
            for i in 1:length(bus["awaiting_ground"])
                bus["awaiting_ground"][i][bus["awaiting_ground"][i].==0] .= neutral
            end

            delete!(bus, "awaiting_ground")
        end
    end
end


function _discover_phases_neutral!(data_eng::Dict{String,<:Any})
    bus_neutral = _find_neutrals(data_eng)
    for (id, bus) in data_eng["bus"]
        terminals = bus["terminals"]
        if haskey(bus_neutral, id)
            bus["neutral"] = bus_neutral[id]
            phases = setdiff(terminals, bus["neutral"])
        else
            phases = terminals
        end
        @assert(length(phases)<=3, "At bus $id, we found $(length(phases))>3 phases; aborting discovery, requires manual inspection.")
    end
end


""
function _find_neutrals(data_eng::Dict{String,<:Any})
    vertices = [(id, t) for (id, bus) in data_eng["bus"] for t in bus["terminals"]]
    neutrals = []
    edges = Set([((eng_obj["f_bus"], eng_obj["f_connections"][c]),(eng_obj["t_bus"], eng_obj["t_connections"][c])) for (id, eng_obj) in data_eng["line"] for c in 1:length(eng_obj["f_connections"])])

    bus_neutrals = [(id,bus["neutral"]) for (id,bus) in data_eng["bus"] if haskey(bus, "neutral")]
    trans_neutrals = []
    for (_, tr) in data_eng["transformer"]
        for w in 1:length(tr["connections"])
            if tr["configuration"][w] == "wye"
                push!(trans_neutrals, (tr["bus"][w], tr["connections"][w][end]))
            end
        end
    end
    load_neutrals = [(eng_obj["bus"],eng_obj["connections"][end]) for (_,eng_obj) in get(data_eng, "load", Dict{String,Any}()) if eng_obj["configuration"]=="wye"]
    neutrals = Set(vcat(bus_neutrals, trans_neutrals, load_neutrals))
    neutrals = Set([(bus,t) for (bus,t) in neutrals if t!=0])
    stack = deepcopy(neutrals)
    while !isempty(stack)
        vertex = pop!(stack)
        candidates_t = [((f,t), t) for (f,t) in edges if f==vertex]
        candidates_f = [((f,t), f) for (f,t) in edges if t==vertex]
        for (edge,next) in [candidates_t..., candidates_f...]
            delete!(edges, edge)
            push!(stack, next)
            push!(neutrals, next)
        end
    end
    bus_neutral = Dict{String, Int}()
    for (bus,t) in neutrals
        bus_neutral[bus] = t
    end
    return bus_neutral
end


"Returns an ordered list of defined conductors. If ground=false, will omit any `0`"
function _get_conductors_ordered_dm(busname::AbstractString; default::Array=[], check_length::Bool=true, pad_ground::Bool=false)::Array
    parts = split(busname, '.'; limit=2)
    ret = []
    if length(parts)==2
        conds_str = split(parts[2], '.')
        ret = [parse(Int, i) for i in conds_str]
    else
        return default
    end

    if pad_ground && length(ret)==length(default)-1
        ret = [ret..., 0]
    end

    if check_length && length(default)!=length(ret)
        # TODO
        Memento.warn(_LOGGER, "An incorrect number of nodes was specified on $(parts[1]); |$(parts[2])|!=$(length(default)).")
    end
    return ret
end


""
function _import_all!(component::Dict{String,<:Any}, defaults::Dict{String,<:Any}, prop_order::Array{<:AbstractString,1})
    component["dss"] = Dict{String,Any}((key, defaults[key]) for key in prop_order)
end


"""
Given a vector and a list of elements to find, this method will return a list
of the positions of the elements in that vector.
"""
function _get_idxs(vec::Vector{<:Any}, els::Vector{<:Any})::Vector{Int}
    ret = Array{Int, 1}(undef, length(els))
    for (i,f) in enumerate(els)
        for (j,l) in enumerate(vec)
            if f==l
                ret[i] = j
            end
        end
    end
    return ret
end


function _get_ilocs(vec::Vector{<:Any}, loc::Any)::Vector{Int}
    return collect(1:length(vec))[vec.==loc]
end


"Discovers all of the buses (not separately defined in OpenDSS), from \"lines\""
function _discover_buses(data_dss::Dict{String,<:Any})::Set
    buses = Set([])
    for obj_type in _dss_edge_components
        for (name, dss_obj) in get(data_dss, obj_type, Dict{String,Any}())
            if obj_type == "transformer"
                dss_obj = _create_transformer(dss_obj["name"]; _to_kwargs(dss_obj)...)
                for bus in dss_obj["buses"]
                    push!(buses, split(bus, '.'; limit=2)[1])
                end
            elseif haskey(dss_obj, "bus2")
                for key in ["bus1", "bus2"]
                    push!(buses, split(dss_obj[key], '.'; limit=2)[1])
                end
            end
        end
    end

    return buses
end


""
function _barrel_roll(x::Vector, shift)
    N = length(x)
    if shift < 0
        shift = shift + ceil(Int, shift/N)*N
    end

    shift = mod(shift, N)

    return x[[(i-1+shift)%N+1 for i in 1:N]]
end


"Parses busnames as defined in OpenDSS, e.g. \"primary.1.2.3.0\""
function _parse_busname(busname::AbstractString)
    parts = split(busname, '.'; limit=2)
    name = parts[1]
    elements = "1.2.3"

    if length(parts) >= 2
        name, elements = split(busname, '.'; limit=2)
    end

    nodes = Array{Bool}([0 0 0 0])

    for num in 1:3
        if occursin("$num", elements)
            nodes[num] = true
        end
    end

    if occursin("0", elements) || sum(nodes[1:3]) == 1
        nodes[4] = true
    end

    return name, nodes
end


"Returns an ordered list of defined conductors. If ground=false, will omit any `0`"
function _get_conductors_ordered(busname::AbstractString; neutral::Bool=true, nconductors::Int=3)::Array
    parts = split(busname, '.'; limit=2)
    ret = []
    if length(parts)==2
        conductors_string = split(parts[2], '.')
        if neutral
            ret = [parse(Int, i) for i in conductors_string]
        else
            ret = [parse(Int, i) for i in conductors_string if i != "0"]
        end
    else
        ret = collect(1:nconductors)
    end

    return ret
end


"converts Dict{String,Any} to Dict{Symbol,Any} for passing as kwargs"
function _to_kwargs(data::Dict{String,Any})::Dict{Symbol,Any}
    return Dict{Symbol,Any}((Symbol(k), v) for (k, v) in data)
end


""
function _apply_ordered_properties(defaults::Dict{String,<:Any}, raw_dss::Dict{String,<:Any}; code_dict::Dict{String,<:Any}=Dict{String,Any}())
    _defaults = deepcopy(defaults)

    for prop in filter(p->p!="like", raw_dss["prop_order"])
        if prop in ["linecode", "loadshape"]
            merge!(defaults, code_dict)
        else
            defaults[prop] = _defaults[prop]
        end
    end

    return defaults
end


"applies `like` to component"
function _apply_like!(raw_dss, data_dss, comp_type)
    links = ["like"]
    if any(link in raw_dss["prop_order"] for link in links)
        new_prop_order = []
        raw_dss_copy = deepcopy(raw_dss)

        for prop in raw_dss["prop_order"]
            push!(new_prop_order, prop)

            if prop in get(_like_exclusions, comp_type, []) || prop in _like_exclusions["all"]
                continue
            end

            if prop in links
                linked_dss = get(get(data_dss, comp_type, Dict{String,Any}()), raw_dss[prop], Dict{String,Any}())
                if isempty(linked_dss)
                    Memento.warn(_LOGGER, "$comp_type.$(raw_dss["name"]): $prop=$(raw_dss[prop]) cannot be found")
                else
                    for linked_prop in linked_dss["prop_order"]
                        if linked_prop in get(_like_exclusions, comp_type, []) || linked_prop in _like_exclusions["all"]
                            continue
                        end

                        push!(new_prop_order, linked_prop)
                        if linked_prop in links
                            _apply_like!(linked_dss, data_dss, comp_type)
                        else
                            raw_dss[linked_prop] = deepcopy(linked_dss[linked_prop])
                        end
                    end
                end
            else
                raw_dss[prop] = deepcopy(raw_dss_copy[prop])
            end
        end

        final_prop_order = []
        while !isempty(new_prop_order)
            prop = popfirst!(new_prop_order)
            if !(prop in new_prop_order)
                push!(final_prop_order, prop)
            end
        end
        raw_dss["prop_order"] = final_prop_order
    end
end


"Parses options defined with the `set` command in OpenDSS"
function parse_dss_options!(data_dss::Dict{String,<:Any})
    if haskey(data_dss, "options")
        for (k,v) in data_dss["options"]
            if haskey(_dss_option_dtypes, k)
                dtype = _dss_option_dtypes[k]
                if _isa_array(v)
                    data_dss["options"][k] = _parse_array(dtype, v)
                elseif _isa_matrix(v)
                    data_dss["options"][k] = _parse_matrix(dtype, v)
                else
                    data_dss["options"][k] = parse(dtype, v)
                end
            end
        end
    end
end


"""
    parse_dss_with_dtypes!(data_dss, to_parse)

Parses the data in keys defined by `to_parse` in `data_dss` using types given by
the default properties from the `get_prop_default` function.
"""
function parse_dss_with_dtypes!(data_dss::Dict{String,<:Any}, to_parse::Array{String}=_dss_supported_components)
    for obj_type in to_parse
        if haskey(data_dss, obj_type)
            dtypes = _dss_parameter_data_types[obj_type]
            if obj_type == "circuit"
                _parse_obj_dtypes!(obj_type, data_dss[obj_type], dtypes)
            else
                for object in values(data_dss[obj_type])
                    _parse_obj_dtypes!(obj_type, object, dtypes)
                end
            end
        end
    end
end


"parses the raw dss values into their expected data types"
function _parse_element_with_dtype(dtype, element)
    if _isa_rpn(element)
        out = _parse_rpn(element, dtype)
    elseif _isa_matrix(element)
        out = _parse_matrix(eltype(dtype), element)
    elseif _isa_array(element)
        out = _parse_array(eltype(dtype), element)
    elseif dtype <: Bool
        if element in ["n", "no"]
            element = "false"
        elseif element in ["y", "yes"]
            element = "true"
        end
        out = parse(dtype, element)
    elseif _isa_rpn(element)
        out = _parse_rpn(element)
    elseif dtype == String
        out = element
    else
        if _isa_conn(element)
            out = _parse_conn(element)
        else
            try
                out = parse(dtype, element)
            catch
                Memento.warn(_LOGGER, "cannot parse $element as $dtype, leaving as String.")
                out = element
            end
        end
    end

    return out
end


""
function _parse_obj_dtypes!(obj_type, object, dtypes)
    for (k, v) in object
        if haskey(dtypes, k)
            if isa(v, Array)
                arrout = []
                for el in v
                    if isa(v, AbstractString)
                        push!(arrout, _parse_element_with_dtype(dtypes[k], el))
                    else
                        push!(arrout, el)
                    end
                end
                object[k] = arrout
            elseif isa(v, AbstractString)
                object[k] = _parse_element_with_dtype(dtypes[k], v)
            end
        end
    end
end


"""
Given a set of addmittances 'y' connected from the conductors 'f_cnds' to the
conductors 't_cnds', this method will return a list of conductors 'cnd' and a
matrix 'Y', which will satisfy I[cnds] = Y*V[cnds].
"""
function _calc_shunt(f_cnds, t_cnds, y)
    # TODO add types
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
function _calc_ground_shunt_admittance_matrix(cnds, Y, ground)
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
function _rm_floating_cnd(cnds, Y, f)
    # TODO add types
    P = setdiff(cnds, f)
    f_inds = _get_idxs(cnds, [f])
    P_inds = _get_idxs(cnds, P)
    Yrm = Y[P_inds,P_inds]-(1/Y[f_inds,f_inds][1])*Y[P_inds,f_inds]*Y[f_inds,P_inds]
    return (P,Yrm)
end


function _register_awaiting_ground!(bus, connections)
    if !haskey(bus, "awaiting_ground")
        bus["awaiting_ground"] = []
    end
    push!(bus["awaiting_ground"], connections)
end
