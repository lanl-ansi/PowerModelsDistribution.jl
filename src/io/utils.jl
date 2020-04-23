import Base.Iterators: flatten

"all node types that can help define buses"
const _dss_node_objects = Vector{String}([
    "isource", "load", "generator", "indmach012", "storage", "pvsystem"
])

"all edge types that can help define buses"
const _dss_edge_objects = Vector{String}([
    "vsource", "fault", "capacitor", "line", "reactor", "transformer", "gictransformer", "gicline"
])

"all data holding objects"
const _dss_data_objects = Vector{String}([
    "options", "xfmrcode", "linecode", "loadshape", "xycurve", "linegeometry",
    "linespacing", "growthshape", "tcc_curve", "cndata", "tsdata", "wiredata"
])

"all objects that define controls"
const _dss_control_objects = Vector{String}([
    "capcontrol", "regcontrol", "swtcontrol", "relay", "recloser", "fuse"
])

"all objects that provide montoring"
const _dss_monitor_objects = Vector{String}([
    "energymeter", "monitor"
])

"components currently supported for automatic data type parsing"
const _dss_supported_components = Vector{String}([
    "line", "linecode", "load", "generator", "capacitor", "reactor",
    "transformer", "pvsystem", "storage", "loadshape", "options",
    "xfmrcode", "vsource", "xycurve"
])

"two number operators for reverse polish notation"
_double_operators = Dict{String,Any}(
    "+" => +,
    "-" => -,
    "*" => *,
    "/" => /,
    "^" => ^,
    "atan2" => (x, y) -> rad2deg(atan(y, x))
)

"single number operators in reverse polish notation"
_single_operators = Dict{String,Any}(
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
const _array_delimiters = Vector{Char}(['\"', '\'', '[', '{', '(', ']', '}', ')'])

"properties that should be excluded from being overwritten during the application of `like`"
const _like_exclusions = Dict{String,Vector{Regex}}(
    "all" => Vector{Regex}([r"name", r"enabled"]),
    "line" => [r"switch"],
)

"Regexes for determining data types"
const _dtype_regex = Dict{Regex, Type}(
    r"^[+-]{0,1}\d*\.{0,1}\d*[eE]{0,1}[+-]{0,1}\d*[+-]\d*\.{0,1}\d*[eE]{0,1}[+-]{0,1}\d*[ij]$" => ComplexF64,
    r"^[+-]{0,1}\d*\.{0,1}\d*[eE]{0,1}[+-]{0,1}\d*$" => Float64,
    r"^\d+$" => Int,
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
function _parse_matrix(dtype::Type, data::AbstractString)::Matrix{dtype}
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
        if occursin(",", clean_data) ||
            (startswith(clean_data, "[") && endswith(clean_data, "]")) ||
            (startswith(clean_data, "\"") && endswith(clean_data, "\"")) ||
            (startswith(clean_data, "\'") && endswith(clean_data, "\'")) ||
            (startswith(clean_data, "(") && endswith(clean_data, ")")) ||
            (startswith(clean_data, "{") && endswith(clean_data, "}"))
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
function _parse_array(dtype::Type, data::AbstractString)::Vector{dtype}
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
        for delim in _array_delimiters
            data = replace(data, delim => "")
        end
        elements = split(data, split_char)
        elements = [strip(el) for el in elements if strip(el) != ""]
    end

    if dtype == String || dtype == AbstractString || dtype == Char
        array = Vector{String}([])
        for el in elements
            push!(array, el)
        end
    else
        array = Vector{dtype}(undef, length(elements))
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
    if haskey(data_eng, "transformer")
        bankable_transformers = Dict()
        for (id, tr) in data_eng["transformer"]
            if haskey(tr, "bank")
                bank = tr["bank"]
                if !haskey(bankable_transformers, bank)
                    bankable_transformers[bank] = ([], [])
                end
                push!(bankable_transformers[bank][1], id)
                push!(bankable_transformers[bank][2], deepcopy(tr))
            end
        end

        for (bank, (ids, trs)) in bankable_transformers
            for tr in trs
                _apply_xfmrcode!(tr, data_eng)
            end
            # across-phase properties should be the same to be eligible for banking
            props = ["bus", "noloadloss", "xsc", "rs", "imag", "vnom", "snom", "polarity", "configuration"]
            btrans = Dict{String, Any}(prop=>trs[1][prop] for prop in props)
            if !all(tr[prop]==btrans[prop] for tr in trs, prop in props)
                Memento.warn(_LOGGER, "Not all across-phase properties match among transfomers identified by bank='$bank', aborting attempt to bank")
                continue
            end
            nrw = length(btrans["bus"])

            # only attempt to bank wye-connected transformers
            if !all(all(tr["configuration"].=="wye") for tr in trs)
                Memento.warn(_LOGGER, "Not all configurations 'wye' on transformers identified by bank='$bank', aborting attempt to bank")
                continue
            end
            neutrals = [conns[end] for conns in trs[1]["connections"]]
            # ensure all windings have the same neutral
            if !all(all(conns[end]==neutrals[w] for (w, conns) in enumerate(tr["connections"])) for tr in trs)
                Memento.warn(_LOGGER, "Not all neutral phases match on transfomers identified by bank='$bank', aborting attempt to bank")
                continue
            end

            # this will merge the per-phase properties in such a way that the
            # f_connections will be sorted from small to large
            f_phases_loc = Dict(hcat([[(c,(i,p)) for (p, c) in enumerate(tr["connections"][1][1:end-1])] for (i, tr) in enumerate(trs)]...))
            locs = [f_phases_loc[x] for x in sort(collect(keys(f_phases_loc)))]
            props_merge = ["connections", "tm_set", "tm_ub", "tm_lb", "tm_step", "tm_fix"]
            for prop in props_merge
                btrans[prop] = [[trs[i][prop][w][p] for (i,p) in locs] for w in 1:nrw]

                # for the connections, also prefix the neutral per winding
                if prop=="connections"
                    for w in 1:nrw
                        push!(btrans[prop][w], neutrals[w])
                    end
                end
            end

            btrans["source_id"] = "transformer.$bank"

            # edit the transformer dict
            for id in ids
                delete!(data_eng["transformer"], id)
            end
            data_eng["transformer"][bank] = btrans
        end
    end
end


"discovers all terminals in the network"
function _discover_terminals!(data_eng::Dict{String,<:Any})
    terminals = Dict{String, Set{Int}}([(name, Set{Int}()) for (name,bus) in data_eng["bus"]])

    if haskey(data_eng, "line")
        for (_,eng_obj) in data_eng["line"]
            # ignore 0 terminal
            push!(terminals[eng_obj["f_bus"]], setdiff(eng_obj["f_connections"], [0])...)
            push!(terminals[eng_obj["t_bus"]], setdiff(eng_obj["t_connections"], [0])...)
        end
    end

    if haskey(data_eng, "switch")
        for (_,eng_obj) in data_eng["switch"]
            # ignore 0 terminal
            push!(terminals[eng_obj["f_bus"]], setdiff(eng_obj["f_connections"], [0])...)
            push!(terminals[eng_obj["t_bus"]], setdiff(eng_obj["t_connections"], [0])...)
        end
    end

    if haskey(data_eng, "transformer")
        for (_,tr) in data_eng["transformer"]
            for w in 1:length(tr["bus"])
                # ignore 0 terminal
                push!(terminals[tr["bus"][w]], setdiff(tr["connections"][w], [0])...)
            end
        end
    end

    for comp_type in [x for x in ["voltage_source", "load", "gen"] if haskey(data_eng, x)]
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


"discovers all phases and neutrals in the network"
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


"Discovers all neutrals in the network"
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
function _get_conductors_ordered(busname::AbstractString; default::Vector{Int}=Vector{Int}([]), check_length::Bool=true, pad_ground::Bool=false)::Vector{Int}
    parts = split(busname, '.'; limit=2)
    ret = Vector{Int}([])
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


"creates a `dss` dict inside `object` that imports all items in `prop_order` from `dss_obj`"
function _import_all!(object::Dict{String,<:Any}, dss_obj::Dict{String,<:Any}, prop_order::Vector{String})
    object["dss"] = Dict{String,Any}((key, dss_obj[key]) for key in prop_order)
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


"Discovers all of the buses (not separately defined in OpenDSS), from \"lines\""
function _discover_buses(data_dss::Dict{String,<:Any})::Set
    buses = Set([])
    for obj_type in _dss_node_objects
        for (name, dss_obj) in get(data_dss, obj_type, Dict{String,Any}())
            _apply_like!(dss_obj, data_dss, obj_type)
            push!(buses, split(dss_obj["bus1"], '.'; limit=2)[1])
        end
    end

    for obj_type in _dss_edge_objects
        for (name, dss_obj) in get(data_dss, obj_type, Dict{String,Any}())
            _apply_like!(dss_obj, data_dss, obj_type)
            if obj_type == "transformer"
                transformer = _create_transformer(name; _to_kwargs(dss_obj)...)
                for bus in transformer["buses"]
                    push!(buses, split(bus, '.'; limit=2)[1])
                end
            elseif obj_type == "gictransformer"
                for key in ["bush", "busx", "busnh", "busnx"]
                    if haskey(dss_obj, key)
                        push!(buses, split(dss_obj[key], '.'; limit=2)[1])
                    end
                end
            elseif obj_type == "vsource"
                push!(buses, split(get(dss_obj, "bus1", "sourcebus"), '.'; limit=2)[1])
                if haskey(dss_obj, "bus2")
                    push!(buses, split(dss_obj["bus2"], '.'; limit=2)[1])
                end
            else
                for key in ["bus1", "bus2"]
                    if haskey(dss_obj, key)
                        push!(buses, split(dss_obj[key], '.'; limit=2)[1])
                    end
                end
            end
        end
    end

    return buses
end


"shifts a vector by `shift` spots to the left"
function _barrel_roll(x::Vector{T}, shift::Int)::Vector{T} where T
    N = length(x)
    if shift < 0
        shift = shift + ceil(Int, shift/N)*N
    end

    shift = mod(shift, N)

    return x[[(i-1+shift)%N+1 for i in 1:N]]
end


"Parses busnames as defined in OpenDSS, e.g. \"primary.1.2.3.0\""
function _parse_busname(busname::AbstractString)::Tuple{String,Vector{Bool}}
    parts = split(busname, '.'; limit=2)
    name = parts[1]
    elements = "1.2.3"

    if length(parts) >= 2
        name, elements = split(busname, '.'; limit=2)
    end

    nodes = Vector{Bool}([0, 0, 0, 0])

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


"converts Dict{String,Any} to Dict{Symbol,Any} for passing as kwargs"
function _to_kwargs(data::Dict{String,Any})::Dict{Symbol,Any}
    return Dict{Symbol,Any}((Symbol(k), v) for (k, v) in data)
end


"apply properties in the order that they are given"
function _apply_ordered_properties(defaults::Dict{String,<:Any}, raw_dss::Dict{String,<:Any}; code_dict::Dict{String,<:Any}=Dict{String,Any}())::Dict{String,Any}
    _defaults = deepcopy(defaults)

    for prop in filter(p->p!="like", raw_dss["prop_order"])
        if prop in ["linecode", "loadshape"]
            merge!(defaults, code_dict)
        else
            if haskey(_defaults, prop)
                defaults[prop] = _defaults[prop]
            end
        end
    end

    return defaults
end


"applies `like` to component"
function _apply_like!(raw_dss::Dict{String,<:Any}, data_dss::Dict{String,<:Any}, comp_type::String)
    links = ["like"]
    if any(link in raw_dss["prop_order"] for link in links)
        new_prop_order = Vector{String}([])
        raw_dss_copy = deepcopy(raw_dss)

        for prop in raw_dss["prop_order"]
            push!(new_prop_order, prop)

            if any(match(key, prop) !== nothing for key in [get(_like_exclusions, comp_type, [])..., _like_exclusions["all"]...])
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

        final_prop_order = Vector{String}([])
        while !isempty(new_prop_order)
            prop = popfirst!(new_prop_order)
            if !(prop in new_prop_order)
                push!(final_prop_order, prop)
            end
        end
        raw_dss["prop_order"] = final_prop_order
    end
end


"""
Parses the data in keys defined by `to_parse` in `data_dss` using types given by
the default properties from the `get_prop_default` function.
"""
function _parse_dss_with_dtypes!(data_dss::Dict{String,<:Any}, to_parse::Vector{String}=_dss_supported_components)
    for obj_type in to_parse
        if haskey(data_dss, obj_type)
            dtypes = _dss_parameter_data_types[obj_type]
            if obj_type == "options"
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
function _parse_element_with_dtype(dtype::Type, element::AbstractString)
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


"parses data type of properties of objects"
function _parse_obj_dtypes!(obj_type::String, object::Dict{String,Any}, dtypes::Dict{String,Type})
    for (k, v) in object
        if isa(v, Vector) && eltype(v) == Any || isa(eltype(v), AbstractString)
            _dtype = get(dtypes, k, _guess_dtype("[$(join(v, ","))]"))
            for i in 1:length(v)
                if isa(v[i], AbstractString)
                    v[i] = _parse_element_with_dtype(_dtype, v[i])
                end
            end
        elseif isa(v, Matrix) && eltype(v) == Any || isa(eltype(v), AbstractString)
            _dtype = get(dtypes, k, _guess_dtype("$(join(collect(flatten(v)), " "))"))
            for i in 1:size(v)[1]
                for j in 1:size(v)[2]
                    if isa(v[i,j], AbstractString)
                        v[i,j] = _parse_element_with_dtype(_dtype, v[i,j])
                    end
                end
            end
        elseif isa(v, AbstractString)
            object[k] = _parse_element_with_dtype(get(dtypes, k, _guess_dtype(v)), v)
        end
    end
end


""
function _register_awaiting_ground!(bus::Dict{String,<:Any}, connections::Vector{Int})
    if !haskey(bus, "awaiting_ground")
        bus["awaiting_ground"] = []
    end

    push!(bus["awaiting_ground"], connections)
end


"checks to see if a property is after linecode"
function _is_after_linecode(prop_order::Vector{String}, property::String)::Bool
    return _is_after(prop_order, property, "linecode")
end


"checks to see if a property is after xfmrcode"
function _is_after_xfmrcode(prop_order::Vector{String}, property::String)::Bool
    return _is_after(prop_order, property, "xfmrcode")
end


"checks to see if property1 is after property2 in the prop_order"
function _is_after(prop_order::Vector{String}, property1::String, property2::String)::Bool
    property1_idx = 0
    property2_idx = 0

    for (i, prop) in enumerate(prop_order)
        if prop == property1
            property1_idx = i
        elseif prop == property2
            property2_idx = i
        end
    end

    return property1_idx > property2_idx
end


"add engineering data object to engineering data model"
function _add_eng_obj!(data_eng::Dict{String,<:Any}, eng_obj_type::String, eng_obj_id::Any, eng_obj::Dict{String,<:Any})
    if !haskey(data_eng, eng_obj_type)
        data_eng[eng_obj_type] = Dict{Any,Any}()
    end

    data_eng[eng_obj_type][eng_obj_id] = eng_obj
end


"guesses the data type of a value using regex, returning Float64, Int, ComplexF64, or String (if number type cannot be determined)"
function _guess_dtype(value::AbstractString)::Type
    if _isa_matrix(value) || _isa_array(value) || _isa_rpn(value)
        for delim in [keys(_double_operators)..., keys(_single_operators)..., _array_delimiters..., "|", ","]
            value = replace(value, delim => " ")
        end
        _dtypes = unique([_guess_dtype(v) for v in split(value)])
        if length(_dtypes) == 1
            return _dtypes[1]
        elseif all(isa(v, Int) for v in _dtypes)
            return Int
        elseif any(isa(v, Complex) for v in _dtypes)
            return ComplexF64
        elseif any(isa(v, Float64) for v in _dtypes)
            return Float64
        else
            return String
        end
    else
        for (re, typ) in _dtype_regex
            if occursin(re, value)
                return typ
            end
        end

        return String
    end
end
