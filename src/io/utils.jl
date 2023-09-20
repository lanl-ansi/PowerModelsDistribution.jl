"all node types that can help define buses"
const _dss_node_objects = String[
    "isource", "load", "generator", "indmach012", "storage", "pvsystem"
]

"all edge types that can help define buses"
const _dss_edge_objects = String[
    "vsource", "fault", "capacitor", "line", "reactor", "transformer", "gictransformer", "gicline"
]

"all data holding objects"
const _dss_data_objects = String[
    "options", "xfmrcode", "linecode", "loadshape", "xycurve", "linegeometry",
    "linespacing", "growthshape", "tcc_curve", "cndata", "tsdata", "wiredata"
]

"all objects that define controls"
const _dss_control_objects = String[
    "capcontrol", "regcontrol", "swtcontrol", "relay", "recloser", "fuse"
]

"all objects that provide montoring"
const _dss_monitor_objects = String[
    "energymeter", "monitor"
]

"components currently supported for automatic data type parsing"
const _dss_supported_components = String[
    "line", "linecode", "load", "generator", "capacitor", "reactor",
    "transformer", "pvsystem", "storage", "loadshape", "options",
    "xfmrcode", "vsource", "xycurve", "spectrum", "capcontrol",
    "regcontrol", "linegeometry", "wiredata", "linespacing",
    "cndata", "tsdata"
]

"two number operators for reverse polish notation"
const _double_operators = Dict{String,Function}(
    "+" => +,
    "-" => -,
    "*" => *,
    "/" => /,
    "^" => ^,
    "atan2" => (x, y) -> rad2deg(atan(y, x))
)

"single number operators in reverse polish notation"
const _single_operators = Dict{String,Function}(
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

"dss to pmd load model"
const _dss2pmd_load_model = Dict{Int,LoadModel}(
    1 => POWER,
    2 => IMPEDANCE,
    5 => CURRENT,
    4 => EXPONENTIAL,
    8 => ZIP,
)

"dss to pmd capcontrol type"
const _dss2pmd_capcontrol_type = Dict{String,CapControlType}(
    "kvar" => CAP_REACTIVE_POWER,
    "current" => CAP_CURRENT,
    "voltage" => CAP_VOLTAGE,
    ""=> CAP_DISABLED,
    "time"=>CAP_TIME,
)

"conversion factors for units to meters"
const _convert_to_meters = Dict{String,Float64}(
    "mi" => 1609.3,
    "km" => 1000.0,
    "kft" => 304.8,
    "m" => 1.0,
    "ft" => 0.3048,
    "in" => 0.0254,
    "cm" => 0.01,
    "mm" => 0.001,
    "none" => 1.0
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


"helper function to parse reverse polish notation"
function _parse_rpn(expr::AbstractString)::Union{Float64, Vector{Float64}, AbstractString}
    if _isa_dss_array(expr)
        parse(Vector{Float64}, "(\"24.9 3 sqrt /\" \"10 2 *\")")
    else
        _parse_rpn(Float64, expr)
    end
end


"helper function to parse reverse polish notation vectors"
function _parse_rpn(::Type{T}, expr::AbstractString)::Union{T,AbstractString} where T <: Vector
    parse(T, expr)
end


"helper function to parse reverse polish notation arrays"
function _parse_array(::Type{T}, expr::AbstractString)::Union{Vector{T},AbstractString} where T
    parse(Vector{T}, expr)
end


"parses Reverse Polish Notation `expr`"
function _parse_rpn(::Type{T}, expr::AbstractString)::Union{T,AbstractString} where T
    clean_expr = strip(expr, _array_delimiters)

    if occursin("rollup", clean_expr) || occursin("rolldn", clean_expr) || occursin("swap", clean_expr)
        @warn "_parse_rpn does not support 'rollup', 'rolldn', or 'swap', leaving as String"
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
                    push!(stack, parse(T, item))
                end
            end
        catch error
            if isa(error, ArgumentError)
                @warn "'$expr' is not valid Reverse Polish Notation, leaving as String"
                return expr
            else
                throw(error)
            end
        end
    end
    if length(stack) > 1
        @warn "'$expr' is not valid Reverse Polish Notation, leaving as String"
        return expr
    else
        return stack[1]
    end
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
            props = ["bus", "noloadloss", "xsc", "rw", "cmag", "vm_nom", "sm_nom", "polarity", "configuration", "sm_ub"]
            btrans = Dict{String, Any}(prop=>trs[1][prop] for prop in props)
            if !all(tr[prop]==btrans[prop] for tr in trs, prop in props)
                @warn "Not all across-phase properties match among transfomers identified by bank='$bank', aborting attempt to bank"
                continue
            end
            nrw = length(btrans["bus"])

            # only attempt to bank wye-connected transformers
            if !all(all(conf==WYE for conf in tr["configuration"]) for tr in trs)
                @warn "Not all configurations 'wye' on transformers identified by bank='$bank', aborting attempt to bank"
                continue
            end
            neutrals = [conns[end] for conns in trs[1]["connections"]]
            # ensure all windings have the same neutral
            if !all(all(conns[end]==neutrals[w] for (w, conns) in enumerate(tr["connections"])) for tr in trs)
                @warn "Not all neutral phases match on transfomers identified by bank='$bank', aborting attempt to bank"
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

            btrans["status"] = all(tr["status"] == ENABLED for tr in trs) ? ENABLED : DISABLED
            btrans["source_id"] = "transformer.$bank"

            # add regulator objects if present
            if any([!isempty(get(tr,"controls", Dict{String,Any}())) for tr in trs])
                btrans["controls"]  = Dict{String,Any}(
                    "vreg" => [[0.0 for (i,p) in locs] for w in 1:nrw],
                    "band" => [[0.0 for (i,p) in locs] for w in 1:nrw],
                    "ptratio" => [[0.0 for (i,p) in locs] for w in 1:nrw],
                    "ctprim" => [[0.0 for (i,p) in locs] for w in 1:nrw],
                    "r" => [[0.0 for (i,p) in locs] for w in 1:nrw],
                    "x" => [[0.0 for (i,p) in locs] for w in 1:nrw]
                )
            end
            for (c,(i,p)) in enumerate(locs)
                if haskey(trs[i],"controls")
                    for w in 1:nrw
                        btrans["controls"]["vreg"][w][c] = trs[i]["controls"]["vreg"][w][p]
                        btrans["controls"]["band"][w][c] = trs[i]["controls"]["band"][w][p]
                        btrans["controls"]["ptratio"][w][c] = trs[i]["controls"]["ptratio"][w][p]
                        btrans["controls"]["ctprim"][w][c] = trs[i]["controls"]["ctprim"][w][p]
                        btrans["controls"]["r"][w][c] = trs[i]["controls"]["r"][w][p]
                        btrans["controls"]["x"][w][c] = trs[i]["controls"]["x"][w][p]
                    end
                end
            end

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
            !all(eng_obj["f_connections"] .== 0) && push!(terminals[eng_obj["f_bus"]], setdiff(eng_obj["f_connections"], [0])...)
            !all(eng_obj["t_connections"] .== 0) && push!(terminals[eng_obj["t_bus"]], setdiff(eng_obj["t_connections"], [0])...)
        end
    end

    if haskey(data_eng, "switch")
        for (_,eng_obj) in data_eng["switch"]
            # ignore 0 terminal
            !all(eng_obj["f_connections"] .== 0) && push!(terminals[eng_obj["f_bus"]], setdiff(eng_obj["f_connections"], [0])...)
            !all(eng_obj["t_connections"] .== 0) && push!(terminals[eng_obj["t_bus"]], setdiff(eng_obj["t_connections"], [0])...)
        end
    end

    if haskey(data_eng, "transformer")
        for (_,tr) in data_eng["transformer"]
            for w in 1:length(tr["bus"])
                # ignore 0 terminal
                tr["connections"][w]!=[0] && push!(terminals[tr["bus"][w]], setdiff(tr["connections"][w], [0])...)
            end
        end
    end

    for comp_type in [x for x in ["voltage_source", "load", "generator", "solar"] if haskey(data_eng, x)]
        for comp in values(data_eng[comp_type])
            !all(comp["connections"] .== 0) && push!(terminals[comp["bus"]], setdiff(comp["connections"], [0])...)
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


"Returns an ordered list of defined conductors. If ground=false, will omit any `0`"
function _get_conductors_ordered(busname::AbstractString; default::Vector{Int}=Int[], check_length::Bool=true, pad_ground::Bool=false)::Vector{Int}
    parts = split(busname, '.'; limit=2)
    ret = Int[]
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
        @info "An inconsistent number of nodes was specified on $(parts[1]); |$(parts[2])|!=$(length(default))."
    end

    return ret
end


"creates a `dss` dict inside `object` that imports all items in `prop_order` from `dss_obj`"
function _import_all!(object::Dict{String,<:Any}, dss_obj::DssObject)
    object["dss"] = Dict{String,Any}((key, property) for (key,property) in dss_obj["raw_dss"])
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


"Discovers all of the buses (not separately defined in OpenDSS), from 'lines'"
function _discover_buses(data_dss::OpenDssDataModel)::Set
    buses = Set([])
    for obj_type in _dss_node_objects
        for (name, dss_obj) in get(data_dss, obj_type, Dict{String,Any}())
            push!(buses, split(dss_obj["bus1"], '.'; limit=2)[1])
        end
    end

    for obj_type in _dss_edge_objects
        for (name, dss_obj) in get(data_dss, obj_type, Dict{String,Any}())
            if obj_type == "transformer"
                for bus in dss_obj["buses"]
                    push!(buses, split(bus, '.'; limit=2)[1])
                end
            elseif obj_type == "gictransformer"
                for key in ["bush", "busx", "busnh", "busnx"]
                    if !isempty(dss_obj[key])
                        push!(buses, split(dss_obj[key], '.'; limit=2)[1])
                    end
                end
            elseif obj_type == "vsource"
                push!(buses, split(get(dss_obj, "bus1", "sourcebus"), '.'; limit=2)[1])
                if !isempty(dss_obj["bus2"])
                    push!(buses, split(dss_obj["bus2"], '.'; limit=2)[1])
                end
            else
                for key in ["bus1", "bus2"]
                    if !isempty(dss_obj[key])
                        push!(buses, split(dss_obj[key], '.'; limit=2)[1])
                    end
                end
            end
        end
    end

    return filter(x->!isempty(x), buses)
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


"Parses busnames as defined in OpenDSS, e.g. 'primary.1.2.3.0'"
function _parse_bus_id(busname::String)::Tuple{String,Vector{Bool}}
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


""
function _register_awaiting_ground!(bus::Dict{String,<:Any}, connections::Vector{Int})
    if !haskey(bus, "awaiting_ground")
        bus["awaiting_ground"] = []
    end

    push!(bus["awaiting_ground"], connections)
end


"add engineering data object to engineering data model"
function _add_eng_obj!(data_eng::Dict{String,<:Any}, eng_obj_type::String, eng_obj_id::Any, eng_obj::Dict{String,<:Any})
    if !haskey(data_eng, eng_obj_type)
        data_eng[eng_obj_type] = Dict{String,Any}()
    end

    if haskey(data_eng[eng_obj_type], eng_obj_id)
        @warn "id '$eng_obj_id' already exists in $eng_obj_type, renaming to '$(eng_obj["source_id"])'"
        eng_obj_id = eng_obj["source_id"]
    end

    data_eng[eng_obj_type][eng_obj_id] = eng_obj
end


"checks if loadshape has both pmult and qmult"
function _is_loadshape_split(dss_obj::Dict{String,<:Any})
    haskey(dss_obj, "pmult") && haskey(dss_obj, "qmult") && all(dss_obj["pmult"] .!= dss_obj["qmult"])
end


"checks if loadshape has both pmult and qmult"
function _is_loadshape_split(dss_obj::DssLoadshape)
    !isempty(dss_obj["pmult"]) && !isempty(dss_obj["qmult"]) && all(dss_obj["pmult"] .!= dss_obj["qmult"])
end


"helper function to properly reference time series variables from opendss"
function _build_time_series_reference!(eng_obj::Dict{String,<:Any}, dss_obj::DssTimeSeriesObjects, data_dss::OpenDssDataModel, time_series::String, active::String, reactive::String)
    if !isempty(dss_obj[time_series]) && !isempty(data_dss["loadshape"]) && haskey(data_dss["loadshape"], dss_obj[time_series])
        eng_obj["time_series"] = get(eng_obj, "time_series", Dict{String,Any}())
        if _is_loadshape_split(data_dss["loadshape"][dss_obj[time_series]])
            eng_obj["time_series"][active] = "$(dss_obj[time_series])_p"
            eng_obj["time_series"][reactive] = "$(dss_obj[time_series])_q"
        else
            eng_obj["time_series"][active] = dss_obj[time_series]
            eng_obj["time_series"][reactive] = dss_obj[time_series]
        end
    end
end
