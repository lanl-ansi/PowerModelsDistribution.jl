"Squares `x`, for parsing Reverse Polish Notation"
function _sqr(x::Float64)
    return x * x
end


_double_operators = Dict("+" => +, "-" => -, "*" => *, "/" => /, "^" => ^, "atan2" => (x, y) -> rad2deg(atan2(y, x)))

_single_operators = Dict("sqr" => _sqr, "sqrt" => sqrt, "inv" => inv, "ln" => log,
                         "exp" => exp, "log10" => log10, "sin" => sind, "cos" => cosd,
                         "tan" => tand, "asin" => asind, "acos" => acosd, "atan" => atand)

const _array_delimiters = ['\"', '\'', '[', '{', '(', ']', '}', ')']

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


"checks is a string is a connection by checking the values"
function _isa_conn(expr::AbstractString)::Bool
    if expr in ["wye", "y", "ln", "delta", "ll"]
        return true
    else
        return false
    end
end


"""
    _get_prop_name(obj_type)

Returns the property names in order for a given component type `obj_type`.
"""
function _get_prop_name(obj_type::AbstractString)::Array
    linecode = ["nphases", "r1", "x1", "r0", "x0", "c1", "c0", "units",
                "rmatrix", "xmatrix", "cmatrix", "basefreq", "normamps",
                "emergamps", "faultrate", "pctperm", "repair", "kron",
                "rg", "xg", "rho", "neutral", "b1", "b0", "like"]

    linegeometry = ["nconds", "nphases", "cond", "wire", "x", "h", "units",
                    "normamps", "emergamps", "reduce", "spacing", "wires",
                    "cncable", "tscable", "cncables", "tscables", "like"]

    linespacing = ["nconds", "nphases", "x", "h", "units"]

    loadshape = ["npts", "interval", "minterval", "sinterval", "pmult",
                 "qmult", "hour", "mean", "stddev", "csvfile", "sngfile",
                 "pqcsvfile", "action", "useactual", "pmax", "qmax",
                 "pbase", "like"]

    growthshape = ["npts", "year", "mult", "csvfile", "sngfile", "dblfile",
                   "like"]

    tcc_curve = ["npts", "c_array", "t_array", "like"]

    cndata = []

    tsdata = []

    wiredata = ["rdc", "rac", "runits", "gmrac", "gmrunits", "radius",
                "radunits", "normamps", "emergamps", "diam", "like"]

    xfmrcode = []

    vsource = ["bus1", "bus2", "basekv", "pu", "angle", "frequency", "phases",
               "mvasc3", "mvasc1", "x1r1", "x0r0", "isc3", "isc1", "r1", "x1",
               "r0", "x0", "scantype", "sequence", "spectrum", "z1", "z2", "z0",
               "puz1", "puz2", "puz0", "basemva", "basefreq", "like", "enabled"]

    isource = ["phases", "bus1", "amps", "angle", "frequency", "scantype",
               "sequence", "spectrum", "basefreq", "enabled", "like"]

    fault = ["phases", "bus1", "bus2", "r", "gmatrix", "minamps", "ontime",
             "pctperm", "temporary", "%stddev", "normamps", "emergamps",
             "basefreq", "faultrate", "repair", "enabled", "like"]

    capacitor = ["bus1", "bus2", "phases", "kvar", "kv", "conn", "cmatrix",
                 "cuf", "r", "xl", "harm", "numsteps", "states", "normamps",
                 "emergamps", "faultrate", "pctperm", "basefreq", "enabled",
                 "like"]

    line = ["bus1", "bus2", "linecode", "length", "phases",
            "r1", "x1", "r0", "x0", "c1", "c0", "b1",
            "b0", "normamps", "emergamps", "faultrate", "pctperm",
            "repair", "basefreq", "rmatrix", "xmatrix", "cmatrix",
            "switch", "rg", "xg", "rho", "geometry", "earthmodel",
            "units", "enabled", "like"]

    reactor = ["phases", "bus1", "bus2", "kv", "kvar", "conn", "parallel",
               "r", "rmatrix", "rp", "x", "xmatrix", "z", "z1", "z2", "z0",
               "rcurve", "lcurve", "lmh", "normamps", "emergamps", "repair",
               "faultrate", "pctperm", "basefreq", "enabled", "like"]

    transformer = ["phases", "windings", "wdg", "bus", "conn", "kv", "kva",
                   "tap", "%r", "rneut", "xneut", "buses", "conns", "kvs",
                   "kvas", "taps", "%rs", "xhl", "xlt", "xht", "xscarray",
                   "thermal", "n", "m", "flrise", "hsrise", "%loadloss",
                   "%noloadloss", "%imag", "ppm_antifloat", "normhkva",
                   "emerghkva", "sub", "maxtap", "mintap", "numtaps",
                   "subname", "bank", "xfmrcode", "xrconst", "leadlag",
                   "faultrate", "basefreq", "enabled", "like"]

    gictransformer = ["basefreq", "bush", "busnh", "busnx", "busx",
                      "emergamps", "enabled", "phases", "r1", "r2", "type",
                      "mva", "kvll1", "kvll2", "%r1", "%r2", "k", "varcurve",
                      "like", "normamps", "emergamps", "pctperm", "repair"]

    gicline = ["angle", "bus1", "bus2", "c", "ee", "en", "frequency", "lat1",
               "lat2", "lon1", "lon2", "phases", "r", "volts", "x", "like",
               "basefreq", "enabled", "spectrum"]

    load = ["phases", "bus1", "kv", "kw", "pf", "model",
            "yearly", "daily", "duty", "growth", "conn", "kvar",
            "rneut", "xneut", "status", "class", "vminpu", "vmaxpu",
            "vminnorm", "vminemerg", "xfkva", "allocationfactor",
            "kva", "%mean", "%stddev", "cvrwatts", "cvrvars", "kwh",
            "kwhdays", "cfactor", "cvrcurve", "numcust", "spectrum",
            "zipv", "%seriesrl", "relweight", "vlowpu", "puxharm",
            "xrharm", "spectrum", "basefreq", "enabled", "like"]

    generator = ["bus1", "phases", "kv", "kw", "pf", "model", "yearly",
                 "daily", "duty", "dispvalue", "conn", "kvar", "rneut",
                 "xneut", "status", "class", "vpu", "maxkvar", "minkvar",
                 "pvfactor", "debugtrace", "vminpu", "vmaxpu", "forceon",
                 "kva", "mva", "xd", "xdp", "xdpp", "h", "d", "usermodel",
                 "userdata", "shaftmodel", "shaftdata", "dutystart", "balanced",
                 "xrdp", "spectrum", "basefreq", "enabled", "like"]

    indmach012 = ["phases", "bus1", "kv", "kw", "pf", "conn", "kva", "h",
                  "d", "purs", "puxs", "purr", "puxr", "puxm", "slip",
                  "maxslip", "slipoption", "spectrum", "enabled"]

    storage = ["phases", "bus1", "%charge", "%discharge", "%effcharge", "%idlingkvar",
               "idlingkw", "%r", "%reserve", "%stored", "%x", "basefreq",
               "chargetrigger", "class", "conn", "daily", "yearly", "debugtrace",
               "dischargetrigger", "dispmode", "duty", "dynadata", "dynadll", "enabled",
               "kv", "kva", "kvar", "kw", "kwhrated", "kwhstored", "kwrated", "like",
               "model", "pf", "spectrum", "state", "timechargetrig", "userdata",
               "usermodel", "vmaxpu", "vminpu", "yearly"]

    capcontrol = ["element", "capacitor", "type", "ctphase", "ctratio", "deadtime",
                  "delay", "delayoff", "eventlog", "offsetting", "onsetting",
                  "ptphase", "ptratio", "terminal", "vbus", "vmax", "vmin",
                  "voltoverride", "enabled"]

    regcontrol = ["transformer", "winding", "vreg", "band", "delay", "ptratio",
                  "ctprim", "r", "x", "pthase", "tapwinding", "bus",
                  "remoteptratio", "debugtrace", "eventlog", "inversetime",
                  "maxtapchange", "revband", "revdelay", "reversible",
                  "revneutral", "revr", "revthreshold", "revvreg", "revx",
                  "tapdelay", "tapnum", "vlimit", "ldc_z", "rev_z", "cogen",
                  "enabled"]

    energymeter = ["element", "terminal", "action", "clear", "save", "take",
                   "option", "kwnorm", "kwemerg", "peakcurrent", "zonelist",
                   "zonelist", "zonelist", "localonly", "mask", "losses",
                   "linelosses", "xfmrlosses", "seqlosses", "3phaselosses",
                   "vbaselosses", "basefreq", "enabled", "like"]

   pvsystem =    ["phases", "bus1", "kv", "irradiance", "pmpp", "temperature",
                  "pf", "conn", "kvar", "kva", "%cutin", "%cutout",
                  "effcurve", "p-tcurve", "%r", "%x", "model", "vminpu",
                  "vmaxpu", "yearly", "daily", "duty", "tyearly", "tduty",
                  "class", "usermodel", "userdata", "debugtrace", "spectrum"]

    obj_types = Dict{String, Array}("linecode" => linecode,
                                 "linegeometry" => linegeometry,
                                 "linespacing" =>linespacing,
                                 "loadshape" => loadshape,
                                 "growthshape" => growthshape,
                                 "tcc_curve" => tcc_curve,
                                 "cndata" => cndata,
                                 "tsdata" => tsdata,
                                 "wiredata" => wiredata,
                                 "xfmrcode" => xfmrcode,
                                 "vsource" => vsource,
                                 "isource" => isource,
                                 "fault" => fault,
                                 "capacitor" => capacitor,
                                 "line" => line,
                                 "reactor" => reactor,
                                 "transformer" => transformer,
                                 "gictransformer" => gictransformer,
                                 "gicline" => gicline,
                                 "load" => load,
                                 "generator" => generator,
                                 "indmach012" => indmach012,
                                 "storage" => storage,
                                 "capcontrol" => capcontrol,
                                 "regcontrol" => regcontrol,
                                 "energymeter" => energymeter,
                                 "pvsystem" => pvsystem,
                                 "circuit" => vsource
                                )

    return obj_types[obj_type]
end


"""
    _parse_matrix(dtype, data)

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


"parse matrices according to active nodes"
function _parse_matrix(data::Array{T}, nodes::Array{Bool}, nph::Int=3, fill_val=0.0)::Array where T
    mat = fill(fill_val, (nph, nph))
    idxs = findall(nodes[1:nph])

    for i in 1:size(idxs)[1]
        mat[idxs[i], idxs[i]] = data[i, i]
        for j in 1:i-1
            mat[idxs[i],idxs[j]] = mat[idxs[j],idxs[i]] = data[i, j]
        end
    end

    return mat
end


"checks if `data` is an opendss-style matrix string"
function _isa_matrix(data::AbstractString)::Bool
    if occursin("|", data)
        return true
    else
        return false
    end
end


"Reorders a `matrix` based on the order that phases are listed in on the from- (`pof`) and to-sides (`pot`)"
function _reorder_matrix(matrix, phase_order)
    mat = zeros(size(matrix))
    for (i, n) in zip(sort(phase_order), phase_order)
        for (j, m) in zip(sort(phase_order), phase_order)
            mat[i, j] = matrix[n, m]
        end
    end
    return mat
end


"""
    _parse_array(dtype, data)

Parses a OpenDSS style array string `data` into a one dimensional array of type
`dtype`. Array strings are capped by either brackets, single quotes, or double
quotes, and elements are separated by spaces.
"""
function _parse_array(dtype::Type, data::AbstractString)
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
        array = []
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


"parse matrices according to active nodes"
function _parse_array(data, nodes::Array{Bool}, nph::Int=3, fill_val=0.0)::Array
    mat = fill(fill_val, nph)
    idxs = findall(nodes[1:nph])

    if length(data) == 1 && nph > 1
        for i in idxs
            mat[i] = data[1]
        end
    else
        for i in 1:length(idxs)
            mat[idxs[i]] = data[i]
        end
    end

    return mat
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


"parses single column load profile files"
function _parse_loadshape_csv(path::AbstractString, type::AbstractString; header::Bool=false, column::Int=1, interval::Bool=false)
    open(path, "r") do f
        lines = readlines(f)
        if header
            lines = lines[2:end]
        end

        if type == "mult"
            return [split(line, ",")[column] for line in lines]
        elseif type == "csvfile"
            if interval
                hour, mult = [], []
                for line in lines
                    d = split(line, ",")
                    push!(hour, d[1])
                    push!(mult, d[2])
                end
                return hour, mult
            else
                return [split(line, ",")[1] for line in lines]
            end
        elseif type == "pqcsvfile"
            if interval
                hour, pmult, qmult = [], [], []
                for line in lines
                    d = split(line, ",")
                    push!(hour, d[1])
                    push!(pmult, d[2])
                    push!(qmult, d[3])
                end

                return hour, pmult, qmult
            else
                pmult, qmult = [], []
                for line in lines
                    d = split(line, ",")
                    push!(pmult, d[1])
                    push!(qmult, d[2])
                end

                return pmult, qmult
            end
        end
    end
end


"parses sng and dbl precision loadshape binary files"
function _parse_binary(path::AbstractString, precision::Type; npts::Union{Int,Nothing}=nothing, interval::Bool=false)
    open(path, "r") do f
        if npts === nothing
            data = precision[]
            while true
                try
                    n = read(f, precision)
                    push!(data, n)
                catch EOFError
                    break
                end
            end

            return data
        else
            data = Array{precision, 1}(undef, interval ? npts * 2 : npts)

            try
                read!(f, data)
            catch EOFError
                error("Error reading binary file: likely npts is wrong")
            end

            if interval
                data = reshape(data, 2, :)
                return data[1, :], data[2, :]
            else
                return data
            end
        end
    end
end


"parses csv or binary loadshape files"
function _parse_loadshape_file(path::AbstractString, type::AbstractString, npts::Union{Int,Nothing}; header::Bool=false, interval::Bool=false, column::Int=1)
    if type in ["csvfile", "mult", "pqcsvfile"]
        return _parse_loadshape_csv(path, type; header=header, column=column, interval=interval)
    elseif type in ["sngfile", "dblfile"]
        return _parse_binary(path, Dict("sngfile" => Float32, "dblfile" => Float64)[type]; npts=npts, interval=interval)
    end
end


"parses pmult and qmult entries on loadshapes"
function _parse_mult(mult_string::AbstractString; path::AbstractString="", npts::Union{Int,Nothing}=nothing)
    if !occursin("=", mult_string)
        return mult_string
    else
        props = Dict{String,Any}(split(p, "=") for p in _parse_properties(mult_string[2:end-1]))
        if haskey(props, "header")
            header = props["header"]
            props["header"] = header == "yes" || header == "true" ? true : false
        end

        if haskey(props, "col")
            props["column"] = pop!(props, "col")
        end

        file_key = [prop for prop in keys(props) if endswith(prop, "file")][1]
        full_path = path == "" ? props[file_key] : join([path, props[file_key]], '/')
        type = file_key == "file" ? "mult" : file_key

        return "($(join(_parse_loadshape_file(full_path, type, npts; header=get(props, "header", false), column=parse(Int, get(props, "column", "1"))), ",")))"
    end
end


"parses loadshape component"
function _parse_loadshape!(current_obj::Dict{String,Any}; path::AbstractString="")
    if any(parse.(Float64, [get(current_obj, "interval", "1.0"), get(current_obj, "minterval", "60.0"), get(current_obj, "sinterval", "3600.0")]) .<= 0.0)
        interval = true
    else
        interval = false
    end

    npts = parse(Int, get(current_obj, "npts", "1"))

    for prop in current_obj["prop_order"]
        if prop in ["pmult", "qmult"]
             current_obj[prop] = _parse_mult(current_obj[prop]; path=path, npts=npts)
        elseif prop in ["csvfile", "pqcsvfile", "sngfile", "dblfile"]
            full_path = path == "" ? current_obj[prop] : join([path, current_obj[prop]], '/')
            data = _parse_loadshape_file(full_path, prop, parse(Int, get(current_obj, "npts", "1")); interval=interval, header=false)
            if prop == "pqcsvfile"
                if interval
                    current_obj["hour"], current_obj["pmult"], current_obj["qmult"] = data
                else
                    current_obj["pmult"], current_obj["qmult"] = data
                end
            else
                if interval
                    current_obj["hour"], current_obj["pmult"] = data
                else
                    current_obj["pmult"] = data
                end
            end
        end
    end

    for prop in ["pmult", "qmult", "hour"]
        if haskey(current_obj, prop) && isa(current_obj[prop], Array)
            current_obj[prop] = "($(join(current_obj[prop], ",")))"
        elseif haskey(current_obj, prop) && isa(current_obj[prop], String) && !_isa_array(current_obj[prop])
            current_obj[prop] = "($(current_obj[prop]))"
        end
    end
end


"strips lines that are either commented (block or single) or empty"
function _strip_lines(lines::Array)::Array
    blockComment = false
    stripped_lines = []
    for line in lines
        if startswith(strip(line), "/*") || endswith(strip(line), "*/")
            blockComment = !blockComment
        elseif !startswith(line, '!') && !startswith(line, "//") && strip(line) != "" && !blockComment
            push!(stripped_lines, line)
        end
    end
    return stripped_lines
end


"""
    _parse_buscoords(file)

Parses a Bus Coordinate `file`, in either "dat" or "csv" formats, where in
"dat", columns are separated by spaces, and in "csv" by commas. File expected
to contain "bus,x,y" on each line.
"""
function _parse_buscoords(file::AbstractString)::Dict{String,Any}
    file_str = read(open(file), String)
    regex = r",\s*"
    if endswith(lowercase(file), "csv") || endswith(lowercase(file), "dss")
        regex = r","
    end

    lines = _strip_lines(split(file_str, '\n'))

    buscoords = Dict{String,Dict{String,Any}}()
    for line in lines
        bus, x, y = split(line, regex; limit=3)
        buscoords[lowercase(strip(bus, [',']))] = Dict{String,Any}("x"=>parse(Float64, strip(x, [','])), "y"=>parse(Float64, strip(y, [',', '\r'])))
    end

    return buscoords
end


"""
    _parse_properties(properties)

Parses a string of `properties` of a component type, character by character
into an array with each element containing (if present) the property name, "=",
and the property value.
"""
function _parse_properties(properties::AbstractString)::Array{<:Any, 1}
    parsed_properties = []
    end_array = true
    end_property = false
    end_equality = false
    end_single_quote = true
    end_double_quote = true
    string_out = ""

    properties = replace(properties, r"\s*=\s*" => "=")
    num_chars = length(properties)

    for (n, char) in enumerate(properties)
        substring_out = split(string_out, "=")
        if length(substring_out) == 2 && substring_out[2] != ""
            end_equality = true
        elseif !occursin("=", string_out) && (char == ' ' || n == num_chars)
            end_equality = true
        elseif occursin("=", string_out) && (char == ' ' || n == num_chars) && end_array
            end_equality = true
        else
            end_equality = false
        end

        if char == ' '
            end_property = true
        else
            end_property = false
        end

        if char in ['[', '(', '{']
            end_array = false
        elseif char in [']', ')', '}']
            end_array = true
        end

        if char == '\"'
            end_double_quote = !end_double_quote
        elseif char == '\''
            end_single_quote = !end_single_quote
        end

        if char != ' ' || !end_array || !end_double_quote || !end_single_quote
            string_out = string(string_out, char)
        end

        if string_out != "" && end_array && end_property && end_equality && end_double_quote && end_single_quote || n == num_chars
            push!(parsed_properties, string_out)
            string_out = ""
        end
    end

    return parsed_properties
end


"""
    _add_component!(data_dss, obj_type_name, object)

Adds a component of type `obj_type_name` with properties given by `object` to
the existing `data_dss` structure. If a component of the same type has already
been added to `data_dss`, the new component is appeneded to the existing array
of components of that type, otherwise a new array is created.
"""
function _add_component!(data_dss::Dict, obj_type_name::AbstractString, object::Dict)
    obj_type = split(obj_type_name, '.'; limit=2)[1]
    if obj_type == "circuit"
        if haskey(data_dss, "circuit")
            Memento.error(_LOGGER, "Cannot have two circuits, invalid dss")
        else
            data_dss[obj_type] = object
        end
    elseif haskey(data_dss, obj_type)
        data_dss[obj_type][object["name"]] = object
    else
        data_dss[obj_type] = Dict{String,Any}(object["name"] => object)
    end
end


"""
    _add_property(object, key, value)

Adds a property to an existing component properties dictionary `object` given
the `key` and `value` of the property. If a property of the same name already
exists inside `object`, the original value is converted to an array, and the
new value is appended to the end.
"""
function _add_property(object::Dict, key::AbstractString, value::Any)::Dict
    if !haskey(object, "prop_order")
        object["prop_order"] = Array{String,1}(["name"])
    end

    current_wdg = "wdg" in object["prop_order"] ? string(filter(p->occursin("wdg", p), object["prop_order"])[end][end]) : ""
    current_wdg = current_wdg == "g" ? "" : current_wdg

    if key in ["wdg", "bus", "conn", "kv", "kva", "tap", "%r", "rneut", "xneut"]
        key = join(filter(p->!isempty(p), [key, current_wdg]), "_")
    end

    if haskey(object, lowercase(key))
        rmatch = match(r"_(\d+)$", key)
        if typeof(rmatch) != Nothing
            end_num = parse(Int, rmatch.captures[1]) + 1
            key = replace(key, r"_(\d+)$" => "_$end_num")
        else
            key = string(key, "_2")
        end
    end

    object[lowercase(key)] = value
    push!(object["prop_order"], lowercase(key))

    return object
end


"""
    _parse_component(component, properies, object=Dict{String,Any}())

Parses a `component` with `properties` into a `object`. If `object` is not
defined, an empty dictionary will be used. Assumes that unnamed properties are
given in order, but named properties can be given anywhere.
"""
function _parse_component(component::AbstractString, properties::AbstractString, object::Dict{String,<:Any}=Dict{String,Any}(); path::AbstractString="")
    Memento.debug(_LOGGER, "Properties: $properties")
    obj_type, name = split(component, '.'; limit=2)

    if !haskey(object, "name")
        object["name"] = name
    end

    property_array = _parse_properties(properties)
    Memento.debug(_LOGGER, "property_array: $property_array")

    property_names = _get_prop_name(obj_type)
    property_idx = 1

    for (n, property) in enumerate(property_array)
        if property == ""
            continue
        elseif !occursin("=", property)
            property = join([property_names[property_idx], property], '=')
            property_idx += 1
        else
            if obj_type == "loadshape" && startswith(property, "mult")
                property = replace(property, "mult" => "pmult")
            elseif obj_type == "transformer"
                prop_name, _ = split(property,'=')
                if prop_name == "ppm"
                    property = replace(property, prop_name => "ppm_antifloat")
                elseif prop_name == "x12"
                    property = replace(property, prop_name => "xhl")
                elseif prop_name == "x23"
                    property = replace(property, prop_name => "xlt")
                elseif prop_name == "x13"
                    property = replace(property, prop_name => "xht")
                end
            end

            property_idxs = findall(e->e==split(property,'=')[1], property_names)
            if length(property_idxs) > 0
                property_idx = findall(e->e==split(property,'=')[1], property_names)[1] + 1
            end
        end

        key, value = split(property, '='; limit=2)

        if occursin(r"\(\s*(sng|dbl)*file=(.+)\)", value)
            value = _parse_mult(value; path=path)
        end

        _add_property(object, key, value)
    end

    return object
end


"""
    _merge_dss!(dss_prime, dss_to_add)

Merges two (partially) parsed OpenDSS files to the same dictionary `dss_prime`.
Used in cases where files are referenced via the "compile" or "redirect"
OpenDSS commands inside the originating file.
"""
function _merge_dss!(dss_prime::Dict{String,<:Any}, dss_to_add::Dict{String,<:Any})
    for (k, v) in dss_to_add
        if k in keys(dss_prime)
            append!(dss_prime[k], v)
        else
            dss_prime[k] = v
        end
    end
end


"""
    _parse_line(elements, current_obj=Dict{String,Any}())

Parses an already separated line given by `elements` (an array) of an OpenDSS
file into `current_obj`. If not defined, `current_obj` is an empty dictionary.
"""
function _parse_line(elements::Array, current_obj::Dict=Dict{String,Any}(); path::AbstractString="")
    current_obj_type = strip(elements[2], ['\"', '\''])
    if startswith(current_obj_type, "object")
        current_obj_type = split(current_obj_type, '=')[2]
        current_obj["name"] = split(current_obj_type, '.')[2]
    else
        if length(elements) != 3
            properties = ""
        else
            properties = elements[3]
        end

        current_obj = _parse_component(current_obj_type, properties; path=path)
    end

    return current_obj_type, current_obj
end


"Strips comments, defined by \"!\" from the ends of lines"
function _strip_comments(line::AbstractString)::String
    return strip(split(line, r"\s*!")[1], ['\r', '\n'])
end


"""
    _assign_property!(data_dss, obj_type, obj_name, property_name, property_value)

Assigns a property with name `property_name` and value `property_value` to the component
of type `obj_type` named `obj_name` in `data_dss`.
"""
function _assign_property!(data_dss::Dict, obj_type::AbstractString, obj_name::AbstractString, property_name::AbstractString, property_value::Any)
    if haskey(data_dss, obj_type)
        for obj in data_dss[obj_type]
            if obj["name"] == obj_name
                obj[property_name] = property_value
            end
        end
    else
        Memento.warn(_LOGGER, "Cannot find $obj_type object $obj_name.")
    end
end


"""
    parse_dss(filename)

Parses a OpenDSS file given by `filename` into a Dict{Array{Dict}}. Only
supports components and options, but not commands, e.g. "plot" or "solve".
Will also parse files defined inside of the originating DSS file via the
"compile", "redirect" or "buscoords" commands.
"""
function parse_dss(io::IOStream)::Dict{String,Any}
    filename = match(r"^<file\s(.+)>$", io.name).captures[1]
    Memento.info(_LOGGER, "Calling parse_dss on $filename")
    current_file = split(filename, "/")[end]
    path = join(split(filename, '/')[1:end-1], '/')
    data_dss = Dict{String,Any}()

    data_dss["filename"] = [string(current_file)]

    current_obj = Dict{String,Any}()
    current_obj_type = ""

    lines = readlines(io)

    stripped_lines = _strip_lines(lines)
    nlines = length(stripped_lines)

    for (n, line) in enumerate(stripped_lines)
        real_line_num = findall(lines .== line)[1]
        Memento.debug(_LOGGER, "LINE $real_line_num: $line")
        line = _strip_comments(line)

        if startswith(strip(line), '~')
            current_obj = _parse_component(current_obj_type, strip(strip(lowercase(line)),  '~'), current_obj)

            if n < nlines && startswith(strip(stripped_lines[n + 1]), '~')
                continue
            else
                _add_component!(data_dss, current_obj_type, current_obj)
            end
        else
            current_obj = Dict{String,Any}()
            line_elements = split(line, r"\s+"; limit=3)
            cmd = lowercase(line_elements[1])

            if cmd == "clear"
                Memento.info(_LOGGER, "`data_dss` has been reset with the \"clear\" command.")
                data_dss = Dict{String,Any}("filename"=>data_dss["filename"])
                continue

            elseif cmd == "redirect"
                file = line_elements[2]
                full_path = path == "" ? file : join([path, file], '/')
                Memento.info(_LOGGER, "Redirecting to file \"$file\"")
                _merge_dss!(data_dss, parse_dss(full_path))
                continue

            elseif cmd == "compile"
                file = split(strip(line_elements[2], ['(',')']), '\\')[end]
                full_path = path == "" ? file : join([path, file], '/')
                Memento.info(_LOGGER, "Compiling file \"$file\"")
                _merge_dss!(data_dss, parse_dss(full_path))
                continue

            elseif cmd == "set"
                Memento.debug(_LOGGER, "set command: $line_elements")
                properties = _parse_properties(join(line_elements[2:end], " "))
                if length(line_elements) == 2
                    property, value = split(lowercase(line_elements[2]), '='; limit=2)
                else
                    property, value = split(lowercase(join(line_elements[2:end], " ")), '='; limit=2)
                end

                if !haskey(data_dss, "options")
                    data_dss["options"] = Dict{String,Any}()
                end

                data_dss["options"]["$(property)"] = value
                continue

            elseif cmd == "buscoords"
                file = line_elements[2]
                full_path = path == "" ? file : join([path, file], '/')
                Memento.debug(_LOGGER, "Buscoords path: $full_path")
                data_dss["buscoords"] = _parse_buscoords(full_path)

            elseif cmd == "new"
                current_obj_type, current_obj = _parse_line([lowercase(line_element) for line_element in line_elements]; path=path)

                if startswith(current_obj_type, "loadshape")
                    _parse_loadshape!(current_obj; path=path)
                end
            else
                try
                    obj_type, obj_name, props = split(lowercase(line), '.'; limit=3)
                    parsed_properties = _parse_properties(props)
                    wdg = ""
                    for prop in parsed_properties
                        property_name, property_value = split(prop, '=')
                        if obj_type == "transformer"
                            wdg = property_name == "wdg" && property_value != "1" ? property_value : property_name == "wdg" && property_value == "1" ? "" : wdg

                            if property_name in ["wdg", "bus", "conn", "kv", "kva", "tap", "%r", "rneut", "xneut"]
                                property_name = join(filter(p->!isempty(p), [property_name, wdg]), "_")
                            end

                            _assign_property!(data_dss, obj_type, obj_name, property_name, property_value)
                        else
                            _assign_property!(data_dss, obj_type, obj_name, property_name, property_value)
                        end
                    end
                catch
                    Memento.warn(_LOGGER, "Command \"$cmd\" on line $real_line_num in \"$current_file\" is not supported, skipping.")
                end
            end

            Memento.debug(_LOGGER, "size current_obj: $(length(current_obj))")
            if n < nlines && startswith(strip(stripped_lines[n + 1]), '~')
                continue
            elseif length(current_obj) > 0
                _add_component!(data_dss, current_obj_type, current_obj)
            else
                continue
            end
        end
    end

    Memento.info(_LOGGER, "Done parsing $filename")
    return data_dss
end


""
function parse_dss(filename::AbstractString)::Dict
    data_dss = open(filename) do io
        parse_dss(io)
    end
    return data_dss
end


"parses the raw dss values into their expected data types"
function _parse_element_with_dtype(dtype, element)
    if _isa_matrix(element)
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


"""
    parse_dss_with_dtypes!(data_dss, to_parse)

Parses the data in keys defined by `to_parse` in `data_dss` using types given by
the default properties from the `get_prop_default` function.
"""
function parse_dss_with_dtypes!(data_dss::Dict{String,<:Any}, to_parse::Array{String}=[])
    for obj_type in to_parse
        if haskey(data_dss, obj_type)
            Memento.debug(_LOGGER, "type: $obj_type")
            dtypes = _get_dtypes(obj_type)
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


""
function _parse_obj_dtypes!(obj_type, object, dtypes)
    for (k, v) in object
        if haskey(dtypes, k)
            Memento.debug(_LOGGER, "key: $k")
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
            else
                Memento.error(_LOGGER, "dtype unknown $obj_type, $k, $v")
            end
        end
    end
end


"""
    _parse_busname(busname)

Parses busnames as defined in OpenDSS, e.g. "primary.1.2.3.0".
"""
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


"""
    _get_conductors_ordered(busname; neutral=true)

Returns an ordered list of defined conductors. If ground=false, will omit any `0`
"""
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
function _to_sym_keys(data::Dict{String,Any})::Dict{Symbol,Any}
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


"properties that should be excluded from being overwritten during the application of `like`"
const _like_exclusions = Dict{String,Array}("all" => ["name", "bus1", "bus2", "phases", "nphases", "enabled"],
                                            "line" => ["switch"],
                                            "transformer" => ["bank", "bus", "bus_2", "bus_3", "buses", "windings", "wdg", "wdg_2", "wdg_3"],
                                            "linegeometry" => ["nconds"]
                                            )
