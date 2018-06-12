"Squares `x`, for parsing Reverse Polish Notation"
function sqr(x::Float64)
    return x * x
end


double_operators = Dict("+" => +, "-" => -, "*" => *, "/" => /, "^" => ^,
                        "atan2" => f(x, y) = rad2deg(atan2(y, x)))

single_operators = Dict("sqr" => sqr, "sqrt" => sqrt, "inv" => inv, "ln" => log,
                        "exp" => exp, "log10" => log10, "sin" => sind, "cos" => cosd,
                        "tan" => tand, "asin" => asind, "acos" => acosd, "atan" => atand)


"parses Reverse Polish Notation `expr`"
function parse_rpn(expr::AbstractString)
    clean_expr = lowercase(strip(expr, ['(', ')']))

    if contains(clean_expr, "rollup") || contains(clean_expr, "rolldn") || contains(clean_expr, "swap")
        warn(LOGGER, "parse_rpn does not support \"rollup\", \"rolldn\", or \"swap\", leaving as String")
        return expr
    end

    stack = []
    split_expr = contains(clean_expr, ",") ? split(clean_expr, ',') : split(clean_expr)

    for item in split_expr
        try
            if haskey(double_operators, item)
                b = pop!(stack)
                a = pop!(stack)
                push!(stack, double_operators[item](a, b))
            elseif haskey(single_operators, item)
                push!(stack, single_operators[item](pop!(stack)))
            else
                if item == "pi"
                    push!(stack, pi)
                else
                    push!(stack, parse(Float64, item))
                end
            end
        catch error
            if isa(error, ArgumentError)
                warn(LOGGER, "\"$expr\" is not valid Reverse Polish Notation, leaving as String")
                return expr
            end
        end
    end
    if length(stack) > 1
        warn(LOGGER, "\"$expr\" is not valid Reverse Polish Notation, leaving as String")
        return expr
    else
        return stack[1]
    end
end


"detects if `expr` is Reverse Polish Notation expression"
function isa_rpn(expr::AbstractString)
    for key in keys(merge(double_operators, single_operators))
        if contains(expr, key)
            return true
        end
    end
    return false
end


components = ["linecode", "linegeometry", "line", "linespacing", "loadshape",
              "growthshape", "tcc_curve", "wiredata", "xfmrcode", "vsource",
              "isource", "fault", "capacitor", "reactor", "transformer",
              "gictransformer", "gicline", "load", "generator", "indmach012",
              "storage", "capcontrol", "regcontrol", "energymeter", "monitor"]


"""
    get_prop_default(ctype)

Returns the default property values, or the expected Types if no default is
known, for a given component type `ctype`.
"""
function get_prop_default(ctype::String)::Dict

    line = Dict{String,Any}("length" => 1.0, "phases" => 3, "r1" => 0.058,
                            "x1" => 0.1206, "r0" => 0.1784, "x0" => 0.4047,
                            "c1" => 0.0, "c0" => 0.0,
                            "rmatrix" => "[0.09813333 |0.04013333 0.09813333 |0.04013333 0.04013333 0.09813333 ]",
                            "xmatrix" => "[0.2153 |0.0947 0.2153 |0.0947 0.0947 0.2153 ]",
                            "cmatrix" => "[0.0 | 0.0 0.0 | 0.0 0.0 0.0 ]",
                            "switch" => false, "rg" => 0.01805, "xg" => 0.155081,
                            "rho" => 100, "units" => "none", "earthmodel" => "Deri",
                            "b1" => 1.28177, "b0" => 0.6031858, "normamps" => 400.0,
                            "emergamps" => 600.0, "faultrate" => 0.1, "pctperm" => 20,
                            "repair" => 3, "basefreq" => 60, "enabled" => true)

    load = Dict{String, Any}("phases" => 3, "kv" => 12.47, "kw" => 10.0, "pf" => 0.88,
                             "model" => 1, "conn" => "wye", "kvar" => 5.39742822138087,
                             "rneut" => -1, "xneut" => 0, "status" => "variable",
                             "class" => 1, "vminpu" => 0.95, "vmaxpu" => 1.05,
                             "vminnorm" => 0.0, "vminemerg" => 0.0, "xfkva" => 0.0,
                             "allocationfactor" => 0.5, "kva" => 11.3636363636364,
                             "%mean" => 50, "%stddev" => 10, "cvrwatts" => 1,
                             "cvrvars" => 2, "kwh" => 0, "kwhdays" => 30,
                             "cfactor" => 4, "numcust" => 1, "%seriesrl" => 50,
                             "refweight" => 1, "vlowpu" => 0.5, "puxharm" => 0,
                             "xrharm" => 6, "spectrum" => "defaultload",
                             "basefreq" => 60, "enabled" => true)

    transformer = Dict{String,Any}("phases" => 3, "windings" => 2, "wdg" => 1, "conn" => "wye",
                                   "kv" => 12.47, "kva" => 1000.0, "tap" => 1.0, "%r" => 0.2,
                                   "rneut" => -1, "xneut" => 0, "conns" => ["wye", "wye"],
                                   "kvs" => [12.47, 12.47], "kvas" => [1000.0, 1000.0],
                                   "taps" => [1.0, 1.0], "xhl" => 7.0, "xht" => 35.0, "xlt" => 30.0,
                                   "xscarray" => [7], "thermal" => 2, "n" => 0.8, "m" => 0.8,
                                   "flrise" => 65, "hsrise" => 15, "%loadloss" => 0.4,
                                   "%noloadloss" => 0.0, "normhkva" => 1100, "emerghkva" => 1500,
                                   "sub" => "n", "maxtap" => 1.1, "mintap" => 0.9, "numtaps" => 32,
                                   "%imag" => 0.0, "ppm_antifloat" => 1, "%rs" => [0.2, 0.2],
                                   "xrconst" => "NO", "x12" => 7, "x13" => 35, "x23" => 30,
                                   "leadlag" => "Lag", "normamps" => 50.929, "emergamps" => 69.449,
                                   "faultrate" => 0.007, "pctperm" => 100, "repair" => 36,
                                   "basefreq" => 60, "enabled" => true)

    transformer3 = Dict{String,Any}()

    gen = Dict{String,Any}("phases" => 3, "kv" => 12.47, "kw" => 1000, "pf" => 0.88,
                           "kvar" => 60, "model" => 1, "vminpu" => 0.90, "vmaxpu" => 1.10,
                           "dispmode" => "Default", "dispvalue" => 0.0, "conn" => "wye",
                           "rneut" => 0, "xneut" => 0, "status" => "variable", "class" => 1,
                           "vpu" => 1.0, "maxkvar" => 120, "minkvar" => -120, "pvfactor" => 0.1,
                           "forceon" => "No", "kva" => 1200, "mva" => 1.2, "xd" => 1,
                           "xdpp" => 0.2, "h" => 1, "d" => 0, "debugtrace" => "no",
                           "balanced" => "No", "xrdp" => 20, "spectrum" => "defaultgen",
                           "basefreq" => 60, "enabled" => true)

    linecode = Dict{String,Any}("nphases" => 3, "r1" => 0.058, "x1" => 0.1206, "r0" => 0.1784,
                                "x0" => 0.4047, "c1" => 3.4, "c0" => 1.6, "units" => "none",
                                "rmatrix" => "[0.09813333 |0.04013333 0.09813333 |0.04013333 0.04013333 0.09813333 ]",
                                "xmatrix" => "[0.2153 |0.0947 0.2153 |0.0947 0.0947 0.2153 ]",
                                "cmatrix" => "[0.0 | 0.0 0.0 | 0.0 0.0 0.0 ]",
                                "basefreq" => 60, "normamps" => 400, "emergamps" => 600,
                                "faultrate" => 0.1, "pctperm" => 20, "repair" => 3,
                                "kron" => "N", "rg" => 0.01805, "xg" => 0.15508,
                                "neutral" => 3, "b1" => 1.2818, "b0" => 0.60319)

    capacitor = Dict{String,Any}("phases" => 3, "kvar" => [1200], "kv" => 12.47, "conn" => "wye",
                                 "cuf" => [20.47], "r" => [0], "xl" => [0], "harm" => [0],
                                 "numsteps" => 1, "states" => [1], "normamps" => 75.0046059412345,
                                 "emergamps" => 100.006141254979, "faultrate" => 0,
                                 "pctperm" => 1e2, "repair" => 3, "basefreq" => 60, "enabled" => true)

    reactor = Dict{String,Any}("phases" => 3, "kvar" => 1200.0, "kv" => 12.47, "conn" => "wye",
                               "parallel" => "NO", "r" => 0.0, "x" => 1555.009, "Rp" => 0.0,
                               "z1" => [0.0, 0.0], "z2" => [0.0, 0.0], "z0" => [0.0, 0.0],
                               "z" => [0.0, 1555.009], "lmh" => 4124.7895, "normamps" => 5.0,
                               "emergamps" => 6.0, "faultrate" => 0, "pctperm" => 1.0e2, "repair" => 3,
                               "basefreq" => 60.0, "enabled" => true)

    circuit = Dict{String,Any}("bus1" => "sourcebus", "basekv" => 115.0, "pu" => 1.0, "angle" => 0.0,
                               "frequency" => 60.0, "phases" => 3, "mvasc3" => 2000.0, "mvasc1" => 2100.0,
                               "x1r1" => 4.0, "x0r0" => 3.0, "lsc3" => 10041.0, "lsc1" => 10543.0,
                               "r1" => 1.6038, "x1" => 6.4151, "r0" => 1.796, "x0" => 5.3881, "scantype" => "Pos",
                               "sequence" => "Pos", "bus2" => "sourcebus.0.0.0", "z1" => [1.6037668, 6.4150673],
                               "z0" => [1.7960358, 5.3881075], "z2" => [1.6037668, 6.4150673],
                               "puz1" => [0.012126781, 0.048507125], "puz0" => [0.013580611, 0.040741834],
                               "puz2" => [0.012126781, 0.048507125], "basemva" => 100.0, "spectrum" => "defaultvsource",
                               "basefreq" => 60.0, "enabled" => true)

    vsource = Dict{String,Any}()

    ctypes = Dict{String, Dict}("line" => line,
                                "load" => load,
                                "linecode" => linecode,
                                "generator" => gen,
                                "capacitor" => capacitor,
                                "transformer" => transformer,
                                "transformer3" => transformer3,
                                "reactor" => reactor,
                                "circuit" => circuit,
                                "linecode" => linecode,
                                "vsource" => vsource)
    try
        return ctypes[ctype]
    catch KeyError
        return Dict{String,Any}()
    end
end


"""
    get_prop_name(ctype)

Returns the property names in order for a given component type `ctype`.
"""
function get_prop_name(ctype::AbstractString)::Array
    line = ["bus1", "bus2", "linecode", "length", "phases",
            "r1", "x1", "r0", "x0", "c1", "c0", "b1",
            "b0", "normamps", "emergamps", "faultrate", "pctperm",
            "repair", "basefreq", "rmatrix", "xmatrix", "cmatrix",
            "switch", "rg", "xg", "rho", "geometry", "earthmodel",
            "units", "like"]

    load = ["bus1", "phases", "kv", "kw", "pf", "model",
            "yearly", "daily", "duty", "growth", "conn", "kvar",
            "rneut", "xneut", "status", "class", "vminpu", "vmaxpu",
            "vminnorm", "vminemerg", "xfkva", "allocationfactor",
            "kva", "%mean", "%stddev", "cvrwatts", "cvrvars", "kwh",
            "kwhdays", "cfactor", "cvrcurve", "numcust", "spectrum",
            "zipv", "%seriesrl", "relweight", "vlowpu", "puxharm",
            "xrharm", "spectrum", "basefreq", "like"]

    gen = ["bus1", "phases", "kv", "kw", "pf", "model", "yearly",
           "daily", "duty", "dispvalue", "conn", "kvar", "rneut",
           "xneut", "status", "class", "maxkvar", "minkvar", "pvfactor",
           "debugtrace", "vminpu", "vmaxpu", "forceon", "kva", "mva",
           "xd", "xdp", "xdpp", "h", "d", "usermodel", "userdata",
           "shaftmodel", "shaftdata", "dutystart", "balanced", "xrdp",
           "spectrum", "basefreq", "like"]

    transformer = ["phases", "windings"]

    energymeter = ["element", "terminal", "action", "clear", "save",
                   "take", "option"]

    linecode = ["nphases", "r1", "x1", "r0", "x0", "c1", "c0", "units",
                "rmatrix", "xmatrix", "cmatrix", "basefreq", "normamps",
                "emergamps", "faultrate", "pctperm", "kron", "rg", "xg",
                "rho", "neutral", "b1", "b0", "like"]

    capacitor = ["bus1", "bus2", "phases", "kvar", "kv", "conn", "cmatrix",
                 "cuf", "r", "xl", "harm", "numsteps", "states", "normamps",
                 "emergamps", "faultrate", "pctperm", "basefreq", "like"]

    ctypes = Dict{String, Array}("line" => line,
                                 "load" => load,
                                 "transformer" => transformer,
                                 "energymeter" => energymeter,
                                 "linecode" => linecode,
                                 "generator" => gen,
                                 "capacitor" => capacitor)

    try
        return ctypes[ctype]
    catch KeyError
        return []
    end
end


"""
    get_prop_name(ctype, i)

Returns the `i`th property name for a given component type `ctype`.
"""
function get_prop_name(ctype::AbstractString, i::Int)::String
    return get_prop_name(ctype)[i]
end


""
function get_linecode(dss_data::Dict, id::AbstractString)
    if haskey(dss_data, "linecode")
        for item in dss_data["linecode"]
            if item["id"] == id
                return item
            end
        end
    end
    return Dict{String,Any}()
end


"""
    parse_matrix(dtype, data)

Parses a OpenDSS style triangular matrix string `data` into a two dimensional
array of type `dtype`. Matrix strings are capped by either parenthesis or
brackets, rows are separated by "|", and columns are separated by spaces.
"""
function parse_matrix(dtype::Type, data::AbstractString, nphases::Int=3)::Array
    rows = []
    for line in split(strip(data, ['[', ']', '(', ')']), '|')
        cols = []
        for item in split(line)
            push!(cols, parse(dtype, item))
        end
        push!(rows, cols)
    end

    matrix = zeros(dtype, nphases, nphases)
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


"""
    parse_array(dtype, data)

Parses a OpenDSS style array string `data` into a one dimensional array of type
`dtype`. Array strings are capped by either brackets, single quotes, or double
quotes, and elements are separated by spaces.
"""
function parse_array(dtype::Type, data::AbstractString)::Array
    if contains(data, ",")
        split_char = ','
    else
        split_char = ' '
    end

    elements = split(strip(data, ['\"', '\'', '[', ']', '(', ')']), split_char)
    elements = [strip(el) for el in elements if strip(el) != ""]

    if dtype == String
        array = []
        for el in elements
            push!(array, el)
        end
    else
        array = zeros(dtype, length(elements))
        for (i, el) in enumerate(elements)
            array[i] = parse(dtype, el)
        end
    end

    return array
end


"strips lines that are either commented (block or single) or empty"
function strip_lines(lines::Array)::Array
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
    parse_buscoords(file)

Parses a Bus Coordinate `file`, in either "dat" or "csv" formats, where in
"dat", columns are separated by spaces, and in "csv" by commas. File expected
to contain "bus,x,y" on each line.
"""
function parse_buscoords(file::AbstractString)::Array
    file_str = readstring(open(file))
    regex = r",\s*"
    if endswith(lowercase(file), "csv") || endswith(lowercase(file), "dss")
        regex = r","
    end

    lines = strip_lines(split(file_str, '\n'))

    coordArray = []
    for line in lines
        bus, x, y = split(line, regex; limit=3)
        push!(coordArray, Dict{String,Any}("bus"=>strip(bus, [',']),
                                            "id"=>strip(bus, [',']),
                                            "x"=>parse(Float64, strip(x, [','])),
                                            "y"=>parse(Float64, strip(y, [',', '\r']))))
    end
    return coordArray
end


"""
    parse_properties(properties)

Parses a string of `properties` of a component type, character by character
into an array with each element containing (if present) the property name, "=",
and the property value.
"""
function parse_properties(properties::AbstractString)::Array
    propsOut = []
    endArray = true
    endProp = false
    endEquality = false
    endSQuot = true
    endDQuot = true
    str_out = ""

    properties = replace(properties, r"\s*=\s*", "=")
    nchars = length(properties)

    for (n, char) in enumerate(properties)
        sstr_out = split(str_out, "=")
        if length(sstr_out) == 2 && sstr_out[2] != ""
            endEquality = true
        elseif !contains(str_out, "=") && (char == ' ' || n == nchars)
            endEquality = true
        else
            endEquality = false
        end

        if char == ' '
            endProp = true
        else
            endProp = false
        end

        if char in ['[', '(']
            endArray = false
        elseif char in [']', ')']
            endArray = true
        end

        if char == '\"'
            endDQuot = !endDQuot
        elseif char == '\''
            endSQuot = !endSQuot
        end

        if char != ' ' || !endArray || !endDQuot || !endSQuot
            str_out = string(str_out, char)
        end

        if str_out != "" && endArray && endProp && endEquality && endDQuot && endSQuot || n == nchars
            push!(propsOut, str_out)
            str_out = ""
        end
    end

    return propsOut
end


"""
    add_component!(dss_data, ctype_name, compDict)

Adds a component of type `ctype_name` with properties given by `compDict` to
the existing `dss_data` structure. If a component of the same type has already
been added to `dss_data`, the new component is appeneded to the existing array
of components of that type, otherwise a new array is created.
"""
function add_component!(dss_data::Dict, ctype_name::AbstractString, compDict::Dict)
    ctype = split(lowercase(ctype_name), '.'; limit=2)[1]
    if haskey(dss_data, ctype)
        push!(dss_data[ctype], compDict)
    else
        dss_data[ctype] = [compDict]
    end
end


"""
    add_property(compDict, key, value)

Adds a property to an existing component properties dictionary `compDict` given
the `key` and `value` of the property. If a property of the same name already
exists inside `compDict`, the `key` is renamed to have `_\d` appended to the
end.
"""
function add_property(compDict::Dict, key::AbstractString, value::Any)::Dict
    if haskey(compDict, lowercase(key))
        rmatch = match(r"_(\d+)$", key)
        if typeof(rmatch) != Void
            endNum = parse(Int, rmatch.captures[1]) + 1
            key = replace(key, r"_(\d+)$", "_$endNum")
        else
            key = string(key, "_2")
        end
    end

    compDict[lowercase(key)] = value

    return compDict
end


"""
    parse_component(component, properies, compDict=Dict{String,Any}())

Parses a `component` with `properties` into a `compDict`. If `compDict` is not
defined, an empty dictionary will be used. Assumes that unnamed properties are
given in order, but named properties can be given anywhere.
"""
function parse_component(component::AbstractString, properties::AbstractString, compDict::Dict=Dict{String,Any}())
    debug(LOGGER, "Properties: $properties")
    ctype, name = split(lowercase(component), '.'; limit=2)

    if !haskey(compDict, "id")
        compDict["id"] = name
    end

    propArray = parse_properties(properties)
    debug(LOGGER, "propArray: $propArray")

    propNames = get_prop_name(lowercase(ctype))

    for (n, property) in enumerate(propArray)
        if property == ""
            continue
        elseif !contains(property, "=")
            property = join([shift!(propNames), property], '=')
        else
            filter!(e->e!=lowercase(split(property,'=')[1]), propNames)
        end

        key, value = split(property, '='; limit=2)

        add_property(compDict, key, value)
    end

    return compDict
end


"""
    merge_dss!(dss_prime, dss_to_add)

Merges two (partially) parsed OpenDSS files to the same dictionary `dss_prime`.
Used in cases where files are referenced via the "compile" or "redirect"
OpenDSS commands inside the originating file.
"""
function merge_dss!(dss_prime::Dict{String,Array}, dss_to_add::Dict{String,Array})
    for (k, v) in dss_to_add
        if k in keys(dss_prime)
            append!(dss_prime[k], v)
        else
            dss_prime[k] = v
        end
    end
end


"""
    parse_line(elements, curCompDict=Dict{String,Any}())

Parses an already separated line given by `elements` (an array) of an OpenDSS
file into `curCompDict`. If not defined, `curCompDict` is an empty dictionary.
"""
function parse_line(elements::Array, curCompDict::Dict=Dict{String,Any}())
    curCtypeName = strip(elements[2], ['\"', '\''])
    if startswith(lowercase(curCtypeName), "object")
        curCtypeName = split(curCtypeName, '=')[2]
        curCompDict["id"] = split(curCtypeName, '.')[2]
    else
        if length(elements) != 3
            properties = ""
        else
            properties = elements[3]
        end

        curCompDict = parse_component(curCtypeName, properties)
    end

    return curCtypeName, curCompDict
end


"Strips comments, defined by \"!\" from the ends of lines"
function strip_comments(line::AbstractString)::String
    return strip(split(line, r"\s*!")[1], ['\r', '\n'])
end


"""
    assign_property!(dss_data, cType, cName, propName, propValue)

Assigns a property with name `propName` and value `propValue` to the component
of type `cType` named `cName` in `dss_data`.
"""
function assign_property!(dss_data::Dict, cType::AbstractString, cName::AbstractString,
                          propName::AbstractString, propValue::Any)
    if haskey(dss_data, cType)
        for obj in dss_data[cType]
            if lowercase(obj["id"]) == cName
                obj[propName] = propValue
            end
        end
    else
        warn(LOGGER, "Cannot find $cType object $cName.")
    end
end


"""
    parse_dss(filename)

Parses a OpenDSS file given by `filename` into a Dict{Array{Dict}}. Only
supports components and options, but not commands, e.g. "plot" or "solve".
Will also parse files defined inside of the originating DSS file via the
"compile", "redirect" or "buscoords" commands.
"""
function parse_dss(filename::AbstractString)::Dict
    info(LOGGER, "Calling parse_dss on $filename")
    currentFile = split(filename, "/")[end]
    path = join(split(filename, '/')[1:end-1], '/')
    dss_str = readstring(open(filename))
    dss_data = Dict{String,Array}()

    dss_data["filename"] = [currentFile]

    curCompDict = Dict{String,Any}()
    curCtypeName = ""

    lines = split(dss_str, '\n')

    stripped_lines = strip_lines(lines)
    nlines = length(stripped_lines)

    for (n, line) in enumerate(stripped_lines)
        real_line_num = find(lines .== line)[1]
        debug(LOGGER, "LINE $real_line_num: $line")
        line = strip_comments(line)

        if contains(line, "{") || contains(line, "}")
            warn(LOGGER, "Line $real_line_num in \"$currentFile\" contains an unsupported symbol, skipping")
            continue
        elseif startswith(strip(line), '~')
            curCompDict = parse_component(curCtypeName, strip(strip(line, '~')), curCompDict)

            if n < nlines && startswith(strip(stripped_lines[n + 1]), '~')
                continue
            else
                add_component!(dss_data, curCtypeName, curCompDict)
            end
        else
            curCompDict = Dict{String,Any}()
            line_elements = split(line, r"\s+"; limit=3)
            cmd = lowercase(line_elements[1])

            if cmd == "clear"
                info(LOGGER, "`dss_data` has been reset with the \"clear\" command.")
                dss_data = Dict{String,Array}()
                continue

            elseif cmd == "redirect"
                file = line_elements[2]
                fullpath = path == "" ? file : join([path, file], '/')
                info(LOGGER, "Redirecting to file \"$file\"")
                merge_dss!(dss_data, parse_dss(fullpath))
                continue

            elseif cmd == "compile"
                file = split(strip(line_elements[2], ['(',')']), '\\')[end]
                fullpath = path == "" ? file : join([path, file], '/')
                info(LOGGER, "Compiling file \"$file\"")
                merge_dss!(dss_data, parse_dss(fullpath))
                continue

            elseif cmd == "set"
                debug(LOGGER, "set command: $line_elements")
                if length(line_elements) == 2
                    property, value = split(line_elements[2], '='; limit=2)
                else
                    property, value = line_elements[2], strip(strip(line_elements[3], '='))
                end

                if !haskey(dss_data, "options")
                    dss_data["options"] = [Dict{String,Any}()]
                end

                dss_data["options"][1]["$(lowercase(property))"] = value
                continue

            elseif cmd == "buscoords"
                file = line_elements[2]
                fullpath = path == "" ? file : join([path, file], '/')
                debug(LOGGER, "Buscoords path: $fullpath")
                dss_data["buscoords"] = parse_buscoords(fullpath)

            elseif cmd == "new"
                curCtypeName, curCompDict = parse_line(line_elements)
            else
                try
                    cType, cName, props = split(lowercase(line), '.'; limit=3)
                    propsOut = parse_properties(props)
                    for prop in propsOut
                        propName, propValue = split(prop, '=')
                        assign_property!(dss_data, cType, cName, propName, propValue)
                    end
                catch
                    warn(LOGGER, "Command \"$cmd\" on line $real_line_num in \"$currentFile\" is not supported, skipping.")
                end
            end

            debug(LOGGER, "size curCompDict: $(length(curCompDict))")
            if n < nlines && startswith(strip(stripped_lines[n + 1]), '~')
                continue
            elseif length(curCompDict) > 0
                add_component!(dss_data, curCtypeName, curCompDict)
            else
                continue
            end
        end
    end

    info(LOGGER, "Done parsing $filename")
    return dss_data
end
