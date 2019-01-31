"Squares `x`, for parsing Reverse Polish Notation"
function sqr(x::Float64)
    return x * x
end


double_operators = Dict("+" => +, "-" => -, "*" => *, "/" => /, "^" => ^,
                        "atan2" => (x, y) -> rad2deg(atan2(y, x)))

single_operators = Dict("sqr" => sqr, "sqrt" => sqrt, "inv" => inv, "ln" => log,
                        "exp" => exp, "log10" => log10, "sin" => sind, "cos" => cosd,
                        "tan" => tand, "asin" => asind, "acos" => acosd, "atan" => atand)

array_delimiters = ['\"', '\'', '[', '{', '(', ']', '}', ')']

"parses Reverse Polish Notation `expr`"
function parse_rpn(expr::AbstractString, dtype::Type=Float64)
    clean_expr = strip(expr, array_delimiters)

    if occursin("rollup", clean_expr) || occursin("rolldn", clean_expr) || occursin("swap", clean_expr)
        Memento.warn(LOGGER, "parse_rpn does not support \"rollup\", \"rolldn\", or \"swap\", leaving as String")
        return expr
    end

    stack = []
    split_expr = occursin(",", clean_expr) ? split(clean_expr, ',') : split(clean_expr)

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
                    push!(stack, parse(dtype, item))
                end
            end
        catch error
            if isa(error, ArgumentError)
                Memento.warn(LOGGER, "\"$expr\" is not valid Reverse Polish Notation, leaving as String")
                return expr
            end
        end
    end
    if length(stack) > 1
        Memento.warn(LOGGER, "\"$expr\" is not valid Reverse Polish Notation, leaving as String")
        return expr
    else
        return stack[1]
    end
end


"detects if `expr` is Reverse Polish Notation expression"
function isa_rpn(expr::AbstractString)::Bool
    expr = split(strip(expr, array_delimiters))
    opkeys = keys(merge(double_operators, single_operators))
    for item in expr
        if item in opkeys
            return true
        end
    end
    return false
end


"parses connection \"conn\" specification reducing to wye or delta"
function parse_conn(conn::String)::String
    if conn in ["wye", "y", "ln"]
        return "wye"
    elseif conn in ["delta", "ll"]
        return "delta"
    else
        Memento.warn(LOGGER, "Unsupported connection $conn, defaulting to \"wye\"")
        return "wye"
    end
end


"checks is a string is a connection by checking the values"
function isa_conn(expr::AbstractString)::Bool
    if expr in ["wye", "y", "ln", "delta", "ll"]
        return true
    else
        return false
    end
end


"""
    get_prop_name(ctype)

Returns the property names in order for a given component type `ctype`.
"""
function get_prop_name(ctype::AbstractString)::Array
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

    ctypes = Dict{String, Array}("linecode" => linecode,
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

    try
        return ctypes[ctype]
    catch e
        throw(e)
    end
end


"""
    get_prop_name(ctype, i)

Returns the `i`th property name for a given component type `ctype`.
"""
function get_prop_name(ctype::AbstractString, i::Int)::String
    return get_prop_name(ctype)[i]
end


"""
    parse_matrix(dtype, data)

Parses a OpenDSS style triangular matrix string `data` into a two dimensional
array of type `dtype`. Matrix strings are capped by either parenthesis or
brackets, rows are separated by "|", and columns are separated by spaces.
"""
function parse_matrix(dtype::Type, data::AbstractString)::Array
    rows = []
    for line in split(strip(data, array_delimiters), '|')
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
function parse_matrix(data::Array{T}, nodes::Array{Bool}, nph::Int=3, fill_val=0.0)::Array where T
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
function isa_matrix(data::AbstractString)::Bool
    if occursin("|", data)
        return true
    else
        return false
    end
end


"""
    parse_array(dtype, data)

Parses a OpenDSS style array string `data` into a one dimensional array of type
`dtype`. Array strings are capped by either brackets, single quotes, or double
quotes, and elements are separated by spaces.
"""
function parse_array(dtype::Type, data::AbstractString)
    if occursin(",", data)
        split_char = ','
    else
        split_char = ' '
    end

    if isa_rpn(data)
        matches = collect((m.match for m = eachmatch(Regex(string("[",join(array_delimiters, '\\'),"]")), data, overlap=false)))
        if length(matches) == 2
            if dtype == String
                return data
            else
                return parse_rpn(data, dtype)
            end

        else
            elements = parse_properties(data[2:end-1])
        end
    else
        elements = split(strip(data, array_delimiters), split_char)
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
            if isa_rpn(data)
                array[i] = parse_rpn(el, dtype)
            else
                array[i] = parse(dtype, el)
            end
        end
    end

    return array
end


"parse matrices according to active nodes"
function parse_array(data, nodes::Array{Bool}, nph::Int=3, fill_val=0.0)::Array
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
function isa_array(data::AbstractString)::Bool
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
        elseif startswith(clean_data, "{") && endswith(clean_data, ")")
            return true
        else
            return false
        end
    else
        return false
    end
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
    file_str = read(open(file), String)
    regex = r",\s*"
    if endswith(lowercase(file), "csv") || endswith(lowercase(file), "dss")
        regex = r","
    end

    lines = strip_lines(split(file_str, '\n'))

    coordArray = []
    for line in lines
        bus, x, y = split(line, regex; limit=3)
        push!(coordArray, Dict{String,Any}("bus"=>strip(bus, [',']),
                                            "name"=>strip(bus, [',']),
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

    properties = replace(properties, r"\s*=\s*" => "=")
    nchars = length(properties)

    for (n, char) in enumerate(properties)
        sstr_out = split(str_out, "=")
        if length(sstr_out) == 2 && sstr_out[2] != ""
            endEquality = true
        elseif !occursin("=", str_out) && (char == ' ' || n == nchars)
            endEquality = true
        else
            endEquality = false
        end

        if char == ' '
            endProp = true
        else
            endProp = false
        end

        if char in ['[', '(', '{']
            endArray = false
        elseif char in [']', ')', '}']
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
    ctype = split(ctype_name, '.'; limit=2)[1]
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
exists inside `compDict`, the original value is converted to an array, and the
new value is appended to the end.
"""
function add_property(compDict::Dict, key::AbstractString, value::Any)::Dict
    if haskey(compDict, lowercase(key))
        rmatch = match(r"_(\d+)$", key)
        if typeof(rmatch) != Nothing
            endNum = parse(Int, rmatch.captures[1]) + 1
            key = replace(key, r"_(\d+)$" => "_$endNum")
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
    Memento.debug(LOGGER, "Properties: $properties")
    ctype, name = split(component, '.'; limit=2)

    if !haskey(compDict, "name")
        compDict["name"] = name
    end

    propArray = parse_properties(properties)
    Memento.debug(LOGGER, "propArray: $propArray")

    propNames = get_prop_name(ctype)
    propIdx = 1

    for (n, property) in enumerate(propArray)
        if property == ""
            continue
        elseif !occursin("=", property)
            property = join([propNames[propIdx], property], '=')
            propIdx += 1
        else
            if split(component,'.')[1] == "loadshape" && startswith(property, "mult")
                property = replace(property, "mult" => "pmult")
            end

            try
                propIdxs = findall(e->e==split(property,'=')[1], propNames)
                if length(propIdxs) > 0
                    propIdx = findall(e->e==split(property,'=')[1], propNames)[1] + 1
                end
            catch e
                if split(component,'.')[1] == "transformer"
                    if split(property,'=')[1] == "ppm"
                        property = replace(property, "ppm" => "ppm_antifloat")
                    elseif split(property,'=')[1] == "x12"
                        property = replace(property, "x12" => "xhl")
                    elseif split(property,'=')[1] == "x23"
                        property = replace(property, "x23" => "xlt")
                    elseif split(property,'=')[1] == "x13"
                        property = replace(property, "x13" => "xht")
                    end
                else
                    throw(e)
                end
                propIdx = findall(e->e==split(property,'=')[1], propNames)[1] + 1
            end
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
    if startswith(curCtypeName, "object")
        curCtypeName = split(curCtypeName, '=')[2]
        curCompDict["name"] = split(curCtypeName, '.')[2]
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
            if obj["name"] == cName
                obj[propName] = propValue
            end
        end
    else
        Memento.warn(LOGGER, "Cannot find $cType object $cName.")
    end
end


"""
    parse_dss(filename)

Parses a OpenDSS file given by `filename` into a Dict{Array{Dict}}. Only
supports components and options, but not commands, e.g. "plot" or "solve".
Will also parse files defined inside of the originating DSS file via the
"compile", "redirect" or "buscoords" commands.
"""
function parse_dss(io::IOStream)::Dict
    filename = match(r"^<file\s(.+)>$", io.name).captures[1]
    Memento.info(LOGGER, "Calling parse_dss on $filename")
    currentFile = split(filename, "/")[end]
    path = join(split(filename, '/')[1:end-1], '/')
    dss_data = Dict{String,Array}()

    dss_data["filename"] = [currentFile]

    curCompDict = Dict{String,Any}()
    curCtypeName = ""

    lines = readlines(io)

    stripped_lines = strip_lines(lines)
    nlines = length(stripped_lines)

    for (n, line) in enumerate(stripped_lines)
        real_line_num = findall(lines .== line)[1]
        Memento.debug(LOGGER, "LINE $real_line_num: $line")
        line = strip_comments(line)

        if startswith(strip(line), '~')
            curCompDict = parse_component(curCtypeName, strip(strip(lowercase(line)),  '~'), curCompDict)

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
                Memento.info(LOGGER, "`dss_data` has been reset with the \"clear\" command.")
                dss_data = Dict{String,Array}("filename"=>dss_data["filename"])
                continue

            elseif cmd == "redirect"
                file = line_elements[2]
                fullpath = path == "" ? file : join([path, file], '/')
                Memento.info(LOGGER, "Redirecting to file \"$file\"")
                merge_dss!(dss_data, parse_dss(fullpath))
                continue

            elseif cmd == "compile"
                file = split(strip(line_elements[2], ['(',')']), '\\')[end]
                fullpath = path == "" ? file : join([path, file], '/')
                Memento.info(LOGGER, "Compiling file \"$file\"")
                merge_dss!(dss_data, parse_dss(fullpath))
                continue

            elseif cmd == "set"
                Memento.debug(LOGGER, "set command: $line_elements")
                if length(line_elements) == 2
                    property, value = split(lowercase(line_elements[2]), '='; limit=2)
                else
                    property, value = lowercase(line_elements[2]), strip(strip(lowercase(line_elements[3]), '='))
                end

                if !haskey(dss_data, "options")
                    dss_data["options"] = [Dict{String,Any}()]
                end

                dss_data["options"][1]["$(property)"] = value
                continue

            elseif cmd == "buscoords"
                file = line_elements[2]
                fullpath = path == "" ? file : join([path, file], '/')
                Memento.debug(LOGGER, "Buscoords path: $fullpath")
                dss_data["buscoords"] = parse_buscoords(fullpath)

            elseif cmd == "new"
                curCtypeName, curCompDict = parse_line([lowercase(line_element) for line_element in line_elements])
            else
                try
                    cType, cName, props = split(lowercase(line), '.'; limit=3)
                    propsOut = parse_properties(props)
                    for prop in propsOut
                        propName, propValue = split(prop, '=')
                        assign_property!(dss_data, cType, cName, propName, propValue)
                    end
                catch
                    Memento.warn(LOGGER, "Command \"$cmd\" on line $real_line_num in \"$currentFile\" is not supported, skipping.")
                end
            end

            Memento.debug(LOGGER, "size curCompDict: $(length(curCompDict))")
            if n < nlines && startswith(strip(stripped_lines[n + 1]), '~')
                continue
            elseif length(curCompDict) > 0
                add_component!(dss_data, curCtypeName, curCompDict)
            else
                continue
            end
        end
    end

    Memento.info(LOGGER, "Done parsing $filename")
    return dss_data
end


""
function parse_dss(filename::AbstractString)::Dict
    dss_data = open(filename) do io
        parse_dss(io)
    end
    return dss_data
end


"parses the raw dss values into their expected data types"
function parse_element_with_dtype(dtype, element)
    if isa_matrix(element)
        out = parse_matrix(eltype(dtype), element)
    elseif isa_array(element)
        out = parse_array(eltype(dtype), element)
    elseif dtype <: Bool
        if element in ["n", "no"]
            element = "false"
        elseif element in ["y", "yes"]
            element = "true"
        end
        out = parse(dtype, element)
    elseif isa_rpn(element)
        out = parse_rpn(element)
    elseif dtype == String
        out = element
    else
        if isa_conn(element)
            out = parse_conn(element)
        else
            try
                out = parse(dtype, element)
            catch
                Memento.warn(LOGGER, "cannot parse $element as $dtype, leaving as String.")
                out = element
            end
        end
    end

    return out
end


"""
    parse_dss_with_dtypes!(dss_data, toParse)

Parses the data in keys defined by `toParse` in `dss_data` using types given by
the default properties from the `get_prop_default` function.
"""
function parse_dss_with_dtypes!(dss_data::Dict, toParse::Array{String}=[])
    for compType in toParse
        if haskey(dss_data, compType)
            Memento.debug(LOGGER, "type: $compType")
            for item in dss_data[compType]
                dtypes = get_dtypes(compType)
                for (k, v) in item
                    if haskey(dtypes, k)
                        Memento.debug(LOGGER, "key: $k")
                        if isa(v, Array)
                            arrout = []
                            for el in v
                                if isa(v, AbstractString)
                                    push!(arrout, parse_element_with_dtype(dtypes[k], el))
                                else
                                    push!(arrout, el)
                                end
                            end
                            item[k] = arrout
                        elseif isa(v, AbstractString)
                            item[k] = parse_element_with_dtype(dtypes[k], v)
                        else
                            Memento.error(LOGGER, "dtype unknown $compType, $k, $v")
                        end
                    end
                end
            end
        end
    end
end


"""
    parse_busname(busname)

Parses busnames as defined in OpenDSS, e.g. "primary.1.2.3.0".
"""
function parse_busname(busname::AbstractString)
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


"converts Dict{String,Any} to Dict{Symbol,Any} for passing as kwargs"
function to_sym_keys(data::Dict{String,Any})::Dict{Symbol,Any}
    return Dict{Symbol,Any}((Symbol(k), v) for (k, v) in data)
end
