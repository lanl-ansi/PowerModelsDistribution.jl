# OpenDSS parser


components = ["linecode", "linegeometry", "line", "linespacing", "loadshape",
              "growthshape", "tcc_curve", "wiredata", "xfmrcode", "vsource",
              "isource", "fault", "capacitor", "reactor", "transformer",
              "gictransformer", "gicline", "load", "generator", "indmach012",
              "storage", "capcontrol", "regcontrol", "energymeter", "monitor"]

function warn_get(data, fname, default; valid = x -> true)
    if haskey(data, fname)
        value = data[fname]
    else
        warn("field $fname isn't present, setting to $default")
        return default
    end

    if isa(value, typeof(default)) || isa(default, String)
        x = value
    elseif isa(value, String)
        x = parse(typeof(default), value)
    else
        warn("Invalid type $(typeof(value)) of $fname, setting to $default")
        return default
    end

    if valid(x)
        return x
    end 

    warn("Invalid value $x of $fname, setting to $default")
    return default
end


function warn_pop!(data, fname, default; dtype=nothing)
    value = default

    if haskey(data, fname)
        value = pop!(data, fname)
    else
        warn("field $fname isn't present, setting to $default")
    end

    if dtype === nothing
        return value
    end

    return parse(dtype, value)
end


"""
    get_prop_default(ctype)

Returns the default property values, or the expected Types if no default is
known, for a given component type `ctype`.
"""
function get_prop_default(ctype::String)::Dict

    line = Dict{String,Any}("length" => 1.0, "phases" => 3, "r1" => 0.058,
                            "x1" => 0.1206, "r0" => 0.1784, "x0" => 0.4047,
                            "c1" => 3.4, "c0" => 1.6,
                            "rmatrix" => "[0.09813333 |0.04013333 0.09813333 |0.04013333 0.04013333 0.09813333 ]",
                            "xmatrix" => "[0.2153 |0.0947 0.2153 |0.0947 0.0947 0.2153 ]",
                            "cmatrix" => "[2.8 |-0.6 2.8 |-0.6 -0.6 2.8 ]",
                            "switch" => false, "rg" => 0.01805, "xg" => 0.155081,
                            "rho" => 100, "units" => "none", "earthmodel" => "Deri",
                            "b1" => 1.28177, "b0" => 0.6031858, "normamps" => 400,
                            "emergamps" => 600, "faultrate" => 0.1, "pctperm" => 20,
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

    transformer2 = Dict{String,Any}("phases" => 3, "windings" => 2, "wdg" => 1, "conn" => "wye",
                                    "kv" => 12.47, "kva" => 1000, "tap" => 1, "%r" => 0.2,
                                    "rneut" => -1, "xneut" => 0, "conns" => ["wye", "wye"],
                                    "kvs" => [12.47, 12.47], "kvas" => [1000, 1000],
                                    "taps" => [1, 1], "xhl" => 7, "xht" => 35, "xlt" => 30,
                                    "xscarray" => [7], "thermal" => 2, "n" => 0.8, "m" => 0.8,
                                    "flrise" => 65, "hsrise" => 15, "%loadloss" => 0.4,
                                    "%noloadloss" => 0, "normhkva" => 1100, "emerghkva" => 1500,
                                    "sub" => "n", "maxtap" => 1.1, "mintap" => 0.9, "numtaps" => 32,
                                    "%imag" => 0, "ppm_antifloat" => 1, "%rs" => [0.2, 0.2],
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
                                "cmatrix" => "[2.8 |-0.6 2.8 |-0.6 -0.6 2.8 ]",
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

    linecode = Dict{String,Any}()

    ctypes = Dict{String, Dict}("line" => line,
                                "load" => load,
                                "linecode" => linecode,
                                "generator" => gen,
                                "capacitor" => capacitor,
                                "transformer2" => transformer2,
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


"""
    parse_matrix(dtype, data)

Parses a OpenDSS style triangular matrix string `data` into a two dimensional
array of type `dtype`. Matrix strings are capped by either parenthesis or
brackets, rows are separated by "|", and columns are separated by spaces.
"""
function parse_matrix(dtype::Type, data::AbstractString)::Array
    rows = []
    for line in split(strip(data, ['[', ']', '(', ')']), '|')
        cols = []
        for item in split(line)
            push!(cols, parse(dtype, item))
        end
        push!(rows, cols)
    end

    matrix = zeros(dtype, length(rows), length(rows[end]))
    for (i, row) in enumerate(rows)
        for (j, col) in enumerate(row)
            matrix[i, j] = matrix[j, i] = col
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
    elements = split(strip(data, ['\"', '\'', '[', ']', '(', ')']))
    array = zeros(dtype, length(elements))

    for (i, el) in enumerate(elements)
        array[i] = parse(dtype, el)
    end

    return array
end


"""
    parse_buscoords(file)

Parses a Bus Coordinate `file`, in either "dat" or "csv" formats, where in
"dat", columns are separated by spaces, and in "csv" by commas. File expected
to contain "bus,x,y" on each line.
"""
function parse_buscoords(file::AbstractString)::Array
    file_str = readstring(open(file))
    regex = r"\s+"
    if endswith(lowercase(file), "csv")
        regex = r","
    end

    coordArray = []
    for line in split(file_str, '\n')
        if line != ""
            bus, x, y = split(line, regex; limit=3)
            push!(coordArray, Dict{String,Any}("bus"=>bus, "x"=>x, "y"=>y))
        end
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
    debug(LOGGER, "add_component! $ctype_name")
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
    curCtypeName = elements[2]
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
    return split(line, r"\s*!")[1]
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
    currentFile = split(filename, "/")[end]
    path = join(split(filename, '/')[1:end-1], '/')
    dss_str = readstring(open(filename))
    dss_data = Dict{String,Array}()

    dss_data["filename"] = [currentFile]

    curCompDict = Dict{String,Any}()
    curCtypeName = ""

    lines = split(dss_str, '\n')

    blockComment = false
    stripped_lines = []
    for line in lines
        if startswith(strip(line), "/*") || endswith(strip(line), "*/")
            blockComment = !blockComment
        elseif !startswith(line, '!') && !startswith(line, "//") && strip(line) != "" && !blockComment
            push!(stripped_lines, line)
        end
    end

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

                dss_data["options"][1]["$property"] = value
                continue

            elseif cmd == "buscoords"
                file = line_elements[2]
                fullpath = path == "" ? file : join([path, file], '/')
                dss_data["buscoords"] = parse_buscoords(fullpath)

            elseif cmd == "new"
                curCtypeName, curCompDict = parse_line(line_elements)
            else
                try
                    cType, cName, prop = split(lowercase(line), '.'; limit=3)
                    propName, propValue = split(prop, '=')
                    assign_property!(dss_data, cType, cName, propName, propValue)
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

    return dss_data
end


"""
"""
function parse_dss_with_dtypes!(dss_data::Dict, toParse::Array{String}=[])
    for compType in toParse
        if haskey(dss_data, compType)
            for item in dss_data[compType]
                defaults = get_prop_default(compType)
                for (k, v) in item
                    if haskey(defaults, k)
                        if isa(defaults[k], Array)
                            try
                                item[k] = parse_array(Float64, v)
                            catch err
                                debug(LOGGER, "$compType $k")
                                error(LOGGER, err)
                            end
                        elseif isa(v, AbstractString) && ~isa(defaults[k], String)
                            try
                                item[k] = parse(typeof(defaults[k]), v)
                            catch err
                                debug(LOGGER, "$compType $k")
                                error(LOGGER, err)
                            end
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
    name, elements = split(busname,'.'; limit=2)
    nodes = Array{Bool}([0 0 0 0])

    for num in range(1,3)
        if contains(elements, "$num")
            nodes[num] = true
        end
    end

    if contains(elements, "0") || sum(nodes[1:3]) == 1
        nodes[4] = true
    end

    return name, nodes
end


"""
    discover_buses(dss_data)

Discovers all of the buses (not separately defined in OpenDSS), from "lines".
"""
function discover_buses(dss_data::Dict)::Array
    bus_names = []
    buses = []
    if haskey(dss_data, "line")
        for branch in dss_data["line"]
            for key in ["bus1", "bus2"]
                name, nodes = parse_busname(branch[key])
                if !(name in bus_names)
                    push!(bus_names, name)
                    push!(buses, (name, nodes))
                end
            end
        end
    else
        error(LOGGER, "dss_data has no branches!")
    end
    return buses
end


"""
    dss2tppm_bus!(tppm_data, dss_data)

Adds PowerModels-style buses to `tppm_data` from `dss_data`.
"""
function dss2tppm_bus!(tppm_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(tppm_data, "bus")
        tppm_data["bus"] = []
    end

    buses = discover_buses(dss_data)

    # TODO: merge with defaults, overriding where necessary
    # TODO: handle "sourcebus" (created with the circuit object)
    for (n, (bus, nodes)) in enumerate(buses)
        busDict = Dict{String,Any}()

        nphases = sum(nodes[1:3])

        busDict["bus_i"] = n  # TODO: how to properly index?
        busDict["index"] = n  # TODO: how to properly index?
        busDict["name"] = bus

        busDict["vm"] = ones(nphases)  # TODO: zeros(numPhases)?
        busDict["va"] = zeros(nphases)  # TODO: zeros(numPhases)?

        busDict["vmin"] = 0.9 * ones(nphases)  # TODO: implicit limits? 10%
        busDict["vmax"] = 1.1 * ones(nphases)  # TODO: implicit limits?

        busDict["base_kv"] = NaN  # TODO:

        push!(tppm_data["bus"], busDict)
    end
end


"""
    find_bus(busname, tppm_data)

Finds the index number of the bus in existing data from the given `busname`.
"""
function find_bus(busname::AbstractString, tppm_data::Dict)::Int
    for bus in tppm_data["bus"]
        if bus["name"] == busname
            return bus["bus_i"]
        else
            return 0
        end
    end
end


"""
    dss2tppm_load!(tppm_data, dss_data)

Adds PowerModels-style loads to `tppm_data` from `dss_data`.
"""
function dss2tppm_load!(tppm_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(tppm_data, "load")
        tppm_data["load"] = []
    end

    defaults = get_prop_default("load")

    if haskey(dss_data, "load")
        for load in dss_data["load"]
            merge!(defaults, load)

            loadDict = Dict{String,Any}()

            nphases = defaults["phases"]

            loadDict["load_bus"] = find_bus(defaults["bus1"], tppm_data)
            loadDict["pd"] = zeros(nphases)  # TODO:
            loadDict["qd"] = zeros(nphases)  # TODO:
            loadDict["status"] = convert(Int, defaults["enabled"])  # TODO: or check bus?

            loadDict["index"] = length(tppm_data["load"]) + 1

            used = ["phases", "bus1", "id"]
            PMs.import_remaining!(loadDict, defaults, import_all; exclude=used)

            push!(tppm_data["load"], loadDict)
        end
    end
end


"""
    dss2tppm_shunt!(tppm_data, dss_data)

Adds PowerModels-style shunts to `tppm_data` from `dss_data`.
"""
function dss2tppm_shunt!(tppm_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(tppm_data, "shunt")
        tppm_data["shunt"] = []
    end

    defaults = get_prop_default("capacitor")

    if haskey(dss_data, "capacitor")
        for shunt in dss_data["capacitor"]
            merge!(defaults, shunt)

            shuntDict = Dict{String,Any}()

            nphases = defaults["phases"]

            shuntDict["shunt_bus"] = find_bus(shunt["bus1"], tppm_data)
            shuntDict["name"] = defaults["id"]
            shuntDict["gs"] = zeros(nphases)  # TODO:
            shuntDict["bs"] = zeros(nphases)  # TODO:
            shuntDict["status"] = 1  # TODO:

            shuntDict["index"] = length(tppm_data["shunt"]) + 1

            used = ["bus1", "phases", "id"]
            PMs.import_remaining!(shuntDict, defaults, import_all; exclude=used)

            push!(tppm_data["shunt"], shuntDict)
        end
    end

    defaults = get_prop_default("reactor")

    if haskey(dss_data, "reactor")
        for shunt in dss_data["reactor"]
            merge!(defaults, shunt)

            shuntDict = Dict{String,Any}()

            nphases = defaults["phases"]

            shuntDict["shunt_bus"] = find_bus(shunt["bus1"], tppm_data)
            shuntDict["name"] = defaults["id"]
            shuntDict["gs"] = zeros(nphases)  # TODO:
            shuntDict["bs"] = zeros(nphases)  # TODO:
            shuntDict["status"] = 1  # TODO:

            shuntDict["index"] = length(tppm_data["shunt"]) + 1

            used = ["bus1", "phases", "id"]
            PMs.import_remaining!(shuntDict, defaults, import_all; exclude=used)

            push!(tppm_data["shunt"], shuntDict)
        end
    end
end


"""
    dss2tppm_gen!(tppm_data, dss_data)

Adds PowerModels-style generators to `tppm_data` from `dss_data`.
"""
function dss2tppm_gen!(tppm_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(tppm_data, "gen")
        tppm_data["gen"] = []
    end

    defaults = get_prop_default("generator")

    if haskey(dss_data, "generator")
        for gen in dss_data["generator"]
            merge!(defaults, gen)

            genDict = Dict{String,Any}()

            nphases = defaults["phases"]

            genDict["gen_bus"] = find_bus(defaults["bus1"])
            genDict["name"] = defaults["id"]
            genDict["gen_status"] = 1  # TODO:
            genDict["pg"] = zeros(nphases)  # TODO:
            genDict["qg"] = zeros(nphases)  # TODO:
            genDict["vg"] = ones(nphases)  # TODO:
            genDict["mbase"] = NaN  # TODO:
            genDict["ramp_agc"] = zeros(nphases)  # TODO:
            genDict["ramp_q"] = zeros(nphases)  # TODO:
            genDict["ramp_10"] = zeros(nphases)  # TODO:
            genDict["ramp_30"] = zeros(nphases)  # TODO:
            genDict["pmin"] = zeros(nphases)  # TODO:
            genDict["pmax"] = zeros(nphases)  # TODO:
            genDict["apf"] = zeros(nphases)  # TODO:
            genDict["qmin"] = zeros(nphases)  # TODO:
            genDict["qmax"] = zeros(nphases)  # TODO:
            genDict["pc1"] = zeros(nphases)  # TODO:
            genDict["pc2"] = zeros(nphases)  # TODO:
            genDict["qc1min"] = zeros(nphases)  # TODO:
            genDict["qc1max"] = zeros(nphases)  # TODO:
            genDict["qc2min"] = zeros(nphases)  # TODO:
            genDict["qc2max"] = zeros(nphases)  # TODO:

            genDict["model"] = 2 * ones(nphases)  # TODO:
            genDict["startup"] = zeros(nphases)  # TODO:
            genDict["shutdown"] = zeros(nphases)  # TODO:
            genDict["ncost"] = 3 * ones(nphases)  # TODO:
            genDict["cost"] = [0.0 1.0 0.0; 0.0 1.0 0.0; 0.0 1.0 0.0]

            genDict["index"] = length(tppm_data["gen"]) + 1

            used = ["id", "phases", "bus1"]
            PMs.import_remaining!(genDict, defaults, import_all; exlude=used)

            push!(tppm_data["gen"], genDict)
        end
    end
end


"""
    dss2tppm_branch!(tppm_data, dss_data)

Adds PowerModels-style branches to `tppm_data` from `dss_data`.
"""
function dss2tppm_branch!(tppm_data::Dict, dss_data::Dict, import_all::Bool)
    if !haskey(tppm_data, "branch")
        tppm_data["branch"] = []
    end

    defaults = get_prop_default("line")

    if haskey(dss_data, "line")
        for line in dss_data["line"]
            merge!(defaults, line)

            branchDict = Dict{String,Any}()

            nphases = warn_get(dss_data, "phases", 3, valid = x -> x > 0)

            branchDict["name"] = defaults["id"]

            branchDict["f_bus"] = find_bus(defaults["bus1"], tppm_data)
            branchDict["t_bus"] = find_bus(defaults["bus2"], tppm_data)

            line_length = warn_get(line, "length", 1.0, valid = x -> x > 0.0)
            branchDict["length"] = line_length
            linecode = warn_get(line, "linecode", "oh_horiz_3ph_4/0ph_1/0neu")
            branchDict["linecode"] = linecode

            # need to read in impedance matrices from linecode

            branchDict["br_r"] = parse_matrix(Float64, defaults["rmatrix"])*line_length
            branchDict["br_x"] = parse_matrix(Float64, defaults["xmatrix"])*line_length

            if haskey(dss_data, "rmatrix")
                branchDict["br_r"] = parse_matrix(Float64, dss_data["rmatrix"])*line_length
            else
                warn("rmatrix not present, setting to default")
            end

            if haskey(dss_data, "xmatrix")
                branchDict["br_"] = parse_matrix(Float64, dss_data["matrix"])*line_length
            else
                warn("matrix not present, setting to default")
            end


            # TODO: cmatrix, from linecode?
            branchDict["g_fr"] = zeros(nphases, nphases)  # TODO:
            branchDict["b_fr"] = zeros(nphases, nphases)  # TODO:
            branchDict["g_to"] = zeros(nphases, nphases)  # TODO:
            branchDict["b_to"] = zeros(nphases, nphases)  # TODO:

            if haskey(dss_data, "cmatrix")
                C_nf = parse_matrix(Float64, dss_data["cmatrix"])*line_length
                C = C/1e9
                G = 1.0/(2.0*pi*60.0*C)
                Gpu = Zbase*G
                branchDict["g_fr"] = Gpu/2.0
                branchDict["g_to"] = Gpu/2.0
            end


            # TODO: find out if this is per phase limit or across all phases
            # TODO: pick a better value for emergamps
            branchDict["rate_a"] = ones(nphases)*warn_get(dss_data, "normamps", 600.0; valid = x -> x > 0.0)
            branchDict["rate_b"] = ones(nphases)*warn_get(dss_data, "emergamps", 600.0; valid = x -> x > 0.0)
            branchDict["rate_c"] = ones(nphases)*warn_get(dss_data, "emergamps", 600.0; valid = x -> x > 0.0)

            branchDict["tap"] = ones(nphases)  
            branchDict["shift"] = zeros(nphases)  

            branchDict["br_status"] = warn_get(dss_data, "active", 1)

            branchDict["angmin"] = zeros(nphases)  
            branchDict["angmax"] = zeros(nphases) 

            branchDict["transformer"] = false

            branchDict["index"] = length(tppm_data["branch"]) + 1

            used = ["id", "phases", "bus1", "bus2", "rmatrix", "xmatrix"]
            PMs.import_remaining!(branchDict, defaults, import_all; exclude=used)

            push!(tppm_data["branch"], branchDict)
        end
    end
end


"""
    update_lookup_structure!(tppm_data)

Updates the Dict lookup structure from Dict{Array{Dict{String,Any}}} to
Dict{String,Dict{String,Any}}, requiring the presence of `"index"` in each
component.
"""
function update_lookup_structure!(tppm_data::Dict)
    for (k, v) in tppm_data
        if isa(v, Array)
            dict = Dict{String,Any}()
            for item in v
                assert("index" in keys(item))
                dict[string(item["index"])] = item
            end
            tppm_data[k] = dict
        end
    end
end


"Parses a Dict resulting from the parsing of a DSS file into a PowerModels usable format."
function parse_opendss(dss_data::Dict; import_all::Bool=false)::Dict
    tppm_data = Dict{String,Any}()

    parse_dss_with_dtypes!(dss_data, ["line", "load", "generator", "capacitor", "reactor"])

    tppm_data["per_unit"] = false
    tppm_data["source_type"] = "dss"
    tppm_data["source_version"] = VersionNumber("0")
    tppm_data["files"] = dss_data["filename"]


    if haskey(dss_data, "circuit")
        defaults = get_prop_default("circuit")
        merge!(defaults, dss_data["circuit"][1])

        tppm_data["name"] = defaults["id"]
        tppm_data["basekv"] = defaults["basekv"]
        tppm_data["baseMVA"] = defaults["basemva"]
        tppm_data["pu"] = defaults["pu"]
    else
        error(LOGGER, "Circuit not defined, not a valid circuit!")
    end

    dss2tppm_bus!(tppm_data, dss_data, import_all)
    dss2tppm_load!(tppm_data, dss_data, import_all)
    dss2tppm_shunt!(tppm_data, dss_data, import_all)
    dss2tppm_branch!(tppm_data, dss_data, import_all)
    dss2tppm_gen!(tppm_data, dss_data, import_all)

    update_lookup_structure!(tppm_data)

    return tppm_data
end


"Parses a DSS file into a PowerModels usable format."
function parse_opendss(filename::String; import_all::Bool=false)::Dict
    dss_data = parse_dss(filename)

    return parse_opendss(dss_data; import_all=import_all)
end
