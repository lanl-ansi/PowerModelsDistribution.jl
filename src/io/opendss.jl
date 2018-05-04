# OpenDSS parser


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
function get_prop_default(ctype::AbstractString)::Array

    line = []
    load = []
    transformer = []

    ctypes = Dict{String, Array}("line" => line,
                                 "load" => load,
                                 "transformer" => transformer)
    try
        return ctypes[ctype]
    catch KeyError
        return []
    end
end


"""
    get_prop_default(ctype, i)

Returns the default property value, or the expected Type if no default is
known, of the `i`th property, for a given component type `ctype`.
"""
function get_prop_default(ctype::AbstractString, i::Int)
    return get_prop_default(ctype)[i]
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

    transformer = ["phases", "windings"]

    ctypes = Dict{String, Array}("line" => line,
                                 "load" => load,
                                 "transformer" => transformer)

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

    propNames = get_prop_name(ctype)

    for (n, property) in enumerate(propArray)
        if !contains(property, "=")
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


"Strips comments, defined by "!" from the ends of lines"
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
    # TODO: parse transformers special case
    path = join(split(filename, '/')[1:end-1], '/')
    dss_str = readstring(open(filename))
    dss_data = Dict{String,Array}()

    dss_data["filename"] = [split(filename, "/")[end]]

    curCompDict = Dict{String,Any}()
    curCtypeName = ""

    lines = split(dss_str, '\n')

    stripped_lines = []
    for line in lines
        if !startswith(line, '!') && line != ""
            push!(stripped_lines, line)
        end
    end

    nlines = length(stripped_lines)

    for (n, line) in enumerate(stripped_lines)
        debug(LOGGER, "LINE $(find(lines .== line)): $line")
        line = strip_comments(line)

        if startswith(strip(line), '~')
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

            if cmd == "redirect"
                file = line_elements[2]
                fullpath = path == "" ? file : join([path, file], '/')
                merge_dss!(dss_data, parse_dss(fullpath))

            elseif cmd == "compile"
                file = split(strip(line_elements[2], ['(',')']), '\\')[end]
                fullpath = path == "" ? file : join([path, file], '/')
                merge_dss!(dss_data, parse_dss(fullpath))

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
                    warn(LOGGER, "Command \"$cmd\" not supported, skipping.")
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
function dss2tppm_bus!(tppm_data::Dict, dss_data::Dict)
    if !haskey(tppm_data, "bus")
        tppm_data["bus"] = []
    end

    buses = discover_buses(dss_data)

    for (n, (bus, nodes)) in enumerate(buses)
        busDict = Dict{String,Any}()

        busDict["bus_i"] = n  # TODO: how to properly index?
        busDict["index"] = n  # TODO: how to properly index?
        busDict["name"] = bus

        busDict["vm"]  # TODO: zeros(numPhases)?
        busDict["va"]  # TODO: zeros(numPhases)?

        busDict["vmin"]  # TODO: implicit limits?
        busDict["vmax"]  # TODO: implicit limits?

        push!(tppm_data["bus"], busDict)
    end
end


"""
    dss2tppm_load!(tppm_data, dss_data)

Adds PowerModels-style loads to `tppm_data` from `dss_data`.
"""
function dss2tppm_load!(tppm_data::Dict, dss_data::Dict)
    if !haskey(tppm_data, "load")
        tppm_data["load"] = []
    end

    if haskey(dss_data, "load")
        for load in dss_data["load"]
            loadDict = Dict{String,Any}()


            push!()
        end
    end
end


"""
    dss2tppm_shunt!(tppm_data, dss_data)

Adds PowerModels-style shunts to `tppm_data` from `dss_data`.
"""
function dss2tppm_shunt!(tppm_data::Dict, dss_data::Dict)

end


"""
    dss2tppm_branch!(tppm_data, dss_data)

Adds PowerModels-style branches to `tppm_data` from `dss_data`.
"""
function dss2tppm_branch!(tppm_data::Dict, dss_data::Dict)

end


"""
    dss2tppm_gen!(tppm_data, dss_data)

Adds PowerModels-style generators to `tppm_data` from `dss_data`.
"""
function dss2tppm_gen!(tppm_data::Dict, dss_data::Dict)

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
            #println("updating $(k)")
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
function parse_opendss(dss_data::Dict)::Dict
    tppm_data = Dict{String,Any}()

    tppm_data["per_unit"] = false
    tppm_data["source_type"] = "dss"
    tppm_data["source_version"] = VersionNumber("0")
    tppm_data["name"] = haskey(dss_data, "circuit") ? dss_data["circuit"]["id"] : dss_data["filename"]

    # TODO: add dss2tppm conversions
    dss2tppm_bus!(tppm_data, dss_data)

    update_lookup_structure!(tppm_data)

    return tppm_data
end


"Parses a DSS file into a PowerModels usable format."
function parse_opendss(filename::String)::Dict
    dss_data = parse_dss(filename)

    return parse_opendss(dss_data)::Dict
end