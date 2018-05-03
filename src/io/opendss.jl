# OpenDSS parser


components = ["linecode", "linegeometry", "line", "linespacing", "loadshape",
              "growthshape", "tcc_curve", "wiredata", "xfmrcode", "vsource",
              "isource", "fault", "capacitor", "reactor", "transformer",
              "gictransformer", "gicline", "load", "generator", "indmach012",
              "storage", "capcontrol", "regcontrol", "energymeter", "monitor"]


""
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


""
function get_prop_default(ctype::AbstractString, i::Int)
    return get_prop_default(ctype)[i]
end


""
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


""
function get_prop_name(ctype::AbstractString, i::Int)::String
    return get_prop_name(ctype)[i]
end


""
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


""
function parse_array(dtype::Type, data::AbstractString)::Array
    elements = split(strip(data, ['\"', '\'', '[', ']', '(', ')']))
    array = zeros(dtype, length(elements))

    for (i, el) in enumerate(elements)
        array[i] = parse(dtype, el)
    end

    return array
end


""
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


""
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


""
function add_component!(dss_data::Dict, ctype_name::AbstractString, compDict::Dict)
    debug(LOGGER, "add_component! $ctype_name")
    ctype = split(lowercase(ctype_name), '.'; limit=2)[1]
    if haskey(dss_data, ctype)
        push!(dss_data[ctype], compDict)
    else
        dss_data[ctype] = [compDict]
    end
end


""
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
        compDict[lowercase(key)] = value
    end

    return compDict
end


""
function merge_dss!(dss_prime::Dict{String,Array}, dss_to_add::Dict{String,Array})
    for (k, v) in dss_to_add
        if k in keys(dss_prime)
            append!(dss_prime[k], v)
        else
            dss_prime[k] = v
        end
    end
end


""
function parse_line(elements::Array, curCtypeName::AbstractString, curCompDict::Dict=Dict{String,Any}())
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


""
function strip_comments(line::AbstractString)::String
    return split(line, r"\s*!")[1]
end


""
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


""
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
                curCtypeName, curCompDict = parse_line(line_elements, curCtypeName)
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


""
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


""
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


""
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


""
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


""
function dss2tppm_shunt!(tppm_data::Dict, dss_data::Dict)
    
end


""
function dss2tppm_branch!(tppm_data::Dict, dss_data::Dict)
    
end


""
function dss2tppm_gen!(tppm_data::Dict, dss_data::Dict)
    
end


""
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


""
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


""
function parse_opendss(filename::String)::Dict
    dss_data = parse_dss(filename)

    return parse_opendss(dss_data)::Dict
end