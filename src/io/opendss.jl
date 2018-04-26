# OpenDSS parser


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
    isArray = false
    isProp = true
    endEquality = false
    str_out = ""

    for char in properties
        sstr_out = split(str_out, r"\s*=\s*")
        if length(sstr_out) == 2 && sstr_out[2] != ""
            endEquality = true
        else
            endEquality = false
        end

        if char == ' ' || char == properties[end]
            isProp = false
        else
            isProp = true
        end

        if char in ['[', '(']
            isArray = true
        elseif char in [']', ')']
            isArray = false
        end

        if char != ' ' || isArray
            str_out = string(str_out, char)
        end

        if str_out != "" && !isArray && !isProp && endEquality
            push!(propsOut, str_out)
            str_out = ""
        end
    end

    return propsOut
end


""
function add_component!(dss_data::Dict, ctype_name::AbstractString, compDict::Dict)
    ctype = split(lowercase(ctype_name), '.'; limit=2)[1]
    if haskey(dss_data, ctype)
        push!(dss_data[ctype], compDict)
    else
        dss_data[ctype] = [compDict]
    end
end


""
function parse_component(component::AbstractString, properties::AbstractString, compDict=Dict{String,Any}())
    ctype, name = split(lowercase(component), '.'; limit=2)

    if !haskey(compDict, "id")
        compDict["id"] = name
    end

    propArray = parse_properties(properties)

    for property in propArray
        key, value = split(property, '=')
        compDict[key] = value
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
        curCompDict = parse_component(curCtypeName, elements[3])
    end

    return curCtypeName, curCompDict
end


""
function parse_dss(filename::AbstractString)::Dict
    path = join(split(filename, '/')[1:end-1], '/')
    dss_str = readstring(open(filename))
    dss_data = Dict{String,Array}()

    curCompDict = Dict{String,Any}()
    curCtypeName = ""

    lines = split(dss_str, '\n')

    for (n, line) in enumerate(lines)
        debug(LOGGER, "$line")
        if startswith(line, '!') || line == ""
            # Is a comment, skip
            # TODO: are !!! lines comments?
            # TODO: strip all comment lines
            continue
        elseif startswith(strip(line, ' '), '~')
            curCompDict = parse_component(curCtypeName, strip(line, '~'), curCompDict)
            if n < length(lines) && startswith(lines[n+1], '~')
                continue
            else
                add_component!(dss_data, curCtypeName, curCompDict)
            end
        else
            curCompDict = Dict{String,Any}()
            line_elements = split(split(line, r"\s*\!")[1], r"\s+"; limit=3)
            cmd = line_elements[1]

            if lowercase(cmd) in ["clear", "solve", "calcvoltagebases", "show", "plot", "help"]
                continue

            elseif lowercase(cmd) == "redirect"
                file = line_elements[2]
                fullpath = path == "" ? file : join([path, file], '/')
                merge_dss!(dss_data, parse_dss(fullpath))

            elseif lowercase(cmd) == "compile"
                file = split(strip(line_elements[2], ['(',')']), '\\')[end]
                fullpath = path == "" ? file : join([path, file], '/')
                merge_dss!(dss_data, parse_dss(fullpath))

            elseif lowercase(cmd) == "set"
                if length(line_elements) == 2
                    property, value = split(line_elements[2], '='; limit=2)
                else
                    property, value = line_elements[2], strip(strip(line_elements[3], '='))
                end

                curCtypeName = "options.$(property)"
                curCompDict["$property"] = value

            elseif lowercase(cmd) == "buscoords"
                file = line_elements[2]
                fullpath = path == "" ? file : join([path, file], '/')
                dss_data["buscoords"] = parse_buscoords(fullpath)

            elseif lowercase(cmd) == "new"
                curCtypeName, curCompDict = parse_line(line_elements, curCtypeName)
            end
        end

        # warn(LOGGER, "$curCtypeName")

        if n < length(lines) && startswith(strip(lines[n + 1], ' '), '~')
            continue
        elseif length(curCompDict) > 0
            add_component!(dss_data, curCtypeName, curCompDict)
            curCompDict = Dict{String,Any}()
        else
            continue
        end

    end

    return dss_data
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
    tppm = Dict{String,Any}()

    tppm["per_unit"] = false
    tppm["source_type"] = "dss"
    tppm["source_version"] = VersionNumber("0")
    
    # TODO: add dss2tppm conversions

    update_lookup_structure!(tppm)

    return tppm
end


""
function parse_opendss(filename::String)::Dict
    dss_data = parse_dss(filename)

    return parse_opendss(dss_data)::Dict
end