"strips lines that are either commented (block or single) or empty"
function _strip_lines!(lines::Vector{String})::Vector{String}
    blockComment = false
    for (i,line) in enumerate(lines)
        if startswith(line, "/*") || endswith(line, "*/")
            blockComment = !blockComment
            lines[i] = ""
        elseif !startswith(line, '!') && !startswith(line, "//") && strip(line) != "" && !blockComment
            continue
        else
            lines[i] = ""
        end
    end
    return lines
end


"""
    _sanatize_line(line::String)::String

Sanitizes lines by stripping them clean of extra space and the beginnging and end, making
everything lowercase, changing `~` or `m` to `more`, and stripping comments
"""
function _sanatize_line(line::String)::String
    # strip comments, spaces, make lowercase
    string(lowercase(replace(replace(strip(split(strip(line), r"(\s*(?:\/\/|!))")[1], ['\r', '\n']), r"^[~m]\s+"=>"more "), r"\s*=\s*"=>"=")))
end


"""
    _parse_command_from_line(line::String)::Tuple{String,String}

Parses the dss command from the line (if present), optionnally making the command `set` if
implicitly used.
"""
function _parse_command_from_line(line::String)::Tuple{String,String}
    line_parts = string.(split(line, r"\s+"; limit=2))
    cmd, elements = "", ""

    if length(line_parts) >= 1 && any(startswith(line_parts[1], String(k)) for k in fieldnames(OpenDssRawDataModel))
        cmd = "edit"
        _elements = split(line_parts[1], '.')
        elements = "$(join(_elements[1:2], '.')) $(join(_elements[3:end], '.'))"
        if length(line_parts) == 2
            elements = "$elements $(line_parts[2])"
        end
    elseif length(line_parts) == 2
        cmd, elements = line_parts
    else
        cmd, elements = line_parts[1], ""
    end

    return (cmd, elements)
end


""
function _parse_dss_obj_type_name(dss_obj_type_name::AbstractString)::Vector{String}
    obj_type_name = string.(split(replace(replace(dss_obj_type_name, "\""=>""), "\'"=>""), '.'))
end


""
function _parse_dss_cmd_new!(data_dss::DssRawModel, line::String, line_number::Int)::Nothing
    orig_dss_obj_name = missing
    orig_dss_obj_type = missing

    rmatches_iterator = eachmatch(_dss_cmd_new_regex, line)

    dss_obj_match = iterate(rmatches_iterator)[1]

    dss_obj_type, dss_obj_name = _parse_dss_obj_type_name(dss_obj_match.captures[3])

    if dss_obj_type == "circuit"
        orig_dss_obj_name = dss_obj_name
        orig_dss_obj_type = dss_obj_type

        dss_obj_name = "source"
        dss_obj_type = "vsource"
    end

    if !(dss_obj_type == data_dss.current_state.active_obj_type && dss_obj_name == data_dss.current_state.active_obj_name)
        data_dss.current_state.active_obj_field = ""
    end
    data_dss.current_state.active_obj_type = dss_obj_type
    data_dss.current_state.active_obj_name = dss_obj_name

    if !haskey(data_dss[dss_obj_type], dss_obj_name)
        data_dss[dss_obj_type][dss_obj_name] = Pair{String,String}["name"=>dss_obj_name]
    end
    if !ismissing(orig_dss_obj_type) && !haskey(data_dss[orig_dss_obj_type], orig_dss_obj_name)
        data_dss[orig_dss_obj_type][orig_dss_obj_name] = Pair{String,String}["name"=>orig_dss_obj_name]
    end

    push!(data_dss[dss_obj_type][dss_obj_name], "__path__" => dirname(data_dss.current_state.current_file))
    !ismissing(orig_dss_obj_type) && push!(data_dss[dss_obj_type][dss_obj_name], "__path__" => dirname(data_dss.current_state.current_file))

    for (i, _match) in enumerate(rmatches_iterator)
        if i > 1
            property_key, property_value = _parse_match_element(_match, dss_obj_type; prev_k=data_dss.current_state.active_obj_field)
            data_dss.current_state.active_obj_field = property_key
            push!(data_dss[dss_obj_type][dss_obj_name], property_key => property_value)
            !ismissing(orig_dss_obj_type) && push!(data_dss[dss_obj_type][dss_obj_name], property_key => property_value)
        end
    end
end


""
function _parse_dss_cmd_more!(data_dss::DssRawModel, line::String, line_number::Int)::Nothing
    for _match in eachmatch(_dss_cmd_more_regex, line)
        property_key, property_value = _parse_match_element(_match, data_dss.current_state.active_obj_type; prev_k=data_dss.current_state.active_obj_field)
        data_dss.current_state.active_obj_field = property_key

        push!(getproperty(data_dss, Symbol(data_dss.current_state.active_obj_type))[data_dss.current_state.active_obj_name], property_key => strip(property_value))
    end
end


""
function _parse_dss_cmd_set!(data_dss::DssRawModel, line::String, line_number::Int)::Nothing
    for _match in eachmatch(_dss_cmd_set_regex, line)
        property_key, property_value = _parse_match_element(_match, "options")
        push!(data_dss["options"], property_key => property_value)
    end
end


""
function _parse_dss_cmd_edit!(data_dss::DssRawModel, line::String, line_number::Int)::Nothing
    rmatches_iterator = eachmatch(_dss_cmd_new_regex, line)

    dss_obj_match = iterate(rmatches_iterator)[1]

    dss_obj_type, dss_obj_name = _parse_dss_obj_type_name(dss_obj_match.captures[3])
    if haskey(data_dss[dss_obj_type], dss_obj_name)
        _parse_dss_cmd_new!(data_dss, line, line_number)
    else
        @warn "Cannot find $(dss_obj_type) object $(dss_obj_name)."
    end
end


""
function _parse_dss_cmd_redirect!(data_dss::DssRawModel, line::String, line_number::Int)::Union{Nothing,DssRawModel}
    file = replace(line, r"(redirect|compile)\s*"=>"")
    file_path = FilePaths.Path(strip(file, ['(',')']))
    current_file = data_dss.filename[end]

    if !(joinpath(file_path...) in data_dss.filename)
        full_path = joinpath(data_dss.current_state.base_path, file_path)
        @info "Redirecting to '$(joinpath(file_path...))' on line $line_number in '$(basename(current_file))'"
        data_dss = parse_raw_dss(full_path, data_dss)
    end
end


""
_parse_dss_cmd_compile!(data_dss::DssRawModel, line::String, line_number::Int) = _parse_dss_cmd_redirect!(data_dss, line, line_number)


""
function _parse_match_element(rmatch::RegexMatch, dss_obj_type::String; prev_k::String="")::Pair{String,String}
    k = rmatch.captures[2]
    if isnothing(k)
        prop_keys = String[x for x in filter(x->x!="name", string.(fieldnames(getfield(PowerModelsDistribution, Symbol("Dss$(titlecase(dss_obj_type))")))))]
        if isempty(prev_k)
            k = prop_keys[1]
        else
            k = circshift(prop_keys, -(findfirst(x->x==prev_k, prop_keys)))[1]
        end
    end

    _sanitize_property_name(string(k), dss_obj_type) => strip(isnothing(rmatch.captures[5]) ? rmatch.captures[6] : rmatch.captures[5])
end


""
function _sanitize_property_name(property_name::String, dss_obj_type::String)::String
    get(get(_dss_property_renames, dss_obj_type, Dict{String,String}()), property_name, property_name)
end


""
function _parse_dss_cmd_clear!(data_dss::DssRawModel, line::String, line_number::Int)::DssRawModel
    @info "Circuit has been reset with the 'clear' on line $(line_number) in '$(basename(data_dss.current_state.current_file))'"
    data_dss = _init_dss_data(; current_file=data_dss.current_state.current_file)
end


""
function _init_dss_data(; current_file::Union{Missing,FilePaths.AbstractPath}=missing, current_command::String="")::DssRawModel
    data_dss = OpenDssRawDataModel()

    if !ismissing(current_file)
        push!(data_dss.filename, current_file)
        data_dss.current_state.base_path = dirname(current_file)
    end

    if !isempty(current_command)
        data_dss.current_state.current_command = current_command
    end

    return data_dss
end


""
function _parse_dss_cmd_buscoords!(data_dss::DssRawModel, file::String, line_number::Int)::Nothing
    file_path = split(strip(file, ['(',')']), r"[\\|\/]")
    full_path = joinpath(data_dss.current_state.base_path, file_path...)

    @info "Reading Buscoords in '$(basename(full_path))' on line $(line_number) in '$(basename(data_dss.current_state.current_file))'"
    open(first(Glob.glob([Glob.FilenameMatch(basename(full_path), "i")], string(dirname(full_path)))), "r") do io
        buscoords = _parse_dss_cmd_buscoords(io)
        for bc in buscoords
            push!(data_dss.buscoordinates, bc)
        end
    end
end


""
function _parse_dss_cmd_buscoords(io::IO)::Vector{DssBuscoords}
    buscoords = DssBuscoords[]
    for line in filter(x->!isempty(x), _strip_lines!(_sanatize_line.(readlines(io))))
        bus, x, y = split(line, _dss_cmd_buscoords_regex; limit=3)
        push!(buscoords, DssBuscoords(bus, parse(Float64, strip(x, [','])), parse(Float64, strip(y, [',', '\r']))))
    end

    return buscoords
end


""
_parse_dss_cmd_latloncoords!(data_dss::DssRawModel, file::String, line_number::Int)::Vector{DssBuscoords} = _parse_dss_cmd_buscoords!(data_dss, file, line_number)


""
function _parse_dss_cmd_setbusxy!(data_dss::DssRawModel, line::String, line_number::Int)::Vector{DssBuscoords}
    properties = Dict{String,String}()
    for _match in eachmatch(_dss_cmd_new_regex, line)
        property_key, property_value = _parse_match_element(_match, "setbusxy"; prev_k=data_dss.current_state.active_obj_field)
        data_dss.current_state.active_obj_field = property_key
        properties[property_key] = property_value
    end

    push!(data_dss.buscoordinates, DssBuscoords([k ∈ ["x", "y"] ? parse(Float64, properties[k]) : properties[k] for k in ["bus", "x", "y"]]...))
end

""
function _parse_dss_cmd_disable!(data_dss::DssRawModel, line::String, line_number::Int)::Nothing
    dss_obj_type, dss_obj_name = _parse_dss_obj_type_name(line)

    push!(data_dss[dss_obj_type][dss_obj_name], "enabled"=>false)
end


""
function _parse_dss_cmd_enable!(data_dss::DssRawModel, line::String, line_number::Int)::Nothing
    dss_obj_type, dss_obj_name = _parse_dss_obj_type_name(line)

    push!(data_dss[dss_obj_type][dss_obj_name], "enabled"=>true)
end


""
parse_raw_dss(filename::AbstractString, data::Missing=missing)::DssRawModel = parse_raw_dss(FilePaths.Path(filename), _init_dss_data())
parse_raw_dss(filename::FilePaths.AbstractPath, data::Missing=missing)::DssRawModel = parse_raw_dss(filename, _init_dss_data())


"""
    parse_raw_dss(filename::String)::DssRawModel

Parses a OpenDSS file given by `filename` into a Dict{Array{Dict}}. Only
supports components and options, but not commands, e.g. "plot" or "solve".
Will also parse files defined inside of the originating DSS file via the
"compile", "redirect" or "buscoords" commands.
"""
function parse_raw_dss(path::FilePaths.AbstractPath, data::DssRawModel)::DssRawModel
    open(first(Glob.glob([Glob.FilenameMatch(basename(path), "i")], string(dirname(path))))) do io
        parse_raw_dss(io, data)
    end
end


""
parse_raw_dss(io::IO, data::Missing=missing)::DssRawModel = parse_raw_dss(io, OpenDssRawDataModel())


"""
    parse_raw_dss(io::IO)::DssRawModel

Parses a OpenDSS file aleady in IO into a Dict{Array{Dict}}. Only
supports components and options, but not commands, e.g. "plot" or "solve".
Will also parse files defined inside of the originating DSS file via the
"compile", "redirect" or "buscoords" commands.
"""
function parse_raw_dss(io::IO, data::DssRawModel)::DssRawModel
    file = isa(io, IOStream) ? match(r"^<file\s(.+)>$", io.name).captures[1] : "GenericIOBuffer"

    push!(data.filename, file)
    data.current_state.current_file = FilePaths.Path(file)
    if ismissing(data.current_state.base_path)
        data.current_state.base_path = FilePaths.Path(dirname(file))
    end

    for (line_number, (cmd, elements)) in enumerate(_parse_command_from_line.(_strip_lines!(_sanatize_line.(readlines(io)))))
        # @warn line_number cmd elements
        if !isempty(cmd)
            if cmd ∉ _dss_supported_commands
                @info "Command '$cmd' on line $line_number in '$(basename(file))' is not supported, skipping."
                continue
            else
                data.current_state.last_command = data.current_state.current_command
                data.current_state.current_command = cmd

                getfield(PowerModelsDistribution, Symbol("_parse_dss_cmd_$(cmd)!"))(data, elements, line_number)
            end
        else
            continue
        end
    end

    return data
end


"""
    parse_dss_voltages_export(file::String)::Dict{String,Any}

Parses a voltages CSV file exported from OpenDSS with the command `Export Voltages [filename.csv]`
"""
function parse_dss_voltages_export(path::Union{String,FilePaths.AbstractPath})::Dict{String,Any}
    open(first(Glob.glob([Glob.FilenameMatch(basename(path), "i")], dirname(path))), "r") do io
        parse_dss_voltages_export(io)
    end
end


"""
    parse_dss_voltages_export(io::IO)::Dict{String,Any}

Parses a voltages CSV file exported from OpenDSS with the command `Export Voltages [filename.csv]`

Units of `vm` are volts, and units of `va` are degrees.
"""
function parse_dss_voltages_export(io::IO)::Dict{String,Any}
    data = Dict{String,Any}()
    for row in CSV.File(IOBuffer(lowercase(read(io, String))); normalizenames=true)
        data[row.bus] = Dict{String,Any}(
            "vbase" => row.basekv * 1e3,  # convert to volts, same units as vm
            "terminals" => Int[],
            "vm" => Real[],
            "va" => Real[],
        )
        nodes = [parse(Int, match(r"(\d+)", string(n)).captures[1]) for n in propertynames(row) if startswith(string(n), "node")]
        for i in nodes
            terminal = getproperty(row, Symbol("node$i"))
            if terminal != 0
                push!(data[row.bus]["terminals"], terminal)
                push!(data[row.bus]["vm"], getproperty(row, Symbol("magnitude$i")))
                push!(data[row.bus]["va"], getproperty(row, Symbol("angle$i")))
            end
        end
    end

    return data
end
