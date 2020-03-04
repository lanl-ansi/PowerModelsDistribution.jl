const _linecode_properties = Vector{String}([
    "nphases", "r1", "x1", "r0","x0", "c1", "c0", "units", "rmatrix",
    "xmatrix", "cmatrix", "basefreq", "normamps", "emergamps", "faultrate",
    "pctperm", "repair", "kron", "rg", "xg", "rho", "neutral", "b1", "b0",
    "like"
])

const _linegeometry_properties = Vector{String}([
    "nconds", "nphases", "cond", "wire", "x", "h", "units", "normamps",
    "emergamps", "reduce", "spacing", "wires", "cncable", "tscable",
    "cncables", "tscables", "like"
])

const _linespacing_properties = Vector{String}([
    "nconds", "nphases", "x", "h", "units"
])

const _loadshape_properties = Vector{String}([
    "npts", "interval", "minterval", "sinterval", "pmult", "qmult", "hour",
    "mean", "stddev", "csvfile", "sngfile", "pqcsvfile", "action", "useactual",
    "pmax", "qmax", "pbase", "like"
])

const _growthshape_properties = Vector{String}([
    "npts", "year", "mult", "csvfile", "sngfile", "dblfile", "like"
])

const _tcc_curve_properties = Vector{String}([
    "npts", "c_array", "t_array", "like"
])

const _cndata_properties = Vector{String}([])

const _tsdata_properties = Vector{String}([])

const _wiredata_properties = Vector{String}([
    "rdc", "rac", "runits", "gmrac", "gmrunits", "radius", "radunits",
    "normamps", "emergamps", "diam", "like"
])

const _xfmrcode_properties = Vector{String}([])

const _vsource_properties = Vector{String}([
    "bus1", "bus2", "basekv", "pu", "angle", "frequency", "phases", "mvasc3",
    "mvasc1", "x1r1", "x0r0", "isc3", "isc1", "r1", "x1", "r0", "x0",
    "scantype", "sequence", "spectrum", "z1", "z2", "z0", "puz1", "puz2",
    "puz0", "basemva", "basefreq", "like", "enabled"
])

const _isource_properties = Vector{String}([
    "phases", "bus1", "amps", "angle", "frequency", "scantype", "sequence",
    "spectrum", "basefreq", "enabled", "like"
])

const _fault_properties = Vector{String}([
    "phases", "bus1", "bus2", "r", "gmatrix", "minamps", "ontime", "pctperm",
    "temporary", "%stddev", "normamps", "emergamps", "basefreq", "faultrate",
    "repair", "enabled", "like"
])

const _capacitor_properties = Vector{String}([
    "bus1", "bus2", "phases", "kvar", "kv", "conn", "cmatrix", "cuf", "r",
    "xl", "harm", "numsteps", "states", "normamps", "emergamps", "faultrate",
    "pctperm", "basefreq", "enabled", "like"
])

const _line_properties = Vector{String}([
    "bus1", "bus2", "linecode", "length", "phases", "r1", "x1", "r0", "x0",
    "c1", "c0", "b1", "b0", "normamps", "emergamps", "faultrate", "pctperm",
    "repair", "basefreq", "rmatrix", "xmatrix", "cmatrix", "switch", "rg",
    "xg", "rho", "geometry", "earthmodel", "units", "enabled", "like"
])

const _reactor_properties = Vector{String}([
    "phases", "bus1", "bus2", "kv", "kvar", "conn", "parallel", "r", "rmatrix",
    "rp", "x", "xmatrix", "z", "z1", "z2", "z0", "rcurve", "lcurve", "lmh",
    "normamps", "emergamps", "repair", "faultrate", "pctperm", "basefreq",
    "enabled", "like"
])

const _transformer_properties = Vector{String}([
    "phases", "windings", "wdg", "bus", "conn", "kv", "kva", "tap", "%r",
    "rneut", "xneut", "buses", "conns", "kvs", "kvas", "taps", "%rs", "xhl",
    "xlt", "xht", "xscarray", "thermal", "n", "m", "flrise", "hsrise",
    "%loadloss", "%noloadloss", "%imag", "ppm_antifloat", "normhkva",
    "emerghkva", "sub", "maxtap", "mintap", "numtaps", "subname", "bank",
    "xfmrcode", "xrconst", "leadlag", "faultrate", "basefreq", "enabled",
    "like"
])

const _gictransformer_properties = Vector{String}([
    "basefreq", "bush", "busnh", "busnx", "busx", "emergamps", "enabled",
    "phases", "r1", "r2", "type", "mva", "kvll1", "kvll2", "%r1", "%r2", "k",
    "varcurve", "like", "normamps", "emergamps", "pctperm", "repair"
])

const _gicline_properties = Vector{String}([
    "angle", "bus1", "bus2", "c", "ee", "en", "frequency", "lat1", "lat2",
    "lon1", "lon2", "phases", "r", "volts", "x", "like", "basefreq",
    "enabled", "spectrum"
])

const _load_properties = Vector{String}([
    "phases", "bus1", "kv", "kw", "pf", "model", "yearly", "daily", "duty",
    "growth", "conn", "kvar", "rneut", "xneut", "status", "class", "vminpu",
    "vmaxpu", "vminnorm", "vminemerg", "xfkva", "allocationfactor", "kva",
    "%mean", "%stddev", "cvrwatts", "cvrvars", "kwh", "kwhdays", "cfactor",
    "cvrcurve", "numcust", "spectrum", "zipv", "%seriesrl", "relweight",
    "vlowpu", "puxharm", "xrharm", "spectrum", "basefreq", "enabled", "like"
])

const _generator_properties = Vector{String}([
    "bus1", "phases", "kv", "kw", "pf", "model", "yearly", "daily", "duty",
    "dispvalue", "conn", "kvar", "rneut", "xneut", "status", "class", "vpu",
    "maxkvar", "minkvar", "pvfactor", "debugtrace", "vminpu", "vmaxpu",
    "forceon", "kva", "mva", "xd", "xdp", "xdpp", "h", "d", "usermodel",
    "userdata", "shaftmodel", "shaftdata", "dutystart", "balanced", "xrdp",
    "spectrum", "basefreq", "enabled", "like"
])

const _indmach012_properties = Vector{String}([
    "phases", "bus1", "kv", "kw", "pf", "conn", "kva", "h", "d", "purs",
    "puxs", "purr", "puxr", "puxm", "slip", "maxslip", "slipoption",
    "spectrum", "enabled"
])

const _storage_properties = Vector{String}([
    "phases", "bus1", "%charge", "%discharge", "%effcharge", "%idlingkvar",
    "idlingkw", "%r", "%reserve", "%stored", "%x", "basefreq", "chargetrigger",
    "class", "conn", "daily", "yearly", "debugtrace", "dischargetrigger",
    "dispmode", "duty", "dynadata", "dynadll", "enabled", "kv", "kva", "kvar",
    "kw", "kwhrated", "kwhstored", "kwrated", "like", "model", "pf",
    "spectrum", "state", "timechargetrig", "userdata", "usermodel", "vmaxpu",
    "vminpu", "yearly"
])

const _capcontrol_properties = Vector{String}([
    "element", "capacitor", "type", "ctphase", "ctratio", "deadtime", "delay",
    "delayoff", "eventlog", "offsetting", "onsetting", "ptphase", "ptratio",
    "terminal", "vbus", "vmax", "vmin", "voltoverride", "enabled"
])

const _regcontrol_properties = Vector{String}([
    "transformer", "winding", "vreg", "band", "delay", "ptratio", "ctprim",
    "r", "x", "pthase", "tapwinding", "bus", "remoteptratio", "debugtrace",
    "eventlog", "inversetime", "maxtapchange", "revband", "revdelay",
    "reversible", "revneutral", "revr", "revthreshold", "revvreg", "revx",
    "tapdelay", "tapnum", "vlimit", "ldc_z", "rev_z", "cogen", "enabled"
])

const _energymeter_properties = Vector{String}([
    "element", "terminal", "action", "clear", "save", "take", "option",
    "kwnorm", "kwemerg", "peakcurrent", "zonelist", "zonelist", "zonelist",
    "localonly", "mask", "losses", "linelosses", "xfmrlosses", "seqlosses",
    "3phaselosses", "vbaselosses", "basefreq", "enabled", "like"
])

const _pvsystem_properties = Vector{String}([
    "phases", "bus1", "kv", "irradiance", "pmpp", "temperature", "pf",
    "conn", "kvar", "kva", "%cutin", "%cutout", "effcurve", "p-tcurve", "%r",
    "%x", "model", "vminpu", "vmaxpu", "yearly", "daily", "duty", "tyearly",
    "tduty", "class", "usermodel", "userdata", "debugtrace", "spectrum"
])

const _relay_properties = Vector{String}([])

const _recloser_properties = Vector{String}([])

const _fuse_properties = Vector{String}([])

const _dss_object_properties = Dict{String,Vector{String}}(
    "linecode" => _linecode_properties,
    "linegeometry" => _linegeometry_properties,
    "linespacing" => _linespacing_properties,
    "loadshape" => _loadshape_properties,
    "growthshape" => _growthshape_properties,
    "tcc_curve" => _tcc_curve_properties,
    "cndata" => _cndata_properties,
    "tsdata" => _tsdata_properties,
    "wiredata" => _wiredata_properties,
    "xfmrcode" => _xfmrcode_properties,
    "vsource" => _vsource_properties,
    "circuit" => _vsource_properties,  # alias circuit to vsource
    "isource" => _isource_properties,
    "fault" => _fault_properties,
    "capacitor" => _capacitor_properties,
    "line" => _line_properties,
    "reactor" => _reactor_properties,
    "transformer" => _transformer_properties,
    "gictransformer" => _gictransformer_properties,
    "gicline" => _gicline_properties,
    "load" => _load_properties,
    "generator" => _generator_properties,
    "indmach012" => _indmach012_properties,
    "storage" => _storage_properties,
    "capcontrol" => _capacitor_properties,
    "regcontrol" => _regcontrol_properties,
    "energymeter" => _energymeter_properties,
    "pvsystem" => _pvsystem_properties,
    "relay" => _relay_properties,
    "recloser" => _reactor_properties,
    "fuse" => _fuse_properties
)


"parses single column load profile files"
function _parse_loadshape_csv_file(path::AbstractString, type::AbstractString; header::Bool=false, column::Int=1, interval::Bool=false)
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
function _parse_binary_file(path::AbstractString, precision::Type; npts::Union{Int,Nothing}=nothing, interval::Bool=false)
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
        return _parse_loadshape_csv_file(path, type; header=header, column=column, interval=interval)
    elseif type in ["sngfile", "dblfile"]
        return _parse_binary_file(path, Dict("sngfile" => Float32, "dblfile" => Float64)[type]; npts=npts, interval=interval)
    end
end


"parses pmult and qmult entries on loadshapes"
function _parse_mult_parameter(mult_string::AbstractString; path::AbstractString="", npts::Union{Int,Nothing}=nothing)
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
             current_obj[prop] = _parse_mult_parameter(current_obj[prop]; path=path, npts=npts)
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
Parses a Bus Coordinate `file`, in either "dat" or "csv" formats, where in
"dat", columns are separated by spaces, and in "csv" by commas. File expected
to contain "bus,x,y" on each line.
"""
function _parse_buscoords_file(file::AbstractString)::Dict{String,Any}
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
Parses a `component` with `properties` into a `object`. If `object` is not
defined, an empty dictionary will be used. Assumes that unnamed properties are
given in order, but named properties can be given anywhere.
"""
function _parse_component(component::AbstractString, properties::AbstractString, object::Dict{String,<:Any}=Dict{String,Any}(); path::AbstractString="")
    obj_type, name = split(component, '.'; limit=2)

    if !haskey(object, "name")
        object["name"] = name
    end

    property_array = _parse_properties(properties)

    property_names = _dss_object_properties[obj_type]
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
            value = _parse_mult_parameter(value; path=path)
        end

        _add_property(object, key, value)
    end

    return object
end


"""
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


""
function parse_dss(filename::AbstractString)::Dict
    data_dss = open(filename) do io
        parse_dss(io)
    end
    return data_dss
end


"""
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
                data_dss["buscoords"] = _parse_buscoords_file(full_path)

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

    parse_dss_with_dtypes!(data_dss)

    parse_dss_options!(data_dss)

    return data_dss
end
