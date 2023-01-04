const _dss_unsupported_commands = String[
    "cleanup", "connect", "disconnect", "_docontrolactions", "_initsnap",
    "_samplecontrols", "_showcontrolqueue", "_solvedirect", "solvenocontrol",
    "_solvepflow", "?", "about", "addbusmarker", "alignfile", "allocateloads",
    "batchedit", "buildy", "calcvoltagebases", "capacity", "cd", "cktlosses",
    "clearbusmarkers", "close", "closedi", "comparecases", "currents",
    "di_plot", "distribute", "doscmd", "dump", "estimate", "export",
    "fileedit", "finishtimestep", "formedit", "get", "guids", "help", "init",
    "interpolate", "losses", "makebuslist", "makeposseq", "nodelist",
    "nodediff", "obfuscate", "open", "phaselosses", "plot", "powers",
    "pstcalc", "reconductor", "reduce", "relcalc", "remove", "rephase",
    "reprocessbuses", "reset", "rotate", "sample", "save", "select",
    "seqcurrents", "seqpowers", "seqvoltages", "setkvbase", "show", "solve",
    "summary", "totals", "updatestorage", "var", "variable", "varnames",
    "vdiff", "visualize", "voltages", "yearlycurves", "ysc", "zsc", "zsc10",
    "zscrefresh"
]

const _linecode_properties = String[
    "nphases", "r1", "x1", "r0","x0", "c1", "c0", "units", "rmatrix",
    "xmatrix", "cmatrix", "basefreq", "normamps", "emergamps", "faultrate",
    "pctperm", "repair", "kron", "rg", "xg", "rho", "neutral", "b1", "b0",
    "like"
]

const _linegeometry_properties = String[
    "nconds", "nphases", "cond", "wire", "x", "h", "units", "normamps",
    "emergamps", "reduce", "spacing", "wires", "cncable", "tscable",
    "cncables", "tscables", "like"
]

const _linespacing_properties = String[
    "nconds", "nphases", "x", "h", "units", "like"
]

const _loadshape_properties = String[
    "npts", "interval", "minterval", "sinterval", "pmult", "qmult", "hour",
    "mean", "stddev", "csvfile", "sngfile", "pqcsvfile", "action", "useactual",
    "pmax", "qmax", "pbase", "like"
]

const _xycurve_properties = String[
    "npts", "points", "yarray", "xarray", "csvfile", "sngfile", "dblfile",
    "x", "y", "xshift", "yshift", "xscale", "yscale", "like"
]

const _growthshape_properties = String[
    "npts", "year", "mult", "csvfile", "sngfile", "dblfile", "like"
]

const _tcc_curve_properties = String[
    "npts", "c_array", "t_array", "like"
]

const _cndata_properties = String[
    "diacable", "diains", "diam", "diastrand", "emergamps", "epsr", "gmrac",
    "gmrstrang", "gmrunits", "inslayer", "k", "like", "normamps", "rac",
    "radius", "radunits", "rdc", "rstrand", "runits"
]

const _tsdata_properties = String[
    "diacable", "diains", "diam", "diashield", "emergamps", "epsr", "gmrac",
    "gmrunits", "inslayer", "like", "normamps", "rac", "radius", "radunits",
    "rdc", "runits", "taplap", "taplayer"
]

const _wiredata_properties = String[
    "rdc", "rac", "runits", "gmrac", "gmrunits", "radius", "radunits",
    "normamps", "emergamps", "diam", "like"
]

const _xfmrcode_properties = String[
    "phases", "windings", "wdg", "conn", "kv", "kva", "tap", "%r", "rneut",
    "xneut", "conns", "kvs", "kvas", "taps", "%rs", "xhl", "xlt", "xht",
    "xscarray", "thermal", "n", "m", "flrise", "hsrise", "%loadloss",
    "%noloadloss", "%imag", "ppm_antifloat", "normhkva", "emerghkva", "sub",
    "maxtap", "mintap", "numtaps", "subname", "xrconst", "leadlag",
    "wdgcurrents", "core", "rdcohms", "faultrate", "basefreq", "enabled", "like"
]

const _spectrum_properties = String[
    "numharm", "harmonic", "%mag", "angle", "csvfile"
]

const _vsource_properties = String[
    "bus1", "bus2", "basekv", "pu", "angle", "frequency", "phases", "mvasc3",
    "mvasc1", "x1r1", "x0r0", "isc3", "isc1", "r1", "x1", "r0", "x0",
    "scantype", "sequence", "spectrum", "z1", "z2", "z0", "puz1", "puz2",
    "puz0", "basemva", "basefreq", "like", "enabled"
]

const _isource_properties = String[
    "phases", "bus1", "amps", "angle", "frequency", "scantype", "sequence",
    "spectrum", "basefreq", "enabled", "like"
]

const _fault_properties = String[
    "phases", "bus1", "bus2", "r", "gmatrix", "minamps", "ontime", "pctperm",
    "temporary", "%stddev", "normamps", "emergamps", "basefreq", "faultrate",
    "repair", "enabled", "like"
]

const _capacitor_properties = String[
    "bus1", "bus2", "phases", "kvar", "kv", "conn", "cmatrix", "cuf", "r",
    "xl", "harm", "numsteps", "states", "normamps", "emergamps", "faultrate",
    "pctperm", "basefreq", "enabled", "like"
]

const _line_properties = String[
    "bus1", "bus2", "linecode", "length", "phases", "r1", "x1", "r0", "x0",
    "c1", "c0", "normamps", "emergamps", "faultrate", "pctperm",
    "repair", "basefreq", "rmatrix", "xmatrix", "cmatrix", "switch", "rg",
    "xg", "rho", "geometry", "units", "spacing", "wires", "earthmodel",
    "cncables", "tscables", "b1", "b0", "seasons", "ratings", "enabled", "like"
]

const _reactor_properties = String[
    "phases", "bus1", "bus2", "kv", "kvar", "conn", "parallel", "r", "rmatrix",
    "rp", "x", "xmatrix", "z", "z1", "z2", "z0", "rcurve", "lcurve", "lmh",
    "normamps", "emergamps", "repair", "faultrate", "pctperm", "basefreq",
    "enabled", "like"
]

const _transformer_properties = String[
    "phases", "windings", "wdg", "bus", "conn", "kv", "kva", "tap", "%r",
    "rneut", "xneut", "buses", "conns", "kvs", "kvas", "taps", "%rs", "xhl",
    "xlt", "xht", "xscarray", "thermal", "n", "m", "flrise", "hsrise",
    "%loadloss", "%noloadloss", "%imag", "ppm_antifloat", "normhkva",
    "emerghkva", "sub", "maxtap", "mintap", "numtaps", "subname", "bank",
    "xfmrcode", "xrconst", "leadlag", "wdgcurrents", "core", "rdcohms",
    "faultrate", "basefreq", "enabled", "like"
]

const _gictransformer_properties = String[
    "basefreq", "bush", "busnh", "busnx", "busx", "emergamps", "enabled",
    "phases", "r1", "r2", "type", "mva", "kvll1", "kvll2", "%r1", "%r2", "k",
    "varcurve", "like", "normamps", "emergamps", "pctperm", "repair"
]

const _gicline_properties = String[
    "angle", "bus1", "bus2", "c", "ee", "en", "frequency", "lat1", "lat2",
    "lon1", "lon2", "phases", "r", "volts", "x", "like", "basefreq",
    "enabled", "spectrum"
]

const _load_properties = String[
    "phases", "bus1", "kv", "kw", "pf", "model", "yearly", "daily", "duty",
    "growth", "conn", "kvar", "rneut", "xneut", "status", "class", "vminpu",
    "vmaxpu", "vminnorm", "vminemerg", "xfkva", "allocationfactor", "kva",
    "%mean", "%stddev", "cvrwatts", "cvrvars", "kwh", "kwhdays", "cfactor",
    "cvrcurve", "numcust", "spectrum", "zipv", "%seriesrl", "relweight",
    "vlowpu", "puxharm", "xrharm", "spectrum", "basefreq", "enabled", "like"
]

const _generator_properties = String[
    "bus1", "phases", "kv", "kw", "pf", "model", "yearly", "daily", "duty",
    "dispvalue", "conn", "kvar", "rneut", "xneut", "status", "class", "vpu",
    "maxkvar", "minkvar", "pvfactor", "debugtrace", "vminpu", "vmaxpu",
    "forceon", "kva", "mva", "xd", "xdp", "xdpp", "h", "d", "usermodel",
    "userdata", "shaftmodel", "shaftdata", "dutystart", "balanced", "xrdp",
    "spectrum", "basefreq", "enabled", "like"
]

const _indmach012_properties = String[
    "phases", "bus1", "kv", "kw", "pf", "conn", "kva", "h", "d", "purs",
    "puxs", "purr", "puxr", "puxm", "slip", "maxslip", "slipoption",
    "spectrum", "enabled"
]

const _storage_properties = String[
    "phases", "bus1", "%charge", "%discharge", "%effcharge", "%idlingkvar",
    "idlingkw", "%r", "%reserve", "%stored", "%x", "basefreq", "chargetrigger",
    "class", "conn", "daily", "yearly", "debugtrace", "dischargetrigger",
    "dispmode", "duty", "dynadata", "dynadll", "enabled", "kv", "kva", "kvar",
    "kw", "kwhrated", "kwhstored", "kwrated", "like", "model", "pf",
    "spectrum", "state", "timechargetrig", "userdata", "usermodel", "vmaxpu",
    "vminpu", "yearly"
]

const _capcontrol_properties = String[
    "element", "capacitor", "type", "ctphase", "ctratio", "deadtime", "delay",
    "delayoff", "eventlog", "offsetting", "onsetting", "ptphase", "ptratio",
    "terminal", "vbus", "vmax", "vmin", "voltoverride", "enabled"
]

const _regcontrol_properties = String[
    "transformer", "winding", "vreg", "band", "delay", "ptratio", "ctprim",
    "r", "x", "pthase", "tapwinding", "bus", "remoteptratio", "debugtrace",
    "eventlog", "inversetime", "maxtapchange", "revband", "revdelay",
    "reversible", "revneutral", "revr", "revthreshold", "revvreg", "revx",
    "tapdelay", "tapnum", "vlimit", "ldc_z", "rev_z", "cogen", "enabled"
]

const _energymeter_properties = String[
    "element", "terminal", "action", "clear", "save", "take", "option",
    "kwnorm", "kwemerg", "peakcurrent", "zonelist", "zonelist", "zonelist",
    "localonly", "mask", "losses", "linelosses", "xfmrlosses", "seqlosses",
    "3phaselosses", "vbaselosses", "basefreq", "enabled", "like"
]

const _monitor_properties = String[
    "element", "terminal", "mode", "action"
]

const _pvsystem_properties = Vector{String}([
    "phases", "bus1", "kv", "irradiance", "pmpp", "%pmpp", "temperature", "pf",
    "conn", "kvar", "kva", "%cutin", "%cutout", "effcurve", "p-tcurve", "%r",
    "%x", "model", "vminpu", "vmaxpu", "yearly", "daily", "duty", "tyearly",
    "tduty", "class", "usermodel", "userdata", "debugtrace", "varfollowinverter",
    "dutystart", "wattpriority", "pfpriority", "%pminnovars", "%pminkvarmax",
    "kvarmax", "kvarmaxabs", "spectrum", "basefreq", "enabled", "like"
])

const _recloser_properties = String[
    "monitoredobj", "monitoredterm", "switchedobj", "switchedterm", "numfast",
    "phasefast", "phasedelayed", "groundfast", "grounddelayed", "phasetrip",
    "groundtrip", "phaseinst", "groundinst", "reset", "shots",
    "recloseintervals", "delay", "action", "tdphfast", "tdgrfast",
    "tdphdelayed", "tdgrdelayed", "basefreq", "enabled", "like"
]

const _relay_properties = String[
    "monitoredobj", "monitoredterm", "switchedobj", "switchedterm", "type",
    "phasecurve", "groundcurve", "phasetrip", "groundtrip", "tdphase",
    "tdground", "phaseinst", "groundinst", "reset", "shots",
    "recloseintervals", "delay", "overvoltcurve", "undervoltcurve",
    "kvbase", "47%pickup", "46baseamps", "46%pickup", "46isqt",
    "variable", "overtrip", "undertrip", "breakertime", "action", "basefreq",
    "enabled"
]

const _fuse_properties = String[
    "monitoredobj", "monitoredterm", "switchedobj", "switchedterm",
    "fusecurve", "ratedcurrent", "delay", "action", "basefreq", "enabled"
]

const _swtcontrol_properties = String[
    "action", "basefreq", "delay", "enabled", "like", "lock", "normal",
    "reset", "state", "switchedobj", "switchedterm"
]

const _dss_object_properties = Dict{String,Vector{String}}(
    "linecode" => _linecode_properties,
    "linegeometry" => _linegeometry_properties,
    "linespacing" => _linespacing_properties,
    "loadshape" => _loadshape_properties,
    "xycurve" => _xycurve_properties,
    "growthshape" => _growthshape_properties,
    "tcc_curve" => _tcc_curve_properties,
    "cndata" => _cndata_properties,
    "tsdata" => _tsdata_properties,
    "wiredata" => _wiredata_properties,
    "xfmrcode" => _xfmrcode_properties,
    "spectrum" => _spectrum_properties,
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
    "monitor" => _monitor_properties,
    "pvsystem" => _pvsystem_properties,
    "relay" => _relay_properties,
    "recloser" => _reactor_properties,
    "fuse" => _fuse_properties,
    "swtcontrol" => _swtcontrol_properties,
)


"parses single column load profile files"
function _parse_csv_file(path::FilePaths.AbstractPath, type::AbstractString; header::Bool=false, column::Int=1, interval::Bool=false)::Union{Vector{String}, Tuple{Vector{String}, Vector{String}, Vector{String}}, Tuple{Vector{String}, Vector{String}}}
    open(first(Glob.glob([Glob.FilenameMatch(basename(path), "i")], dirname(path))), "r") do f
        lines = readlines(f)
        if header
            lines = lines[2:end]
        end

        if type == "mult"
            return Vector{String}([split(line, ",")[column] for line in lines])
        elseif type == "csvfile"
            if interval
                hour, mult = String[], String[]
                for line in lines
                    d = split(line, ",")
                    push!(hour, string(d[1]))
                    push!(mult, string(d[2]))
                end
                return hour, mult
            else
                return Vector{String}([split(line, ",")[1] for line in lines])
            end
        elseif type == "pqcsvfile"
            if interval
                hour, pmult, qmult = String[], String[], String[]
                for line in lines
                    d = split(line, ",")
                    push!(hour, string(d[1]))
                    push!(pmult, string(d[2]))
                    push!(qmult, string(d[3]))
                end

                return hour, pmult, qmult
            else
                pmult, qmult = String[], String[]
                for line in lines
                    d = split(line, ",")
                    push!(pmult, string(d[1]))
                    push!(qmult, string(d[2]))
                end

                return pmult, qmult
            end
        end
    end
end


"parses sng and dbl precision loadshape binary files"
function _parse_binary_file(path::FilePaths.AbstractPath, precision::Type; npts::Union{Int,Nothing}=nothing, interval::Bool=false)::Union{Vector{precision}, Tuple{Vector{precision}, Vector{precision}}}
    open(first(Glob.glob([Glob.FilenameMatch(basename(path), "i")], dirname(path))), "r") do f
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
        else
            data = Array{precision, 1}(undef, interval ? npts * 2 : npts)

            try
                read!(f, data)
            catch EOFError
                error("Error reading binary file: likely npts is wrong")
            end
        end

        if interval
            data = reshape(data, 2, :)
            return data[1, :], data[2, :]
        else
            return data
        end
    end
end


"parses csv or binary loadshape files"
function _parse_data_file(path::FilePaths.AbstractPath, type::AbstractString, npts::Union{Int,Nothing}; header::Bool=false, interval::Bool=false, column::Int=1)
    if type in ["csvfile", "mult", "pqcsvfile"]
        return _parse_csv_file(path, type; header=header, column=column, interval=interval)
    elseif type in ["sngfile", "dblfile"]
        return _parse_binary_file(path, Dict("sngfile" => Float32, "dblfile" => Float64)[type]; npts=npts, interval=interval)
    end
end


"parses pmult and qmult entries on loadshapes"
function _parse_mult_parameter(mult_string::AbstractString; path::FilePaths.AbstractPath=FilePaths.p".", npts::Union{Int,Nothing}=nothing)::String
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
        full_path = joinpath(path, props[file_key])
        type = file_key == "file" ? "mult" : file_key

        return "($(join(_parse_data_file(full_path, type, npts; header=get(props, "header", false), column=parse(Int, get(props, "column", "1"))), ",")))"
    end
end


"parses loadshape component"
function _parse_loadshape!(current_obj::Dict{String,<:Any}; path::FilePaths.AbstractPath=FilePaths.p".")
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
            full_path = joinpath(path, current_obj[prop])
            data = _parse_data_file(full_path, prop, parse(Int, get(current_obj, "npts", "1")); interval=interval, header=false)
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

    _clean_arrays!(current_obj, ["pmult", "qmult", "hour"])
end


"parse xycurve component"
function _parse_xycurve!(current_obj::Dict{String,<:Any}; path::FilePaths.AbstractPath=FilePaths.p".")
    for prop in current_obj["prop_order"]
        if prop in ["csvfile", "sngfile", "dblfile"]
            full_path = joinpath(path, current_obj[prop])
            data = _parse_data_file(full_path, prop, nothing; interval=true, header=false)
            current_obj["xarray"], current_obj["yarray"] = data
        end

    end

    _clean_arrays!(current_obj, ["xarray", "yarray", "points"])
end


"parse spectrum component"
function _parse_spectrum!(current_obj::Dict{String,<:Any}; path::FilePaths.AbstractPath=FilePaths.p".")
    for prop in current_obj["prop_order"]
        if prop == "csvfile"
            full_path = joinpath(path, current_obj[prop])
            current_obj["harmonic"], current_obj["%mag"], current_obj["angle"] = _parse_csv_file(full_path, "pqcsvfile"; header=false, interval=true)
        end
    end

    _clean_arrays!(current_obj, ["harmonic", "%mag", "angle"])

end


"cleans up array properties back into strings for later conversion"
function _clean_arrays!(current_obj::Dict{String,<:Any}, properties::Vector{String})
    for prop in properties
        if haskey(current_obj, prop) && isa(current_obj[prop], Array)
            current_obj[prop] = "($(join(current_obj[prop], ",")))"
        elseif haskey(current_obj, prop) && isa(current_obj[prop], String) && !_isa_array(current_obj[prop])
            current_obj[prop] = "($(current_obj[prop]))"
        end
    end
end

"strips lines that are either commented (block or single) or empty"
function _strip_lines(lines::Vector{<:AbstractString})::Vector{String}
    blockComment = false
    stripped_lines = String[]
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
function _parse_buscoords_file(path::FilePaths.AbstractPath)::Dict{String,Any}
    file_str = read(open(first(Glob.glob([Glob.FilenameMatch(basename(path), "i")], dirname(path)))), String)
    regex = r"[\s,\t]+"

    lines = _strip_lines(split(file_str, '\n'))

    buscoords = Dict{String,Dict{String,Any}}()
    for line in lines
        line = _strip_comments(line)

        bus, x, y = split(line, regex; limit=3)
        buscoords[lowercase(strip(bus, [',']))] = Dict{String,Any}("x"=>parse(Float64, strip(x, [','])), "y"=>parse(Float64, strip(y, [',', '\r'])))
    end

    return buscoords
end


"parses the setbusxy command"
function _parse_setbusxy!(data_dss::Dict{String,<:Any}, properties_str::AbstractString)
    properties = _parse_properties(properties_str)
    object = Dict{String,Any}()
    for (prop_name, prop) in zip(["bus", "x", "y"], properties)
        if occursin("=", prop)
            prop_name, prop = split(prop, "=")
        end
        object[prop_name] = prop
    end
    bus = pop!(object, "bus")
    for (k,v) in object
        object[k] = parse(Float64, v)
    end

    if !haskey(data_dss, "buscoords")
        data_dss["buscoords"] = Dict{String,Any}(bus => object)
    else
        data_dss["buscoords"][bus] = object
    end
end


"""
Parses a string of `properties` of a component type, character by character
into an array with each element containing (if present) the property name, "=",
and the property value.
"""
function _parse_properties(properties::AbstractString)::Vector{String}
    parsed_properties = Vector{SubString{String}}([])
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
        elseif !occursin("=", string_out) && (isspace(char) || n == num_chars)
            end_equality = true
        elseif occursin("=", string_out) && (isspace(char) || n == num_chars) && end_array
            end_equality = true
        else
            end_equality = false
        end

        if isspace(char)
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
function _add_component!(data_dss::Dict{String,<:Any}, obj_type_name::AbstractString, object::Dict{String,<:Any})
    obj_type = split(obj_type_name, '.'; limit=2)[1]
    if obj_type == "circuit"
        if haskey(data_dss, "circuit")
            error("Cannot have two circuits, invalid dss")
        else
            data_dss[obj_type] = object
        end
    elseif haskey(data_dss, obj_type)
        data_dss[obj_type][object["name"]] = object
    else
        data_dss[obj_type] = Dict{String,Any}(object["name"] => object)
    end
end


function _add_component_edits!(data_dss::Dict{String,<:Any}, obj_type_name::AbstractString, object::Dict{String,<:Any})
    obj_type = split(obj_type_name, '.'; limit=2)[1]
    if !haskey(data_dss, obj_type)
        data_dss[obj_type] = Dict{String,Any}(
            object["name"] => object
        )
    else
        if !haskey(data_dss[obj_type], object["name"])
            data_dss[obj_type][object["name"]] = object
        else
            filtered_prop_order = filter(x->x!="name", object["prop_order"])
            object["prop_order"] = vcat(filter(x->!(x in filtered_prop_order), data_dss[obj_type][object["name"]]["prop_order"]), filtered_prop_order)
            merge!(data_dss[obj_type][object["name"]], object)
        end
    end
end


"""
Adds a property to an existing component properties dictionary `object` given
the `key` and `value` of the property. If a property of the same name already
exists inside `object`, the original value is converted to an array, and the
new value is appended to the end.
"""
function _add_property(object::Dict{String,<:Any}, key::SubString{String}, value::Any)::Dict{String,Any}
    if !haskey(object, "prop_order")
        object["prop_order"] = Vector{String}(["name"])
    end

    current_wdg = key == "wdg" ? value == "1" ? "" : "$value" : any(occursin("wdg", prop) for prop in object["prop_order"]) ? replace(split(filter(x->occursin("wdg", x), object["prop_order"])[end], "_")[end], "wdg"=>"") : ""

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
    push!(object["prop_order"], string(lowercase(key)))

    return object
end


"""
Parses a `component` with `properties` into a `object`. If `object` is not
defined, an empty dictionary will be used. Assumes that unnamed properties are
given in order, but named properties can be given anywhere.
"""
function _parse_component(component::AbstractString, properties::AbstractString, object::Dict{String,<:Any}=Dict{String,Any}(); path::FilePaths.AbstractPath=FilePaths.p".")::Dict{String,Any}
    obj_type, name = split(component, '.'; limit=2)

    if !haskey(object, "prop_order")
        object["prop_order"] = String[]
    end

    if !haskey(object, "name")
        object["name"] = string(name)
    end

    push!(object["prop_order"], "name")

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
        if k in keys(dss_prime) && (isa(v, Dict) || isa(v, Set))
            if isa(v, Dict)
                if haskey(v, "prop_order")
                    v["prop_order"] = vcat(get(dss_prime[k], "prop_order", []), filter(x->x!="name", v["prop_order"]))
                end
                _merge_dss!(dss_prime[k], v)
            elseif isa(v, Set)
                union!(dss_prime[k], v)
            end
        else
            dss_prime[k] = v
        end
    end
end


"""
Parses an already separated line given by `elements` (an array) of an OpenDSS
file into `current_obj`. If not defined, `current_obj` is an empty dictionary.
"""
function _parse_line(elements::Vector{String}; current_obj::Dict{String,<:Any}=Dict{String,Any}(), path::FilePaths.AbstractPath=FilePaths.p".")::Tuple{SubString{String}, Dict{String,Any}}
    current_obj_type = strip(elements[2], ['\"', '\''])
    if startswith(current_obj_type, "object")
        current_obj_type = split(current_obj_type, '=')[2]
        current_obj["name"] = split(current_obj_type, '.')[2]
    end

    if length(elements) != 3
        properties = ""
    else
        properties = elements[3]
    end

    current_obj = _parse_component(current_obj_type, properties; path=path)

    return current_obj_type, current_obj
end


"Strips comments, defined by '!' from the ends of lines"
function _strip_comments(line::AbstractString)::String
    return strip(split(line, r"\s*!")[1], ['\r', '\n'])
end


"""
Assigns a property with name `property_name` and value `property_value` to the component
of type `obj_type` named `obj_name` in `data_dss`.
"""
function _assign_property!(data_dss::Dict{String,<:Any}, obj_type::AbstractString, obj_name::AbstractString, property_name::AbstractString, property_value::Any)
    if haskey(data_dss, obj_type) && haskey(data_dss[obj_type], obj_name)
        data_dss[obj_type][obj_name][property_name] = property_value
    else
        @warn "Cannot find $obj_type object $obj_name."
    end
end


"""
    parse_dss(filename::String; data_dss::Union{Missing,Dict{String,Any}}=missing)::Dict{String,Any}

Parses a OpenDSS file given by `filename` into a Dict{Array{Dict}}. Only
supports components and options, but not commands, e.g. "plot" or "solve".
Will also parse files defined inside of the originating DSS file via the
"compile", "redirect" or "buscoords" commands.
"""
function parse_dss(path::Union{AbstractString,FilePaths.AbstractPath}; data_dss::Union{Missing,Dict{String,Any}}=missing)::Dict{String,Any}
    path = isa(path, FilePaths.AbstractPath) ? path : FilePaths.Path(path)
    data_dss = open(first(Glob.glob([Glob.FilenameMatch(basename(path), "i")], dirname(path)))) do io
        parse_dss(io; data_dss=data_dss)
    end
    return data_dss
end


"""
    parse_dss(io::IO; data_dss::Union{Missing,Dict{String,Any}}=missing)::Dict{String,Any}

Parses a OpenDSS file aleady in IO into a Dict{Array{Dict}}. Only
supports components and options, but not commands, e.g. "plot" or "solve".
Will also parse files defined inside of the originating DSS file via the
"compile", "redirect" or "buscoords" commands.
"""
function parse_dss(io::IO; data_dss::Union{Missing,Dict{String,Any}}=missing)::Dict{String,Any}
    filename = isa(io, IOStream) ? match(r"^<file\s(.+)>$", io.name).captures[1] : "GenericIOBuffer"
    current_file = basename(FilePaths.Path(filename))
    path = dirname(FilePaths.Path(filename))
    data_dss = ismissing(data_dss) ? Dict{String,Any}() : data_dss

    data_dss["filename"] = haskey(data_dss, "filename") ? union(data_dss["filename"], Set{String}([string(filename)])) : Set{String}([string(filename)])

    current_obj = Dict{String,Any}()
    current_obj_type = ""

    lines = readlines(io)

    stripped_lines = _strip_lines(lines)
    nlines = length(stripped_lines)

    for (n, line) in enumerate(stripped_lines)
        real_line_num = findall(lines .== line)[1]
        line = _strip_comments(line)

        if startswith(strip(line), '~') || startswith(strip(lowercase(line)), "more")
            if startswith(strip(lowercase(line)), "more")
                line = lowercase(strip(line)[5:end])
            end

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

            if cmd in _dss_unsupported_commands
                @info "Command '$cmd' on line $real_line_num in '$current_file' is not supported, skipping."
                continue

            elseif cmd == "clear"
                @info "Circuit has been reset with the 'clear' on line $real_line_num in '$current_file'"
                data_dss = Dict{String,Any}("filename"=>data_dss["filename"])
                continue

            elseif cmd in ["redirect", "compile"]
                file_path = split(strip(line_elements[2], ['(',')']), r"[\\|\/]")

                if !(joinpath(file_path...) in data_dss["filename"])
                    full_path = joinpath(path, file_path...)
                    @info "Redirecting to '$(joinpath(file_path...))' on line $real_line_num in '$current_file'"
                    data_dss = parse_dss(full_path; data_dss=data_dss)
                end

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

            elseif cmd == "edit"
                current_obj_type, current_obj = _parse_line([lowercase(line_element) for line_element in line_elements]; path=path)
                if startswith(current_obj_type, "circuit")
                    current_obj_type = "vsource.source"
                    current_obj["name"] = "source"
                end

                _add_component_edits!(data_dss, current_obj_type, current_obj)
                continue

            elseif cmd in ["disable", "enable"]
                current_obj_type, current_obj_name = split(join(line_elements[2:end], ""), ".")

                enabled = cmd == "enable" ? "true" : "false"

                _add_component_edits!(data_dss, current_obj_type, Dict{String,Any}("name"=>current_obj_name, "enabled"=>enabled))
                continue

            elseif cmd in ["buscoords", "latloncoords"]
                file = FilePaths.Path(line_elements[2])
                full_path = joinpath(path, file)
                @info "Reading Buscoords in '$file' on line $real_line_num in '$current_file'"
                data_dss["buscoords"] = _parse_buscoords_file(full_path)

            elseif cmd == "setbusxy"
                _parse_setbusxy!(data_dss, join(line_elements[2:end], " "))
                continue

            elseif cmd == "new"
                current_obj_type, current_obj = _parse_line([lowercase(line_element) for line_element in line_elements]; path=path)

                if startswith(current_obj_type, "loadshape")
                    _parse_loadshape!(current_obj; path=path)
                elseif startswith(current_obj_type, "xycurve")
                    _parse_xycurve!(current_obj; path=path)
                elseif startswith(current_obj_type, "spectrum")
                    _parse_spectrum!(current_obj; path=path)
                end

                if startswith(current_obj_type, "circuit")
                    current_obj_type = "vsource.source"

                    data_dss["circuit"] = current_obj["name"]

                    current_obj["name"] = "source"
                end
            elseif split(cmd, '.')[1] in keys(_dss_object_properties)
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
                    end
                    _assign_property!(data_dss, obj_type, obj_name, property_name, property_value)
                end
            else
                @warn "Command '$cmd' on line $real_line_num in '$current_file' is not recognized, skipping."
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

    _parse_dss_with_dtypes!(data_dss)

    data_dss["data_model"] = DSS

    return data_dss
end


"""
    parse_dss_voltages_export(file::String)::Dict{String,Any}

Parses a voltages CSV file exported from OpenDSS with the command `Export Voltages [filename.csv]`
"""
function parse_dss_voltages_export(path::String)::Dict{String,Any}
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
