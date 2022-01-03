const _dss_supported_commands = String[
    "buscoords", "clear", "clearall", "clone", "close",
    "compile", "redirect", "disable", "edit", "enable", "get",
    "more", "m", "new", "open", "remove", "set", "setbusxy"
]

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

const _dss_cmd_new_regex = r"(([\w%_]+)=){0,1}(([\[\(\{\"\'](?<=[\[\(\{\"\'])(.+?)(?=[\]\)\}\"\'])[\]\)\}\"\'])|(.+?))(?:\s+|$)"
const _dss_cmd_more_regex = _dss_cmd_new_regex
const _dss_cmd_set_regex = _dss_cmd_new_regex
const _dss_cmd_buscoords_regex = r"[\s,\t]+"
const _dss_matrix_regex = r"\s*,\s*|\s+"
const _dss_array_regex = r"\s*,\s*|\s+"
const _dss_rpn_array_sep_regex = Regex(string("[",join(_array_delimiters, '\\'),"]"))

const _dss_property_renames = Dict{String,Dict{String,String}}(
    "loadshape" => Dict{String,String}(
        "mult" => "pmult",
    ),
    "transformer" => Dict{String,String}(
        "ppm" => "ppm_antifloat",
        "x12" => "xhl",
        "x23" => "xlt",
        "x13" => "xht",
    ),
    "linecode" => Dict{String,String}(
        "phases" => "nphases",
    )
)
