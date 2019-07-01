"""
    parse_file(io)

Parses the IOStream of a file into a Three-Phase PowerModels data structure.
"""
function parse_file(io::IO; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1, filetype="json")
    if filetype == "m"
        pmd_data = PowerModelsDistribution.parse_matlab(io)
    elseif filetype == "dss"
        Memento.warn(_LOGGER, "Not all OpenDSS features are supported, currently only minimal support for lines, loads, generators, and capacitors as shunts. Transformers and reactors as transformer branches are included, but value translation is not fully supported.")
        tppm_data = PowerModelsDistribution.parse_opendss(io; import_all=import_all, vmin=vmin, vmax=vmax)
    elseif filetype == "json"
        pmd_data = PowerModels.parse_json(io; validate=false)
    else
        Memento.error(_LOGGER, "only .m and .dss files are supported")
    end

    correct_network_data!(pmd_data)

    return pmd_data
end


""
function parse_file(file::String; kwargs...)
    pmd_data = open(file) do io
        parse_file(io; filetype=split(lowercase(file), '.')[end], kwargs...)
    end
    return pmd_data
end


""
function correct_network_data!(data::Dict{String,Any})
    _PMs.make_per_unit!(data)

    _PMs.check_connectivity(data)
    _PMs.correct_transformer_parameters!(data)
    _PMs.correct_voltage_angle_differences!(data)
    _PMs.correct_thermal_limits!(data)
    _PMs.correct_branch_directions!(data)
    _PMs.check_branch_loops(data)
    _PMs.correct_bus_types!(data)
    _PMs.correct_dcline_limits!(data)
    _PMs.correct_cost_functions!(data)
    _PMs.standardize_cost_terms!(data)
end
