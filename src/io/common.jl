"""
    parse_file(io)

Parses the IOStream of a file into a Three-Phase PowerModels data structure.
"""
function parse_file(io::IO; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1, filetype="json")
    if filetype == "m"
        tppm_data = ThreePhasePowerModels.parse_matlab(io)
    elseif filetype == "dss"
        Memento.warn(_LOGGER, "Not all OpenDSS features are supported, currently only minimal support for lines, loads, generators, and capacitors as shunts. Transformers and reactors as transformer branches are included, but value translation is not fully supported.")
        tppm_data = ThreePhasePowerModels.parse_opendss(io; import_all=import_all, vmin=vmin, vmax=vmax)
    elseif filetype == "json"
        tppm_data = PowerModels.parse_json(io; validate=false)
    else
        Memento.error(_LOGGER, "only .m and .dss files are supported")
    end

    correct_network_data!(tppm_data)

    return tppm_data
end


""
function parse_file(file::String; kwargs...)
    tppm_data = open(file) do io
        parse_file(io; filetype=split(lowercase(file), '.')[end], kwargs...)
    end
    return tppm_data
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
    # _PMs.check_voltage_setpoints(data)
    _PMs.correct_cost_functions!(data)
    _PMs.standardize_cost_terms!(data)
end

function _wrap_to_180(degrees)
    return degrees - 360*floor.((degrees .+ 180)/360)
end

function _wrap_to_pi(radians)
    return radians - 2*pi*floor.((radians .+ pi)/(2*pi))
end
