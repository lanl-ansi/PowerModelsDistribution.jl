"""
    parse_file(io)

Parses the IOStream of a file into a Three-Phase PowerModels data structure.
"""
function parse_file(io::IO; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1, filetype="json")
    if filetype == "m"
        tppm_data = ThreePhasePowerModels.parse_matlab(io)
    elseif filetype == "dss"
        Memento.warn(LOGGER, "Not all OpenDSS features are supported, currently only minimal support for lines, loads, generators, and capacitors as shunts. Transformers and reactors as transformer branches are included, but value translation is not fully supported.")
        tppm_data = ThreePhasePowerModels.parse_opendss(io; import_all=import_all, vmin=vmin, vmax=vmax)
    elseif filetype == "json"
        tppm_data = PowerModels.parse_json(io; validate=false)
    else
        Memento.error(LOGGER, "only .m and .dss files are supported")
    end

    check_network_data(tppm_data)

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
function check_network_data(data::Dict{String,Any})
    PMs.make_per_unit(data)

    PMs.check_connectivity(data)
    PMs.check_transformer_parameters(data)
    PMs.check_voltage_angle_differences(data)
    PMs.check_thermal_limits(data)
    PMs.check_branch_directions(data)
    PMs.check_branch_loops(data)
    PMs.check_bus_types(data)
    PMs.check_dcline_limits(data)
    # PMs.check_voltage_setpoints(data)
    PMs.check_cost_functions(data)
    PMs.standardize_cost_terms(data)
end

function wrapto180(degrees)
    return degrees - 360*floor.((degrees .+ 180)/360)
end

function wraptopi(radians)
    return radians - 2*pi*floor.((radians .+ pi)/(2*pi))
end
