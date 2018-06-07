"""
    parse_file(file)

Parses a matlab .m `file` into a Three Phase PowerModels data structure.
"""
function parse_file(file::String; import_all::Bool=false, vmin::Float64=0.9, vmax::Float64=1.1)
    if endswith(file, ".m")
        tppm_data = ThreePhasePowerModels.parse_matlab(file)
    elseif endswith(lowercase(file), ".dss")
        warn(LOGGER, "Not all OpenDSS features are supported, currently only minimal support for \
                      lines, loads, generators, and capacitors as shunts. Transformers and reactors \
                      as transformer branches are included, but value translation is not fully supported.")
        tppm_data = ThreePhasePowerModels.parse_opendss(file; import_all=import_all, vmin=vmin, vmax=vmax)
    else
        error(LOGGER, "only .m and .dss files are supported")
    end

    check_network_data(tppm_data)

    return tppm_data
end


""
function check_network_data(data::Dict{String,Any})
    data["version"] = Pkg.installed("ThreePhasePowerModels")
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
end
