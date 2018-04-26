"""
    parse_file(file)

Parses a matlab .m `file` into a Three Phase PowerModels data structure.
"""
function parse_file(file::String)
    if endswith(file, ".m")
        tppm_data = ThreePhasePowerModels.parse_matlab(file)
    elseif endswith(lowercase(file), ".dss")
        warn(LOGGER, "Not all OpenDSS features are supported, currently only raw data is loaded.")
        tppm_data = ThreePhasePowerModels.parse_opendss(file)
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

    PMs.check_voltage_angle_differences(data)
    PMs.check_connectivity(data)
    PMs.check_bus_types(data)
    # PMs.check_voltage_setpoints(data)  # CHECK: no "vg" in gens?
    PMs.check_transformer_parameters(data)
    PMs.check_dcline_limits(data)
    PMs.check_cost_functions(data)
    PMs.check_branch_loops(data)
    PMs.check_branch_directions(data)
    PMs.check_thermal_limits(data)
end
