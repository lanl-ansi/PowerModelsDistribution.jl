"""
    parse_file(file)

Parses a matlab .m `file` into a Three Phase PowerModels data structure.
"""
function parse_file(file::String)
    if endswith(file, ".m")
        tppm_data = ThreePhasePowerModels.parse_matlab(file)
    else
        error(LOGGER, "only .m files are supported")
    end

    check_network_data(tppm_data)

    return tppm_data
end


""
function check_network_data(data::Dict{String,Any})
    data["version"] = Pkg.installed("ThreePhasePowerModels")
    PMs.make_per_unit(data)

    # TODO see which of these can be reused
    #=
    PMs.check_connectivity(data)
    PMs.check_transformer_parameters(data)
    PMs.check_voltage_angle_differences(data)
    PMs.check_thermal_limits(data)
    PMs.check_branch_directions(data)
    PMs.check_branch_loops(data)
    PMs.check_bus_types(data)
    PMs.check_dcline_limits(data)
    PMs.check_voltage_setpoints(data)
    PMs.check_cost_functions(data)
    =#
end