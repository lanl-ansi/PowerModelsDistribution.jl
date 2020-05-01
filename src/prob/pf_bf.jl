"Branch Flow Power Flow Problem"
function run_mc_pf_bf(data::Union{Dict{String,<:Any},String}, model_type::DataType, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mc_pf_bf; kwargs...)
end


"Constructor for Branch Flow Power Flow"
function build_mc_pf_bf(pm::_PM.AbstractPowerModel)
    # Variables
    variable_mc_bus_voltage(pm; bounded=false)
    variable_mc_branch_current(pm)
    variable_mc_branch_power(pm)
    variable_mc_gen_power_setpoint(pm; bounded=false)

    # Constraints
    constraint_mc_model_current(pm)

    for (i,bus) in ref(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)

        @assert bus["bus_type"] == 3
        constraint_mc_voltage_magnitude_only(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_mc_power_balance(pm, i)

        # PV Bus Constraints
        if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            constraint_mc_voltage_magnitude_only(pm, i)
            for j in ref(pm, :bus_gens, i)
                constraint_mc_gen_power_setpoint_real(pm, j)
            end
        end
    end

    for i in ids(pm, :branch)
        constraint_mc_power_losses(pm, i)
        constraint_mc_model_voltage_magnitude_difference(pm, i)
        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    # Objective
    _PM.objective_min_fuel_cost(pm)
end
