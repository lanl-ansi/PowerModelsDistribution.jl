# This problem includes load models beyond simple constant power ones; this is
# handled in variable_mc_load_setpoint and constraint_mc_load_setpoint.


"branch flow with loadmodels"
function run_mc_opf_bf_lm(data::Union{Dict{String,<:Any},String}, model_type::DataType, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mc_opf_bf_lm; kwargs...)
end


"constructor for branch flow with loadmodels"
function build_mc_opf_bf_lm(pm::_PM.AbstractPowerModel)
    # Variables
    variable_mc_bus_voltage(pm)
    variable_mc_branch_current(pm)
    variable_mc_branch_power(pm)

    variable_mc_gen_power_setpoint(pm)
    variable_mc_load_setpoint(pm)

    # Constraints
    constraint_mc_model_current(pm)

    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_power_losses(pm, i)
        constraint_mc_model_voltage_magnitude_difference(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :load)
        constraint_mc_load_setpoint(pm, i)
    end

    for i in ids(pm, :gen)
        constraint_mc_gen_setpoint(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_mc_power_balance(pm, i)
    end

    # Objective
    _PM.objective_min_fuel_cost(pm)
end
