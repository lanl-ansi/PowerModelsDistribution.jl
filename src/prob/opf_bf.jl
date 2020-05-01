"branch flow opf"
function run_mc_opf_bf(data::Union{Dict{String,<:Any},String}, model_type::DataType, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mc_opf_bf; kwargs...)
end


"constructor for branch flow opf"
function build_mc_opf_bf(pm::_PM.AbstractPowerModel)
    # Variables
    variable_mc_bus_voltage(pm)
    variable_mc_branch_current(pm)
    variable_mc_branch_power(pm)
    variable_mc_transformer_power(pm)
    variable_mc_gen_power_setpoint(pm)

    # Constraints
    constraint_mc_model_current(pm)

    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_mc_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_power_losses(pm, i)
        constraint_mc_model_voltage_magnitude_difference(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    # Objective
    _PM.objective_min_fuel_cost(pm)
end
