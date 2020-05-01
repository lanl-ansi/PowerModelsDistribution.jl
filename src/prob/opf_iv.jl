"OPF in current-voltage variable space"
function run_mc_opf_iv(data::Union{Dict{String,<:Any},String}, model_type::DataType, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mc_opf_iv; kwargs...)
end


"constructor for OPF in current-voltage variable space"
function build_mc_opf_iv(pm::_PM.AbstractPowerModel)
    # Variables
    variable_mc_bus_voltage(pm)
    variable_mc_branch_current(pm)
    variable_mc_transformer_current(pm)
    variable_mc_gen_power_setpoint(pm)
    variable_mc_load_setpoint(pm)

    # Constraints
    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_gen_setpoint(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_setpoint(pm, id)
    end

    for i in ids(pm, :bus)
        constraint_mc_load_current_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_current_from(pm, i)
        constraint_mc_current_to(pm, i)

        constraint_mc_bus_voltage_drop(pm, i)

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
