"Solve optimal switching problem"
function _solve_mc_osw(data::Union{Dict{String,<:Any}, String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, _build_mc_osw; kwargs...)
end


"Solve mixed-integer optimal switching problem"
function _solve_mc_osw_mi(data::Union{Dict{String,<:Any}, String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, _build_mc_osw_mi; kwargs...)
end


"Constructor for Optimal Switching"
function _build_mc_osw(pm::_PM.AbstractPowerModel)
    variable_mc_bus_voltage(pm)
    variable_mc_branch_power(pm)
    variable_mc_transformer_power(pm)

    variable_mc_switch_power(pm)
    variable_mc_switch_state(pm; relax=true)

    variable_mc_generator_power(pm)
    variable_mc_load_power(pm)
    variable_mc_storage_power(pm)

    constraint_mc_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    # generators should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for i in ids(pm, :bus)
        constraint_mc_power_balance(pm, i)
    end

    for i in ids(pm, :storage)
        _PM.constraint_storage_state(pm, i)
        _PM.constraint_storage_complementarity_nl(pm, i)
        constraint_mc_storage_losses(pm, i)
        constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state_on_off(pm, i; relax=true)
        constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    objective_mc_min_fuel_cost_switch(pm)
end


"constructor for branch flow osw"
function _build_mc_osw(pm::AbstractUBFModels)
    # Variables
    variable_mc_bus_voltage(pm)
    variable_mc_branch_current(pm)
    variable_mc_branch_power(pm)

    variable_mc_switch_power(pm)
    variable_mc_switch_state(pm; relax=true)

    variable_mc_transformer_power(pm)
    variable_mc_generator_power(pm)
    variable_mc_load_power(pm)
    variable_mc_storage_power(pm)

    # Constraints
    constraint_mc_model_current(pm)

    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for i in ids(pm, :bus)
        constraint_mc_power_balance(pm, i)
    end

    for i in ids(pm, :storage)
        _PM.constraint_storage_state(pm, i)
        _PM.constraint_storage_complementarity_nl(pm, i)
        constraint_mc_storage_losses(pm, i)
        constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_power_losses(pm, i)
        constraint_mc_model_voltage_magnitude_difference(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state_on_off(pm, i; relax=true)
        constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    # Objective
    objective_mc_min_fuel_cost_switch(pm)
end


"constructor for OSW in current-voltage variable space"
function _build_mc_osw(pm::_PM.AbstractIVRModel)
    # Variables
    variable_mc_bus_voltage(pm)
    variable_mc_branch_current(pm)
    variable_mc_switch_current(pm)
    variable_mc_switch_state(pm; relax=true)
    variable_mc_transformer_current(pm)
    variable_mc_generator_current(pm)
    variable_mc_load_current(pm)

    # Constraints
    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for i in ids(pm, :bus)
        constraint_mc_current_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_current_from(pm, i)
        constraint_mc_current_to(pm, i)

        constraint_mc_bus_voltage_drop(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state_on_off(pm, i; relax=true)
        constraint_mc_switch_current_limit(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    # Objective
    objective_mc_min_fuel_cost_switch(pm)
end


"Constructor for Optimal Switching"
function _build_mc_osw_mi(pm::_PM.AbstractPowerModel)
    variable_mc_bus_voltage(pm)
    variable_mc_branch_power(pm)
    variable_mc_transformer_power(pm)

    variable_mc_switch_power(pm)
    variable_mc_switch_state(pm; relax=false)

    variable_mc_generator_power(pm)
    variable_mc_load_power(pm)
    variable_mc_storage_power_mi(pm)

    constraint_mc_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    # generators should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for i in ids(pm, :bus)
        constraint_mc_power_balance(pm, i)
    end

    for i in ids(pm, :storage)
        _PM.constraint_storage_state(pm, i)
        _PM.constraint_storage_complementarity_mi(pm, i)
        constraint_mc_storage_losses(pm, i)
        constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state_on_off(pm, i; relax=false)
        constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    objective_mc_min_fuel_cost_switch(pm)
end


"constructor for mixed-integer branch flow osw"
function _build_mc_osw_mi(pm::AbstractUBFModels)
    # Variables
    variable_mc_bus_voltage(pm)
    variable_mc_branch_current(pm)
    variable_mc_branch_power(pm)

    variable_mc_switch_power(pm)
    variable_mc_switch_state(pm; relax=false)

    variable_mc_transformer_power(pm)
    variable_mc_generator_power(pm)
    variable_mc_load_power(pm)
    variable_mc_storage_power_mi(pm; relax=false)

    # Constraints
    constraint_mc_model_current(pm)

    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for i in ids(pm, :bus)
        constraint_mc_power_balance(pm, i)
    end

    for i in ids(pm, :storage)
        _PM.constraint_storage_state(pm, i)
        _PM.constraint_storage_complementarity_mi(pm, i)
        constraint_mc_storage_losses(pm, i)
        constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_power_losses(pm, i)
        constraint_mc_model_voltage_magnitude_difference(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state_on_off(pm, i; relax=false)
        constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    # Objective
    objective_mc_min_fuel_cost_switch(pm)
end

# Depreciated run_ functions (remove after ~4-6 months)

"depreciation warning for _run_mc_osw"
function _run_mc_osw(data::Union{Dict{String,<:Any}, String}, model_type::Type, solver; kwargs...)
    @warn "_run_mc_osw is being depreciated in favor of _solve_mc_osw, please update your code accordingly"
    return _solve_mc_osw(data, model_type, solver; kwargs...)
end


"depreciation warning for _run_mc_osw_mi"
function _run_mc_osw_mi(data::Union{Dict{String,<:Any}, String}, model_type::Type, solver; kwargs...)
    @warn "_run_mc_osw_mi is being depreciated in favor of _solve_mc_osw_mi, please update your code accordingly"
    return _solve_mc_osw_mi(data, model_type, solver; kwargs...)
end
