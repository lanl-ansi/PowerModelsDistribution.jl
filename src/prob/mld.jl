"Run load shedding problem with storage"
function run_mc_mld(data::Dict{String,Any}, model_type, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mc_mld; kwargs...)
end


""
function run_mc_mld(file::String, model_type, solver; kwargs...)
    return run_mc_mld(PowerModelsDistribution.parse_file(file), model_type, solver; kwargs...)
end


"Run Branch Flow Model Load Shedding Problem"
function run_mc_mld_bf(data::Dict{String,Any}, model_type, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mc_mld_bf; kwargs...)
end


""
function run_mc_mld_bf(file::String, model_type, solver; kwargs...)
    return run_mc_mld_bf(PowerModelsDistribution.parse_file(file), model_type, solver; kwargs...)
end


"Run unit commitment load shedding problem (!relaxed)"
function run_mc_mld_uc(data::Dict{String,Any}, model_type, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mc_mld_uc; kwargs...)
end


""
function run_mc_mld_uc(file::String, model_type, solver; kwargs...)
    return run_mc_mld(PowerModelsDistribution.parse_file(file), model_type, solver; kwargs...)
end


"Load shedding problem including storage (snap-shot)"
function build_mc_mld(pm::_PM.AbstractPowerModel)
    variable_mc_bus_voltage_indicator(pm; relax=true)
    variable_mc_bus_voltage_on_off(pm)

    variable_mc_branch_power(pm)
    variable_mc_transformer_power(pm)

    variable_mc_gen_indicator(pm; relax=true)
    variable_mc_gen_power_setpoint_on_off(pm)

    variable_mc_storage_power_mi(pm; relax=true)

    variable_mc_load_indicator(pm; relax=true)
    variable_mc_shunt_indicator(pm; relax=true)

    constraint_mc_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    constraint_mc_bus_voltage_on_off(pm)

    for i in ids(pm, :gen)
        constraint_mc_gen_power_on_off(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_mc_shed_power_balance(pm, i)
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

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    objective_mc_min_load_setpoint_delta(pm)
end


"Load shedding problem for Branch Flow model"
function build_mc_mld_bf(pm::_PM.AbstractPowerModel)
    variable_mc_bus_voltage_indicator(pm; relax=true)
    variable_mc_bus_voltage_on_off(pm)

    variable_mc_branch_current(pm)
    variable_mc_branch_power(pm)
    variable_mc_transformer_power(pm)

    variable_mc_gen_indicator(pm; relax=true)
    variable_mc_gen_power_setpoint_on_off(pm)

    variable_mc_load_indicator(pm; relax=true)
    variable_mc_shunt_indicator(pm; relax=true)

    constraint_mc_model_current(pm)

    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    constraint_mc_bus_voltage_on_off(pm)

    for i in ids(pm, :gen)
        constraint_mc_gen_power_on_off(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_mc_shed_power_balance(pm, i)
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

    objective_mc_min_load_setpoint_delta(pm)
end


"Standard unit commitment (!relaxed) load shedding problem"
function build_mc_mld_uc(pm::_PM.AbstractPowerModel)
    variable_mc_bus_voltage_indicator(pm; relax=false)
    variable_mc_bus_voltage_on_off(pm)

    variable_mc_branch_power(pm)
    variable_mc_transformer_power(pm)

    variable_mc_gen_indicator(pm; relax=false)
    variable_mc_gen_power_setpoint_on_off(pm)

    variable_mc_storage_power(pm)
    variable_mc_storage_indicator(pm; relax=false)
    variable_mc_storage_power_on_off(pm)

    variable_mc_load_indicator(pm; relax=false)
    variable_mc_shunt_indicator(pm; relax=false)

    constraint_mc_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    constraint_mc_bus_voltage_on_off(pm)

    for i in ids(pm, :gen)
        constraint_mc_gen_power_on_off(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_mc_shed_power_balance(pm, i)
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

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    objective_mc_min_load_setpoint_delta(pm)
end
