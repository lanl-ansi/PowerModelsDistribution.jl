"Run load shedding problem with storage"
function run_mc_mld(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mc_mld; kwargs...)
end


"Run multinetwork load shedding problem with storage"
function run_mn_mc_mld(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mn_mc_mld; kwargs...)
end


"Run Branch Flow Model Load Shedding Problem"
function run_mc_mld_bf(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    Memento.info(_LOGGER, "We recommend using run_mc_mld, which will attempt to select appropriate variables and constraints based on the specified formulation, instead of run_mc_mld_bf")
    return run_mc_model(data, model_type, solver, build_mc_mld_bf; kwargs...)
end


"Run unit commitment load shedding problem (!relaxed)"
function run_mc_mld_uc(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mc_mld_uc; kwargs...)
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


"Multinetwork load shedding problem including storage (snap-shot)"
function build_mn_mc_mld(pm::_PM.AbstractPowerModel)
    for (n, network) in _PM.nws(pm)
        variable_mc_bus_voltage_indicator(pm; nw=n, relax=true)
        variable_mc_bus_voltage_on_off(pm; nw=n)

        variable_mc_branch_power(pm; nw=n)
        variable_mc_transformer_power(pm; nw=n)

        variable_mc_gen_indicator(pm; nw=n, relax=true)
        variable_mc_gen_power_setpoint_on_off(pm; nw=n)

        variable_mc_storage_power_mi(pm; nw=n, relax=true)

        variable_mc_load_indicator(pm; nw=n, relax=true)
        variable_mc_shunt_indicator(pm; nw=n, relax=true)

        constraint_mc_model_voltage(pm; nw=n)

        for i in ids(pm, n, :ref_buses)
            constraint_mc_theta_ref(pm, i; nw=n)
        end

        constraint_mc_bus_voltage_on_off(pm; nw=n)

        for i in ids(pm, n, :gen)
            constraint_mc_gen_power_on_off(pm, i; nw=n)
        end

        for i in ids(pm, n, :bus)
            constraint_mc_shed_power_balance(pm, i; nw=n)
        end

        for i in ids(pm, n, :storage)
            _PM.constraint_storage_complementarity_nl(pm, i; nw=n)
            constraint_mc_storage_losses(pm, i; nw=n)
            constraint_mc_storage_thermal_limit(pm, i; nw=n)
        end

        for i in ids(pm, n, :branch)
            constraint_mc_ohms_yt_from(pm, i; nw=n)
            constraint_mc_ohms_yt_to(pm, i; nw=n)

            constraint_mc_voltage_angle_difference(pm, i; nw=n)

            constraint_mc_thermal_limit_from(pm, i; nw=n)
            constraint_mc_thermal_limit_to(pm, i; nw=n)
        end

        for i in ids(pm, n, :transformer)
            constraint_mc_transformer_power(pm, i; nw=n)
        end
    end

    network_ids = sort(collect(_PM.nw_ids(pm)))

    n_1 = network_ids[1]

    for i in _PM.ids(pm, :storage; nw=n_1)
        _PM.constraint_storage_state(pm, i; nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in _PM.ids(pm, :storage; nw=n_2)
            _PM.constraint_storage_state(pm, i, n_1, n_2)
        end

        n_1 = n_2
    end

    objective_mc_min_load_setpoint_delta(pm)
end


"Load shedding problem for Branch Flow model"
function build_mc_mld(pm::AbstractUBFModels)
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


"Multinetwork load shedding problem for Branch Flow model"
function build_mn_mc_mld(pm::AbstractUBFModels)
    for (n, network) in _PM.nws(pm)
        variable_mc_bus_voltage_indicator(pm; nw=n, relax=true)
        variable_mc_bus_voltage_on_off(pm; nw=n)

        variable_mc_branch_current(pm; nw=n)
        variable_mc_branch_power(pm; nw=n)
        variable_mc_transformer_power(pm; nw=n)

        variable_mc_gen_indicator(pm; nw=n, relax=true)
        variable_mc_gen_power_setpoint_on_off(pm; nw=n)

        variable_mc_load_indicator(pm; nw=n, relax=true)
        variable_mc_shunt_indicator(pm; nw=n, relax=true)

        constraint_mc_model_current(pm; nw=n)

        for i in ids(pm, n, :ref_buses)
            constraint_mc_theta_ref(pm, i; nw=n)
        end

        constraint_mc_bus_voltage_on_off(pm; nw=n)

        for i in ids(pm, n, :gen)
            constraint_mc_gen_power_on_off(pm, i; nw=n)
        end

        for i in ids(pm, n, :bus)
            constraint_mc_shed_power_balance(pm, i; nw=n)
        end

        for i in ids(pm, n, :branch)
            constraint_mc_power_losses(pm, i; nw=n)
            constraint_mc_model_voltage_magnitude_difference(pm, i; nw=n)

            constraint_mc_voltage_angle_difference(pm, i; nw=n)

            constraint_mc_thermal_limit_from(pm, i; nw=n)
            constraint_mc_thermal_limit_to(pm, i; nw=n)
        end

        for i in ids(pm, n, :transformer)
            constraint_mc_transformer_power(pm, i; nw=n)
        end
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
