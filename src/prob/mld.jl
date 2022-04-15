"Solve load shedding problem with storage"
function solve_mc_mld(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mc_mld; kwargs...)
end


"Solve multinetwork load shedding problem with storage"
function solve_mn_mc_mld_simple(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mn_mc_mld_simple; multinetwork=true, kwargs...)
end


"Solve unit commitment load shedding problem (!relaxed)"
function solve_mc_mld_uc(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mc_mld_uc; kwargs...)
end


"Load shedding problem including storage (snap-shot)"
function build_mc_mld(pm::AbstractUnbalancedPowerModel)
    variable_mc_bus_voltage_indicator(pm; relax=true)
    variable_mc_bus_voltage_on_off(pm)

    variable_mc_branch_power(pm)
    variable_mc_switch_power(pm)
    variable_mc_transformer_power(pm)

    variable_mc_gen_indicator(pm; relax=true)
    variable_mc_generator_power_on_off(pm)

    variable_mc_storage_indicator(pm, relax=true)
    variable_mc_storage_power_mi_on_off(pm, relax=true)

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
        constraint_mc_power_balance_shed(pm, i)
    end

    for i in ids(pm, :storage)
        constraint_storage_state(pm, i)
        constraint_storage_complementarity_mi(pm, i)
        constraint_mc_storage_on_off(pm, i)
        constraint_mc_storage_losses(pm, i)
        constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
        constraint_mc_ampacity_from(pm, i)
        constraint_mc_ampacity_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
        constraint_mc_switch_thermal_limit(pm, i)
        constraint_mc_switch_ampacity(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    objective_mc_min_load_setpoint_delta(pm)
end


""
function build_mc_mld(pm::AbstractUnbalancedIVRModel)
    error("IVRUPowerModel is not yet supported in the MLD problem space")
    # TODO
end


"Multinetwork load shedding problem including storage"
function build_mn_mc_mld_simple(pm::AbstractUnbalancedPowerModel)
    for (n, network) in nws(pm)
        variable_mc_branch_power(pm; nw=n)
        variable_mc_switch_power(pm; nw=n)
        variable_mc_transformer_power(pm; nw=n)
        variable_mc_generator_power(pm; nw=n)
        variable_mc_bus_voltage(pm; nw=n)

        variable_mc_load_indicator(pm; nw=n, relax=true)
        variable_mc_shunt_indicator(pm; nw=n, relax=true)
        variable_mc_storage_power_mi(pm; nw=n, relax=true)
        variable_mc_storage_power_mi_on_off_ne(pm; nw=n)

        constraint_mc_model_voltage(pm; nw=n)

        for i in ids(pm, n, :ref_buses)
            constraint_mc_theta_ref(pm, i; nw=n)
        end

        for i in ids(pm, n, :gen)
            constraint_mc_generator_power(pm, i; nw=n)
        end

        for i in ids(pm, n, :bus)
            constraint_mc_power_balance_shed(pm, i; nw=n)
        end

        for i in ids(pm, n, :storage)
            constraint_mc_storage_losses(pm, i; nw=n)
            constraint_mc_storage_thermal_limit(pm, i; nw=n)
            constraint_storage_complementarity_mi(pm, i; nw=n)
        end

        for i in ids(pm, n, :storage_ne)
            constraint_mc_storage_losses_ne(pm, i; nw=n)
            constraint_mc_storage_thermal_limit_ne(pm, i; nw=n)
            constraint_storage_complementarity_mi_ne(pm, i; nw=n)
            constraint_storage_indicator_expand_ne(pm, i; nw=n)
        end

        for i in ids(pm, n, :branch)
            constraint_mc_ohms_yt_from(pm, i; nw=n)
            constraint_mc_ohms_yt_to(pm, i; nw=n)
            constraint_mc_voltage_angle_difference(pm, i; nw=n)
            constraint_mc_thermal_limit_from(pm, i; nw=n)
            constraint_mc_thermal_limit_to(pm, i; nw=n)
            constraint_mc_ampacity_from(pm, i; nw=n)
            constraint_mc_ampacity_to(pm, i; nw=n)
        end

        for i in ids(pm, n, :switch)
            constraint_mc_switch_state(pm, i; nw=n)
            constraint_mc_switch_thermal_limit(pm, i; nw=n)
            constraint_mc_switch_ampacity(pm, i; nw=n)
        end

        for i in ids(pm, n, :transformer)
            constraint_mc_transformer_power(pm, i; nw=n)
        end
    end

    network_ids = sort(collect(nw_ids(pm)))

    n_1 = network_ids[1]

    for i in ids(pm, :storage; nw=n_1)
        constraint_storage_state(pm, i; nw=n_1)
    end

    for i in ids(pm, :storage_ne; nw=n_1)
        constraint_storage_state_ne(pm, i; nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in ids(pm, :storage; nw=n_2)
            constraint_storage_state(pm, i, n_1, n_2)
        end

        for i in ids(pm, :storage_ne; nw=n_2)
            constraint_storage_state_ne(pm, i, n_1, n_2)
        end

        n_1 = n_2
    end
    
    objective_mc_min_load_setpoint_delta_simple(pm)
end


"Load shedding problem for Branch Flow model"
function build_mc_mld(pm::AbstractUBFModels)
    variable_mc_bus_voltage_indicator(pm; relax=true)
    variable_mc_bus_voltage_on_off(pm)

    variable_mc_branch_current(pm)
    variable_mc_branch_power(pm)
    variable_mc_switch_power(pm)
    variable_mc_transformer_power(pm)

    variable_mc_gen_indicator(pm; relax=true)
    variable_mc_generator_power_on_off(pm)

    variable_mc_storage_indicator(pm, relax=true)
    variable_mc_storage_power_mi_on_off(pm, relax=true)

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
        constraint_mc_power_balance_shed(pm, i)
    end

    for i in ids(pm, :storage)
        constraint_storage_state(pm, i)
        constraint_storage_complementarity_mi(pm, i)
        constraint_mc_storage_on_off(pm, i)
        constraint_mc_storage_losses(pm, i)
        constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_power_losses(pm, i)
        constraint_mc_model_voltage_magnitude_difference(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
        constraint_mc_ampacity_from(pm, i)
        constraint_mc_ampacity_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
        constraint_mc_switch_thermal_limit(pm, i)
        constraint_mc_switch_ampacity(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    objective_mc_min_load_setpoint_delta(pm)
end


"Multinetwork load shedding problem for Branch Flow model"
function build_mn_mc_mld_simple(pm::AbstractUBFModels)
    for (n, network) in nws(pm)
        variable_mc_branch_current(pm; nw=n)
        variable_mc_branch_power(pm; nw=n)
        variable_mc_switch_power(pm; nw=n)
        variable_mc_transformer_power(pm; nw=n)
        variable_mc_generator_power(pm; nw=n)
        variable_mc_bus_voltage(pm; nw=n)

        variable_mc_load_indicator(pm; nw=n, relax=true)
        variable_mc_shunt_indicator(pm; nw=n, relax=true)
        variable_mc_storage_power_mi(pm; nw=n, relax=true)
        variable_mc_storage_power_mi_on_off_ne(pm; nw=n)

        constraint_mc_model_current(pm; nw=n)

        for i in ids(pm, n, :ref_buses)
            constraint_mc_theta_ref(pm, i; nw=n)
        end

        for i in ids(pm, n, :gen)
            constraint_mc_generator_power(pm, i; nw=n)
        end

        for i in ids(pm, n, :bus)
            constraint_mc_power_balance_shed(pm, i; nw=n)
        end

        for i in ids(pm, n, :storage)
            constraint_mc_storage_losses(pm, i; nw=n)
            constraint_mc_storage_thermal_limit(pm, i; nw=n)
            constraint_storage_complementarity_mi(pm, i; nw=n)
        end

        for i in ids(pm, n, :storage_ne)
            constraint_mc_storage_losses_ne(pm, i; nw=n)
            constraint_mc_storage_thermal_limit_ne(pm, i; nw=n)
            constraint_storage_complementarity_mi_ne(pm, i; nw=n)
            constraint_storage_indicator_expand_ne(pm, i; nw=n)
        end

        for i in ids(pm, n, :branch)
            constraint_mc_power_losses(pm, i; nw=n)
            constraint_mc_model_voltage_magnitude_difference(pm, i; nw=n)
            constraint_mc_voltage_angle_difference(pm, i; nw=n)
            constraint_mc_thermal_limit_from(pm, i; nw=n)
            constraint_mc_thermal_limit_to(pm, i; nw=n)
            constraint_mc_ampacity_from(pm, i; nw=n)
            constraint_mc_ampacity_to(pm, i; nw=n)
        end

        for i in ids(pm, n, :switch)
            constraint_mc_switch_state(pm, i; nw=n)
            constraint_mc_switch_thermal_limit(pm, i; nw=n)
            constraint_mc_switch_ampacity(pm, i; nw=n)
        end

        for i in ids(pm, n, :transformer)
            constraint_mc_transformer_power(pm, i; nw=n)
        end
    end

    network_ids = sort(collect(nw_ids(pm)))

    n_1 = network_ids[1]

    for i in ids(pm, :storage; nw=n_1)
        constraint_storage_state(pm, i; nw=n_1)
    end

    for i in ids(pm, :storage_ne; nw=n_1)
        constraint_storage_state_ne(pm, i; nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in ids(pm, :storage; nw=n_2)
            constraint_storage_state(pm, i, n_1, n_2)
        end

        for i in ids(pm, :storage_ne; nw=n_2)
            constraint_storage_state_ne(pm, i, n_1, n_2)
        end

        n_1 = n_2
    end

    objective_mc_min_load_setpoint_delta_simple(pm)
end


"Load shedding problem for Branch Flow model"
function build_mc_mld_bf(pm::AbstractUnbalancedPowerModel)
    build_mc_mld(pm)

    variable_mc_bus_voltage_indicator(pm; relax=true)
    variable_mc_bus_voltage_on_off(pm)

    variable_mc_branch_current(pm)
    variable_mc_branch_power(pm)
    variable_mc_switch_power(pm)
    variable_mc_transformer_power(pm)

    variable_mc_gen_indicator(pm; relax=true)
    variable_mc_generator_power_on_off(pm)

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
        constraint_mc_power_balance_shed(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_power_losses(pm, i)
        constraint_mc_model_voltage_magnitude_difference(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
        constraint_mc_ampacity_from(pm, i)
        constraint_mc_ampacity_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
        constraint_mc_switch_thermal_limit(pm, i)
        constraint_mc_switch_ampacity(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    objective_mc_min_load_setpoint_delta(pm)
end


"Standard unit commitment (!relaxed) load shedding problem"
function build_mc_mld_uc(pm::AbstractUnbalancedPowerModel)
    variable_mc_bus_voltage_indicator(pm; relax=false)
    variable_mc_bus_voltage_on_off(pm)

    variable_mc_branch_power(pm)
    variable_mc_switch_power(pm)
    variable_mc_transformer_power(pm)

    variable_mc_gen_indicator(pm; relax=false)
    variable_mc_generator_power_on_off(pm)

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
        constraint_mc_power_balance_shed(pm, i)
    end

    for i in ids(pm, :storage)
        constraint_storage_state(pm, i)
        constraint_storage_complementarity_nl(pm, i)
        constraint_mc_storage_losses(pm, i)
        constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
        constraint_mc_ampacity_from(pm, i)
        constraint_mc_ampacity_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
        constraint_mc_switch_thermal_limit(pm, i)
        constraint_mc_switch_ampacity(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    objective_mc_min_load_setpoint_delta(pm)
end
