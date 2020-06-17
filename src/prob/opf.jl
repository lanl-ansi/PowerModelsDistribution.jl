"OPF with ACPPowerModel"
function run_ac_mc_opf(data::Union{Dict{String,<:Any},String}, solver; kwargs...)
    return run_mc_opf(data, ACPPowerModel, solver; kwargs...)
end


"Optimal Power Flow"
function run_mc_opf(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mc_opf; kwargs...)
end


"Run multinetwork optimal power flow problem"
function run_mn_mc_opf(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mn_mc_opf; kwargs...)
end


"Constructor for Optimal Power Flow"
function build_mc_opf(pm::_PM.AbstractPowerModel)
    variable_mc_bus_voltage(pm)
    variable_mc_branch_power(pm)
    variable_mc_transformer_power(pm)
    variable_mc_gen_power_setpoint(pm)
    variable_mc_load_setpoint(pm)
    variable_mc_storage_power(pm)

    constraint_mc_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    # generators should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_gen_setpoint(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_setpoint(pm, id)
    end

    for i in ids(pm, :bus)
        constraint_mc_load_power_balance(pm, i)
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

    _PM.objective_min_fuel_cost(pm)
end


"constructor for OPF in current-voltage variable space"
function build_mc_opf(pm::_PM.AbstractIVRModel)
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


"constructor for branch flow opf"
function build_mc_opf(pm::AbstractUBFModels)
    # Variables
    variable_mc_bus_voltage(pm)
    variable_mc_branch_current(pm)
    variable_mc_branch_power(pm)
    variable_mc_transformer_power(pm)
    variable_mc_gen_power_setpoint(pm)
    variable_mc_load_setpoint(pm)

    # Constraints
    constraint_mc_model_current(pm)

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
        constraint_mc_load_power_balance(pm, i)
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


"Multinetwork optimal power flow problem"
function build_mn_mc_opf(pm::_PM.AbstractPowerModel)
    for (n, network) in _PM.nws(pm)
        variable_mc_bus_voltage(pm; nw=n)
        variable_mc_branch_power(pm; nw=n)
        variable_mc_transformer_power(pm; nw=n)
        variable_mc_gen_power_setpoint(pm; nw=n)
        variable_mc_load_setpoint(pm; nw=n)
        variable_mc_storage_power(pm; nw=n)

        constraint_mc_model_voltage(pm; nw=n)

        for i in ids(pm, n, :ref_buses)
            constraint_mc_theta_ref(pm, i; nw=n)
        end

        for id in ids(pm, n, :gen)
            constraint_mc_gen_setpoint(pm, id; nw=n)
        end

        for id in ids(pm, n, :load)
            constraint_mc_load_setpoint(pm, id; nw=n)
        end

        for i in ids(pm, n, :bus)
            constraint_mc_load_power_balance(pm, i; nw=n)
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

    _PM.objective_min_fuel_cost(pm)
end

"Multinetwork current-voltage optimal power flow problem"
function build_mn_mc_opf(pm::_PM.AbstractIVRModel)
    for (n, network) in _PM.nws(pm)
        variable_mc_bus_voltage(pm; nw=n)
        variable_mc_branch_current(pm; nw=n)
        variable_mc_transformer_current(pm; nw=n)
        variable_mc_gen_power_setpoint(pm; nw=n)
        variable_mc_load_setpoint(pm; nw=n)

        for i in ids(pm, n, :ref_buses)
            constraint_mc_theta_ref(pm, i; nw=n)
        end

        for id in ids(pm, n, :gen)
            constraint_mc_gen_setpoint(pm, id; nw=n)
        end

        for id in ids(pm, n, :load)
            constraint_mc_load_setpoint(pm, id; nw=n)
        end

        for i in ids(pm, n, :bus)
            constraint_mc_load_current_balance(pm, i; nw=n)
        end

        for i in ids(pm, n, :branch)
            constraint_mc_current_from(pm, i; nw=n)
            constraint_mc_current_to(pm, i; nw=n)
            constraint_mc_bus_voltage_drop(pm, i; nw=n)
            constraint_mc_voltage_angle_difference(pm, i; nw=n)
            constraint_mc_thermal_limit_from(pm, i; nw=n)
            constraint_mc_thermal_limit_to(pm, i; nw=n)
        end

        for i in ids(pm, n, :transformer)
            constraint_mc_transformer_power(pm, i; nw=n)
        end
    end

    _PM.objective_min_fuel_cost(pm)
end

"Multinetwork branch flow optimal power flow problem"
function build_mn_mc_opf(pm::AbstractUBFModels)
    for (n, network) in _PM.nws(pm)
        variable_mc_bus_voltage(pm; nw=n)
        variable_mc_branch_current(pm; nw=n)
        variable_mc_branch_power(pm; nw=n)
        variable_mc_transformer_power(pm; nw=n)
        variable_mc_gen_power_setpoint(pm; nw=n)
        variable_mc_load_setpoint(pm; nw=n)

        constraint_mc_model_current(pm; nw=n)

        for i in ids(pm, n, :ref_buses)
            constraint_mc_theta_ref(pm, i; nw=n)
        end

        for id in ids(pm, n, :gen)
            constraint_mc_gen_setpoint(pm, id; nw=n)
        end

        for id in ids(pm, n, :load)
            constraint_mc_load_setpoint(pm, id; nw=n)
        end

        for i in ids(pm, n, :bus)
            constraint_mc_load_power_balance(pm, i; nw=n)
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

    _PM.objective_min_fuel_cost(pm)
end
