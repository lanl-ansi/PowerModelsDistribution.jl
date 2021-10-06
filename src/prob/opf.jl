"""
	function solve_mc_opf(
		data::Union{Dict{String,<:Any},String},
		model_type::Type,
		solver;
		kwargs...
	)

Solve Optimal Power Flow
"""
function solve_mc_opf(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mc_opf; kwargs...)
end


"""
	function solve_mn_mc_opf(
		data::Union{Dict{String,<:Any},String},
		model_type::Type,
		solver;
		kwargs...
	)

Solve multinetwork optimal power flow problem
"""
function solve_mn_mc_opf(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mn_mc_opf; multinetwork=true, kwargs...)
end


"""
	function build_mc_opf(
		pm::AbstractUnbalancedPowerModel
	)

Constructor for Optimal Power Flow
"""
function build_mc_opf(pm::AbstractUnbalancedPowerModel)
    variable_mc_bus_voltage(pm)
    variable_mc_branch_power(pm)
    variable_mc_transformer_power(pm)
    variable_mc_switch_power(pm)
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
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
        constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    objective_mc_min_fuel_cost(pm)
end


"""
	function build_mc_opf(
		pm::AbstractUnbalancedIVRModel
	)

constructor for OPF in current-voltage variable space
"""
function build_mc_opf(pm::AbstractUnbalancedIVRModel)
    # Variables
    variable_mc_bus_voltage(pm)
    variable_mc_branch_current(pm)
    variable_mc_switch_current(pm)
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
        constraint_mc_switch_state(pm, i)
        constraint_mc_switch_current_limit(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    # Objective
    objective_mc_min_fuel_cost(pm)
end


"""
	function build_mc_opf(
		pm::AbstractUBFModels
	)

constructor for branch flow opf
"""
function build_mc_opf(pm::AbstractUBFModels)
    # Variables
    variable_mc_bus_voltage(pm)
    variable_mc_branch_current(pm)
    variable_mc_branch_power(pm)
    variable_mc_switch_power(pm)
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
        constraint_storage_state(pm, i)
        constraint_storage_complementarity_nl(pm, i)
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
        constraint_mc_switch_state(pm, i)
        constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    # Objective
    objective_mc_min_fuel_cost(pm)
end


"""
	function build_mc_opf(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for OPF in current-voltage variable space with explicit neutrals
"""
function build_mc_opf(pm::AbstractExplicitNeutralIVRModel)
    # Variables
    variable_mc_bus_voltage(pm)
    variable_mc_branch_current(pm)
    variable_mc_load_current(pm)
    variable_mc_load_power(pm)
    variable_mc_generator_current(pm)
    variable_mc_generator_power(pm)
    variable_mc_transformer_current(pm)
    variable_mc_transformer_power(pm)
    variable_mc_switch_current(pm)

    # Constraints
    for i in ids(pm, :bus)

        if i in ids(pm, :ref_buses)
            constraint_mc_voltage_reference(pm, i)
        end

        constraint_mc_voltage_absolute(pm, i)
        constraint_mc_voltage_pairwise(pm, i)
    end

    # components should be constrained before KCL, or the bus current variables might be undefined

    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
        constraint_mc_generator_current(pm, id)
    end

    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
        constraint_mc_load_current(pm, id)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_voltage(pm, i)
        constraint_mc_transformer_current(pm, i)

        constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_current_from(pm, i)
        constraint_mc_current_to(pm, i)
        constraint_mc_bus_voltage_drop(pm, i)

        constraint_mc_branch_current_limit(pm, i)
        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_current(pm, i)
        constraint_mc_switch_state(pm, i)

        constraint_mc_switch_current_limit(pm, i)
        constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_mc_current_balance(pm, i)
    end

    # Objective
    objective_mc_min_fuel_cost(pm)
end


"""
	function build_mc_opf(
		pm::AbstractExplicitNeutralACRModel
	)

constructor for OPF in power-voltage variable space with explicit neutrals
"""
function build_mc_opf(pm::AbstractExplicitNeutralACRModel)
    # Variables
    variable_mc_bus_voltage(pm)
    variable_mc_branch_power(pm)
    variable_mc_load_power(pm)
    variable_mc_generator_power(pm)
    variable_mc_transformer_power(pm)
    variable_mc_switch_power(pm)

    # Constraints
    for i in ids(pm, :bus)

        if i in ids(pm, :ref_buses)
            constraint_mc_voltage_reference(pm, i)
        end

        constraint_mc_voltage_absolute(pm, i)
        constraint_mc_voltage_pairwise(pm, i)
    end

    # components should be constrained before KCL, or the bus current variables might be undefined

    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
        # constraint_mc_generator_current(pm, id)
    end

    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_voltage(pm, i)
        constraint_mc_transformer_power(pm, i)

        constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)

        constraint_mc_branch_current_limit(pm, i)
        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_power(pm, i)
        constraint_mc_switch_state(pm, i)

        constraint_mc_switch_current_limit(pm, i)
        constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_mc_power_balance(pm, i)
    end

    # Objective
    objective_mc_min_fuel_cost(pm)
end


"""
	function build_mn_mc_opf(
		pm::AbstractUnbalancedPowerModel
	)

Multinetwork optimal power flow problem
"""
function build_mn_mc_opf(pm::AbstractUnbalancedPowerModel)
    for (n, network) in nws(pm)
        variable_mc_bus_voltage(pm; nw=n)
        variable_mc_branch_power(pm; nw=n)
        variable_mc_switch_power(pm; nw=n)
        variable_mc_transformer_power(pm; nw=n)
        variable_mc_generator_power(pm; nw=n)
        variable_mc_load_power(pm; nw=n)
        variable_mc_storage_power(pm; nw=n)

        constraint_mc_model_voltage(pm; nw=n)

        for i in ids(pm, n, :ref_buses)
            constraint_mc_theta_ref(pm, i; nw=n)
        end

        for id in ids(pm, n, :gen)
            constraint_mc_generator_power(pm, id; nw=n)
        end

        for id in ids(pm, n, :load)
            constraint_mc_load_power(pm, id; nw=n)
        end

        for i in ids(pm, n, :bus)
            constraint_mc_power_balance(pm, i; nw=n)
        end

        for i in ids(pm, n, :storage)
            constraint_storage_complementarity_nl(pm, i; nw=n)
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

        for i in ids(pm, n, :switch)
            constraint_mc_switch_state(pm, i; nw=n)
            constraint_mc_switch_thermal_limit(pm, i; nw=n)
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

    for n_2 in network_ids[2:end]
        for i in ids(pm, :storage; nw=n_2)
            constraint_storage_state(pm, i, n_1, n_2)
        end

        n_1 = n_2
    end

    objective_mc_min_fuel_cost(pm)
end

"""
	function build_mn_mc_opf(
		pm::AbstractUnbalancedIVRModel
	)

Multinetwork current-voltage optimal power flow problem
"""
function build_mn_mc_opf(pm::AbstractUnbalancedIVRModel)
    for (n, network) in nws(pm)
        variable_mc_bus_voltage(pm; nw=n)
        variable_mc_branch_current(pm; nw=n)
        variable_mc_switch_current(pm; nw=n)
        variable_mc_transformer_current(pm; nw=n)
        variable_mc_generator_current(pm; nw=n)
        variable_mc_load_current(pm; nw=n)
        variable_mc_storage_power(pm; nw=n)

        for i in ids(pm, n, :ref_buses)
            constraint_mc_theta_ref(pm, i; nw=n)
        end

        for id in ids(pm, n, :gen)
            constraint_mc_generator_power(pm, id; nw=n)
        end

        for id in ids(pm, n, :load)
            constraint_mc_load_power(pm, id; nw=n)
        end

        for i in ids(pm, n, :bus)
            constraint_mc_current_balance(pm, i; nw=n)
        end

        for i in ids(pm, n, :storage)
            constraint_storage_complementarity_nl(pm, i; nw=n)
            constraint_mc_storage_losses(pm, i; nw=n)
            constraint_mc_storage_thermal_limit(pm, i; nw=n)
        end

        for i in ids(pm, n, :branch)
            constraint_mc_current_from(pm, i; nw=n)
            constraint_mc_current_to(pm, i; nw=n)
            constraint_mc_bus_voltage_drop(pm, i; nw=n)
            constraint_mc_voltage_angle_difference(pm, i; nw=n)
            constraint_mc_thermal_limit_from(pm, i; nw=n)
            constraint_mc_thermal_limit_to(pm, i; nw=n)
        end

        for i in ids(pm, n, :switch)
            constraint_mc_switch_state(pm, i; nw=n)
            constraint_mc_switch_current_limit(pm, i; nw=n)
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

    for n_2 in network_ids[2:end]
        for i in ids(pm, :storage; nw=n_2)
            constraint_storage_state(pm, i, n_1, n_2)
        end

        n_1 = n_2
    end

    objective_mc_min_fuel_cost(pm)
end

"""
	function build_mn_mc_opf(
		pm::AbstractUBFModels
	)

Multinetwork branch flow optimal power flow problem
"""
function build_mn_mc_opf(pm::AbstractUBFModels)
    for (n, network) in nws(pm)
        variable_mc_bus_voltage(pm; nw=n)
        variable_mc_branch_current(pm; nw=n)
        variable_mc_branch_power(pm; nw=n)
        variable_mc_switch_power(pm; nw=n)
        variable_mc_transformer_power(pm; nw=n)
        variable_mc_generator_power(pm; nw=n)
        variable_mc_load_power(pm; nw=n)
        variable_mc_storage_power(pm; nw=n)

        constraint_mc_model_current(pm; nw=n)

        for i in ids(pm, n, :ref_buses)
            constraint_mc_theta_ref(pm, i; nw=n)
        end

        for id in ids(pm, n, :gen)
            constraint_mc_generator_power(pm, id; nw=n)
        end

        for id in ids(pm, n, :load)
            constraint_mc_load_power(pm, id; nw=n)
        end

        for i in ids(pm, n, :bus)
            constraint_mc_power_balance(pm, i; nw=n)
        end

        for i in ids(pm, n, :storage)
            constraint_storage_complementarity_nl(pm, i; nw=n)
            constraint_mc_storage_losses(pm, i; nw=n)
            constraint_mc_storage_thermal_limit(pm, i; nw=n)
        end

        for i in ids(pm, n, :branch)
            constraint_mc_power_losses(pm, i; nw=n)
            constraint_mc_model_voltage_magnitude_difference(pm, i; nw=n)
            constraint_mc_voltage_angle_difference(pm, i; nw=n)
            constraint_mc_thermal_limit_from(pm, i; nw=n)
            constraint_mc_thermal_limit_to(pm, i; nw=n)
        end

        for i in ids(pm, n, :switch)
            constraint_mc_switch_state(pm, i; nw=n)
            constraint_mc_switch_thermal_limit(pm, i; nw=n)
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

    for n_2 in network_ids[2:end]
        for i in ids(pm, :storage; nw=n_2)
            constraint_storage_state(pm, i, n_1, n_2)
        end

        n_1 = n_2
    end

    objective_mc_min_fuel_cost(pm)
end
