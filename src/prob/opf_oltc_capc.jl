"""
    solve_mc_opf_capc(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)

Solve OPF with capacitor control
"""
function solve_mc_opf_oltc_capc(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mc_opf_oltc_capc; kwargs...)
end


"""
    build_mc_opf_capc(pm::AbstractUnbalancedPowerModel)

Constructor for capcontrol OPF
"""
function build_mc_opf_oltc_capc(pm::AbstractUnbalancedPowerModel)
    variable_mc_bus_voltage(pm)

    variable_mc_branch_power(pm)
    variable_mc_switch_power(pm)
    variable_mc_transformer_power(pm)

    variable_mc_oltc_transformer_tap(pm)

    variable_mc_generator_power(pm)
    variable_mc_storage_power(pm)
    variable_mc_load_power(pm)

    variable_mc_capcontrol(pm; relax=true)

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
        constraint_mc_power_balance_capc(pm, i)
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
        constraint_mc_transformer_power(pm, i; fix_taps=false)
    end

    objective_mc_min_fuel_cost(pm)
end


"constructor for branch flow opf"
function build_mc_opf_oltc_capc(pm::AbstractUBFModels)
    # Variables
    variable_mc_bus_voltage(pm)

    variable_mc_branch_current(pm)
    variable_mc_branch_power(pm)
    variable_mc_switch_power(pm)
    variable_mc_transformer_power(pm)

    variable_mc_oltc_transformer_tap(pm)

    variable_mc_generator_power(pm)
    variable_mc_storage_power(pm)
    variable_mc_load_power(pm)

    variable_mc_capcontrol(pm; relax=true)

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
        constraint_mc_power_balance_capc(pm, i)
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

        constraint_mc_ampacity_from(pm, i)
        constraint_mc_ampacity_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
        constraint_mc_switch_thermal_limit(pm, i)
        constraint_mc_switch_ampacity(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i; fix_taps=false)
    end

    # Objective
    objective_mc_min_fuel_cost(pm)
end
