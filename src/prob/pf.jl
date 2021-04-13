"Power Flow Problem"
function solve_mc_pf(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mc_pf; kwargs...)
end


"Constructor for Power Flow Problem"
function build_mc_pf(pm::AbstractUnbalancedPowerModel)
    variable_mc_bus_voltage(pm; bounded=false)
    variable_mc_branch_power(pm; bounded=false)
    variable_mc_switch_power(pm; bounded=false)
    variable_mc_transformer_power(pm; bounded=false)
    variable_mc_generator_power(pm; bounded=false)
    variable_mc_load_power(pm; bounded=false)
    variable_mc_storage_power(pm; bounded=false)

    constraint_mc_model_voltage(pm)

    for (i,bus) in ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3

        constraint_mc_theta_ref(pm, i)
        constraint_mc_voltage_magnitude_only(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for (i,bus) in ref(pm, :bus)
        constraint_mc_power_balance(pm, i)

        # PV Bus Constraints
        if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            constraint_mc_voltage_magnitude_only(pm, i)
            for j in ref(pm, :bus_gens, i)
                constraint_mc_gen_power_setpoint_real(pm, j)
            end
        end
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
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
        constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end
end


"Constructor for Power Flow in current-voltage variable space"
function build_mc_pf(pm::AbstractUnbalancedIVRModel)
    # Variables
    variable_mc_bus_voltage(pm, bounded = false)
    variable_mc_branch_current(pm, bounded = false)
    variable_mc_switch_current(pm, bounded=false)
    variable_mc_transformer_current(pm, bounded = false)
    variable_mc_generator_current(pm, bounded = false)
    variable_mc_load_current(pm, bounded = false)

    # Constraints
    for (i,bus) in ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3
        constraint_mc_theta_ref(pm, i)
        constraint_mc_voltage_magnitude_only(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for (i,bus) in ref(pm, :bus)
        constraint_mc_current_balance(pm, i)

        # PV Bus Constraints
        if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2
            constraint_mc_voltage_magnitude_only(pm, i)
            for j in ref(pm, :bus_gens, i)
                constraint_mc_gen_power_setpoint_real(pm, j)
            end
        end
    end

    for i in ids(pm, :branch)
        constraint_mc_current_from(pm, i)
        constraint_mc_current_to(pm, i)

        constraint_mc_bus_voltage_drop(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
        constraint_mc_switch_current_limit(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end
end


"Constructor for Branch Flow Power Flow"
function build_mc_pf(pm::AbstractUBFModels)
    # Variables
    variable_mc_bus_voltage(pm; bounded=true) # TODO should be false
    variable_mc_branch_current(pm)
    variable_mc_branch_power(pm)
    variable_mc_switch_power(pm)
    variable_mc_transformer_power(pm; bounded=false)
    variable_mc_generator_power(pm; bounded=false)
    variable_mc_load_power(pm)
    variable_mc_storage_power(pm; bounded=false)

    # Constraints
    constraint_mc_model_current(pm)

    for (i,bus) in ref(pm, :ref_buses)
        if !(typeof(pm)<:LPUBFDiagPowerModel)
            constraint_mc_theta_ref(pm, i)
        end

        @assert bus["bus_type"] == 3
        constraint_mc_voltage_magnitude_only(pm, i)
    end

    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
    end

    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for (i,bus) in ref(pm, :bus)
        constraint_mc_power_balance(pm, i)

        # PV Bus Constraints
        if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            constraint_mc_voltage_magnitude_only(pm, i)
            for j in ref(pm, :bus_gens, i)
                constraint_mc_gen_power_setpoint_real(pm, j)
            end
        end
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
end

# Deprecated run_ functions (remove in ~4-6 months)

"Power Flow problem with ACPUPowerModel"
function run_ac_mc_pf(data::Union{Dict{String,<:Any},String}, solver; kwargs...)
    @warn "run_ac_mc_pf is being depreciated in favor of solve_mc_pf(data, ACPUPowerModel, solver; kwargs...), please update your code accordingly"
    return solve_mc_pf(data, ACPUPowerModel, solver; kwargs...)
end


"Power Flow Problem"
function run_mc_pf(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    @warn "run_mc_pf is being depreciated in favor of solve_mc_pf, please update your code accordingly"
    return solve_mc_pf(data, model_type, solver; kwargs...)
end
