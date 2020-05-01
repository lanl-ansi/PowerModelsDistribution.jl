"Power Flow in current-voltage variable space"
function run_mc_pf_iv(data::Union{Dict{String,<:Any},String}, model_type::DataType, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mc_pf_iv; kwargs...)
end


"Constructor for Power Flow in current-voltage variable space"
function build_mc_pf_iv(pm::_PM.AbstractPowerModel)
    # Variables
    variable_mc_bus_voltage(pm, bounded = false)
    variable_mc_branch_current(pm, bounded = false)
    variable_mc_transformer_current(pm, bounded = false)
    variable_mc_gen_power_setpoint(pm, bounded = false)
    variable_mc_load_setpoint(pm, bounded = false)

    # Constraints
    for (i,bus) in ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3
        constraint_mc_theta_ref(pm, i)
        constraint_mc_voltage_magnitude_only(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_gen_setpoint(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_setpoint(pm, id)
    end

    for (i,bus) in ref(pm, :bus)
        constraint_mc_load_current_balance(pm, i)

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

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end
end
