# This problem includes load models beyond simple constant power ones; this is
# handled in variable_mc_load_mx and constraint_mc_load_mx.
# Also, KCL now links the full matrix variables, and not only the diagonals.
# This is what 'mx' refers to. Long-term, this might become a formulation and
# not a separate problem specification.


""
function run_mc_opf_bf_lm_mx(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, post_mc_opf_bf_lm_mx; solution_builder=solution_tp!, multiconductor=true, kwargs...)
end


""
function run_mc_opf_bf_lm_mx(file::String, model_constructor, solver; kwargs...)
    return run_mc_opf_bf_lm_mx(PowerModelsDistribution.parse_file(file), model_constructor, solver; kwargs...)
end


""
function post_mc_opf_bf_lm_mx(pm::_PMs.AbstractPowerModel)
    # Variables
    variable_mc_voltage(pm)
    variable_mc_branch_current(pm)
    variable_mc_branch_flow(pm)

    variable_mc_generation_mx(pm)
    variable_mc_load_mx(pm)

    # Constraints
    constraint_mc_model_current(pm)

    for i in _PMs.ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    for i in _PMs.ids(pm, :branch)
        constraint_mc_flow_losses(pm, i)
        constraint_mc_model_voltage_magnitude_difference(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    for i in _PMs.ids(pm, :load)
        constraint_mc_load_mx(pm, i)
    end

    for i in _PMs.ids(pm, :gen)
        constraint_mc_generation_mx(pm, i)
    end

    # This has to happen after constraining the load now!
    # Otherwise incorrect values on the diagonal, initialized at zero!
    for i in _PMs.ids(pm, :bus)
        constraint_mc_power_balance_mx_shunt(pm, i)
    end

    # Objective
    _PMs.objective_min_fuel_cost(pm)
end
