# This problem includes load models beyond simple constant power ones; this is
# handled in variable_mc_load and constraint_mc_load.


""
function run_mc_opf_bf_lm(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, build_mc_opf_bf_lm;  multiconductor=true, kwargs...)
end


""
function run_mc_opf_bf_lm(file::String, model_constructor, solver; kwargs...)
    return run_mc_opf_bf_lm(PowerModelsDistribution.parse_file(file), model_constructor, solver; kwargs...)
end


""
function build_mc_opf_bf_lm(pm::_PMs.AbstractPowerModel)
    # Variables
    variable_mc_voltage(pm)
    variable_mc_branch_current(pm)
    variable_mc_branch_flow(pm)

    variable_mc_generation(pm)
    variable_mc_load(pm)

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
        constraint_mc_load(pm, i)
    end

    for i in _PMs.ids(pm, :gen)
        constraint_mc_generation(pm, i)
    end

    for i in _PMs.ids(pm, :bus)
        constraint_mc_power_balance(pm, i)
    end

    # Objective
    _PMs.objective_min_fuel_cost(pm)
end
