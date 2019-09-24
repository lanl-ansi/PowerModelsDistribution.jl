""
function run_mc_pf_bf(data::Dict{String,Any}, model_type, solver; kwargs...)
    if model_type != SDPUBFPowerModel && model_type != SOCNLPUBFPowerModel && model_type != SOCConicUBFPowerModel && model_type != LPLinUBFPowerModel && model_type != LPdiagUBFPowerModel && model_type !=  SOCBFPowerModel
        Memento.error(_LOGGER, "The problem type mc_opf_bf at the moment only supports a limited set of formulations")
    end
    return _PMs.run_model(data, model_type, solver, post_mc_pf_bf; solution_builder=solution_bf!, ref_extensions=[ref_add_arcs_trans!], multiconductor=true, kwargs...)
end


""
function run_mc_pf_bf(file::String, model_type, solver; kwargs...)
    return run_mc_pf_bf(PowerModelsDistribution.parse_file(file), model_type, solver; kwargs...)
end


""
function post_mc_pf_bf(pm::_PMs.AbstractPowerModel)
    # Variables
    variable_mc_voltage(pm; bounded=false)
    variable_mc_branch_current(pm)
    variable_mc_branch_flow(pm)
    variable_mc_generation(pm; bounded=false)

    # Constraints
    constraint_mc_model_current(pm)

    for (i,bus) in _PMs.ref(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)

        @assert bus["bus_type"] == 3
        constraint_mc_voltage_magnitude_setpoint(pm, i)
    end

    for i in _PMs.ids(pm, :bus)
        constraint_mc_power_balance(pm, i)

        # PV Bus Constraints
        if length(_PMs.ref(pm, :bus_gens, i)) > 0 && !(i in _PMs.ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            constraint_mc_voltage_magnitude_setpoint(pm, i)
            for j in _PMs.ref(pm, :bus_gens, i)
                constraint_mc_active_gen_setpoint(pm, j)
            end
        end
    end

    for i in _PMs.ids(pm, :branch)
        constraint_mc_flow_losses(pm, i)
        constraint_mc_model_voltage_magnitude_difference(pm, i)
        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    # Objective
    _PMs.objective_min_fuel_cost(pm)
end
