export run_tp_opf_bf_v2

""
function run_tp_opf_bf_v2(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, post_tp_opf_bf_v2; solution_builder=solution_tp!, multiconductor=true, kwargs...)
end


""
function run_tp_opf_bf_v2(file::String, model_constructor, solver; kwargs...)
    return run_tp_opf_bf_v2(PowerModelsDistribution.parse_file(file), model_constructor, solver; kwargs...)
end


""
function post_tp_opf_bf_v2(pm::_PMs.GenericPowerModel)
    shunts_diag2mat!(pm)

    # Variables
    variable_tp_voltage_v2(pm)
    variable_tp_branch_current_v2(pm)
    variable_tp_branch_flow_v2(pm)

    for c in _PMs.conductor_ids(pm)
        _PMs.variable_generation(pm, cnd=c)
        _PMs.variable_dcline_flow(pm, cnd=c)
    end

    # Constraints
    constraint_tp_model_current(pm)

    for i in _PMs.ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    for i in _PMs.ids(pm, :bus), c in _PMs.conductor_ids(pm)
        _PMs.constraint_power_balance_shunt(pm, i, cnd=c)
    end

    for i in _PMs.ids(pm, :branch)
        constraint_tp_flow_losses(pm, i)

        constraint_tp_model_voltage_magnitude_difference(pm, i)

        for c in _PMs.conductor_ids(pm)
            _PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

            _PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            _PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end

    for i in _PMs.ids(pm, :dcline), c in _PMs.conductor_ids(pm)
        _PMs.constraint_dcline(pm, i, cnd=c)
    end

    # Objective
    _PMs.objective_min_fuel_cost(pm)
end
