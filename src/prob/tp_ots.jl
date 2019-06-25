""
function run_tp_ots(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, post_tp_ots; multiconductor=true, solution_builder=_PMs.get_ots_solution, kwargs...)
end


""
function run_tp_ots(file::String, model_constructor, solver; kwargs...)
    data = PowerModelsDistribution.parse_file(file)
    return _PMs.run_model(data, model_constructor, solver, post_tp_ots; multiconductor=true, solution_builder=_PMs.get_ots_solution, kwargs...)
end


""
function post_tp_ots(pm::_PMs.GenericPowerModel)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_branch_indicator(pm, cnd=c)
        _PMs.variable_voltage_on_off(pm, cnd=c)
        _PMs.variable_generation(pm, cnd=c)
        _PMs.variable_branch_flow(pm, cnd=c)
        _PMs.variable_dcline_flow(pm, cnd=c)
    end

    for c in _PMs.conductor_ids(pm)
        _PMs.constraint_voltage_on_off(pm, cnd=c)

        for i in _PMs.ids(pm, :ref_buses)
            constraint_tp_theta_ref(pm, i, cnd=c)
        end

        for i in _PMs.ids(pm, :bus)
            _PMs.constraint_power_balance_shunt(pm, i, cnd=c)
        end

        for i in _PMs.ids(pm, :branch)
            constraint_tp_ohms_yt_from_on_off(pm, i, cnd=c)
            constraint_tp_ohms_yt_to_on_off(pm, i, cnd=c)

            _PMs.constraint_voltage_angle_difference_on_off(pm, i, cnd=c)

            _PMs.constraint_thermal_limit_from_on_off(pm, i, cnd=c)
            _PMs.constraint_thermal_limit_to_on_off(pm, i, cnd=c)
        end

        for i in _PMs.ids(pm, :dcline)
            _PMs.constraint_dcline(pm, i, cnd=c)
        end
    end

    _PMs.objective_min_fuel_cost(pm)
end
