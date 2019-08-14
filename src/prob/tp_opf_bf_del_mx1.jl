""
function run_tp_opf_bf_del_mx1(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, post_tp_opf_bf_del_mx1; solution_builder=solution_tp!, multiconductor=true, kwargs...)
end


""
function run_tp_opf_bf_del_mx1(file::String, model_constructor, solver; kwargs...)
    return run_tp_opf_bf_del_mx1(PowerModelsDistribution.parse_file(file), model_constructor, solver; kwargs...)
end


""
function post_tp_opf_bf_del_mx1(pm::_PMs.GenericPowerModel)
    # Variables
    variable_tp_voltage(pm)
    variable_tp_branch_current(pm)
    variable_tp_branch_flow(pm)

    variable_tp_generation_mx(pm)
    variable_tp_load_mx(pm)

    for c in _PMs.conductor_ids(pm)
        _PMs.variable_dcline_flow(pm, cnd=c)
    end

    # Constraints
    constraint_tp_model_current(pm)

    for i in _PMs.ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    for i in _PMs.ids(pm, :branch)
        constraint_tp_flow_losses(pm, i)

        constraint_tp_model_voltage_magnitude_difference(pm, i)

        for c in _PMs.conductor_ids(pm)
            #_PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

            _PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            _PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end

    for i in _PMs.ids(pm, :load)
        constraint_tp_load_mx(pm, i)
        constraint_tp_load_mx_SWL_only_delta(pm, i)
    end

    for i in _PMs.ids(pm, :gen)
        #constraint_tp_generation_mx_SWL(pm, i)
    end

    for i in _PMs.ids(pm, :dcline), c in _PMs.conductor_ids(pm)
        #_PMs.constraint_dcline(pm, i, cnd=c)
    end

    # This has to happen after constraining the load now!
    # Otherwise incorrect values on the diagonal, initialized at zero!
    for i in _PMs.ids(pm, :bus)
        constraint_tp_power_balance_mx_shunt(pm, i)
    end

    # Objective
    _PMs.objective_min_fuel_cost(pm)
end
