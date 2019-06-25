""
function run_tp_pf_bf(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    if model_constructor != SDPUBFPowerModel && model_constructor != SOCNLPUBFPowerModel && model_constructor != SOCConicUBFPowerModel && model_constructor != LPUBFPowerModel && model_constructor != LPdiagUBFPowerModel && model_constructor !=  SOCBFPowerModel
        Memento.error(_LOGGER, "The problem type tp_opf_bf at the moment only supports a limited set of formulations")
    end
    return _PMs.run_model(data, model_constructor, solver, post_tp_pf_bf; solution_builder=get_solution_tp, multiconductor=true, kwargs...)
end


""
function run_tp_pf_bf(file::String, model_constructor, solver; kwargs...)
    if model_constructor != SDPUBFPowerModel && model_constructor != SOCNLPUBFPowerModel && model_constructor != SOCConicUBFPowerModel && model_constructor != LPUBFPowerModel && model_constructor !=  SOCBFPowerModel
        Memento.error(_LOGGER, "The problem type tp_opf_bf at the moment only supports a limited set of formulations")
    end
    data = PowerModelsDistribution.parse_file(file)
    return _PMs.run_model(data, model_constructor, solver, post_tp_pf_bf; solution_builder=get_solution_tp, multiconductor=true, kwargs...)
end


""
function post_tp_pf_bf(pm::_PMs.GenericPowerModel)
    # Variables
    variable_tp_voltage(pm, bounded=false)
    variable_tp_branch_current(pm)
    variable_tp_branch_flow(pm)

    for c in _PMs.conductor_ids(pm)
        _PMs.variable_generation(pm, bounded=false, cnd=c)
        _PMs.variable_dcline_flow(pm, bounded=false, cnd=c)
    end

    # Constraints
    constraint_tp_model_current(pm)

    for (i,bus) in _PMs.ref(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)

        for c in _PMs.conductor_ids(pm)
            @assert bus["bus_type"] == 3
            # _PMs.constraint_voltage_magnitude_setpoint(pm, i, cnd=c) #TODO add back
        end
    end

    for i in _PMs.ids(pm, :bus), c in _PMs.conductor_ids(pm)
        _PMs.constraint_power_balance_shunt(pm, i, cnd=c)

        # PV Bus Constraints
        if length(_PMs.ref(pm, :bus_gens, i)) > 0 && !(i in _PMs.ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            _PMs.constraint_voltage_magnitude_setpoint(pm, i, cnd=c)
            for j in _PMs.ref(pm, :bus_gens, i)
                _PMs.constraint_active_gen_setpoint(pm, j, cnd=c)
            end
        end
    end

    for i in _PMs.ids(pm, :branch)
        constraint_tp_flow_losses(pm, i)

        constraint_tp_model_voltage_magnitude_difference(pm, i)

        for c in _PMs.conductor_ids(pm)
            constraint_voltage_angle_difference(pm, i, cnd=c)

            _PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            _PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end

    for (i,dcline) in _PMs.ref(pm, :dcline)
        for c in _PMs.conductor_ids(pm)

            _PMs.constraint_active_dcline_setpoint(pm, i, cnd=c)

            f_bus = _PMs.ref(pm, :bus)[dcline["f_bus"]]
            if f_bus["bus_type"] == 1
                _PMs.constraint_voltage_magnitude_setpoint(pm, f_bus["index"], cnd=c)
            end

            t_bus = _PMs.ref(pm, :bus)[dcline["t_bus"]]
            if t_bus["bus_type"] == 1
                _PMs.constraint_voltage_magnitude_setpoint(pm, t_bus["index"], cnd=c)
            end
        end
    end

     _PMs.objective_min_fuel_cost(pm)
end
