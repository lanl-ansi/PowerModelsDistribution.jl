# These problem formulations are used to debug Three Phase datasets
# that do not converge using the standard formulations
""
function run_tp_opf_pbs(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, post_tp_opf_pbs; multiconductor=true, solution_builder = get_pbs_solution, kwargs...)
end


""
function run_tp_opf_pbs(file::String, model_constructor, solver; kwargs...)
    data = PowerModelsDistribution.parse_file(file)
    return _PMs.run_model(data, model_constructor, solver, post_tp_opf_pbs; multiconductor=true, solution_builder = get_pbs_solution, kwargs...)
end


""
function post_tp_opf_pbs(pm::_PMs.GenericPowerModel)
    variable_tp_voltage(pm)
    variable_tp_branch_flow(pm)

    for c in _PMs.conductor_ids(pm)
        variable_tp_bus_power_slack(pm, cnd=c)
        _PMs.variable_generation(pm, cnd=c)
        _PMs.variable_dcline_flow(pm, cnd=c)
    end

    constraint_tp_model_voltage(pm)

    for i in _PMs.ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    for i in _PMs.ids(pm, :bus), c in _PMs.conductor_ids(pm)
        constraint_tp_power_balance_shunt_slack(pm, i, cnd=c)
    end

    for i in _PMs.ids(pm, :branch), c in _PMs.conductor_ids(pm)
        constraint_tp_ohms_yt_from(pm, i, cnd=c)
        constraint_tp_ohms_yt_to(pm, i, cnd=c)

        _PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

        _PMs.constraint_thermal_limit_from(pm, i, cnd=c)
        _PMs.constraint_thermal_limit_to(pm, i, cnd=c)
    end

    for i in _PMs.ids(pm, :dcline), c in _PMs.conductor_ids(pm)
        _PMs.constraint_dcline(pm, i, cnd=c)
    end

    objective_min_bus_power_slack(pm)
end



""
function run_tp_pf_pbs(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, post_tp_pf_pbs; multiconductor=true, solution_builder = get_pbs_solution, kwargs...)
end


""
function run_tp_pf_pbs(file::String, model_constructor, solver; kwargs...)
    data = PowerModelsDistribution.parse_file(file)
    return _PMs.run_model(data, model_constructor, solver, post_tp_pf_pbs; multiconductor=true, solution_builder = get_pbs_solution, kwargs...)
end

""
function post_tp_pf_pbs(pm::_PMs.GenericPowerModel)
    variable_tp_voltage(pm, bounded=false)
    variable_tp_branch_flow(pm, bounded=false)

    for c in _PMs.conductor_ids(pm)
        variable_tp_bus_power_slack(pm, cnd=c)
        _PMs.variable_generation(pm, bounded=false, cnd=c)
        _PMs.variable_dcline_flow(pm, bounded=false, cnd=c)
    end

    constraint_tp_model_voltage(pm)

    for (i,bus) in _PMs.ref(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
        for c in _PMs.conductor_ids(pm)

            @assert bus["bus_type"] == 3
            _PMs.constraint_voltage_magnitude_setpoint(pm, i, cnd=c)
        end
    end

    for (i,bus) in _PMs.ref(pm, :bus), c in _PMs.conductor_ids(pm)
        constraint_tp_power_balance_shunt_slack(pm, i, cnd=c)

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

    for i in _PMs.ids(pm, :branch), c in _PMs.conductor_ids(pm)
        constraint_tp_ohms_yt_from(pm, i, cnd=c)
        constraint_tp_ohms_yt_to(pm, i, cnd=c)
    end

    for (i,dcline) in _PMs.ref(pm, :dcline), c in _PMs.conductor_ids(pm)
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

    objective_min_bus_power_slack(pm)
end


""
function get_pbs_solution(pm::_PMs.GenericPowerModel, sol::Dict{String,Any})
    _PMs.add_setpoint_bus_voltage!(sol, pm)
    _PMs.add_setpoint_generator_power!(sol, pm)
    _PMs.add_setpoint_branch_flow!(sol, pm)
    add_bus_slack_setpoint(sol, pm)
end


""
function add_bus_slack_setpoint(sol, pm::_PMs.GenericPowerModel)
    _PMs.add_setpoint!(sol, pm, "bus", "p_slack", :p_slack, status_name="bus_type", inactive_status_value=4)
    _PMs.add_setpoint!(sol, pm, "bus", "q_slack", :q_slack, status_name="bus_type", inactive_status_value=4)
end
