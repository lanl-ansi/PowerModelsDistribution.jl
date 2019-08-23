# These problem formulations are used to debug Three Phase datasets
# that do not converge using the standard formulations
""
function run_mc_opf_pbs(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, post_mc_opf_pbs; multiconductor=true, solution_builder=solution_pbs!, kwargs...)
end


""
function run_mc_opf_pbs(file::String, model_constructor, solver; kwargs...)
    return run_mc_opf_pbs(PowerModelsDistribution.parse_file(file), model_constructor, solver; kwargs...)
end


""
function post_mc_opf_pbs(pm::_PMs.AbstractPowerModel)
    variable_mc_voltage(pm)
    variable_mc_branch_flow(pm)

    for c in _PMs.conductor_ids(pm)
        variable_mc_bus_power_slack(pm, cnd=c)
        _PMs.variable_generation(pm, cnd=c)
        _PMs.variable_dcline_flow(pm, cnd=c)
    end

    constraint_mc_model_voltage(pm)

    for i in _PMs.ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    for i in _PMs.ids(pm, :bus), c in _PMs.conductor_ids(pm)
        constraint_mc_power_balance_shunt_slack(pm, i, cnd=c)
    end

    for i in _PMs.ids(pm, :branch), c in _PMs.conductor_ids(pm)
        constraint_mc_ohms_yt_from(pm, i, cnd=c)
        constraint_mc_ohms_yt_to(pm, i, cnd=c)

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
function run_mc_pf_pbs(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, post_mc_pf_pbs; multiconductor=true, solution_builder=solution_pbs!, kwargs...)
end


""
function run_mc_pf_pbs(file::String, model_constructor, solver; kwargs...)
    return run_mc_pf_pbs(PowerModelsDistribution.parse_file(file), model_constructor, solver; kwargs...)
end


""
function post_mc_pf_pbs(pm::_PMs.AbstractPowerModel)
    variable_mc_voltage(pm, bounded=false)
    variable_mc_branch_flow(pm, bounded=false)

    for c in _PMs.conductor_ids(pm)
        variable_mc_bus_power_slack(pm, cnd=c)
        _PMs.variable_generation(pm, bounded=false, cnd=c)
        _PMs.variable_dcline_flow(pm, bounded=false, cnd=c)
    end

    constraint_mc_model_voltage(pm)

    for (i,bus) in _PMs.ref(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
        for c in _PMs.conductor_ids(pm)

            @assert bus["bus_type"] == 3
            _PMs.constraint_voltage_magnitude_setpoint(pm, i, cnd=c)
        end
    end

    for (i,bus) in _PMs.ref(pm, :bus), c in _PMs.conductor_ids(pm)
        constraint_mc_power_balance_shunt_slack(pm, i, cnd=c)

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
        constraint_mc_ohms_yt_from(pm, i, cnd=c)
        constraint_mc_ohms_yt_to(pm, i, cnd=c)
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
