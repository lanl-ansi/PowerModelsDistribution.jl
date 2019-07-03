""
function run_tp_mld(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, post_tp_mld; multiconductor=true, ref_extensions=[ref_add_arcs_trans!], solution_builder=solution_mld!, kwargs...)
end


""
function run_tp_mld(file::String, model_constructor, solver; kwargs...)
    return run_tp_mld(PowerModelsDistribution.parse_file(file), model_constructor, solver; kwargs...)
end


""
function run_tp_mld_strg(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, post_tp_mld_strg; multiconductor=true, ref_extensions=[ref_add_arcs_trans!], solution_builder=solution_mld!, kwargs...)
end


""
function run_tp_mld_strg(file::String, model_constructor, solver; kwargs...)
    return run_tp_mld_strg(PowerModelsDistribution.parse_file(file), model_constructor, solver; kwargs...)
end


""
function run_tp_mld_bf(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, post_tp_mld_bf; multiconductor=true, solution_builder=solution_mld_bf!, kwargs...)
end


""
function run_tp_mld_bf(file::String, model_constructor, solver; kwargs...)
    return run_tp_mld_bf(PowerModelsDistribution.parse_file(file), model_constructor, solver; kwargs...)
end


""
function post_tp_mld(pm::_PMs.GenericPowerModel)
    variable_tp_indicator_bus_voltage(pm; relax=true)
    variable_tp_bus_voltage_on_off(pm)

    variable_tp_trans_flow(pm)
    variable_tp_branch_flow(pm)

    _PMs.variable_generation_indicator(pm; relax=true)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_generation_on_off(pm, cnd=c)
        _PMs.variable_dcline_flow(pm, cnd=c)
    end

    variable_tp_demand_factor(pm)
    variable_tp_shunt_factor(pm)

    constraint_tp_model_voltage(pm)

    for i in _PMs.ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    for i in _PMs.ids(pm, :bus), c in _PMs.conductor_ids(pm)
        constraint_tp_power_balance_shunt_trans_shed(pm, i, cnd=c)
    end

    for i in _PMs.ids(pm, :branch)
        for c in _PMs.conductor_ids(pm)
            constraint_tp_ohms_yt_from(pm, i, cnd=c)
            constraint_tp_ohms_yt_to(pm, i, cnd=c)

            _PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

            _PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            _PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end

    for i in _PMs.ids(pm, :dcline), c in _PMs.conductor_ids(pm)
        _PMs.constraint_dcline(pm, i, cnd=c)
    end

    for i in _PMs.ids(pm, :trans)
        constraint_tp_trans(pm, i)
    end

    objective_tp_min_load_delta(pm)
end


""
function post_tp_mld_strg(pm::_PMs.GenericPowerModel)
    variable_tp_indicator_bus_voltage(pm; relax=true)
    variable_tp_bus_voltage_on_off(pm)

    variable_tp_trans_flow(pm)
    variable_tp_branch_flow(pm)

    variable_tp_storage(pm)

    _PMs.variable_generation_indicator(pm; relax=true)
    variable_tp_indicator_storage(pm; relax=true)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_generation_on_off(pm, cnd=c)
        variable_tp_storage_on_off(pm, cnd=c)
        _PMs.variable_dcline_flow(pm, cnd=c)
    end

    variable_tp_demand_factor(pm)
    variable_tp_shunt_factor(pm)

    constraint_tp_model_voltage(pm)

    for i in _PMs.ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    for i in _PMs.ids(pm, :bus), c in _PMs.conductor_ids(pm)
        constraint_tp_power_balance_shunt_trans_shed(pm, i, cnd=c)
    end

    for i in _PMs.ids(pm, :storage)
        _PMs.constraint_storage_state(pm, i)
        constraint_tp_storage_exchange(pm, i)
        for c in _PMs.conductor_ids(pm)
            _PMs.constraint_storage_thermal_limit(pm, i, cnd=c)
        end
    end

    for i in _PMs.ids(pm, :branch)
        for c in _PMs.conductor_ids(pm)
            constraint_tp_ohms_yt_from(pm, i, cnd=c)
            constraint_tp_ohms_yt_to(pm, i, cnd=c)

            _PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

            _PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            _PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end

    for i in _PMs.ids(pm, :dcline), c in _PMs.conductor_ids(pm)
        _PMs.constraint_dcline(pm, i, cnd=c)
    end

    for i in _PMs.ids(pm, :trans)
        constraint_tp_trans(pm, i)
    end

    objective_tp_min_load_delta_strg(pm)
end


""
function post_tp_mld_bf(pm::_PMs.GenericPowerModel)
    variable_tp_indicator_bus_voltage(pm; relax=true)
    variable_tp_bus_voltage_on_off(pm)

    variable_tp_branch_flow(pm)
    variable_tp_branch_current(pm)

    _PMs.variable_generation_indicator(pm; relax=true)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_generation_on_off(pm, cnd=c)
        _PMs.variable_dcline_flow(pm, cnd=c)
    end

    variable_tp_demand_factor(pm)
    variable_tp_shunt_factor(pm)

    constraint_tp_model_current(pm)

    for i in _PMs.ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    for i in _PMs.ids(pm, :bus), c in _PMs.conductor_ids(pm)
        constraint_tp_power_balance_shunt_shed(pm, i, cnd=c)
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

    objective_tp_max_loadability(pm)
end