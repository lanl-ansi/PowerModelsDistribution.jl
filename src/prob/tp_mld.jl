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
function post_tp_mld(pm::_PMs.GenericPowerModel)
    variable_tp_bus_voltage_indicator(pm; relax=true)
    variable_tp_bus_voltage_on_off(pm)

    variable_tp_trans_flow(pm)
    variable_tp_branch_flow(pm)

    variable_tp_generation_indicator(pm; relax=true)
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
    variable_tp_bus_voltage_indicator(pm; relax=true)
    variable_tp_bus_voltage_on_off(pm)

    variable_tp_trans_flow(pm)
    variable_tp_branch_flow(pm)

    variable_tp_storage(pm)

    variable_tp_generation_indicator(pm; relax=true)
    variable_tp_storage_indicator(pm; relax=true)
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
function solution_mld!(pm::_PMs.GenericPowerModel{T}, sol::Dict{String,Any}) where T
    _PMs.add_setpoint_bus_voltage!(sol, pm)
    _PMs.add_setpoint_generator_power!(sol, pm)
    _PMs.add_setpoint_storage!(sol, pm)
    _PMs.add_setpoint_branch_flow!(sol, pm)
    add_setpoint_bus_status!(sol, pm)
    add_setpoint_load!(sol, pm)
    add_setpoint_shunt!(sol, pm)
    add_setpoint_generator_status!(sol, pm)
    add_setpoint_storage_status!(sol, pm)
end


""
function add_setpoint_load!(sol, pm::_PMs.GenericPowerModel{T}) where T
    _PMs.add_setpoint!(sol, pm, "load", "pd", :z_demand; scale = (x,item,cnd) -> x*item["pd"][cnd])
    _PMs.add_setpoint!(sol, pm, "load", "qd", :z_demand; scale = (x,item,cnd) -> x*item["qd"][cnd])
    _PMs.add_setpoint!(sol, pm, "load", "status", :z_demand; default_value = (item) -> if (item["status"] == 0) 0.0 else 1.0 end, conductorless=true)
end


""
function add_setpoint_shunt!(sol, pm::_PMs.GenericPowerModel{T}) where T
    _PMs.add_setpoint!(sol, pm, "shunt", "gs", :z_shunt; scale = (x,item,cnd) -> x*item["gs"][cnd])
    _PMs.add_setpoint!(sol, pm, "shunt", "bs", :z_shunt; scale = (x,item,cnd) -> x*item["bs"][cnd])
    _PMs.add_setpoint!(sol, pm, "shunt", "status", :z_shunt; default_value = (item) -> if (item["status"] == 0) 0.0 else 1.0 end, conductorless=true)
end


""
function add_setpoint_bus_status!(sol, pm::_PMs.GenericPowerModel{T}) where T
   _PMs.add_setpoint!(sol, pm, "bus", "status", :z_voltage; status_name="bus_type", inactive_status_value=4, default_value = (item) -> if item["bus_type"] == 4 0.0 else 1.0 end, conductorless=true)
end


""
function add_setpoint_generator_status!(sol, pm::_PMs.GenericPowerModel{T}) where T
   _PMs.add_setpoint!(sol, pm, "gen", "gen_status", :z_gen; status_name="gen_status", default_value = (item) -> item["gen_status"]*1.0, conductorless=true)
end


""
function add_setpoint_storage_status!(sol, pm::_PMs.GenericPowerModel{T}) where T
    _PMs.add_setpoint!(sol, pm, "storage", "status", :z_storage; default_value = (item) -> item["status"]*1.0, conductorless=true)
end
