export run_tp_opf, run_ac_tp_opf

""
function run_ac_tp_opf(file, solver; kwargs...)
    return run_tp_opf(file, _PMs.ACPPowerModel, solver; multiconductor=true, kwargs...)
end


""
function run_tp_opf(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, post_tp_opf; multiconductor=true, kwargs...)
end


""
function run_tp_opf(file::String, model_constructor, solver; kwargs...)
    data = PowerModelsDistribution.parse_file(file)
    return _PMs.run_model(data, model_constructor, solver, post_tp_opf; multiconductor=true, kwargs...)
end


""
function post_tp_opf(pm::_PMs.GenericPowerModel)
    add_arcs_trans!(pm)

    variable_tp_voltage(pm)
    variable_tp_branch_flow(pm)

    for c in _PMs.conductor_ids(pm)
        _PMs.variable_generation(pm, cnd=c)
        _PMs.variable_dcline_flow(pm, cnd=c)
    end
    variable_tp_trans_flow(pm)

    constraint_tp_model_voltage(pm)

    for i in _PMs.ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    for i in _PMs.ids(pm, :bus), c in _PMs.conductor_ids(pm)
        constraint_tp_power_balance_shunt_trans(pm, i, cnd=c)
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

    _PMs.objective_min_fuel_cost(pm)
end
