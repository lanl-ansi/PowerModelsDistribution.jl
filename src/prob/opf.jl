""
function run_ac_mc_opf(file, solver; kwargs...)
    return run_mc_opf(file, _PMs.ACPPowerModel, solver; kwargs...)
end


""
function run_mc_opf(data::Dict{String,Any}, model_type, solver; kwargs...)
    return _PMs.run_model(data, model_type, solver, post_mc_opf; multiconductor=true, ref_extensions=[ref_add_arcs_trans!], kwargs...)
end


""
function run_mc_opf(file::String, model_type, solver; kwargs...)
    return run_mc_opf(PowerModelsDistribution.parse_file(file), model_type, solver; kwargs...)
end


""
function post_mc_opf(pm::_PMs.AbstractPowerModel)
    variable_mc_voltage(pm)
    variable_mc_branch_flow(pm)
    variable_mc_transformer_flow(pm)

    for c in _PMs.conductor_ids(pm)
        _PMs.variable_generation(pm, cnd=c)
        _PMs.variable_dcline_flow(pm, cnd=c)
    end

    constraint_mc_model_voltage(pm)

    for i in _PMs.ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    for i in _PMs.ids(pm, :bus)
        constraint_mc_power_balance(pm, i)
    end

    for i in _PMs.ids(pm, :branch)
        constraint_mc_voltage_angle_difference(pm, i)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)

        for c in _PMs.conductor_ids(pm)
            _PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            _PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end

    for i in _PMs.ids(pm, :dcline), c in _PMs.conductor_ids(pm)
        _PMs.constraint_dcline(pm, i, cnd=c)
    end

    for i in _PMs.ids(pm, :transformer)
        constraint_mc_trans(pm, i)
    end

    _PMs.objective_min_fuel_cost(pm)
end
