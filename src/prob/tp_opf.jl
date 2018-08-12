export run_tp_opf, run_ac_tp_opf

""
function run_ac_tp_opf(file, solver; kwargs...)
    return run_tp_opf(file, PMs.ACPPowerModel, solver; multiconductor=true, kwargs...)
end


""
function run_tp_opf(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf; multiconductor=true, kwargs...)
end


""
function run_tp_opf(file::String, model_constructor, solver; kwargs...)
    data = ThreePhasePowerModels.parse_file(file)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf; multiconductor=true, kwargs...)
end


""
function post_tp_opf(pm::GenericPowerModel)
    variable_tp_voltage(pm)
    variable_tp_branch_flow(pm)

    for c in PMs.conductor_ids(pm)
        PMs.variable_generation(pm, cnd=c)
        PMs.variable_dcline_flow(pm, cnd=c)
    end

    constraint_tp_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    for i in ids(pm, :bus), c in PMs.conductor_ids(pm)
        PMs.constraint_kcl_shunt(pm, i, cnd=c)
    end

    for i in ids(pm, :branch)
        for c in PMs.conductor_ids(pm)
            constraint_ohms_tp_yt_from(pm, i, cnd=c)
            constraint_ohms_tp_yt_to(pm, i, cnd=c)

            PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

            PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end

    for i in ids(pm, :dcline), c in PMs.conductor_ids(pm)
        PMs.constraint_dcline(pm, i, cnd=c)
    end

    PMs.objective_min_fuel_cost(pm)
end
