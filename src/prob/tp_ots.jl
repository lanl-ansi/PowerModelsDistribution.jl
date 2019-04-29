export run_tp_ots

""
function run_tp_ots(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_ots; multiconductor=true, solution_builder=PMs.get_ots_solution, kwargs...)
end


""
function run_tp_ots(file::String, model_constructor, solver; kwargs...)
    data = ThreePhasePowerModels.parse_file(file)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_ots; multiconductor=true, solution_builder=PMs.get_ots_solution, kwargs...)
end


""
function post_tp_ots(pm::PMs.GenericPowerModel)
    for c in PMs.conductor_ids(pm)
        PMs.variable_branch_indicator(pm, cnd=c)
        PMs.variable_voltage_on_off(pm, cnd=c)
        PMs.variable_generation(pm, cnd=c)
        PMs.variable_branch_flow(pm, cnd=c)
        PMs.variable_dcline_flow(pm, cnd=c)
    end

    for c in PMs.conductor_ids(pm)
        PMs.constraint_voltage_on_off(pm, cnd=c)

        for i in PMs.ids(pm, :ref_buses)
            constraint_tp_theta_ref(pm, i, cnd=c)
        end

        for i in PMs.ids(pm, :bus)
            PMs.constraint_kcl_shunt(pm, i, cnd=c)
        end

        for i in PMs.ids(pm, :branch)
            constraint_ohms_tp_yt_from_on_off(pm, i, cnd=c)
            constraint_ohms_tp_yt_to_on_off(pm, i, cnd=c)

            PMs.constraint_voltage_angle_difference_on_off(pm, i, cnd=c)

            PMs.constraint_thermal_limit_from_on_off(pm, i, cnd=c)
            PMs.constraint_thermal_limit_to_on_off(pm, i, cnd=c)
        end

        for i in PMs.ids(pm, :dcline)
            PMs.constraint_dcline(pm, i, cnd=c)
        end
    end

    PMs.objective_min_fuel_cost(pm)
end
