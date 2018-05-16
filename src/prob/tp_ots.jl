export run_tp_ots

""
function run_tp_ots(file, model_constructor, solver; kwargs...)
    return run_generic_tp_model(file, model_constructor, solver, post_tp_ots; solution_builder=get_ots_solution, kwargs...)
end


""
function post_tp_ots(pm::GenericPowerModel)
    for h in PMs.phase_ids(pm)
        PMs.variable_branch_indicator(pm, ph=h)
        PMs.variable_voltage_on_off(pm, ph=h)
        PMs.variable_generation(pm, ph=h)
        PMs.variable_branch_flow(pm, ph=h)
        PMs.variable_dcline_flow(pm, ph=h)

        PMs.constraint_voltage_on_off(pm, ph=h)

        for i in ids(pm, :ref_buses)
            PMs.constraint_theta_ref(pm, i, ph=h)
        end

        for i in ids(pm, :bus)
            PMs.constraint_kcl_shunt(pm, i, ph=h)
        end

        for i in ids(pm, :branch)
            constraint_ohms_yt_from_on_off(pm, i, ph=h)
            constraint_ohms_yt_to_on_off(pm, i, ph=h)

            PMs.constraint_voltage_angle_difference_on_off(pm, i, ph=h)

            PMs.constraint_thermal_limit_from_on_off(pm, i, ph=h)
            PMs.constraint_thermal_limit_to_on_off(pm, i, ph=h)
        end

        for i in ids(pm, :dcline)
            PMs.constraint_dcline(pm, i, ph=h)
        end
    end

    PMs.objective_min_fuel_cost(pm)
end


""
function get_ots_solution(pm::GenericPowerModel, sol::Dict{String,Any})
    PMs.add_bus_voltage_setpoint(sol, pm)
    PMs.add_generator_power_setpoint(sol, pm)
    PMs.add_branch_flow_setpoint(sol, pm)
    PMs.add_branch_status_setpoint(sol, pm)
end