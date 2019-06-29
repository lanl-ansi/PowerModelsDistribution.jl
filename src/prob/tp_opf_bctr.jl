export run_tp_opf_bctr, run_ac_tp_opf_bctr


""
function run_tp_opf_bctr(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf_bctr; multiconductor=true, solution_builder=solution_bctr!, kwargs...)
end


""
function run_tp_opf_bctr(file::String, model_constructor, solver; kwargs...)
    data = ThreePhasePowerModels.parse_file(file)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf_bctr; multiconductor=true, solution_builder=solution_bctr!,  kwargs...)
end


""
function post_tp_opf_bctr(pm::PMs.GenericPowerModel)
    add_arcs_trans!(pm)

    variable_tp_voltage(pm)
    variable_tp_branch_flow(pm)

    for c in PMs.conductor_ids(pm)
        PMs.variable_generation(pm, cnd=c)
        PMs.variable_dcline_flow(pm, cnd=c)
    end
    variable_tp_trans_flow(pm)

    constraint_tp_voltage(pm)

    for i in PMs.ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    for i in PMs.ids(pm, :bus)
        bus = PMs.ref(pm, pm.cnw, :bus, i)
        constraint_tp_voltage_balance(pm, i)
        # KCL
        for c in PMs.conductor_ids(pm)
            PMs.constraint_kcl_shunt(pm, i, cnd=c)
        end
    end

    for i in PMs.ids(pm, :branch)
        for c in PMs.conductor_ids(pm)
            constraint_ohms_tp_yt_from(pm, i, cnd=c)
            constraint_ohms_tp_yt_to(pm, i, cnd=c)

            PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

            PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end

    for i in PMs.ids(pm, :dcline), c in PMs.conductor_ids(pm)
        PMs.constraint_dcline(pm, i, cnd=c)
    end

    for i in PMs.ids(pm, :trans)
        constraint_tp_trans(pm, i)
    end

    PMs.objective_min_fuel_cost(pm)
end

""
function solution_bctr!(pm::PMs.GenericPowerModel, sol::Dict{String,<:Any})
    PMs.add_bus_voltage_setpoint(sol, pm)
    add_setpoint_bus_voltage_balance_indicators!(pm, sol)
    PMs.add_generator_power_setpoint(sol, pm)
end
