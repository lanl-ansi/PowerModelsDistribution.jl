export run_tp_opf_bf_mx

""
function run_tp_opf_bf_mx(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf_bf_mx; solution_builder=get_solution_tp, multiconductor=true, kwargs...)
end


""
function run_tp_opf_bf_mx(file::String, model_constructor, solver; kwargs...)
    data = ThreePhasePowerModels.parse_file(file)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf_bf_mx; solution_builder=get_solution_tp, multiconductor=true, kwargs...)
end


""
function post_tp_opf_bf_mx(pm::PMs.GenericPowerModel)
    shunts_diag2mat!(pm)

    # Variables
    variable_tp_voltage_v2(pm)
    variable_tp_branch_current_v2(pm)
    variable_tp_branch_flow_v2(pm)

    variable_tp_generation_power_mx(pm)
    variable_tp_generation_current_mx(pm)
    variable_tp_load_power_mx(pm)
    variable_tp_load_current_mx(pm)

    for c in PMs.conductor_ids(pm)
        PMs.variable_dcline_flow(pm, cnd=c)
    end

    # Constraints
    for i in PMs.ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    for i in PMs.ids(pm, :bus), c in PMs.conductor_ids(pm)
        PMs.constraint_kcl_shunt(pm, i, cnd=c)
    end
    constraint_tp_voltage_psd(pm)

    for i in PMs.ids(pm, :branch)
        constraint_tp_flow_losses_v2(pm, i)

        constraint_tp_voltage_magnitude_difference_v2(pm, i)

        constraint_tp_branch_current_v2(pm, i)


        for c in PMs.conductor_ids(pm)
            PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

            PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end

    for i in PMs.conductor_ids(pm, :load)
        constraint_tp_load_current(pm, i)
    end

    for i in PMs.conductor_ids(pm, :gen)
        constraint_tp_generation_current(pm, i)
    end

    for i in PMs.ids(pm, :dcline), c in PMs.conductor_ids(pm)
        PMs.constraint_dcline(pm, i, cnd=c)
    end

    PMs.objective_min_fuel_cost(pm)
end
