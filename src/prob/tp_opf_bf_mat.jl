export run_tp_opf_bf_mat

""
function run_tp_opf_bf_mat(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    if model_constructor != SDPUBFPowerModel && model_constructor != SOCNLPUBFPowerModel && model_constructor != SOCConicUBFPowerModel && model_constructor != LPUBFPowerModel && model_constructor != LPdiagUBFPowerModel
        error(LOGGER, "The problem type tp_opf_bf_mat at the moment only supports the SDPUBFPowerModel formulation")
    end
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf_bf_mat; solution_builder=get_solution_tp, multiconductor=true, kwargs...)
end


""
function run_tp_opf_bf_mat(file::String, model_constructor, solver; kwargs...)
    if model_constructor != SDPUBFPowerModel && model_constructor != SOCNLPUBFPowerModel && model_constructor != SOCConicUBFPowerModel && model_constructor != LPUBFPowerModel
        error(LOGGER, "The problem type tp_opf_bf_mat at the moment only supports the SDPUBFPowerModel formulation")
    end
    data = ThreePhasePowerModels.parse_file(file)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf_bf_mat; solution_builder=get_solution_tp, multiconductor=true, kwargs...)
end


""
function post_tp_opf_bf_mat(pm::GenericPowerModel)
    # Variables
    n_cond = 3
    variable_tp_voltage(pm, n_cond = n_cond)
    variable_branch_current(pm, n_cond = n_cond)
    variable_branch_flow(pm, n_cond = n_cond)

    for c in PMs.conductor_ids(pm)
        PMs.variable_generation(pm, cnd=c)

        # PMs.variable_dcline_flow(pm, cnd=c)
    end

    # Constraints
    for i in ids(pm, :ref_buses)
        constraint_tp_voltage_ref(pm, i)
    end

    for i in ids(pm, :bus)
        for c in PMs.conductor_ids(pm)
            constraint_tp_kcl_shunt(pm, i, cnd=c)
        end
        constraint_tp_valid_inquality_bus(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_tp_flow_losses_mat(pm, i)
        constraint_tp_branch_current_mat(pm, i)

        constraint_tp_voltage_magnitude_difference_mat(pm, i)
        constraint_tp_valid_inquality_branch(pm, i)

        for c in PMs.conductor_ids(pm)
            # PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

            PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end
    #
    # for i in ids(pm, :dcline)
    #     for h in PMs.phase_ids(pm)
    #         PMs.constraint_dcline(pm, i, ph=h)
    #     end
    # end

    #Objective
     # PMs.objective_min_fuel_cost(pm)

    #objective_min_losses(pm)
    objective_min_slack(pm)
end
