export run_tp_opf_bf

""
function run_tp_opf_bf(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    if model_constructor != SDPUBFPowerModel && model_constructor != SOCNLPUBFPowerModel && model_constructor != SOCConicUBFPowerModel && model_constructor != LPUBFPowerModel && model_constructor != LPdiagUBFPowerModel && model_constructor !=  SOCBFPowerModel
        error(LOGGER, "The problem type tp_opf_bf_mat at the moment only supports a limited set of formulations")
    end
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf_bf; solution_builder=get_solution_tp, multiconductor=true, kwargs...)
end


""
function run_tp_opf_bf(file::String, model_constructor, solver; kwargs...)
    if model_constructor != SDPUBFPowerModel && model_constructor != SOCNLPUBFPowerModel && model_constructor != SOCConicUBFPowerModel && model_constructor != LPUBFPowerModel && model_constructor !=  SOCBFPowerModel
        error(LOGGER, "The problem type tp_opf_bf_mat at the moment only supports a limited set of formulations")
    end
    data = ThreePhasePowerModels.parse_file(file)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf_bf; solution_builder=get_solution_tp, multiconductor=true, kwargs...)
end


""
function post_tp_opf_bf(pm::GenericPowerModel)
    # Variables
    variable_tp_voltage(pm)
    variable_tp_branch_current(pm)
    variable_tp_branch_flow(pm)

    for c in PMs.conductor_ids(pm)
        PMs.variable_generation(pm, cnd=c)
        PMs.variable_dcline_flow(pm, cnd=c)
    end

    # Constraints
    for i in ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_tp_kcl_shunt(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_tp_flow_losses(pm, i)

        constraint_tp_voltage_magnitude_difference(pm, i)

        constraint_tp_branch_current(pm, i)


        for c in PMs.conductor_ids(pm)
            constraint_voltage_angle_difference(pm, i, cnd=c)

            PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end

    for i in ids(pm, :dcline)
        for h in PMs.phase_ids(pm)
            PMs.constraint_dcline(pm, i, ph=h)
        end
    end

     PMs.objective_min_fuel_cost(pm)
end
