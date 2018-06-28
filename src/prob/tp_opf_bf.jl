export run_tp_opf_bf

""
function run_tp_opf_bf(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    if model_constructor != PMs.SOCBFPowerModel
        error(LOGGER, "The problem type opf_bf at the moment only supports the SOCBFForm formulation")
    end
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf_bf; multiphase=true, kwargs...)
end


""
function run_tp_opf_bf(file::String, model_constructor, solver; kwargs...)
    if model_constructor != PMs.SOCBFPowerModel
        error(LOGGER, "The problem type opf_bf at the moment only supports the SOCBFForm formulation")
    end
    data = ThreePhasePowerModels.parse_file(file)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf_bf; multiphase=true, kwargs...)
end


""
function post_tp_opf_bf(pm::GenericPowerModel)
    for h in PMs.phase_ids(pm)
        variable_tp_voltage(pm, ph=h)
        PMs.variable_generation(pm, ph=h)
        PMs.variable_branch_flow(pm, ph=h)
        PMs.variable_branch_current(pm, ph=h)
        PMs.variable_dcline_flow(pm, ph=h)
    end

    for h in PMs.phase_ids(pm)
        for i in ids(pm, :ref_buses)
            constraint_tp_theta_ref(pm, i, ph=h)
        end

        for i in ids(pm, :bus)
            PMs.constraint_kcl_shunt(pm, i, ph=h)
        end

        for i in ids(pm, :branch)
            PMs.constraint_flow_losses(pm, i, ph=h)

            constraint_tp_voltage_magnitude_difference(pm, i, ph=h)

            PMs.constraint_branch_current(pm, i, ph=h)

            PMs.constraint_voltage_angle_difference(pm, i, ph=h)

            PMs.constraint_thermal_limit_from(pm, i, ph=h)
            PMs.constraint_thermal_limit_to(pm, i, ph=h)
        end

        for i in ids(pm, :dcline)
            PMs.constraint_dcline(pm, i, ph=h)
        end
    end

    PMs.objective_min_fuel_cost(pm)
end