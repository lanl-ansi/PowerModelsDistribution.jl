export run_tp_opf, run_ac_tp_opf

""
function run_ac_tp_opf(file, solver; kwargs...)
    return run_tp_opf(file, PMs.ACPPowerModel, solver; multiphase=true, kwargs...)
end


""
function run_tp_opf(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf; multiphase=true, kwargs...)
end


""
function run_tp_opf(file::String, model_constructor, solver; kwargs...)
    data = ThreePhasePowerModels.parse_file(file)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf; multiphase=true, kwargs...)
end


""
function post_tp_opf(pm::GenericPowerModel)
    for h in PMs.phase_ids(pm)
        variable_tp_voltage(pm, ph=h)
        PMs.variable_generation(pm, ph=h)
        PMs.variable_branch_flow(pm, ph=h)
        PMs.variable_dcline_flow(pm, ph=h)
    end

    for h in PMs.phase_ids(pm)
        constraint_tp_voltage(pm, ph=h)

        for i in ids(pm, :ref_buses)
            constraint_tp_theta_ref(pm, i, ph=h)
        end

        for i in ids(pm, :bus)
            PMs.constraint_kcl_shunt(pm, i, ph=h)
        end

        for i in ids(pm, :branch)
            constraint_ohms_tp_yt_from(pm, i, ph=h)
            constraint_ohms_tp_yt_to(pm, i, ph=h)

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
