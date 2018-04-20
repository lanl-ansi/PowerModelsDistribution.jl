export run_tp_opf, run_ac_tp_opf

""
function run_ac_tp_opf(file, solver; kwargs...)
    return run_tp_opf(file, PMs.ACPPowerModel, solver; kwargs...)
end

""
function run_tp_opf(file, model_constructor, solver; kwargs...)
    return PMs.run_generic_model(file, model_constructor, solver, post_tp_opf; kwargs...)
end

""
function post_tp_opf(pm::GenericPowerModel)
    for (n, nw_ref) in pm.ref[:nw]
        for (h, ph_ref) in nw_ref[:ph]
            PMs.variable_voltage(pm, n, h)
            PMs.variable_generation(pm, n, h)
            PMs.variable_branch_flow(pm, n, h)
            PMs.variable_dcline_flow(pm, n, h)

            PMs.constraint_voltage(pm, n, h)

            for i in ids(pm, n, p, :ref_buses)
                PMs.constraint_theta_ref(pm, n, h, i)
            end

            for i in ids(pm, n, p, :bus)
                PMs.constraint_kcl_shunt(pm, n, h, i)
            end

            for i in ids(pm, n, p, :branch)
                PMs.constraint_ohms_yt_from(pm, n, h, i)
                PMs.constraint_ohms_yt_to(pm, n, h, i)

                PMs.constraint_voltage_angle_difference(pm, n, h, i)

                PMs.constraint_thermal_limit_from(pm, n, h, i)
                PMs.constraint_thermal_limit_to(pm, n, h, i)
            end

            for i in ids(pm, n, h, :dcline)
                PMs.constraint_tp_dcline(pm, n, h, i)
            end
        end
    end

    PMs.objective_min_fuel_cost(pm)
end
