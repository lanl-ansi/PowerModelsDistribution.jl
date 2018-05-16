export run_ac_tp_pf, run_tp_pf


""
function run_ac_tp_pf(file, solver; kwargs...)
    return run_tp_pf(file, PMs.ACPPowerModel, solver; kwargs...)
end


""
function run_dc_tp_pf(file, solver; kwargs...)
    return run_tp_pf(file, PMs.DCPPowerModel, solver; kwargs...)
end


""
function run_tp_pf(file, model_constructor, solver; kwargs...)
    return run_generic_tp_model(file, model_constructor, solver, post_tp_pf; multiphase=true, kwargs...)
end


""
function post_tp_pf(pm::GenericPowerModel)
    for h in PMs.phase_ids(pm)
        PMs.variable_voltage(pm, bounded=false, ph=h)
        PMs.variable_generation(pm, bounded=false, ph=h)
        PMs.variable_branch_flow(pm, bounded=false, ph=h)
        PMs.variable_dcline_flow(pm, bounded=false, ph=h)

        PMs.constraint_voltage(pm, ph=h)

        for (i,bus) in ref(pm, :ref_buses)
            @assert bus["bus_type"] == 3
            PMs.constraint_theta_ref(pm, i, ph=h)
            PMs.constraint_voltage_magnitude_setpoint(pm, i, ph=h)
        end

        for (i,bus) in ref(pm, :bus)
            PMs.constraint_kcl_shunt(pm, i)

            # PV Bus Constraints
            if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
                # this assumes inactive generators are filtered out of bus_gens
                @assert bus["bus_type"] == 2

                PMs.constraint_voltage_magnitude_setpoint(pm, i, ph=h)
                for j in ref(pm, :bus_gens, i)
                    PMs.constraint_active_gen_setpoint(pm, j, ph=h)
                end
            end
        end

        for i in ids(pm, :branch)
            constraint_ohms_yt_from(pm, i, ph=h)
            constraint_ohms_yt_to(pm, i, ph=h)
        end

        for (i,dcline) in ref(pm, :dcline)
            PMs.constraint_active_dcline_setpoint(pm, i, ph=h)

            f_bus = ref(pm, :bus)[dcline["f_bus"]]
            if f_bus["bus_type"] == 1
                PMs.constraint_voltage_magnitude_setpoint(pm, f_bus["index"], ph=h)
            end

            t_bus = ref(pm, :bus)[dcline["t_bus"]]
            if t_bus["bus_type"] == 1
                PMs.constraint_voltage_magnitude_setpoint(pm, t_bus["index"], ph=h)
            end
        end

    end
end