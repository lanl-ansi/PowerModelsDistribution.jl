export run_tp_pf_bf

""
function run_tp_pf_bf(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    if model_constructor != SDPUBFPowerModel && model_constructor != SOCNLPUBFPowerModel && model_constructor != SOCConicUBFPowerModel && model_constructor != LPUBFPowerModel && model_constructor != LPdiagUBFPowerModel && model_constructor !=  SOCBFPowerModel
        Memento.error(LOGGER, "The problem type tp_opf_bf at the moment only supports a limited set of formulations")
    end
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_pf_bf; solution_builder=get_solution_tp, multiconductor=true, kwargs...)
end


""
function run_tp_pf_bf(file::String, model_constructor, solver; kwargs...)
    if model_constructor != SDPUBFPowerModel && model_constructor != SOCNLPUBFPowerModel && model_constructor != SOCConicUBFPowerModel && model_constructor != LPUBFPowerModel && model_constructor !=  SOCBFPowerModel
        Memento.error(LOGGER, "The problem type tp_opf_bf at the moment only supports a limited set of formulations")
    end
    data = ThreePhasePowerModels.parse_file(file)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_pf_bf; solution_builder=get_solution_tp, multiconductor=true, kwargs...)
end


""
function post_tp_pf_bf(pm::PMs.GenericPowerModel)
    # Variables
    variable_tp_voltage(pm, bounded=false)
    variable_tp_branch_current(pm)
    variable_tp_branch_flow(pm)

    for c in PMs.conductor_ids(pm)
        PMs.variable_generation(pm, bounded=false, cnd=c)
        PMs.variable_dcline_flow(pm, bounded=false, cnd=c)
    end

    # Constraints
    for (i,bus) in PMs.ref(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)

        for c in PMs.conductor_ids(pm)
            @assert bus["bus_type"] == 3
            # PMs.constraint_voltage_magnitude_setpoint(pm, i, cnd=c) #TODO add back
        end
    end

    for i in PMs.ids(pm, :bus), c in PMs.conductor_ids(pm)
        PMs.constraint_kcl_shunt(pm, i, cnd=c)

        # PV Bus Constraints
        if length(PMs.ref(pm, :bus_gens, i)) > 0 && !(i in PMs.ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            PMs.constraint_voltage_magnitude_setpoint(pm, i, cnd=c)
            for j in PMs.ref(pm, :bus_gens, i)
                PMs.constraint_active_gen_setpoint(pm, j, cnd=c)
            end
        end
    end

    for i in PMs.ids(pm, :branch)
        constraint_tp_flow_losses(pm, i)

        constraint_tp_voltage_magnitude_difference(pm, i)

        constraint_tp_branch_current(pm, i)


        for c in PMs.conductor_ids(pm)
            constraint_voltage_angle_difference(pm, i, cnd=c)

            PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end

    for (i,dcline) in PMs.ref(pm, :dcline)
        for c in PMs.conductor_ids(pm)

            PMs.constraint_active_dcline_setpoint(pm, i, cnd=c)

            f_bus = PMs.ref(pm, :bus)[dcline["f_bus"]]
            if f_bus["bus_type"] == 1
                PMs.constraint_voltage_magnitude_setpoint(pm, f_bus["index"], cnd=c)
            end

            t_bus = PMs.ref(pm, :bus)[dcline["t_bus"]]
            if t_bus["bus_type"] == 1
                PMs.constraint_voltage_magnitude_setpoint(pm, t_bus["index"], cnd=c)
            end
        end
    end

     PMs.objective_min_fuel_cost(pm)
end
