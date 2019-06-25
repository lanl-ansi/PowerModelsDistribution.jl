export run_ac_tp_pf_lm, run_tp_pf_lm


""
function run_ac_tp_pf_lm(file, solver; kwargs...)
    return run_tp_pf_lm(file, PMs.ACPPowerModel, solver; kwargs...)
end


""
function run_dc_tp_opf_lm(file, solver; kwargs...)
    return run_tp_pf_lm(file, PMs.DCPPowerModel, solver; kwargs...)
end


""
function run_tp_pf_lm(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_pf; multiconductor=true, kwargs...)
end


""
function run_tp_pf_lm(file::String, model_constructor, solver; kwargs...)
    data = ThreePhasePowerModels.parse_file(file)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_pf; multiconductor=true, kwargs...)
end


""
function post_tp_pf_lm(pm::PMs.GenericPowerModel)
    add_arcs_trans!(pm)

    variable_tp_voltage(pm, bounded=false)
    variable_tp_branch_flow(pm, bounded=false)

    for c in PMs.conductor_ids(pm)
        PMs.variable_generation(pm, bounded=false, cnd=c)
        variable_load(pm, bounded=false, cnd=c)
        PMs.variable_dcline_flow(pm, bounded=false, cnd=c)
    end

    variable_tp_trans_flow(pm, bounded=false)

    constraint_tp_voltage(pm)

    for (i,bus) in PMs.ref(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)

        for c in PMs.conductor_ids(pm)
            @assert bus["bus_type"] == 3
            PMs.constraint_voltage_magnitude_setpoint(pm, i, cnd=c)
        end
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in PMs.ids(pm, :load)
        constraint_tp_load(pm, id)
    end

    for (i,bus) in PMs.ref(pm, :bus), c in PMs.conductor_ids(pm)
        constraint_kcl_shunt_trans_load(pm, i, cnd=c)

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

    for i in PMs.ids(pm, :branch), c in PMs.conductor_ids(pm)
        constraint_ohms_tp_yt_from(pm, i, cnd=c)
        constraint_ohms_tp_yt_to(pm, i, cnd=c)
        # PMs.constraint_ohms_yt_from(pm, i, cnd=c)
        # PMs.constraint_ohms_yt_to(pm, i, cnd=c)
    end

    for (i,dcline) in PMs.ref(pm, :dcline), c in PMs.conductor_ids(pm)

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

    for i in PMs.ids(pm, :trans)
        constraint_tp_trans(pm, i)
    end

end
