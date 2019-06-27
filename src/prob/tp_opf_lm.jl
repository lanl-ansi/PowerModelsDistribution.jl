export run_tp_opf_lm, run_ac_tp_opf_lm

""
function run_ac_tp_opf_lm(file, solver; kwargs...)
    return run_tp_opf_lm(file, PMs.ACPPowerModel, solver; multiconductor=true, kwargs...)
end


""
function run_tp_opf_lm(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf; multiconductor=true, kwargs...)
end


""
function run_tp_opf_lm(file::String, model_constructor, solver; kwargs...)
    data = ThreePhasePowerModels.parse_file(file)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf; multiconductor=true, kwargs...)
end


"""
This problem specification includes advanced load models, including
constant power, constant current and constabt impedance
delta-connected and wye-connected
"""
function post_tp_opf_lm(pm::PMs.GenericPowerModel)
    add_arcs_trans!(pm)

    variable_tp_voltage(pm)
    variable_tp_branch_flow(pm)

    for c in PMs.conductor_ids(pm)
        PMs.variable_generation(pm, cnd=c)
        variable_load(pm, cnd=c)
        PMs.variable_dcline_flow(pm, cnd=c)
    end
    variable_tp_trans_flow(pm)

    constraint_tp_voltage(pm)

    for i in PMs.ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in PMs.ids(pm, :load)
        constraint_tp_load(pm, id)
    end

    for i in PMs.ids(pm, :bus), c in PMs.conductor_ids(pm)
        constraint_kcl_shunt_trans_load(pm, i, cnd=c)
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
