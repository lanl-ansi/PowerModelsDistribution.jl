######
#
# These are toy problem formulations used to test advanced features
# such as storage devices
#
######

"opf with storage"
function run_tp_strg_opf(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_strg_opf; multiconductor=true, kwargs...)
end

""
function post_tp_strg_opf(pm::GenericPowerModel)
    variable_tp_voltage(pm)
    variable_tp_branch_flow(pm)
    variable_tp_storage(pm)

    for c in PMs.conductor_ids(pm)
        PMs.variable_generation(pm, cnd=c)
        PMs.variable_dcline_flow(pm, cnd=c)
    end

    constraint_tp_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    for i in ids(pm, :bus), c in PMs.conductor_ids(pm)
        PMs.constraint_kcl_shunt_storage(pm, i, cnd=c)
    end

    for i in ids(pm, :storage)
        PMs.constraint_storage_state(pm, i)
        constraint_tp_storage_exchange(pm, i)
        for c in PMs.conductor_ids(pm)
            PMs.constraint_storage_thermal_limit(pm, i, cnd=c)
        end
    end

    for i in ids(pm, :branch)
        for c in PMs.conductor_ids(pm)
            constraint_ohms_tp_yt_from(pm, i, cnd=c)
            constraint_ohms_tp_yt_to(pm, i, cnd=c)

            PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

            PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end

    for i in ids(pm, :dcline), c in PMs.conductor_ids(pm)
        PMs.constraint_dcline(pm, i, cnd=c)
    end

    PMs.objective_min_fuel_cost(pm)
end
