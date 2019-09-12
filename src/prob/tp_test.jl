######
#
# These are toy problem formulations used to test advanced features
# such as storage devices
#
######

"opf with storage"
function run_tp_strg_opf(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, post_tp_strg_opf; multiconductor=true, ref_extensions=[ref_add_arcs_trans!], kwargs...)
end


""
function run_tp_strg_opf(file::String, model_constructor, solver; kwargs...)
    return run_tp_strg_opf(PowerModelsDistribution.parse_file(file), model_constructor, solver; kwargs...)
end


""
function post_tp_strg_opf(pm::_PMs.GenericPowerModel)
    variable_tp_voltage(pm)
    variable_tp_branch_flow(pm)
    variable_tp_storage(pm)

    for c in _PMs.conductor_ids(pm)
        _PMs.variable_generation(pm, cnd=c)
        _PMs.variable_dcline_flow(pm, cnd=c)
    end

    variable_tp_trans_flow(pm)
    constraint_tp_model_voltage(pm)

    for i in _PMs.ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    for i in _PMs.ids(pm, :bus), c in _PMs.conductor_ids(pm)
        constraint_tp_power_balance_shunt_storage_trans(pm, i, cnd=c)
    end

    for i in _PMs.ids(pm, :storage)
        _PMs.constraint_storage_state(pm, i)
        constraint_tp_storage_exchange(pm, i)
        for c in _PMs.conductor_ids(pm)
            _PMs.constraint_storage_thermal_limit(pm, i, cnd=c)
        end
    end

    for i in _PMs.ids(pm, :branch)
        constraint_mc_voltage_angle_difference(pm, i)

        for c in _PMs.conductor_ids(pm)
            constraint_tp_ohms_yt_from(pm, i, cnd=c)
            constraint_tp_ohms_yt_to(pm, i, cnd=c)

            _PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            _PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end

    for i in _PMs.ids(pm, :dcline), c in _PMs.conductor_ids(pm)
        _PMs.constraint_dcline(pm, i, cnd=c)
    end

    for i in _PMs.ids(pm, :trans)
        constraint_tp_trans(pm, i)
    end

    _PMs.objective_min_fuel_cost(pm)
end




"multi-network opf with storage"
function run_mn_tp_strg_opf(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return _PMs.run_model(data, model_constructor, solver, post_mn_tp_strg_opf; multiconductor=true, multinetwork=true, kwargs...)
end


""
function run_mn_tp_strg_opf(file::String, model_constructor, solver; kwargs...)
    return run_mn_tp_strg_opf(PowerModelsDistribution.parse_file(file), model_constructor, solver; kwargs...)
end


""
function post_mn_tp_strg_opf(pm::_PMs.GenericPowerModel)
    for (n, network) in _PMs.nws(pm)
        variable_tp_voltage(pm, nw=n)
        variable_tp_branch_flow(pm, nw=n)
        variable_tp_storage(pm, nw=n)

        for c in _PMs.conductor_ids(pm, nw=n)
            _PMs.variable_generation(pm, cnd=c, nw=n)
            _PMs.variable_dcline_flow(pm, cnd=c, nw=n)
        end

        constraint_tp_model_voltage(pm, nw=n)

        for i in _PMs.ids(pm, :ref_buses, nw=n)
            constraint_tp_theta_ref(pm, i, nw=n)
        end

        for i in _PMs.ids(pm, :bus, nw=n), c in _PMs.conductor_ids(pm, nw=n)
            _PMs.constraint_power_balance_shunt_storage(pm, i, cnd=c, nw=n)
        end

        for i in _PMs.ids(pm, :storage, nw=n)
            constraint_tp_storage_exchange(pm, i, nw=n)
            for c in _PMs.conductor_ids(pm, nw=n)
                _PMs.constraint_storage_thermal_limit(pm, i, cnd=c, nw=n)
            end
        end

        for i in _PMs.ids(pm, :branch, nw=n)
            for c in _PMs.conductor_ids(pm, nw=n)
                constraint_tp_ohms_yt_from(pm, i, cnd=c, nw=n)
                constraint_tp_ohms_yt_to(pm, i, cnd=c, nw=n)

                _PMs.constraint_voltage_angle_difference(pm, i, cnd=c, nw=n)

                _PMs.constraint_thermal_limit_from(pm, i, cnd=c, nw=n)
                _PMs.constraint_thermal_limit_to(pm, i, cnd=c, nw=n)
            end
        end

        for i in _PMs.ids(pm, :dcline, nw=n), c in _PMs.conductor_ids(pm, nw=n)
            _PMs.constraint_dcline(pm, i, cnd=c, nw=n)
        end
    end

    network_ids = sort(collect(_PMs.nw_ids(pm)))

    n_1 = network_ids[1]
    for i in _PMs.ids(pm, :storage, nw=n_1)
        _PMs.constraint_storage_state(pm, i, nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in _PMs.ids(pm, :storage, nw=n_2)
            _PMs.constraint_storage_state(pm, i, n_1, n_2)
        end
        n_1 = n_2
    end

    _PMs.objective_min_fuel_cost(pm)
end
