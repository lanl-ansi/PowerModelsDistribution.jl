######
#
# These are toy problem formulations used to test advanced features
# such as storage devices
#
######
"multi-network opf with storage"
function run_mn_mc_opf(data::Dict{String,Any}, model_type, solver; kwargs...)
    return _PMs.run_model(data, model_type, solver, post_mn_mc_strg_opf; multiconductor=true, multinetwork=true, kwargs...)
end


""
function run_mn_mc_opf(file::String, model_type, solver; kwargs...)
    return run_mn_mc_opf(PowerModelsDistribution.parse_file(file), model_type, solver; kwargs...)
end


""
function post_mn_mc_strg_opf(pm::_PMs.AbstractPowerModel)
    for (n, network) in _PMs.nws(pm)
        variable_mc_voltage(pm, nw=n)
        constraint_mc_model_voltage(pm, nw=n)
        variable_mc_branch_flow(pm, nw=n)
        variable_mc_storage(pm, nw=n)

        for c in _PMs.conductor_ids(pm, nw=n)
            _PMs.variable_generation(pm, cnd=c, nw=n)
            _PMs.variable_dcline_flow(pm, cnd=c, nw=n)
        end

        for i in _PMs.ids(pm, :ref_buses, nw=n)
            constraint_mc_theta_ref(pm, i, nw=n)
        end

        for i in _PMs.ids(pm, :bus, nw=n), c in _PMs.conductor_ids(pm, nw=n)
            _PMs.constraint_power_balance(pm, i, cnd=c, nw=n)
        end

        for i in _PMs.ids(pm, :storage, nw=n)
            _PMs.constraint_storage_state(pm, i, nw=n)
            _PMs.constraint_storage_complementarity_nl(pm, i, nw=n)
            _PMs.constraint_storage_loss(pm, i, conductors=_PMs.conductor_ids(pm), nw=n)
            for c in _PMs.conductor_ids(pm, nw=n)
                _PMs.constraint_storage_thermal_limit(pm, i, cnd=c, nw=n)
            end
        end

        for i in _PMs.ids(pm, :branch, nw=n)
            constraint_mc_ohms_yt_from(pm, i, nw=n)
            constraint_mc_ohms_yt_to(pm, i, nw=n)

            for c in _PMs.conductor_ids(pm, nw=n)
                _PMs.constraint_voltage_angle_difference(pm, i, cnd=c, nw=n)

                _PMs.constraint_thermal_limit_from(pm, i, cnd=c, nw=n)
                _PMs.constraint_thermal_limit_to(pm, i, cnd=c, nw=n)
            end
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
