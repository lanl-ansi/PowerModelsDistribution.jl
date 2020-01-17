######
#
# These are toy problem formulations used to test advanced features
# such as storage devices
#
######
"multi-network opf with storage"
function _run_mn_mc_opf(data::Dict{String,Any}, model_type, solver; kwargs...)
    return _PMs.run_model(data, model_type, solver, _build_mn_mc_strg_opf; ref_extensions=[ref_add_arcs_trans!], multiconductor=true, multinetwork=true, kwargs...)
end


"multi-network opf with storage"
function _run_mn_mc_opf(file::String, model_type, solver; kwargs...)
    return run_mn_mc_opf(PowerModelsDistribution.parse_file(file), model_type, solver; kwargs...)
end


"multi-network opf with storage"
function _build_mn_mc_strg_opf(pm::_PMs.AbstractPowerModel)
    for (n, network) in _PMs.nws(pm)
        variable_mc_voltage(pm; nw=n)
        constraint_mc_model_voltage(pm; nw=n)
        variable_mc_branch_flow(pm; nw=n)
        variable_mc_generation(pm; nw=n)
        variable_mc_storage(pm; nw=n)

        for i in _PMs.ids(pm, :ref_buses; nw=n)
            constraint_mc_theta_ref(pm, i; nw=n)
        end

        for i in _PMs.ids(pm, :bus; nw=n)
            constraint_mc_power_balance(pm, i; nw=n)
        end

        for i in _PMs.ids(pm, :storage; nw=n)
            _PMs.constraint_storage_state(pm, i; nw=n)
            _PMs.constraint_storage_complementarity_nl(pm, i; nw=n)
            constraint_mc_storage_loss(pm, i; nw=n)
            constraint_mc_storage_thermal_limit(pm, i; nw=n)
        end

        for i in _PMs.ids(pm, :branch; nw=n)
            constraint_mc_ohms_yt_from(pm, i; nw=n)
            constraint_mc_ohms_yt_to(pm, i; nw=n)

            constraint_mc_voltage_angle_difference(pm, i; nw=n)

            constraint_mc_thermal_limit_from(pm, i; nw=n)
            constraint_mc_thermal_limit_to(pm, i; nw=n)
        end

        for i in _PMs.ids(pm, :transformer; nw=n)
            constraint_mc_trans(pm, i; nw=n)
        end
    end

    network_ids = sort(collect(_PMs.nw_ids(pm)))

    n_1 = network_ids[1]
    for i in _PMs.ids(pm, :storage; nw=n_1)
        _PMs.constraint_storage_state(pm, i; nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in _PMs.ids(pm, :storage; nw=n_2)
            _PMs.constraint_storage_state(pm, i, n_1, n_2)
        end
        n_1 = n_2
    end

    _PMs.objective_min_fuel_cost(pm)
end


######
#
# Formulations from PowerModels
#
######

""
function _run_mc_ucopf(file, model_type::Type, solver; kwargs...)
    return run_model(file, model_type, solver, _build_mc_ucopf; solution_builder = _solution_uc!, multiconductor=true, kwargs...)
end

""
function _build_mc_ucopf(pm::_PMs.AbstractPowerModel)
    variable_generation_indicator(pm)

    variable_storage_indicator(pm)
    variable_storage_energy(pm)
    variable_storage_charge(pm)
    variable_storage_discharge(pm)
    variable_storage_complementary_indicator(pm)

    for c in conductor_ids(pm)
        variable_voltage(pm, cnd=c)

        variable_generation_on_off(pm, cnd=c)

        variable_active_storage(pm, cnd=c)
        variable_reactive_storage(pm, cnd=c)

        variable_branch_flow(pm, cnd=c)
        variable_dcline_flow(pm, cnd=c)

        constraint_model_voltage(pm, cnd=c)

        for i in ids(pm, :ref_buses)
            constraint_theta_ref(pm, i, cnd=c)
        end

        for i in ids(pm, :gen)
            constraint_generation_on_off(pm, i, cnd=c)
        end

        for i in ids(pm, :bus)
            constraint_power_balance(pm, i, cnd=c)
        end

        for i in ids(pm, :branch)
            constraint_ohms_yt_from(pm, i, cnd=c)
            constraint_ohms_yt_to(pm, i, cnd=c)

            constraint_voltage_angle_difference(pm, i, cnd=c)

            constraint_thermal_limit_from(pm, i, cnd=c)
            constraint_thermal_limit_to(pm, i, cnd=c)
        end

        for i in ids(pm, :dcline)
            constraint_dcline(pm, i, cnd=c)
        end
    end

    for i in ids(pm, :storage)
        constraint_storage_state(pm, i)
        constraint_storage_complementarity_mi(pm, i)
        constraint_storage_loss(pm, i, conductors=conductor_ids(pm))

        for c in conductor_ids(pm)
            constraint_storage_on_off(pm, i, cnd=c)
            constraint_storage_thermal_limit(pm, i, cnd=c)
        end
    end

    objective_min_fuel_and_flow_cost(pm)
end



""
function _run_mc_opf(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, _build_mc_opf; multiconductor=true, kwargs...)
end

""
function _build_mc_opf(pm::_PMs.AbstractPowerModel)
    for c in conductor_ids(pm)
        variable_voltage(pm, cnd=c)
        variable_generation(pm, cnd=c)
        variable_branch_flow(pm, cnd=c)
        variable_dcline_flow(pm, cnd=c)

        constraint_model_voltage(pm, cnd=c)

        for i in ids(pm, :ref_buses)
            constraint_theta_ref(pm, i, cnd=c)
        end

        for i in ids(pm, :bus)
            constraint_power_balance(pm, i, cnd=c)
        end

        for i in ids(pm, :branch)
            constraint_ohms_yt_from(pm, i, cnd=c)
            constraint_ohms_yt_to(pm, i, cnd=c)

            constraint_voltage_angle_difference(pm, i, cnd=c)

            constraint_thermal_limit_from(pm, i, cnd=c)
            constraint_thermal_limit_to(pm, i, cnd=c)
        end

        for i in ids(pm, :dcline)
            constraint_dcline(pm, i, cnd=c)
        end
    end

    objective_min_fuel_and_flow_cost(pm)
end


""
function _run_mc_opf_iv(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, _build_mc_opf_iv; multiconductor=true, kwargs...)
end

""
function _build_mc_opf_iv(pm::_PMs.AbstractPowerModel)
    for c in conductor_ids(pm)
        variable_voltage(pm, cnd=c)
        variable_branch_current(pm, cnd=c)

        variable_gen(pm, cnd=c)
        variable_dcline(pm, cnd=c)


        for i in ids(pm, :ref_buses)
            constraint_theta_ref(pm, i, cnd=c)
        end

        for i in ids(pm, :bus)
            constraint_current_balance(pm, i, cnd=c)
        end

        for i in ids(pm, :branch)
            constraint_current_from(pm, i, cnd=c)
            constraint_current_to(pm, i, cnd=c)

            constraint_voltage_drop(pm, i, cnd=c)
            constraint_voltage_angle_difference(pm, i, cnd=c)

            constraint_thermal_limit_from(pm, i, cnd=c)
            constraint_thermal_limit_to(pm, i, cnd=c)
        end

        for i in ids(pm, :dcline)
            constraint_dcline(pm, i, cnd=c)
        end
    end

    objective_min_fuel_and_flow_cost(pm)
end


""
function _run_mn_mc_opf(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, _build_mn_mc_opf; multinetwork=true, multiconductor=true, kwargs...)
end

""
function _build_mn_mc_opf(pm::_PMs.AbstractPowerModel)
    for (n, network) in nws(pm)
        for c in conductor_ids(pm, nw=n)
            variable_voltage(pm, nw=n, cnd=c)
            variable_generation(pm, nw=n, cnd=c)
            variable_branch_flow(pm, nw=n, cnd=c)
            variable_dcline_flow(pm, nw=n, cnd=c)

            constraint_model_voltage(pm, nw=n, cnd=c)

            for i in ids(pm, :ref_buses, nw=n)
                constraint_theta_ref(pm, i, nw=n, cnd=c)
            end

            for i in ids(pm, :bus, nw=n)
                constraint_power_balance(pm, i, nw=n, cnd=c)
            end

            for i in ids(pm, :branch, nw=n)
                constraint_ohms_yt_from(pm, i, nw=n, cnd=c)
                constraint_ohms_yt_to(pm, i, nw=n, cnd=c)

                constraint_voltage_angle_difference(pm, i, nw=n, cnd=c)

                constraint_thermal_limit_from(pm, i, nw=n, cnd=c)
                constraint_thermal_limit_to(pm, i, nw=n, cnd=c)
            end

            for i in ids(pm, :dcline, nw=n)
                constraint_dcline(pm, i, nw=n, cnd=c)
            end
        end
    end

    objective_min_fuel_and_flow_cost(pm)
end


""
function _run_mn_mc_opf_strg(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, _build_mn_mc_opf_strg; multinetwork=true, multiconductor=true, kwargs...)
end

"warning: this model is not realistic or physically reasonable, it is only for test coverage"
function _build_mn_mc_opf_strg(pm::_PMs.AbstractPowerModel)
    for (n, network) in nws(pm)
        variable_storage_energy(pm, nw=n)
        variable_storage_charge(pm, nw=n)
        variable_storage_discharge(pm, nw=n)

        for c in conductor_ids(pm, nw=n)
            variable_voltage(pm, nw=n, cnd=c)
            variable_generation(pm, nw=n, cnd=c)
            variable_active_storage(pm, nw=n, cnd=c)
            variable_reactive_storage(pm, nw=n, cnd=c)
            variable_branch_flow(pm, nw=n, cnd=c)
            variable_dcline_flow(pm, nw=n, cnd=c)

            constraint_model_voltage(pm, nw=n, cnd=c)

            for i in ids(pm, :ref_buses, nw=n)
                constraint_theta_ref(pm, i, nw=n, cnd=c)
            end

            for i in ids(pm, :bus, nw=n)
                constraint_power_balance(pm, i, nw=n, cnd=c)
            end

            for i in ids(pm, :storage, nw=n)
                constraint_storage_thermal_limit(pm, i, nw=n, cnd=c)
            end

            for i in ids(pm, :branch, nw=n)
                constraint_ohms_yt_from(pm, i, nw=n, cnd=c)
                constraint_ohms_yt_to(pm, i, nw=n, cnd=c)

                constraint_voltage_angle_difference(pm, i, nw=n, cnd=c)

                constraint_thermal_limit_from(pm, i, nw=n, cnd=c)
                constraint_thermal_limit_to(pm, i, nw=n, cnd=c)
            end

            for i in ids(pm, :dcline, nw=n)
                constraint_dcline(pm, i, nw=n, cnd=c)
            end
        end

        for i in ids(pm, :storage, nw=n)
            constraint_storage_complementarity_nl(pm, i, nw=n)
            constraint_storage_loss(pm, i, nw=n, conductors=conductor_ids(pm, nw=n))
        end
    end

    network_ids = sort(collect(nw_ids(pm)))

    n_1 = network_ids[1]
    for i in ids(pm, :storage, nw=n_1)
        constraint_storage_state(pm, i, nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in ids(pm, :storage, nw=n_2)
            constraint_storage_state(pm, i, n_1, n_2)
        end
        n_1 = n_2
    end

    objective_min_fuel_and_flow_cost(pm)
end
