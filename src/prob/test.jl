######
#
# Formulations from PowerModels
#
######

""
function _run_mc_ucopf(file, model_type::Type, solver; kwargs...)
    return run_mc_model(file, model_type, solver, _build_mc_ucopf; kwargs...)
end

""
function _build_mc_ucopf(pm::_PM.AbstractPowerModel)
    for (n, network) in nws(pm)
        variable_mc_bus_voltage(pm, nw=n)
        variable_mc_branch_power(pm, nw=n)
        variable_mc_transformer_power(pm, nw=n)

        variable_mc_gen_indicator(pm, nw=n)
        variable_mc_gen_power_setpoint_on_off(pm, nw=n)


        constraint_mc_model_voltage(pm, nw=n)



        variable_mc_storage_power_on_off(pm, nw=n)
        _PM.variable_storage_energy(pm, nw=n)
        _PM.variable_storage_charge(pm, nw=n)
        _PM.variable_storage_discharge(pm, nw=n)
        _PM.variable_storage_indicator(pm, nw=n)
        _PM.variable_storage_complementary_indicator(pm, nw=n)


        for i in ids(pm, :ref_buses, nw=n)
            constraint_mc_theta_ref(pm, i, nw=n)
        end

        for i in ids(pm, :bus, nw=n)
            constraint_mc_power_balance(pm, i, nw=n)

        end

        for i in ids(pm, :branch, nw=n)
            constraint_mc_ohms_yt_from(pm, i, nw=n)
            constraint_mc_ohms_yt_to(pm, i, nw=n)

            constraint_mc_voltage_angle_difference(pm, i, nw=n)

            constraint_mc_thermal_limit_from(pm, i, nw=n)
            constraint_mc_thermal_limit_to(pm, i, nw=n)
        end


        for i in ids(pm, :transformer, nw=n)
            constraint_mc_transformer_power(pm, i, nw=n)
        end

        # for i in ids(pm, :dcline, nw=n)
        #     constraint_mc_dcline(pm, i, nw=n)
        # end

        for i in ids(pm, :storage; nw=n)
            # _PM.constraint_storage_state(pm, i; nw=n)
            _PM.constraint_storage_complementarity_mi(pm, i; nw=n)
            constraint_mc_storage_losses(pm, i; nw=n)
            constraint_mc_storage_thermal_limit(pm, i; nw=n)

            constraint_mc_storage_on_off(pm, i; nw=n)

        end
    end

    network_ids = sort(collect(nw_ids(pm)))

    n_1 = network_ids[1]
    for i in ids(pm, :storage, nw=n_1)
        _PM.constraint_storage_state(pm, i, nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in ids(pm, :storage, nw=n_2)
            _PM.constraint_storage_state(pm, i, n_1, n_2)
        end
        n_1 = n_2
    end
    _PM.objective_min_fuel_cost(pm)
end


""
function _run_mn_mc_opf(file, model_type::Type, optimizer; kwargs...)
    return run_mc_model(file, model_type, optimizer, _build_mn_mc_opf; multinetwork=true, kwargs...)
end

""
function _build_mn_mc_opf(pm::_PM.AbstractPowerModel)
    for (n, network) in nws(pm)
        variable_mc_bus_voltage(pm, nw=n)
        variable_mc_branch_power(pm, nw=n)
        variable_mc_transformer_power(pm, nw=n)
        variable_mc_gen_power_setpoint(pm, nw=n)

        constraint_mc_model_voltage(pm, nw=n)

        for i in ids(pm, :ref_buses, nw=n)
            constraint_mc_theta_ref(pm, i, nw=n)
        end

        for i in ids(pm, :bus, nw=n)
            constraint_mc_power_balance(pm, i, nw=n)

        end

        for i in ids(pm, :branch, nw=n)
            constraint_mc_ohms_yt_from(pm, i, nw=n)
            constraint_mc_ohms_yt_to(pm, i, nw=n)

            constraint_mc_voltage_angle_difference(pm, i, nw=n)

            constraint_mc_thermal_limit_from(pm, i, nw=n)
            constraint_mc_thermal_limit_to(pm, i, nw=n)
        end


        for i in ids(pm, :transformer, nw=n)
            constraint_mc_transformer_power(pm, i, nw=n)
        end

        # for i in ids(pm, :dcline, nw=n)
        #     constraint_mc_dcline(pm, i, nw=n)
        # end
    end
    _PM.objective_min_fuel_cost(pm)
end


""
function _run_mn_mc_opf_strg(file, model_type::Type, optimizer; kwargs...)
    return run_mc_model(file, model_type, optimizer, _build_mn_mc_opf_strg; multinetwork=true, kwargs...)
end

"warning: this model is not realistic or physically reasonable, it is only for test coverage"
function _build_mn_mc_opf_strg(pm::_PM.AbstractPowerModel)

    for (n, network) in nws(pm)
        variable_mc_bus_voltage(pm, nw=n)
        variable_mc_branch_power(pm, nw=n)
        variable_mc_transformer_power(pm, nw=n)
        variable_mc_gen_power_setpoint(pm, nw=n)
        variable_mc_load_setpoint(pm, nw=n)

        variable_mc_storage_power(pm, nw=n)

        constraint_mc_model_voltage(pm, nw=n)

        for i in ids(pm, :ref_buses, nw=n)
            constraint_mc_theta_ref(pm, i, nw=n)
        end

        # generators should be constrained before KCL, or Pd/Qd undefined
        for id in ids(pm, :gen)
            constraint_mc_gen_setpoint(pm, id, nw=n)
        end

        # loads should be constrained before KCL, or Pd/Qd undefined
        for id in ids(pm, :load)
            constraint_mc_load_setpoint(pm, id, nw=n)
        end

        for i in ids(pm, :bus, nw=n)
            constraint_mc_load_power_balance(pm, i, nw=n)

        end

        for i in ids(pm, :branch, nw=n)
            constraint_mc_ohms_yt_from(pm, i, nw=n)
            constraint_mc_ohms_yt_to(pm, i, nw=n)

            constraint_mc_voltage_angle_difference(pm, i, nw=n)

            constraint_mc_thermal_limit_from(pm, i, nw=n)
            constraint_mc_thermal_limit_to(pm, i, nw=n)
        end


        for i in ids(pm, :transformer, nw=n)
            constraint_mc_transformer_power(pm, i, nw=n)
        end

        for i in ids(pm, :storage, nw=n)
            _PM.constraint_storage_complementarity_nl(pm, i; nw=n)
            constraint_mc_storage_losses(pm, i; nw=n)
            constraint_mc_storage_thermal_limit(pm, i, nw=n)
        end

        # for i in ids(pm, :dcline, nw=n)
        #     constraint_mc_dcline(pm, i, nw=n)
        # end

    end
    network_ids = sort(collect(nw_ids(pm)))

    n_1 = network_ids[1]
    for i in ids(pm, :storage, nw=n_1)
        _PM.constraint_storage_state(pm, i, nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in ids(pm, :storage, nw=n_2)
            _PM.constraint_storage_state(pm, i, n_1, n_2)
        end
        n_1 = n_2
    end
    _PM.objective_min_fuel_cost(pm)
end
