######
#
# These are toy problem formulations used to test advanced features
# such as storage devices
#
######
# "multi-network opf with storage"
# function _run_mn_mc_opf(data::Dict{String,Any}, model_type, solver; kwargs...)
#     return _PMs.run_model(data, model_type, solver, _build_mn_mc_strg_opf; ref_extensions=[ref_add_arcs_trans!], multiconductor=true, multinetwork=true, kwargs...)
# end
#
#
# "multi-network opf with storage"
# function _run_mn_mc_opf(file::String, model_type, solver; kwargs...)
#     return run_mn_mc_opf(PowerModelsDistribution.parse_file(file), model_type, solver; kwargs...)
# end
#
#
# "multi-network opf with storage"
# function _build_mn_mc_strg_opf(pm::_PMs.AbstractPowerModel)
#     for (n, network) in _PMs.nws(pm)
#         variable_mc_voltage(pm; nw=n)
#         constraint_mc_model_voltage(pm; nw=n)
#         variable_mc_branch_flow(pm; nw=n)
#         variable_mc_generation(pm; nw=n)
#         variable_mc_storage(pm; nw=n)
#
#         for i in _PMs.ids(pm, :ref_buses; nw=n)
#             constraint_mc_theta_ref(pm, i; nw=n)
#         end
#
#         for i in _PMs.ids(pm, :bus; nw=n)
#             constraint_mc_power_balance(pm, i; nw=n)
#         end
#
#         for i in _PMs.ids(pm, :storage; nw=n)
#             _PMs.constraint_storage_state(pm, i; nw=n)
#             _PMs.constraint_storage_complementarity_nl(pm, i; nw=n)
#             constraint_mc_storage_loss(pm, i; nw=n)
#             constraint_mc_storage_thermal_limit(pm, i; nw=n)
#         end
#
#         for i in _PMs.ids(pm, :branch; nw=n)
#             constraint_mc_ohms_yt_from(pm, i; nw=n)
#             constraint_mc_ohms_yt_to(pm, i; nw=n)
#
#             constraint_mc_voltage_angle_difference(pm, i; nw=n)
#
#             constraint_mc_thermal_limit_from(pm, i; nw=n)
#             constraint_mc_thermal_limit_to(pm, i; nw=n)
#         end
#
#         for i in _PMs.ids(pm, :transformer; nw=n)
#             constraint_mc_trans(pm, i; nw=n)
#         end
#     end
#
#     network_ids = sort(collect(_PMs.nw_ids(pm)))
#
#     n_1 = network_ids[1]
#     for i in _PMs.ids(pm, :storage; nw=n_1)
#         _PMs.constraint_storage_state(pm, i; nw=n_1)
#     end
#
#     for n_2 in network_ids[2:end]
#         for i in _PMs.ids(pm, :storage; nw=n_2)
#             _PMs.constraint_storage_state(pm, i, n_1, n_2)
#         end
#         n_1 = n_2
#     end
#
#     _PMs.objective_min_fuel_cost(pm)
# end


######
#
# Formulations from PowerModels
#
######

""
function _run_mc_ucopf(file, model_type::Type, solver; kwargs...)
    return _PMs.run_model(file, model_type, solver, _build_mc_ucopf; ref_extensions=[ref_add_arcs_trans!], multiconductor=true, kwargs...)
end

""
function _build_mc_ucopf(pm::_PMs.AbstractPowerModel)
    for (n, network) in _PMs.nws(pm)
        variable_mc_voltage(pm, nw=n)
        variable_mc_branch_flow(pm, nw=n)
        variable_mc_transformer_flow(pm, nw=n)

        variable_mc_indicator_generation(pm, nw=n)
        variable_mc_generation_on_off(pm, nw=n)


        constraint_mc_model_voltage(pm, nw=n)



        variable_mc_on_off_storage(pm, nw=n)
        _PMs.variable_storage_energy(pm, nw=n)
        _PMs.variable_storage_charge(pm, nw=n)
        _PMs.variable_storage_discharge(pm, nw=n)
        _PMs.variable_storage_indicator(pm, nw=n)
        _PMs.variable_storage_complementary_indicator(pm, nw=n)


        for i in _PMs.ids(pm, :ref_buses, nw=n)
            constraint_mc_theta_ref(pm, i, nw=n)
        end

        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_mc_power_balance(pm, i, nw=n)

        end

        for i in _PMs.ids(pm, :branch, nw=n)
            constraint_mc_ohms_yt_from(pm, i, nw=n)
            constraint_mc_ohms_yt_to(pm, i, nw=n)

            constraint_mc_voltage_angle_difference(pm, i, nw=n)

            constraint_mc_thermal_limit_from(pm, i, nw=n)
            constraint_mc_thermal_limit_to(pm, i, nw=n)
        end


        for i in _PMs.ids(pm, :transformer, nw=n)
            constraint_mc_trans(pm, i, nw=n)
        end

        # for i in ids(pm, :dcline, nw=n)
        #     constraint_mc_dcline(pm, i, nw=n)
        # end

        for i in _PMs.ids(pm, :storage; nw=n)
            # _PMs.constraint_storage_state(pm, i; nw=n)
            _PMs.constraint_storage_complementarity_mi(pm, i; nw=n)
            constraint_mc_storage_loss(pm, i; nw=n)
            constraint_mc_storage_thermal_limit(pm, i; nw=n)

            constraint_mc_storage_on_off(pm, i; nw=n)

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


""
function _run_mn_mc_opf(file, model_type::Type, optimizer; kwargs...)
    return _PMs.run_model(file, model_type, optimizer, _build_mn_mc_opf; ref_extensions=[ref_add_arcs_trans!], multinetwork=true, multiconductor=true, kwargs...)
end

""
function _build_mn_mc_opf(pm::_PMs.AbstractPowerModel)
    for (n, network) in _PMs.nws(pm)
        variable_mc_voltage(pm, nw=n)
        variable_mc_branch_flow(pm, nw=n)
        variable_mc_transformer_flow(pm, nw=n)
        variable_mc_generation(pm, nw=n)

        constraint_mc_model_voltage(pm, nw=n)

        for i in _PMs.ids(pm, :ref_buses, nw=n)
            constraint_mc_theta_ref(pm, i, nw=n)
        end

        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_mc_power_balance(pm, i, nw=n)

        end

        for i in _PMs.ids(pm, :branch, nw=n)
            constraint_mc_ohms_yt_from(pm, i, nw=n)
            constraint_mc_ohms_yt_to(pm, i, nw=n)

            constraint_mc_voltage_angle_difference(pm, i, nw=n)

            constraint_mc_thermal_limit_from(pm, i, nw=n)
            constraint_mc_thermal_limit_to(pm, i, nw=n)
        end


        for i in _PMs.ids(pm, :transformer, nw=n)
            constraint_mc_trans(pm, i, nw=n)
        end

        # for i in ids(pm, :dcline, nw=n)
        #     constraint_mc_dcline(pm, i, nw=n)
        # end
    end
    _PMs.objective_min_fuel_cost(pm)
end


""
function _run_mn_mc_opf_strg(file, model_type::Type, optimizer; kwargs...)
    return _PMs.run_model(file, model_type, optimizer, _build_mn_mc_opf_strg; ref_extensions=[ref_add_arcs_trans!], multinetwork=true, multiconductor=true, kwargs...)
end

"warning: this model is not realistic or physically reasonable, it is only for test coverage"
function _build_mn_mc_opf_strg(pm::_PMs.AbstractPowerModel)

    for (n, network) in _PMs.nws(pm)
        variable_mc_voltage(pm, nw=n)
        variable_mc_branch_flow(pm, nw=n)
        variable_mc_transformer_flow(pm, nw=n)
        variable_mc_generation(pm, nw=n)

        variable_mc_storage(pm, nw=n)

        constraint_mc_model_voltage(pm, nw=n)

        for i in _PMs.ids(pm, :ref_buses, nw=n)
            constraint_mc_theta_ref(pm, i, nw=n)
        end

        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_mc_power_balance(pm, i, nw=n)

        end

        for i in _PMs.ids(pm, :branch, nw=n)
            constraint_mc_ohms_yt_from(pm, i, nw=n)
            constraint_mc_ohms_yt_to(pm, i, nw=n)

            constraint_mc_voltage_angle_difference(pm, i, nw=n)

            constraint_mc_thermal_limit_from(pm, i, nw=n)
            constraint_mc_thermal_limit_to(pm, i, nw=n)
        end


        for i in _PMs.ids(pm, :transformer, nw=n)
            constraint_mc_trans(pm, i, nw=n)
        end

        for i in _PMs.ids(pm, :storage, nw=n)
            _PMs.constraint_storage_complementarity_nl(pm, i; nw=n)
            constraint_mc_storage_loss(pm, i; nw=n)
            constraint_mc_storage_thermal_limit(pm, i, nw=n)
        end

        # for i in ids(pm, :dcline, nw=n)
        #     constraint_mc_dcline(pm, i, nw=n)
        # end

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
