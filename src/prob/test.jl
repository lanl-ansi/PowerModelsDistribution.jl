"test mc mn"
function _run_mc_mn_opb(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return run_mc_model(data, model_type, solver, _build_mc_mn_opb; ref_extensions=[_PM.ref_add_connected_components!], multinetwork=true, kwargs...)
end


"Constructor for Optimal Power Flow"
function _build_mc_mn_opb(pm::_PM.AbstractPowerModel)
    for (n, network) in nws(pm)
        variable_mc_gen_power_setpoint(pm; nw=n)

        for i in ids(pm, :components, nw=n)
            constraint_mc_network_power_balance(pm, i; nw=n)
        end
    end

    _PM.objective_min_fuel_cost(pm)
end
