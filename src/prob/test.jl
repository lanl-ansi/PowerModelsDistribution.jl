"solve test mn mc problem"
function _solve_mn_mc_opb(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, _build_mc_mn_opb; ref_extensions=[_PM.ref_add_connected_components!], multinetwork=true, kwargs...)
end


"Constructor for Optimal Power Balance"
function _build_mc_mn_opb(pm::_PM.AbstractPowerModel)
    for (n, network) in nws(pm)
        variable_mc_generator_power(pm; nw=n)

        for i in ids(pm, :components, nw=n)
            constraint_mc_network_power_balance(pm, i; nw=n)
        end
    end

    objective_mc_min_fuel_cost(pm)
end

# Depreciated run_ functions (remove in ~4-6 months)

"depreciation warning for _run_mc_mn_opb"
function _run_mc_mn_opb(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    @warn "_run_mc_mn_opb is being depreciated in favor of _solve_mn_mc_opb, please update your code accordingly"
    return _solve_mn_mc_opb(data, model_type, solver; kwargs...)
end
