"solve test mn mc problem"
function _solve_mn_mc_opb(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, _build_mc_mn_opb; ref_extensions=[ref_add_connected_components!], multinetwork=true, kwargs...)
end


"Constructor for Optimal Power Balance"
function _build_mc_mn_opb(pm::AbstractUnbalancedPowerModel)
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


""
function ref_add_connected_components!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    apply_pmd!(_ref_add_connected_components!, ref, data)
end


""
function _ref_add_connected_components!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    component_sets = calc_connected_components(data)
    ref[:components] = Dict(i => c for (i,c) in enumerate(sort(collect(component_sets); by = length)))
end
