"Solve multinetwork network expansion problem"
function solve_mn_mc_mld_simple_ne(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mn_mc_mld_simple_ne; multinetwork=true, kwargs...)
end


"Multinetwork network expansion problem including storage"
function build_mn_mc_mld_simple_ne(pm::AbstractUnbalancedPowerModel)
    build_mn_mc_mld_simple(pm)
    objective_ne(pm)
end

"Multinetwork network expansion problem including storage"
function build_mn_mc_mld_multi_scenario_ne(pm::AbstractUnbalancedPowerModel)
    build_mn_mc_mld_multi_scenario(pm)
    objective_ne(pm)
end