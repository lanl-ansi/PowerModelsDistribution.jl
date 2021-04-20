""
function _solve_mc_model(data::Dict{String,<:Any}, model_type::Type, optimizer, build_method::Function;
        ref_extensions::Vector{<:Function}=Function[], solution_processors::Vector{<:Function}=Function[], relax_integrality::Bool=false,
        multinetwork::Bool=false, kwargs...)

    if multinetwork != ismultinetwork(data)
        model_requirement = multinetwork ? "multi-network" : "single-network"
        data_type = ismultinetwork(data) ? "multi-network" : "single-network"
        error("attempted to build a $(model_requirement) model with $(data_type) data")
    end

    start_time = time()
    pm = instantiate_mc_model(data, model_type, build_method; ref_extensions=ref_extensions, kwargs...)
    @debug "pm model build time: $(time() - start_time)"

    start_time = time()
    result = optimize_model!(pm, relax_integrality=relax_integrality, optimizer=optimizer, solution_processors=solution_processors)
    @debug "pm model solve and solution time: $(time() - start_time)"

    return result
end


""
function instantiate_mc_model(data::Dict{String,<:Any}, model_type::Type, build_method::Function; ref_extensions::Vector{<:Function}=Function[], multinetwork::Bool=false, kwargs...)
    if get(data, "data_model", MATHEMATICAL) == ENGINEERING
        @info "Converting ENGINEERING data model to MATHEMATICAL first to build JuMP model"
        data = transform_data_model(data; multinetwork=multinetwork)
    end

    return _IM.instantiate_model(
        data, model_type, build_method, ref_add_core!, _pmd_global_keys,
        pmd_it_sym; ref_extensions = ref_extensions, kwargs...)
end


""
function solve_mc_model(data::Dict{String,<:Any}, model_type::Type, solver, build_mc::Function; ref_extensions::Vector{<:Function}=Function[], make_si::Bool=!get(data, "per_unit", false), multinetwork::Bool=false, kwargs...)::Dict{String,Any}
    if get(data, "data_model", MATHEMATICAL) == ENGINEERING
        data_math = transform_data_model(data; multinetwork=multinetwork)

        result = _solve_mc_model(data_math, model_type, solver, build_mc; ref_extensions=ref_extensions, multinetwork=multinetwork, kwargs...)

        result["solution"] = transform_solution(result["solution"], data_math; make_si=make_si)
    elseif get(data, "data_model", MATHEMATICAL) == MATHEMATICAL
        result = _solve_mc_model(data, model_type, solver, build_mc; ref_extensions=ref_extensions, multinetwork=multinetwork, kwargs...)
    end

    return result
end


""
function solve_mc_model(file::String, model_type::Type, solver, build_mc::Function; ref_extensions::Vector{<:Function}=Function[], kwargs...)::Dict{String,Any}
    return solve_mc_model(parse_file(file), model_type, solver, build_mc; ref_extensions=ref_extensions, kwargs...)
end


"depreciation message for run_mc_model"
function run_mc_model(data::Union{String,Dict{String,<:Any}}, model_type::Type, solver, build_mc::Function; kwargs...)::Dict{String,Any}
    @warn "run_mc_model is being depreciated in favor of solve_mc_model, please update your code accordingly"
    return solve_mc_model(data, model_type, solver, build_mc; kwargs...)
end
