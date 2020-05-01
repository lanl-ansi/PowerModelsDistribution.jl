"alias to run_model in PMs with multiconductor=true, and transformer ref extensions added by default"
function run_mc_model(data::Dict{String,<:Any}, model_type::DataType, solver, build_mc::Function; ref_extensions::Vector{<:Function}=Vector{Function}([]), make_si=!get(data, "per_unit", false), kwargs...)::Dict{String,Any}

    if get(data, "data_model", MATHEMATICAL) == ENGINEERING
        data_math = transform_data_model(data)

        result = _PM.run_model(data_math, model_type, solver, build_mc; ref_extensions=[ref_add_arcs_transformer!, ref_extensions...], multiconductor=true, kwargs...)

        result["solution"] = transform_solution(result["solution"], data_math; make_si=make_si)
    elseif get(data, "data_model", MATHEMATICAL) == MATHEMATICAL
        result = _PM.run_model(data, model_type, solver, build_mc; ref_extensions=[ref_add_arcs_transformer!, ref_extensions...], multiconductor=true, kwargs...)
    end

    return result
end


"alias to run_model in PMs with multiconductor=true, and transformer ref extensions added by default"
function run_mc_model(file::String, model_type::DataType, solver, build_mc::Function; ref_extensions::Vector{<:Function}=Vector{Function}([]), kwargs...)::Dict{String,Any}
    return run_mc_model(parse_file(file), model_type, solver, build_mc; ref_extensions=ref_extensions, kwargs...)
end
