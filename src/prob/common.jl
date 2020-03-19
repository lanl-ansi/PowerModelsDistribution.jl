"alias to run_model in PMs with multiconductor=true, and transformer ref extensions added by default"
function run_mc_model(data::Dict{String,<:Any}, model_type, solver, build_mc; ref_extensions=[], kwargs...)
    if get(data, "data_model", "mathematical") == "engineering"
        data_math = transform_data_model(data)

        result = _PMs.run_model(data_math, model_type, solver, build_mc; ref_extensions=[ref_add_arcs_trans!, ref_extensions...], multiconductor=true, kwargs...)

        result["solution"] = solution_math2eng(result["solution"], data_math; make_si=get(data, "per_unit", false), make_deg=get(data, "per_unit", false))
    else
        result = _PMs.run_model(data, model_type, solver, build_mc; ref_extensions=[ref_add_arcs_trans!, ref_extensions...], multiconductor=true, kwargs...)
    end

    return result
end
