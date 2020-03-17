"alias to run_model in PMs with multiconductor=true, and transformer ref extensions added by default"
function run_mc_model(data::Dict{String,<:Any}, model_type, solver, build_mc; ref_extensions=[], kwargs...)
    _PMs.run_model(data, model_type, solver, build_mc; ref_extensions=[ref_add_arcs_trans!, ref_extensions...], multiconductor=true, kwargs...)
end
