export run_generic_tp_model

""
function run_generic_tp_model(file::String, model_constructor, solver, post_method; solution_builder=PMs.get_solution, kwargs...)
    data = ThreePhasePowerModels.parse_file(file)
    return PMs.run_generic_model(data, model_constructor, solver, post_method; multiphase=true, solution_builder=solution_builder, kwargs...)
end

""
function run_generic_tp_model(data::Dict{String,Any}, model_constructor, solver, post_method; solution_builder=PMs.get_solution, kwargs...)
    return PMs.run_generic_model(data, model_constructor, solver, post_method; multiphase=true, solution_builder=solution_builder, kwargs...)
end