"multiconductor version of instantiate model from PowerModels"
function instantiate_mc_model(data::Dict{String,<:Any}, model_type::Type, build_method::Function; ref_extensions::Vector{<:Function}=Vector{Function}([]), kwargs...)
    if get(data, "data_model", MATHEMATICAL) == ENGINEERING
        Memento.info(_LOGGER, "Converting ENGINEERING data model to MATHEMATICAL first to build JuMP model")
        data = transform_data_model(data)
    end

    return _PM.instantiate_model(data, model_type, build_method; ref_extensions=[ref_extensions..., ref_add_arcs_transformer!, ref_add_connections!], kwargs...)
end
