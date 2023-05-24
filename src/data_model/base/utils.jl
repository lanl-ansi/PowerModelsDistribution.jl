"Helper function to convert InfrastructureModel,InfrastructureObject into Dict{String,Any}"
function _convert_model_to_dict(data::Union{InfrastructureModel,InfrastructureObject})::Dict{String,Any}
    out = Dict{String,Any}(
    )

    for property in propertynames(data)
        item = getproperty(data, property)

        if isa(item, Dict)
            out["$property"] = Dict{String,Any}()
            for (id, obj) in item
                if isa(obj, InfrastructureObject)
                    out["$property"]["$id"] = _convert_model_to_dict(obj)
                else
                    out["$property"]["$id"] = obj
                end
            end
        elseif isa(item, InfrastructureObject)
            out["$property"] = _convert_model_to_dict(item)
        else
            out["$property"] = item
        end
    end

    return filter(x->!ismissing(x.second), out)
end


"Helper function to pull the specified properties from dss property pairs"
function _get_raw_fields(property_pairs::Vector{Pair{String,String}})::Vector{Symbol}
    collect(Symbol(x.first) for x in property_pairs)
end
