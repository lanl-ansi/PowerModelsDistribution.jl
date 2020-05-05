# This file contains useful transformation functions for the engineering data model

const _loss_model_objects = Dict{String,Vector{String}}(
    "switch" => Vector{String}(["linecode", "rs", "xs", "g_fr", "b_fr", "g_to", "b_to"]),
    "voltage_source" => Vector{String}(["rs", "xs"]),
    "transformer" => Vector{String}(["rs", "xsc", "imag", "noloadloss"])
)


"remove parameters from objects with loss models to make them lossless"
function make_lossless!(data_eng::Dict{String,<:Any})
    for (object_type, parameters) in _loss_model_objects
        if haskey(data_eng, object_type)
            for (id, eng_obj) in data_eng[object_type]
                for parameter in parameters
                    if haskey(eng_obj, parameter)
                        if object_type == "transformer"
                            eng_obj[parameter] = 0 .* eng_obj[parameter]
                        else
                            delete!(eng_obj, parameter)
                        end
                    end
                end
            end
        end
    end
end
