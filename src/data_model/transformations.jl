# This file contains useful transformation functions for the engineering data model

const _loss_model_objects = Dict{String,Vector{String}}(
    "switch" => Vector{String}(["linecode", "rs", "xs"]),
    "voltage_source" => Vector{String}(["rs", "xs"]),
    "transformer" => Vector{String}(["rw", "xsc", "imag", "noloadloss"])
)


"remove parameters from objects with loss models to make them lossless"
function make_lossless!(data_eng::Dict{String,<:Any})
    @assert data_eng["data_model"] == ENGINEERING "incorrect data model type"

    for (object_type, parameters) in _loss_model_objects
        if haskey(data_eng, object_type)
            for (id, eng_obj) in data_eng[object_type]
                for parameter in parameters
                    if haskey(eng_obj, parameter)
                        if parameter == "linecode"
                            delete!(eng_obj, parameter)
                        else
                            eng_obj[parameter] = 0 .* eng_obj[parameter]
                        end
                    end
                end
            end
        end
    end
end


"add voltage bounds"
function apply_voltage_bounds!(data_eng::Dict{String,<:Any}; vm_lb::Union{Real,Missing}=0.9, vm_ub::Union{Real,Missing}=1.1)
    @assert data_eng["data_model"] == ENGINEERING "incorrect data model type"

    (bus_vbases, edge_vbases) = calc_voltage_bases(data_eng, data_eng["settings"]["vbases_default"])
    if haskey(data_eng, "bus")
        for (id,bus) in data_eng["bus"]
            vbase = bus_vbases[id]
            if !ismissing(vm_lb) && !haskey(bus, "vm_lb")
                bus["vm_lb"] = vbase .* fill(vm_lb, length(bus["terminals"]))
            end

            if !ismissing(vm_ub) && !haskey(bus, "vm_ub")
                bus["vm_ub"] = vbase .* fill(vm_ub, length(bus["terminals"]))
            end
        end
    end
end
