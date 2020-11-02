function _IM.build_solution_values(var::JuMP.Containers.DenseAxisArray{<:JuMP.VariableRef,1})
    return JuMP.value.(var.data)
end


function _IM.build_solution_values(var::JuMP.Containers.DenseAxisArray{<:JuMP.NonlinearExpression,1})
    return JuMP.value.(var.data)
end


function _IM.build_solution_values(var::JuMP.Containers.DenseAxisArray{<:JuMP.GenericAffExpr,1})
    return JuMP.value.(var.data)
end


function _IM.build_solution_values(var::JuMP.Containers.DenseAxisArray{<:Number,1})
    return var.data
end


""
function sol_data_model!(pm::_PM.AbstractWModels, solution::Dict{String,<:Any})
    if haskey(solution, "nw")
        nws_data = solution["nw"]
    else
        nws_data = Dict("0" => solution)
    end

    for (n, nw_data) in nws_data
        if haskey(nw_data, "bus")
            for (i,bus) in nw_data["bus"]
                if haskey(bus, "w")
                    bus["vm"] = sqrt.(bus["w"])
                    delete!(bus, "w")
                end
            end
        end
    end
end


""
function sol_data_model!(pm::_PM.AbstractACRModel, solution::Dict{String,<:Any})
    if haskey(solution, "nw")
        nws_data = solution["nw"]
    else
        nws_data = Dict("0" => solution)
    end

    for (n, nw_data) in nws_data
        if haskey(nw_data, "bus")
            for (i,bus) in nw_data["bus"]
                if haskey(bus, "vr") && haskey(bus, "vi")
                    bus["vm"] = sqrt.(bus["vr"].^2 + bus["vi"].^2)
                    bus["va"] = atan.(bus["vi"], bus["vr"])

                    delete!(bus, "vr")
                    delete!(bus, "vi")
                end
            end
        end
    end
end


"converts the solution data into the data model's standard space, polar voltages and rectangular power"
function sol_data_model!(pm::_PM.AbstractPowerModel, solution::Dict{String,<:Any})
    Memento.warn(_LOGGER, "sol_data_model! not defined for power model of type $(typeof(pm))")
end
