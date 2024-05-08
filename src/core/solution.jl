"""
Definition of the default solution preprocessor for PowerModelsDistribution
"""
function _IM.solution_preprocessor(pm::AbstractUnbalancedPowerModel, solution::Dict)
    per_unit = _IM.get_data(x -> x["per_unit"], pm.data, pmd_it_name; apply_to_subnetworks=true)

    for (nw_id, nw_ref) in nws(pm)
        solution["it"][pmd_it_name]["nw"]["$(nw_id)"]["settings"] = nw_ref[:settings]
        solution["it"][pmd_it_name]["nw"]["$(nw_id)"]["per_unit"] = ismultinetwork(pm) ? per_unit["$(nw_id)"] : per_unit
    end
end


"custom `build_solution_values` for multiconductor (vector) JuMP expressions"
function _IM.build_solution_values(var::JuMP.AbstractJuMPScalar)
    return JuMP.value(var)
end


"custom `build_solution_values` for multiconductor (vector) constants"
function _IM.build_solution_values(var::JuMP.Containers.DenseAxisArray{<:Number,1})
    return var.data
end


"custom `build_solution_values` for generic dense axis arrays"
function _IM.build_solution_values(var::JuMP.Containers.DenseAxisArray{<:Any,1})
    return [_IM.build_solution_values(x) for x in var.data]
end


"custom `build_solution_values` for multiconductor (vector) constants"
function _IM.build_solution_values(var::LinearAlgebra.Symmetric{JuMP.VariableRef, Matrix{JuMP.VariableRef}})
    return JuMP.value.(var.data)
end


"converts w models voltages to standard voltage magnitude (sqrt)"
function _sol_data_model_w!(solution::Dict{String,<:Any})
    if haskey(solution, "nw")
        nws_data = solution["nw"]
    else
        nws_data = Dict("0" => solution)
    end

    for (n, nw_data) in nws_data
        if haskey(nw_data, "bus")
            for (i,bus) in nw_data["bus"]
                if haskey(bus, "w")
                    if any(bus["w"] .< 0) # e.g., as allowed by constraint violation settings
                        bus["vm"] = zeros(length(bus["w"]))
                        bus["vm"][bus["w"] .>= 0.0] .= sqrt.(bus["w"][bus["w"] .>= 0.0])
                    else
                        bus["vm"] = sqrt.(bus["w"])
                    end
                    delete!(bus, "w")
                end
                if haskey(bus, "Wr")
                    w = LinearAlgebra.diag(bus["Wr"])
                    if any(w .< 0) # e.g., as allowed by constraint violation settings
                        bus["vm"] = zeros(length(w))
                        bus["vm"][w .>= 0.0] .= sqrt.(w[w .>= 0.0])
                    else
                        bus["vm"] = sqrt.(w)
                        if length(w) == 3
                            t = [-1 1 0; -1 0 1; 0 -1 1]
                            va = LinearAlgebra.pinv(t)*[atan(bus["Wi"][2,1], bus["Wr"][2,1]);
                                                        atan(bus["Wi"][3,1], bus["Wr"][3,1]);
                                                        atan(bus["Wi"][3,2], bus["Wr"][3,2])]
                            bus["va"] = [va[findmin(abs.(va .- 0))[2]],
                                         va[findmin(abs.(va .+ 2*pi/3))[2]],
                                         va[findmin(abs.(va .- 2*pi/3))[2]]] # TODO: better way to get angles in order
                        end
                    end
                    delete!(bus, "Wr")
                    delete!(bus, "Wi")
                end
            end
        end
    end
end


"""
    sol_data_model!(pm::AbstractUnbalancedWModels, solution::Dict{String,<:Any})

solution_processor, see [`solve_mc_model`](@ref solve_mc_model), to convert W variables
back into polar representation (default data model voltage form)
"""
function sol_data_model!(pm::AbstractUnbalancedWModels, solution::Dict{String,<:Any})
    apply_pmd!(_sol_data_model_w!, solution)
end


"""
    _sol_data_model_acr!(solution::Dict{String,<:Any})

solution_processor, see [`solve_mc_model`](@ref solve_mc_model),  to convert ACR variables
back into polar representation (default data model voltage form)
"""
function _sol_data_model_acr!(solution::Dict{String,<:Any})
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


"""
    sol_data_model!(pm::AbstractUnbalancedACRModel, solution::Dict{String,<:Any})

solution_processor, see [`solve_mc_model`](@ref solve_mc_model),  to convert ACR variables
back into polar representation (default data model voltage form)
"""
function sol_data_model!(pm::AbstractUnbalancedACRModel, solution::Dict{String,<:Any})
    apply_pmd!(_sol_data_model_acr!, solution)
end


"""
    sol_data_model!(pm::FBSUBFPowerModel, solution::Dict{String,<:Any})

solution_processor, to convert FBS variables
back into polar representation (default data model voltage form)
"""
function sol_data_model!(pm::FBSUBFPowerModel, solution::Dict{String,<:Any})
    apply_pmd!(_sol_data_model_acr!, solution)
end


"""
    sol_data_model!(pm::FOTRUPowerModel, solution::Dict{String,<:Any})

solution_processor, to convert FOT rectangular variables
back into polar representation (default data model voltage form)
"""
function sol_data_model!(pm::FOTRUPowerModel, solution::Dict{String,<:Any})
    apply_pmd!(_sol_data_model_acr!, solution)
end


"""
    sol_data_model!(pm::AbstractUnbalancedPowerModel, solution::Dict{String,<:Any})

does nothing (no `sol_data_model!` exists for the formulation attempting to be converted)
"""
function sol_data_model!(pm::AbstractUnbalancedPowerModel, solution::Dict{String,<:Any})
    @info "sol_data_model! not defined for power model of type $(typeof(pm))"
end
