# core InfrastructureModels types
abstract type MultinetworkModel end
abstract type NetworkModel end
abstract type InfrastructureModel{T} end
abstract type InfrastructureObject end

# PowerModelsDistribution
abstract type DistributionModel{T} <: InfrastructureModel{T} end

# Core Data Model Types
abstract type DssRawModel <: DistributionModel{NetworkModel} end
abstract type DssModel <: DistributionModel{NetworkModel} end

abstract type EngineeringModel{T} <: DistributionModel{T} end
abstract type MathematicalModel{T} <: DistributionModel{T} end

# OpenDss Model Abstract Object Types
abstract type DssObject <: InfrastructureObject end
abstract type DssNodeObject <: DssObject end
abstract type DssEdgeObject <: DssObject end
abstract type DssDataObject <: DssObject end
abstract type DssControlObject <: DssObject end

# Engineering Model Abstract Object Types
abstract type EngObject <: InfrastructureObject end
abstract type EngNodeObject <: EngObject end
abstract type EngEdgeObject <: EngObject end
abstract type EngDataObject <: EngObject end
abstract type EngControlObject <: EngObject end

# Mathematical Model Abstract Object Types
abstract type MathObject <: InfrastructureObject end
abstract type MathNodeObject <: MathObject end
abstract type MathEdgeObject <: MathObject end
abstract type MathDataObject <: MathObject end
abstract type MathControlObject <: MathObject end


"""
"""
macro extend(base, name, fields)
    base_type = Core.eval(@__MODULE__, base)

    base_fieldnames = fieldnames(base_type)
    base_types = [t for t in base_type.types]

    base_fields = [:($f::$T) for (f,T) in zip(base_fieldnames, base_types)]

    res = :(mutable struct $name <: $(supertype(base_type)) end)

    push!(res.args[end].args, base_fields...)
    push!(res.args[end].args, fields.args...)

    return res
end

"default `transform_data_model` ErrorException for unsupported combinations"
transform_data_model(type::Type{T}, data::U; kwargs...) where {V, T <: InfrastructureModel{V}, U <: DistributionModel{V}} = throw(ErrorException("Transforming from $(typeof(data)) to $type is not supported"))
