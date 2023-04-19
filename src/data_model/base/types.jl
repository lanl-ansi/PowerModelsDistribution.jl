# TODO: Move to InfrastructureModels
abstract type MultinetworkModel end
abstract type SubnetworkModel end
abstract type InfrastructureDataModel{T} end
abstract type GenericInfrastructureObject end

# PowerModelsDistribution
abstract type DistributionDataModel{T} <: InfrastructureDataModel{T} end

# Core Data Model Types
abstract type DssRawModel <: DistributionDataModel{SubnetworkModel} end
abstract type DssModel <: DistributionDataModel{SubnetworkModel} end

abstract type EngineeringModel{T} <: DistributionDataModel{T} end
abstract type MathematicalModel{T} <: DistributionDataModel{T} end

# OpenDss Model Abstract Object Types
abstract type DssObject <: GenericInfrastructureObject end
abstract type DssNodeObject <: DssObject end
abstract type DssEdgeObject <: DssObject end
abstract type DssDataObject <: DssObject end
abstract type DssControlObject <: DssObject end

# Engineering Model Abstract Object Types
abstract type EngObject <: GenericInfrastructureObject end
abstract type EngNodeObject <: EngObject end
abstract type EngEdgeObject <: EngObject end
abstract type EngDataObject <: EngObject end
abstract type EngControlObject <: EngObject end

# Mathematical Model Abstract Object Types
abstract type MathObject <: GenericInfrastructureObject end
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
