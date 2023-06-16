# core InfrastructureModels types
"Multinetwork subtype"
abstract type MultinetworkModel end

"Single network (non-multinetwork) subtype"
abstract type NetworkModel end

"Base Generic Infrastructure Model type"
abstract type InfrastructureModel{T} end

"Base Generic Infrastructur Object type"
abstract type InfrastructureObject end

# PowerModelsDistribution
"Base Model type for power distribution networks"
abstract type DistributionModel{T} <: InfrastructureModel{T} end

# Core Data Model Types
"Base Model type for the raw dss model (property pairs)"
abstract type DssRawModel <: DistributionModel{NetworkModel} end

"Base Model type of the parsed dss model (structs)"
abstract type DssModel <: DistributionModel{NetworkModel} end

# OpenDss Model Abstract Object Types
"Generic dss object"
abstract type DssObject <: InfrastructureObject end

"Generic dss node-type object"
abstract type DssNodeObject <: DssObject end

"Generic dss edge-type object"
abstract type DssEdgeObject <: DssObject end

"Generic dss data-type object"
abstract type DssDataObject <: DssObject end

"Generic dss control-type object"
abstract type DssControlObject <: DssObject end

"default `transform_data_model` ErrorException for unsupported combinations"
transform_data_model(type::Type{T}, data::U; kwargs...) where {V, T <: InfrastructureModel{V}, U <: DistributionModel{V}} = throw(ErrorException("Transforming from $(typeof(data)) to $type is not supported"))
