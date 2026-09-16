# core InfrastructureModels types
"Multinetwork subtype"
abstract type MultinetworkModel end

"Single network (non-multinetwork) subtype"
abstract type NetworkModel end

"Base Generic Infrastructure Model type"
abstract type InfrastructureModel{T} end

"Base Generic Infrastructure Object type"
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

# EngineeringModel Abstract Object Types
abstract type EngineeringObject <: InfrastructureObject end
abstract type EngineeringEdgeObject <: EngineeringObject end
abstract type EngineeringNodeObject <: EngineeringObject end
abstract type EngineeringDataObject <: EngineeringObject end

# MathematicalModel Abstract Object Types
abstract type MathObject <: InfrastructureObject end
abstract type MathEdgeObject <: MathObject end
abstract type MathNodeObject <: MathObject end

"default `transform_data_model` ErrorException for unsupported combinations"
function transform_data_model(type::Type{T}, data::U; kwargs...) where {V,T<:InfrastructureModel{V},U<:DistributionModel{V}}
    if type == U
        return data
    else
        throw(ErrorException("Transforming from $(typeof(data)) to $type is not supported"))
    end
end

get_model_type(::DistributionModel{T}) where {T} = T
get_model_type(::Type{<:DistributionModel{T}}) where {T} = T
get_model_type(::Type{<:DistributionModel}) = Nothing

mutable struct MathematicalModel{T<:Union{NetworkModel,MultinetworkModel}} <: DistributionModel{T}
    data::Dict{String,Any}

    function MathematicalModel(data::Union{Dict{String,Any},DistributionModel,Nothing}=nothing; multinetwork::Bool=false)
        if data isa Dict
            if get(data, "multinetwork", multinetwork)
                obj = new{MultinetworkModel}(data)
            else
                obj = new{NetworkModel}(data)
            end
        elseif data isa MathematicalModel
            if get_model_type(data) == NetworkModel && multinetwork
                obj = make_multinetwork(data)
            else
                obj = data
            end
        elseif data isa EngineeringModel
            obj = transform_data_model(MathematicalModel{get_model_type(data)}, data)
        elseif data isa DssModel
            if multinetwork
                obj = transform_data_model(MathematicalModel{MultinetworkModel}, data)
            else
                obj = transform_data_model(MathematicalModel{NetworkModel}, data)
            end
        else
            error("Unrecognized data format, cannot instantiate MathematicalModel")
        end

        if obj isa MathematicalModel{MultinetworkModel}
            nw_dict = get!(() -> Dict{String,Any}(), obj.data, "nw")

            for (k, v) in copy(nw_dict)
                if !(v isa MathematicalModel{NetworkModel})
                    nw_dict[k] = MathematicalModel(v)
                end
            end
        end

        return obj
    end
end

MathematicalModel{T}(data::Dict{String,Any}) where {T<:Union{NetworkModel,MultinetworkModel}} = MathematicalModel(data; multinetwork=T === MultinetworkModel)

mutable struct EngineeringModel{T<:Union{NetworkModel,MultinetworkModel}} <: DistributionModel{T}
    data::Dict{String,Any}

    function EngineeringModel(data::Union{Dict{String,<:Any},DistributionModel,Nothing}=nothing; multinetwork::Bool=false)
        if data isa Dict
            if get(data, "multinetwork", multinetwork)
                obj = new{MultinetworkModel}(data)
            else
                obj = new{NetworkModel}(data)
            end
        elseif data isa MathematicalModel
            error("Cannot construct an EngineeringModel from a MathematicalModel")
        elseif data isa EngineeringModel
            if get_model_type(data) == NetworkModel && multinetwork
                obj = make_multinetwork(data)
            else
                obj = data
            end
        elseif data isa DssModel
            if multinetwork
                obj = transform_data_model(EngineeringModel{MultinetworkModel}, data)
            else
                obj = transform_data_model(EngineeringModel{NetworkModel}, data)
            end
        else
            error("Unrecognized data format, cannot instantiate EngineeringModel")
        end

        if obj isa EngineeringModel{MultinetworkModel}
            nw_dict = get!(() -> Dict{String,Any}(), obj.data, "nw")

            for (k, v) in copy(nw_dict)
                if !(v isa EngineeringModel{NetworkModel})
                    nw_dict[k] = EngineeringModel(v)
                end
            end
        end

        return obj
    end
end

EngineeringModel{T}(data::Dict{String,Any}) where {T<:Union{NetworkModel,MultinetworkModel}} = EngineeringModel(data; multinetwork=T === MultinetworkModel)

function Base.getproperty(m::InfrastructureModel{MultinetworkModel}, name::Symbol)
    if name === :nw
        return get!(() -> Dict{String,MathematicalModel{NetworkModel}}(), m.data, "nw")
    else
        return getfield(m, name)
    end
end

function Base.setproperty!(m::InfrastructureModel{MultinetworkModel}, name::Symbol, value)
    if name === :nw
        m.data["nw"] = value
    else
        setfield!(m, name, value)
    end
end


mutable struct MathematicalSolution{T<:Union{NetworkModel,MultinetworkModel}} <: DistributionModel{T}
    data::Dict{String,<:Any}

    function MathematicalSolution(data::Union{Dict{String,<:Any},Nothing}=nothing; multinetwork::Bool=false)
        if data isa Dict
            if get(data, "multinetwork", multinetwork)
                obj = new{MultinetworkModel}(data)
            else
                obj = new{NetworkModel}(data)
            end
        elseif data isa MathematicalSolution
            obj = data
        else
            error("Unrecognized data format, cannot instantiate MathematicalSolution")
        end

        if obj isa MathematicalSolution{MultinetworkModel}
            nw_dict = get!(() -> Dict{String,Any}(), obj.data, "nw")

            for (k, v) in copy(nw_dict)
                if !(v isa MathematicalSolution{NetworkModel})
                    nw_dict[k] = MathematicalSolution(v)
                end
            end
        end

        return obj
    end
end


mutable struct EngineeringSolution{T<:Union{NetworkModel,MultinetworkModel}} <: DistributionModel{T}
    data::Dict{String,<:Any}

    function EngineeringSolution(data::Union{Dict{String,<:Any},Nothing}=nothing; multinetwork::Bool=false)
        if data isa Dict
            if get(data, "multinetwork", multinetwork)
                obj = new{MultinetworkModel}(data)
            else
                obj = new{NetworkModel}(data)
            end
        elseif data isa MathematicalSolution
            obj = data
        else
            error("Unrecognized data format, cannot instantiate MathematicalSolution")
        end

        if obj isa MathematicalSolution{MultinetworkModel}
            nw_dict = get!(() -> Dict{String,Any}(), obj.data, "nw")

            for (k, v) in copy(nw_dict)
                if !(v isa MathematicalSolution{NetworkModel})
                    nw_dict[k] = MathematicalSolution(v)
                end
            end
        end

        return obj
    end
end

"""
    solution_to_dict(solution)::Dict{String,Any}

Recursively converts a PMD `MathematicalSolution` or `EngineeringSolution`
to a `Dict{String,Any}`.

Nested solution objects, such as subnetworks under `"nw"`, are also
converted to dictionaries.
"""
function solution_to_dict(
    solution::Union{
        MathematicalSolution,
        EngineeringSolution,
    },
)::Dict{String,Any}
    return _solution_to_dict(solution.data)
end


function _solution_to_dict(data::AbstractDict)::Dict{String,Any}
    return Dict{String,Any}(
        string(k) => _solution_to_dict(v)
        for (k, v) in data
    )
end


function _solution_to_dict(
    solution::Union{
        MathematicalSolution,
        EngineeringSolution,
    },
)
    return _solution_to_dict(solution.data)
end


function _solution_to_dict(data::AbstractVector)
    return [_solution_to_dict(v) for v in data]
end


_solution_to_dict(x) = deepcopy(x)
StandardDataModels = Union{EngineeringModel,MathematicalModel,MathematicalSolution,EngineeringSolution}

Base.getindex(rd::T, key::String) where T<:(StandardDataModels) = rd.data[key]
Base.setindex!(rd::T, value, key::String) where T<:(StandardDataModels) = (rd.data[key] = value)
Base.haskey(rd::T, key::String) where T<:(StandardDataModels) = haskey(rd.data, key)
Base.get(rd::T, key::String) where T<:(StandardDataModels) = get(rd.data, key, nothing)
Base.get(rd::T, key::String, default) where T<:(StandardDataModels) = get(rd.data, key, default)
Base.keys(rd::T) where T<:(StandardDataModels) = keys(rd.data)
Base.values(rd::T) where T<:(StandardDataModels) = values(rd.data)
Base.iterate(rd::T, state...) where T<:(StandardDataModels) = iterate(rd.data, state...)
Base.iterate(rd::T, i::Int=1) where T<:StandardDataModels = iterate(rd.data, i)
Base.length(rd::T) where T<:(StandardDataModels) = length(rd.data)
Base.delete!(rd::T, key::String) where T<:(StandardDataModels) = delete!(rd.data, key)
Base.pop!(rd::T, key) where T<:(StandardDataModels) = pop!(rd.data, key)
# Base.convert(rd::T, ::U) where {T<:StandardDataModels,U<:Dict{String,<:Any}} = rd.data

transform_data_model(::Type{T}, data::EngineeringModel; kwargs...) where T<:MathematicalModel = transform_data_model(data; kwargs...)

