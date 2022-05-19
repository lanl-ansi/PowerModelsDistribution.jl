abstract type InfrastructureDataModel{T} end
abstract type MultinetworkDataModel end
abstract type SimpleDataModel end
abstract type GenericObject end

abstract type NodeObject <: GenericObject end
abstract type EdgeObject <: GenericObject end
abstract type DataObject <: GenericObject end

abstract type DistributionDataModel{T} <: InfrastructureDataModel{T} end
abstract type EngineeringModel{T} <: DistributionDataModel{T} end

abstract type EngNodeObject <: NodeObject end
abstract type EngEdgeObject <: EdgeObject end
abstract type EngDataObject <: DataObject end

EngObject = Union{EngNodeObject,EngEdgeObject,EngDataObject}

abstract type EngBus <: EngNodeObject end

Base.@kwdef mutable struct EngBusObj <: EngBus
    name::String = "default_name"
end

function get_name(s::GenericObject)
    return s.name
end

function create_object(::Type{T}) where T <: EngBus
    EngBusObj("name")
end


macro kwdef(expr)
    expr = macroexpand(__module__, expr) # to expand @static
    expr isa Expr && expr.head === :struct || error("Invalid usage of @kwdef")
    expr = expr::Expr
    T = expr.args[2]
    if T isa Expr && T.head === :<:
        T = T.args[1]
    end

    params_ex = Expr(:parameters)
    call_args = Any[]

    Base._kwdef!(expr.args[3], params_ex.args, call_args)
    # Only define a constructor if the type has fields, otherwise we'll get a stack
    # overflow on construction
    if !isempty(params_ex.args)
        if T isa Symbol
            kwdefs = :(($(esc(T)))($params_ex) = ($(esc(T)))($(call_args...)))
        elseif T isa Expr && T.head === :curly
            T = T::Expr
            # if T == S{A<:AA,B<:BB}, define two methods
            #   S(...) = ...
            #   S{A,B}(...) where {A<:AA,B<:BB} = ...
            S = T.args[1]
            P = T.args[2:end]
            Q = Any[U isa Expr && U.head === :<: ? U.args[1] : U for U in P]
            SQ = :($S{$(Q...)})
            kwdefs = quote
                ($(esc(S)))($params_ex) =($(esc(S)))($(call_args...))
                ($(esc(SQ)))($params_ex) where {$(esc.(P)...)} =
                    ($(esc(SQ)))($(call_args...))
            end
        else
            error("Invalid usage of @kwdef")
        end
    else
        kwdefs = nothing
    end
    quote
        Base.@__doc__($(esc(expr)))
        $kwdefs
    end
end


macro extend(base, name, fields)
    for arg in fields.args
        println(dump(arg))
    end
    base_type = Core.eval(@__MODULE__, base)

    base_fieldnames = fieldnames(base_type)
    base_types = [t for t in base_type.types]

    base_fields = [:($f::$T) for (f,T) in zip(base_fieldnames, base_types)]

    res = :(mutable struct $name <: $(supertype(base_type)) end)

    push!(res.args[end].args, base_fields...)
    push!(res.args[end].args, fields.args...)

    return res
end

@extend EngBusObj TestBusObj begin
    extra::Vector{<:Real} = Real[]
    source_id::String = "default_source_id"
end

function get_source_id(s::TestBusObj)
    return s.source_id
end

extbus = TestBusObj("name", [1, 2, 3], "sourceid")

get_source_id(extbus)

typeof(extbus)

extbus <: EngBus

get_name(extbus)

params = Dict("name" => Dict("type"=>"string","default"=>"default_name"))

macro gen(sname, params)
    println(dump($(esc(params))))
end

@gen Test params
