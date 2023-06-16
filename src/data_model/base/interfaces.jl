"Base.get for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.get(@nospecialize(data::T), k::String, @nospecialize(v::Any)) where {T<:Union{InfrastructureModel,InfrastructureObject}} = hasproperty(data, Symbol(k)) ? getproperty(data, Symbol(k)) : v

"Base.getindex for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.getindex(@nospecialize(s::T), k::String) where {T<:Union{InfrastructureModel,InfrastructureObject}} = getproperty(s, Symbol(k))
Base.getindex(@nospecialize(s::T), k::Symbol) where {T<:Union{InfrastructureModel,InfrastructureObject}} = getproperty(s, k)

"Base.setindex! for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.setindex!(@nospecialize(s::T), v::U, k::String) where {U,T<:Union{InfrastructureModel,InfrastructureObject}} = setproperty!(s, Symbol(k), v)
Base.setindex!(@nospecialize(s::T), v::U, k::Symbol) where {U,T<:Union{InfrastructureModel,InfrastructureObject}} = setproperty!(s, k, v)

"Base.zero for String"
Base.zero(::Type{String})::String = ""

"Base.zero for Char"
Base.zero(::Type{Char})::Char = ' '

"Base.zero for ConnConfig enum"
Base.zero(::Type{ConnConfig})::ConnConfig = WYE

"Base.zero for SwitchState enum"
Base.zero(::Type{SwitchState})::SwitchState = OPEN

"Base.zero for Status enum"
Base.zero(::Type{Status})::Status = DISABLED

"Base.zero for Dispatchable enum"
Base.zero(::Type{Dispatchable})::Dispatchable = NO

"Base.zero for LoadModel enum"
Base.zero(::Type{LoadModel})::LoadModel = POWER

"Base.zero for ShuntModel enum"
Base.zero(::Type{ShuntModel})::ShuntModel = CAPACITOR

"Base.iterate for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.@propagate_inbounds function Base.iterate(@nospecialize(itr::T), i::Int=1) where {T<:Union{InfrastructureModel,InfrastructureObject}}
    pn = propertynames(itr)
    i > length(pn) && return nothing

    val = itr[pn[i]]
    i <= length(pn) ? (Base.@inbounds Pair{String,typeof(val)}(string(pn[i]), val), i + 1) : nothing
end

"Base.haskey for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.haskey(@nospecialize(h::T), key::String) where {T<:Union{InfrastructureModel,InfrastructureObject}} = (Symbol(key) ∈ [pn for pn in propertynames(h) if !isempty(getproperty(h, pn))])
Base.haskey(@nospecialize(h::T), key::Symbol) where {T<:Union{InfrastructureModel,InfrastructureObject}} = (key ∈ [pn for pn in propertynames(h) if !isempty(getproperty(h, pn))])

"Base.isempty for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.isempty(@nospecialize(h::T)) where {T<:Union{InfrastructureModel,InfrastructureObject}} = all(isempty(getproperty(h, pn)) for pn in propertynames(h))
Base.isempty(@nospecialize(h::Missing)) = true
Base.isempty(::Status) = false

"Base.keytype for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.keytype(@nospecialize(::InfrastructureModel)) = String
Base.keytype(@nospecialize(::InfrastructureObject)) = String

"Base.valtype for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.valtype(@nospecialize(h::T)) where {T<:Union{InfrastructureModel,InfrastructureObject}} = typeof(h)

"Base.eltype for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.eltype(@nospecialize(h::T)) where {T<:Union{InfrastructureModel,InfrastructureObject}} = typeof(h)

"Base.length for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.length(@nospecialize(X::T)) where {T<:Union{InfrastructureModel,InfrastructureObject}} = length(propertynames(X))

"Base.summary for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
function Base.summary(@nospecialize(io::IO), @nospecialize(t::T)) where {T<:Union{InfrastructureModel,InfrastructureObject}}
    Base.showarg(io, t, true)
    if Base.IteratorSize(t) isa Base.HasLength
        n = length(t)
        print(io, " with ", n, (n == 1 ? " entry" : " entries"))
    else
        print(io, "(...)")
    end
end

"Base.merge! for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
function Base.merge!(@nospecialize(x::T), @nospecialize(y::T)) where {T<:DssObject}
    for pn in propertynames(y)
        if pn ∉ [:switch, :name, :enabled]  # global exclusions for applying "like"
            setproperty!(x, pn, getproperty(y, pn))
        end
    end
end

"Base.getproperty for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.getproperty(@nospecialize(x::Any), k::String) = getproperty(x, Symbol(k))

"Base.setproperty! for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.setproperty!(@nospecialize(x::Any), k::String, @nospecialize(v::Any)) = setproperty!(x, Symbol(k), v)

"Base.keys for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.keys(@nospecialize(x::T)) where {T<:Union{InfrastructureObject,InfrastructureModel}} = string.(collect(propertynames(x)))

"Base.delete! for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.delete!(@nospecialize(x::T), k::String) where {T<:InfrastructureObject} = delete!(x, Symbol(k))

"Base.delete! for InfrastructureModel,InfrastructureObject to give them Dict-like behavior"
Base.delete!(@nospecialize(x::T), k::Symbol) where {T<:InfrastructureObject} = setproperty!(x, k, missing)
