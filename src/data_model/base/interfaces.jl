""
Base.getindex(@nospecialize(s::Union{InfrastructureDataModel,GenericInfrastructureObject}), k::String) = getproperty(s, Symbol(k))
Base.getindex(@nospecialize(s::Union{InfrastructureDataModel,GenericInfrastructureObject}), k::Symbol) = getproperty(s, k)

Base.setindex!(@nospecialize(s::T), v::U, k::String) where {U, T <: Union{DssObject, EngObject, MathObject} } = setproperty!(s, Symbol(k), v)
Base.setindex!(@nospecialize(s::T), v::U, k::Symbol) where {U, T <: Union{DssObject, EngObject, MathObject} } = setproperty!(s, k, v)

Base.zero(::Type{String})::String = ""
Base.zero(::Type{Char})::Char = ' '
Base.zero(::Type{ConnConfig})::ConnConfig = WYE
Base.zero(::Type{SwitchState})::SwitchState = OPEN

Base.@propagate_inbounds function Base.iterate(@nospecialize(itr::Union{InfrastructureDataModel,GenericInfrastructureObject}), i::Int=1)
    pn = propertynames(itr)
    i > length(pn) && return nothing

    val = itr[pn[i]]
    i <= length(pn) ? (Base.@inbounds Pair{String, typeof(val)}(string(pn[i]), val), i+1) : nothing
end

Base.haskey(@nospecialize(h::Union{InfrastructureDataModel,GenericInfrastructureObject}), key::String) = (Symbol(key) ∈ [pn for pn in propertynames(h) if !isempty(getproperty(h, pn))])
Base.haskey(@nospecialize(h::Union{InfrastructureDataModel,GenericInfrastructureObject}), key::Symbol) = (key ∈ [pn for pn in propertynames(h) if !isempty(getproperty(h, pn))])

Base.isempty(@nospecialize(h::Union{InfrastructureDataModel,GenericInfrastructureObject})) = all(isempty(getproperty(h, pn)) for pn in propertynames(h))
Base.isempty(@nospecialize(h::Missing)) = true
Base.isempty(::Status) = false

Base.keytype(@nospecialize(::InfrastructureDataModel)) = String
Base.keytype(@nospecialize(::GenericInfrastructureObject)) = String

Base.valtype(@nospecialize(h::Union{InfrastructureDataModel,GenericInfrastructureObject})) = typeof(h)

Base.eltype(@nospecialize(h::Union{InfrastructureDataModel,GenericInfrastructureObject})) = typeof(h)

Base.length(@nospecialize(X::T)) where T <: Union{InfrastructureDataModel,GenericInfrastructureObject} = length(propertynames(X))


function Base.summary(io::IO, @nospecialize(t::Union{InfrastructureDataModel,GenericInfrastructureObject}))
    Base.showarg(io, t, true)
    if Base.IteratorSize(t) isa Base.HasLength
        n = length(t)
        print(io, " with ", n, (n==1 ? " entry" : " entries"))
    else
        print(io, "(...)")
    end
end


function Base.merge!(a::EngShuntControls, b::EngShuntControls)
    for (pn,a_property) in a
        b_property = getproperty(b, pn)
        if ismissing(a_property)
            setproperty!(a, pn, b_property)
        elseif isa(a_property, Vector)
            push!(a_property, b_property...)
        end
    end
end


function Base.merge!(a::EngTransformerControls, b::EngTransformerControls)
    if ismissing(a.dss) && !ismissing(b.dss)
        a.dss = Vector{Union{Missing,DssRegcontrol}}[Union{Missing,DssRegcontrol}[missing]]
    end

    for (w, wdg) in enumerate(b.windings)
        if wdg ∉ a.windings
            for (pn, prop) in a
                if pn == "dss" && !ismissing(b.dss)
                    push!(prop, getproperty(b, Symbol(pn))[w])
                end
            end
        else
            wdg_idx = findfirst(isequal(wdg), a.windings)
            for (t, term) in enumerate(b.terminals[w])
                if term ∉ a.terminals[wdg_idx]
                    for (pn, prop) in a
                        if pn ∉ ["windings","dss"]
                            push!(prop[wdg_idx], getproperty(b, Symbol(pn))[w][t])
                        elseif pn == "dss"
                            if !ismissing(getproperty(b, Symbol(pn)))
                                push!(prop[wdg_idx], getproperty(b, Symbol(pn))[w][t])
                            end
                        end
                    end
                else
                    term_idx = findfirst(isequal(term), a.terminals[wdg_idx])
                    for (pn, prop) in a
                        if pn ∉ ["windings", "terminals", "dss"]
                            prop[wdg_idx][term_idx] = getproperty(b, Symbol(pn))[w][t]
                        elseif pn == "dss"
                            if !ismissing(getproperty(b, Symbol(pn)))
                                prop[wdg_idx][term_idx] = getproperty(b, Symbol(pn))[w][t]
                            end
                        end
                    end
                end
            end
        end
    end
end


"""
"""
function Base.merge!(@nospecialize(x::T), @nospecialize(y::T)) where T <: DssObject
    for pn in propertynames(y)
        if pn ∉ [:switch, :name, :enabled]  # global exclusions for applying "like"
            setproperty!(x, pn, getproperty(y, pn))
        end
    end
end
