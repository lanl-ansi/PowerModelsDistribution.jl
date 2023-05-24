"Base.show variant for InfrastructureModel,InfrastructureObject to make them Dict-like"
function Base.show(@nospecialize(io::IO), @nospecialize(t::T)) where {T<:Union{InfrastructureModel,InfrastructureObject}}
    recur_io = IOContext(io, :SHOWN_SET => t,
        :typeinfo => eltype(t))

    limit = get(io, :limit, false)::Bool
    # show in a Julia-syntax-like form: Dict(k=>v, ...)
    print(io, Base.typeinfo_prefix(io, t)[1])
    print(io, '(')
    if !isempty(t) && !Base.show_circular(io, t)
        first = true
        n = 0
        for pair in t
            first || print(io, ", ")
            first = false
            show(recur_io, pair)
            n += 1
            limit && n >= 10 && (print(io, "…"); break)
        end
    end
    print(io, ')')
end

"Base.show variant for InfrastructureModel,InfrastructureObject to make them Dict-like"
function Base.show(@nospecialize(io::IO), ::MIME"text/plain", @nospecialize(t::T)) where {T<:Union{InfrastructureModel,InfrastructureObject}}
    isempty(t) && return show(io, t)
    # show more descriptively, with one line per key/value pair
    recur_io = IOContext(io, :SHOWN_SET => t)
    limit = get(io, :limit, false)::Bool
    if !haskey(io, :compact)
        recur_io = IOContext(recur_io, :compact => true)
    end
    recur_io_k = IOContext(recur_io, :typeinfo => keytype(t))
    recur_io_v = IOContext(recur_io, :typeinfo => valtype(t))

    summary(io, t)
    isempty(t) && return
    print(io, ":")
    Base.show_circular(io, t) && return
    if limit
        sz = displaysize(io)
        rows, cols = sz[1] - 3, sz[2]
        rows < 2 && (print(io, " …"); return)
        cols < 12 && (cols = 12) # Minimum widths of 2 for key, 4 for value
        cols -= 6 # Subtract the widths of prefix "  " separator " => "
        rows -= 1 # Subtract the summary

        # determine max key width to align the output, caching the strings
        ks = Vector{String}(undef, min(rows, length(t)))
        vs = Vector{String}(undef, min(rows, length(t)))
        keylen = 0
        vallen = 0
        for (i, (k, v)) in enumerate(t)
            i > rows && break
            ks[i] = sprint(show, k, context=recur_io_k, sizehint=0)
            vs[i] = sprint(show, v, context=recur_io_v, sizehint=0)
            keylen = clamp(length(ks[i]), keylen, cols)
            vallen = clamp(length(vs[i]), vallen, cols)
        end
        if keylen > max(div(cols, 2), cols - vallen)
            keylen = max(cld(cols, 3), cols - vallen)
        end
    else
        rows = cols = typemax(Int)
    end

    for (i, (k, v)) in enumerate(t)
        print(io, "\n  ")
        if i == rows < length(t)
            print(io, rpad("⋮", keylen), " => ⋮")
            break
        end

        if limit
            if Base.VERSION >= v"1.9.0"
                key = rpad(Base._truncate_at_width_or_chars(false, ks[i], keylen, false, "\r\n"), keylen)
            else
                key = rpad(Base._truncate_at_width_or_chars(ks[i], keylen), keylen)
            end
        else
            key = sprint(show, k, context=recur_io_k, sizehint=0)
        end
        print(recur_io, key)
        print(io, " => ")

        if limit
            if Base.VERSION >= v"1.9.0"
                val = Base._truncate_at_width_or_chars(false, vs[i], cols - keylen, false, "\r\n")
            else
                val = Base._truncate_at_width_or_chars(vs[i], cols - keylen)
            end
            print(io, val)
        else
            show(recur_io_v, v)
        end
    end
end

"Helper functions for Base.show for InfrastructureModel,InfrastructureObject"
_show(@nospecialize(io::IO), @nospecialize(t::T)) where {T<:Union{InfrastructureModel,InfrastructureObject}} = Base.show(io, t)
_show(@nospecialize(io::IO), m::MIME"text/plain", @nospecialize(t::T)) where {T<:Union{InfrastructureModel,InfrastructureObject}} = Base.show(io, m, t)
