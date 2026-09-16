struct ReferenceProxy
    ref_string::String
    target::Dict{String,Any}
end

Base.getindex(p::ReferenceProxy, key) = p.target[key]
Base.setindex!(p::ReferenceProxy, value, key) = (p.target[key] = value)
Base.get(p::ReferenceProxy, key, default) = get(p.target, key, default)
Base.keys(p::ReferenceProxy) = keys(p.target)
Base.values(p::ReferenceProxy) = values(p.target)
Base.iterate(p::ReferenceProxy, state...) = iterate(p.target, state...)
Base.show(io::IO, p::ReferenceProxy) = print(io, p.target)
Base.show(io::IO, ::MIME"text/plain", p::ReferenceProxy) = show(io, MIME("text/plain"), p.target)

mutable struct RavensModel{T<:NetworkModel} <: DistributionModel{T}
    data::Dict{String,Any}
    unraveled::Dict{String,Any}
    paths::Dict{String,Any}
    iter::Dict{String,Dict{String,Any}}
    containers::Dict{String,Any}  # Maps container names to all objects underneath
    container_paths::Dict{String,Vector{String}}  # Maps container names to their path in the tree

    function RavensModel(network_profile=nothing)
        obj = new{NetworkModel}(
            Dict{String,Any}(),
            Dict{String,Any}(),
            Dict{String,Any}(),
            Dict{String,Dict{String,Any}}(),
            Dict{String,Dict{String,Any}}(),
            Dict{String,Vector{String}}()
        )

        if network_profile isa Dict
            obj.data = copy(network_profile)
        elseif network_profile isa String
            open(network_profile, "r") do f
                obj.data = JSON.parse(f)
            end
        elseif network_profile isa IO
            obj.data = JSON.parse(network_profile)
        end

        build_paths!(obj)
        build_containers!(obj)
        build_iterates!(obj)
        update_refs!(obj)

        return obj
    end
end


"""
    safe_get_container(rd::RavensModel, path::Vector{String}, default=Dict{String,Any}())

Safely navigate nested container paths, returning default if any part doesn't exist.
"""
function safe_get_container(rd::RavensModel, path::Vector{String}, default=Dict{String,Any}())
    result = rd.data
    for key in path
        if result isa AbstractDict && haskey(result, key)
            result = result[key]
        else
            return default
        end
    end
    return result
end

"""
    safe_get_container(rd::RavensModel, container::String, default=Dict{String,Any}())

Get a container by name, returning default if it doesn't exist.
"""
function safe_get_container(rd::RavensModel, container::String, default=Dict{String,Any}())
    get(rd.containers, container, default)
end


Base.show(io::IO, rd::RavensModel) = print(io, "RavensModel with $(length(rd.containers)) collections and $(length(rd.unraveled)) objects")
Base.show(io::IO, ::MIME"text/plain", rd::RavensModel) = begin
    println(io, "RavensModel:")
    println(io, "  Collections: $(length(rd.containers))")
    println(io, "  Objects: $(length(rd.unraveled))")
    if !isempty(rd.containers)
        println(io, "  Available collections:")
        for (name, objs) in sort(collect(rd.containers), by=x -> x[1])
            println(io, "    $name: $(length(objs)) objects")
        end
    end
end

# Enable tab completion for RavensModel - show containers, not data keys
function Base.propertynames(rd::RavensModel, private::Bool=false)
    public_props = tuple(Symbol.(keys(rd.containers))...)
    if private
        return (public_props..., :data, :unraveled, :paths, :iter, :containers, :container_paths)
    else
        return public_props
    end
end

function Base.getproperty(rd::RavensModel, name::Symbol)
    if name in (:data, :unraveled, :paths, :iter, :containers, :container_paths)
        return getfield(rd, name)
    else
        key = String(name)
        if haskey(getfield(rd, :containers), key)
            return getfield(rd, :containers)[key]
        else
            throw(KeyError(key))
        end
    end
end

# Enable tab completion for bracket notation
if isdefined(Base, :key_completion_hint)
    Base.key_completion_hint(rd::RavensModel) = keys(rd.containers)
end

function add_path!(rd::RavensModel, data::AbstractDict)
    for (k, v) in data
        if v isa AbstractDict
            cim_obj_type = get(v, "Ravens.cimObjectType", nothing)
            name = get(v, "IdentifiedObject.name", nothing)
            mrid = get(v, "IdentifiedObject.mRID", nothing)

            if !isnothing(cim_obj_type)
                isnothing(mrid) && (mrid = string(UUIDs.uuid4()))
                isnothing(name) && (name = mrid)

                rd.unraveled[mrid] = v
                rd.paths["$cim_obj_type::'$name'"] = v
                name != mrid && (rd.paths["$cim_obj_type::'$mrid'"] = v)
            end

            add_path!(rd, v)
        elseif v isa AbstractVector
            for item in v
                if item isa AbstractDict
                    cim_obj_type = get(item, "Ravens.cimObjectType", nothing)
                    name = get(item, "IdentifiedObject.name", nothing)
                    mrid = get(item, "IdentifiedObject.mRID", nothing)

                    if !isnothing(cim_obj_type)
                        isnothing(mrid) && (mrid = string(UUIDs.uuid4()))
                        isnothing(name) && (name = mrid)

                        rd.unraveled[mrid] = item
                        rd.paths["$cim_obj_type::'$name'"] = item
                        name != mrid && (rd.paths["$cim_obj_type::'$mrid'"] = item)
                    end

                    add_path!(rd, item)
                end
            end
        end
    end
end

function build_paths!(rd::RavensModel)
    add_path!(rd, rd.data)
end

function collect_objects_from_container!(container_name::String, data::AbstractDict, objects::Dict{String,Any})
    """Recursively collect all objects with Ravens.cimObjectType from a container dict"""
    for (k, v) in data
        if v isa AbstractDict
            if haskey(v, "Ravens.cimObjectType")
                # This is an object with a type - add it
                obj_name = get(v, "IdentifiedObject.name", k)
                objects[obj_name] = v
                @debug "Added $(v["Ravens.cimObjectType"])::$obj_name to container $container_name"
            else
                # No type - keep recursing to find objects
                collect_objects_from_container!(container_name, v, objects)
            end
        end
    end
end

function any_typed_descendants(data::AbstractDict)
    """Check if this dict or any of its descendants have Ravens.cimObjectType"""
    for (k, v) in data
        if v isa AbstractDict
            if haskey(v, "Ravens.cimObjectType")
                return true
            elseif any_typed_descendants(v)
                return true
            end
        end
    end
    return false
end

function build_containers!(rd::RavensModel, data=nothing, path::Vector{String}=String[])
    """Build containers dict by identifying container nodes and collecting all objects beneath them"""
    isnothing(data) && (data = rd.data)

    if data isa AbstractDict
        has_cim_type = haskey(data, "Ravens.cimObjectType")

        if !has_cim_type
            # Check if any descendants have Ravens.cimObjectType
            has_typed_descendants = false
            for (k, v) in data
                if v isa AbstractDict
                    if haskey(v, "Ravens.cimObjectType")
                        has_typed_descendants = true
                        break
                    else
                        # Recurse to check descendants
                        if any_typed_descendants(v)
                            has_typed_descendants = true
                            break
                        end
                    end
                end
            end

            # If this dict doesn't have a type but has descendants with types, it's a container
            if has_typed_descendants && !isempty(path)
                container_name = path[end]
                if !haskey(rd.containers, container_name)
                    rd.containers[container_name] = Dict{String,Any}()
                    rd.container_paths[container_name] = copy(path)
                end
                collect_objects_from_container!(container_name, data, rd.containers[container_name])
                @debug "Built container: $container_name with $(length(rd.containers[container_name])) objects at path $(join(path, " > "))"
            end

            # Recurse into children
            for (k, v) in data
                if v isa AbstractDict
                    build_containers!(rd, v, vcat(path, [k]))
                end
            end
        end
    end
end

function build_iterates!(rd::RavensModel)
    ref_pattern = r"(\w+)::'(.+)'"

    for (k, v) in rd.paths
        m = match(ref_pattern, k)
        if !isnothing(m)
            uobj = m.captures[1]
            obj_id = m.captures[2]

            if !haskey(rd.iter, uobj)
                rd.iter[uobj] = Dict{String,Any}()
            end

            rd.iter[uobj][obj_id] = v
        end
    end
end

function update_refs!(rd::RavensModel, obj=nothing)
    ref_pattern = r"(\w+)::'(.+)'"
    isnothing(obj) && (obj = rd.data)

    if obj isa AbstractDict
        for (k, v) in collect(obj)  # collect to avoid mutation during iteration
            if v isa AbstractString && !isnothing(match(ref_pattern, v))
                if haskey(rd.paths, v)
                    obj[k] = ReferenceProxy(v, rd.paths[v])
                    @debug "Created proxy for $v"
                else
                    @warn "Reference not found: $v"
                end
            elseif v isa Union{AbstractDict,AbstractVector}
                update_refs!(rd, v)
            end
        end
    elseif obj isa AbstractVector
        for i in eachindex(obj)
            item = obj[i]
            if item isa AbstractString && !isnothing(match(ref_pattern, item))
                if haskey(rd.paths, item)
                    obj[i] = ReferenceProxy(item, rd.paths[item])
                    @debug "Created proxy for $item"
                else
                    @warn "Reference not found: $item"
                end
            elseif item isa Union{AbstractDict,AbstractVector}
                update_refs!(rd, item)
            end
        end
    end
end

function serialize_refs(obj)
    if obj isa ReferenceProxy
        return obj.ref_string
    elseif obj isa Dict
        return Dict(k => serialize_refs(v) for (k, v) in obj)
    elseif obj isa Vector
        return [serialize_refs(item) for item in obj]
    else
        return obj
    end
end

# Collections-first interface
Base.getindex(rd::RavensModel, key::String) = rd.containers[key]
Base.setindex!(rd::RavensModel, value, key::String) = (rd.containers[key] = value)
Base.haskey(rd::RavensModel, key::String) = haskey(rd.containers, key)
Base.get(rd::RavensModel, key::String) = get(rd.containers, key, nothing)
Base.get(rd::RavensModel, key::String, default) = get(rd.containers, key, default)
Base.keys(rd::RavensModel) = keys(rd.containers)
Base.values(rd::RavensModel) = values(rd.containers)
Base.iterate(rd::RavensModel, state...) = iterate(rd.containers, state...)
Base.length(rd::RavensModel) = length(rd.containers)

ref(rd::RavensModel, key::String) = rd.paths[key]

function dump_ravens(rd::RavensModel, file_path::AbstractString; indent::Union{Int,Nothing}=nothing)
    open(file_path, "w") do f
        serialized = serialize_refs(rd.data)
        JSON.print(f, serialized, isnothing(indent) ? 0 : indent)
    end
end

function _deep_merge!(dst::Dict{String,Any}, src::Dict{String,Any})
    for (k, v_src) in src
        if haskey(dst, k)
            v_dst = dst[k]
            if v_dst isa AbstractDict && v_src isa AbstractDict
                # Recurse into the nested dicts
                _deep_merge!(v_dst, v_src)
            else
                # Overwrite with the source value (including scalars, vectors, etc.)
                dst[k] = v_src
            end
        else
            # New key – just copy it over
            dst[k] = v_src
        end
    end
    return dst
end

"""
    Base.merge(rd::RavensModel, extra::Dict)

Return a new `RavensModel` whose `data` is the result of recursively merging
`rd.data` with `extra`.  All internal tables (`paths`, `containers`,
`iter`, reference proxies) are rebuilt so the model behaves exactly like one
constructed from the merged dictionary.
"""
function Base.merge(rd::RavensModel, extra::Dict)
    merged_data = deepcopy(rd.data) # work on a copy to keep `rd` unchanged
    _deep_merge!(merged_data, extra)

    new_rd = RavensModel() # empty model
    new_rd.data = merged_data # replace its data with the merged result

    empty!(new_rd.unraveled)
    empty!(new_rd.paths)
    empty!(new_rd.containers)
    empty!(new_rd.container_paths)
    empty!(new_rd.iter)

    build_paths!(new_rd)
    build_containers!(new_rd)
    build_iterates!(new_rd)
    update_refs!(new_rd)

    return new_rd
end

transform_data_model(::Type{T}, data::RavensModel; kwargs...) where T <: MathematicalModel = transform_data_model(data; kwargs...)
