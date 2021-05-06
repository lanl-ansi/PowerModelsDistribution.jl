const _pmd_eng_global_keys = Set{String}([
    "files", "name", "data_model", "dss_options"
])


"list of multinetwork keys that belong at the root level"
const _pmd_math_global_keys = Set{String}([
    "data_model", "per_unit", "name", "map", "bus_lookup", "dss_options"
])


"""
    ismultinetwork(pm::AbstractUnbalancedPowerModel)

Checks if power model struct is multinetwork
"""
ismultinetwork(pm::AbstractUnbalancedPowerModel) = ismultinetwork(pm, pmd_it_sym)


"""
    make_multinetwork(
        data::Dict{String,<:Any};
        sparse::Bool=false,
        time_elapsed::Union{Missing,Real,Vector{<:Real}}=missing,
        global_keys::Set{String}=Set{String}(),
    )::Dict{String,Any}

Expands a data structure into a multinetwork assuming there are `time_series` objects defined and assigned
to some components.

If `global_keys::Set{String}` is defined, the global keys that will appear at the top-level of the data
structure will include both the default global keys for that data model type, and additionally the keys
defined in `global_keys`.

`time_elapsed` defines the time elapsed between subnetworkin __hours__, and can either be a single Real value, and
thus a constant time between each time step, or a Vector with the same length as the number of time steps,
or can be left missing, in which case time elapsed will attempt to be discovered, with a fallback on 1 hour default.
Time elapsed can be adjusted later using [`set_time_elapsed!`](@ref set_time_elapsed!)

`make_multinetwork` assumes all `"time"` values in "time_series" objects to be in the same units, and will attempt
to automatically sort multinetworks in the correct order. [`sort_multinetwork!`](@ref sort_multinetwork!) can be
used after the fact to re-sort the subnetworks.

`sparse` is currently unsupported, and is only included for future compatibility
"""
function make_multinetwork(
    data::Dict{String,<:Any};
    sparse::Bool=false,
    time_elapsed::Union{Missing,Real,Vector{<:Real}}=missing,
    global_keys::Set{String}=Set{String}(),
    )::Dict{String,Any}

    if ismultinetwork(data)
        @info "data is already multinetwork, returning"
        return data
    elseif iseng(data)
        return _make_multinetwork_eng(data; sparse=sparse, time_elapsed=time_elapsed, global_keys=global_keys)
    elseif ismath(data)
        return _make_multinetwork_math(data; global_keys=global_keys)
    end
end


"Expands an ENGINEERING data structure into a multinetwork, see [`make_multinetwork`](@ref make_multinetwork)"
function _make_multinetwork_eng(
    data_eng::Dict{String,<:Any};
    sparse::Bool=false,
    time_elapsed::Union{Missing,Real,Vector{<:Real}}=missing,
    global_keys::Set{String}=Set{String}(),
    )::Dict{String,Any}

    @assert iseng(data_eng) "wrong data model for _make_multinetwork_eng"

    mn_data = Dict{String,Any}(
        k => data_eng[k] for k in union(_pmd_eng_global_keys, global_keys) if haskey(data_eng, k)
    )
    mn_data["multinetwork"] = true

    ts_lookup = Dict(
        type => Dict(
            name => obj["time_series"] for (name, obj) in get(data_eng, type, Dict()) if !isempty(get(obj, "time_series", Dict()))
        ) for type in pmd_eng_asset_types
    )

    times = Set()
    for (type, ts_info) in ts_lookup
        for ts_dict in values(ts_info)
            for ts_id in values(ts_dict)
                if haskey(get(data_eng, "time_series", Dict()), ts_id)
                    push!(times, data_eng["time_series"][ts_id]["time"]...)
                end
            end
        end
    end

    if eltype(times) == Any && length(unique([typeof(el) for el in times])) > 1
        @warn "Time steps in this data model are of mixed type, multinetwork frames may be out of order. You may want to manually correct this with sort_multinetwork!"
        times = [t for t in times]
        if ismissing(time_elapsed)
            @warn "Cannot infer time_elapsed, setting to 1.0. You can manually adjust with set_time_elapsed!"
            time_elapsed = 1.0
        end
    else
        times = sort([t for t in times])
        if ismissing(time_elapsed)
            if eltype(times) <: Real
                @info "assuming time is in hours for time_elapsed inference. if this is incorrect, manually adjust with set_time_elapsed!"
                time_elapsed = diff(times)
                push!(time_elapsed, time_elapsed[end])
            else
                try
                    date_times = [Dates.DateFormat(t) for t in times]
                    delta = diff(date_times)
                    push!(delta, delta[end])
                    time_elapsed = [Dates.toms(dt) / 3.6e6 for dt in delta]
                catch err
                    @warn "Cannot infer time_elapsed, setting to 1.0. You can manually adjust with set_time_elapsed!"
                    time_elapsed = 1.0
                end
            end
        end
    end

    if length(times) == 0
        mn_data["nw"] = Dict{String,Any}(
            "0" => Dict{String,Any}(
                k => v for (k,v) in data_eng if !(k in _pmd_eng_global_keys)
            )
        )
        mn_data["nw"]["0"]["time"] = 0
        mn_data["mn_lookup"] = Dict{String,Any}("0" => 0)
    else
        mn_lookup = Dict("$i" => time for (i,time) in enumerate(times))
        mn_data["nw"] = Dict{String,Any}(
            i => Dict{String,Any}(
                "time" => time
            ) for (i,time) in mn_lookup
        )
        mn_data["mn_lookup"] = mn_lookup

        _nw = Dict{String,Any}(
            k => v for (k,v) in data_eng if !(k in _pmd_eng_global_keys)
        )
        for n in sort([n for n in keys(mn_data["nw"])])
            time = mn_lookup["$n"]
            _nw["time"] = time

            for (type, ts_info) in ts_lookup
                for (obj_name, ts_dict) in ts_info
                    for (property_name, ts_name) in ts_dict
                        if haskey(get(data_eng, "time_series", Dict()), ts_name)
                            ts = data_eng["time_series"][ts_name]
                            if time in ts["time"]
                                idx = findfirst(x->x==time, ts["time"])
                                if ts["replace"]
                                    _nw[type][obj_name][property_name] = zeros(size(_nw[type][obj_name][property_name])) .+ ts["values"][idx]
                                else
                                    _nw[type][obj_name][property_name] = deepcopy(data_eng[type][obj_name][property_name]) .* ts["values"][idx]
                                end
                            end
                        end
                    end
                end
            end
            merge!(mn_data["nw"]["$n"], deepcopy(_nw))
        end
    end

    return mn_data
end


"Expands an MATHEMATICAL data structure into a multinetwork, see [`make_multinetwork`](@ref make_multinetwork)"
function _make_multinetwork_math(
    data_math::Dict{String,<:Any};
    global_keys::Set{String}=Set{String}(),
    )::Dict{String,Any}

    @assert ismath(data_math) "wrong data model for _make_multinetwork_math"
    @info "Converting a MATHEMATICAL data model into multinetwork assumes the InfrastructureModels data structure ( see https://github.com/lanl-ansi/InfrastructureModels.jl/blob/master/src/core/data.jl#L135 )"

    return _IM.make_multinetwork(data_math, pmd_it_name, union(_pmd_math_global_keys, global_keys))
end


"""
    sort_multinetwork!(mn_data::Dict{String,<:Any}, times::Vector{<:Any})

Helper function to manually sort your multinetwork frames, given some pre-sorted vector of time values `times`
"""
function sort_multinetwork!(mn_data::Dict{String,<:Any}, times::Vector{<:Any})
    @assert ismultinetwork(mn_data) "Data is not multinetwork, cannot sort"
    mn_lookup = Dict("$i" => time for (i,time) in enumerate(times))
    mn_lookup_reverse = Dict(nw["time"] => n for (n,nw) in mn_data["nw"])

    _nw_data = Dict{String,Any}(
        new_n => deepcopy(mn_data["nw"][mn_lookup_reverse[time]]) for (new_n, time) in mn_lookup
    )

    mn_data["nw"] = _nw_data
    mn_data["mn_lookup"] = mn_lookup
end


"""
    set_time_elapsed!(data::Dict{String,<:Any}, time_elapsed::Union{Real,Vector{<:Real}})

Helper function to set time_elapsed in multinetwork data, given either some constant value of time elapsed
or a Vector of time elapsed values of the same length as the number of subnetworks.
"""
function set_time_elapsed!(data::Dict{String,<:Any}, time_elapsed::Union{Real,Vector{<:Real}})
    if ismultinetwork(data)
        nw_data = data["nw"]
    else
        nw_data = Dict("0" => data)
    end

    if isa(time_elapsed, Vector)
        @assert length(time_elapsed) == length(nw_data) "provided time_elapsed length doesn't match number of multinetwork frames"
    end

    for (n,nw) in nw_data
        if isa(time_elapsed, Vector)
            nw["time_elapsed"] = popfirst!(time_elapsed)
        else
            nw["time_elapsed"] = time_elapsed
        end
    end
end
