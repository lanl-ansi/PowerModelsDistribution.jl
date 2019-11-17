"wraps angles in degrees to 180"
function _wrap_to_180(degrees)
    return degrees - 360*floor.((degrees .+ 180)/360)
end


"wraps angles in radians to pi"
function _wrap_to_pi(radians)
    return radians - 2*pi*floor.((radians .+ pi)/(2*pi))
end


"rolls a 1d array left or right by idx"
function _roll(array::Array{T, 1}, idx::Int; right=true) where T <: Number
    out = Array{T}(undef, size(array))
    pos = idx % length(out)

    if right
        out[1+pos:end] = array[1:end-pos]
        out[1:pos] = array[end-(pos-1):end]
    else
        out[1:end-pos] = array[1+pos:end]
        out[end-(pos-1):end] = array[1:pos]
    end

    return out
end


"Corrects the shunts from vectors to matrices after the call to PMs."
function make_multiconductor!(mp_data, n_conductors::Int)
    PowerModels.make_multiconductor!(mp_data, n_conductors)
    # replace matrix shunts by matrices instead of vectors
    for (_, br) in mp_data["branch"]
        for key in ["b_fr", "b_to", "g_fr", "g_to"]
            br[key] = _PMs.MultiConductorMatrix(LinearAlgebra.diagm(0=>br[key].values))
        end
    end
end


"Replaces NaN values with zeros"
_replace_nan(v) = map(x -> isnan(x) ? zero(x) : x, v)


"Counts number of nodes in network"
function count_nodes(dss_data::Dict{String,Array})::Int
    sourcebus = get(dss_data["circuit"][1], "bus1", "sourcebus")
    all_nodes = Dict()
    for comp_type in values(dss_data)
        for comp in values(comp_type)
            if isa(comp, Dict) && haskey(comp, "buses")
                for busname in values(_parse_array(String, comp["buses"]))
                    name, nodes = _parse_busname(busname)

                    if !haskey(all_nodes, name)
                        all_nodes[name] = Set([])
                    end

                    for (n, node) in enumerate(nodes[1:3])
                        if node
                            push!(all_nodes[name], n)
                        end
                    end
                end
            elseif isa(comp, Dict)
                for (prop, val) in comp
                    if startswith(prop, "bus") && prop != "buses"
                        name, nodes = _parse_busname(val)

                        if !haskey(all_nodes, name)
                            all_nodes[name] = Set([])
                        end

                        for (n, node) in enumerate(nodes[1:3])
                            if node
                                push!(all_nodes[name], n)
                            end
                        end
                    end
                end
            end
        end
    end

    n_nodes = 0
    for (name, phases) in all_nodes
        if name != sourcebus
            n_nodes += length(phases)
        end
    end

    return n_nodes
end


"Counts number of nodes in network"
function count_nodes(pmd_data::Dict{String,Any})::Int
    if pmd_data["source_type"] == "dss"
        Memento.info(_LOGGER, "counting nodes from PowerModelsDistribution structure may not be as accurate as directly from `parse_dss` data due to virtual buses, etc.")
    end

    n_nodes = 0
    for bus in values(pmd_data["bus"])
        if pmd_data["source_type"] == "matlab"
            n_nodes += sum(bus["vm"] .> 0.0)
        elseif pmd_data["source_type"] == "dss"
            if !(pmd_data["source_type"] == "dss" && bus["name"] in ["virtual_sourcebus", pmd_data["sourcebus"]]) && !(pmd_data["source_type"] == "dss" && startswith(bus["name"], "tr") && endswith(bus["name"], r"_b\d"))
                n_nodes += sum(bus["vm"] .> 0.0)
            end
        end
    end

    return n_nodes
end


"Calculates the tap scale factor for the non-dimensionalized equations."
function calculate_tm_scale(trans::Dict{String,Any}, bus_fr::Dict{String,Any}, bus_to::Dict{String,Any})
    f_vnom = trans["config_fr"]["vm_nom"]
    t_vnom = trans["config_to"]["vm_nom"]
    f_vbase = bus_fr["base_kv"]
    t_vbase = bus_to["base_kv"]
    f_type = trans["config_fr"]["type"]
    t_type = trans["config_to"]["type"]

    tm_scale = (f_vnom/t_vnom)*(t_vbase/f_vbase)
    if f_type == "delta"
        tm_scale *= sqrt(3)
    end
    if t_type == "delta"
        tm_scale *= 1/sqrt(3)
    end
    if f_type == "zig-zag"
        Memento.error("Zig-zag not yet supported.")
    end
    if t_type == "zig-zag"
        Memento.error("Zig-zag not yet supported.")
    end

    return tm_scale
end
