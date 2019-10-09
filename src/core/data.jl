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

"""
Returns bounds in LN base.
"""
function _bus_vm_ll_bounds(bus::Dict; eps=0.1)
    vmax = bus["vmax"].values
    vmin = bus["vmin"].values
    if haskey(bus, "vm_ll_max")
        vdmax = bus["vm_ll_max"].values*sqrt(3)
    else
        vdmax = 2*vmax
    end
    if haskey(bus, "vm_ll_min")
        vdmin = bus["vm_ll_min"].values*sqrt(3)
    else
        vdmin = ones(3)*eps*sqrt(3)
    end
    return (vdmin, vdmax)
end


function _load_pq_bounds(load::Dict, bus::Dict)
    a, α, b, β = _load_expmodel_params(load, bus)
    vmin, vmax = _load_vbounds(load, bus)
    # get bounds
    pmin = min.(a.*vmin.^α, a.*vmax.^α)
    pmax = max.(a.*vmin.^α, a.*vmax.^α)
    qmin = min.(b.*vmin.^β, b.*vmax.^β)
    qmax = max.(b.*vmin.^β, b.*vmax.^β)
    return (pmin, pmax, qmin, qmax)
end

function _load_curr_max(load::Dict, bus::Dict)
    pmin, pmax, qmin, qmax = _load_pq_bounds(load, bus)
    pabsmax = max.(abs.(pmin), abs.(pmax))
    qabsmax = max.(abs.(qmin), abs.(qmax))
    smax = sqrt.(pabsmax.^2 + qabsmax.^2)

    vmin, vmax = _load_vbounds(load, bus)

    return smax./vmin
end

function _load_curr_mag_bounds(load::Dict, bus::Dict)
    a, α, b, β = _load_expmodel_params(load, bus)
    vmin, vmax = _load_vbounds(load, bus)
    cb1 = sqrt.(a.^(2).*vmin.^(2*α.-2) + b.^(2).*vmin.^(2*β.-2))
    cb2 = sqrt.(a.^(2).*vmax.^(2*α.-2) + b.^(2).*vmax.^(2*β.-2))
    cmin = min.(cb1, cb2)
    cmax = max.(cb1, cb2)
    return cmin, cmax
end


function _load_expmodel_params(load::Dict, bus::Dict)
    pd = load["pd"].values
    qd = load["qd"].values
    ncnds = length(pd)
    if load["model"]=="constant_power"
        return (pd, zeros(ncnds), qd, zeros(ncnds))
    else
        # get exponents
        if load["model"]=="constant_current"
            α = ones(ncnds)
            β  =ones(ncnds)
        elseif load["model"]=="constant_impedance"
            α = ones(ncnds)*2
            β  =ones(ncnds)*2
        elseif load["model"]=="exponential"
            α = load["alpha"].values
            β = load["beta"].values
        end
        # calculate proportionality constants
        v0 = load["vnom_kv"]/(bus["base_kv"]/sqrt(3))
        a = pd./v0.^α
        b = qd./v0.^β
        # get bounds
        return (a, α, b, β)
    end
end

function _load_vbounds(load::Dict, bus::Dict)
    if load["conn"]=="wye"
        vmin = bus["vmin"].values
        vmax = bus["vmax"].values
    elseif load["conn"]=="delta"
        vmin, vmax = _bus_vm_ll_bounds(bus)
    end
    return vmin, vmax
end

"""
Returns a Bool, indicating whether the convex hull of the voltage-dependent
relationship needs a cone inclusion constraint.
"""
function _load_needs_cone(load::Dict)
    if load["model"]=="constant_current"
        return true
    elseif load["model"]=="exponential"
        return true
    else
        return false
    end
end


function _gen_curr_max(gen::Dict, bus::Dict)
    pabsmax = max.(abs.(gen["pmin"].values), abs.(gen["pmax"].values))
    qabsmax = max.(abs.(gen["qmax"].values), abs.(gen["qmax"].values))
    smax = sqrt.(pabsmax.^2 + qabsmax.^2)

    vmin = bus["vmin"].values

    return smax./vmin
end
