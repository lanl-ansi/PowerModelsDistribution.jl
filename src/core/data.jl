"wraps angles in degrees to 180"
function _wrap_to_180(degrees)
    return degrees - 360*floor.((degrees .+ 180)/360)
end


"wraps angles in radians to pi"
function _wrap_to_pi(radians)
    return radians - 2*pi*floor.((radians .+ pi)/(2*pi))
end


"creates a delta transformation matrix"
function _get_delta_transformation_matrix(n_phases::Int)
    @assert(n_phases>2, "We only define delta transforms for three and more conductors.")
    Md = LinearAlgebra.diagm(0=>fill(1, n_phases), 1=>fill(-1, n_phases-1))
    Md[end,1] = -1
    return Md
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
        Memento.error(_LOGGER, "Zig-zag not yet supported.")
    end
    if t_type == "zig-zag"
        Memento.error(_LOGGER, "Zig-zag not yet supported.")
    end

    return tm_scale
end


"""
Returns bounds in line-to-line bounds on the voltage magnitude.
If these are not part of the problem specification, then a valid upper bound is
implied by the line-to-neutral bounds, but a lower bound (greater than zero) is
not. Therefore, a default lower bound is then used, specified by the keyword
argument vdmin_eps.
The returned bounds are for the pairs 1->2, 2->3, 3->1
"""
function _calc_bus_vm_ll_bounds(bus::Dict; vdmin_eps=0.1)
    vmax = bus["vmax"].values
    vmin = bus["vmin"].values
    if haskey(bus, "vm_ll_max")
        vdmax = bus["vm_ll_max"].values*sqrt(3)
    else
        # implied valid upper bound
        vdmax = [1 1 0; 0 1 1; 1 0 1]*vmax
        id = bus["index"]
    end
    if haskey(bus, "vm_ll_min")
        vdmin = bus["vm_ll_min"].values*sqrt(3)
    else
        vdmin = ones(3)*vdmin_eps*sqrt(3)
        id = bus["index"]
        Memento.info(_LOGGER, "Bus $id has no phase-to-phase vm upper bound; instead, $vdmin_eps was used as a valid upper bound.")
    end
    return (vdmin, vdmax)
end


"""
Calculates lower and upper bounds for the loads themselves (not the power
withdrawn at the bus).
"""
function _calc_load_pq_bounds(load::Dict, bus::Dict)
    a, alpha, b, beta = _load_expmodel_params(load, bus)
    vmin, vmax = _calc_load_vbounds(load, bus)
    # get bounds
    pmin = min.(a.*vmin.^alpha, a.*vmax.^alpha)
    pmax = max.(a.*vmin.^alpha, a.*vmax.^alpha)
    qmin = min.(b.*vmin.^beta, b.*vmax.^beta)
    qmax = max.(b.*vmin.^beta, b.*vmax.^beta)
    return (pmin, pmax, qmin, qmax)
end


"Returns a magnitude bound for the current going through the load."
function _calc_load_current_max(load::Dict, bus::Dict)
    pmin, pmax, qmin, qmax = _calc_load_pq_bounds(load, bus)
    pabsmax = max.(abs.(pmin), abs.(pmax))
    qabsmax = max.(abs.(qmin), abs.(qmax))
    smax = sqrt.(pabsmax.^2 + qabsmax.^2)

    vmin, vmax = _calc_load_vbounds(load, bus)

    return smax./vmin
end


"""
Returns magnitude bounds for the current going through the load.
"""
function _calc_load_current_magnitude_bounds(load::Dict, bus::Dict)
    a, alpha, b, beta = _load_expmodel_params(load, bus)
    vmin, vmax = _calc_load_vbounds(load, bus)
    cb1 = sqrt.(a.^(2).*vmin.^(2*alpha.-2) + b.^(2).*vmin.^(2*beta.-2))
    cb2 = sqrt.(a.^(2).*vmax.^(2*alpha.-2) + b.^(2).*vmax.^(2*beta.-2))
    cmin = min.(cb1, cb2)
    cmax = max.(cb1, cb2)
    return cmin, cmax
end


"""
Returns the exponential load model parameters for a load.
For an exponential load it simply returns certain data model properties, whilst
for constant_power, constant_current and constant_impedance it returns the
equivalent exponential model parameters.
"""
function _load_expmodel_params(load::Dict, bus::Dict)
    pd = load["pd"].values
    qd = load["qd"].values
    ncnds = length(pd)
    if load["model"]=="constant_power"
        return (pd, zeros(ncnds), qd, zeros(ncnds))
    else
        # get exponents
        if load["model"]=="constant_current"
            alpha = ones(ncnds)
            beta  =ones(ncnds)
        elseif load["model"]=="constant_impedance"
            alpha = ones(ncnds)*2
            beta  =ones(ncnds)*2
        elseif load["model"]=="exponential"
            alpha = load["alpha"].values
            @assert(all(alpha.>=0))
            beta = load["beta"].values
            @assert(all(beta.>=0))
        end
        # calculate proportionality constants
        v0 = load["vnom_kv"]/(bus["base_kv"]/sqrt(3))
        a = pd./v0.^alpha
        b = qd./v0.^beta
        # get bounds
        return (a, alpha, b, beta)
    end
end


"""
Returns the voltage magnitude bounds for the individual load elements in a
multiphase load. These are inferred from vmin/vmax for wye loads and from
_calc_bus_vm_ll_bounds for delta loads.
"""
function _calc_load_vbounds(load::Dict, bus::Dict)
    if load["conn"]=="wye"
        vmin = bus["vmin"].values
        vmax = bus["vmax"].values
    elseif load["conn"]=="delta"
        vmin, vmax = _calc_bus_vm_ll_bounds(bus)
    end
    return vmin, vmax
end

"""
Returns a Bool, indicating whether the convex hull of the voltage-dependent
relationship needs a cone inclusion constraint.
"""
function _check_load_needs_cone(load::Dict)
    if load["model"]=="constant_current"
        return true
    elseif load["model"]=="exponential"
        return true
    else
        return false
    end
end


"""
Returns a current magnitude bound for the generators.
"""
function _calc_gen_current_max(gen::Dict, bus::Dict)
    pabsmax = max.(abs.(gen["pmin"].values), abs.(gen["pmax"].values))
    qabsmax = max.(abs.(gen["qmax"].values), abs.(gen["qmax"].values))
    smax = sqrt.(pabsmax.^2 + qabsmax.^2)

    vmin = bus["vmin"].values

    return smax./vmin
end


"""
Returns a total (shunt+series) current magnitude bound for the from and to side
of a branch. The total power rating also implies a current bound through the
lower bound on the voltage magnitude of the connected buses.
"""
function _calc_branch_current_max_frto(branch::Dict, bus_fr::Dict, bus_to::Dict)
    bounds_fr = []
    bounds_to = []
    if haskey(branch, "c_rating_a")
        push!(bounds_fr, branch["c_rating_a"].values)
        push!(bounds_to, branch["c_rating_a"].values)
    end
    if haskey(branch, "rate_a")
        push!(bounds_fr, branch["rate_a"].values./bus_fr["vmin"].values)
        push!(bounds_to, branch["rate_a"].values./bus_to["vmin"].values)
    end
    @assert(length(bounds_fr)>=0, "no (implied/valid) current bounds defined")
    return min.(bounds_fr...), min.(bounds_to...)
end


"""
Returns a total (shunt+series) power magnitude bound for the from and to side
of a branch. The total current rating also implies a current bound through the
upper bound on the voltage magnitude of the connected buses.
"""
function _calc_branch_power_ub_frto(branch::Dict, bus_fr::Dict, bus_to::Dict)
    bounds_fr = []
    bounds_to = []
    if haskey(branch, "c_rating_a")
        push!(bounds_fr, branch["c_rating_a"].values.*bus_fr["vmax"].values)
        push!(bounds_to, branch["c_rating_a"].values.*bus_to["vmax"].values)
    end
    if haskey(branch, "rate_a")
        push!(bounds_fr, branch["rate_a"].values)
        push!(bounds_to, branch["rate_a"].values)
    end
    @assert(length(bounds_fr)>=0, "no (implied/valid) current bounds defined")
    return min.(bounds_fr...), min.(bounds_to...)
end


"""
Returns a valid series current magnitude bound for a branch.
"""
function _calc_branch_series_current_ub(branch::Dict, bus_fr::Dict, bus_to::Dict)
    vmin_fr = bus_fr["vmin"].values
    vmin_to = bus_to["vmin"].values

    vmax_fr = bus_fr["vmax"].values
    vmax_to = bus_to["vmax"].values

    # assumed to be matrices already
    # temportary fix by shunts_diag2mat!

    # get valid bounds on total current
    c_max_fr_tot, c_max_to_tot = _calc_branch_current_max_frto(branch, bus_fr, bus_to)

    # get valid bounds on shunt current
    y_fr = branch["g_fr"].values + im* branch["b_fr"].values
    y_to = branch["g_to"].values + im* branch["b_to"].values
    c_max_fr_sh = abs.(y_fr)*vmax_fr
    c_max_to_sh = abs.(y_to)*vmax_to

    # now select element-wise lowest valid bound between fr and to
    return min.(c_max_fr_sh.+c_max_fr_tot, c_max_to_sh.+c_max_to_tot)
end


# from PowerModels
"Transforms single-conductor network data into multi-conductor data"
function make_multiconductor!(data::Dict{String,<:Any}, conductors::Int)
    if InfrastructureModels.ismultinetwork(data)
        for (i,nw_data) in data["nw"]
            _make_multiconductor!(nw_data, conductors)
        end
    else
         _make_multiconductor!(data, conductors)
    end
end


"field names that should not be multi-conductor values"
const _conductorless = Set(["index", "bus_i", "bus_type", "status", "gen_status",
    "br_status", "gen_bus", "load_bus", "shunt_bus", "storage_bus", "f_bus", "t_bus",
    "transformer", "area", "zone", "base_kv", "energy", "energy_rating", "charge_rating",
    "discharge_rating", "charge_efficiency", "discharge_efficiency", "p_loss", "q_loss",
    "model", "ncost", "cost", "startup", "shutdown", "name", "source_id", "active_phases"])

"field names that should become multi-conductor matrix not arrays"
const _conductor_matrix = Set(["br_r", "br_x", "b_fr", "b_to", "g_fr", "g_to", "gs", "bs"])


""
function _make_multiconductor!(data::Dict{String,<:Any}, conductors::Real)
    if haskey(data, "conductors")
        Memento.warn(_LOGGER, "skipping network that is already multiconductor")
        return
    end

    data["conductors"] = conductors

    for (key, item) in data
        if isa(item, Dict{String,Any})
            for (item_id, item_data) in item
                if isa(item_data, Dict{String,Any})
                    item_ref_data = Dict{String,Any}()
                    for (param, value) in item_data
                        if param in _conductorless
                            item_ref_data[param] = value
                        else
                            if param in _conductor_matrix
                                item_ref_data[param] = MultiConductorMatrix(value, conductors)
                            else
                                item_ref_data[param] = MultiConductorVector(value, conductors)
                            end
                        end
                    end
                    item[item_id] = item_ref_data
                end
            end
        else
            #root non-dict items
        end
    end
end





""
function InfrastructureModels._value2string(v, float_precision::Int)
    if typeof(v) <: AbstractFloat
        return InfrastructureModels._float2string(v, float_precision)
    end
    if typeof(v) <: Array
        return "[$(join([InfrastructureModels._value2string(val, float_precision) for val in v], ", "))]"
    end
    if typeof(v) <: Dict
        return "{($(length(v)))}"
    end

    return "$(v)"
end
