import LinearAlgebra: Adjoint, pinv


"field names that should not be multi-conductor values"
const _conductorless = Set(["index", "bus_i", "bus_type", "status", "gen_status",
    "br_status", "gen_bus", "load_bus", "shunt_bus", "storage_bus", "f_bus", "t_bus",
    "transformer", "area", "zone", "base_kv", "energy", "energy_rating", "charge_rating",
    "discharge_rating", "charge_efficiency", "discharge_efficiency", "p_loss", "q_loss",
    "model", "ncost", "cost", "startup", "shutdown", "name", "source_id", "active_phases"])

"field names that should become multi-conductor matrix not arrays"
const _conductor_matrix = Set(["br_r", "br_x", "b_fr", "b_to", "g_fr", "g_to", "gs", "bs"])

const _excluded_count_busname_patterns = Vector{Regex}([
    r"^_virtual.*",
])

const _pmd_math_component_status_parameters = Set(["status", "gen_status", "br_status"])

"maps component types to status parameters"
const pmd_math_component_status = Dict(
    "bus" => "bus_type",
    "load" => "status",
    "shunt" => "status",
    "gen" => "gen_status",
    "storage" => "status",
    "switch" => "status",
    "branch" => "br_status",
    "transformer" => "status",
)

"maps component types to inactive status values"
const pmd_math_component_status_inactive = Dict(
    "bus" => 4,
    "load" => 0,
    "shunt" => 0,
    "gen" => 0,
    "storage" => 0,
    "switch" => 0,
    "branch" => 0,
    "transformer" => 0,
)

"""
    apply_pmd!(func!::Function, data::Dict{String, <:Any}; apply_to_subnetworks::Bool = true)

PowerModelsDistribution wrapper for the InfrastructureModels `apply!` function, working only on data
"""
function apply_pmd!(func!::Function, data::Dict{String, <:Any}; apply_to_subnetworks::Bool = true)
    _IM.apply!(func!, data, pmd_it_name; apply_to_subnetworks = apply_to_subnetworks)
end


"""
    get_pmd_data(data::Dict{String, <:Any})

Convenience function for retrieving the power-distribution-only portion of network data
"""
function get_pmd_data(data::Dict{String, <:Any})
    return _IM.ismultiinfrastructure(data) ? data["it"][pmd_it_name] : data
end


"BOUND manipulation methods (0*Inf->0 is often desired)"
_sum_rm_nan(X::Vector) = sum([X[(!).(isnan.(X))]..., 0.0])


""
function _mat_mult_rm_nan(A::Matrix, B::Union{Matrix, Adjoint}) where T
    N, A_ncols = size(A)
    B_nrows, M = size(B)
    @assert(A_ncols==B_nrows)
    return [_sum_rm_nan(A[n,:].*B[:,m]) for n in 1:N, m in 1:M]
end


_mat_mult_rm_nan(A::Union{Matrix, Adjoint}, b::Vector) = dropdims(_mat_mult_rm_nan(A, reshape(b, length(b), 1)), dims=2)
_mat_mult_rm_nan(a::Vector, B::Union{Matrix, Adjoint}) = _mat_mult_rm_nan(reshape(a, length(a), 1), B)


"Replaces NaN values with zeros"
function _nan2zero(b, a; val=0)
    and(x, y) = x && y
    b[and.(isnan.(b), a.==val)] .= 0.0
    return b
end


"Replaces NaN values with zeros"
_replace_nan(v) = map(x -> isnan(x) ? zero(x) : x, v)


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


"""
    count_nodes(data::Dict{String,<:Any})::Int

Counts number of nodes in network
"""
function count_nodes(data::Dict{String,<:Any})::Int
    n_nodes = 0

    if !ismissing(get(data, "data_model", missing)) && data["data_model"] == DSS
        all_nodes = Dict()
        for obj_type in values(data)
            if isa(obj_type, Dict)
                for object in values(obj_type)
                    if isa(object, Dict)
                        if haskey(object, "buses")
                            for busname in values(object["buses"])
                                name, nodes = _parse_bus_id(busname)

                                if !haskey(all_nodes, name)
                                    all_nodes[name] = Set([])
                                end

                                for (n, node) in enumerate(nodes[1:3])
                                    if node
                                        push!(all_nodes[name], n)
                                    end
                                end
                            end
                        else
                            for (prop, val) in object
                                if startswith(prop, "bus") && prop != "buses"
                                    name, nodes = _parse_bus_id(val)

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
            end
        end

        for (name, phases) in all_nodes
            n_nodes += length(phases)
        end
    elseif !ismissing(get(data, "data_model", missing)) && data["data_model"] in [MATHEMATICAL, ENGINEERING] || (haskey(data, "source_type") && data["source_type"] == "matlab")
        n_nodes = 0
        for (name, bus) in data["bus"]
            if get(data, "data_model", missing) == MATPOWER || (haskey(data, "source_type") && data["source_type"] == "matlab")
                n_nodes += sum(bus["vm"] .> 0.0)
            else
                if data["data_model"] == MATHEMATICAL
                    name = bus["name"]
                end

                if all(!occursin(pattern, name) for pattern in [_excluded_count_busname_patterns...])
                    if data["data_model"] == MATHEMATICAL
                        if get(data, "is_projected", false)
                            n_nodes += count(i->i>0, get(bus, "vmax", []))
                        else
                            n_nodes += length(bus["terminals"][.!get(bus, "grounded", zeros(length(bus["terminals"])))])
                        end
                    else
                        n_nodes += length([n for n in bus["terminals"] if !(n in get(bus, "grounded", []))])
                    end
                end
            end
        end
    else
        error("Origin of data structure not recognized, cannot count nodes reliably")
    end

    return n_nodes
end

"""
    count_active_connections(data::Dict{String,<:Any})

Counts active ungrounded connections on edge components
"""
function count_active_connections(data::Dict{String,<:Any})
    data_model = get(data, "data_model", MATHEMATICAL)
    edge_elements = data_model == MATHEMATICAL ? PowerModelsDistribution._math_edge_elements : PowerModelsDistribution._eng_edge_elements
    # bus_connections = Dict(id => [] for (id, _) in data["bus"])
    active_connections = 0

    for edge_type in edge_elements
        for (_, component) in get(data, edge_type, Dict())
            counted_connections = Set([])
            if edge_type == "transformer" && !haskey(component, "f_connections") && data_model == ENGINEERING
                for (wdg, connections) in enumerate(component["connections"])
                    for terminal in connections
                        if !(terminal in counted_connections)
                            if !(terminal in data["bus"][component["bus"][wdg]]["grounded"])
                                push!(counted_connections, terminal)
                                active_connections += 1
                            end
                        end
                    end
                end
            else
                for (bus, connections) in [(component["f_bus"], component["f_connections"]), (component["t_bus"], component["t_connections"])]
                    for (i, terminal) in enumerate(connections)
                        if !(terminal in counted_connections)
                            if data_model == ENGINEERING
                                if edge_type == "transformer" && (get(data, "is_kron_reduced", false) || (component["configuration"] == WYE && terminal != connections[end]))
                                    push!(counted_connections, terminal)
                                    active_connections += 1
                                elseif !(terminal in data["bus"][bus]["grounded"])
                                    push!(counted_connections, terminal)
                                    active_connections += 1
                                end
                            else
                                if edge_type == "transformer"
                                    if get(data, "is_kron_reduced", false) || component["configuration"] == DELTA || (component["configuration"] == WYE && terminal != connections[end])
                                        push!(counted_connections, terminal)
                                        active_connections += 1
                                    end
                                elseif !get(data["bus"]["$bus"]["grounded"], i, false)
                                    if get(data, "is_projected", false) && get(data["bus"]["$bus"]["vmax"], i, Inf) > 0
                                        push!(counted_connections, terminal)
                                        active_connections += 1
                                    else
                                        push!(counted_connections, terminal)
                                        active_connections += 1
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return active_connections
end


"""
    count_active_terminals(data::Dict{String,<:Any}; count_grounded::Bool=false)

Counts active ungrounded terminals on buses
"""
function count_active_terminals(data::Dict{String,<:Any}; count_grounded::Bool=false)
    data_model = get(data, "data_model", MATHEMATICAL)
    active_terminal_count = 0
    for (_,bus) in data["bus"]
        counted_terminals = []
        for (i, terminal) in enumerate(bus["terminals"])
            if !(terminal in counted_terminals)
                if count_grounded
                    push!(counted_terminals, terminal)
                    active_terminal_count += 1
                else
                    if data_model == ENGINEERING
                        if !(terminal in bus["grounded"])
                            push!(counted_terminals, terminal)
                            active_terminal_count += 1
                        end
                    else
                        if !get(bus["grounded"], i, false)
                            push!(counted_terminals, terminal)
                            active_terminal_count += 1
                        end
                    end
                end
            end
        end
    end
    return active_terminal_count
end


"returns the conductor indicator for a ENGINEERING component"
function _get_conductor_indicator(comp::Dict{String,<:Any})::String
    if haskey(comp, "terminals")
        return "terminals"
    elseif haskey(comp, "connections")
        return "connections"
    elseif haskey(comp, "f_connections")
        return "f_connections"
    else
        return ""
    end
end


"creates a delta transformation matrix"
function _get_delta_transformation_matrix(n_phases::Int)
    @assert(n_phases>2, "We only define delta transforms for three and more conductors.")
    Md = LinearAlgebra.diagm(0=>fill(1, n_phases), 1=>fill(-1, n_phases-1))
    Md[end,1] = -1
    return Md
end


"Calculates the tap scale factor for the non-dimensionalized equations."
function calculate_tm_scale(trans::Dict{String,Any}, bus_fr::Dict{String,Any}, bus_to::Dict{String,Any})
    tm_nom = trans["tm_nom"]

    f_vbase = haskey(bus_fr, "vbase") ? bus_fr["vbase"] : bus_fr["base_kv"]
    t_vbase = haskey(bus_to, "vbase") ? bus_to["vbase"] : bus_to["base_kv"]
    config = trans["configuration"]

    tm_scale = tm_nom*(t_vbase/f_vbase)
    if config == DELTA
        #TODO is this still needed?
        tm_scale *= sqrt(3)
    elseif config == "zig-zag"
        error("Zig-zag not yet supported.")
    end

    return tm_nom
end


"""
Returns a power magnitude bound for the from and to side of a transformer.
The total current rating also implies a current bound through the
upper bound on the voltage magnitude of the connected buses.
"""
function _calc_transformer_power_ub_frto(trans::Dict{String,<:Any}, bus_fr::Dict{String,<:Any}, bus_to::Dict{String,<:Any})
    bounds_fr = []
    bounds_to = []
    if haskey(trans, "c_rating_a")
        push!(bounds_fr, trans["c_rating_a"].*bus_fr["vmax"].*bus_fr["vbase"])
        push!(bounds_to, trans["c_rating_a"].*bus_to["vmax"].*bus_to["vbase"])
    end
    if haskey(trans, "rate_a")
        push!(bounds_fr, trans["rate_a"])
        push!(bounds_to, trans["rate_a"])
    end

    N = length(trans["f_connections"])
    return min.(fill(Inf, N), bounds_fr...), min.(fill(Inf, N), bounds_to...)
end


"""
Returns a current magnitude bound for the from and to side of a transformer.
The total power rating also implies a current bound through the lower bound on
the voltage magnitude of the connected buses.
"""
function _calc_transformer_current_max_frto(trans::Dict{String,<:Any}, bus_fr::Dict{String,<:Any}, bus_to::Dict{String,<:Any})
    bounds_fr = []
    bounds_to = []
    if haskey(trans, "c_rating_a")
        push!(bounds_fr, trans["c_rating_a"])
        push!(bounds_to, trans["c_rating_a"])
    end
    if haskey(trans, "rate_a")
        push!(bounds_fr, trans["rate_a"]./(bus_fr["vmax"].*bus_fr["vbase"]))
        push!(bounds_to, trans["rate_a"]./(bus_to["vmax"].*bus_to["vbase"]))
    end

    N = length(trans["f_connections"])
    return min.(fill(Inf, N), bounds_fr...), min.(fill(Inf, N), bounds_to...)
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
    vmax = bus["vmax"]
    vmin = bus["vmin"]
    if haskey(bus, "vm_ll_max")
        vdmax = bus["vm_ll_max"]
    else
        # implied valid upper bound
        vdmax = _mat_mult_rm_nan([1 1 0; 0 1 1; 1 0 1], vmax)
    end
    if haskey(bus, "vm_ll_min")
        vdmin = bus["vm_ll_min"]
    else
        vdmin = fill(0.0, length(vmin))
    end

    return (vdmin, vdmax)
end


"""
Calculates lower and upper bounds for the loads themselves (not the power
withdrawn at the bus).
"""
function _calc_load_pq_bounds(load::Dict{String,<:Any}, bus::Dict{String,<:Any})
    a, alpha, b, beta = _load_expmodel_params(load, bus)
    vmin, vmax = _calc_load_vbounds(load, bus)
    # get bounds
    pmin = _nan2zero(min.(a.*vmin.^alpha, a.*vmax.^alpha), a)
    pmax = _nan2zero(max.(a.*vmin.^alpha, a.*vmax.^alpha), a)
    qmin = _nan2zero(min.(b.*vmin.^beta, b.*vmax.^beta), b)
    qmax = _nan2zero(max.(b.*vmin.^beta, b.*vmax.^beta), b)
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
    cb1 = sqrt.(_nan2zero(a.^(2).*vmin.^(2*alpha.-2), a) + _nan2zero(b.^(2).*vmin.^(2*beta.-2), b))
    cb2 = sqrt.(_nan2zero(a.^(2).*vmax.^(2*alpha.-2), a) + _nan2zero(b.^(2).*vmax.^(2*beta.-2), b))
    cmin = min.(cb1, cb2)
    cmax = max.(cb1, cb2)
    return cmin, cmax
end


"""
    _load_expmodel_params(load::Dict{String,<:Any}, bus::Dict{String,<:Any})

Returns the exponential load model parameters for a load.
For an exponential load it simply returns certain data model properties, whilst
for constant_power, constant_current and constant_impedance it returns the
equivalent exponential model parameters.
"""
function _load_expmodel_params(load::Dict, bus::Dict)
    pd = load["pd"]
    qd = load["qd"]
    ncnds = length(pd)
    if load["model"]==POWER
        return (pd, zeros(ncnds), qd, zeros(ncnds))
    else
        # get exponents
        if load["model"]==CURRENT
            alpha = ones(ncnds)
            beta  =ones(ncnds)
        elseif load["model"]==IMPEDANCE
            alpha = ones(ncnds)*2
            beta  =ones(ncnds)*2
        elseif load["model"]==EXPONENTIAL
            alpha = load["alpha"]
            @assert(all(alpha.>=0))
            beta = load["beta"]
            @assert(all(beta.>=0))
        else
            error("load model '$(load["model"])' on 'load.$(load["name"])' unrecongized")
        end
        # calculate proportionality constants
        v0 = load["vnom_kv"]
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
function _calc_load_vbounds(load::Dict{String,<:Any}, bus::Dict{String,<:Any})
    terminals = bus["terminals"]
    connections = [findfirst(isequal(c), terminals) for c in load["connections"]]

    if load["configuration"]==WYE
        vmin = bus["vmin"]
        vmax = bus["vmax"]
    elseif load["configuration"]==DELTA
        vmin, vmax = _calc_bus_vm_ll_bounds(bus)
    end
    return vmin[connections], vmax[connections]
end


"""
Returns a Bool, indicating whether the convex hull of the voltage-dependent
relationship needs a cone inclusion constraint.
"""
function _check_load_needs_cone(load::Dict)
    if load["model"]==CURRENT
        return true
    elseif load["model"]==EXPONENTIAL
        return true
    else
        return false
    end
end


"""
Returns a current magnitude bound for the generators.
"""
function _calc_gen_current_max(gen::Dict{String,<:Any}, bus::Dict{String,<:Any})
    if all([haskey(gen, prop) for prop in ["pmax", "pmin", "qmax", "qmin"]]) && haskey(bus, "vmin")
        pabsmax = max.(abs.(gen["pmin"]), abs.(gen["pmax"]))
        qabsmax = max.(abs.(gen["qmin"]), abs.(gen["qmax"]))
        smax = sqrt.(pabsmax.^2 + qabsmax.^2)

        vmin = bus["vmin"][[findfirst(isequal(c), bus["terminals"]) for c in gen["connections"]]]

        return smax./vmin
    else
        return fill(Inf, length(gen["connections"]))
    end
end


"""
Returns a total (shunt+series) current magnitude bound for the from and to side
of a branch. The total power rating also implies a current bound through the
lower bound on the voltage magnitude of the connected buses.
"""
function _calc_branch_current_max(branch::Dict{String,<:Any}, bus::Dict{String,<:Any})
    bounds = []

    if haskey(branch, "c_rating_a")
        push!(bounds, branch["c_rating_a"])
    end
    if haskey(branch, "rate_a") && haskey(bus, "vmin")
        push!(bounds, branch["rate_a"]./bus["vmin"][[findfirst(isequal(c), bus["terminals"]) for c in branch["f_connections"]]])
    end

    return min.(fill(Inf, length(branch["f_connections"])), bounds...)
end


"""
Returns a total (shunt+series) current magnitude bound for the from and to side
of a branch. The total power rating also implies a current bound through the
lower bound on the voltage magnitude of the connected buses.
"""
function _calc_branch_current_max_frto(branch::Dict{String,<:Any}, bus_fr::Dict{String,<:Any}, bus_to::Dict{String,<:Any})
    bounds_fr = []
    bounds_to = []

    if haskey(branch, "c_rating_a")
        push!(bounds_fr, branch["c_rating_a"])
        push!(bounds_to, branch["c_rating_a"])
    end
    if haskey(branch, "rate_a")
        push!(bounds_fr, branch["rate_a"]./bus_fr["vmin"][[findfirst(isequal(c), bus_fr["terminals"]) for c in branch["f_connections"]]])
        push!(bounds_to, branch["rate_a"]./bus_to["vmin"][[findfirst(isequal(c), bus_to["terminals"]) for c in branch["t_connections"]]])
    end

    return min.(fill(Inf, length(branch["f_connections"])), bounds_fr...), min.(fill(Inf, length(branch["t_connections"])), bounds_to...)
end


"""
Returns a total (shunt+series) power magnitude bound for the from and to side
of a branch. The total current rating also implies a current bound through the
upper bound on the voltage magnitude of the connected buses.
"""
function _calc_branch_power_max(branch::Dict{String,<:Any}, bus::Dict{String,<:Any})
    bounds = []

    terminals = bus["terminals"]
    connections = bus["bus_i"] == branch["f_bus"] ? branch["f_connections"] : branch["t_connections"]
    connections = [findfirst(isequal(cnd), terminals) for cnd in connections]

    if haskey(bus, "vmax") && (haskey(branch, "c_rating_b") || haskey(branch, "c_rating_a"))
        c_rating = haskey(branch, "c_rating_b") ? branch["c_rating_b"] : branch["c_rating_a"]
        push!(bounds, c_rating .* bus["vmax"][connections] .* bus["vbase"])
    end
    if haskey(branch, "rate_a")
        push!(bounds, branch["rate_a"])
    end

    N = length(connections)
    return min.(fill(Inf, N), bounds...)
end


"""
Returns a total (shunt+series) power magnitude bound for the from and to side
of a branch. The total current rating also implies a current bound through the
upper bound on the voltage magnitude of the connected buses.
"""
function _calc_branch_power_max_frto(branch::Dict{String,<:Any}, bus_fr::Dict{String,<:Any}, bus_to::Dict{String,<:Any})
    return _calc_branch_power_max(branch, bus_fr), _calc_branch_power_max(branch, bus_to)
end


"""
Returns a valid series current magnitude bound for a branch.
"""
function _calc_branch_series_current_max(branch::Dict{String,<:Any}, bus_fr::Dict{String,<:Any}, bus_to::Dict{String,<:Any})
    ncnds = length(branch["f_connections"])
    vmin_fr = haskey(bus_fr, "vmin") ? bus_fr["vmin"][[findfirst(isequal(c), bus_fr["terminals"]) for c in branch["f_connections"]]] : fill(0.0, ncnds)
    vmin_to = haskey(bus_to, "vmin") ? bus_fr["vmin"][[findfirst(isequal(c), bus_to["terminals"]) for c in branch["t_connections"]]] : fill(0.0, ncnds)

    vmax_fr = haskey(bus_fr, "vmax") ? bus_fr["vmax"][[findfirst(isequal(c), bus_fr["terminals"]) for c in branch["f_connections"]]] : fill(Inf, ncnds)
    vmax_to = haskey(bus_to, "vmax") ? bus_to["vmax"][[findfirst(isequal(c), bus_to["terminals"]) for c in branch["t_connections"]]] : fill(Inf, ncnds)

    # assumed to be matrices already
    # temportary fix by shunts_diag2mat!

    # get valid bounds on total current
    c_max_fr_tot = _calc_branch_current_max(branch, bus_fr)
    c_max_to_tot = _calc_branch_current_max(branch, bus_to)

    # get valid bounds on shunt current
    y_fr = branch["g_fr"] + im* branch["b_fr"]
    y_to = branch["g_to"] + im* branch["b_to"]
    c_max_fr_sh = _mat_mult_rm_nan(abs.(y_fr), vmax_fr)
    c_max_to_sh = _mat_mult_rm_nan(abs.(y_to), vmax_to)

    # now select element-wise lowest valid bound between fr and to
    N = length(branch["f_connections"])
    return min.(fill(Inf, N), c_max_fr_sh.+c_max_fr_tot, c_max_to_sh.+c_max_to_tot)
end


"https://stackoverflow.com/questions/39039553/lower-triangular-matrix-in-julia"
function _vec2utri!(v::Vector{T}) where T
    d = length(v)
    n = Int((sqrt(8d+1)+1)/2)
    n*(n-1)/2 == d || error("vec2utri: length of vector is not triangular")
    [ i<j ? v[Int((j-1)*(j-2)/2)+i] : 0 for i=1:n, j=1:n ]
end


"vector to lower triangular"
function _vec2ltri!(v::Vector{T}) where T
    _vec2utri!(v)'
end


"matrix to upper triangular vector"
function _mat2utrivec!(m::Union{Matrix{T}, LinearAlgebra.Symmetric{T}}) where T
    @assert size(m,1) == size(m,2)
    n = size(m,1)
    [m[i,j] for i=1:n, j=1:n if i < j]
end


"matrix to lower triangular vector"
function _mat2ltrivec!(m::Union{Matrix{T}, LinearAlgebra.Symmetric{T}}) where T
    @assert size(m,1) == size(m,2)
    n = size(m,1)
    [m[j,i] for i=1:n, j=1:n if i < j]
end


"makes a hermitian matrix variable from diagonal, and lower and upper triangular vectors"
function _make_hermitian_matrix_variable(diag, lowertrianglereal, lowertriangleimag)
    #TODO clean up
    matrixreal = []
    if length(diag) == 3
        matrixreal = [
        diag[1]                 lowertrianglereal[1]    lowertrianglereal[2];
        lowertrianglereal[1]    diag[2]                 lowertrianglereal[3];
        lowertrianglereal[2]    lowertrianglereal[3]    diag[3]
        ]
    elseif length(diag) == 4
        matrixreal = [
        diag[1]                 lowertrianglereal[1]    lowertrianglereal[2]    lowertrianglereal[4];
        lowertrianglereal[1]    diag[2]                 lowertrianglereal[3]    lowertrianglereal[5];
        lowertrianglereal[2]    lowertrianglereal[3]    diag[3]                 lowertrianglereal[6];
        lowertrianglereal[4]    lowertrianglereal[5]    lowertrianglereal[6]    diag[4]
        ]
    elseif length(diag) == 5
        matrixreal = [
        diag[1]                 lowertrianglereal[1]    lowertrianglereal[2]    lowertrianglereal[4]    lowertrianglereal[7];
        lowertrianglereal[1]    diag[2]                 lowertrianglereal[3]    lowertrianglereal[5]    lowertrianglereal[8];
        lowertrianglereal[2]    lowertrianglereal[3]    diag[3]                 lowertrianglereal[6]    lowertrianglereal[9];
        lowertrianglereal[4]    lowertrianglereal[5]    lowertrianglereal[6]    diag[4]                 lowertrianglereal[10];
        lowertrianglereal[7]    lowertrianglereal[8]    lowertrianglereal[9]    lowertrianglereal[10]    diag[5]
        ]
    end

    matriximag = []
    if length(diag) == 3
        matriximag = [
        0                       -lowertriangleimag[1]   -lowertriangleimag[2];
        lowertriangleimag[1]    0                       -lowertriangleimag[3];
        lowertriangleimag[2]    lowertriangleimag[3]    0
        ]
    elseif length(diag) == 4
        matriximag = [
        0                       -lowertriangleimag[1]   -lowertriangleimag[2]   -lowertriangleimag[4];
        lowertriangleimag[1]    0                       -lowertriangleimag[3]   -lowertriangleimag[5];
        lowertriangleimag[2]    lowertriangleimag[3]    0                       -lowertriangleimag[6];
        lowertriangleimag[4]    lowertriangleimag[5]    lowertriangleimag[6]    0
        ]
    elseif length(diag) == 5
        matriximag = [
        0                       -lowertriangleimag[1]   -lowertriangleimag[2]   -lowertriangleimag[4]   -lowertriangleimag[7];
        lowertriangleimag[1]    0                       -lowertriangleimag[3]   -lowertriangleimag[5]   -lowertriangleimag[8];
        lowertriangleimag[2]    lowertriangleimag[3]    0                       -lowertriangleimag[6]   -lowertriangleimag[9];
        lowertriangleimag[4]    lowertriangleimag[5]    lowertriangleimag[6]    0                       -lowertriangleimag[10];
        lowertriangleimag[7]    lowertriangleimag[8]    lowertriangleimag[9]    lowertriangleimag[10]    0
        ]
    end
    return matrixreal, matriximag
end


"makes a full matrix variable from a diagonal, and lower and upper triangular vectors"
function _make_full_matrix_variable(diag, lowertriangle, uppertriangle)
    #TODO clean up
    matrix = []
    if length(diag) == 3
        matrix = [
        diag[1]             uppertriangle[1]    uppertriangle[2];
        lowertriangle[1]    diag[2]             uppertriangle[3];
        lowertriangle[2]    lowertriangle[3]    diag[3]
        ]
    elseif length(diag) == 4
        matrix = [
        diag[1]             uppertriangle[1]    uppertriangle[2]    uppertriangle[4];
        lowertriangle[1]    diag[2]             uppertriangle[3]    uppertriangle[5];
        lowertriangle[2]    lowertriangle[3]    diag[3]             uppertriangle[6];
        lowertriangle[4]    lowertriangle[5]    lowertriangle[6]    diag[4]
        ]
    elseif length(diag) == 5
        matrix = [
        diag[1]             uppertriangle[1]    uppertriangle[2]    uppertriangle[4]    uppertriangle[7];
        lowertriangle[1]    diag[2]             uppertriangle[3]    uppertriangle[5]    uppertriangle[8];
        lowertriangle[2]    lowertriangle[3]    diag[3]             uppertriangle[6]    uppertriangle[9];
        lowertriangle[4]    lowertriangle[5]    lowertriangle[6]    diag[4]             uppertriangle[10]
        lowertriangle[7]    lowertriangle[8]    lowertriangle[9]    lowertriangle[10]    diag[5]
        ]
    end
    # matrix = diagm(0 => diag) + _vec2ltri!(lowertriangle) + _vec2utri!(uppertriangle)
    return matrix
end


"helper to determine if expession has any Nonlinear terms"
function _has_nl_expression(x)::Bool
    if isa(x, JuMP.NonlinearExpression)
        return true
    elseif isa(x, Vector)
        for i in x
            if _has_nl_expression(i)
                return true
            end
        end
    elseif isa(x, Dict)
        for i in values(x)
            if _has_nl_expression(i)
                return true
            end
        end
    elseif isa(x, JuMP.Containers.DenseAxisArray)
        for i in values(x.data)
            if _has_nl_expression(i)
                return true
            end
        end
    end
    return false
end


"""
    correct_mc_voltage_angle_differences!(data::Dict{String,<:Any}, default_pad::Real=deg2rad(10.0))

checks that voltage angle differences are within 90 deg., if not tightens to a default of 10deg (adjustable)
"""
function correct_mc_voltage_angle_differences!(data::Dict{String,<:Any}, default_pad::Real=deg2rad(10.0))
    if ismultinetwork(data)
        nw_data = data["nw"]
    else
        nw_data = Dict("0" => data)
    end

    @assert("per_unit" in keys(data) && data["per_unit"])
    default_pad_deg = round(rad2deg(default_pad), digits=2)

    modified = Set{Int}()

    for (n,nw) in nw_data
        if ismultinetwork(data)
            subnet = ", subnetwork $n"
        else
            subnet = ""
        end

        for (i, branch) in nw["branch"]
            angmin = branch["angmin"]
            angmax = branch["angmax"]

            if any(angmin .<= -pi/2)
                @warn "this code only supports angmin values in -90 deg. to 90 deg., tightening the value on branch $(i)$(subnet) from $(rad2deg(angmin)) to -$(default_pad_deg) deg."
                branch["angmin"][angmin .<= -pi/2] .= -default_pad

                push!(modified, branch["index"])
            end

            if any(angmax .>= pi/2)
                @warn "this code only supports angmax values in -90 deg. to 90 deg., tightening the value on branch $(i)$(subnet) from $(rad2deg(angmax)) to $(default_pad_deg) deg."
                branch["angmax"][angmax .>= pi/2] .= default_pad

                push!(modified, branch["index"])
            end

            if any((angmin .== 0.0) .& (angmax .== 0.0))
                @warn "angmin and angmax values are 0, widening these values on branch $(i)$(subnet) to +/- $(default_pad_deg) deg."
                branch["angmin"][(angmin .== 0.0) .& (angmax .== 0.0)] .= -default_pad
                branch["angmax"][(angmin .== 0.0) .& (angmax .== 0.0)] .=  default_pad

                push!(modified, branch["index"])
            end
        end
    end

    return modified
end

"""
    correct_mc_thermal_limits!(data::Dict{String,<:Any})

checks that each branch has non-negative thermal ratings and removes zero thermal ratings
"""
function correct_mc_thermal_limits!(data::Dict{String,<:Any})
    if ismultinetwork(data)
        nw_data = data["nw"]
    else
        nw_data = Dict("0" => data)
    end

    modified = Set{Int}()

    for (n,nw) in nw_data
        if ismultinetwork(data)
            subnet = ", subnetwork $n"
        else
            subnet = ""
        end

        branches = [branch for branch in values(nw["branch"])]
        if haskey(data, "ne_branch")
            append!(branches, values(data["ne_branch"]))
        end

        for branch in branches
            for rate_key in ["rate_a", "rate_b", "rate_c"]
                if haskey(branch, rate_key)
                    rate_value = branch[rate_key]

                    if any(rate_value .< 0.0)
                        error("negative $(rate_key) value on branch $(branch["index"])$(subnet), this code only supports non-negative $(rate_key) values")
                    end

                    if all(isapprox.(rate_value, 0.0))
                        delete!(branch, rate_key)
                        @warn "removing zero $(rate_key) limit on branch $(branch["index"])$(subnet)"

                        push!(modified, branch["index"])
                    end
                end
            end
        end
    end

    return modified
end


"helper function to build bus shunt matrices for power balance constraints"
function _build_bus_shunt_matrices(pm::AbstractUnbalancedPowerModel, nw::Int, terminals::Vector{Int}, bus_shunts::Vector{<:Tuple{Int,Vector{Int}}})::Tuple{Matrix{<:Real},Matrix{<:Real}}
    ncnds = length(terminals)
    Gs = fill(0.0, ncnds, ncnds)
    Bs = fill(0.0, ncnds, ncnds)
    for (i, connections) in bus_shunts
        shunt = ref(pm,nw,:shunt,i)
        for (idx,c) in enumerate(connections)
            for (jdx,d) in enumerate(connections)
                Gs[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] += shunt["gs"][idx,jdx]
                Bs[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] += shunt["bs"][idx,jdx]
            end
        end
    end

    return (Gs, Bs)
end


"""
    calc_branch_y(branch::Dict{String,<:Any})

computes branch admittance matrices
"""
function calc_branch_y(branch::Dict{String,<:Any})
    y = pinv(branch["br_r"] + im * branch["br_x"])
    g, b = real(y), imag(y)
    return g, b
end


"""
    identify_load_blocks(data::Dict{String,<:Any})

computes load blocks based on switch locations
"""
identify_load_blocks(data::Dict{String,<:Any}) = calc_connected_components(data; type="load_blocks")


"""
    identify_blocks(data::Dict{String,<:Any})

computes connected blocks currently in the model based on switch states
"""
identify_blocks(data::Dict{String,<:Any}) = calc_connected_components(data; type="blocks")


"""
    identify_islands(data::Dict{String,<:Any})

computes component islands base only on edge and bus status
"""
identify_islands(data::Dict{String,<:Any}) = calc_connected_components(data)


"""
    calc_connected_components(data::Dict{String,<:Any}; edges::Union{Missing, Vector{<:String}}=missing, type::Union{Missing,String}=missing, check_enabled::Bool=true)::Set

computes the connected components of the network graph
returns a set of sets of bus ids, each set is a connected component
"""
function calc_connected_components(data::Dict{String,<:Any}; edges::Union{Missing, Vector{<:String}}=missing, type::Union{Missing,String}=missing, check_enabled::Bool=true)::Set
    pmd_data = get_pmd_data(data)

    if ismultinetwork(pmd_data)
        error("multinetwork data is not yet supported, recommend to use on each subnetwork independently")
    end

    if get(pmd_data, "data_model", MATHEMATICAL) == ENGINEERING
        return _calc_connected_components_eng(pmd_data; edges=ismissing(edges) ? _eng_edge_elements : edges, type=type, check_enabled=check_enabled)
    elseif get(pmd_data, "data_model", MATHEMATICAL) == MATHEMATICAL
        return _calc_connected_components_math(pmd_data; edges=ismissing(edges) ? _math_edge_elements : edges, type=type, check_enabled=check_enabled)
    else
        error("data_model `$(get(pmd_data, "data_model", MATHEMATICAL))` is unrecongized")
    end
end


"""
computes the connected components of the network graph
returns a set of sets of bus ids, each set is a connected component
"""
function _calc_connected_components_eng(data; edges::Vector{<:String}=_eng_edge_elements, type::Union{Missing,String}=missing, check_enabled::Bool=true)::Set
    @assert get(data, "data_model", MATHEMATICAL) == ENGINEERING

    active_bus = Dict{Any,Dict{String,Any}}(x for x in data["bus"] if x.second["status"] == ENABLED || !check_enabled)
    active_bus_ids = Set{Any}([i for (i,bus) in active_bus])

    neighbors = Dict{Any,Vector{Any}}(i => [] for i in active_bus_ids)
    for edge_type in edges
        for (id, edge_obj) in get(data, edge_type, Dict{Any,Dict{String,Any}}())
            if edge_obj["status"] == ENABLED || !check_enabled
                if edge_type == "transformer" && haskey(edge_obj, "bus")
                    for f_bus in edge_obj["bus"]
                        for t_bus in edge_obj["bus"]
                            if f_bus != t_bus
                                push!(neighbors[f_bus], t_bus)
                                push!(neighbors[t_bus], f_bus)
                            end
                        end
                    end
                else
                    if edge_type == "switch" && !ismissing(type)
                        if type == "load_blocks"
                            if edge_obj["dispatchable"] == NO && edge_obj["state"] == CLOSED
                                push!(neighbors[edge_obj["f_bus"]], edge_obj["t_bus"])
                                push!(neighbors[edge_obj["t_bus"]], edge_obj["f_bus"])
                            end
                        elseif type == "blocks"
                            if edge_obj["state"] == CLOSED
                                push!(neighbors[edge_obj["f_bus"]], edge_obj["t_bus"])
                                push!(neighbors[edge_obj["t_bus"]], edge_obj["f_bus"])
                            end
                        end
                    else
                        push!(neighbors[edge_obj["f_bus"]], edge_obj["t_bus"])
                        push!(neighbors[edge_obj["t_bus"]], edge_obj["f_bus"])
                    end
                end
            end
        end
    end

    component_lookup = Dict(i => Set{String}([i]) for i in active_bus_ids)
    touched = Set{String}()

    for i in active_bus_ids
        if !(i in touched)
            _cc_dfs(i, neighbors, component_lookup, touched)
        end
    end

    ccs = (Set(values(component_lookup)))

    return ccs
end


"""
computes the connected components of the network graph
returns a set of sets of bus ids, each set is a connected component
"""
function _calc_connected_components_math(data::Dict{String,<:Any}; edges::Vector{<:String}=_math_edge_elements, type::Union{Missing,String}=missing, check_enabled::Bool=true)::Set
    @assert get(data, "data_model", MATHEMATICAL) == MATHEMATICAL

    active_bus = Dict{Any,Dict{String,Any}}(x for x in data["bus"] if x.second[pmd_math_component_status["bus"]] != pmd_math_component_status_inactive["bus"] || !check_enabled)
    active_bus_ids = Set{Any}([parse(Int,i) for (i,bus) in active_bus])

    neighbors = Dict{Any,Vector{Any}}(i => [] for i in active_bus_ids)
    for edge_type in edges
        for (id, edge_obj) in get(data, edge_type, Dict{Any,Dict{String,Any}}())
            if edge_obj[pmd_math_component_status[edge_type]] != pmd_math_component_status_inactive[edge_type] || !check_enabled
                if edge_type == "switch" && !ismissing(type)
                    if type == "load_blocks"
                        if edge_obj["dispatchable"] != 1 && edge_obj["state"] == 1
                            push!(neighbors[edge_obj["f_bus"]], edge_obj["t_bus"])
                            push!(neighbors[edge_obj["t_bus"]], edge_obj["f_bus"])
                        end
                    elseif type == "blocks"
                        if edge_obj["state"] != 0
                            push!(neighbors[edge_obj["f_bus"]], edge_obj["t_bus"])
                            push!(neighbors[edge_obj["t_bus"]], edge_obj["f_bus"])
                        end
                    end
                else
                    push!(neighbors[edge_obj["f_bus"]], edge_obj["t_bus"])
                    push!(neighbors[edge_obj["t_bus"]], edge_obj["f_bus"])
                end
            end
        end
    end

    component_lookup = Dict(i => Set{Int}([i]) for i in active_bus_ids)
    touched = Set{Int64}()

    for i in active_bus_ids
        if !(i in touched)
            _cc_dfs(i, neighbors, component_lookup, touched)
        end
    end

    ccs = (Set(values(component_lookup)))

    return ccs
end


"DFS on a graph"
function _cc_dfs(i, neighbors, component_lookup, touched)
    push!(touched, i)
    for j in neighbors[i]
        if !(j in touched)
            for k in  component_lookup[j]
                push!(component_lookup[i], k)
            end
            for k in component_lookup[j]
                component_lookup[k] = component_lookup[i]
            end
            _cc_dfs(j, neighbors, component_lookup, touched)
        end
    end
end


"rescales the cost model terms"
function _rescale_cost_model!(comp::Dict{String,<:Any}, scale::Real)
    if "model" in keys(comp) && "cost" in keys(comp)
        if comp["model"] == 1
            for i in 1:2:length(comp["cost"])
                comp["cost"][i] = comp["cost"][i]/scale
            end
        elseif comp["model"] == 2
            degree = length(comp["cost"])
            for (i, item) in enumerate(comp["cost"])
                comp["cost"][i] = item*(scale^(degree-i))
            end
        else
            @warn "Skipping cost model of type $(comp["model"]) in per unit transformation"
        end
    end
end


"""
    correct_branch_directions!(data::Dict{String,<:Any})

checks that all parallel branches have the same orientation
"""
function correct_branch_directions!(data::Dict{String,<:Any})
    apply_pmd!(_correct_branch_directions!, data)
end


"checks that all parallel branches have the same orientation"
function _correct_branch_directions!(pm_data::Dict{String,<:Any})
    orientations = Set()
    for (i, branch) in pm_data["branch"]
        orientation = (branch["f_bus"], branch["t_bus"])
        orientation_rev = (branch["t_bus"], branch["f_bus"])

        if in(orientation_rev, orientations)
            @warn "reversing the orientation of branch $(i) $(orientation) to be consistent with other parallel branches"
            branch_orginal = copy(branch)
            branch["f_bus"] = branch_orginal["t_bus"]
            branch["t_bus"] = branch_orginal["f_bus"]
            branch["f_connections"] = branch_orginal["t_connections"]
            branch["t_connections"] = branch_orginal["f_connections"]
            branch["g_to"] = branch_orginal["g_fr"]
            branch["b_to"] = branch_orginal["b_fr"]
            branch["g_fr"] = branch_orginal["g_to"]
            branch["b_fr"] = branch_orginal["b_to"]
            branch["br_r"] = branch_orginal["br_r"]
            branch["br_x"] = branch_orginal["br_x"]
            branch["angmin"] = -branch_orginal["angmax"]
            branch["angmax"] = -branch_orginal["angmin"]

        else
            push!(orientations, orientation)
        end

    end

end


"""
    check_branch_loops(data::Dict{String,<:Any})

checks that all branches connect two distinct buses
"""
function check_branch_loops(data::Dict{String,<:Any})
    apply_pmd!(_check_branch_loops, data)
end


"checks that all branches connect two distinct buses"
function _check_branch_loops(pm_data::Dict{String, <:Any})
    for (i, branch) in pm_data["branch"]
        if branch["f_bus"] == branch["t_bus"]
            error("both sides of branch $(i) connect to bus $(branch["f_bus"])")
        end
    end
end


"""
    check_connectivity(data::Dict{String,<:Any})

checks that all buses are unique and other components link to valid buses
"""
function check_connectivity(data::Dict{String,<:Any})
    apply_pmd!(_check_connectivity, data)
end


"checks that all buses are unique and other components link to valid buses"
function _check_connectivity(data::Dict{String,<:Any})
    bus_ids = Set(bus["index"] for (i,bus) in data["bus"])
    @assert(length(bus_ids) == length(data["bus"])) # if this is not true something very bad is going on

    for (i, load) in data["load"]
        if !(load["load_bus"] in bus_ids)
            error("bus $(load["load_bus"]) in load $(i) is not defined")
        end
    end

    for (i, shunt) in data["shunt"]
        if !(shunt["shunt_bus"] in bus_ids)
            error("bus $(shunt["shunt_bus"]) in shunt $(i) is not defined")
        end
    end

    for (i, gen) in data["gen"]
        if !(gen["gen_bus"] in bus_ids)
            error("bus $(gen["gen_bus"]) in generator $(i) is not defined")
        end
    end

    for (i, strg) in data["storage"]
        if !(strg["storage_bus"] in bus_ids)
            error("bus $(strg["storage_bus"]) in storage unit $(i) is not defined")
        end
    end

    for (i, switch) in data["switch"]
        if !(switch["f_bus"] in bus_ids)
            error("from bus $(switch["f_bus"]) in switch $(i) is not defined")
        end

        if !(switch["t_bus"] in bus_ids)
            error("to bus $(switch["t_bus"]) in switch $(i) is not defined")
        end
    end

    for (i, branch) in data["branch"]
        if !(branch["f_bus"] in bus_ids)
            error("from bus $(branch["f_bus"]) in branch $(i) is not defined")
        end

        if !(branch["t_bus"] in bus_ids)
            error("to bus $(branch["t_bus"]) in branch $(i) is not defined")
        end
    end
end


"""
checks bus types are suitable for a power flow study, if not, fixes them.
the primary checks are that all type 2 buses (i.e., PV) have a connected and
active generator and there is a single type 3 bus (i.e., slack bus) with an
active connected generator.
assumes that the network is a single connected component
"""
function correct_bus_types!(data::Dict{String,<:Any})
    apply_pmd!(_correct_bus_types!, data)
end


"""
checks bus types are suitable for a power flow study, if not, fixes them.
the primary checks are that all type 2 buses (i.e., PV) have a connected and
active generator and there is a single type 3 bus (i.e., slack bus) with an
active connected generator. Assumes that the network is a single connected component
"""
function _correct_bus_types!(pm_data::Dict{String,<:Any})
    bus_gens = Dict{String,Vector{String}}(i => String[] for (i,bus) in pm_data["bus"])

    for (i,gen) in pm_data["gen"]
        if gen[pmd_math_component_status["gen"]] != pmd_math_component_status_inactive["gen"]
            push!(bus_gens["$(gen["gen_bus"])"], i)
        end
    end

    slack_found = false
    for (i, bus) in pm_data["bus"]
        if bus["bus_type"] == 1
            if !isempty(bus_gens[i]) # PQ
                @warn "active generators found at bus $(bus["bus_i"]), updating to bus type from $(bus["bus_type"]) to 2"
                bus["bus_type"] = 2
            end
        elseif bus["bus_type"] == 2 # PV
            if isempty(bus_gens[i])
                @warn "no active generators found at bus $(bus["bus_i"]), updating to bus type from $(bus["bus_type"]) to 1"
                bus["bus_type"] = 1
            end
        elseif bus["bus_type"] == 3 # Slack
            if !isempty(bus_gens[i])
                slack_found = true
            else
                @warn "no active generators found at bus $(bus["bus_i"]), updating to bus type from $(bus["bus_type"]) to 1"
                bus["bus_type"] = 1
            end
        elseif bus["bus_type"] == 4 # inactive bus
            # do nothing
        else  # unknown bus type
            new_bus_type = 1
            if length(bus_gens[i]) != 0
                new_bus_type = 2
            end
            @warn "bus $(bus["bus_i"]) has an unrecongized bus_type $(bus["bus_type"]), updating to bus_type $(new_bus_type)"
            bus["bus_type"] = new_bus_type
        end
    end

    if !slack_found
        gen = _biggest_generator(pm_data["gen"])
        if length(gen) > 0
            gen_bus = gen["gen_bus"]
            ref_bus = pm_data["bus"]["$(gen_bus)"]
            ref_bus["bus_type"] = 3
            @warn "no reference bus found, setting bus $(gen_bus) as reference based on generator $(gen["index"])"
        else
            error("no generators found in the given network data, correct_bus_types! requires at least one generator at the reference bus")
        end
    end

end


"find the largest active generator in a collection of generators"
function _biggest_generator(gens::Dict)::Dict
    if length(gens) == 0
        error("generator list passed to _biggest_generator was empty.  please report this bug.")
    end

    biggest_gen = Dict{String,Any}()
    biggest_value = -Inf

    for (k,gen) in gens
        if gen["gen_status"] != 0
            pmax = maximum(get(gen, "pmax", fill(Inf, length(gen["connections"]))))
            if pmax > biggest_value
                biggest_gen = gen
                biggest_value = pmax
            end
        end
    end

    return biggest_gen
end


"""
    correct_cost_functions!(data::Dict{String,<:Any})

throws warnings if cost functions are malformed
"""
function correct_cost_functions!(data::Dict{String,<:Any})
    apply_pmd!(_correct_cost_functions!, data)
end


"throws warnings if cost functions are malformed"
function _correct_cost_functions!(pm_data::Dict{String,<:Any})
    for (i,gen) in pm_data["gen"]
        _correct_cost_function!(i, gen, "generator", "pmin", "pmax")
    end
end


"throws warnings if cost functions are malformed"
function _correct_cost_function!(id, comp, type_name, pmin_key, pmax_key)
    if "model" in keys(comp) && "cost" in keys(comp)
        if comp["model"] == 1
            if length(comp["cost"]) != 2*comp["ncost"]
                error("ncost of $(comp["ncost"]) not consistent with $(length(comp["cost"])) cost values on $(type_name) $(id)")
            end
            if length(comp["cost"]) < 4
                error("cost includes $(comp["ncost"]) points, but at least two points are required on $(type_name) $(id)")
            end

            _remove_pwl_cost_duplicates!(id, comp, type_name)

            for i in 3:2:length(comp["cost"])
                if comp["cost"][i-2] >= comp["cost"][i]
                    error("non-increasing x values in pwl cost model on $(type_name) $(id)")
                end
            end

            _simplify_pwl_cost!(id, comp, type_name)
        elseif comp["model"] == 2
            if length(comp["cost"]) != comp["ncost"]
                error("ncost of $(comp["ncost"]) not consistent with $(length(comp["cost"])) cost values on $(type_name) $(id)")
            end
        else
            @warn "Unknown cost model of type $(comp["model"]) on $(type_name) $(id)"
        end
    end

end


"checks that each point in the a pwl function is unqiue, simplifies the function if duplicates appear"
function _remove_pwl_cost_duplicates!(id, comp, type_name; tolerance=1e-2)
    @assert comp["model"] == 1

    unique_costs = Float64[comp["cost"][1], comp["cost"][2]]
    for i in 3:2:length(comp["cost"])
        x1 = unique_costs[end-1]
        y1 = unique_costs[end]
        x2 = comp["cost"][i+0]
        y2 = comp["cost"][i+1]
        if !(isapprox(x1, x2) && isapprox(y1, y2))
            push!(unique_costs, x2)
            push!(unique_costs, y2)
        end
    end

    if length(unique_costs) < length(comp["cost"])
        @warn "removing duplicate points from pwl cost on $(type_name) $(id), $(comp["cost"]) -> $(unique_costs)"
        comp["cost"] = unique_costs
        comp["ncost"] = div(length(unique_costs), 2)
        return true
    end
    return false
end


"checks the slope of each segment in a pwl function, simplifies the function if the slope changes is below a tolerance"
function _simplify_pwl_cost!(id, comp, type_name; tolerance=1e-2)
    @assert comp["model"] == 1

    slopes = Float64[]
    smpl_cost = Float64[]
    prev_slope = nothing

    x2, y2 = 0.0, 0.0

    for i in 3:2:length(comp["cost"])
        x1 = comp["cost"][i-2]
        y1 = comp["cost"][i-1]
        x2 = comp["cost"][i-0]
        y2 = comp["cost"][i+1]

        m = (y2 - y1)/(x2 - x1)

        if prev_slope == nothing || (abs(prev_slope - m) > tolerance)
            push!(smpl_cost, x1)
            push!(smpl_cost, y1)
            prev_slope = m
        end

        push!(slopes, m)
    end

    push!(smpl_cost, x2)
    push!(smpl_cost, y2)

    if length(smpl_cost) < length(comp["cost"])
        @warn "simplifying pwl cost on $(type_name) $(id), $(comp["cost"]) -> $(smpl_cost)"
        comp["cost"] = smpl_cost
        comp["ncost"] = div(length(smpl_cost), 2)
        return true
    end
    return false
end


"""
    simplify_cost_terms!(data::Dict{String,<:Any})

trims zeros from higher order cost terms
"""
function simplify_cost_terms!(data::Dict{String,<:Any})
    apply_pmd!(_simplify_cost_terms!, data)
end


""
function _simplify_cost_terms!(pm_data::Dict{String,<:Any})

    if haskey(pm_data, "gen")
        for (i, gen) in pm_data["gen"]
            if haskey(gen, "model") && gen["model"] == 2
                ncost = length(gen["cost"])
                for j in 1:ncost
                    if gen["cost"][1] == 0.0
                        gen["cost"] = gen["cost"][2:end]
                    else
                        break
                    end
                end
                if length(gen["cost"]) != ncost
                    gen["ncost"] = length(gen["cost"])
                    @info "removing $(ncost - gen["ncost"]) cost terms from generator $(i): $(gen["cost"])"
                end
            end
        end
    end
end


"""
    standardize_cost_terms!(data::Dict{String,<:Any}; order=-1)

ensures all polynomial costs functions have the same number of terms
"""
function standardize_cost_terms!(data::Dict{String,<:Any}; order=-1)
    pm_data = get_pmd_data(data)

    comp_max_order = 1

    if ismultinetwork(pm_data)
        networks = pm_data["nw"]
    else
        networks = [("0", pm_data)]
    end

    for (i, network) in networks
        if haskey(network, "gen")
            for (i, gen) in network["gen"]
                if haskey(gen, "model") && gen["model"] == 2
                    max_nonzero_index = 1
                    for i in 1:length(gen["cost"])
                        max_nonzero_index = i
                        if gen["cost"][i] != 0.0
                            break
                        end
                    end

                    max_oder = length(gen["cost"]) - max_nonzero_index + 1

                    comp_max_order = max(comp_max_order, max_oder)
                end
            end
        end
    end

    if comp_max_order <= order+1
        comp_max_order = order+1
    else
        if order != -1 # if not the default
            @warn "a standard cost order of $(order) was requested but the given data requires an order of at least $(comp_max_order-1)"
        end
    end

    for (i, network) in networks
        if haskey(network, "gen")
            _standardize_cost_terms!(network["gen"], comp_max_order, "generator"; nw=ismultinetwork(data) ? i : "")
        end
    end

end


"ensures all polynomial costs functions have at exactly comp_order terms"
function _standardize_cost_terms!(components::Dict{String,<:Any}, comp_order::Int, cost_comp_name::String; nw::String="")
    modified = Set{Int}()
    for (i, comp) in components
        if haskey(comp, "model") && comp["model"] == 2 && length(comp["cost"]) != comp_order
            std_cost = [0.0 for i in 1:comp_order]
            current_cost = reverse(comp["cost"])
            #println("gen cost: $(comp["cost"])")
            for i in 1:min(comp_order, length(current_cost))
                std_cost[i] = current_cost[i]
            end
            comp["cost"] = reverse(std_cost)
            comp["ncost"] = comp_order
            #println("std gen cost: $(comp["cost"])")

            subnet = !isempty(nw) ? " on subnetwork $nw" : ""
            @info "updated $(cost_comp_name) $(comp["index"])$(subnet) cost function with order $(length(current_cost)) to a function of order $(comp_order): $(comp["cost"])"
            push!(modified, comp["index"])
        end
    end
    return modified
end


"infer the internal dimension of a winding, load or generator based on the connections and the configuration"
function _infer_int_dim(connections::Vector, configuration::ConnConfig, kron_reduced)
    if configuration==WYE
        if kron_reduced
            return length(connections)
        else
            return length(connections)-1
        end
    else # DELTA
        if length(connections)==2
            return 1
        elseif length(connections)==3
            return 3
        else
            error("Only 1 and 3 phase delta-connections are supported.")
        end
    end
end


"infer the internal dimension for a unit, i.e. any one-port component with `connections` and `configuration` properties"
function _infer_int_dim_unit(unit::Dict{String,<:Any}, kron_reduced)
    return _infer_int_dim(unit["connections"], unit["configuration"], kron_reduced)
end


"infer the internal dimension for a transformer (only in the MATHEMATICAL data model format)"
function _infer_int_dim_transformer(trans::Dict{String,<:Any}, kron_reduced)
    return _infer_int_dim(trans["f_connections"], trans["configuration"], kron_reduced)
end
