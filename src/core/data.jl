import LinearAlgebra: Adjoint

const _excluded_count_busname_patterns = Vector{Regex}([
    r"^_virtual.*",
])

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
function count_nodes(data::Dict{String,<:Any})::Int
    n_nodes = 0

    if get(data, "data_model", missing) == DSS
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
    elseif get(data, "data_model", missing) in [MATHEMATICAL, ENGINEERING] || (haskey(data, "source_type") && data["source_type"] == "matlab")
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
        Memento.error(_LOGGER, "Origin of data structure not recognized, cannot count nodes reliably")
    end

    return n_nodes
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
        Memento.error(_LOGGER, "Zig-zag not yet supported.")
    end

    return tm_nom
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
        id = bus["index"]
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
function _calc_load_pq_bounds(load::Dict, bus::Dict)
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
function _calc_load_vbounds(load::Dict, bus::Dict)
    if load["configuration"]==WYE
        vmin = bus["vmin"]
        vmax = bus["vmax"]
    elseif load["configuration"]==DELTA
        vmin, vmax = _calc_bus_vm_ll_bounds(bus)
    end
    return vmin, vmax
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
function _calc_gen_current_max(gen::Dict, bus::Dict)
    if all([haskey(gen, prop) for prop in ["pmax", "pmin", "qmax", "qmin"]]) && haskey(bus, "vmin")
        pabsmax = max.(abs.(gen["pmin"]), abs.(gen["pmax"]))
        qabsmax = max.(abs.(gen["qmin"]), abs.(gen["qmax"]))
        smax = sqrt.(pabsmax.^2 + qabsmax.^2)

        vmin = bus["vmin"]

        return smax./vmin
    else
        N = 3 #TODO update for 4-wire
        return fill(Inf, N)
    end
end


"""
Returns a total (shunt+series) current magnitude bound for the from and to side
of a branch. The total power rating also implies a current bound through the
lower bound on the voltage magnitude of the connected buses.
"""
function _calc_branch_current_max(branch::Dict, bus::Dict)
    bounds = []

    if haskey(branch, "c_rating_a")
        push!(bounds, branch["c_rating_a"])
    end
    if haskey(branch, "rate_a") && haskey(bus, "vmin")
        push!(bounds, branch["rate_a"]./bus["vmin"])
    end

    N = 3 #TODO update for 4-wire
    return min.(fill(Inf, N), bounds...)
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
        push!(bounds_fr, branch["c_rating_a"])
        push!(bounds_to, branch["c_rating_a"])
    end
    if haskey(branch, "rate_a")
        push!(bounds_fr, branch["rate_a"]./bus_fr["vmin"])
        push!(bounds_to, branch["rate_a"]./bus_to["vmin"])
    end

    N = 3 #TODO update for 4-wire
    return min.(fill(Inf, N), bounds_fr...), min.(fill(Inf, N), bounds_to...)
end


"""
Returns a power magnitude bound for the from and to side of a transformer.
The total current rating also implies a current bound through the
upper bound on the voltage magnitude of the connected buses.
"""
function _calc_transformer_power_ub_frto(trans::Dict, bus_fr::Dict, bus_to::Dict)
    bounds_fr = []
    bounds_to = []
    #TODO redefine transformer bounds
    # if haskey(trans, "c_rating_a")
    #     push!(bounds_fr, trans["c_rating_a"].*bus_fr["vmax"])
    #     push!(bounds_to, trans["c_rating_a"].*bus_to["vmax"])
    # end
    # if haskey(trans, "rate_a")
    #     push!(bounds_fr, trans["rate_a"])
    #     push!(bounds_to, trans["rate_a"])
    # end


    N = 3 #TODO update for 4-wire
    return min.(fill(Inf, N), bounds_fr...), min.(fill(Inf, N), bounds_to...)
end


"""
Returns a current magnitude bound for the from and to side of a transformer.
The total power rating also implies a current bound through the lower bound on
the voltage magnitude of the connected buses.
"""
function _calc_transformer_current_max_frto(trans::Dict, bus_fr::Dict, bus_to::Dict)
    bounds_fr = []
    bounds_to = []
    #TODO redefine transformer bounds
    # if haskey(trans, "c_rating_a")
    #     push!(bounds_fr, trans["c_rating_a"].*bus_fr["vmax"])
    #     push!(bounds_to, trans["c_rating_a"].*bus_to["vmax"])
    # end
    # if haskey(trans, "rate_a")
    #     push!(bounds_fr, trans["rate_a"])
    #     push!(bounds_to, trans["rate_a"])
    # end


    N = 3 #TODO update for 4-wire
    return min.(fill(Inf, N), bounds_fr...), min.(fill(Inf, N), bounds_to...)
end


"""
Returns a total (shunt+series) power magnitude bound for the from and to side
of a branch. The total current rating also implies a current bound through the
upper bound on the voltage magnitude of the connected buses.
"""
function _calc_branch_power_max(branch::Dict, bus::Dict)
    bounds = []

    if haskey(branch, "c_rating_a") && haskey(bus, "vmax")
        push!(bounds, branch["c_rating_a"].*bus["vmax"])
    end
    if haskey(branch, "rate_a")
        push!(bounds, branch["rate_a"])
    end

    N = 3 #TODO update for 4-wire
    return min.(fill(Inf, N), bounds...)
end


"""
Returns a total (shunt+series) power magnitude bound for the from and to side
of a branch. The total current rating also implies a current bound through the
upper bound on the voltage magnitude of the connected buses.
"""
function _calc_branch_power_max_frto(branch::Dict, bus_fr::Dict, bus_to::Dict)
    return _calc_branch_power_max(branch, bus_fr), _calc_branch_power_max(branch, bus_to)
end


"""
Returns a valid series current magnitude bound for a branch.
"""
function _calc_branch_series_current_max(branch::Dict, bus_fr::Dict, bus_to::Dict)
    ncnds = 3 #TODO update for four-wire
    vmin_fr = get(bus_fr, "vmin", fill(0.0, ncnds))
    vmin_to = get(bus_to, "vmin", fill(0.0, ncnds))

    vmax_fr = get(bus_fr, "vmax", fill(Inf, ncnds))
    vmax_to = get(bus_to, "vmax", fill(Inf, ncnds))

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
    N = 3 #TODO update for 4-wire
    return min.(fill(Inf, N), c_max_fr_sh.+c_max_fr_tot, c_max_to_sh.+c_max_to_tot)
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
                                item_ref_data[param] = LinearAlgebra.diagm(0=>fill(value, conductors))
                            else
                                item_ref_data[param] = fill(value, conductors)
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

    for (_, load) in data["load"]
        load["model"] = POWER
        load["configuration"] = WYE
    end

    for (_, load) in data["gen"]
        load["configuration"] = WYE
    end
end


"https://stackoverflow.com/questions/39039553/lower-triangular-matrix-in-julia"
function _vec2utri!(v::Vector{T}) where T
    d = length(v)
    n = Int((sqrt(8d+1)+1)/2)
    n*(n-1)/2 == d || error("vec2utri: length of vector is not triangular")
    [ i<j ? v[Int((j-1)*(j-2)/2)+i] : 0 for i=1:n, j=1:n ]
end


""
function _vec2ltri!(v::Vector{T}) where T
    _vec2utri!(v)'
end


""
function _mat2utrivec!(m::Union{Matrix{T}, LinearAlgebra.Symmetric{T}}) where T
    @assert size(m,1) == size(m,2)
    n = size(m,1)
    [m[i,j] for i=1:n, j=1:n if i < j]
end


""
function _mat2ltrivec!(m::Union{Matrix{T}, LinearAlgebra.Symmetric{T}}) where T
    @assert size(m,1) == size(m,2)
    n = size(m,1)
    [m[j,i] for i=1:n, j=1:n if i < j]
end


""
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


""
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
    end
    return false
end


macro smart_constraint(model, vars, expr)
    esc(quote
        if _has_nl_expression($vars)
            JuMP.@NLconstraint($model, $expr)
        else
            JuMP.@constraint($model, $expr)
        end
    end)
end


"Local wrapper method for JuMP.set_lower_bound, which skips NaN and infinite (-Inf only)"
function set_lower_bound(x::JuMP.VariableRef, bound; loose_bounds::Bool=false, pm=missing, category=:default)
    if !(isnan(bound) || bound==-Inf)
        JuMP.set_lower_bound(x, bound)
    elseif loose_bounds
        lbs = pm.ext[:loose_bounds]
        JuMP.set_lower_bound(x, -lbs.bound_values[category])
        push!(lbs.loose_lb_vars, x)
    end
end


"Local wrapper method for JuMP.set_upper_bound, which skips NaN and infinite (+Inf only)"
function set_upper_bound(x::JuMP.VariableRef, bound; loose_bounds::Bool=false, pm=missing, category=:default)
    if !(isnan(bound) || bound==Inf)
        JuMP.set_upper_bound(x, bound)
    elseif loose_bounds
        lbs = pm.ext[:loose_bounds]
        JuMP.set_upper_bound(x, lbs.bound_values[category])
        push!(lbs.loose_ub_vars, x)
    end
end


""
function sol_polar_voltage!(pm::_PM.AbstractPowerModel, solution::Dict)
    if haskey(solution, "nw")
        nws_data = solution["nw"]
    else
        nws_data = Dict("0" => solution)
    end

    for (n, nw_data) in nws_data
        if haskey(nw_data, "bus")
            for (i,bus) in nw_data["bus"]
                if haskey(bus, "vr") && haskey(bus, "vi")
                    vr = bus["vr"]
                    vi = bus["vi"]
                    if isa(vr, Dict)
                        bus["vm"] = Dict(t=>abs(vr[t]+im*vi[t]) for t in keys(vr))
                        bus["va"] = Dict(t=>_wrap_to_pi(atan(vi[t], vr[t])) for t in keys(vr))
                    else
                        bus["vm"] = abs.(vr .+ im*vi)
                        bus["va"] = _wrap_to_pi(atan.(vi, vr))
                    end
                end
            end
        end
    end
end

# BOUND manipulation methods (0*Inf->0 is often desired)
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


""
function _nan2zero(b, a; val=0)
    and(x, y) = x && y
    b[and.(isnan.(b), a.==val)] .= 0.0
    return b
end


"Counts active ungrounded connections on edge components"
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
                                if edge_type == "transformer" && component["configuration"] == WYE && terminal != connections[end]
                                    push!(counted_connections, terminal)
                                    active_connections += 1
                                elseif !(terminal in data["bus"][bus]["grounded"])
                                    push!(counted_connections, terminal)
                                    active_connections += 1
                                end
                            else
                                if edge_type == "transformer"
                                    if component["configuration"] == DELTA || (component["configuration"] == WYE && terminal != connections[end])
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


"Counts active ungrounded terminals on buses"
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
