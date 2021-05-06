"""
    objective_mc_min_slack_bus_power(pm::AbstractUnbalancedPowerModel)

a quadratic penalty for bus power slack variables
"""
function objective_mc_min_slack_bus_power(pm::AbstractUnbalancedPowerModel)
    return JuMP.@objective(pm.model, Min,
        sum(
            sum(
                sum( var(pm, n, :p_slack, i)[t]^2 + var(pm, n, :q_slack, i)[t]^2 for t in ref(pm, n, :bus, i, "terminals")
                ) for (i,bus) in nw_ref[:bus]
            ) for (n, nw_ref) in nws(pm))
        )
end


"""
    objective_mc_min_load_setpoint_delta(pm::AbstractUnbalancedPowerModel)

minimum load delta objective with storage
"""
function objective_mc_min_load_setpoint_delta(pm::AbstractUnbalancedPowerModel)
    for (n, nw_ref) in nws(pm)
        var(pm, n)[:delta_pg] = Dict(i => JuMP.@variable(pm.model,
                [c in ref(pm, n, :gen, i)["connections"]], base_name="$(n)_$(i)_delta_pg",
                start = 0.0) for i in ids(pm, n, :gen))

        var(pm, n)[:delta_ps] = Dict(i => JuMP.@variable(pm.model,
                [c in ref(pm, n, :storage, i)["connections"]], base_name="$(n)_$(i)_delta_ps",
                start = 0.0) for i in ids(pm, n, :storage))

        for (i, gen) in nw_ref[:gen]
            for (idx, c) in enumerate(gen["connections"])
                JuMP.@constraint(pm.model, var(pm, n, :delta_pg, i)[c] >=  (gen["pg"][idx] - var(pm, n, :pg, i)[c]))
                JuMP.@constraint(pm.model, var(pm, n, :delta_pg, i)[c] >= -(gen["pg"][idx] - var(pm, n, :pg, i)[c]))
            end
        end

        for (i, strg) in nw_ref[:storage]
            for (idx, c) in enumerate(strg["connections"])
                JuMP.@constraint(pm.model, var(pm, n, :delta_ps, i)[c] >=  (strg["ps"][idx] - var(pm, n, :ps, i)[c]))
                JuMP.@constraint(pm.model, var(pm, n, :delta_ps, i)[c] >= -(strg["ps"][idx] - var(pm, n, :ps, i)[c]))
            end
        end
    end

    w = Dict(n => Dict(i => 10*get(load, "weight", 1.0) for (i,load) in ref(pm, n, :load)) for n in nw_ids(pm))

    JuMP.@objective(pm.model, Min,
        sum(
            sum(                      10*(1 - var(pm, n, :z_voltage, i)) for (i,bus) in nw_ref[:bus]) +
            sum( w[n][i]*sum(load["pd"])*(1 - var(pm, n, :z_demand, i)) for (i,load) in nw_ref[:load]) +
            sum(        sum(shunt["gs"])*(1 - var(pm, n, :z_shunt, i)) for (i,shunt) in nw_ref[:shunt]) +
            sum( sum(                         var(pm, n, :delta_pg, i)[c] for (idx,c) in enumerate(gen["connections"])) for (i,gen)  in nw_ref[:gen]) +
            sum( sum(                         var(pm, n, :delta_ps, i)[c] for (idx,c) in enumerate(strg["connections"])) for (i,strg) in nw_ref[:storage])
        for (n, nw_ref) in nws(pm))
    )
end


"""
    objective_mc_min_load_setpoint_delta_simple(pm::AbstractUnbalancedPowerModel)

simplified minimum load delta objective (continuous load shed)
"""
function objective_mc_min_load_setpoint_delta_simple(pm::AbstractUnbalancedPowerModel)
    JuMP.@objective(pm.model, Min,
        sum(
            sum( ((1 - var(pm, n, :z_demand, i))) for i in keys(nw_ref[:load])) +
            sum( ((1 - var(pm, n, :z_shunt, i))) for (i,shunt) in nw_ref[:shunt])
        for (n, nw_ref) in nws(pm))
    )
end


"""
    objective_mc_min_load_setpoint_delta_simple_switch(pm::AbstractUnbalancedPowerModel)

simplified minimum load delta objective (continuous load shed) including a switch state term
"""
function objective_mc_min_load_setpoint_delta_simple_switch(pm::AbstractUnbalancedPowerModel)
    JuMP.@objective(pm.model, Min,
        sum(
            sum( ((1 - var(pm, n, :z_demand, i))) for i in keys(nw_ref[:load])) +
            sum( ((1 - var(pm, n, :z_shunt, i))) for (i,shunt) in nw_ref[:shunt]) +
            sum( var(pm, n, :switch_state, l) for l in ids(pm, n, :switch_dispatchable))
        for (n, nw_ref) in nws(pm))
    )
end

"""
    objective_mc_max_load_setpoint(pm::AbstractUnbalancedPowerModel)

maximum loadability objective (continuous load shed) with storage
"""
function objective_mc_max_load_setpoint(pm::AbstractUnbalancedPowerModel)
    w = Dict(n => Dict(i => get(load, "weight", 1.0) for (i,load) in ref(pm, n, :load)) for n in nw_ids(pm))

    JuMP.@objective(pm.model, Max,
        sum(
            sum( (              10*var(pm, n, :z_voltage, i) for (i, bus) in nw_ref[:bus])) +
            sum( (         w[n][i]*var(pm, n, :z_demand, i)) for (i,load) in nw_ref[:load]) +
            sum( (                 var(pm, n, :z_shunt, i)) for (i,shunt) in nw_ref[:shunt]) +
            sum( (                 var(pm, n, :z_gen, i) for (i,gen) in nw_ref[:gen])) +
            sum( (                 var(pm, n, :z_storage, i) for (i,storage) in nw_ref[:storage]))
        for (n, nw_ref) in nws(pm))
    )
end


"""
    objective_mc_min_fuel_cost(pm::AbstractUnbalancedPowerModel)

Standard fuel cost minimization objective
"""
function objective_mc_min_fuel_cost(pm::AbstractUnbalancedPowerModel; kwargs...)
    model = check_gen_cost_models(pm)

    if model == 1
        return objective_mc_min_fuel_cost_pwl(pm; kwargs...)
    elseif model == 2
        return objective_mc_min_fuel_cost_polynomial(pm; kwargs...)
    else
        error("Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end


"""
    objective_mc_min_fuel_cost_switch(pm::AbstractUnbalancedPowerModel)

Standard fuel cost minimization objective including switches
"""
function objective_mc_min_fuel_cost_switch(pm::AbstractUnbalancedPowerModel; kwargs...)
    model = check_gen_cost_models(pm)

    if model == 1
        return objective_mc_min_fuel_cost_pwl_switch(pm; kwargs...)
    elseif model == 2
        return objective_mc_min_fuel_cost_polynomial_switch(pm; kwargs...)
    else
        error("Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end



"""
    objective_mc_min_fuel_cost_pwl(pm::AbstractUnbalancedPowerModel)

Fuel cost minimization objective with piecewise linear terms
"""
function objective_mc_min_fuel_cost_pwl(pm::AbstractUnbalancedPowerModel; kwargs...)
    objective_mc_variable_pg_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n, :pg_cost, i) for (i,gen) in nw_ref[:gen])
        for (n, nw_ref) in nws(pm))
    )
end


"""
    objective_mc_min_fuel_cost_pwl_switch(pm::AbstractUnbalancedPowerModel)

Fuel cost minimization objective with piecewise linear terms including switches
"""
function objective_mc_min_fuel_cost_pwl_switch(pm::AbstractUnbalancedPowerModel; kwargs...)
    objective_mc_variable_pg_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n, :pg_cost, i) for (i,gen) in nw_ref[:gen]) +
            sum( var(pm, n, :switch_state, l) for l in ids(pm, n, :switch_dispatchable))
        for (n, nw_ref) in nws(pm))
    )
end


"""
    objective_mc_variable_pg_cost(pm::AbstractUnbalancedPowerModel)

adds pg_cost variables and constraints
"""
function objective_mc_variable_pg_cost(pm::AbstractUnbalancedPowerModel; report::Bool=true)
    for (n, nw_ref) in nws(pm)
        pg_cost = var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in ref(pm, n, :gen)
            points = calc_pwl_points(gen["ncost"], gen["cost"], gen["pmin"], gen["pmax"])

            pg_cost_lambda = JuMP.@variable(pm.model,
                [i in 1:length(points)], base_name="$(n)_pg_cost_lambda",
                lower_bound = 0.0,
                upper_bound = 1.0
            )
            JuMP.@constraint(pm.model, sum(pg_cost_lambda) == 1.0)

            pg_expr = 0.0
            pg_cost_expr = 0.0
            for (i,point) in enumerate(points)
                pg_expr += point.mw*pg_cost_lambda[i]
                pg_cost_expr += point.cost*pg_cost_lambda[i]
            end
            JuMP.@constraint(pm.model, pg_expr == sum(var(pm, n, :pg, i)[c] for c in gen["connections"]))
            pg_cost[i] = pg_cost_expr
        end

        report && _IM.sol_component_value(pm, pmd_it_sym, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)
    end
end


"""
    objective_mc_min_fuel_cost_polynomial(pm::AbstractUnbalancedPowerModel)

Fuel cost minimization objective for polynomial terms
"""
function objective_mc_min_fuel_cost_polynomial(pm::AbstractUnbalancedPowerModel; kwargs...)
    order = calc_max_cost_index(pm.data)-1

    if order <= 2
        return _objective_mc_min_fuel_cost_polynomial_linquad(pm; kwargs...)
    else
        return _objective_mc_min_fuel_cost_polynomial_nl(pm; kwargs...)
    end
end


"""
    objective_mc_min_fuel_cost_polynomial_switch(pm::AbstractUnbalancedPowerModel)

Fuel cost minimization objective for polynomial terms including switches
"""
function objective_mc_min_fuel_cost_polynomial_switch(pm::AbstractUnbalancedPowerModel; kwargs...)
    order = calc_max_cost_index(pm.data)-1

    if order <= 2
        return _objective_mc_min_fuel_cost_polynomial_linquad_switch(pm; kwargs...)
    else
        return _objective_mc_min_fuel_cost_polynomial_nl_switch(pm; kwargs...)
    end
end


"gen connections adaptation of min fuel cost polynomial linquad objective"
function _objective_mc_min_fuel_cost_polynomial_linquad(pm::AbstractUnbalancedPowerModel; report::Bool=true)
    gen_cost = Dict()
    for (n, nw_ref) in nws(pm)
        for (i,gen) in nw_ref[:gen]
            pg = sum( var(pm, n, :pg, i)[c] for c in gen["connections"] )

            if length(gen["cost"]) == 1
                gen_cost[(n,i)] = gen["cost"][1]
            elseif length(gen["cost"]) == 2
                gen_cost[(n,i)] = gen["cost"][1]*pg + gen["cost"][2]
            elseif length(gen["cost"]) == 3
                gen_cost[(n,i)] = gen["cost"][1]*pg^2 + gen["cost"][2]*pg + gen["cost"][3]
            else
                gen_cost[(n,i)] = 0.0
            end
        end
    end

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] )
        for (n, nw_ref) in nws(pm))
    )
end


"gen connections adaptation of min fuel cost polynomial linquad objective"
function _objective_mc_min_fuel_cost_polynomial_linquad_switch(pm::AbstractUnbalancedPowerModel; report::Bool=true)
    gen_cost = Dict()
    for (n, nw_ref) in nws(pm)
        for (i,gen) in nw_ref[:gen]
            pg = sum( var(pm, n, :pg, i)[c] for c in gen["connections"] )

            if length(gen["cost"]) == 1
                gen_cost[(n,i)] = gen["cost"][1]
            elseif length(gen["cost"]) == 2
                gen_cost[(n,i)] = gen["cost"][1]*pg + gen["cost"][2]
            elseif length(gen["cost"]) == 3
                gen_cost[(n,i)] = gen["cost"][1]*pg^2 + gen["cost"][2]*pg + gen["cost"][3]
            else
                gen_cost[(n,i)] = 0.0
            end
        end
    end

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] ) +
            sum( var(pm, n, :switch_state, l) for l in ids(pm, n, :switch_dispatchable))
        for (n, nw_ref) in nws(pm))
    )
end


"Multiconductor adaptation of min fuel cost polynomial linquad objective"
function _objective_mc_min_fuel_cost_polynomial_linquad(pm::AbstractUnbalancedIVRModel; report::Bool=true)
    gen_cost = Dict()

    for (n, nw_ref) in nws(pm)
        for (i,gen) in nw_ref[:gen]
            bus = gen["gen_bus"]

            #to avoid function calls inside of @NLconstraint:
            pg = var(pm, n, :pg, i)
            if length(gen["cost"]) == 1
                gen_cost[(n,i)] = gen["cost"][1]
            elseif length(gen["cost"]) == 2
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*sum(pg[c] for c in gen["connections"]) + gen["cost"][2])
            elseif length(gen["cost"]) == 3
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*sum(pg[c] for c in gen["connections"])^2 + gen["cost"][2]*sum(pg[c] for c in gen["connections"]) + gen["cost"][3])
            else
                gen_cost[(n,i)] = 0.0
            end
        end
    end

    return JuMP.@NLobjective(pm.model, Min,
        sum(
            sum(    gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] )
        for (n, nw_ref) in nws(pm))
    )
end


"Multiconductor adaptation of min fuel cost polynomial linquad objective"
function _objective_mc_min_fuel_cost_polynomial_linquad_switch(pm::AbstractUnbalancedIVRModel; report::Bool=true)
    gen_cost = Dict()
    switch_state = Dict()

    for (n, nw_ref) in nws(pm)
        for (i,gen) in nw_ref[:gen]
            bus = gen["gen_bus"]

            #to avoid function calls inside of @NLconstraint:
            pg = var(pm, n, :pg, i)
            if length(gen["cost"]) == 1
                gen_cost[(n,i)] = gen["cost"][1]
            elseif length(gen["cost"]) == 2
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*sum(pg[c] for c in gen["connections"]) + gen["cost"][2])
            elseif length(gen["cost"]) == 3
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*sum(pg[c] for c in gen["connections"])^2 + gen["cost"][2]*sum(pg[c] for c in gen["connections"]) + gen["cost"][3])
            else
                gen_cost[(n,i)] = 0.0
            end
        end
        for i in ids(pm, n, :switch_dispatchable)
            switch_state[(n,i)] = var(pm, n, :switch_state, i)
        end
    end

    return JuMP.@NLobjective(pm.model, Min,
        sum(
            sum(    gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] )
            + sum( switch_state[(n,i)] for i in ids(pm, n, :switch_dispatchable))
        for (n, nw_ref) in nws(pm))
    )
end


""
function _objective_mc_min_fuel_cost_polynomial_nl(pm::AbstractUnbalancedPowerModel; report::Bool=true)
    gen_cost = Dict()
    for (n, nw_ref) in nws(pm)
        for (i,gen) in nw_ref[:gen]
            pg = sum( var(pm, n, :pg, i)[c] for c in gen["connections"] )

            cost_rev = reverse(gen["cost"])
            if length(cost_rev) == 1
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1])
            elseif length(cost_rev) == 2
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*pg)
            elseif length(cost_rev) == 3
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*pg + cost_rev[3]*pg^2)
            elseif length(cost_rev) >= 4
                cost_rev_nl = cost_rev[4:end]
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*pg + cost_rev[3]*pg^2 + sum( v*pg^(d+3) for (d,v) in enumerate(cost_rev_nl)) )
            else
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, 0.0)
            end
        end
    end

    return JuMP.@NLobjective(pm.model, Min,
        sum(
            sum( gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] )
        for (n, nw_ref) in nws(pm))
    )
end


""
function _objective_mc_min_fuel_cost_polynomial_nl_switch(pm::AbstractUnbalancedPowerModel)
    gen_cost = Dict()
    for (n, nw_ref) in nws(pm)
        for (i,gen) in nw_ref[:gen]
            pg = sum( var(pm, n, :pg, i)[c] for c in gen["connections"] )

            cost_rev = reverse(gen["cost"])
            if length(cost_rev) == 1
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1])
            elseif length(cost_rev) == 2
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*pg)
            elseif length(cost_rev) == 3
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*pg + cost_rev[3]*pg^2)
            elseif length(cost_rev) >= 4
                cost_rev_nl = cost_rev[4:end]
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*pg + cost_rev[3]*pg^2 + sum( v*pg^(d+3) for (d,v) in enumerate(cost_rev_nl)) )
            else
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, 0.0)
            end
        end
    end

    return JuMP.@NLobjective(pm.model, Min,
        sum(
            sum( gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] ) +
            sum( var(pm, n, :switch_state, l) for l in ids(pm, n, :switch_dispatchable))
        for (n, nw_ref) in nws(pm))
    )
end


"""
    objective_variable_pg_cost(pm::AbstractUnbalancedIVRModel)

adds pg_cost variables and constraints for the IVR formulation
"""
function objective_variable_pg_cost(pm::AbstractUnbalancedIVRModel)
    for (n, nw_ref) in nws(pm)
        gen_lines = calc_cost_pwl_lines(nw_ref[:gen])

        #to avoid function calls inside of @NLconstraint
        pg_cost = var(pm, n)[:pg_cost] = JuMP.@variable(pm.model,
            [i in ids(pm, n, :gen)], base_name="$(n)_pg_cost",
        )
        report && _IM.sol_component_value(pm, pmd_it_sym, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)

        # gen pwl cost
        for (i, gen) in nw_ref[:gen]
            pg = var(pm, n, :pg, i)
            for line in gen_lines[i]
                JuMP.@NLconstraint(pm.model, pg_cost[i] >= line.slope*sum(pg[c] for c in gen["connections"]) + line.intercept)
            end
        end
    end
end


"""
    check_cost_models(pm::AbstractUnbalancedPowerModel)

Checks that all cost models are of the same type
"""
function check_cost_models(pm::AbstractUnbalancedPowerModel)
    return check_gen_cost_models(pm)
end


"""
    check_gen_cost_models(pm::AbstractUnbalancedPowerModel)

Checks that all generator cost models are of the same type
"""
function check_gen_cost_models(pm::AbstractUnbalancedPowerModel)
    model = nothing

    for (n, nw_ref) in nws(pm)
        for (i,gen) in nw_ref[:gen]
            if haskey(gen, "cost")
                if model === nothing
                    model = gen["model"]
                else
                    if gen["model"] != model
                        error("cost models are inconsistent, the typical model is $(model) however model $(gen["model"]) is given on generator $(i)")
                    end
                end
            else
                error("no cost given for generator $(i)")
            end
        end
    end

    return model
end


"""
    calc_pwl_points(ncost::Int, cost::Vector{<:Real}, pmin::Real, pmax::Real; tolerance=1e-2)

cleans up raw pwl cost points in preparation for building a mathamatical model.
The key mathematical properties,
- the first and last points are strickly outside of the pmin-to-pmax range
- pmin and pmax occur in the first and last line segments.
"""
function calc_pwl_points(ncost::Int, cost::Vector{<:Real}, pmin::Real, pmax::Real; tolerance=1e-2)
    @assert ncost >= 1 && length(cost) >= 2
    @assert 2*ncost == length(cost)
    @assert pmin <= pmax

    if isinf(pmin) || isinf(pmax)
        error("a bounded operating range is required for modeling pwl costs.  Given active power range in $(pmin) - $(pmax)")
    end

    points = []
    for i in 1:ncost
        push!(points, (mw=cost[2*i-1], cost=cost[2*i]))
    end

    first_active = 0
    for i in 1:(ncost-1)
        #mw_0 = points[i].mw
        mw_1 = points[i+1].mw
        first_active = i
        if pmin <= mw_1
            break
        end
    end

    last_active = 0
    for i in 1:(ncost-1)
        mw_0 = points[end - i].mw
        #mw_1 = points[end - i + 1].mw
        last_active = ncost - i + 1
        if pmax >= mw_0
            break
        end
    end

    points = points[first_active : last_active]


    x1 = points[1].mw
    y1 = points[1].cost
    x2 = points[2].mw
    y2 = points[2].cost

    if x1 > pmin
        x0 = pmin - tolerance

        m = (y2 - y1)/(x2 - x1)

        if !isnan(m)
            y0 = y2 - m*(x2 - x0)
            points[1] = (mw=x0, cost=y0)
        else
            points[1] = (mw=x0, cost=y1)
        end

        modified = true
    end


    x1 = points[end-1].mw
    y1 = points[end-1].cost
    x2 = points[end].mw
    y2 = points[end].cost

    if x2 < pmax
        x3 = pmax + tolerance

        m = (y2 - y1)/(x2 - x1)

        if !isnan(m)
            y3 = m*(x3 - x1) + y1

            points[end] = (mw=x3, cost=y3)
        else
            points[end] = (mw=x3, cost=y2)
        end
    end

    return points
end


"""
    calc_max_cost_index(data::Dict{String,<:Any})

Computes maximum cost index
"""
function calc_max_cost_index(data::Dict{String,<:Any})
    if ismultinetwork(data)
        max_index = 0
        for (i,nw_data) in data["nw"]
            nw_max_index = _calc_max_cost_index(nw_data)
            max_index = max(max_index, nw_max_index)
        end
        return max_index
    else
        return _calc_max_cost_index(data)
    end
end


"""
    _calc_max_cost_index(data::Dict{String,<:Any})

Computes maximum cost index of subnetworks
"""
function _calc_max_cost_index(data::Dict{String,<:Any})
    max_index = 0

    for (i,gen) in data["gen"]
        if haskey(gen, "model")
            if gen["model"] == 2
                if haskey(gen, "cost")
                    max_index = max(max_index, length(gen["cost"]))
                end
            else
                @warn "skipping cost generator $(i) cost model in calc_cost_order, only model 2 is supported."
            end
        end
    end

    return max_index
end


"""
    calc_cost_pwl_lines(comp_dict::Dict)

compute lines in m and b from from pwl cost models data is a list of components.
Can be run on data or ref data structures
"""
function calc_cost_pwl_lines(comp_dict::Dict)
    lines = Dict()
    for (i,comp) in comp_dict
        lines[i] = _calc_comp_lines(comp)
    end
    return lines
end


"""
    _calc_comp_lines(component::Dict{String,<:Any})

compute lines in m and b from from pwl cost models
"""
function _calc_comp_lines(component::Dict{String,<:Any})
    @assert component["model"] == 1
    points = component["cost"]

    line_data = []
    for i in 3:2:length(points)
        x1 = points[i-2]
        y1 = points[i-1]
        x2 = points[i-0]
        y2 = points[i+1]

        m = (y2 - y1)/(x2 - x1)
        b = y1 - m * x1

        push!(line_data, (slope=m, intercept=b))
    end

    for i in 2:length(line_data)
        if line_data[i-1].slope > line_data[i].slope
            error("non-convex pwl function found in points $(component["cost"])\nlines: $(line_data)")
        end
    end

    return line_data
end
