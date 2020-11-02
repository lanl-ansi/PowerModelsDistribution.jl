"a quadratic penalty for bus power slack variables"
function objective_mc_min_slack_bus_power(pm::_PM.AbstractPowerModel)
    return JuMP.@objective(pm.model, Min,
        sum(
            sum(
                sum( var(pm, n, :p_slack, i)[c]^2 + var(pm, n, :q_slack, i)[c]^2 for (i,bus) in nw_ref[:bus])
            for c in conductor_ids(pm, n))
        for (n, nw_ref) in nws(pm))
    )
end


"minimum load delta objective (continuous load shed) with storage"
function objective_mc_min_load_setpoint_delta(pm::_PM.AbstractPowerModel)
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


"simplified minimum load delta objective (continuous load shed)"
function objective_mc_min_load_setpoint_delta_simple(pm::_PM.AbstractPowerModel)
    JuMP.@objective(pm.model, Min,
        sum(
            sum( ((1 - var(pm, n, :z_demand, i))) for i in keys(nw_ref[:load])) +
            sum( ((1 - var(pm, n, :z_shunt, i))) for (i,shunt) in nw_ref[:shunt])
        for (n, nw_ref) in nws(pm))
    )
end


"maximum loadability objective (continuous load shed) with storage"
function objective_mc_max_load_setpoint(pm::_PM.AbstractPowerModel)
    load_weight = Dict(n => Dict(i => get(load, "weight", 1.0) for (i,load) in ref(pm, n, :load)) for n in nw_ids(pm))

    JuMP.@objective(pm.model, Max,
        sum(
            sum(
                sum( (              10*var(pm, n, :z_voltage, i) for (i, bus) in nw_ref[:bus])) +
                sum( (         M[n][c]*var(pm, n, :z_demand, i)) for (i,load) in nw_ref[:load]) +
                sum( (abs(shunt["gs"])*var(pm, n, :z_shunt, i)) for (i,shunt) in nw_ref[:shunt]) +
                sum( (                 var(pm, n, :z_gen, i) for (i,gen) in nw_ref[:gen])) +
                sum( (                 var(pm, n, :z_storage, i) for (i,storage) in nw_ref[:storage]))
            for c in conductor_ids(pm, n))
        for (n, nw_ref) in nws(pm))
    )
end


"gen connections adaptation of min fuel cost polynomial linquad objective"
function _PM._objective_min_fuel_cost_polynomial_linquad(pm::_PM.AbstractPowerModel; report::Bool=true)
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


"Multiconductor adaptation of min fuel cost polynomial linquad objective"
function _PM._objective_min_fuel_cost_polynomial_linquad(pm::_PM.AbstractIVRModel; report::Bool=true)
    gen_cost = Dict()
    dcline_cost = Dict()

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
            + sum( dcline_cost[(n,i)] for (i,dcline) in nw_ref[:dcline] )
        for (n, nw_ref) in nws(pm))
    )
end


"adds pg_cost variables and constraints"
function objective_variable_pg_cost(pm::_PM.AbstractIVRModel; report::Bool=true)
    for (n, nw_ref) in nws(pm)
        gen_lines = calc_cost_pwl_lines(nw_ref[:gen])

        #to avoid function calls inside of @NLconstraint
        pg_cost = var(pm, n)[:pg_cost] = JuMP.@variable(pm.model,
            [i in ids(pm, n, :gen)], base_name="$(n)_pg_cost",
        )
        report && _IM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)

        # gen pwl cost
        for (i, gen) in nw_ref[:gen]
            pg = var(pm, n, :pg, i)
            for line in gen_lines[i]
                JuMP.@NLconstraint(pm.model, pg_cost[i] >= line.slope*sum(pg[c] for c in gen["connections"]) + line.intercept)
            end
        end
    end
end
