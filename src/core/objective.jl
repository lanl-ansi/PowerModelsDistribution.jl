"a quadratic penalty for bus power slack variables"
function objective_min_bus_power_slack(pm::_PM.AbstractPowerModel)
    return JuMP.@objective(pm.model, Min,
        sum(
            sum(
                sum( var(pm, n, :p_slack, i)[c]^2 + var(pm, n, :q_slack, i)[c]^2 for (i,bus) in nw_ref[:bus])
            for c in conductor_ids(pm, n))
        for (n, nw_ref) in nws(pm))
    )
end


"minimum load delta objective (continuous load shed) with storage"
function objective_mc_min_load_delta(pm::_PM.AbstractPowerModel)
    for (n, nw_ref) in nws(pm)
        var(pm, n)[:delta_pg] = Dict(i => JuMP.@variable(pm.model,
                [c in conductor_ids(pm, n)], base_name="$(n)_$(i)_delta_pg",
                start = 0.0) for i in ids(pm, n, :gen))

        var(pm, n)[:delta_ps] = Dict(i => JuMP.@variable(pm.model,
                [c in conductor_ids(pm, n)], base_name="$(n)_$(i)_delta_ps",
                start = 0.0) for i in ids(pm, n, :storage))

        for c in conductor_ids(pm, n)
            for (i, gen) in nw_ref[:gen]
                JuMP.@constraint(pm.model, var(pm, n, :delta_pg, i)[c] >=  (gen["pg"][c] - var(pm, n, :pg, i)[c]))
                JuMP.@constraint(pm.model, var(pm, n, :delta_pg, i)[c] >= -(gen["pg"][c] - var(pm, n, :pg, i)[c]))
            end

            for (i, strg) in nw_ref[:storage]
                JuMP.@constraint(pm.model, var(pm, n, :delta_ps, i)[c] >=  (strg["ps"][c] - var(pm, n, :ps, i)[c]))
                JuMP.@constraint(pm.model, var(pm, n, :delta_ps, i)[c] >= -(strg["ps"][c] - var(pm, n, :ps, i)[c]))
            end

        end
    end

    load_weight = Dict(n => Dict(i => get(load, "weight", 1.0) for (i,load) in ref(pm, n, :load)) for n in nw_ids(pm))
    M = Dict(n => Dict(c => 10*maximum([load_weight[n][i]*abs(load["pd"][c]) for (i,load) in ref(pm, n, :load)]) for c in conductor_ids(pm, n)) for n in nw_ids(pm))

    JuMP.@objective(pm.model, Min,
        sum(
            sum(
                sum( (                 10*(1 - var(pm, n, :z_voltage, i)) for i in keys(nw_ref[:bus]))) +
                sum( (            M[n][c]*(1 - var(pm, n, :z_demand, i))) for i in keys(nw_ref[:load])) +
                sum( (abs(shunt["gs"][c])*(1 - var(pm, n, :z_shunt, i))) for (i,shunt) in nw_ref[:shunt]) +
                sum( (                     var(pm, n, :delta_pg, i)[c] for i in keys(nw_ref[:gen]))) +
                sum( (                     var(pm, n, :delta_ps, i)[c] for i in keys(nw_ref[:storage])))
            for c in conductor_ids(pm, n))
        for (n, nw_ref) in nws(pm))
    )
end



"maximum loadability objective (continuous load shed) with storage"
function objective_mc_max_loadability(pm::_PM.AbstractPowerModel)
    load_weight = Dict(n => Dict(i => get(load, "weight", 1.0) for (i,load) in ref(pm, n, :load)) for n in nw_ids(pm))
    M = Dict(n => Dict(c => 10*maximum([load_weight[n][i]*abs(load["pd"][c]) for (i,load) in ref(pm, n, :load)]) for c in conductor_ids(pm, n)) for n in nw_ids(pm))

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
