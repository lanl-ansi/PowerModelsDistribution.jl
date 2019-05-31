"a quadratic penalty for bus power slack variables"
function objective_min_bus_power_slack(pm::PMs.GenericPowerModel)
    return JuMP.@objective(pm.model, Min,
        sum(
            sum(
                sum( PMs.var(pm, n, c, :p_slack, i)^2 + PMs.var(pm, n, c, :q_slack, i)^2 for (i,bus) in nw_ref[:bus])
            for c in PMs.conductor_ids(pm, n))
        for (n, nw_ref) in PMs.nws(pm))
    )
end
