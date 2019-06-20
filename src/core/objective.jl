"a quadratic penalty for bus power slack variables"
function objective_min_bus_power_slack(pm::_PMs.GenericPowerModel)
    return JuMP.@objective(pm.model, Min,
        sum(
            sum(
                sum( _PMs.var(pm, n, c, :p_slack, i)^2 + _PMs.var(pm, n, c, :q_slack, i)^2 for (i,bus) in nw_ref[:bus])
            for c in _PMs.conductor_ids(pm, n))
        for (n, nw_ref) in _PMs.nws(pm))
    )
end
