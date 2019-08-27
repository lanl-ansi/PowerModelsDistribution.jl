"a quadratic penalty for bus power slack variables"
function objective_min_bus_power_slack(pm::_PMs.AbstractPowerModel)
    return JuMP.@objective(pm.model, Min,
        sum(
            sum(
                sum( _PMs.var(pm, n, c, :p_slack, i)^2 + _PMs.var(pm, n, c, :q_slack, i)^2 for (i,bus) in nw_ref[:bus])
            for c in _PMs.conductor_ids(pm, n))
        for (n, nw_ref) in _PMs.nws(pm))
    )
end


"minimum load delta objective (continuous load shed)"
function objective_mc_min_load_delta(pm::_PMs.AbstractPowerModel)
    for (n, nw_ref) in _PMs.nws(pm)
        for c in _PMs.conductor_ids(pm, n)
            _PMs.var(pm, n, c)[:delta_pg] = JuMP.@variable(pm.model,
                [i in _PMs.ids(pm, n, :gen)], base_name="$(n)_$(c)_delta_pg",
                start = 0.0
            )
            for (i, gen) in nw_ref[:gen]
                JuMP.@constraint(pm.model, _PMs.var(pm, n, c, :delta_pg, i) >=  (gen["pg"][c] - _PMs.var(pm, n, c, :pg, i)))
                JuMP.@constraint(pm.model, _PMs.var(pm, n, c, :delta_pg, i) >= -(gen["pg"][c] - _PMs.var(pm, n, c, :pg, i)))
            end
        end
    end

    load_weight = Dict(n => Dict(i => get(load, "weight", 1.0) for (i,load) in _PMs.ref(pm, n, :load)) for n in _PMs.nw_ids(pm))
    M = Dict(n => Dict(c => 10*maximum([load_weight[n][i]*abs(load["pd"][c]) for (i,load) in _PMs.ref(pm, n, :load)]) for c in _PMs.conductor_ids(pm, n)) for n in _PMs.nw_ids(pm))

    JuMP.@objective(pm.model, Min,
        sum(
            sum(
                sum( (                 10*(1 - _PMs.var(pm, n, :z_voltage, i)) for i in keys(nw_ref[:bus]))) +
                sum( (            M[n][c]*(1 - _PMs.var(pm, n, :z_demand, i))) for i in keys(nw_ref[:load])) +
                sum( (abs(shunt["gs"][c])*(1 - _PMs.var(pm, n, :z_shunt, i))) for (i,shunt) in nw_ref[:shunt]) +
                sum( (                     _PMs.var(pm, n, c, :delta_pg, i) for i in keys(nw_ref[:gen])))
            for c in _PMs.conductor_ids(pm, n))
        for (n, nw_ref) in _PMs.nws(pm))
    )
end


"minimum load delta objective (continuous load shed) with storage"
function objective_mc_min_load_delta_strg(pm::_PMs.AbstractPowerModel)
    for (n, nw_ref) in _PMs.nws(pm)
        for c in _PMs.conductor_ids(pm, n)
            _PMs.var(pm, n, c)[:delta_pg] = JuMP.@variable(pm.model,
                [i in _PMs.ids(pm, n, :gen)], base_name="$(n)_$(c)_delta_pg",
                start = 0.0
            )
            for (i, gen) in nw_ref[:gen]
                JuMP.@constraint(pm.model, _PMs.var(pm, n, c, :delta_pg, i) >=  (gen["pg"][c] - _PMs.var(pm, n, c, :pg, i)))
                JuMP.@constraint(pm.model, _PMs.var(pm, n, c, :delta_pg, i) >= -(gen["pg"][c] - _PMs.var(pm, n, c, :pg, i)))
            end

            _PMs.var(pm, n, c)[:delta_ps] = JuMP.@variable(pm.model,
                [i in _PMs.ids(pm, n, :storage)], base_name="$(n)_$(c)_delta_ps",
                start = 0.0
            )
            for (i, strg) in nw_ref[:storage]
                JuMP.@constraint(pm.model, _PMs.var(pm, n, c, :delta_ps, i) >=  (strg["ps"][c] - _PMs.var(pm, n, c, :ps, i)))
                JuMP.@constraint(pm.model, _PMs.var(pm, n, c, :delta_ps, i) >= -(strg["ps"][c] - _PMs.var(pm, n, c, :ps, i)))
            end

        end
    end

    load_weight = Dict(n => Dict(i => get(load, "weight", 1.0) for (i,load) in _PMs.ref(pm, n, :load)) for n in _PMs.nw_ids(pm))
    M = Dict(n => Dict(c => 10*maximum([load_weight[n][i]*abs(load["pd"][c]) for (i,load) in _PMs.ref(pm, n, :load)]) for c in _PMs.conductor_ids(pm, n)) for n in _PMs.nw_ids(pm))

    JuMP.@objective(pm.model, Min,
        sum(
            sum(
                sum( (                 10*(1 - _PMs.var(pm, n, :z_voltage, i)) for i in keys(nw_ref[:bus]))) +
                sum( (            M[n][c]*(1 - _PMs.var(pm, n, :z_demand, i))) for i in keys(nw_ref[:load])) +
                sum( (abs(shunt["gs"][c])*(1 - _PMs.var(pm, n, :z_shunt, i))) for (i,shunt) in nw_ref[:shunt]) +
                sum( (                     _PMs.var(pm, n, c, :delta_pg, i) for i in keys(nw_ref[:gen]))) +
                sum( (                     _PMs.var(pm, n, c, :delta_ps, i) for i in keys(nw_ref[:storage])))
            for c in _PMs.conductor_ids(pm, n))
        for (n, nw_ref) in _PMs.nws(pm))
    )
end



"maximum loadability objective (continuous load shed) with storage"
function objective_mc_max_loadability_strg(pm::_PMs.AbstractPowerModel)
    load_weight = Dict(n => Dict(i => get(load, "weight", 1.0) for (i,load) in _PMs.ref(pm, n, :load)) for n in _PMs.nw_ids(pm))
    M = Dict(n => Dict(c => 10*maximum([load_weight[n][i]*abs(load["pd"][c]) for (i,load) in _PMs.ref(pm, n, :load)]) for c in _PMs.conductor_ids(pm, n)) for n in _PMs.nw_ids(pm))

    JuMP.@objective(pm.model, Max,
        sum(
            sum(
                sum( (              10*_PMs.var(pm, n, :z_voltage, i) for (i, bus) in nw_ref[:bus])) +
                sum( (         M[n][c]*_PMs.var(pm, n, :z_demand, i)) for (i,load) in nw_ref[:load]) +
                sum( (abs(shunt["gs"])*_PMs.var(pm, n, :z_shunt, i)) for (i,shunt) in nw_ref[:shunt]) +
                sum( (                 _PMs.var(pm, n, :z_gen, i) for (i,gen) in nw_ref[:gen])) +
                sum( (                 _PMs.var(pm, n, :z_storage, i) for (i,storage) in nw_ref[:storage]))
            for c in _PMs.conductor_ids(pm, n))
        for (n, nw_ref) in _PMs.nws(pm))
    )
end
