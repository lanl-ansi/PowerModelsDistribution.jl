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


"minimum load delta objective (continuous load shed)"
function objective_tp_min_load_delta(pm::_PMs.GenericPowerModel{T}) where T
    load_weight = Dict(n => Dict(i => get(load, "weight", 1.0) for (i,load) in _PMs.ref(pm, n, :load)) for n in _PMs.nw_ids(pm))
    M = Dict(n => Dict(c => 10*maximum([load_weight[n][i]*abs(load["pd"][c]) for (i,load) in _PMs.ref(pm, n, :load)]) for c in _PMs.conductor_ids(pm, n)) for n in _PMs.nw_ids(pm))

    JuMP.@objective(pm.model, Min,
        sum(
            sum(
                sum( (M[n][c]*10*(1 - _PMs.var(pm, n, :z_voltage, i)) for (i, bus) in nw_ref[:bus])) +
                sum( (load_weight[n][i]*load["pd"][c]*(1 - _PMs.var(pm, n, :z_demand, i))) for (i,load) in nw_ref[:load]) +
                sum( (M[n][c]*(1 - _PMs.var(pm, n, :z_shunt, i))) for (i,shunt) in nw_ref[:shunt]) +
                sum( (M[n][c]*(1 - _PMs.var(pm, n, :z_gen, i)) for (i,gen) in nw_ref[:gen]))
            for c in _PMs.conductor_ids(pm, n))
        for (n, nw_ref) in _PMs.nws(pm))
    )
end


"minimum load delta objective (continuous load shed) with storage"
function objective_tp_min_load_delta_strg(pm::_PMs.GenericPowerModel{T}) where T
    load_weight = Dict(n => Dict(i => get(load, "weight", 1.0) for (i,load) in _PMs.ref(pm, n, :load)) for n in _PMs.nw_ids(pm))
    M = Dict(n => Dict(c => 10*maximum([load_weight[n][i]*abs(load["pd"][c]) for (i,load) in _PMs.ref(pm, n, :load)]) for c in _PMs.conductor_ids(pm, n)) for n in _PMs.nw_ids(pm))

    JuMP.@objective(pm.model, Min,
        sum(
            sum(
                sum( (M[n][c]*10*(1 - _PMs.var(pm, n, :z_voltage, i)) for (i, bus) in nw_ref[:bus])) +
                sum( (load_weight[n][i]*load["pd"][c]*(1 - _PMs.var(pm, n, :z_demand, i))) for (i,load) in nw_ref[:load]) +
                sum( (M[n][c]*(1 - _PMs.var(pm, n, :z_shunt, i))) for (i,shunt) in nw_ref[:shunt]) +
                sum( (M[n][c]*(1 - _PMs.var(pm, n, :z_gen, i)) for (i,gen) in nw_ref[:gen])) +
                sum( (M[n][c]*(1 - _PMs.var(pm, n, :z_storage, i)) for (i,storage) in nw_ref[:storage]))
            for c in _PMs.conductor_ids(pm, n))
        for (n, nw_ref) in _PMs.nws(pm))
    )
end


"maximum loadability objective (continuous)"
function objective_tp_max_loadability(pm::_PMs.GenericPowerModel{T}) where T
    load_weight = Dict(n => Dict(i => get(load, "weight", 1.0) for (i,load) in _PMs.ref(pm, n, :load)) for n in _PMs.nw_ids(pm))
    M = Dict(n => Dict(c => 10*maximum([load_weight[n][i]*abs(load["pd"][c]) for (i,load) in _PMs.ref(pm, n, :load)]) for c in _PMs.conductor_ids(pm, n)) for n in _PMs.nw_ids(pm))

    JuMP.@objective(pm.model, Max,
        sum(
            sum(
                sum( (M[n][c]*10*_PMs.var(pm, n, :z_voltage, i) for (i, bus) in nw_ref[:bus])) +
                sum( (load_weight[n][i]*abs(load["pd"][c])*_PMs.var(pm, n, :z_demand, i)) for (i,load) in nw_ref[:load]) +
                sum( (M[n][c]*_PMs.var(pm, n, :z_shunt, i)) for (i,shunt) in nw_ref[:shunt]) +
                sum( (M[n][c]*_PMs.var(pm, n, :z_gen, i) for (i,gen) in nw_ref[:gen]))
            for c in _PMs.conductor_ids(pm, n))
        for (n, nw_ref) in _PMs.nws(pm))
    )
end


"maximum loadability objective (continuous) with storage"
function objective_tp_max_loadability_strg(pm::_PMs.GenericPowerModel{T}) where T
    load_weight = Dict(n => Dict(i => get(load, "weight", 1.0) for (i,load) in _PMs.ref(pm, n, :load)) for n in _PMs.nw_ids(pm))
    M = Dict(n => Dict(c => 10*maximum([load_weight[n][i]*abs(load["pd"][c]) for (i,load) in _PMs.ref(pm, n, :load)]) for c in _PMs.conductor_ids(pm, n)) for n in _PMs.nw_ids(pm))

    JuMP.@objective(pm.model, Max,
        sum(
            sum(
                sum( (M[n][c]*10*_PMs.var(pm, n, :z_voltage, i) for (i, bus) in nw_ref[:bus])) +
                sum( (load_weight[n][i]*abs(load["pd"][c])*_PMs.var(pm, n, :z_demand, i)) for (i,load) in nw_ref[:load]) +
                sum( (M[n][c]*_PMs.var(pm, n, :z_shunt, i)) for (i,shunt) in nw_ref[:shunt]) +
                sum( (M[n][c]*_PMs.var(pm, n, :z_gen, i) for (i,gen) in nw_ref[:gen])) +
                sum( (M[n][c]*_PMs.var(pm, n, :z_storage, i) for (i,storage) in nw_ref[:storage]))
            for c in _PMs.conductor_ids(pm, n))
        for (n, nw_ref) in _PMs.nws(pm))
    )
end
