""
function constraint_mc_model_voltage(pm::_PMs.AbstractWRModel, n::Int, c::Int)
    w  = _PMs.var(pm, n,  :w)
    wr = _PMs.var(pm, n, :wr)
    wi = _PMs.var(pm, n, :wi)

    for d in c:length(_PMs.conductor_ids(pm))
        for (i,j) in _PMs.ids(pm, n, :buspairs)
            InfrastructureModels.relaxation_complex_product(pm.model, w[(i,d)], w[(j,c)], wr[(i,j,c,d)], wi[(i,j,c,d)])
        end

        if d != c
            for i in _PMs.ids(pm, n, :bus)
                InfrastructureModels.relaxation_complex_product(pm.model, w[(i,d)], w[(i,c)], wr[(i,i,c,d)], wi[(i,i,c,d)])
            end
        end
    end
end


"power balance constraint with line shunts and transformers for relaxed WR forms"
function constraint_mc_power_balance(pm::_PMs.AbstractWRModel, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    w    = _PMs.var(pm, nw, c, :w, i)
    p    = get(_PMs.var(pm, nw, c),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(_PMs.var(pm, nw, c),    :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(_PMs.var(pm, nw, c),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(_PMs.var(pm, nw, c),   :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PMs.var(pm, nw, c),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PMs.var(pm, nw, c),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(_PMs.var(pm, nw, c),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(_PMs.var(pm, nw, c),  :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt   = get(_PMs.var(pm, nw, c),   :pt, Dict()); _PMs._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt   = get(_PMs.var(pm, nw, c),   :qt, Dict()); _PMs._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")


    _PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        - sum(gs for gs in values(bus_gs))*w
    )
    _PMs.con(pm, nw, c, :kcl_q)[i] = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(qt[a_trans] for a_trans in bus_arcs_trans)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd for qd in values(bus_qd))
        + sum(bs for bs in values(bus_bs))*w
    )
end
