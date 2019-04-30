# Three-phase specific constraints


""
function constraint_tp_voltage(pm::PMs.GenericPowerModel{T}, n::Int, c::Int) where T <: PMs.AbstractWRForm
    w  = PMs.var(pm, n,  :w)
    wr = PMs.var(pm, n, :wr)
    wi = PMs.var(pm, n, :wi)

    for d in c:length(PMs.conductor_ids(pm))
        for (i,j) in PMs.ids(pm, n, :buspairs)
            InfrastructureModels.relaxation_complex_product(pm.model, w[(i,d)], w[(j,c)], wr[(i,j,c,d)], wi[(i,j,c,d)])
        end

        if d != c
            for i in PMs.ids(pm, n, :bus)
                InfrastructureModels.relaxation_complex_product(pm.model, w[(i,d)], w[(i,c)], wr[(i,i,c,d)], wi[(i,i,c,d)])
            end
        end
    end
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] ==        g/tm*w_fr[i] + (-g*tr+b*ti)/tm*(wr[i]) + (-b*tr-g*ti)/tm*(wi[i])
q[f_idx] == -(b+c/2)/tm*w_fr[i] - (-b*tr-g*ti)/tm*(wr[i]) + (-g*tr+b*ti)/tm*(wi[i])
```
"""
function constraint_ohms_tp_yt_from_on_off(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractWRForm
    p_fr = PMs.var(pm, n, c, :p, f_idx)
    q_fr = PMs.var(pm, n, c, :q, f_idx)
    w    = PMs.var(pm, n, :w)
    wr   = PMs.var(pm, n, :wr)
    wi   = PMs.var(pm, n, :wi)

    JuMP.@constraint(pm.model, p_fr ==  ( g_fr[c]+g[c,c]) * w[(f_bus, c)] +
                                sum( g[c,d] * wr[(f_bus, f_bus, c, d)] +
                                     b[c,d] * wi[(f_bus, f_bus, c, d)] for d in PMs.conductor_ids(pm) if d != c) +
                                sum(-g[c,d] * wr[(f_bus, t_bus, c, d)] +
                                    -b[c,d] * wi[(f_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm)) )
    JuMP.@constraint(pm.model, q_fr == -( b_fr[c]+b[c,c]) * w[(f_bus, c)] -
                                sum( b[c,d] * wr[(f_bus, f_bus, c, d)] -
                                     g[c,d] * wi[(f_bus, f_bus, c, d)] for d in PMs.conductor_ids(pm) if d != c) -
                                sum(-b[c,d] * wr[(f_bus, t_bus, c, d)] +
                                     g[c,d] * wi[(f_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm)) )
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==        g*w_to[i] + (-g*tr-b*ti)/tm*(wr[i]) + (-b*tr+g*ti)/tm*(-wi[i])
q[t_idx] == -(b+c/2)*w_to[i] - (-b*tr+g*ti)/tm*(wr[i]) + (-g*tr-b*ti)/tm*(-wi[i])
```
"""
function constraint_ohms_tp_yt_to_on_off(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractWRForm
    p_to = PMs.var(pm, n, c, :p, t_idx)
    q_to = PMs.var(pm, n, c, :q, t_idx)
    w    = PMs.var(pm, n, :w)
    wr   = PMs.var(pm, n, :wr)
    wi   = PMs.var(pm, n, :wi)

    JuMP.@constraint(pm.model, p_to ==  ( g_to[c]+g[c,c]) * w[(t_bus, c)] +
                                sum( g[c,d] * wr[(t_bus, t_bus, c, d)] +
                                     b[c,d] *-wi[(t_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm) if d != c) +
                                sum(-g[c,d] * wr[(f_bus, t_bus, c, d)] +
                                    -b[c,d] *-wi[(f_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm)) )
    JuMP.@constraint(pm.model, q_to == -( b_to[c]+b[c,c]) * w[(t_bus, c)] -
                                sum( b[c,d] * wr[(t_bus, t_bus, c, d)] -
                                     g[c,d] *-wi[(t_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm) if d != c) -
                                sum(-b[c,d] * wr[(f_bus, t_bus, c, d)] +
                                     g[c,d] *-wi[(f_bus, t_bus, c, d)] for d in PMs.conductor_ids(pm)) )
end


function constraint_kcl_shunt_trans(pm::PMs.GenericPowerModel{T}, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: PMs.AbstractWRForm
    w    = PMs.var(pm, nw, c, :w, i)
    pg   = PMs.var(pm, nw, c, :pg)
    qg   = PMs.var(pm, nw, c, :qg)
    p    = PMs.var(pm, nw, c, :p)
    q    = PMs.var(pm, nw, c, :q)
    p_dc = PMs.var(pm, nw, c, :p_dc)
    q_dc = PMs.var(pm, nw, c, :q_dc)
    p_trans = PMs.var(pm, nw, c, :pt)
    q_trans = PMs.var(pm,  nw, c, :qt)

    JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) + sum(p_trans[a_trans] for a_trans in bus_arcs_trans) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*w)
    JuMP.@constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) + sum(q_trans[a_trans] for a_trans in bus_arcs_trans) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*w)
end
