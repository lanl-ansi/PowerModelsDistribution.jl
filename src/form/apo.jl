### generic features that apply to all active-power-only (apo) approximations

"do nothing, no reactive power in this model"
function variable_tp_trans_reactive_flow(pm::PMs.GenericPowerModel{T}; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded=true) where T <: PMs.AbstractActivePowerFormulation
end

""
function constraint_kcl_shunt_trans(pm::PMs.GenericPowerModel{T}, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: PMs.AbstractActivePowerFormulation
    pg   = PMs.var(pm, nw, c, :pg)
    p    = PMs.var(pm, nw, c, :p)
    p_dc = PMs.var(pm, nw, c, :p_dc)
    p_trans = PMs.var(pm, nw, c, :pt)

    PMs.con(pm, nw, c, :kcl_p)[i] = JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) + sum(p_trans[a_trans] for a_trans in bus_arcs_trans) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*1.0^2)
    # omit reactive constraint
end


function constraint_tp_storage_loss(pm::PMs.GenericPowerModel{T}, n::Int, i, bus, r, x, standby_loss) where T <: PMs.AbstractActivePowerFormulation
    conductors = PMs.conductor_ids(pm)
    ps = [PMs.var(pm, n, c, :ps, i) for c in conductors]
    sc = PMs.var(pm, n, :sc, i)
    sd = PMs.var(pm, n, :sd, i)

    JuMP.@NLconstraint(pm.model, sum(ps[c] for c in conductors) + (sd - sc) == standby_loss + sum( r[c]*ps[c]^2 for c in conductors) )
end
