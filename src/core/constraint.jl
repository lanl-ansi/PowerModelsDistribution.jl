
function constraint_tp_branch_current(pm::PMs.GenericPowerModel, i::Int; kwargs...)
    for c in PMs.conductor_ids(pm)
        PMs.constraint_branch_current(pm, i, cnd=c; kwargs...)
    end
end


function constraint_tp_theta_ref(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    for cnd in PMs.conductor_ids(pm)
        constraint_tp_theta_ref(pm, nw, cnd, i)
    end
end


function constraint_tp_storage_loss(pm::PMs.GenericPowerModel, n::Int, i, bus, r, x, standby_loss)
    conductors = PMs.conductor_ids(pm)
    vm = [PMs.var(pm, n, c, :vm, bus) for c in conductors]
    ps = [PMs.var(pm, n, c, :ps, i) for c in conductors]
    qs = [PMs.var(pm, n, c, :qs, i) for c in conductors]
    sc = PMs.var(pm, n, :sc, i)
    sd = PMs.var(pm, n, :sd, i)

    JuMP.@NLconstraint(pm.model, sum(ps[c] for c in conductors) + (sd - sc) == standby_loss + sum( r[c]*(ps[c]^2 + qs[c]^2)/vm[c]^2 for c in conductors))
end
