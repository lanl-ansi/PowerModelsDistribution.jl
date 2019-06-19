
function constraint_tp_model_current(pm::_PMs.GenericPowerModel; kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.constraint_model_current(pm; cnd=c, kwargs...)
    end
end


function constraint_tp_theta_ref(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    for cnd in _PMs.conductor_ids(pm)
        constraint_tp_theta_ref(pm, nw, cnd, i)
    end
end


function constraint_tp_storage_loss(pm::_PMs.GenericPowerModel, n::Int, i, bus, r, x, standby_loss)
    conductors = _PMs.conductor_ids(pm)
    vm = [_PMs.var(pm, n, c, :vm, bus) for c in conductors]
    ps = [_PMs.var(pm, n, c, :ps, i) for c in conductors]
    qs = [_PMs.var(pm, n, c, :qs, i) for c in conductors]
    sc = _PMs.var(pm, n, :sc, i)
    sd = _PMs.var(pm, n, :sd, i)

    JuMP.@NLconstraint(pm.model, sum(ps[c] for c in conductors) + (sd - sc) == standby_loss + sum( r[c]*(ps[c]^2 + qs[c]^2)/vm[c]^2 for c in conductors))
end

