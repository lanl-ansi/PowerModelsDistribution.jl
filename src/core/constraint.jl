
function constraint_tp_branch_current(pm::GenericPowerModel, i::Int; kwargs...)
    for c in PMs.conductor_ids(pm)
        PMs.constraint_branch_current(pm, i, cnd=c; kwargs...)
    end
end


function constraint_tp_theta_ref(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw)
    for cnd in PMs.conductor_ids(pm)
        constraint_tp_theta_ref(pm, nw, cnd, i)
    end
end


function constraint_tp_storage_loss(pm::GenericPowerModel, n::Int, i, bus, r, x, standby_loss)
    conductors = PMs.conductor_ids(pm)
    vm = [var(pm, n, c, :vm, bus) for c in conductors]
    ps = [var(pm, n, c, :ps, i) for c in conductors]
    qs = [var(pm, n, c, :qs, i) for c in conductors]
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)

    @NLconstraint(pm.model, sum(ps[c] for c in conductors) + (sd - sc) == standby_loss + sum( r[c]*(ps[c]^2 + qs[c]^2)/vm[c]^2 for c in conductors))
end

function constraint_tp_vuf(pm::GenericPowerModel; nw::Int=pm.cnw, default=Inf)
    for id in PMs.ids(pm, nw, :bus)
        bus = ref(pm, nw, :bus, id)
        if haskey(bus, "vufmax")
            constraint_tp_vuf(pm, nw, id, bus["vufmax"])
        elseif default < Inf
            constraint_tp_vuf(pm, nw, id, default)
        end
    end
end

function constraint_tp_vmneg(pm::GenericPowerModel; nw::Int=pm.cnw, default=Inf)
    for id in PMs.ids(pm, nw, :bus)
        bus = ref(pm, nw, :bus, id)
        if haskey(bus, "vmnegmax")
            constraint_tp_vmneg(pm, nw, id, bus["vmnegmax"])
        elseif default < Inf
            constraint_tp_vmneg(pm, nw, id, default)
        end
    end
end

function constraint_tp_vmpos(pm::GenericPowerModel; nw::Int=pm.cnw, default=Inf)
    for id in PMs.ids(pm, nw, :bus)
        bus = ref(pm, nw, :bus, id)
        if haskey(bus, "vmnegmax")
            constraint_tp_vmneg(pm, nw, id, bus["vmposmax"])
        elseif default < Inf
            constraint_tp_vmneg(pm, nw, id, default)
        end
    end
end

function constraint_tp_vmzero(pm::GenericPowerModel; nw::Int=pm.cnw, default=Inf)
    for id in PMs.ids(pm, nw, :bus)
        bus = ref(pm, nw, :bus, id)
        if haskey(bus, "vmnegmax")
            constraint_tp_vmneg(pm, nw, id, bus["vmzeromax"])
        elseif default < Inf
            constraint_tp_vmneg(pm, nw, id, default)
        end
    end
end
