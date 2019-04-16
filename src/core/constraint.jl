
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


function constraint_tp_load_flow_setpoint(pm::GenericPowerModel, load_id::Int; nw=pm.cnw)
    load = ref(pm, nw, :load, load_id)
    pd = load["pd"]
    qd = load["qd"]
    conn = load["conn"]
    if conn=="wye"
        for c in PMs.conductor_ids(pm)
            constraint_load_flow_setpoint_wye(pm, nw, c, load_id, pd[c], qd[c])
        end
    elseif conn=="delta"
        constraint_tp_load_flow_setpoint_delta(pm, nw, load_id, load["load_bus"], pd, qd)
    else
        error(LOGGER, "Unknown load connection type $conn.")
    end
end


"Sets va_starts on every bus if the 'va_start' key is not present."
function set_tp_va_start_if_unset(pm::GenericPowerModel; nw::Int=pm.cnw, va1=0.0, va2=-2*pi/3, va3=2*pi/3, va_offset=0.0)
    va1 = va1+va_offset
    va2 = va2+va_offset
    va3 = va3+va_offset

    for bus_id in ids(pm, pm.cnw, :bus)
        if !haskey(ref(pm, nw, :bus, bus_id), "va_start")
            ref(pm, nw, :bus, bus_id)["va_start"] = MultiConductorVector([va1,va2, va3])
        end
    end
end
