
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


"""
Also referred to as a constant power load.
Fixes the load power sd.
sd = [sd_1, sd_2, sd_3]
What is actually fixed, depends on whether the load is connected in delta or wye.
When connected in wye, the load power equals the per-phase power sn drawn at the
bus to which the load is connected.
sd_1 = v_a.conj(i_a) = sn_a
When connected in delta, the load power gives the reference in the delta reference
frame. This means
sd_1 = v_ab.conj(i_ab) = (v_a-v_b).conj(i_ab)
We can relate this to the per-phase power by
sn_a = v_a.conj(i_a)
    = v_a.conj(i_ab-i_ca)
    = v_a.conj(conj(s_ab/v_ab) - conj(s_ca/v_ca))
    = v_a.(s_ab/(v_a-v_b) - s_ca/(v_c-v_a))
So for delta, sn is constrained indirectly.
"""
function constraint_tp_load_power(pm::PMs.GenericPowerModel, load_id::Int; nw=pm.cnw)
    load = PMs.ref(pm, nw, :load, load_id)
    pd = load["pd"]
    qd = load["qd"]
    conn = load["conn"]
    if conn=="wye"
        for c in PMs.conductor_ids(pm)
            constraint_load_power_wye(pm, nw, c, load_id, pd[c], qd[c])
        end
    elseif conn=="delta"
        @assert(PMs.ref(pm, 0, :conductors)==3)
        constraint_tp_load_power_delta(pm, nw, load_id, load["load_bus"], pd, qd)
    else
        Memento.error(LOGGER, "Unknown load connection type $conn.")
    end
end


"""
Also referred to as a constant current load.
Sets the active and reactive load power sd to be proportional to
the the voltage magnitude.
pd = cp.|vm|
qd = cq.|vm|
sd = cp.|vm| + j.cq.|vm|
The same remark applies on delta/wye as for the fixed setpoint.
"""
function constraint_tp_load_current(pm::PMs.GenericPowerModel, load_id::Int; nw=pm.cnw)
    load = PMs.ref(pm, nw, :load, load_id)

    vnom_kv = load["vnom_kv"]
    vbase_kv_LL = PMs.ref(pm, nw, :bus, load["load_bus"])["base_kv"]
    vbase_kv_LN = vbase_kv_LL/sqrt(3)

    pd = load["pd"]
    qd = load["qd"]
    cp = pd/(vnom_kv/vbase_kv_LN)
    cq = qd/(vnom_kv/vbase_kv_LN)

    conn = load["conn"]
    if conn=="wye"
        for c in PMs.conductor_ids(pm)
            constraint_load_current_wye(pm, nw, c, load_id, load["load_bus"], cp[c], cq[c])
        end
    elseif conn=="delta"
        @assert(PMs.ref(pm, 0, :conductors)==3)
        constraint_tp_load_current_delta(pm, nw, load_id, load["load_bus"], cp, cq)
    else
        Memento.error(LOGGER, "Unknown load connection type $conn.")
    end
end


"""
Also referred to as a constant impedance load.
Sets the active and reactive power drawn by the load to be proportional to
the square of the voltage magnitude.
pd = cp.|vm|^2
qd = cq.|vm|^2
sd = cp.|vm|^2 + j.cq.|vm|^2
The same remark applies on delta/wye as for the fixed setpoint.
"""
function constraint_tp_load_impedance(pm::PMs.GenericPowerModel, load_id::Int; nw=pm.cnw)
    load = PMs.ref(pm, nw, :load, load_id)

    vnom_kv = load["vnom_kv"]
    vbase_kv_LL = PMs.ref(pm, nw, :bus, load["load_bus"])["base_kv"]
    vbase_kv_LN = vbase_kv_LL/sqrt(3)

    pd = load["pd"]
    qd = load["qd"]
    cp = pd/(vnom_kv/vbase_kv_LN)^2
    cq = qd/(vnom_kv/vbase_kv_LN)^2

    conn = load["conn"]
    if conn=="wye"
        for c in PMs.conductor_ids(pm)
            constraint_load_impedance_wye(pm, nw, c, load_id, load["load_bus"], cp[c], cq[c])
        end
    elseif conn=="delta"
        @assert(PMs.ref(pm, 0, :conductors)==3)
        constraint_tp_load_impedance_delta(pm, nw, load_id, load["load_bus"], cp, cq)
    else
        Memento.error(LOGGER, "Unknown load connection type $conn.")
    end
end
