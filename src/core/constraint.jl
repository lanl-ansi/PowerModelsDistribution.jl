
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
function constraint_tp_load_power_setpoint(pm::GenericPowerModel, load_id::Int; nw=pm.cnw)
    load = ref(pm, nw, :load, load_id)
    pd = load["pd"]
    qd = load["qd"]
    conn = load["conn"]
    if conn=="wye"
        for c in PMs.conductor_ids(pm)
            constraint_load_power_setpoint_wye(pm, nw, c, load_id, pd[c], qd[c])
        end
    elseif conn=="delta"
        constraint_tp_load_power_setpoint_delta(pm, nw, load_id, load["load_bus"], pd, qd)
    else
        error(LOGGER, "Unknown load connection type $conn.")
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
function constraint_tp_load_power_prop_vm(pm::GenericPowerModel, load_id::Int; nw=pm.cnw)
    load = ref(pm, nw, :load, load_id)

    vnom_kv = load["vnom_kv"]
    vbase_kv_LL = ref(pm, nw, :bus, load["load_bus"])["base_kv"]
    vbase_kv_LN = vbase_kv_LL/sqrt(3)

    pd = load["pd"]
    qd = load["qd"]
    cp = pd/(vnom_kv/vbase_kv_LN)
    cq = qd/(vnom_kv/vbase_kv_LN)

    conn = load["conn"]
    if conn=="wye"
        for c in PMs.conductor_ids(pm)
            constraint_load_power_prop_vm_wye(pm, nw, c, load_id, load["load_bus"], cp[c], cq[c])
        end
    elseif conn=="delta"
        constraint_tp_load_power_prop_vm_delta(pm, nw, load_id, load["load_bus"], cp, cq)
    else
        error(LOGGER, "Unknown load connection type $conn.")
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
function constraint_tp_load_power_prop_vmsqr(pm::GenericPowerModel, load_id::Int; nw=pm.cnw)
    load = ref(pm, nw, :load, load_id)

    vnom_kv = load["vnom_kv"]
    vbase_kv_LL = ref(pm, nw, :bus, load["load_bus"])["base_kv"]
    vbase_kv_LN = vbase_kv_LL/sqrt(3)

    pd = load["pd"]
    qd = load["qd"]
    cp = pd/(vnom_kv/vbase_kv_LN)^2
    cq = qd/(vnom_kv/vbase_kv_LN)^2

    conn = load["conn"]
    if conn=="wye"
        for c in PMs.conductor_ids(pm)
            constraint_load_power_prop_vmsqr_wye(pm, nw, c, load_id, load["load_bus"], cp[c], cq[c])
        end
    elseif conn=="delta"
        constraint_tp_load_power_prop_vmsqr_delta(pm, nw, load_id, load["load_bus"], cp, cq)
    else
        error(LOGGER, "Unknown load connection type $conn.")
    end
end


"""
Sets va_start on every bus if the 'va_start' key is not present. This is needed
for delta loads, where division occurs by the difference of voltage phasors.
If the voltage phasors at one bus are initialized in the same point, this would
lead to division by zero.
"""
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
