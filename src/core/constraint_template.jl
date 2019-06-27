# Three-phase specific constraints

""
function constraint_kcl_shunt_slack(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    if !haskey(PMs.con(pm, nw, cnd), :kcl_p)
        PMs.con(pm, nw, cnd)[:kcl_p] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(PMs.con(pm, nw, cnd), :kcl_q)
        PMs.con(pm, nw, cnd)[:kcl_q] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = PMs.ref(pm, nw, :bus, i)
    bus_arcs = PMs.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = PMs.ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = PMs.ref(pm, nw, :bus_gens, i)
    bus_loads = PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts = PMs.ref(pm, nw, :bus_shunts, i)

    bus_pd = Dict(k => PMs.ref(pm, nw, :load, k, "pd", cnd) for k in bus_loads)
    bus_qd = Dict(k => PMs.ref(pm, nw, :load, k, "qd", cnd) for k in bus_loads)

    bus_gs = Dict(k => PMs.ref(pm, nw, :shunt, k, "gs", cnd) for k in bus_shunts)
    bus_bs = Dict(k => PMs.ref(pm, nw, :shunt, k, "bs", cnd) for k in bus_shunts)

    constraint_kcl_shunt_slack(pm, nw, cnd, i, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
end


""
function constraint_tp_voltage(pm::PMs.GenericPowerModel; nw::Int=pm.cnw)
    for c in PMs.conductor_ids(pm)
        constraint_tp_voltage(pm, nw, c)
    end
end


"delegate back to PowerModels by default"
function constraint_tp_voltage(pm::PMs.GenericPowerModel, n::Int, c::Int)
        PMs.constraint_voltage(pm, n, c)
end


""
function constraint_ohms_tp_yt_from(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = PMs.calc_branch_y(branch)
    tr, ti = PMs.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    constraint_ohms_tp_yt_from(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


""
function constraint_ohms_tp_yt_to(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = PMs.calc_branch_y(branch)
    tr, ti = PMs.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    constraint_ohms_tp_yt_to(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end


""
function constraint_ohms_tp_yt_from_on_off(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = PMs.calc_branch_y(branch)
    tr, ti = PMs.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    vad_min = [PMs.ref(pm, nw, :off_angmin, c) for c in PMs.conductor_ids(pm)]
    vad_max = [PMs.ref(pm, nw, :off_angmax, c) for c in PMs.conductor_ids(pm)]

    constraint_ohms_tp_yt_from_on_off(pm, nw, cnd, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
end


""
function constraint_ohms_tp_yt_to_on_off(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = PMs.calc_branch_y(branch)
    tr, ti = PMs.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    vad_min = [PMs.ref(pm, nw, :off_angmin, c) for c in PMs.conductor_ids(pm)]
    vad_max = [PMs.ref(pm, nw, :off_angmax, c) for c in PMs.conductor_ids(pm)]

    constraint_ohms_tp_yt_to_on_off(pm, nw, cnd, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
end


""
function constraint_voltage_magnitude_difference(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    r = branch["br_r"]
    x = branch["br_x"]
    g_sh_fr = branch["g_fr"]
    b_sh_fr = branch["b_fr"]
    tm = branch["tap"]

    constraint_voltage_magnitude_difference(pm, nw, cnd, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
end

function constraint_tp_voltage_magnitude_difference(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    branch = PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    r = branch["br_r"].values
    x = branch["br_x"].values
    g_sh_fr = diagm(0 => branch["g_fr"].values)
    b_sh_fr = diagm(0 => branch["b_fr"].values)
    tm = branch["tap"].values

    constraint_tp_voltage_magnitude_difference(pm, nw, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
end


function constraint_tp_branch_current(pm::PMs.GenericPowerModel{T}, i::Int; nw::Int=pm.cnw) where T <: AbstractUBFForm
    branch = PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    g_sh_fr = diagm(0 => branch["g_fr"].values)
    b_sh_fr = diagm(0 => branch["b_fr"].values)

    constraint_tp_branch_current(pm, nw, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
end


function constraint_flow_losses(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    r = branch["br_r"][cnd]
    x = branch["br_x"][cnd]
    tm = branch["tap"][cnd]
    g_sh_fr = branch["g_fr"][cnd]
    g_sh_to = branch["g_to"][cnd]
    b_sh_fr = branch["b_fr"][cnd]
    b_sh_to = branch["b_to"][cnd]

    constraint_flow_losses(pm::PMs.GenericPowerModel, nw, cnd, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to, tm)
end

function constraint_tp_flow_losses(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    branch = PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    r = branch["br_r"].values
    x = branch["br_x"].values
    g_sh_fr = diagm(0 => branch["g_fr"]).values
    g_sh_to = diagm(0 => branch["g_to"]).values
    b_sh_fr = diagm(0 => branch["b_fr"]).values
    b_sh_to = diagm(0 => branch["b_to"]).values

    constraint_tp_flow_losses(pm::PMs.GenericPowerModel, nw, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to)
end


""
function constraint_tp_storage_exchange(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    storage = PMs.ref(pm, nw, :storage, i)

    PMs.constraint_storage_complementarity(pm, nw, i)
    constraint_tp_storage_loss(pm, nw, i, storage["storage_bus"], storage["r"], storage["x"], storage["standby_loss"])
end


function constraint_tp_trans(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    if PMs.ref(pm, pm.cnw, :conductors)!=3
        Memento.error(LOGGER, "Transformers only work with networks with three conductors.")
    end
    (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to) = calc_tp_trans_Tvi(pm, i)
    f_bus = PMs.ref(pm, :trans, i)["f_bus"]
    t_bus = PMs.ref(pm, :trans, i)["t_bus"]
    tm = PMs.ref(pm, :trans, i)["tm"]
    constraint_tp_trans_voltage(pm, nw, i, f_bus, t_bus, tm, Tv_fr, Tv_im, Cv_to)
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    constraint_tp_trans_flow(pm, nw, i, f_bus, t_bus, f_idx, t_idx, tm, Ti_fr, Ti_im, Cv_to)
end


function constraint_tp_oltc(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    if PMs.ref(pm, pm.cnw, :conductors)!=3
        Memento.error(LOGGER, "Transformers only work with networks with three conductors.")
    end
    (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to) = calc_tp_trans_Tvi(pm, i)
    trans = PMs.ref(pm, :trans, i)
    f_bus = trans["f_bus"]
    t_bus = trans["t_bus"]
    constraint_tp_oltc_voltage(pm, nw, i, f_bus, t_bus, Tv_fr, Tv_im, Cv_to)
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    constraint_tp_oltc_flow(pm, nw, i, f_bus, t_bus, f_idx, t_idx, Ti_fr, Ti_im, Cv_to)
    # fix the taps with a constraint which are not free
    trans = PMs.ref(pm, :trans, i)
    fixed = trans["fixed"]
    tm = trans["tm"]
    constraint_tp_oltc_tap_fix(pm, i, fixed, tm)
end


"KCL including transformer arcs."
function constraint_kcl_shunt_trans(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    if !haskey(PMs.con(pm, nw, cnd), :kcl_p)
        PMs.con(pm, nw, cnd)[:kcl_p] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(PMs.con(pm, nw, cnd), :kcl_q)
        PMs.con(pm, nw, cnd)[:kcl_q] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = PMs.ref(pm, nw, :bus, i)
    bus_arcs = PMs.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = PMs.ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = PMs.ref(pm, nw, :bus_gens, i)
    bus_loads = PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts = PMs.ref(pm, nw, :bus_shunts, i)
    bus_arcs_trans = PMs.ref(pm, nw, :bus_arcs_trans, i)

    bus_pd = Dict(k => PMs.ref(pm, nw, :load, k, "pd", cnd) for k in bus_loads)
    bus_qd = Dict(k => PMs.ref(pm, nw, :load, k, "qd", cnd) for k in bus_loads)

    bus_gs = Dict(k => PMs.ref(pm, nw, :shunt, k, "gs", cnd) for k in bus_shunts)
    bus_bs = Dict(k => PMs.ref(pm, nw, :shunt, k, "bs", cnd) for k in bus_shunts)

    constraint_kcl_shunt_trans(pm, nw, cnd, i, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
end


"""
Impose all balance related constraints for which key present in data model of bus.
For a discussion of sequence components and voltage unbalance factor (VUF), see
@INPROCEEDINGS{girigoudar_molzahn_roald-2019,
	author={K. Girigoudar and D. K. Molzahn and L. A. Roald},
	booktitle={submitted},
	title={{Analytical and Empirical Comparisons of Voltage Unbalance Definitions}},
	year={2019},
	month={},
    url={https://molzahn.github.io/pubs/girigoudar_molzahn_roald-2019.pdf}
}
"""
function constraint_tp_voltage_balance(pm::PMs.GenericPowerModel, bus_id::Int; nw=pm.cnw)
    @assert(PMs.ref(pm, nw, :conductors)==3)

    bus = PMs.ref(pm, nw, :bus, bus_id)

    if haskey(bus, "vm_vuf_max")
        constraint_tp_vm_vuf(pm, nw, bus_id, bus["vm_vuf_max"])
    end

    if haskey(bus, "vm_seq_neg_max")
        constraint_tp_vm_neg_seq(pm, nw, bus_id, bus["vm_seq_neg_max"])
    end

    if haskey(bus, "vm_seq_pos_max")
        constraint_tp_vm_pos_seq(pm, nw, bus_id, bus["vm_seq_pos_max"])
    end

    if haskey(bus, "vm_seq_zero_max")
        constraint_tp_vm_zero_seq(pm, nw, bus_id, bus["vm_seq_zero_max"])
    end

    if haskey(bus, "vm_ll_min")|| haskey(bus, "vm_ll_max")
        vm_ll_min = haskey(bus, "vm_ll_min") ? bus["vm_ll_min"] : PMs.MultiConductorVector(fill(0, 3))
        vm_ll_max = haskey(bus, "vm_ll_max") ? bus["vm_ll_max"] : PMs.MultiConductorVector(fill(Inf, 3))
        constraint_tp_vm_ll(pm, nw, bus_id, vm_ll_min, vm_ll_max)
    end
end


"KCL including transformer arcs and load variables."
function constraint_kcl_shunt_trans_load(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    if !haskey(PMs.con(pm, nw, cnd), :kcl_p)
        PMs.con(pm, nw, cnd)[:kcl_p] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(PMs.con(pm, nw, cnd), :kcl_q)
        PMs.con(pm, nw, cnd)[:kcl_q] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = PMs.ref(pm, nw, :bus, i)
    bus_arcs = PMs.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = PMs.ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = PMs.ref(pm, nw, :bus_gens, i)
    bus_loads = PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts = PMs.ref(pm, nw, :bus_shunts, i)
    bus_arcs_trans = PMs.ref(pm, nw, :bus_arcs_trans, i)

    bus_gs = Dict(k => PMs.ref(pm, nw, :shunt, k, "gs", cnd) for k in bus_shunts)
    bus_bs = Dict(k => PMs.ref(pm, nw, :shunt, k, "bs", cnd) for k in bus_shunts)

    constraint_kcl_shunt_trans_load(pm, nw, cnd, i, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_loads, bus_gs, bus_bs)
end


"""
CONSTANT POWER
Fixes the load power sd.
sd = [sd_1, sd_2, sd_3]
What is actually fixed, depends on whether the load is connected in delta or wye.
When connected in wye, the load power equals the per-phase power sn drawn at the
bus to which the load is connected.
sd_1 = v_a.conj(i_a) = sn_a

CONSTANT CURRENT
Sets the active and reactive load power sd to be proportional to
the the voltage magnitude.
pd = cp.|vm|
qd = cq.|vm|
sd = cp.|vm| + j.cq.|vm|

CONSTANT IMPEDANCE
Sets the active and reactive power drawn by the load to be proportional to
the square of the voltage magnitude.
pd = cp.|vm|^2
qd = cq.|vm|^2
sd = cp.|vm|^2 + j.cq.|vm|^2

DELTA
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
function constraint_tp_load(pm::PMs.GenericPowerModel, id::Int; nw=pm.cnw)
    load = PMs.ref(pm, nw, :load, id)
    model = load["model"]
    conn = PMs.ref(pm, nw, :load, id, "conn")
    @assert(conn in ["delta", "wye"])

    if model=="constant_power"
        pd = load["pd"]
        qd = load["qd"]

        if conn=="wye"
            for c in PMs.conductor_ids(pm)
                constraint_load_power_wye(pm, nw, c, id, pd[c], qd[c])
            end
        elseif conn=="delta"
            @assert(PMs.ref(pm, 0, :conductors)==3)
            constraint_tp_load_power_delta(pm, nw, id, load["load_bus"], pd, qd)
        end

    elseif model=="constant_current"
        vnom_kv = load["vnom_kv"]
        vbase_kv_LL = PMs.ref(pm, nw, :bus, load["load_bus"])["base_kv"]
        vbase_kv_LN = vbase_kv_LL/sqrt(3)

        pd = load["pd"]
        qd = load["qd"]
        cp = pd/(vnom_kv/vbase_kv_LN)
        cq = qd/(vnom_kv/vbase_kv_LN)

        if conn=="wye"
            for c in PMs.conductor_ids(pm)
                constraint_load_current_wye(pm, nw, c, id, load["load_bus"], cp[c], cq[c])
            end
        elseif conn=="delta"
            @assert(PMs.ref(pm, 0, :conductors)==3)
            constraint_tp_load_current_delta(pm, nw, id, load["load_bus"], cp, cq)
        end

    elseif model=="constant_impedance"
        vnom_kv = load["vnom_kv"]
        vbase_kv_LL = PMs.ref(pm, nw, :bus, load["load_bus"])["base_kv"]
        vbase_kv_LN = vbase_kv_LL/sqrt(3)

        pd = load["pd"]
        qd = load["qd"]
        cp = pd/(vnom_kv/vbase_kv_LN)^2
        cq = qd/(vnom_kv/vbase_kv_LN)^2

        if conn=="wye"
            for c in PMs.conductor_ids(pm)
                constraint_load_impedance_wye(pm, nw, c, id, load["load_bus"], cp[c], cq[c])
            end
        elseif conn=="delta"
            @assert(PMs.ref(pm, 0, :conductors)==3)
            constraint_tp_load_impedance_delta(pm, nw, id, load["load_bus"], cp, cq)
        end

    else
        Memento.@error(LOGGER, "Unknown model $model for load $id.")
    end
end
