# Three-phase specific constraints

""
function constraint_kcl_shunt_slack(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    if !haskey(con(pm, nw, cnd), :kcl_p)
        con(pm, nw, cnd)[:kcl_p] = Dict{Int,ConstraintRef}()
    end
    if !haskey(con(pm, nw, cnd), :kcl_q)
        con(pm, nw, cnd)[:kcl_q] = Dict{Int,ConstraintRef}()
    end

    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts = ref(pm, nw, :bus_shunts, i)

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd", cnd) for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd", cnd) for k in bus_loads)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs", cnd) for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs", cnd) for k in bus_shunts)

    constraint_kcl_shunt_slack(pm, nw, cnd, i, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
end


""
function constraint_tp_voltage(pm::GenericPowerModel; nw::Int=pm.cnw)
    for c in PMs.conductor_ids(pm)
        constraint_tp_voltage(pm, nw, c)
    end
end


"Set loose voltage bounds; they will only be set if the voltage magnitude was unconstrained before."
function constraint_tp_voltage_mag_unbound(pm::GenericPowerModel; nw::Int=pm.cnw, vmin=0.5, vmax=Inf)
    for bus_id in ids(pm, pm.cnw, :bus)
        constraint_tp_voltage_mag_unbound(pm, bus_id, vmin, vmax, nw=nw)
    end
end


""
function constraint_ohms_tp_yt_from(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd, atol_impzero=1E-13)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    tr, ti = PMs.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    if all(abs.(branch["br_r"]).<=atol_impzero) && all(abs.(branch["br_x"]).<=atol_impzero)
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        constraint_ohms_tp_yt_from_impzero(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g_fr, b_fr, g_to, b_to, tr, ti, tm)
    else
        g, b = PMs.calc_branch_y(branch)
        constraint_ohms_tp_yt_from(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    end
end


""
function constraint_ohms_tp_yt_to(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd, atol_impzero=1E-13)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = PMs.calc_branch_y(branch)
    tr, ti = PMs.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    if all(abs.(branch["br_r"]).<=atol_impzero) && all(abs.(branch["br_x"]).<=atol_impzero)
        # do nothing, constraint_ohms_tp_yt_from_impzero covers both sides
    else
        g, b = PMs.calc_branch_y(branch)
        constraint_ohms_tp_yt_to(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    end
end


""
function constraint_ohms_tp_yt_from_on_off(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = PMs.calc_branch_y(branch)
    tr, ti = PMs.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    vad_min = [ref(pm, nw, :off_angmin, c) for c in PMs.conductor_ids(pm)]
    vad_max = [ref(pm, nw, :off_angmax, c) for c in PMs.conductor_ids(pm)]

    constraint_ohms_tp_yt_from_on_off(pm, nw, cnd, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
end


""
function constraint_ohms_tp_yt_to_on_off(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = PMs.calc_branch_y(branch)
    tr, ti = PMs.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    vad_min = [ref(pm, nw, :off_angmin, c) for c in PMs.conductor_ids(pm)]
    vad_max = [ref(pm, nw, :off_angmax, c) for c in PMs.conductor_ids(pm)]

    constraint_ohms_tp_yt_to_on_off(pm, nw, cnd, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
end

""
function constraint_voltage_magnitude_difference(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = ref(pm, nw, :branch, i)
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

function constraint_tp_voltage_magnitude_difference(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw)
    branch = ref(pm, nw, :branch, i)
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


function constraint_tp_branch_current(pm::GenericPowerModel{T}, i::Int; nw::Int=pm.cnw) where T <: AbstractUBFForm
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    g_sh_fr = diagm(0 => branch["g_fr"].values)
    b_sh_fr = diagm(0 => branch["b_fr"].values)

    constraint_tp_branch_current(pm, nw, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
end


function constraint_flow_losses(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = ref(pm, nw, :branch, i)
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

    constraint_flow_losses(pm::GenericPowerModel, nw, cnd, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to, tm)
end

function constraint_tp_flow_losses(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw)
    branch = ref(pm, nw, :branch, i)
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

    constraint_tp_flow_losses(pm::GenericPowerModel, nw, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to)
end


""
function constraint_tp_storage_exchange(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw)
    storage = ref(pm, nw, :storage, i)

    PMs.constraint_storage_complementarity(pm, nw, i)
    constraint_tp_storage_loss(pm, nw, i, storage["storage_bus"], storage["r"], storage["x"], storage["standby_loss"])
end

""
function constraint_tp_trans_voltage(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw)
    (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to) = calc_tp_trans_Tvi(pm, i)
    # TODO add checks to simplify
    f_bus = ref(pm, :trans, i)["f_bus"]
    t_bus = ref(pm, :trans, i)["t_bus"]
    tapset = ref(pm, :trans, i)["tapset"]
    constraint_tp_trans_voltage(pm, i, f_bus, t_bus, tapset, Tv_fr, Tv_im, Cv_to)
end

""
function constraint_tp_trans_voltage_var(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw)
    (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to) = calc_tp_trans_Tvi(pm, i)
    # TODO add checks to simplify
    trans = ref(pm, :trans, i)
    f_bus = trans["f_bus"]
    t_bus = trans["t_bus"]
    constraint_tp_trans_voltage_var(pm, i, f_bus, t_bus, Tv_fr, Tv_im, Cv_to)
    # fix the taps with a constraint which are not free
    # at least one tap will be free
    tapfix = trans["tapfix"]
    tapset = trans["tapset"]
    constraint_tp_trans_tap_fix(pm, i, tapfix, tapset)
end

""
function constraint_tp_trans_flow(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw)
    (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to) = calc_tp_trans_Tvi(pm, i)
    # TODO add checks to simplify
    f_bus = ref(pm, :trans, i)["f_bus"]
    t_bus = ref(pm, :trans, i)["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    constraint_tp_trans_flow(pm, i, f_bus, t_bus, f_idx, t_idx, Ti_fr, Ti_im)
end

""
function constraint_tp_trans_vartap_fix(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw)
    tapset = ref(pm, nw, :trans, i)["tapset"]
    tapfix = ref(pm, nw, :trans, i)["tapfix"]
    constraint_tp_trans_vartap_fix(pm, i, tapset, tapfix, nw=nw)
end

""
function constraint_tp_trans_voltage_vartap(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw)
    f_bus = ref(pm, :trans, i)["f_bus"]
    t_bus = ref(pm, :trans, i)["t_bus"]
    bkv_fr = ref(pm, :bus, f_bus)["base_kv"]
    bkv_to = ref(pm, :bus, t_bus)["base_kv"]
    vnom_kv_fr = ref(pm, :trans, i)["vnom_kv"][1]
    vnom_kv_to = ref(pm, :trans, i)["vnom_kv"][2]
    constraint_tp_trans_voltage_vartap(pm, i, f_bus, t_bus, bkv_fr, bkv_to, vnom_kv_fr, vnom_kv_to)
end

""
function constraint_tp_trans_power_vartap(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw)
    f_bus = ref(pm, :trans, i)["f_bus"]
    t_bus = ref(pm, :trans, i)["t_bus"]
    bkv_fr = ref(pm, :bus, f_bus)["base_kv"]
    bkv_to = ref(pm, :bus, t_bus)["base_kv"]
    f_idx = (i, f_bus,  t_bus)
    t_idx = (i, t_bus,  f_bus)
    f_vnom_kv = ref(pm, :trans, i)["vnom_kv"][1]
    t_vnom_kv = ref(pm, :trans, i)["vnom_kv"][2]
    constraint_tp_trans_power_vartap(pm, i, f_bus, t_bus, f_idx, t_idx, bkv_fr, bkv_to, f_vnom_kv, t_vnom_kv)
end

""
function constraint_kcl_shunt_trans(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    if !haskey(con(pm, nw, cnd), :kcl_p)
        con(pm, nw, cnd)[:kcl_p] = Dict{Int,ConstraintRef}()
    end
    if !haskey(con(pm, nw, cnd), :kcl_q)
        con(pm, nw, cnd)[:kcl_q] = Dict{Int,ConstraintRef}()
    end

    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts = ref(pm, nw, :bus_shunts, i)
    bus_arcs_trans = ref(pm, nw, :bus_arcs_trans, i)

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd", cnd) for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd", cnd) for k in bus_loads)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs", cnd) for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs", cnd) for k in bus_shunts)

    constraint_kcl_shunt_trans(pm, nw, cnd, i, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
end

""
function constraint_kcl_shunt_trans(pm::GenericPowerModel, nw::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
    vm = var(pm, nw, c, :vm, i)
    p = var(pm, nw, c, :p)
    q = var(pm, nw, c, :q)
    pg = var(pm, nw, c, :pg)
    qg = var(pm, nw, c, :qg)
    p_dc = var(pm, nw, c, :p_dc)
    q_dc = var(pm, nw, c, :q_dc)
    p_trans = var(pm, nw, c, :p_trans)
    q_trans = var(pm,  nw, c, :q_trans)
    con(pm, nw, c, :kcl_p)[i] = @constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) + sum(p_trans[a_trans] for a_trans in bus_arcs_trans) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*vm^2)
    con(pm, nw, c, :kcl_q)[i] = @constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) + sum(q_trans[a_trans] for a_trans in bus_arcs_trans) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*vm^2)
end
