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


"""
In the general case, do not add the bounds;
expand this on a per formulation case as deemed necessary.
"""
function constraint_tp_voltage_mag_unbound(pm::GenericPowerModel{T}, i::Int, vmin::Float64, vmax::Float64; nw::Int=pm.cnw) where T
    # do nothing
end


""
function constraint_ohms_tp_yt_from(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
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

    constraint_ohms_tp_yt_from(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


""
function constraint_ohms_tp_yt_to(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
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

    constraint_ohms_tp_yt_to(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
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


"Links the voltage at both windings of a fixed tap transformer."
function constraint_tp_trans_voltage(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw)
    (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to) = calc_tp_trans_Tvi(pm, i)
    # TODO add checks to simplify
    f_bus = ref(pm, :trans, i)["f_bus"]
    t_bus = ref(pm, :trans, i)["t_bus"]
    tapset = ref(pm, :trans, i)["tapset"]
    constraint_tp_trans_voltage(pm, i, f_bus, t_bus, tapset, Tv_fr, Tv_im, Cv_to)
end


"Links the voltage at both windings of a variable tap transformer."
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


"Links the power flowing into both windings of a fixed tap transformer."
function constraint_tp_trans_flow(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw)
    (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to) = calc_tp_trans_Tvi(pm, i)
    # TODO add checks to simplify
    f_bus = ref(pm, :trans, i)["f_bus"]
    t_bus = ref(pm, :trans, i)["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    constraint_tp_trans_flow(pm, i, f_bus, t_bus, f_idx, t_idx, Ti_fr, Ti_im)
end


"Links the power flowing into both windings of a variable tap transformer."
function constraint_tp_trans_flow_var(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw)
    (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to) = calc_tp_trans_Tvi(pm, i)
    # TODO add checks to simplify
    f_bus = ref(pm, :trans, i)["f_bus"]
    t_bus = ref(pm, :trans, i)["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    constraint_tp_trans_flow(pm, i, f_bus, t_bus, f_idx, t_idx, Ti_fr, Ti_im)
end


"KCL including transformer arcs."
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


"KCL including transformer arcs and load variables."
function constraint_kcl_shunt_trans_load(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
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

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs", cnd) for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs", cnd) for k in bus_shunts)

    constraint_kcl_shunt_trans_load(pm, nw, cnd, i, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_loads, bus_gs, bus_bs)
end
