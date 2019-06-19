# Three-phase specific constraints

""
function constraint_tp_power_balance_shunt_slack(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    if !haskey(_PMs.con(pm, nw, cnd), :kcl_p)
        _PMs.con(pm, nw, cnd)[:kcl_p] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(_PMs.con(pm, nw, cnd), :kcl_q)
        _PMs.con(pm, nw, cnd)[:kcl_q] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = _PMs.ref(pm, nw, :bus, i)
    bus_arcs = _PMs.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = _PMs.ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = _PMs.ref(pm, nw, :bus_gens, i)
    bus_loads = _PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PMs.ref(pm, nw, :bus_shunts, i)

    bus_pd = Dict(k => _PMs.ref(pm, nw, :load, k, "pd", cnd) for k in bus_loads)
    bus_qd = Dict(k => _PMs.ref(pm, nw, :load, k, "qd", cnd) for k in bus_loads)

    bus_gs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "gs", cnd) for k in bus_shunts)
    bus_bs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "bs", cnd) for k in bus_shunts)

    constraint_tp_power_balance_shunt_slack(pm, nw, cnd, i, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
end


""
function constraint_tp_model_voltage(pm::_PMs.GenericPowerModel; nw::Int=pm.cnw)
    for c in _PMs.conductor_ids(pm)
        constraint_tp_model_voltage(pm, nw, c)
    end
end


""
function constraint_tp_ohms_yt_from(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PMs.calc_branch_y(branch)
    tr, ti = _PMs.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    constraint_tp_ohms_yt_from(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


""
function constraint_tp_ohms_yt_to(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PMs.calc_branch_y(branch)
    tr, ti = _PMs.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    constraint_tp_ohms_yt_to(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end


""
function constraint_tp_ohms_yt_from_on_off(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PMs.calc_branch_y(branch)
    tr, ti = _PMs.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    vad_min = [_PMs.ref(pm, nw, :off_angmin, c) for c in _PMs.conductor_ids(pm)]
    vad_max = [_PMs.ref(pm, nw, :off_angmax, c) for c in _PMs.conductor_ids(pm)]

    constraint_tp_ohms_yt_from_on_off(pm, nw, cnd, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
end


""
function constraint_tp_ohms_yt_to_on_off(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PMs.calc_branch_y(branch)
    tr, ti = _PMs.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    vad_min = [_PMs.ref(pm, nw, :off_angmin, c) for c in _PMs.conductor_ids(pm)]
    vad_max = [_PMs.ref(pm, nw, :off_angmax, c) for c in _PMs.conductor_ids(pm)]

    constraint_tp_ohms_yt_to_on_off(pm, nw, cnd, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
end


""
function constraint_tp_voltage_magnitude_difference(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    r = branch["br_r"]
    x = branch["br_x"]
    g_sh_fr = branch["g_fr"]
    b_sh_fr = branch["b_fr"]
    tm = branch["tap"]

    constraint_tp_voltage_magnitude_difference(pm, nw, cnd, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
end

function constraint_tp_model_voltage_magnitude_difference(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    r = branch["br_r"].values
    x = branch["br_x"].values
    g_sh_fr = diagm(0 => branch["g_fr"].values)
    b_sh_fr = diagm(0 => branch["b_fr"].values)
    tm = branch["tap"].values

    constraint_tp_model_voltage_magnitude_difference(pm, nw, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
end


function constraint_tp_model_current(pm::_PMs.GenericPowerModel{T}; nw::Int=pm.cnw) where T <: AbstractUBFForm
    for (i,branch) in _PMs.ref(pm, nw, :branch)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_idx = (i, f_bus, t_bus)

        g_sh_fr = diagm(0 => branch["g_fr"].values)
        b_sh_fr = diagm(0 => branch["b_fr"].values)

        constraint_tp_model_current(pm, nw, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
    end
end

#= ??? DEPRECIATED ???
function constraint_tp_flow_losses(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = _PMs.ref(pm, nw, :branch, i)
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

    constraint_tp_flow_losses(pm::_PMs.GenericPowerModel, nw, cnd, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to, tm)
end
=#

function constraint_tp_flow_losses(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PMs.ref(pm, nw, :branch, i)
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

    constraint_tp_flow_losses(pm::_PMs.GenericPowerModel, nw, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to)
end


""
function constraint_tp_storage_exchange(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    storage = _PMs.ref(pm, nw, :storage, i)

    _PMs.constraint_storage_complementarity(pm, nw, i)
    constraint_tp_storage_loss(pm, nw, i, storage["storage_bus"], storage["r"], storage["x"], storage["standby_loss"])
end


function constraint_tp_trans(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    if _PMs.ref(pm, pm.cnw, :conductors)!=3
        Memento.error(_LOGGER, "Transformers only work with networks with three conductors.")
    end
    (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to) = _calc_tp_trans_Tvi(pm, i)
    f_bus = _PMs.ref(pm, :trans, i)["f_bus"]
    t_bus = _PMs.ref(pm, :trans, i)["t_bus"]
    tm = _PMs.ref(pm, :trans, i)["tm"]
    constraint_tp_trans_voltage(pm, nw, i, f_bus, t_bus, tm, Tv_fr, Tv_im, Cv_to)
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    constraint_tp_trans_flow(pm, nw, i, f_bus, t_bus, f_idx, t_idx, tm, Ti_fr, Ti_im, Cv_to)
end


function constraint_tp_oltc(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    if _PMs.ref(pm, pm.cnw, :conductors)!=3
        Memento.error(_LOGGER, "Transformers only work with networks with three conductors.")
    end
    (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to) = _calc_tp_trans_Tvi(pm, i)
    trans = _PMs.ref(pm, :trans, i)
    f_bus = trans["f_bus"]
    t_bus = trans["t_bus"]
    constraint_tp_oltc_voltage(pm, nw, i, f_bus, t_bus, Tv_fr, Tv_im, Cv_to)
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    constraint_tp_oltc_flow(pm, nw, i, f_bus, t_bus, f_idx, t_idx, Ti_fr, Ti_im, Cv_to)
    # fix the taps with a constraint which are not free
    trans = _PMs.ref(pm, :trans, i)
    fixed = trans["fixed"]
    tm = trans["tm"]
    constraint_tp_oltc_tap_fix(pm, i, fixed, tm)
end


"KCL including transformer arcs."
function constraint_tp_power_balance_shunt_trans(pm::_PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    if !haskey(_PMs.con(pm, nw, cnd), :kcl_p)
        _PMs.con(pm, nw, cnd)[:kcl_p] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(_PMs.con(pm, nw, cnd), :kcl_q)
        _PMs.con(pm, nw, cnd)[:kcl_q] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = _PMs.ref(pm, nw, :bus, i)
    bus_arcs = _PMs.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = _PMs.ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = _PMs.ref(pm, nw, :bus_gens, i)
    bus_loads = _PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PMs.ref(pm, nw, :bus_shunts, i)
    bus_arcs_trans = _PMs.ref(pm, nw, :bus_arcs_trans, i)

    bus_pd = Dict(k => _PMs.ref(pm, nw, :load, k, "pd", cnd) for k in bus_loads)
    bus_qd = Dict(k => _PMs.ref(pm, nw, :load, k, "qd", cnd) for k in bus_loads)

    bus_gs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "gs", cnd) for k in bus_shunts)
    bus_bs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "bs", cnd) for k in bus_shunts)

    constraint_tp_power_balance_shunt_trans(pm, nw, cnd, i, bus_arcs, bus_arcs_dc, bus_arcs_trans, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
end
