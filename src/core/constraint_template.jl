# Three-phase specific constraints


""
function constraint_ohms_yt_from(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, ph::Int=pm.cph)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = PMs.calc_branch_y(branch)
    g, b = g.values, b.values

    tr, ti = PMs.calc_branch_t(branch)
    tr, ti = tr.values, ti.values

    g_fr = branch["g_fr"].values
    b_fr = branch["b_fr"].values
    tm = branch["tap"].values

    constraint_ohms_yt_from(pm, nw, ph, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


""
function constraint_ohms_yt_to(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, ph::Int=pm.cph)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = PMs.calc_branch_y(branch)
    g, b = g.values, b.values

    tr, ti = PMs.calc_branch_t(branch)
    tr, ti = tr.values, ti.values

    g_to = branch["g_to"].values
    b_to = branch["b_to"].values
    tm = branch["tap"].values

    constraint_ohms_yt_to(pm, nw, ph, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end


""
function constraint_ohms_yt_from_on_off(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, ph::Int=pm.cph)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = PMs.calc_branch_y(branch)
    g, b = g.values, b.values

    tr, ti = PMs.calc_branch_t(branch)
    tr, ti = tr.values, ti.values

    g_fr = branch["g_fr"].values
    b_fr = branch["b_fr"].values
    tm = branch["tap"].values

    vad_min = [ref(pm, nw, :off_angmin, i) for i in PMs.phase_ids(pm)]
    vad_max = [ref(pm, nw, :off_angmax, i) for i in PMs.phase_ids(pm)]

    constraint_ohms_yt_from_on_off(pm, nw, ph, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
end


""
function constraint_ohms_yt_to_on_off(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, ph::Int=pm.cph)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = PMs.calc_branch_y(branch)
    g, b = g.values, b.values

    tr, ti = PMs.calc_branch_t(branch)
    tr, ti = tr.values, ti.values

    g_to = branch["g_to"].values
    b_to = branch["b_to"].values
    tm = branch["tap"].values

    vad_min = ref(pm, nw, :off_angmin, ph)
    vad_max = ref(pm, nw, :off_angmax, ph)

    constraint_ohms_yt_to_on_off(pm, nw, ph, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
end
