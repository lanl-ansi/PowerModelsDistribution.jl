# Three-phase specific constraints


""
function constraint_ohms_yt_from(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, ph::Int=pm.cph)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_fr = branch["g_fr"][ph]
    b_fr = branch["b_fr"][ph]
    tm = branch["tap"][ph]

    constraint_ohms_yt_from(pm, nw, ph, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


""
function constraint_ohms_yt_to(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, ph::Int=pm.cph)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_to = branch["g_to"][ph]
    b_to = branch["b_to"][ph]
    tm = branch["tap"][ph]

    constraint_ohms_yt_to(pm, nw, ph, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end


""
function constraint_ohms_y_from(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, ph::Int=pm.cph)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    g_fr = branch["g_fr"][ph]
    b_fr = branch["b_fr"][ph]
    tm = branch["tap"][ph]
    ta = branch["shift"][ph]

    constraint_ohms_y_from(pm, nw, ph, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tm, ta)
end


""
function constraint_ohms_y_to(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, ph::Int=pm.cph)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    g_to = branch["g_to"][ph]
    b_to = branch["b_to"][ph]
    tm = branch["tap"][ph]
    ta = branch["shift"][ph]

    constraint_ohms_y_to(pm, nw, ph, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tm, ta)
end


""
function constraint_ohms_yt_from_on_off(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, ph::Int=pm.cph)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_fr = branch["g_fr"][ph]
    b_fr = branch["b_fr"][ph]
    tm = branch["tap"][ph]

    vad_min = ref(pm, nw, :off_angmin, ph)
    vad_max = ref(pm, nw, :off_angmax, ph)

    constraint_ohms_yt_from_on_off(pm, nw, ph, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
end


""
function constraint_ohms_yt_to_on_off(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, ph::Int=pm.cph)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_to = branch["g_to"][ph]
    b_to = branch["b_to"][ph]
    tm = branch["tap"][ph]

    vad_min = ref(pm, nw, :off_angmin, ph)
    vad_max = ref(pm, nw, :off_angmax, ph)

    constraint_ohms_yt_to_on_off(pm, nw, ph, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
end


""
function constraint_ohms_yt_from_ne(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, ph::Int=pm.cph)
    branch = ref(pm, nw, :ne_branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_fr = branch["g_fr"][ph]
    b_fr = branch["b_fr"][ph]
    tm = branch["tap"][ph]

    vad_min = ref(pm, nw, :off_angmin, ph)
    vad_max = ref(pm, nw, :off_angmax, ph)

    constraint_ohms_yt_from_ne(pm, nw, ph, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
end


""
function constraint_ohms_yt_to_ne(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, ph::Int=pm.cph)
    branch = ref(pm, nw, :ne_branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_to = branch["g_to"][ph]
    b_to = branch["b_to"][ph]
    tm = branch["tap"][ph]

    vad_min = ref(pm, nw, :off_angmin, ph)
    vad_max = ref(pm, nw, :off_angmax, ph)

    constraint_ohms_yt_to_ne(pm, nw, ph, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
end