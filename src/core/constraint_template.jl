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
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]


    g_frout = Array{Float64}(length(g_fr))
    for i in range(1, length(g_fr))
        g_frout[i] = g_fr[i]
    end

    b_frout = Array{Float64}(length(b_fr))
    for i in range(1, length(b_fr))
        b_frout[i] = b_fr[i]
    end

    tmout = Array{Float64}(length(tm))
    for i in range(1, length(tm))
        tmout[i] = tm[i]
    end


    constraint_ohms_yt_from(pm, nw, ph, f_bus, t_bus, f_idx, t_idx, g, b, g_frout, b_frout, tr, ti, tmout)
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
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    g_toout = Array{Float64}(length(g_to))
    for i in range(1, length(g_to))
        g_toout[i] = g_to[i]
    end

    b_toout = Array{Float64}(length(b_to))
    for i in range(1, length(b_to))
        b_toout[i] = b_to[i]
    end

    tmout = Array{Float64}(length(tm))
    for i in range(1, length(tm))
        tmout[i] = tm[i]
    end

    constraint_ohms_yt_to(pm, nw, ph, f_bus, t_bus, f_idx, t_idx, g, b, g_toout, b_toout, tr, ti, tmout)
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
