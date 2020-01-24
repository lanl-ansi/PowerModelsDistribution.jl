"reference angle constraints"
function constraint_mc_theta_ref(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    constraint_mc_theta_ref(pm, nw, i)
end


""
function constraint_mc_power_balance_slack(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    for cnd in _PMs.conductor_ids(pm; nw=nw)
        if !haskey(_PMs.con(pm, nw, cnd), :kcl_p)
            _PMs.con(pm, nw, cnd)[:kcl_p] = Dict{Int,JuMP.ConstraintRef}()
        end
        if !haskey(_PMs.con(pm, nw, cnd), :kcl_q)
            _PMs.con(pm, nw, cnd)[:kcl_q] = Dict{Int,JuMP.ConstraintRef}()
        end

        bus = _PMs.ref(pm, nw, :bus, i)
        bus_arcs = _PMs.ref(pm, nw, :bus_arcs, i)
        bus_arcs_sw = _PMs.ref(pm, nw, :bus_arcs_sw, i)
        bus_arcs_trans = _PMs.ref(pm, nw, :bus_arcs_trans, i)
        bus_gens = _PMs.ref(pm, nw, :bus_gens, i)
        bus_storage = _PMs.ref(pm, nw, :bus_storage, i)
        bus_loads = _PMs.ref(pm, nw, :bus_loads, i)
        bus_shunts = _PMs.ref(pm, nw, :bus_shunts, i)

        bus_pd = Dict(k => _PMs.ref(pm, nw, :load, k, "pd", cnd) for k in bus_loads)
        bus_qd = Dict(k => _PMs.ref(pm, nw, :load, k, "qd", cnd) for k in bus_loads)

        bus_gs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "gs", cnd) for k in bus_shunts)
        bus_bs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "bs", cnd) for k in bus_shunts)

        constraint_mc_power_balance_slack(pm, nw, cnd, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    end
end


""
function constraint_mc_model_voltage(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    constraint_mc_model_voltage(pm, nw)
end


"ohms constraint for branches on the from-side"
function constraint_mc_ohms_yt_from(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
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

    #TODO why was this not required before?
    if length(size(g_fr)) == 1
        tmp = MultiConductorMatrix(0.0, length(g_fr))
        for c in 1:length(g_fr)
            tmp[c,c] = g_fr[c]
        end
        g_fr = tmp
    end

    if length(size(b_fr)) == 1
        tmp = MultiConductorMatrix(0.0, length(b_fr))
        for c in 1:length(b_fr)
            tmp[c,c] = b_fr[c]
        end
        b_fr = tmp
    end

    constraint_mc_ohms_yt_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


"ohms constraint for branches on the to-side"
function constraint_mc_ohms_yt_to(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
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

    #TODO why was this not required before?
    if length(size(g_to)) == 1
        tmp = MultiConductorMatrix(0.0, length(g_to))
        for c in 1:length(g_to)
            tmp[c,c] = g_to[c]
        end
        g_to = tmp
    end

    if length(size(b_to)) == 1
        tmp = MultiConductorMatrix(0.0, length(b_to))
        for c in 1:length(b_to)
            tmp[c,c] = b_to[c]
        end
        b_to = tmp
    end

    constraint_mc_ohms_yt_to(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end


""
function constraint_mc_model_voltage_magnitude_difference(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    r = branch["br_r"].values
    x = branch["br_x"].values
    g_sh_fr = branch["g_fr"].values
    b_sh_fr = branch["b_fr"].values
    tm = branch["tap"].values

    constraint_mc_model_voltage_magnitude_difference(pm, nw, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
end


""
function constraint_mc_model_current(pm::Union{AbstractUBFModels,LPLinUBFPowerModel}; nw::Int=pm.cnw)
    for (i,branch) in _PMs.ref(pm, nw, :branch)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_idx = (i, f_bus, t_bus)

        g_sh_fr = branch["g_fr"].values
        b_sh_fr = branch["b_fr"].values

        constraint_mc_model_current(pm, nw, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
    end
end


""
function constraint_mc_flow_losses(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    r = branch["br_r"].values
    x = branch["br_x"].values
    g_sh_fr = branch["g_fr"].values
    g_sh_to = branch["g_to"].values
    b_sh_fr = branch["b_fr"].values
    b_sh_to = branch["b_to"].values

    tm = [1, 1, 1] #TODO

    constraint_mc_flow_losses(pm, nw, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to, tm)
end


"Transformer constraints, considering winding type, conductor order, polarity and tap settings."
function constraint_mc_trans(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw, fix_taps::Bool=true)
    if _PMs.ref(pm, pm.cnw, :conductors)!=3
        Memento.error(_LOGGER, "Transformers only work with networks with three conductors.")
    end

    trans = _PMs.ref(pm, :transformer, i)
    f_bus = trans["f_bus"]
    t_bus = trans["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    f_type = trans["config_fr"]["type"]
    t_type = trans["config_to"]["type"]
    f_cnd = trans["config_fr"]["cnd"]
    t_cnd = trans["config_to"]["cnd"]
    tm_set = trans["tm"]
    tm_fixed = fix_taps ? ones(Bool, length(tm_set)) : trans["fixed"]
    tm_scale = calculate_tm_scale(trans, _PMs.ref(pm, nw, :bus, f_bus), _PMs.ref(pm, nw, :bus, t_bus))

    #TODO change data model
    # there is redundancy in specifying polarity seperately on from and to side
    f_pol = trans["config_fr"]["polarity"]=='+' ? 1 : -1
    t_pol = trans["config_to"]["polarity"]=='+' ? 1 : -1
    pol = f_pol*t_pol

    if f_type=="wye" && t_type=="wye"
        constraint_mc_trans_yy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
    elseif f_type=="delta" && t_type=="wye"
        constraint_mc_trans_dy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
    elseif f_type=="wye" && t_type=="delta"
        constraint_mc_trans_dy(pm, nw, i, t_bus, f_bus, t_idx, f_idx, t_cnd, f_cnd, pol, tm_set, tm_fixed, (tm_scale)^-1)
    elseif f_type=="delta" && t_type=="delta"
        Memento.error(_LOGGER, "Dd transformers are not supported at the low-level data format. This can be cast as a combo of two dy transformers.")
    end
    if f_type=="zig-zag" || t_type=="zig-zag"
        Memento.error(_LOGGER, "Zig-zag not yet supported.")
    end
end


"KCL including transformer arcs"
function constraint_mc_power_balance(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    bus = _PMs.ref(pm, nw, :bus, i)
    bus_arcs = _PMs.ref(pm, nw, :bus_arcs, i)
    bus_arcs_sw = _PMs.ref(pm, nw, :bus_arcs_sw, i)
    bus_arcs_trans = _PMs.ref(pm, nw, :bus_arcs_trans, i)
    bus_gens = _PMs.ref(pm, nw, :bus_gens, i)
    bus_storage = _PMs.ref(pm, nw, :bus_storage, i)
    bus_loads = _PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PMs.ref(pm, nw, :bus_shunts, i)

    bus_pd = Dict(k => _PMs.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => _PMs.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_mc_power_balance(pm, nw, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
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
function constraint_mc_voltage_balance(pm::_PMs.AbstractPowerModel, bus_id::Int; nw=pm.cnw)
    @assert(_PMs.ref(pm, nw, :conductors)==3)

    bus = _PMs.ref(pm, nw, :bus, bus_id)

    if haskey(bus, "vm_vuf_max")
        constraint_mc_vm_vuf(pm, nw, bus_id, bus["vm_vuf_max"])
    end

    if haskey(bus, "vm_seq_neg_max")
        constraint_mc_vm_neg_seq(pm, nw, bus_id, bus["vm_seq_neg_max"])
    end

    if haskey(bus, "vm_seq_pos_max")
        constraint_mc_vm_pos_seq(pm, nw, bus_id, bus["vm_seq_pos_max"])
    end

    if haskey(bus, "vm_seq_zero_max")
        constraint_mc_vm_zero_seq(pm, nw, bus_id, bus["vm_seq_zero_max"])
    end

    if haskey(bus, "vm_ll_min")|| haskey(bus, "vm_ll_max")
        vm_ll_min = haskey(bus, "vm_ll_min") ? bus["vm_ll_min"] : MultiConductorVector(fill(0, 3))
        vm_ll_max = haskey(bus, "vm_ll_max") ? bus["vm_ll_max"] : MultiConductorVector(fill(Inf, 3))
        constraint_mc_vm_ll(pm, nw, bus_id, vm_ll_min, vm_ll_max)
    end
end


"KCL including transformer arcs and load variables."
function constraint_mc_power_balance_load(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    bus = _PMs.ref(pm, nw, :bus, i)
    bus_arcs = _PMs.ref(pm, nw, :bus_arcs, i)
    bus_arcs_sw = _PMs.ref(pm, nw, :bus_arcs_sw, i)
    bus_arcs_trans = _PMs.ref(pm, nw, :bus_arcs_trans, i)
    bus_gens = _PMs.ref(pm, nw, :bus_gens, i)
    bus_storage = _PMs.ref(pm, nw, :bus_storage, i)
    bus_loads = _PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PMs.ref(pm, nw, :bus_shunts, i)

    bus_gs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_mc_power_balance_load(pm, nw, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_gs, bus_bs)
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
function constraint_mc_load(pm::_PMs.AbstractPowerModel, id::Int; nw::Int=pm.cnw)
    load = _PMs.ref(pm, nw, :load, id)
    bus = _PMs.ref(pm, nw,:bus, load["load_bus"])
    model = load["model"]
    conn = _PMs.ref(pm, nw, :load, id, "conn")
    @assert(conn in ["delta", "wye"])

    if model=="constant_power"
        pd = load["pd"]
        qd = load["qd"]

        if conn=="wye"
            for c in _PMs.conductor_ids(pm; nw=nw)
                constraint_mc_load_power_wye(pm, nw, id, pd, qd)
            end
        elseif conn=="delta"
            @assert(_PMs.ref(pm, 0, :conductors)==3)
            constraint_mc_load_power_delta(pm, nw, id, load["load_bus"], pd, qd)
        end

    elseif model=="constant_current"
        vnom_kv = load["vnom_kv"]
        vbase_kv_LL = _PMs.ref(pm, nw, :bus, load["load_bus"])["base_kv"]
        vbase_kv_LN = vbase_kv_LL/sqrt(3)

        pd = load["pd"]
        qd = load["qd"]
        cp = pd/(vnom_kv/vbase_kv_LN)
        cq = qd/(vnom_kv/vbase_kv_LN)

        if conn=="wye"
            for c in _PMs.conductor_ids(pm; nw=nw)
                constraint_load_current_wye(pm, nw, c, id, load["load_bus"], cp[c], cq[c])
            end
        elseif conn=="delta"
            @assert(_PMs.ref(pm, 0, :conductors)==3)
            constraint_mc_load_current_delta(pm, nw, id, load["load_bus"], cp, cq)
        end

    elseif model=="constant_impedance"
        vnom_kv = load["vnom_kv"]
        vbase_kv_LL = _PMs.ref(pm, nw, :bus, load["load_bus"])["base_kv"]
        vbase_kv_LN = vbase_kv_LL/sqrt(3)

        pd = load["pd"]
        qd = load["qd"]
        cp = pd/(vnom_kv/vbase_kv_LN)^2
        cq = qd/(vnom_kv/vbase_kv_LN)^2

        if conn=="wye"
            for c in _PMs.conductor_ids(pm; nw=nw)
                constraint_load_impedance_wye(pm, nw, c, id, load["load_bus"], cp[c], cq[c])
            end
        elseif conn=="delta"
            @assert(_PMs.ref(pm, 0, :conductors)==3)
            constraint_mc_load_impedance_delta(pm, nw, id, load["load_bus"], cp, cq)
        end
    elseif model=="exponential"
        a, alpha, b, beta = _load_expmodel_params(load, bus)

        if conn=="wye"
            for c in _PMs.conductor_ids(pm)
                constraint_load_exponential_wye(pm, nw, c, id, load["load_bus"], a[c], alpha[c], b[c], beta[c])
            end
        elseif conn=="delta"
            @assert(_PMs.ref(pm, 0, :conductors)==3)
            constraint_mc_load_exponential_delta(pm, nw, id, load["load_bus"], a, alpha, b, beta)
        end

    else
        Memento.@error(_LOGGER, "Unknown model $model for load $id.")
    end
end


"KCL for load shed problem with transformers"
function constraint_mc_power_balance_shed(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    for cnd in _PMs.conductor_ids(pm; nw=nw)
        if !haskey(_PMs.con(pm, nw), :kcl_p)
            _PMs.con(pm, nw)[:kcl_p] = Dict{Int,Array{JuMP.ConstraintRef,1}}()
        end
        if !haskey(_PMs.con(pm, nw), :kcl_q)
            _PMs.con(pm, nw)[:kcl_q] = Dict{Int,Array{JuMP.ConstraintRef,1}}()
        end

        bus = _PMs.ref(pm, nw, :bus, i)
        bus_arcs = _PMs.ref(pm, nw, :bus_arcs, i)
        bus_arcs_sw = _PMs.ref(pm, nw, :bus_arcs_sw, i)
        bus_arcs_trans = _PMs.ref(pm, nw, :bus_arcs_trans, i)
        bus_gens = _PMs.ref(pm, nw, :bus_gens, i)
        bus_storage = _PMs.ref(pm, nw, :bus_storage, i)
        bus_loads = _PMs.ref(pm, nw, :bus_loads, i)
        bus_shunts = _PMs.ref(pm, nw, :bus_shunts, i)

        bus_pd = Dict(k => _PMs.ref(pm, nw, :load, k, "pd") for k in bus_loads)
        bus_qd = Dict(k => _PMs.ref(pm, nw, :load, k, "qd") for k in bus_loads)

        bus_gs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
        bus_bs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

        constraint_mc_power_balance_shed(pm, nw, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    end
end


"on/off constraint for bus voltages"
function constraint_mc_bus_voltage_on_off(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, kwargs...)
    constraint_mc_bus_voltage_on_off(pm, nw; kwargs...)
end


"on/off voltage magnitude constraint"
function constraint_mc_voltage_magnitude_on_off(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    bus = _PMs.ref(pm, nw, :bus, i)

    constraint_mc_voltage_magnitude_on_off(pm, nw, i, bus["vmin"], bus["vmax"])
end


"on/off voltage magnitude squared constraint for relaxed formulations"
function constraint_mc_voltage_magnitude_sqr_on_off(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    bus = _PMs.ref(pm, nw, :bus, i)

    constraint_mc_voltage_magnitude_sqr_on_off(pm, nw, i, bus["vmin"], bus["vmax"])
end


"This is duplicated at PMD level to correctly handle the indexing of the shunts."
function constraint_mc_voltage_angle_difference(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    pair = (f_bus, t_bus)
    buspair = _PMs.ref(pm, nw, :buspairs, pair)

    if buspair["branch"] == i
        constraint_mc_voltage_angle_difference(pm, nw, f_idx, buspair["angmin"], buspair["angmax"])
    end
end


"storage loss constraints, delegate to PowerModels"
function constraint_mc_storage_loss(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw, kwargs...)
    storage = ref(pm, nw, :storage, i)

    _PMs.constraint_storage_loss(pm, nw, i, storage["storage_bus"], storage["r"], storage["x"], storage["p_loss"], storage["q_loss"];
        conductors = _PMs.conductor_ids(pm, nw)
    )
end


"storage thermal limit constraints, delegate to PowerModels per conductor"
function constraint_mc_storage_thermal_limit(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm; nw=nw)
        _PMs.constraint_storage_thermal_limit(pm, i; cnd=c, nw=nw, kwargs...)
    end
end


"branch thermal constraints from"
function constraint_mc_thermal_limit_from(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if haskey(branch, "rate_a")
        constraint_mc_thermal_limit_from(pm, nw, f_idx, branch["rate_a"])
    end
end


"branch thermal constraints to"
function constraint_mc_thermal_limit_to(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    if haskey(branch, "rate_a")
        constraint_mc_thermal_limit_to(pm, nw, t_idx, branch["rate_a"])
    end
end


"voltage magnitude setpoint constraint"
function constraint_mc_voltage_magnitude_setpoint(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw, kwargs...)
    bus = _PMs.ref(pm, nw, :bus, i)
    vmref = bus["vm"].values #Not sure why this is needed
    constraint_mc_voltage_magnitude_setpoint(pm, nw, i, vmref)
end


"generator active power setpoint constraint"
function constraint_mc_active_gen_setpoint(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw, kwargs...)
    pg_set = _PMs.ref(pm, nw, :gen, i)["pg"].values
    constraint_mc_active_gen_setpoint(pm, nw, i, pg_set)
end


"""
This constraint captures problem agnostic constraints that define limits for
voltage magnitudes (where variable bounds cannot be used)
Notable examples include IVRPowerModel and ACRPowerModel
"""
function constraint_mc_voltage_magnitude_bounds(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    bus = _PMs.ref(pm, nw, :bus, i)
    constraint_mc_voltage_magnitude_bounds(pm, nw, i, bus["vmin"], bus["vmax"])
end


function constraint_mc_generation_on_off(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    gen = _PMs.ref(pm, nw, :gen, i)

    constraint_mc_generation_on_off(pm, nw, i, gen["pmin"], gen["pmax"], gen["qmin"], gen["qmax"])
end
