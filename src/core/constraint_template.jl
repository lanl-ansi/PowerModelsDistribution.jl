"reference angle constraints"
function constraint_mc_theta_ref(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    va_ref = _PMs.ref(pm, nw, :bus, i, "va")
    constraint_mc_theta_ref(pm, nw, i, va_ref)
end


""
function constraint_mc_power_balance_slack(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
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

    constraint_mc_power_balance_slack(pm, nw, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
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

    constraint_mc_ohms_yt_to(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end


""
function constraint_mc_model_voltage_magnitude_difference(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
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

    constraint_mc_model_voltage_magnitude_difference(pm, nw, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
end


""
function constraint_mc_model_current(pm::AbstractUBFModels; nw::Int=pm.cnw)
    for (i,branch) in _PMs.ref(pm, nw, :branch)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_idx = (i, f_bus, t_bus)

        g_sh_fr = branch["g_fr"]
        b_sh_fr = branch["b_fr"]

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

    r = branch["br_r"]
    x = branch["br_x"]
    g_sh_fr = branch["g_fr"]
    g_sh_to = branch["g_to"]
    b_sh_fr = branch["b_fr"]
    b_sh_to = branch["b_to"]

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
    config = trans["configuration"]
    type = trans["configuration"]
    f_cnd = trans["f_connections"][1:3]
    t_cnd = trans["t_connections"][1:3]
    tm_set = trans["tm"]
    tm_fixed = fix_taps ? ones(Bool, length(tm_set)) : trans["fixed"]
    tm_scale = calculate_tm_scale(trans, _PMs.ref(pm, nw, :bus, f_bus), _PMs.ref(pm, nw, :bus, t_bus))

    #TODO change data model
    # there is redundancy in specifying polarity seperately on from and to side
    #TODO change this once migrated to new data model
    pol = haskey(trans, "poalrity") ? trans["polarity"] : trans["config_fr"]["polarity"]

    if config=="wye"
        constraint_mc_trans_yy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
    elseif config=="delta"
        constraint_mc_trans_dy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_cnd, t_cnd, pol, tm_set, tm_fixed, tm_scale)
    end
    if type=="zig-zag"
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
        vm_ll_min = haskey(bus, "vm_ll_min") ? bus["vm_ll_min"] : fill(0, 3)
        vm_ll_max = haskey(bus, "vm_ll_max") ? bus["vm_ll_max"] : fill(Inf, 3)
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
function constraint_mc_load(pm::_PMs.AbstractPowerModel, id::Int; nw::Int=pm.cnw, report::Bool=true)
    load = _PMs.ref(pm, nw, :load, id)
    bus = _PMs.ref(pm, nw,:bus, load["load_bus"])

    conn = haskey(load, "conn") ? load["conn"] : "wye"

    a, alpha, b, beta = _load_expmodel_params(load, bus)

    if conn=="wye"
        constraint_mc_load_wye(pm, nw, id, load["load_bus"], a, alpha, b, beta)
    else
        constraint_mc_load_delta(pm, nw, id, load["load_bus"], a, alpha, b, beta)
    end
end


"""
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
function constraint_mc_generation(pm::_PMs.AbstractPowerModel, id::Int; nw::Int=pm.cnw, report::Bool=true, bounded::Bool=true)
    generator = _PMs.ref(pm, nw, :gen, id)
    bus = _PMs.ref(pm, nw,:bus, generator["gen_bus"])

    N = 3
    pmin = get(generator, "pmin", fill(-Inf, N))
    pmax = get(generator, "pmax", fill( Inf, N))
    qmin = get(generator, "qmin", fill(-Inf, N))
    qmax = get(generator, "qmax", fill( Inf, N))

    if generator["conn"]=="wye"
        constraint_mc_generation_wye(pm, nw, id, bus["index"], pmin, pmax, qmin, qmax; report=report, bounded=bounded)
    else
        constraint_mc_generation_delta(pm, nw, id, bus["index"], pmin, pmax, qmin, qmax; report=report, bounded=bounded)
    end
end


"KCL for load shed problem with transformers"
function constraint_mc_power_balance_shed(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
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

    constraint_mc_voltage_angle_difference(pm, nw, f_idx, branch["angmin"], branch["angmax"])
end


"storage loss constraints, delegate to PowerModels"
function constraint_mc_storage_loss(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw, kwargs...)
    storage = _PMs.ref(pm, nw, :storage, i)

    _PMs.constraint_storage_loss(pm, nw, i, storage["storage_bus"], storage["r"], storage["x"], storage["p_loss"], storage["q_loss"];
        conductors = _PMs.conductor_ids(pm, nw)
    )
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
    vmref = bus["vm"] #Not sure why this is needed
    constraint_mc_voltage_magnitude_setpoint(pm, nw, i, vmref)
end


"generator active power setpoint constraint"
function constraint_mc_active_gen_setpoint(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw, kwargs...)
    pg_set = _PMs.ref(pm, nw, :gen, i)["pg"]
    constraint_mc_active_gen_setpoint(pm, nw, i, pg_set)
end


"""
This constraint captures problem agnostic constraints that define limits for
voltage magnitudes (where variable bounds cannot be used)
Notable examples include IVRPowerModel and ACRPowerModel
"""
function constraint_mc_voltage_magnitude_bounds(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    bus = _PMs.ref(pm, nw, :bus, i)
    vmin = get(bus, "vmin", fill(0.0, 3)) #TODO update for four-wire
    vmax = get(bus, "vmax", fill(Inf, 3)) #TODO update for four-wire
    constraint_mc_voltage_magnitude_bounds(pm, nw, i, vmin, vmax)
end


function constraint_mc_generation_on_off(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    gen = _PMs.ref(pm, nw, :gen, i)

    constraint_mc_generation_on_off(pm, nw, i, gen["pmin"], gen["pmax"], gen["qmin"], gen["qmax"])
end

""
function constraint_mc_storage_thermal_limit(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    storage = _PMs.ref(pm, nw, :storage, i)
    constraint_mc_storage_thermal_limit(pm, nw, i, storage["thermal_rating"])
end

""
function constraint_mc_storage_current_limit(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    storage = _PMs.ref(pm, nw, :storage, i)
    constraint_mc_storage_current_limit(pm, nw, i, storage["storage_bus"], storage["current_rating"])
end

""
function constraint_mc_storage_on_off(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    storage = _PMs.ref(pm, nw, :storage, i)
    charge_ub = storage["charge_rating"]
    discharge_ub = storage["discharge_rating"]

    cnds = _PMs.conductor_ids(pm, nw)
    ncnds = length(cnds)
    pmin = zeros(ncnds)
    pmax = zeros(ncnds)
    qmin = zeros(ncnds)
    qmax = zeros(ncnds)

    for c in 1:ncnds
        inj_lb, inj_ub = _PMs.ref_calc_storage_injection_bounds(_PMs.ref(pm, nw, :storage), _PMs.ref(pm, nw, :bus), c)
        pmin[c] = inj_lb[i]
        pmax[c] = inj_ub[i]
        qmin[c] = max(inj_lb[i], _PMs.ref(pm, nw, :storage, i, "qmin")[c])
        qmax[c] = min(inj_ub[i], _PMs.ref(pm, nw, :storage, i, "qmax")[c])
    end

    constraint_mc_storage_on_off(pm, nw, i, pmin, pmax, qmin, qmax, charge_ub, discharge_ub)
end

"defines limits on active power output of a generator where bounds can't be used"
function constraint_mc_generation_active_power_limits(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    gen = _PMs.ref(pm, nw, :gen, i)
    bus = gen["gen_bus"]
    constraint_mc_generation_active_power_limits(pm, nw, i, bus, gen["pmax"], gen["pmin"])
end

"defines limits on reactive power output of a generator where bounds can't be used"
function constraint_mc_generation_reactive_power_limits(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    gen = _PMs.ref(pm, nw, :gen, i)
    bus = gen["gen_bus"]
    constraint_mc_generation_reactive_power_limits(pm, nw, i, bus, gen["qmax"], gen["qmin"])
end
""
function constraint_mc_current_balance_load(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
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

    constraint_mc_current_balance_load(pm, nw, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_gs, bus_bs)
end


""
function constraint_mc_current_from(pm::_PMs.AbstractIVRModel, i::Int; nw::Int=pm.cnw)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    tr, ti = _PMs.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    constraint_mc_current_from(pm, nw, f_bus, f_idx, g_fr, b_fr, tr, ti, tm)
end

""
function constraint_mc_current_to(pm::_PMs.AbstractIVRModel, i::Int; nw::Int=pm.cnw)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    tr, ti = _PMs.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    constraint_mc_current_to(pm, nw, t_bus, f_idx, t_idx, g_to, b_to)
end

""
function constraint_mc_voltage_drop(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    tr, ti = _PMs.calc_branch_t(branch)
    r = branch["br_r"]
    x = branch["br_x"]
    tm = branch["tap"]

    constraint_mc_voltage_drop(pm, nw, i, f_bus, t_bus, f_idx, r, x, tr, ti, tm)
end
