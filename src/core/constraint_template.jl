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
    for c in _PMs.conductor_ids(pm; nw=nw)
        constraint_mc_model_voltage(pm, nw, c)
    end
end


"delegate back to PowerModels by default"
function constraint_mc_model_voltage(pm::_PMs.AbstractPowerModel, n::Int, c::Int)
    _PMs.constraint_model_voltage(pm, n, c)
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

    for cnd in _PMs.conductor_ids(pm)
        constraint_mc_ohms_yt_from(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    end
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

    for cnd in _PMs.conductor_ids(pm; nw=nw)
        constraint_mc_ohms_yt_to(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    end
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
function constraint_mc_model_current(pm::AbstractUBFModels; nw::Int=pm.cnw)
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

    constraint_mc_flow_losses(pm::_PMs.AbstractPowerModel, nw, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to)
end


""
function constraint_mc_trans(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    if _PMs.ref(pm, pm.cnw, :conductors)!=3
        Memento.error(_LOGGER, "Transformers only work with networks with three conductors.")
    end
    (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to) = _calc_mc_transformer_Tvi(pm, i)
    f_bus = _PMs.ref(pm, :transformer, i)["f_bus"]
    t_bus = _PMs.ref(pm, :transformer, i)["t_bus"]
    tm = _PMs.ref(pm, :transformer, i)["tm"]
    constraint_mc_transformer_voltage(pm, nw, i, f_bus, t_bus, tm, Tv_fr, Tv_im, Cv_to)
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    constraint_mc_transformer_flow(pm, nw, i, f_bus, t_bus, f_idx, t_idx, tm, Ti_fr, Ti_im, Cv_to)
end


""
function constraint_mc_oltc(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    if _PMs.ref(pm, pm.cnw, :conductors)!=3
        Memento.error(_LOGGER, "Transformers only work with networks with three conductors.")
    end
    (Tv_fr,Tv_im,Ti_fr,Ti_im,Cv_to) = _calc_mc_transformer_Tvi(pm, i)
    trans = _PMs.ref(pm, :transformer, i)
    f_bus = trans["f_bus"]
    t_bus = trans["t_bus"]
    constraint_mc_oltc_voltage(pm, nw, i, f_bus, t_bus, Tv_fr, Tv_im, Cv_to)
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    constraint_mc_oltc_flow(pm, nw, i, f_bus, t_bus, f_idx, t_idx, Ti_fr, Ti_im, Cv_to)
    # fix the taps with a constraint which are not free
    trans = _PMs.ref(pm, :transformer, i)
    fixed = trans["fixed"]
    tm = trans["tm"]
    constraint_mc_oltc_tap_fix(pm, i, fixed, tm)
end


"KCL including transformer arcs"
function constraint_mc_power_balance(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
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

        constraint_mc_power_balance(pm, nw, cnd, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    end
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
        vm_ll_min = haskey(bus, "vm_ll_min") ? bus["vm_ll_min"] : _PMs.MultiConductorVector(fill(0, 3))
        vm_ll_max = haskey(bus, "vm_ll_max") ? bus["vm_ll_max"] : _PMs.MultiConductorVector(fill(Inf, 3))
        constraint_mc_vm_ll(pm, nw, bus_id, vm_ll_min, vm_ll_max)
    end
end


"KCL including transformer arcs and load variables."
function constraint_mc_power_balance_load(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
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

        bus_gs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "gs", cnd) for k in bus_shunts)
        bus_bs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "bs", cnd) for k in bus_shunts)

        constraint_mc_power_balance_load(pm, nw, cnd, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_gs, bus_bs)
    end
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
                constraint_load_power_wye(pm, nw, c, id, pd[c], qd[c])
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
        a, α, b, β = _load_expmodel_params(load, bus)

        if conn=="wye"
            for c in _PMs.conductor_ids(pm)
                constraint_load_exponential_wye(pm, nw, c, id, load["load_bus"], a[c], α[c], b[c], β[c])
            end
        elseif conn=="delta"
            @assert(_PMs.ref(pm, 0, :conductors)==3)
            constraint_mc_load_exponential_delta(pm, nw, id, load["load_bus"], a, α, b, β)
        end

    else
        Memento.@error(LOGGER, "Unknown model $model for load $id.")
    end
end


"KCL for load shed problem with transformers"
function constraint_mc_power_balance_shed(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
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

        constraint_mc_power_balance_shed(pm, nw, cnd, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    end
end


"on/off constraint for bus voltages"
function constraint_mc_bus_voltage_on_off(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm; nw=nw)
        constraint_mc_bus_voltage_on_off(pm, nw, c; kwargs...)
    end
end


"on/off voltage magnitude constraint"
function constraint_mc_voltage_magnitude_on_off(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    bus = _PMs.ref(pm, nw, :bus, i)

    constraint_mc_voltage_magnitude_on_off(pm, nw, cnd, i, bus["vmin"][cnd], bus["vmax"][cnd])
end


"on/off voltage magnitude squared constraint for relaxed formulations"
function constraint_mc_voltage_magnitude_sqr_on_off(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    bus = _PMs.ref(pm, nw, :bus, i)

    constraint_mc_voltage_magnitude_sqr_on_off(pm, nw, cnd, i, bus["vmin"][cnd], bus["vmax"][cnd])
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
    _PMs.constraint_storage_loss(pm, i; conductors=_PMs.conductor_ids(pm; nw=nw), nw=nw, kwargs...)
end


"storage thermal limit constraints, delegate to PowerModels per conductor"
function constraint_mc_storage_thermal_limit(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm; nw=nw)
        _PMs.constraint_storage_thermal_limit(pm, i; cnd=c, nw=nw, kwargs...)
    end
end


"branch thermal constraints from, delegate to PowerModels per conductor"
function constraint_mc_thermal_limit_from(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm; nw=nw)
        _PMs.constraint_thermal_limit_from(pm, i; cnd=c, nw=nw, kwargs...)
    end
end


"branch thermal constraints to, delegate to PowerModels per conductor"
function constraint_mc_thermal_limit_to(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm; nw=nw)
        _PMs.constraint_thermal_limit_to(pm, i; cnd=c, nw=nw, kwargs...)
    end
end


"voltage magnitude setpoint constraint, delegate to PowerModels per conductor"
function constraint_mc_voltage_magnitude_setpoint(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm; nw=nw)
        _PMs.constraint_voltage_magnitude_setpoint(pm, i; nw=nw, cnd=c, kwargs...)
    end
end


"generator active power setpoint constraint, delegate to PowerModels"
function constraint_mc_active_gen_setpoint(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm; nw=nw)
        _PMs.constraint_active_gen_setpoint(pm, i; nw=nw, cnd=c, kwargs...)
    end
end
