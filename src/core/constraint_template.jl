# Voltage constraints

""
function constraint_mc_model_voltage(pm::AbstractMCPowerModel; nw::Int=nw_id_default)
    constraint_mc_model_voltage(pm, nw)
end


"reference angle constraints"
function constraint_mc_theta_ref(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    va_ref = ref(pm, nw, :bus, i, "va")
    constraint_mc_theta_ref(pm, nw, i, va_ref)
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
function constraint_mc_bus_voltage_balance(pm::AbstractMCPowerModel, bus_id::Int; nw=nw_id_default)
    @assert(ref(pm, nw, :conductors)==3)

    bus = ref(pm, nw, :bus, bus_id)

    if haskey(bus, "vm_vuf_max")
        constraint_mc_bus_voltage_magnitude_vuf(pm, nw, bus_id, bus["vm_vuf_max"])
    end

    if haskey(bus, "vm_seq_neg_max")
        constraint_mc_bus_voltage_magnitude_negative_sequence(pm, nw, bus_id, bus["vm_seq_neg_max"])
    end

    if haskey(bus, "vm_seq_pos_max")
        constraint_mc_bus_voltage_magnitude_positive_sequence(pm, nw, bus_id, bus["vm_seq_pos_max"])
    end

    if haskey(bus, "vm_seq_zero_max")
        constraint_mc_bus_voltage_magnitude_zero_sequence(pm, nw, bus_id, bus["vm_seq_zero_max"])
    end

    if haskey(bus, "vm_ll_min")|| haskey(bus, "vm_ll_max")
        vm_ll_min = haskey(bus, "vm_ll_min") ? bus["vm_ll_min"] : fill(0, 3)
        vm_ll_max = haskey(bus, "vm_ll_max") ? bus["vm_ll_max"] : fill(Inf, 3)
        constraint_mc_bus_voltage_magnitude_ll(pm, nw, bus_id, vm_ll_min, vm_ll_max)
    end
end


"voltage magnitude setpoint constraint"
function constraint_mc_voltage_magnitude_only(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default, kwargs...)
    bus = ref(pm, nw, :bus, i)
    if haskey(bus, "vm")
        constraint_mc_voltage_magnitude_only(pm, nw, i, bus["vm"])
    end
end


"""
This constraint captures problem agnostic constraints that define limits for
voltage magnitudes (where variable bounds cannot be used)
Notable examples include IVRPowerModel and ACRPowerModel
"""
function constraint_mc_voltage_magnitude_bounds(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    vmin = get(bus, "vmin", fill(0.0, 3)) #TODO update for four-wire
    vmax = get(bus, "vmax", fill(Inf, 3)) #TODO update for four-wire
    constraint_mc_voltage_magnitude_bounds(pm, nw, i, vmin, vmax)
end

## Voltage on/off constraints

"on/off constraint for bus voltages"
function constraint_mc_bus_voltage_on_off(pm::AbstractMCPowerModel; nw::Int=nw_id_default, kwargs...)
    constraint_mc_bus_voltage_on_off(pm, nw; kwargs...)
end


"on/off voltage magnitude constraint"
function constraint_mc_bus_voltage_magnitude_on_off(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)

    constraint_mc_bus_voltage_magnitude_on_off(pm, nw, i, bus["vmin"], bus["vmax"])
end


"on/off voltage magnitude squared constraint for relaxed formulations"
function constraint_mc_bus_voltage_magnitude_sqr_on_off(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)

    constraint_mc_bus_voltage_magnitude_sqr_on_off(pm, nw, i, bus["vmin"], bus["vmax"])
end


# Switch constraints

""
function constraint_mc_switch_state(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    switch = ref(pm, nw, :switch, i)
    f_bus = switch["f_bus"]
    t_bus = switch["t_bus"]

    f_idx = (i, f_bus, t_bus)

    if switch["state"] != 0
        constraint_mc_switch_state_closed(pm, nw, f_bus, t_bus, switch["f_connections"], switch["t_connections"])
    else
        constraint_mc_switch_state_open(pm, nw, f_idx)
    end
end


function constraint_mc_switch_thermal_limit(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    switch = ref(pm, nw, :switch, i)

    if haskey(switch, "thermal_rating")
        f_idx = (i, switch["f_bus"], switch["t_bus"])
        constraint_mc_switch_thermal_limit(pm, nw, f_idx, switch["fr_connections"], switch["thermal_rating"])
    end
end


function constraint_mc_switch_current_limit(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    switch = ref(pm, nw, :switch, i)

    if haskey(switch, "current_rating")
        f_idx = (i, switch["f_bus"], switch["t_bus"])
        constraint_mc_switch_current_limit(pm, nw, f_idx, switch["fr_connections"], switch["current_rating"])
    end
end

## switch on/off constraints

""
function constraint_mc_switch_state_on_off(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default, relax::Bool=false)
    switch = ref(pm, nw, :switch, i)
    f_bus = switch["f_bus"]
    t_bus = switch["t_bus"]

    f_idx = (i, f_bus, t_bus)

    if switch["dispatchable"] != 0
        constraint_mc_switch_state_on_off(pm, nw, i, f_bus, t_bus, switch["f_connections"], switch["t_connections"]; relax=relax)
        constraint_mc_switch_power_on_off(pm, nw, f_idx; relax=relax)
    else
        if switch["state"] != 0
            constraint_mc_switch_state_closed(pm, nw, f_bus, t_bus, switch["f_connections"], switch["t_connections"])
        else
            constraint_mc_switch_state_open(pm, nw, f_idx)
        end
    end
end


# Balance constraints

"KCL including transformer arcs and load variables."
function constraint_mc_power_balance(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs_conns_branch, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_conns_switch, i)
    bus_arcs_trans = ref(pm, nw, :bus_arcs_conns_transformer, i)
    bus_gens = ref(pm, nw, :bus_conns_gen, i)
    bus_storage = ref(pm, nw, :bus_conns_storage, i)
    bus_loads = ref(pm, nw, :bus_conns_load, i)
    bus_shunts = ref(pm, nw, :bus_conns_shunt, i)

    if !haskey(con(pm, nw), :lam_kcl_r)
        con(pm, nw)[:lam_kcl_r] = Dict{Int,Array{JuMP.ConstraintRef}}()
    end

    if !haskey(con(pm, nw), :lam_kcl_i)
        con(pm, nw)[:lam_kcl_i] = Dict{Int,Array{JuMP.ConstraintRef}}()
    end

    constraint_mc_power_balance(pm, nw, i, bus["terminals"], bus["grounded"], bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_shunts)
end


""
function constraint_mc_power_balance_slack(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs_conns_branch, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_conns_switch, i)
    bus_arcs_trans = ref(pm, nw, :bus_arcs_conns_transformer, i)
    bus_gens = ref(pm, nw, :bus_conns_gen, i)
    bus_storage = ref(pm, nw, :bus_conns_storage, i)
    bus_loads = ref(pm, nw, :bus_conns_load, i)
    bus_shunts = ref(pm, nw, :bus_conns_shunt, i)

    if !haskey(con(pm, nw), :lam_kcl_r)
        con(pm, nw)[:lam_kcl_r] = Dict{Int,Array{JuMP.ConstraintRef}}()
    end

    if !haskey(con(pm, nw), :lam_kcl_i)
        con(pm, nw)[:lam_kcl_i] = Dict{Int,Array{JuMP.ConstraintRef}}()
    end

    constraint_mc_power_balance_slack(pm, nw, i, bus["terminals"], bus["grounded"], bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_shunts)
end


"KCL including transformer arcs"
function constraint_mc_power_balance_simple(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs_conns_branch, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_conns_switch, i)
    bus_arcs_trans = ref(pm, nw, :bus_arcs_conns_transformer, i)
    bus_gens = ref(pm, nw, :bus_conns_gen, i)
    bus_storage = ref(pm, nw, :bus_conns_storage, i)
    bus_loads = ref(pm, nw, :bus_conns_load, i)
    bus_shunts = ref(pm, nw, :bus_conns_shunt, i)

    if !haskey(con(pm, nw), :lam_kcl_r)
        con(pm, nw)[:lam_kcl_r] = Dict{Int,Array{JuMP.ConstraintRef}}()
    end

    if !haskey(con(pm, nw), :lam_kcl_i)
        con(pm, nw)[:lam_kcl_i] = Dict{Int,Array{JuMP.ConstraintRef}}()
    end

    constraint_mc_power_balance_simple(pm, nw, i, bus["terminals"], bus["grounded"], bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_shunts)
end


"KCL for load shed problem with transformers"
function constraint_mc_power_balance_shed(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs_conns_branch, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_conns_switch, i)
    bus_arcs_trans = ref(pm, nw, :bus_arcs_conns_transformer, i)
    bus_gens = ref(pm, nw, :bus_conns_gen, i)
    bus_storage = ref(pm, nw, :bus_conns_storage, i)
    bus_loads = ref(pm, nw, :bus_conns_load, i)
    bus_shunts = ref(pm, nw, :bus_conns_shunt, i)

    if !haskey(con(pm, nw), :lam_kcl_r)
        con(pm, nw)[:lam_kcl_r] = Dict{Int,Array{JuMP.ConstraintRef}}()
    end

    if !haskey(con(pm, nw), :lam_kcl_i)
        con(pm, nw)[:lam_kcl_i] = Dict{Int,Array{JuMP.ConstraintRef}}()
    end

    constraint_mc_power_balance_shed(pm, nw, i, bus["terminals"], bus["grounded"], bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_shunts)
end


""
function constraint_mc_current_balance(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs_conns_branch, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_conns_switch, i)
    bus_arcs_trans = ref(pm, nw, :bus_arcs_conns_transformer, i)
    bus_gens = ref(pm, nw, :bus_conns_gen, i)
    bus_storage = ref(pm, nw, :bus_conns_storage, i)
    bus_loads = ref(pm, nw, :bus_conns_load, i)
    bus_shunts = ref(pm, nw, :bus_conns_shunt, i)

    constraint_mc_current_balance(pm, nw, i, bus["terminals"], bus["grounded"], bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_shunts)
end


"ensures that power generation and demand are balanced"
function constraint_mc_network_power_balance(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    comp_bus_ids = ref(pm, nw, :components, i)

    comp_gen_ids = Set{Tuple{Int,Vector{Int}}}()
    for bus_id in comp_bus_ids, gen_id in ref(pm, nw, :bus_gens, bus_id)
        push!(comp_gen_ids, (gen_id, ref(pm, nw, :gen, gen_id, "connections")))
    end

    comp_loads = Set()
    for bus_id in comp_bus_ids, load_id in ref(pm, nw, :bus_loads, bus_id)
        push!(comp_loads, ref(pm, nw, :load, load_id))
    end

    comp_shunts = Set()
    for bus_id in comp_bus_ids, shunt_id in ref(pm, nw, :bus_shunts, bus_id)
        push!(comp_shunts, ref(pm, nw, :shunt, shunt_id))
    end

    comp_branches = Set()
    for (branch_id, branch) in ref(pm, nw, :branch)
        if in(branch["f_bus"], comp_bus_ids) && in(branch["t_bus"], comp_bus_ids)
            push!(comp_branches, branch)
        end
    end

    comp_pd = Dict(load["index"] => (load["load_bus"], load["connections"], load["pd"]) for load in comp_loads)
    comp_qd = Dict(load["index"] => (load["load_bus"], load["connections"], load["qd"]) for load in comp_loads)

    comp_gs = Dict(shunt["index"] => (shunt["shunt_bus"], shunt["connections"], shunt["gs"]) for shunt in comp_shunts)
    comp_bs = Dict(shunt["index"] => (shunt["shunt_bus"], shunt["connections"], shunt["bs"]) for shunt in comp_shunts)

    comp_branch_g = Dict(branch["index"] => (branch["f_bus"], branch["t_bus"], branch["f_connections"], branch["t_connections"], branch["br_r"], branch["br_x"], fill(1.0, size(branch["br_r"])[1]), branch["g_fr"], branch["g_to"]) for branch in comp_branches)
    comp_branch_b = Dict(branch["index"] => (branch["f_bus"], branch["t_bus"], branch["f_connections"], branch["t_connections"], branch["br_r"], branch["br_x"], fill(1.0, size(branch["br_r"])[1]), branch["b_fr"], branch["b_to"]) for branch in comp_branches)

    constraint_mc_network_power_balance(pm, nw, i, comp_gen_ids, comp_pd, comp_qd, comp_gs, comp_bs, comp_branch_g, comp_branch_b)
end


# Branch constraints

"ohms constraint for branches on the from-side"
function constraint_mc_ohms_yt_from(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    G, B = calc_branch_y(branch)

    constraint_mc_ohms_yt_from(pm, nw, f_bus, t_bus, f_idx, t_idx, branch["f_connections"], branch["t_connections"], G, B, branch["g_fr"], branch["b_fr"])
end


"ohms constraint for branches on the to-side"
function constraint_mc_ohms_yt_to(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    G, B = calc_branch_y(branch)

    constraint_mc_ohms_yt_to(pm, nw, f_bus, t_bus, f_idx, t_idx, branch["f_connections"], branch["t_connections"], G, B, branch["g_to"], branch["b_to"])
end


""
function constraint_mc_model_voltage_magnitude_difference(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    r = branch["br_r"]
    x = branch["br_x"]
    g_sh_fr = branch["g_fr"]
    b_sh_fr = branch["b_fr"]

    constraint_mc_model_voltage_magnitude_difference(pm, nw, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr)
end


""
function constraint_mc_model_current(pm::AbstractUBFModels; nw::Int=nw_id_default)
    for (i,branch) in ref(pm, nw, :branch)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_idx = (i, f_bus, t_bus)

        g_sh_fr = branch["g_fr"]
        b_sh_fr = branch["b_fr"]

        constraint_mc_model_current(pm, nw, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
    end
end


""
function constraint_mc_power_losses(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
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

    constraint_mc_power_losses(pm, nw, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to)
end


"This is duplicated at PowerModelsDistribution level to correctly handle the indexing of the shunts."
function constraint_mc_voltage_angle_difference(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    constraint_mc_voltage_angle_difference(pm, nw, f_idx, branch["f_connections"], branch["t_connections"], branch["angmin"], branch["angmax"])
end


"branch thermal constraints from"
function constraint_mc_thermal_limit_from(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if haskey(branch, "rate_a")
        constraint_mc_thermal_limit_from(pm, nw, f_idx, branch["f_connections"], branch["rate_a"])
    end
end


"branch thermal constraints to"
function constraint_mc_thermal_limit_to(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    if haskey(branch, "rate_a")
        constraint_mc_thermal_limit_to(pm, nw, t_idx, branch["t_connections"], branch["rate_a"])
    end
end


""
function constraint_mc_current_from(pm::AbstractIVRModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]

    constraint_mc_current_from(pm, nw, f_bus, f_idx, branch["f_connections"], g_fr, b_fr)
end


""
function constraint_mc_current_to(pm::AbstractIVRModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g_to = branch["g_to"]
    b_to = branch["b_to"]

    constraint_mc_current_to(pm, nw, t_bus, f_idx, t_idx, branch["f_connections"], branch["t_connections"], g_to, b_to)
end


""
function constraint_mc_bus_voltage_drop(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    r = branch["br_r"]
    x = branch["br_x"]

    constraint_mc_bus_voltage_drop(pm, nw, i, f_bus, t_bus, f_idx, branch["f_connections"], branch["t_connections"], r, x)
end


# Transformer constraints

"Transformer constraints, considering winding type, conductor order, polarity and tap settings."
function constraint_mc_transformer_power(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default, fix_taps::Bool=true)
    # if ref(pm, nw_id_default, :conductors)!=3
    #     Memento.error(_LOGGER, "Transformers only work with networks with three conductors.")
    # end

    transformer = ref(pm, :transformer, i)
    f_bus = transformer["f_bus"]
    t_bus = transformer["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    configuration = transformer["configuration"]
    f_connections = transformer["f_connections"]
    t_connections = transformer["t_connections"]
    tm_set = transformer["tm_set"]
    tm_fixed = fix_taps ? ones(Bool, length(tm_set)) : transformer["tm_fix"]
    tm_scale = calculate_tm_scale(transformer, ref(pm, nw, :bus, f_bus), ref(pm, nw, :bus, t_bus))

    #TODO change data model
    # there is redundancy in specifying polarity seperately on from and to side
    #TODO change this once migrated to new data model
    pol = transformer["polarity"]

    if configuration == WYE
        constraint_mc_transformer_power_yy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == DELTA
        constraint_mc_transformer_power_dy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == "zig-zag"
        Memento.error(_LOGGER, "Zig-zag not yet supported.")
    end
end


# Load constraints

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
function constraint_mc_load_power(pm::AbstractMCPowerModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
    load = ref(pm, nw, :load, id)
    bus = ref(pm, nw,:bus, load["load_bus"])

    configuration = load["configuration"]

    a, alpha, b, beta = _load_expmodel_params(load, bus)

    if configuration==WYE
        constraint_mc_load_power_wye(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    else
        constraint_mc_load_power_delta(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    end
end


# Gernerator constraints

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
function constraint_mc_generator_power(pm::AbstractMCPowerModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)
    generator = ref(pm, nw, :gen, id)
    bus = ref(pm, nw,:bus, generator["gen_bus"])

    N = length(generator["connections"])
    pmin = get(generator, "pmin", fill(-Inf, N))
    pmax = get(generator, "pmax", fill( Inf, N))
    qmin = get(generator, "qmin", fill(-Inf, N))
    qmax = get(generator, "qmax", fill( Inf, N))

    if get(generator, "configuration", WYE) == WYE
        constraint_mc_generator_power_wye(pm, nw, id, bus["index"], generator["connections"], pmin, pmax, qmin, qmax; report=report, bounded=bounded)
    else
        constraint_mc_generator_power_delta(pm, nw, id, bus["index"], generator["connections"], pmin, pmax, qmin, qmax; report=report, bounded=bounded)
    end
end


"generator active power setpoint constraint"
function constraint_mc_gen_power_setpoint_real(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default, kwargs...)
    pg_set = ref(pm, nw, :gen, i)["pg"]
    constraint_mc_gen_power_setpoint_real(pm, nw, i, pg_set)
end


""
function constraint_mc_gen_power_on_off(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    gen = ref(pm, nw, :gen, i)
    ncnds = length(gen["connections"])

    pmin = get(gen, "pmin", fill(-Inf, ncnds))
    pmax = get(gen, "pmax", fill( Inf, ncnds))
    qmin = get(gen, "qmin", fill(-Inf, ncnds))
    qmax = get(gen, "qmax", fill( Inf, ncnds))

    constraint_mc_gen_power_on_off(pm, nw, i, gen["connections"], pmin, pmax, qmin, qmax)
end


"defines limits on active power output of a generator where bounds can't be used"
function constraint_mc_gen_active_bounds(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    gen = ref(pm, nw, :gen, i)
    bus = gen["gen_bus"]
    constraint_mc_gen_active_bounds(pm, nw, i, bus, gen["connections"], gen["pmax"], gen["pmin"])
end


"defines limits on reactive power output of a generator where bounds can't be used"
function constraint_mc_gen_reactive_bounds(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    gen = ref(pm, nw, :gen, i)
    bus = gen["gen_bus"]
    constraint_mc_gen_reactive_bounds(pm, nw, i, bus, gen["connections"], gen["qmax"], gen["qmin"])
end


# Storage constraints

"storage loss constraints, delegate to PowerModels"
function constraint_mc_storage_losses(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default, kwargs...)
    storage = ref(pm, nw, :storage, i)

    constraint_storage_losses(pm, nw, i, storage["storage_bus"], storage["r"], storage["x"], storage["p_loss"], storage["q_loss"];
        conductors = storage["connections"]
    )
end


""
function constraint_mc_storage_thermal_limit(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    storage = ref(pm, nw, :storage, i)
    constraint_mc_storage_thermal_limit(pm, nw, i, storage["connections"], storage["thermal_rating"])
end


""
function constraint_mc_storage_current_limit(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    storage = ref(pm, nw, :storage, i)
    constraint_mc_storage_current_limit(pm, nw, i, storage["storage_bus"], storage["connections"], storage["current_rating"])
end


""
function constraint_mc_storage_on_off(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    storage = ref(pm, nw, :storage, i)
    charge_ub = storage["charge_rating"]
    discharge_ub = storage["discharge_rating"]

    ncnds = length(storage["connections"])
    pmin = zeros(ncnds)
    pmax = zeros(ncnds)
    qmin = zeros(ncnds)
    qmax = zeros(ncnds)

    inj_lb, inj_ub = ref_calc_storage_injection_bounds(ref(pm, nw, :storage), ref(pm, nw, :bus))
    for (idx,c) in enumerate(storage["connections"])
        pmin[idx] = inj_lb[i][idx]
        pmax[idx] = inj_ub[i][idx]
        qmin[idx] = max(inj_lb[i][idx], ref(pm, nw, :storage, i, "qmin")[idx])
        qmax[idx] = min(inj_ub[i][idx], ref(pm, nw, :storage, i, "qmax")[idx])
    end

    constraint_mc_storage_on_off(pm, nw, i, storage["connections"], pmin, pmax, qmin, qmax, charge_ub, discharge_ub)
end


""
function constraint_storage_state(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    storage = ref(pm, nw, :storage, i)

    if haskey(ref(pm, nw), :time_elapsed)
        time_elapsed = ref(pm, nw, :time_elapsed)
    else
        Memento.warn(_LOGGER, "network data should specify time_elapsed, using 1.0 as a default")
        time_elapsed = 1.0
    end

    constraint_storage_state_initial(pm, nw, i, storage["energy"], storage["charge_efficiency"], storage["discharge_efficiency"], time_elapsed)
end


""
function constraint_storage_state(pm::AbstractMCPowerModel, i::Int, nw_1::Int, nw_2::Int)
    storage = ref(pm, nw_2, :storage, i)

    if haskey(ref(pm, nw_2), :time_elapsed)
        time_elapsed = ref(pm, nw_2, :time_elapsed)
    else
        Memento.warn(_LOGGER, "network $(nw_2) should specify time_elapsed, using 1.0 as a default")
        time_elapsed = 1.0
    end

    if haskey(ref(pm, nw_1, :storage), i)
        constraint_storage_state(pm, nw_1, nw_2, i, storage["charge_efficiency"], storage["discharge_efficiency"], time_elapsed)
    else
        # if the storage device has status=0 in nw_1, then the stored energy variable will not exist. Initialize storage from data model instead.
        Memento.warn(_LOGGER, "storage component $(i) was not found in network $(nw_1) while building constraint_storage_state between networks $(nw_1) and $(nw_2). Using the energy value from the storage component in network $(nw_2) instead")
        constraint_storage_state_initial(pm, nw_2, i, storage["energy"], storage["charge_efficiency"], storage["discharge_efficiency"], time_elapsed)
    end
end


""
function constraint_storage_complementarity_nl(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    constraint_storage_complementarity_nl(pm, nw, i)
end


""
function constraint_storage_complementarity_mi(pm::AbstractMCPowerModel, i::Int; nw::Int=nw_id_default)
    storage = ref(pm, nw, :storage, i)
    charge_ub = storage["charge_rating"]
    discharge_ub = storage["discharge_rating"]

    constraint_storage_complementarity_mi(pm, nw, i, charge_ub, discharge_ub)
end
