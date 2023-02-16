# Voltage constraints

"""
    constraint_mc_model_voltage(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default)::Nothing

Template function for model voltage constraints.
"""
function constraint_mc_model_voltage(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default)::Nothing
    constraint_mc_model_voltage(pm, nw)
    nothing
end


"""
    constraint_mc_theta_ref(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for reference angle constraints.
"""
function constraint_mc_theta_ref(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    bus = ref(pm, nw, :bus, i)
    terminals = bus["terminals"]
    if haskey(bus, "va")
        va_ref = get(ref(pm, nw, :bus, i), "va", [deg2rad.([0.0, -120.0, 120.0])..., zeros(length(terminals))...][terminals])
        constraint_mc_theta_ref(pm, nw, i, va_ref)
    end
    nothing
end


"""
    constraint_mc_bus_voltage_balance(pm::AbstractUnbalancedPowerModel, bus_id::Int; nw=nw_id_default)::Nothing

Template function for bus voltage balance constraints.

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
function constraint_mc_bus_voltage_balance(pm::AbstractUnbalancedPowerModel, bus_id::Int; nw=nw_id_default)::Nothing
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
    nothing
end


"""
    constraint_mc_voltage_magnitude_only(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for voltage magnitude setpoint constraint
"""
function constraint_mc_voltage_magnitude_only(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    bus = ref(pm, nw, :bus, i)
    if haskey(bus, "vm")
        vm_ref = get(bus, "vm", ones(length(bus["terminals"])))
        constraint_mc_voltage_magnitude_only(pm, nw, i, vm_ref)
    end
    nothing
end


"""
    constraint_mc_voltage_magnitude_bounds(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for voltage magnitude bounds constraints.

This constraint captures problem agnostic constraints that define limits for
voltage magnitudes (where variable bounds cannot be used).
Notable examples include IVRUPowerModel and ACRUPowerModel.
"""
function constraint_mc_voltage_magnitude_bounds(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    bus = ref(pm, nw, :bus, i)
    vmin = get(bus, "vmin", fill(0.0, length(bus["terminals"])))
    vmax = get(bus, "vmax", fill(Inf, length(bus["terminals"])))
    constraint_mc_voltage_magnitude_bounds(pm, nw, i, vmin, vmax)
    nothing
end

## Voltage on/off constraints
"""
    constraint_mc_bus_voltage_on_off(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default)::Nothing

Template function for on/off constraint for bus voltages"
"""
function constraint_mc_bus_voltage_on_off(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default)::Nothing
    constraint_mc_bus_voltage_on_off(pm, nw)
    nothing
end


"""
    constraint_mc_bus_voltage_magnitude_on_off(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for on/off voltage magnitude constraint
"""
function constraint_mc_bus_voltage_magnitude_on_off(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    bus = ref(pm, nw, :bus, i)

    constraint_mc_bus_voltage_magnitude_on_off(pm, nw, i, bus["vmin"], bus["vmax"])
    nothing
end


"""
    constraint_mc_bus_voltage_magnitude_sqr_on_off(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for on/off voltage magnitude squared constraint for relaxed formulations
"""
function constraint_mc_bus_voltage_magnitude_sqr_on_off(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    bus = ref(pm, nw, :bus, i)

    constraint_mc_bus_voltage_magnitude_sqr_on_off(pm, nw, i, bus["vmin"], bus["vmax"])
    nothing
end


# Switch constraints

"""
    constraint_mc_switch_state(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for switch state constraints
"""
function constraint_mc_switch_state(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    switch = ref(pm, nw, :switch, i)
    f_bus = switch["f_bus"]
    t_bus = switch["t_bus"]

    f_idx = (i, f_bus, t_bus)

    if switch["state"] != 0
        constraint_mc_switch_state_closed(pm, nw, f_bus, t_bus, switch["f_connections"], switch["t_connections"])
    else
        constraint_mc_switch_state_open(pm, nw, f_idx)
    end
    nothing
end


"""
    constraint_mc_switch_thermal_limit(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for switch thermal limit constraint
"""
function constraint_mc_switch_thermal_limit(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    switch = ref(pm, nw, :switch, i)
    f_idx = (i, switch["f_bus"], switch["t_bus"])

    if !haskey(con(pm, nw), :mu_sm_switch)
        con(pm, nw)[:mu_sm_switch] = Dict{Tuple{Int,Int,Int},Vector{JuMP.ConstraintRef}}()
    end

    if haskey(switch, "thermal_rating") && any(switch["thermal_rating"] .< Inf)
        constraint_mc_switch_thermal_limit(pm, nw, f_idx, switch["f_connections"], switch["thermal_rating"])
    end
    nothing
end


"""
    constraint_mc_switch_current_limit(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for switch current limit constraints
"""
function constraint_mc_switch_current_limit(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    switch = ref(pm, nw, :switch, i)

    if !haskey(con(pm, nw), :mu_cm_switch)
        con(pm, nw)[:mu_cm_switch] = Dict{Tuple{Int,Int,Int},Vector{JuMP.ConstraintRef}}()
    end

    if haskey(switch, "current_rating") && any(switch["current_rating"] .< Inf)
        f_idx = (i, switch["f_bus"], switch["t_bus"])
        constraint_mc_switch_current_limit(pm, nw, f_idx, switch["f_connections"], switch["current_rating"])
    end
    nothing
end

## switch on/off constraints

"""
    constraint_mc_switch_state_on_off(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default, relax::Bool=false)::Nothing

Template function for switch state on/off constraints (MLD problems)
"""
function constraint_mc_switch_state_on_off(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default, relax::Bool=false)::Nothing
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
    nothing
end


# Balance constraints

"""
    constraint_mc_power_balance(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for KCL constraints.
"""
function constraint_mc_power_balance(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
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
    nothing
end


"""
    constraint_mc_power_balance_slack(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for KCL constraints which include a slack power at every bus
"""
function constraint_mc_power_balance_slack(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
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
    nothing
end


"""
    constraint_mc_power_balance_simple(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for KCL constraints for simple load shedding
"""
function constraint_mc_power_balance_simple(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
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
    nothing
end


"""
    constraint_mc_power_balance_shed(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for KCL constraints for load shed problem
"""
function constraint_mc_power_balance_shed(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs_conns_branch, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_conns_switch, i)
    bus_arcs_trans = ref(pm, nw, :bus_arcs_conns_transformer, i)
    bus_gens = ref(pm, nw, :bus_conns_gen, i)
    bus_storage = ref(pm, nw, :bus_conns_storage, i)
    bus_storage_ne = ref(pm, nw, :bus_conns_storage_ne, i)
    bus_loads = ref(pm, nw, :bus_conns_load, i)
    bus_shunts = ref(pm, nw, :bus_conns_shunt, i)

    if !haskey(con(pm, nw), :lam_kcl_r)
        con(pm, nw)[:lam_kcl_r] = Dict{Int,Array{JuMP.ConstraintRef}}()
    end

    if !haskey(con(pm, nw), :lam_kcl_i)
        con(pm, nw)[:lam_kcl_i] = Dict{Int,Array{JuMP.ConstraintRef}}()
    end

    constraint_mc_power_balance_shed(pm, nw, i, bus["terminals"], bus["grounded"], bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_storage_ne, bus_loads, bus_shunts)
    nothing
end


"""
    constraint_mc_power_balance_capc(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)

Template function for KCL constraints with capacitor control variables.
"""
function constraint_mc_power_balance_capc(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
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

    constraint_mc_power_balance_capc(pm, nw, i, bus["terminals"], bus["grounded"], bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_shunts)
    nothing
end


"""
    constraint_mc_current_balance(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for KCL constraints in current-voltage variable space
"""
function constraint_mc_current_balance(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs_conns_branch, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_conns_switch, i)
    bus_arcs_trans = ref(pm, nw, :bus_arcs_conns_transformer, i)
    bus_gens = ref(pm, nw, :bus_conns_gen, i)
    bus_storage = ref(pm, nw, :bus_conns_storage, i)
    bus_loads = ref(pm, nw, :bus_conns_load, i)
    bus_shunts = ref(pm, nw, :bus_conns_shunt, i)

    constraint_mc_current_balance(pm, nw, i, bus["terminals"], bus["grounded"], bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_shunts)
    nothing
end


"""
    constraint_mc_network_power_balance(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for constraints that ensures that power generation and demand are balanced in OBF problem
"""
function constraint_mc_network_power_balance(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
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
    nothing
end


# Branch constraints

"""
    constraint_mc_ohms_yt_from(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)

Template function for ohms constraint for branches on the from-side
"""
function constraint_mc_ohms_yt_from(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    if all(all(isapprox.(branch[k], 0.0)) for k in ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to"])
        @info "branch $(branch["source_id"]) being treated as superconducting (effective zero impedance)"
        if !haskey(con(pm, nw), :branch_flow)
            con(pm, nw)[:branch_flow] = Dict{Int,Vector{Vector{<:JuMP.ConstraintRef}}}()
        end
        constraint_mc_branch_flow(pm, nw, f_idx, t_idx, branch["f_connections"], branch["t_connections"])
    else
        if !haskey(con(pm, nw), :ohms_yt_from)
            con(pm, nw)[:ohms_yt] = Dict{Tuple{Int,Int,Int},Vector{Vector{<:JuMP.ConstraintRef}}}()
        end
        G, B = calc_branch_y(branch)
        constraint_mc_ohms_yt_from(pm, nw, f_bus, t_bus, f_idx, t_idx, branch["f_connections"], branch["t_connections"], G, B, branch["g_fr"], branch["b_fr"])
    end
end
function constraint_fixed_load(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    constraint_fixed_load(pm, nw, i)
end


"""
    constraint_mc_ohms_yt_to(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)

Template function for ohms constraint for branches on the to-side
"""
function constraint_mc_ohms_yt_to(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    if all(all(isapprox.(branch[k], 0.0)) for k in ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to"])
        @info "branch $(branch["source_id"]) being treated as superconducting (effective zero impedance)"
        if !haskey(con(pm, nw), :branch_flow)
            con(pm, nw)[:branch_flow] = Dict{Int,Vector{Vector{<:JuMP.ConstraintRef}}}()
        end
        constraint_mc_branch_flow(pm, nw, f_idx, t_idx, branch["f_connections"], branch["t_connections"])
    else
        if !haskey(con(pm, nw), :ohms_yt_to)
            con(pm, nw)[:ohms_yt] = Dict{Tuple{Int,Int,Int},Vector{Vector{<:JuMP.ConstraintRef}}}()
        end
        G, B = calc_branch_y(branch)
        constraint_mc_ohms_yt_to(pm, nw, f_bus, t_bus, f_idx, t_idx, branch["f_connections"], branch["t_connections"], G, B, branch["g_to"], branch["b_to"])
    end
end


"""
    constraint_mc_model_voltage_magnitude_difference(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for constraints for modeling voltage magnitude difference across branches
"""
function constraint_mc_model_voltage_magnitude_difference(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
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
    nothing
end


"""
    constraint_mc_model_current(pm::AbstractUBFModels; nw::Int=nw_id_default)::Nothing

Template function for constraints for model current
"""
function constraint_mc_model_current(pm::AbstractUBFModels; nw::Int=nw_id_default)::Nothing
    for (i,branch) in ref(pm, nw, :branch)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_idx = (i, f_bus, t_bus)

        g_sh_fr = branch["g_fr"]
        b_sh_fr = branch["b_fr"]

        constraint_mc_model_current(pm, nw, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
    end
    nothing
end


"""
    constraint_mc_power_losses(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for constraints for modeling power losses across branches
"""
function constraint_mc_power_losses(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
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
    nothing
end


"""
    constraint_mc_voltage_angle_difference(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for constraints of the voltage angle difference across branches
"""
function constraint_mc_voltage_angle_difference(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    constraint_mc_voltage_angle_difference(pm, nw, f_idx, branch["f_connections"], branch["t_connections"], branch["angmin"], branch["angmax"])
    nothing
end


"""
    constraint_mc_thermal_limit_from(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch thermal constraints (from-side)
"""
function constraint_mc_thermal_limit_from(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = ref(pm, nw, :branch, i)
    f_idx = (i, branch["f_bus"], branch["t_bus"])

    if !haskey(con(pm, nw), :mu_sm_branch)
        con(pm, nw)[:mu_sm_branch] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    if haskey(branch, "rate_a") && any(branch["rate_a"] .< Inf)
        constraint_mc_thermal_limit_from(pm, nw, f_idx, branch["f_connections"], branch["rate_a"])
    end
    nothing
end


"""
    constraint_mc_thermal_limit_to(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch thermal constraints (to-side)
"""
function constraint_mc_thermal_limit_to(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = ref(pm, nw, :branch, i)
    t_idx = (i, branch["t_bus"], branch["f_bus"])

    if !haskey(con(pm, nw), :mu_sm_branch)
        con(pm, nw)[:mu_sm_branch] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    if haskey(branch, "rate_a") && any(branch["rate_a"] .< Inf)
        constraint_mc_thermal_limit_to(pm, nw, t_idx, branch["t_connections"], branch["rate_a"])
    end
    nothing
end


"""
    constraint_mc_current_from(pm::AbstractUnbalancedIVRModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for current constraints on branches (from-side)
"""
function constraint_mc_current_from(pm::AbstractUnbalancedIVRModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]

    constraint_mc_current_from(pm, nw, f_bus, f_idx, branch["f_connections"], g_fr, b_fr)
    nothing
end


"""
    constraint_mc_current_to(pm::AbstractUnbalancedIVRModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for current constraints on branches (to-side)
"""
function constraint_mc_current_to(pm::AbstractUnbalancedIVRModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g_to = branch["g_to"]
    b_to = branch["b_to"]

    constraint_mc_current_to(pm, nw, t_bus, f_idx, t_idx, branch["f_connections"], branch["t_connections"], g_to, b_to)
    nothing
end


"""
    constraint_mc_bus_voltage_drop(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for bus voltage drop constraints
"""
function constraint_mc_bus_voltage_drop(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    r = branch["br_r"]
    x = branch["br_x"]

    constraint_mc_bus_voltage_drop(pm, nw, i, f_bus, t_bus, f_idx, branch["f_connections"], branch["t_connections"], r, x)
    nothing
end


# Transformer constraints

"""
    constraint_mc_transformer_power(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default, fix_taps::Bool=true)::Nothing

Template function for Transformer constraints in Power-voltage space, considering winding type, conductor order, polarity and tap settings.
"""
function constraint_mc_transformer_power(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default, fix_taps::Bool=true)::Nothing
    transformer = ref(pm, nw, :transformer, i)
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
        error("Zig-zag not yet supported.")
    end
    nothing
end


# Load constraints

@doc raw"""
    constraint_mc_load_power(pm::AbstractUnbalancedPowerModel, id::Int; nw::Int=nw_id_default, report::Bool=true)::Nothing

Template function for load constraints.

## CONSTANT POWER

Fixes the load power sd.

```math
sd = [sd_1, sd_2, sd_3]
```

What is actually fixed, depends on whether the load is connected in delta or wye.
When connected in wye, the load power equals the per-phase power sn drawn at the
bus to which the load is connected.

```math
sd_1 = v_a.conj(i_a) = sn_a
```

## CONSTANT CURRENT

Sets the active and reactive load power sd to be proportional to
the the voltage magnitude.

```math
pd = cp.|vm|
qd = cq.|vm|
sd = cp.|vm| + j.cq.|vm|
```

## CONSTANT IMPEDANCE

Sets the active and reactive power drawn by the load to be proportional to
the square of the voltage magnitude.

```math
pd = cp.|vm|^2
qd = cq.|vm|^2
sd = cp.|vm|^2 + j.cq.|vm|^2
```

## DELTA

When connected in delta, the load power gives the reference in the delta reference
frame. This means

```math
sd_1 = v_ab.conj(i_ab) = (v_a-v_b).conj(i_ab)
```

We can relate this to the per-phase power by

```math
sn_a = v_a.conj(i_a)
    = v_a.conj(i_ab-i_ca)
    = v_a.conj(conj(s_ab/v_ab) - conj(s_ca/v_ca))
    = v_a.(s_ab/(v_a-v_b) - s_ca/(v_c-v_a))
```

So for delta, sn is constrained indirectly.
"""
function constraint_mc_load_power(pm::AbstractUnbalancedPowerModel, id::Int; nw::Int=nw_id_default, report::Bool=true)::Nothing
    load = ref(pm, nw, :load, id)
    bus = ref(pm, nw,:bus, load["load_bus"])

    configuration = load["configuration"]

    a, alpha, b, beta = _load_expmodel_params(load, bus)

    if configuration==WYE
        constraint_mc_load_power_wye(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    else
        constraint_mc_load_power_delta(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    end
    nothing
end


# Gernerator constraints

@doc raw"""
    constraint_mc_generator_power(pm::AbstractUnbalancedPowerModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)::Nothing

Template function for generator power constraints

## DELTA

When connected in delta, the load power gives the reference in the delta reference
frame. This means

```math
sd_1 = v_ab.conj(i_ab) = (v_a-v_b).conj(i_ab)
```

We can relate this to the per-phase power by

```math
sn_a = v_a.conj(i_a)
    = v_a.conj(i_ab-i_ca)
    = v_a.conj(conj(s_ab/v_ab) - conj(s_ca/v_ca))
    = v_a.(s_ab/(v_a-v_b) - s_ca/(v_c-v_a))
```

So for delta, sn is constrained indirectly.
"""
function constraint_mc_generator_power(pm::AbstractUnbalancedPowerModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)::Nothing
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
    nothing
end


function constraint_mc_generator_power_ne(pm::AbstractUnbalancedPowerModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)::Nothing
    generator = ref(pm, nw, :gen_ne, id)
    bus = ref(pm, nw,:bus_ne, generator["gen_ne_bus"])

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
    nothing
end

"""
    constraint_mc_gen_power_setpoint_real(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for generator active power setpoint constraint, for power flow problems
"""
function constraint_mc_gen_power_setpoint_real(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    pg_set = ref(pm, nw, :gen, i)["pg"]
    constraint_mc_gen_power_setpoint_real(pm, nw, i, pg_set)
    nothing
end


"""
    constraint_mc_storage_power_setpoint_real(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for storage active power setpoint constraint, for power flow problems
"""
function constraint_mc_storage_power_setpoint_real(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    ps_set = ref(pm, nw, :storage, i)["ps"]
    constraint_mc_storage_power_setpoint_real(pm, nw, i, ps_set)
    nothing
end


"""
    constraint_mc_gen_power_on_off(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for generator power on/off constraints (MLD problems)
"""
function constraint_mc_gen_power_on_off(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    gen = ref(pm, nw, :gen, i)
    ncnds = length(gen["connections"])

    pmin = get(gen, "pmin", fill(-Inf, ncnds))
    pmax = get(gen, "pmax", fill( Inf, ncnds))
    qmin = get(gen, "qmin", fill(-Inf, ncnds))
    qmax = get(gen, "qmax", fill( Inf, ncnds))

    constraint_mc_gen_power_on_off(pm, nw, i, gen["connections"], pmin, pmax, qmin, qmax)
    nothing
end


"""
    constraint_mc_gen_active_bounds(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for defining limits on active power output of a generator where bounds can't be used.
"""
function constraint_mc_gen_active_bounds(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    gen = ref(pm, nw, :gen, i)
    bus = gen["gen_bus"]
    constraint_mc_gen_active_bounds(pm, nw, i, bus, gen["connections"], gen["pmax"], gen["pmin"])
    nothing
end


"""
    constraint_mc_gen_reactive_bounds(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for defines limits on reactive power output of a generator where bounds can't be used.
"""
function constraint_mc_gen_reactive_bounds(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    gen = ref(pm, nw, :gen, i)
    bus = gen["gen_bus"]
    constraint_mc_gen_reactive_bounds(pm, nw, i, bus, gen["connections"], gen["qmax"], gen["qmin"])
    nothing
end


# Storage constraints

"""
    constraint_mc_storage_losses(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for storage loss constraints
"""
function constraint_mc_storage_losses(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    storage = ref(pm, nw, :storage, i)

    constraint_mc_storage_losses(pm, nw, i, storage["storage_bus"], storage["connections"], storage["r"], storage["x"], storage["p_loss"], storage["q_loss"])
    nothing
end

"""
    constraint_mc_storage_losses_ne(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for storage loss constraints for network expansion storage
"""
function constraint_mc_storage_losses_ne(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    storage = ref(pm, nw, :storage_ne, i)

    constraint_mc_storage_losses_ne(pm, nw, i, storage["storage_ne_bus"], storage["connections"], storage["r"], storage["x"], storage["p_loss"], storage["q_loss"])
    nothing
end

""" 
            constraint_mc_storage_ne_power_on_off(
"""

"""
    constraint_mc_storage_thermal_limit(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for storage thermal limit constraints
"""
function constraint_mc_storage_thermal_limit(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)
    storage = ref(pm, nw, :storage, i)
    constraint_mc_storage_thermal_limit(pm, nw, i, storage["connections"], storage["thermal_rating"])
end



"""
    constraint_mc_storage_thermal_limit(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for storage thermal limit constraints (network expansion)
"""
function constraint_mc_storage_thermal_limit_ne(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)
    storage = ref(pm, nw, :storage_ne, i)
    constraint_mc_storage_thermal_limit_ne(pm, nw, i, storage["connections"], storage["thermal_rating"])
end

"""
    constraint_mc_storage_current_limit(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for storage current limit constraints
"""
function constraint_mc_storage_current_limit(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    storage = ref(pm, nw, :storage, i)
    constraint_mc_storage_current_limit(pm, nw, i, storage["storage_bus"], storage["connections"], storage["current_rating"])
    nothing
end


"""
    constraint_mc_storage_on_off(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Tempate function for storage on/off constraints for MLD problems
"""
function constraint_mc_storage_on_off(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
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
        qmin[idx] = max(inj_lb[i][idx], ref(pm, nw, :storage, i, "qmin"))
        qmax[idx] = min(inj_ub[i][idx], ref(pm, nw, :storage, i, "qmax"))
    end

    constraint_mc_storage_on_off(pm, nw, i, storage["connections"], maximum(pmin), minimum(pmax), maximum(qmin), minimum(qmax), charge_ub, discharge_ub)
    nothing
end


"""
    constraint_storage_state(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for storage state constraints (non multinetwork)
"""
function constraint_storage_state(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    storage = ref(pm, nw, :storage, i)

    if haskey(ref(pm, nw), :time_elapsed)
        time_elapsed = ref(pm, nw, :time_elapsed)
    else
        @warn "network data should specify time_elapsed in hours, using 1.0 as a default"
        time_elapsed = 1.0
    end

    constraint_storage_state_initial(pm, nw, i, storage["energy"], storage["charge_efficiency"], storage["discharge_efficiency"], time_elapsed)
    nothing
end

"""
    constraint_storage_state(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for storage state constraints (Network Expansion version)
"""
function constraint_storage_state_ne(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    storage = ref(pm, nw, :storage_ne, i)

    if haskey(ref(pm, nw), :time_elapsed)
        time_elapsed = ref(pm, nw, :time_elapsed)
    else
        @warn "network data should specify time_elapsed in hours, using 1.0 as a default"
        time_elapsed = 1.0
    end

    constraint_storage_state_initial_ne(pm, nw, i, storage["energy"], storage["charge_efficiency"], storage["discharge_efficiency"], time_elapsed)
    nothing
end

"""
    constraint_storage_state(pm::AbstractUnbalancedPowerModel, i::Int, nw_1::Int, nw_2::Int)::Nothing

Template function for storage state constraints for multinetwork problems
"""
function constraint_storage_state(pm::AbstractUnbalancedPowerModel, i::Int, nw_1::Int, nw_2::Int)::Nothing
    storage = ref(pm, nw_2, :storage, i)

    if haskey(ref(pm, nw_2), :time_elapsed)
        time_elapsed = ref(pm, nw_2, :time_elapsed)
    else
        @warn "network $(nw_2) should specify time_elapsed in hours, using 1.0 as a default"
        time_elapsed = 1.0
    end

    if haskey(ref(pm, nw_1, :storage), i)
        constraint_storage_state(pm, nw_1, nw_2, i, storage["charge_efficiency"], storage["discharge_efficiency"], time_elapsed)
    else
        # if the storage device has status=0 in nw_1, then the stored energy variable will not exist. Initialize storage from data model instead.
        @warn "storage component $(i) was not found in network $(nw_1) while building constraint_storage_state between networks $(nw_1) and $(nw_2). Using the energy value from the storage component in network $(nw_2) instead"
        constraint_storage_state_initial(pm, nw_2, i, storage["energy"], storage["charge_efficiency"], storage["discharge_efficiency"], time_elapsed)
    end
    nothing
end


"""
    constraint_round_trip_storage(pm::AbstractUnbalancedPowerModel, i::Int, nw_1::Int, nw_2::Int)::Nothing

Template function for storage round trip constraint
"""
function constraint_round_trip_storage(pm::AbstractUnbalancedPowerModel, i::Int, nw_1::Int, nw_2::Int)::Nothing
    storage = ref(pm, nw_2, :storage, i)

    if haskey(ref(pm, nw_2), :time_elapsed)
        time_elapsed = ref(pm, nw_2, :time_elapsed)
    else
        @warn "network $(nw_2) should specify time_elapsed in hours, using 1.0 as a default"
        time_elapsed = 1.0
    end

    if haskey(ref(pm, nw_1, :storage), i)
        constraint_storage_round_trip(pm, nw_1, nw_2, i, time_elapsed)
    else
        # if the storage device has status=0 in nw_1, then the stored energy variable will not exist. Initialize storage from data model instead.
        @warn "storage component $(i) was not found in network $(nw_1) while building constraint_storage_state between networks $(nw_1) and $(nw_2). Using the energy value from the storage component in network $(nw_2) instead"
    end
    nothing
end


"""
    constraint_round_trip_storage_ne(pm::AbstractUnbalancedPowerModel, i::Int, nw_1::Int, nw_2::Int)::Nothing

Template function for storage_ne round trip constraint
"""
function constraint_round_trip_storage_ne(pm::AbstractUnbalancedPowerModel, i::Int, nw_1::Int, nw_2::Int)::Nothing
    storage_ne = ref(pm, nw_2, :storage_ne, i)

    if haskey(ref(pm, nw_2), :time_elapsed)
        time_elapsed = ref(pm, nw_2, :time_elapsed)
    else
        @warn "network $(nw_2) should specify time_elapsed in hours, using 1.0 as a default"
        time_elapsed = 1.0
    end

    if haskey(ref(pm, nw_1, :storage_ne), i)
        constraint_storage_round_trip_ne(pm, nw_1, nw_2, i, time_elapsed)
    else
        # if the storage device has status=0 in nw_1, then the stored energy variable will not exist. Initialize storage from data model instead.
        @warn "storage component $(i) was not found in network $(nw_1) while building constraint_storage_state between networks $(nw_1) and $(nw_2). Using the energy value from the storage component in network $(nw_2) instead"
    end
    nothing
end

"""
    constraint_storage_state(pm::AbstractUnbalancedPowerModel, i::Int, nw_1::Int, nw_2::Int)::Nothing

Template function for storage state constraints for multinetwork problems (Network Expansion)
"""
function constraint_storage_state_ne(pm::AbstractUnbalancedPowerModel, i::Int, nw_1::Int, nw_2::Int)::Nothing
    storage = ref(pm, nw_2, :storage_ne, i)

    if haskey(ref(pm, nw_2), :time_elapsed)
        time_elapsed = ref(pm, nw_2, :time_elapsed)
    else
        @warn "network $(nw_2) should specify time_elapsed in hours, using 1.0 as a default"
        time_elapsed = 1.0
    end

    if haskey(ref(pm, nw_1, :storage_ne), i)
        constraint_storage_state_ne(pm, nw_1, nw_2, i, storage["charge_efficiency"], storage["discharge_efficiency"], time_elapsed)
    else
        # if the storage device has status=0 in nw_1, then the stored energy variable will not exist. Initialize storage from data model instead.
        @warn "storage component $(i) was not found in network $(nw_1) while building constraint_storage_state between networks $(nw_1) and $(nw_2). Using the energy value from the storage component in network $(nw_2) instead"
        constraint_storage_state_initial_ne(pm, nw_2, i, storage["energy"], storage["charge_efficiency"], storage["discharge_efficiency"], time_elapsed)
    end
    nothing
end

"""
    constraint_storage_complementarity_nl(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)

Template function for nonlinear storage complementarity constraints
"""
function constraint_storage_complementarity_nl(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    constraint_storage_complementarity_nl(pm, nw, i)
    nothing
end


"""
    constraint_storage_complementarity_mi(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)

Template function for mixed-integer storage complementarity constraints
"""
function constraint_storage_complementarity_mi(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    storage = ref(pm, nw, :storage, i)
    charge_ub = storage["charge_rating"]
    discharge_ub = storage["discharge_rating"]

    constraint_storage_complementarity_mi(pm, nw, i, charge_ub, discharge_ub)
    nothing
end

"""
    constraint_storage_complementarity_mi(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)

Template function for mixed-integer storage complementarity constraints
"""
function constraint_storage_complementarity_mi_ne(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    storage = ref(pm, nw, :storage_ne, i)
    charge_ub = storage["charge_rating"]
    discharge_ub = storage["discharge_rating"]

    constraint_storage_complementarity_mi_ne(pm, nw, i, charge_ub, discharge_ub)
    nothing
end

""" 
    constraint_storage_indication_expand_ne(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default
Template function for allow energy storage indicator on only if expansion variable is on
"""
function constraint_storage_indicator_expand_ne(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)
    constraint_storage_indicator_expand_ne(pm, nw, i)
    nothing
end


"""
    constraint_mc_ampacity_from(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch current limit constraint from-side
"""
function constraint_mc_ampacity_from(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = ref(pm, nw, :branch, i)
    f_idx = (i, branch["f_bus"], branch["t_bus"])

    if !haskey(con(pm, nw), :mu_cm_branch)
        con(pm, nw)[:mu_cm_branch] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    if haskey(branch, "c_rating_a") && any(branch["c_rating_a"] .< Inf)
        constraint_mc_ampacity_from(pm, nw, f_idx, branch["f_connections"], branch["c_rating_a"])
    end
    nothing
end


"""
    constraint_mc_ampacity_to(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothin

Template function for branch current limit constraint to-side
"""
function constraint_mc_ampacity_to(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = ref(pm, nw, :branch, i)
    t_idx = (i, branch["t_bus"], branch["f_bus"])

    if !haskey(con(pm, nw), :mu_cm_branch)
        con(pm, nw)[:mu_cm_branch] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    if haskey(branch, "c_rating_a") && any(branch["c_rating_a"] .< Inf)
        constraint_mc_ampacity_to(pm, nw, t_idx, branch["t_connections"], branch["c_rating_a"])
    end
    nothing
end


"""
    constraint_mc_switch_ampacity(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for switch current limit constraint from-side
"""
function constraint_mc_switch_ampacity(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    switch = ref(pm, nw, :switch, i)
    f_idx = (i, switch["f_bus"], switch["t_bus"])

    if !haskey(con(pm, nw), :mu_cm_switch)
        con(pm, nw)[:mu_cm_switch] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    if haskey(switch, "current_rating") && any(switch["current_rating"] .< Inf)
        constraint_mc_switch_ampacity(pm, nw, f_idx, switch["f_connections"], switch["current_rating"])
    end
    nothing
end
