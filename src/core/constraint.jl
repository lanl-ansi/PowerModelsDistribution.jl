"do nothing by default"
function constraint_mc_model_voltage(pm::AbstractUnbalancedPowerModel, nw::Int)::Nothing
    nothing
end


"""
    constraint_mc_thermal_limit_from(pm::AbstractUnbalancedPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})::Nothing

Generic thermal limit constraint for branches (from-side)
"""
function constraint_mc_thermal_limit_from(pm::AbstractUnbalancedPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})::Nothing
    p_fr = [var(pm, nw, :p, f_idx)[c] for c in f_connections]
    q_fr = [var(pm, nw, :q, f_idx)[c] for c in f_connections]

    con(pm, nw, :mu_sm_branch)[f_idx] = mu_sm_fr = [JuMP.@constraint(pm.model, p_fr[idx]^2 + q_fr[idx]^2 <= rate_a[idx]^2) for idx in findall(rate_a .< Inf)]

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, f_idx[1])[:mu_sm_fr] = mu_sm_fr
    end
    nothing
end


"""
    constraint_mc_thermal_limit_to(pm::AbstractUnbalancedPowerModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})::Nothing

Generic thermal limit constraint for branches (to-side)
"""
function constraint_mc_thermal_limit_to(pm::AbstractUnbalancedPowerModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})::Nothing
    p_to = [var(pm, nw, :p, t_idx)[c] for c in t_connections]
    q_to = [var(pm, nw, :q, t_idx)[c] for c in t_connections]

    con(pm, nw, :mu_sm_branch)[t_idx] = mu_sm_to = [JuMP.@constraint(pm.model, p_to[idx]^2 + q_to[idx]^2 <= rate_a[idx]^2) for idx in findall(rate_a .< Inf)]

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, t_idx[1])[:mu_sm_to] = mu_sm_to
    end
    nothing
end


"""
    constraint_mc_bus_voltage_magnitude_on_off(pm::AbstractUnbalancedPowerModel, nw::Int, i::Int, vmin::Vector{<:Real}, vmax::Vector{<:Real})::Nothing

Generic on/off bus voltage magnitude constraint
"""
function constraint_mc_bus_voltage_magnitude_on_off(pm::AbstractUnbalancedPowerModel, nw::Int, i::Int, vmin::Vector{<:Real}, vmax::Vector{<:Real})::Nothing
    vm = var(pm, nw, :vm, i)
    z_voltage = var(pm, nw, :z_voltage, i)

    terminals = ref(pm, nw, :bus, i)["terminals"]
    grounded = ref(pm, nw, :bus, i)["grounded"]

    for (idx, t) in [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]
        if isfinite(vmax[idx])
            JuMP.@constraint(pm.model, vm[t] <= vmax[idx]*z_voltage)
        end

        if isfinite(vmin[idx])
            JuMP.@constraint(pm.model, vm[t] >= vmin[idx]*z_voltage)
        end
    end
    nothing
end


"on/off bus voltage magnitude squared constraint for relaxed formulations"
function constraint_mc_bus_voltage_magnitude_sqr_on_off(pm::AbstractUnbalancedPowerModel, nw::Int, i::Int, vmin::Vector{<:Real}, vmax::Vector{<:Real})::Nothing
    w = var(pm, nw, :w, i)
    z_voltage = var(pm, nw, :z_voltage, i)

    terminals = ref(pm, nw, :bus, i)["terminals"]
    grounded = ref(pm, nw, :bus, i)["grounded"]

    for (idx,t) in [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]
        if isfinite(vmax[idx])
            JuMP.@constraint(pm.model, w[t] <= vmax[idx]^2*z_voltage)
        end

        if isfinite(vmin[idx])
            JuMP.@constraint(pm.model, w[t] >= vmin[idx]^2*z_voltage)
        end
    end
    nothing
end


@doc raw"""
    constraint_mc_gen_power_setpoint_real(pm::AbstractUnbalancedPowerModel, nw::Int, i::Int, pg::Vector{<:Real})::Nothing

Generic generator real power setpoint constraint

```math
P_g == P_g^{setpoint}
```
"""
function constraint_mc_gen_power_setpoint_real(pm::AbstractUnbalancedPowerModel, nw::Int, i::Int, pg::Vector{<:Real})::Nothing
    pg_var = [var(pm, nw, :pg, i)[c] for c in ref(pm, nw, :gen, i)["connections"]]
    JuMP.@constraint(pm.model, pg_var .== pg)

    nothing
end


@doc raw"""
    constraint_mc_storage_power_setpoint_real(pm::AbstractUnbalancedPowerModel, nw::Int, i::Int, ps::Real)::Nothing

Generic storage real power setpoint constraint

```math
P_s == P_s^{setpoint}
```
"""
function constraint_mc_storage_power_setpoint_real(pm::AbstractUnbalancedPowerModel, nw::Int, i::Int, ps::Real)
    ps_var = [var(pm, nw, :ps, i)[c] for c in ref(pm, nw, :storage, i)["connections"]]

    JuMP.@constraint(pm.model, sum(ps_var) == ps)
end


"on/off constraint for generators"
function constraint_mc_gen_power_on_off(pm::AbstractUnbalancedPowerModel, nw::Int, i::Int, connections::Vector{<:Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real})
    pg = var(pm, nw, :pg, i)
    qg = var(pm, nw, :qg, i)
    z = var(pm, nw, :z_gen, i)

    for (idx, c) in enumerate(connections)
        if isfinite(pmax[idx])
            JuMP.@constraint(pm.model, pg[c] .<= pmax[idx].*z)
        end

        if isfinite(pmin[idx])
            JuMP.@constraint(pm.model, pg[c] .>= pmin[idx].*z)
        end

        if isfinite(qmax[idx])
            JuMP.@constraint(pm.model, qg[c] .<= qmax[idx].*z)
        end

        if isfinite(qmin[idx])
            JuMP.@constraint(pm.model, qg[c] .>= qmin[idx].*z)
        end
    end
    nothing
end


""
function constraint_mc_storage_thermal_limit(pm::AbstractUnbalancedPowerModel, nw::Int, i::Int, connections::Vector{Int}, rating::Real)
    ps = [var(pm, nw, :ps, i)[c] for c in connections]
    qs = [var(pm, nw, :qs, i)[c] for c in connections]

    JuMP.@constraint(pm.model, sum(ps.^2 .+ qs.^2) <= rating^2)
end


""
function constraint_mc_switch_state_open(pm::AbstractUnbalancedPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int})::Nothing
    psw = var(pm, nw, :psw, f_idx)
    qsw = var(pm, nw, :qsw, f_idx)

    JuMP.@constraint(pm.model, psw .== 0.0)
    JuMP.@constraint(pm.model, qsw .== 0.0)
    nothing
end


""
function constraint_mc_switch_power_on_off(pm::AbstractUnbalancedPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}; relax::Bool=false)::Nothing
    i, f_bus, t_bus = f_idx

    psw = var(pm, nw, :psw, f_idx)
    qsw = var(pm, nw, :qsw, f_idx)

    z = var(pm, nw, :switch_state, i)

    connections = ref(pm, nw, :switch, i)["f_connections"]

    rating = get(ref(pm, nw, :switch, i), "rate_a", fill(1e20, length(connections)))

    for (idx, c) in enumerate(connections)
        if relax
            JuMP.@constraint(pm.model, psw[c] <=  rating[idx] * z)
            JuMP.@constraint(pm.model, psw[c] >= -rating[idx] * z)
            JuMP.@constraint(pm.model, qsw[c] <=  rating[idx] * z)
            JuMP.@constraint(pm.model, qsw[c] >= -rating[idx] * z)
        else
            JuMP.@constraint(pm.model, !z => {psw[c] == 0.0})
            JuMP.@constraint(pm.model, !z => {qsw[c] == 0.0})
        end
    end
    nothing
end


@doc raw"""
    constraint_switch_thermal_limit(pm::AbstractUnbalancedPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{<:Real})::Nothing

Generic thermal limit constraint for switches (from-side)

math```
p_{fr}^2 + q_{fr}^2 \leq S_{max}^2
```
"""
function constraint_switch_thermal_limit(pm::AbstractUnbalancedPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{<:Real})::Nothing
    psw_fr = [var(pm, nw, :psw, f_idx)[c] for c in f_connections]
    qsw_fr = [var(pm, nw, :qsw, f_idx)[c] for c in f_connections]

    con(pm, nw, :mu_sm_switch_fr)[f_idx] = mu_sm_fr = [JuMP.@constraint(pm.model, psw_fr[idx]^2 + qsw_fr[idx]^2 <= rating[idx]^2) for idx in findall(rating .< Inf)]

    if _IM.report_duals(pm)
        sol(pm, nw, :switch, f_idx[1])[:mu_sm_fr] = mu_sm_fr
    end

    nothing
end


""
function constraint_storage_state_initial(pm::AbstractUnbalancedPowerModel, n::Int, i::Int, energy::Real, charge_eff::Real, discharge_eff::Real, time_elapsed::Real)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    se = var(pm, n, :se, i)

    JuMP.@constraint(pm.model, se - energy == time_elapsed*(charge_eff*sc - sd/discharge_eff))
    nothing
end


""
function constraint_storage_state(pm::AbstractUnbalancedPowerModel, n_1::Int, n_2::Int, i::Int, charge_eff::Real, discharge_eff::Real, time_elapsed::Real)
    sc_2 = var(pm, n_2, :sc, i)
    sd_2 = var(pm, n_2, :sd, i)
    se_2 = var(pm, n_2, :se, i)
    se_1 = var(pm, n_1, :se, i)

    JuMP.@constraint(pm.model, se_2 - se_1 == time_elapsed*(charge_eff*sc_2 - sd_2/discharge_eff))
    nothing
end


""
function constraint_storage_complementarity_nl(pm::AbstractUnbalancedPowerModel, n::Int, i::Int)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)

    JuMP.@constraint(pm.model, sc*sd == 0.0)
    nothing
end


""
function constraint_storage_complementarity_mi(pm::AbstractUnbalancedPowerModel, n::Int, i::Int, charge_ub::Real, discharge_ub::Real)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    sc_on = var(pm, n, :sc_on, i)
    sd_on = var(pm, n, :sd_on, i)

    JuMP.@constraint(pm.model, sc_on + sd_on == 1)
    JuMP.@constraint(pm.model, sc_on*charge_ub >= sc)
    JuMP.@constraint(pm.model, sd_on*discharge_ub >= sd)
    nothing
end


"""
    constraint_mc_branch_flow(pm::AbstractUnbalancedPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int})

For superconducting branch flow (br_r and br_x all zeros)
"""
function constraint_mc_branch_flow(pm::AbstractUnbalancedPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int})
    p_fr = [var(pm, nw, :p, f_idx)[c] for c in f_connections]
    p_to = [var(pm, nw, :p, t_idx)[c] for c in t_connections]

    q_fr = [var(pm, nw, :q, f_idx)[c] for c in f_connections]
    q_to = [var(pm, nw, :q, t_idx)[c] for c in t_connections]

    if !haskey(con(pm, nw, :branch_flow), f_idx[1])
        con(pm, nw, :branch_flow)[f_idx[1]] = [
            JuMP.@constraint(pm.model, p_fr .+ p_to .== 0),
            JuMP.@constraint(pm.model, q_fr .+ q_to .== 0)
        ]
    end
end
