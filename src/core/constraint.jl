"checks if a sufficient number of variables exist for the given keys collection"
function _check_var_keys(vars, keys, var_name, comp_name)
    if length(vars) < length(keys)
        error(_LOGGER, "$(var_name) decision variables appear to be missing for $(comp_name) components")
    end
end


"do nothing by default"
function constraint_mc_model_voltage(pm::AbstractMCPowerModel, nw::Int)
end


"Generic thermal limit constraint from-side"
function constraint_mc_thermal_limit_from(pm::AbstractMCPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})
    p_fr = [var(pm, nw, :p, f_idx)[c] for c in f_connections]
    q_fr = [var(pm, nw, :q, f_idx)[c] for c in f_connections]

    mu_sm_fr = JuMP.@constraint(pm.model, p_fr.^2 + q_fr.^2 .<= rate_a.^2)

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, f_idx[1])[:mu_sm_fr] = mu_sm_fr
    end
end


"Generic thermal limit constraint to-side"
function constraint_mc_thermal_limit_to(pm::AbstractMCPowerModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})
    p_to = [var(pm, nw, :p, t_idx)[c] for c in t_connections]
    q_to = [var(pm, nw, :q, t_idx)[c] for c in t_connections]

    mu_sm_to = JuMP.@constraint(pm.model, p_to.^2 + q_to.^2 .<= rate_a.^2)

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, t_idx[1])[:mu_sm_to] = mu_sm_to
    end
end


"on/off bus voltage magnitude constraint"
function constraint_mc_bus_voltage_magnitude_on_off(pm::AbstractMCPowerModel, nw::Int, i::Int, vmin::Vector{<:Real}, vmax::Vector{<:Real})
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
end


"on/off bus voltage magnitude squared constraint for relaxed formulations"
function constraint_mc_bus_voltage_magnitude_sqr_on_off(pm::AbstractMCPowerModel, nw::Int, i::Int, vmin::Vector{<:Real}, vmax::Vector{<:Real})
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
end


function constraint_mc_gen_power_setpoint_real(pm::AbstractMCPowerModel, nw::Int, i::Int, pg::Vector{<:Real})
    pg_var = [var(pm, nw, :pg, i)[c] for c in ref(pm, nw, :gen, i)["connections"]]
    JuMP.@constraint(pm.model, pg_var .== pg)
end


"on/off constraint for generators"
function constraint_mc_gen_power_on_off(pm::AbstractMCPowerModel, nw::Int, i::Int, connections::Vector{<:Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real})
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
end


""
function constraint_mc_storage_thermal_limit(pm::AbstractMCPowerModel, nw::Int, i::Int, connections::Vector{Int}, rating::Vector{<:Real})
    ps = [var(pm, nw, :ps, i)[c] for c in connections]
    qs = [var(pm, nw, :qs, i)[c] for c in connections]

    JuMP.@constraint(pm.model, ps.^2 + qs.^2 .<= rating.^2)
end


""
function constraint_mc_switch_state_open(pm::AbstractMCPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int})
    psw = var(pm, nw, :psw, f_idx)
    qsw = var(pm, nw, :qsw, f_idx)

    JuMP.@constraint(pm.model, psw .== 0.0)
    JuMP.@constraint(pm.model, qsw .== 0.0)
end


""
function constraint_mc_switch_power_on_off(pm::AbstractMCPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}; relax::Bool=false)
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
end


""
function constraint_switch_thermal_limit(pm::AbstractMCPowerModel, n::Int, f_idx::Tuple{Int,Int,Int}, connections::Vector{Int}, rating::Vector{<:Real})
    psw = var(pm, n, :psw, f_idx)
    qsw = var(pm, n, :qsw, f_idx)

    for (idx, c) in enumerate(connections)
        JuMP.@constraint(pm.model, psw[c]^2 + qsw[c]^2 <= rating[idx]^2)
    end
end


""
function constraint_storage_state_initial(pm::AbstractMCPowerModel, n::Int, i::Int, energy, charge_eff, discharge_eff, time_elapsed)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    se = var(pm, n, :se, i)

    JuMP.@constraint(pm.model, se - energy == time_elapsed*(charge_eff*sc - sd/discharge_eff))
end


""
function constraint_storage_state(pm::AbstractMCPowerModel, n_1::Int, n_2::Int, i::Int, charge_eff, discharge_eff, time_elapsed)
    sc_2 = var(pm, n_2, :sc, i)
    sd_2 = var(pm, n_2, :sd, i)
    se_2 = var(pm, n_2, :se, i)
    se_1 = var(pm, n_1, :se, i)

    JuMP.@constraint(pm.model, se_2 - se_1 == time_elapsed*(charge_eff*sc_2 - sd_2/discharge_eff))
end


""
function constraint_storage_complementarity_nl(pm::AbstractMCPowerModel, n::Int, i)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)

    JuMP.@constraint(pm.model, sc*sd == 0.0)
end


""
function constraint_storage_complementarity_mi(pm::AbstractMCPowerModel, n::Int, i, charge_ub, discharge_ub)
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    sc_on = var(pm, n, :sc_on, i)
    sd_on = var(pm, n, :sd_on, i)

    JuMP.@constraint(pm.model, sc_on + sd_on == 1)
    JuMP.@constraint(pm.model, sc_on*charge_ub >= sc)
    JuMP.@constraint(pm.model, sd_on*discharge_ub >= sd)
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
