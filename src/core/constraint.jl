"do nothing by default"
function constraint_mc_model_voltage(pm::_PM.AbstractPowerModel, nw::Int)
end


"Generic thermal limit constraint from-side"
function constraint_mc_thermal_limit_from(pm::_PM.AbstractPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})
    p_fr = [var(pm, nw, :p, f_idx)[c] for c in f_connections]
    q_fr = [var(pm, nw, :q, f_idx)[c] for c in f_connections]

    mu_sm_fr = JuMP.@constraint(pm.model, p_fr.^2 + q_fr.^2 .<= rate_a.^2)

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, f_idx[1])[:mu_sm_fr] = mu_sm_fr
    end
end


"Generic thermal limit constraint to-side"
function constraint_mc_thermal_limit_to(pm::_PM.AbstractPowerModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})
    p_to = [var(pm, nw, :p, t_idx)[c] for c in t_connections]
    q_to = [var(pm, nw, :q, t_idx)[c] for c in t_connections]

    mu_sm_to = JuMP.@constraint(pm.model, p_to.^2 + q_to.^2 .<= rate_a.^2)

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, t_idx[1])[:mu_sm_to] = mu_sm_to
    end
end


"on/off bus voltage magnitude constraint"
function constraint_mc_bus_voltage_magnitude_on_off(pm::_PM.AbstractPowerModel, nw::Int, i::Int, vmin::Vector{<:Real}, vmax::Vector{<:Real})
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
function constraint_mc_bus_voltage_magnitude_sqr_on_off(pm::_PM.AbstractPowerModel, nw::Int, i::Int, vmin::Vector{<:Real}, vmax::Vector{<:Real})
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


function constraint_mc_gen_power_setpoint_real(pm::_PM.AbstractPowerModel, nw::Int, i::Int, pg::Vector{<:Real})
    pg_var = [var(pm, nw, :pg, i)[c] for c in ref(pm, nw, :gen, i)["connections"]]
    JuMP.@constraint(pm.model, pg_var .== pg)
end


"on/off constraint for generators"
function constraint_mc_gen_power_on_off(pm::_PM.AbstractPowerModel, nw::Int, i::Int, connections::Vector{<:Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real})
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
function constraint_mc_storage_thermal_limit(pm::_PM.AbstractPowerModel, nw::Int, i::Int, connections::Vector{Int}, rating::Vector{<:Real})
    ps = [var(pm, nw, :ps, i)[c] for c in connections]
    qs = [var(pm, nw, :qs, i)[c] for c in connections]

    JuMP.@constraint(pm.model, ps.^2 + qs.^2 .<= rating.^2)
end


""
function constraint_mc_switch_state_open(pm::_PM.AbstractPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int})
    psw = var(pm, nw, :psw, f_idx)
    qsw = var(pm, nw, :qsw, f_idx)

    JuMP.@constraint(pm.model, psw .== 0.0)
    JuMP.@constraint(pm.model, qsw .== 0.0)
end
