"do nothing by default"
function constraint_mc_model_voltage(pm::_PM.AbstractPowerModel, n::Int)
end


"Generic thermal limit constraint from-side"
function constraint_mc_thermal_limit_from(pm::_PM.AbstractPowerModel, n::Int, f_idx, rate_a)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)

    mu_sm_fr = JuMP.@constraint(pm.model, p_fr.^2 + q_fr.^2 .<= rate_a.^2)

    if _IM.report_duals(pm)
        sol(pm, n, :branch, f_idx[1])[:mu_sm_fr] = mu_sm_fr
    end
end


"Generic thermal limit constraint to-side"
function constraint_mc_thermal_limit_to(pm::_PM.AbstractPowerModel, n::Int, t_idx, rate_a)
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)

    mu_sm_to = JuMP.@constraint(pm.model, p_to.^2 + q_to.^2 .<= rate_a.^2)

    if _IM.report_duals(pm)
        sol(pm, n, :branch, t_idx[1])[:mu_sm_to] = mu_sm_to
    end
end


"on/off bus voltage magnitude constraint"
function constraint_mc_bus_voltage_magnitude_on_off(pm::_PM.AbstractPowerModel, n::Int, i::Int, vmin, vmax)
    vm = var(pm, n, :vm, i)
    z_voltage = var(pm, n, :z_voltage, i)

    for c in conductor_ids(pm, n)
        if isfinite(vmax[c])
            JuMP.@constraint(pm.model, vm[c] <= vmax[c]*z_voltage)
        end

        if isfinite(vmin[c])
            JuMP.@constraint(pm.model, vm[c] >= vmin[c]*z_voltage)
        end
    end
end


"on/off bus voltage magnitude squared constraint for relaxed formulations"
function constraint_mc_bus_voltage_magnitude_sqr_on_off(pm::_PM.AbstractPowerModel, n::Int, i::Int, vmin, vmax)
    w = var(pm, n, :w, i)
    z_voltage = var(pm, n, :z_voltage, i)

    for c in conductor_ids(pm, n)
        if isfinite(vmax[c])
            JuMP.@constraint(pm.model, w[c] <= vmax[c]^2*z_voltage)
        end

        if isfinite(vmin[c])
            JuMP.@constraint(pm.model, w[c] >= vmin[c]^2*z_voltage)
        end
    end
end


function constraint_mc_gen_power_setpoint_real(pm::_PM.AbstractPowerModel, n::Int, i, pg)
    pg_var = var(pm, n, :pg, i)
    JuMP.@constraint(pm.model, pg_var .== pg)
end


"on/off constraint for generators"
function constraint_mc_gen_power_on_off(pm::_PM.AbstractPowerModel, n::Int, i::Int, pmin, pmax, qmin, qmax)
    pg = var(pm, n, :pg, i)
    qg = var(pm, n, :qg, i)
    z = var(pm, n, :z_gen, i)

    for c in conductor_ids(pm, n)
        if isfinite(pmax[c])
            JuMP.@constraint(pm.model, pg[c] .<= pmax[c].*z)
        end

        if isfinite(pmin[c])
            JuMP.@constraint(pm.model, pg[c] .>= pmin[c].*z)
        end

        if isfinite(qmax[c])
            JuMP.@constraint(pm.model, qg[c] .<= qmax[c].*z)
        end

        if isfinite(qmin[c])
            JuMP.@constraint(pm.model, qg[c] .>= qmin[c].*z)
        end
    end
end


""
function constraint_mc_converter_thermal_limit(pm::_PM.AbstractPowerModel, n::Int, i, rating)
    ps = var(pm, n, :ps, i)
    qs = var(pm, n, :qs, i)

    JuMP.@constraint(pm.model, ps.^2 + qs.^2 .<= rating.^2)
end


"define balance between converter and storage subsystem, store as expression in variable dict"
function constraint_converter_storage_balance(pm::_PM.AbstractPowerModel, n::Int, i::Int, converter_storage)
    sc = [var(pm, n, :sc, c) for c in converter_storage]
    sd = [var(pm, n, :sd, c) for c in converter_storage]
    #store expression to inject into the converter model
    var(pm, n, :pdc)[i] = sum(sd) - sum(sc)
end
