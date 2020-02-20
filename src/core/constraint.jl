"do nothing by default"
function constraint_mc_model_voltage(pm::_PMs.AbstractPowerModel, n::Int)
end

# Generic thermal limit constraint
""
function constraint_mc_thermal_limit_from(pm::_PMs.AbstractPowerModel, n::Int, f_idx, s_rating)
    p_fr = _PMs.var(pm, n, :p, f_idx)
    q_fr = _PMs.var(pm, n, :q, f_idx)

    mu_sm_fr = JuMP.@constraint(pm.model, p_fr.^2 + q_fr.^2 .<= s_rating.^2)

    if _PMs.report_duals(pm)
        _PMs.sol(pm, n, :branch, f_idx[1])[:mu_sm_fr] = mu_sm_fr
    end
end

""
function constraint_mc_thermal_limit_to(pm::_PMs.AbstractPowerModel, n::Int, t_idx, s_rating)
    p_to = _PMs.var(pm, n, :p, t_idx)
    q_to = _PMs.var(pm, n, :q, t_idx)

    mu_sm_to = JuMP.@constraint(pm.model, p_to.^2 + q_to.^2 .<= s_rating.^2)

    if _PMs.report_duals(pm)
        _PMs.sol(pm, n, :branch, t_idx[1])[:mu_sm_to] = mu_sm_to
    end
end


"on/off bus voltage magnitude constraint"
function constraint_mc_voltage_magnitude_on_off(pm::_PMs.AbstractPowerModel, n::Int, i::Int, vmin, vmax)
    vm = _PMs.var(pm, n, :vm, i)
    z_voltage = _PMs.var(pm, n, :z_voltage, i)

    for c in _PMs.conductor_ids(pm, n)
        if isfinite(vmax[c])
            JuMP.@constraint(pm.model, vm[c] <= vmax[c]*z_voltage)
        end

        if isfinite(vmin[c])
            JuMP.@constraint(pm.model, vm[c] >= vmin[c]*z_voltage)
        end
    end
end


"on/off bus voltage magnitude squared constraint for relaxed formulations"
function constraint_mc_voltage_magnitude_sqr_on_off(pm::_PMs.AbstractPowerModel, n::Int, i::Int, vmin, vmax)
    w = _PMs.var(pm, n, :w, i)
    z_voltage = _PMs.var(pm, n, :z_voltage, i)

    for c in _PMs.conductor_ids(pm, n)
        if isfinite(vmax[c])
            JuMP.@constraint(pm.model, w[c] <= vmax[c]^2*z_voltage)
        end

        if isfinite(vmin[c])
            JuMP.@constraint(pm.model, w[c] >= vmin[c]^2*z_voltage)
        end
    end
end


function constraint_mc_active_gen_setpoint(pm::_PMs.AbstractPowerModel, n::Int, i, pg)
    pg_var = _PMs.var(pm, n, :pg, i)
    JuMP.@constraint(pm.model, pg_var .== pg)
end


"on/off constraint for generators"
function constraint_mc_generation_on_off(pm::_PMs.AbstractPowerModel, n::Int, i::Int, pmin, pmax, qmin, qmax)
    pg = _PMs.var(pm, n, :pg, i)
    qg = _PMs.var(pm, n, :qg, i)
    z = _PMs.var(pm, n, :z_gen, i)

    JuMP.@constraint(pm.model, pg .<= pmax.*z)
    JuMP.@constraint(pm.model, pg .>= pmin.*z)
    JuMP.@constraint(pm.model, qg .<= qmax.*z)
    JuMP.@constraint(pm.model, qg .>= qmin.*z)
end

""
function constraint_mc_storage_thermal_limit(pm::_PMs.AbstractPowerModel, n::Int, i, rating)
    ps = _PMs.var(pm, n, :ps, i)
    qs = _PMs.var(pm, n, :qs, i)

    JuMP.@constraint(pm.model, ps.^2 + qs.^2 .<= rating.^2)
end
