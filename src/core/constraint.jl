"do nothing by default"
function constraint_mc_model_voltage(pm::_PMs.AbstractPowerModel, n::Int)
end

# Generic thermal limit constraint
""
function constraint_mc_thermal_limit_from(pm::_PMs.AbstractPowerModel, n::Int, f_idx, rate_a)
    p_fr = _PMs.var(pm, n, :p, f_idx)
    q_fr = _PMs.var(pm, n, :q, f_idx)

    for c in _PMs.conductor_ids(pm; nw=n)
        JuMP.@constraint(pm.model, p_fr[c]^2 + q_fr[c]^2 <= rate_a[c]^2)
    end
end

""
function constraint_mc_thermal_limit_to(pm::_PMs.AbstractPowerModel, n::Int, t_idx, rate_a)
    p_to = _PMs.var(pm, n, :p, t_idx)
    q_to = _PMs.var(pm, n, :q, t_idx)

    for c in _PMs.conductor_ids(pm; nw=n)
        JuMP.@constraint(pm.model, p_to[c]^2 + q_to[c]^2 <= rate_a[c]^2)
    end
end




"model current constraints"
function constraint_mc_model_current(pm::_PMs.AbstractPowerModel; kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.constraint_model_current(pm; cnd=c, kwargs...)
    end
end

"on/off bus voltage magnitude constraint"
function constraint_mc_voltage_magnitude_on_off(pm::_PMs.AbstractPowerModel, n::Int, c::Int, i::Int, vmin, vmax)
    vm = _PMs.var(pm, n, c, :vm, i)
    z_voltage = _PMs.var(pm, n, :z_voltage, i)

    JuMP.@constraint(pm.model, vm <= vmax*z_voltage)
    JuMP.@constraint(pm.model, vm >= vmin*z_voltage)
end


"on/off bus voltage magnitude squared constraint for relaxed formulations"
function constraint_mc_voltage_magnitude_sqr_on_off(pm::_PMs.AbstractPowerModel, n::Int, c::Int, i::Int, vmin, vmax)
    w = _PMs.var(pm, n, c, :w, i)
    z_voltage = _PMs.var(pm, n, :z_voltage, i)

    if isfinite(vmax)
        JuMP.@constraint(pm.model, w <= vmax^2*z_voltage)
    end

    if isfinite(vmin)
        JuMP.@constraint(pm.model, w >= vmin^2*z_voltage)
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
