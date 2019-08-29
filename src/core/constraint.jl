"model current constraints"
function constraint_mc_model_current(pm::_PMs.AbstractPowerModel; kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.constraint_model_current(pm; cnd=c, kwargs...)
    end
end


"reference angle constraints"
function constraint_mc_theta_ref(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    for cnd in _PMs.conductor_ids(pm)
        constraint_mc_theta_ref(pm, nw, cnd, i)
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
