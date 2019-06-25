""
function variable_tp_voltage(pm::_PMs.GenericPowerModel; kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_voltage(pm, cnd=c; kwargs...)
    end
end


""
function variable_tp_branch_flow(pm::_PMs.GenericPowerModel; kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_branch_flow(pm, cnd=c; kwargs...)
    end
end


""
function variable_tp_voltage(pm::_PMs.GenericPowerModel{T}; kwargs...) where T <: _PMs.AbstractWRForm
    for c in _PMs.conductor_ids(pm)
        variable_tp_voltage_magnitude_sqr(pm, cnd=c; kwargs...)
        variable_tp_voltage_product(pm, cnd=c; kwargs...)
    end
end


"variable: `w[i] >= 0` for `i` in `bus`es"
function variable_tp_voltage_magnitude_sqr(pm::_PMs.GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded=true)
    bus_cnd = [(i, c) for i in _PMs.ids(pm, nw, :bus) for c in _PMs.conductor_ids(pm)]

    if bounded
        W = _PMs.var(pm, nw)[:w] = JuMP.@variable(pm.model,
            [i in bus_cnd], base_name="$(nw)_w",
            lower_bound = _PMs.ref(pm, nw, :bus, i[1], "vmin", i[2])^2,
            upper_bound = _PMs.ref(pm, nw, :bus, i[1], "vmax", i[2])^2,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i[1]), "w_start", i[2], 1.001)
        )
    else
        W = _PMs.var(pm, nw)[:w] = JuMP.@variable(pm.model,
            [i in bus_cnd], base_name="$(nw)_w",
            lower_bound = 0,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i[1]), "w_start", i[2], 1.001)
        )
    end

    _PMs.var(pm, nw, cnd)[:w] = Dict{Int,Any}()
    for i in _PMs.ids(pm, nw, :bus)
        _PMs.var(pm, nw, cnd, :w)[i] = W[(i, cnd)]
    end
end


""
function variable_tp_voltage_product(pm::_PMs.GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded=true)
    bp_cndf_cndt = [(i, j, c, d) for (i,j) in keys(_PMs.ref(pm, nw, :buspairs)) for c in _PMs.conductor_ids(pm) for d in _PMs.conductor_ids(pm)]
    bus_cnd = [(i, i, c, d) for i in _PMs.ids(pm, nw, :bus) for c in _PMs.conductor_ids(pm) for d in _PMs.conductor_ids(pm) if c != d]
    append!(bus_cnd, bp_cndf_cndt)

    WR = _PMs.var(pm, nw)[:wr] = JuMP.@variable(pm.model,
        [b in bus_cnd], base_name="$(nw)_wr",
        start = _PMs.comp_start_value(b[1] != b[2] ? _PMs.ref(pm, nw, :buspairs, b[1:2]) : _PMs.ref(pm, nw, :bus, b[1]), "wr_start", b[3], 1.0)
    )

    WI = _PMs.var(pm, nw)[:wi] = JuMP.@variable(pm.model,
        [b in bus_cnd], base_name="$(nw)_wi",
        start = _PMs.comp_start_value(b[1] != b[2] ? _PMs.ref(pm, nw, :buspairs, b[1:2]) : _PMs.ref(pm, nw, :bus, b[1]), "wi_start", b[3])
    )

    if bounded
        # Diagonal bounds
        wr_min, wr_max, wi_min, wi_max = _PMs.ref_calc_voltage_product_bounds(_PMs.ref(pm, nw, :buspairs), cnd)
        for (i, j) in _PMs.ids(pm, nw, :buspairs)
            JuMP.set_upper_bound(WR[(i, j, cnd, cnd)], wr_max[(i,j)])
            JuMP.set_upper_bound(WI[(i, j, cnd, cnd)], wi_max[(i,j)])

            JuMP.set_lower_bound(WR[(i, j, cnd, cnd)], wr_min[(i,j)])
            JuMP.set_lower_bound(WI[(i, j, cnd, cnd)], wi_min[(i,j)])
        end

        # Off-diagonal bounds
        for c in _PMs.conductor_ids(pm)
            if c != cnd
                wr_min, wr_max, wi_min, wi_max = calc_tp_voltage_product_bounds(pm, bus_cnd)
                for k in bus_cnd
                    JuMP.set_upper_bound(WR[k], wr_max[k])
                    JuMP.set_upper_bound(WI[k], wi_max[k])

                    JuMP.set_lower_bound(WR[k], wr_min[k])
                    JuMP.set_lower_bound(WI[k], wi_min[k])
                end
            end
        end
    end

    _PMs.var(pm, nw, cnd)[:wr] = Dict{Tuple{Int,Int},Any}()
    _PMs.var(pm, nw, cnd)[:wi] = Dict{Tuple{Int,Int},Any}()
    for (i, j) in _PMs.ids(pm, nw, :buspairs)
        _PMs.var(pm, nw, cnd, :wr)[(i,j)] = WR[(i, j, cnd, cnd)]
        _PMs.var(pm, nw, cnd, :wi)[(i,j)] = WI[(i, j, cnd, cnd)]
    end
end


"variables for modeling storage units, includes grid injection and internal variables"
function variable_tp_storage(pm::_PMs.GenericPowerModel; kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_active_storage(pm, cnd=c; kwargs...)
        _PMs.variable_reactive_storage(pm, cnd=c; kwargs...)
    end
    _PMs.variable_storage_energy(pm; kwargs...)
    _PMs.variable_storage_charge(pm; kwargs...)
    _PMs.variable_storage_discharge(pm; kwargs...)
end


"generates variables for both `active` and `reactive` slack at each bus"
function variable_tp_bus_power_slack(pm::_PMs.GenericPowerModel; kwargs...)
    variable_tp_active_bus_power_slack(pm; kwargs...)
    variable_tp_reactive_bus_power_slack(pm; kwargs...)
end


""
function variable_tp_active_bus_power_slack(pm::_PMs.GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    _PMs.var(pm, nw, cnd)[:p_slack] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_$(cnd)_p_slack",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i), "p_slack_start", cnd)
    )
end


""
function variable_tp_reactive_bus_power_slack(pm::_PMs.GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    _PMs.var(pm, nw, cnd)[:q_slack] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_$(cnd)_q_slack",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i), "q_slack_start", cnd)
    )
end

"Creates variables for both `active` and `reactive` power flow at each transformer."
function variable_tp_trans_flow(pm::_PMs.GenericPowerModel; kwargs...)
    variable_tp_trans_active_flow(pm; kwargs...)
    variable_tp_trans_reactive_flow(pm; kwargs...)
end


"Create variables for the active power flowing into all transformer windings."
function variable_tp_trans_active_flow(pm::_PMs.GenericPowerModel; nw::Int=pm.cnw, bounded=true)
    for cnd in _PMs.conductor_ids(pm)
        _PMs.var(pm, nw, cnd)[:pt] = JuMP.@variable(pm.model,
            [(l,i,j) in _PMs.ref(pm, nw, :arcs_trans)],
            base_name="$(nw)_$(cnd)_p_trans",
            start=0
        )
        if bounded
            for arc in _PMs.ref(pm, nw, :arcs_trans)
                tr_id = arc[1]
                flow_lb  = -_PMs.ref(pm, nw, :trans, tr_id, "rate_a")[cnd]
                flow_ub  =  _PMs.ref(pm, nw, :trans, tr_id, "rate_a")[cnd]
                JuMP.set_lower_bound(_PMs.var(pm, nw, cnd, :pt, arc), flow_lb)
                JuMP.set_upper_bound(_PMs.var(pm, nw, cnd, :pt, arc), flow_ub)
            end
        end
    end
end


"Create variables for the reactive power flowing into all transformer windings."
function variable_tp_trans_reactive_flow(pm::_PMs.GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded=true)
    for cnd in _PMs.conductor_ids(pm)
        _PMs.var(pm, nw, cnd)[:qt] = JuMP.@variable(pm.model,
            [(l,i,j) in _PMs.ref(pm, nw, :arcs_trans)],
            base_name="$(nw)_$(cnd)_q_trans",
            start=0
        )
        if bounded
            for arc in _PMs.ref(pm, nw, :arcs_trans)
                tr_id = arc[1]
                flow_lb  = -_PMs.ref(pm, nw, :trans, tr_id, "rate_a")[cnd]
                flow_ub  = _PMs.ref(pm, nw, :trans, tr_id, "rate_a")[cnd]
                JuMP.set_lower_bound(_PMs.var(pm, nw, cnd, :qt, arc), flow_lb)
                JuMP.set_upper_bound(_PMs.var(pm, nw, cnd, :qt, arc), flow_ub)
            end
        end
    end
end


"Create tap variables."
function variable_tp_oltc_tap(pm::_PMs.GenericPowerModel; nw::Int=pm.cnw, bounded=true)
    nphases = 3
    oltc_ids = _PMs.ids(pm, pm.cnw, :trans)
    for c in 1:nphases
        _PMs.var(pm, nw, c)[:tap] = JuMP.@variable(pm.model,
            [i in oltc_ids],
            base_name="$(nw)_tm",
            start=_PMs.ref(pm, nw, :trans, i, "tm")[c]
        )
        if bounded
            for tr_id in oltc_ids
                JuMP.set_lower_bound(_PMs.var(pm, nw, c)[:tap][tr_id], _PMs.ref(pm, nw, :trans, tr_id, "tm_min")[c])
                JuMP.set_upper_bound(_PMs.var(pm, nw, c)[:tap][tr_id], _PMs.ref(pm, nw, :trans, tr_id, "tm_max")[c])
            end
        end
    end
end


"""
Create a dictionary with values of type Any for the load.
Depending on the load model, this can be a parameter or a NLexpression.
These will be inserted into KCL.
"""
function variable_load(pm::_PMs.GenericPowerModel; nw=pm.cnw, cnd::Int=pm.ccnd, bounded=true)
    _PMs.var(pm, nw, cnd)[:pd] = Dict{Int, Any}()
    _PMs.var(pm, nw, cnd)[:qd] = Dict{Int, Any}()
end
