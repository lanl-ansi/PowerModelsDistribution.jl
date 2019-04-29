""
function variable_tp_voltage(pm::PMs.GenericPowerModel; kwargs...)
    for c in PMs.conductor_ids(pm)
        PMs.variable_voltage(pm, cnd=c; kwargs...)
    end
end


""
function variable_tp_branch_flow(pm::PMs.GenericPowerModel; kwargs...)
    for c in PMs.conductor_ids(pm)
        PMs.variable_branch_flow(pm, cnd=c; kwargs...)
    end
end



""
function variable_tp_voltage(pm::PMs.GenericPowerModel{T}; kwargs...) where T <: PMs.AbstractWRForm
    for c in PMs.conductor_ids(pm)
        variable_tp_voltage_magnitude_sqr(pm, cnd=c; kwargs...)
        variable_tp_voltage_product(pm, cnd=c; kwargs...)
    end
end


"variable: `w[i] >= 0` for `i` in `bus`es"
function variable_tp_voltage_magnitude_sqr(pm::PMs.GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded=true)
    bus_cnd = [(i, c) for i in PMs.ids(pm, nw, :bus) for c in PMs.conductor_ids(pm)]

    if bounded
        W = PMs.var(pm, nw)[:w] = JuMP.@variable(pm.model,
            [i in bus_cnd], basename="$(nw)_w",
            lowerbound = PMs.ref(pm, nw, :bus, i[1], "vmin", i[2])^2,
            upperbound = PMs.ref(pm, nw, :bus, i[1], "vmax", i[2])^2,
            start = PMs.getval(PMs.ref(pm, nw, :bus, i[1]), "w_start", i[2], 1.001)
        )
    else
        W = PMs.var(pm, nw)[:w] = JuMP.@variable(pm.model,
            [i in bus_cnd], basename="$(nw)_w",
            lowerbound = 0,
            start = PMs.getval(PMs.ref(pm, nw, :bus, i[1]), "w_start", i[2], 1.001)
        )
    end

    PMs.var(pm, nw, cnd)[:w] = Dict{Int,Any}()
    for i in PMs.ids(pm, nw, :bus)
        PMs.var(pm, nw, cnd, :w)[i] = W[(i, cnd)]
    end
end


""
function variable_tp_voltage_product(pm::PMs.GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded=true)
    bp_cndf_cndt = [(i, j, c, d) for (i,j) in keys(PMs.ref(pm, nw, :buspairs)) for c in PMs.conductor_ids(pm) for d in PMs.conductor_ids(pm)]
    bus_cnd = [(i, i, c, d) for i in PMs.ids(pm, nw, :bus) for c in PMs.conductor_ids(pm) for d in PMs.conductor_ids(pm) if c != d]
    append!(bus_cnd, bp_cndf_cndt)

    WR = PMs.var(pm, nw)[:wr] = JuMP.@variable(pm.model,
        [b in bus_cnd], basename="$(nw)_wr",
        start = PMs.getval(b[1] != b[2] ? PMs.ref(pm, nw, :buspairs, b[1:2]) : PMs.ref(pm, nw, :bus, b[1]), "wr_start", b[3], 1.0)
    )

    WI = PMs.var(pm, nw)[:wi] = JuMP.@variable(pm.model,
        [b in bus_cnd], basename="$(nw)_wi",
        start = PMs.getval(b[1] != b[2] ? PMs.ref(pm, nw, :buspairs, b[1:2]) : PMs.ref(pm, nw, :bus, b[1]), "wi_start", b[3])
    )

    if bounded
        # Diagonal bounds
        wr_min, wr_max, wi_min, wi_max = PMs.calc_voltage_product_bounds(PMs.ref(pm, nw, :buspairs), cnd)
        for (i, j) in PMs.ids(pm, nw, :buspairs)
            PMs.JuMP.setupperbound(WR[(i, j, cnd, cnd)], wr_max[(i,j)])
            PMs.JuMP.setupperbound(WI[(i, j, cnd, cnd)], wi_max[(i,j)])

            JuMP.setlowerbound(WR[(i, j, cnd, cnd)], wr_min[(i,j)])
            JuMP.setlowerbound(WI[(i, j, cnd, cnd)], wi_min[(i,j)])
        end

        # Off-diagonal bounds
        for c in PMs.conductor_ids(pm)
            if c != cnd
                wr_min, wr_max, wi_min, wi_max = calc_tp_voltage_product_bounds(pm, bus_cnd)
                for k in bus_cnd
                    PMs.JuMP.setupperbound(WR[k], wr_max[k])
                    PMs.JuMP.setupperbound(WI[k], wi_max[k])

                    JuMP.setlowerbound(WR[k], wr_min[k])
                    JuMP.setlowerbound(WI[k], wi_min[k])
                end
            end
        end
    end

    PMs.var(pm, nw, cnd)[:wr] = Dict{Tuple{Int,Int},Any}()
    PMs.var(pm, nw, cnd)[:wi] = Dict{Tuple{Int,Int},Any}()
    for (i, j) in PMs.ids(pm, nw, :buspairs)
        PMs.var(pm, nw, cnd, :wr)[(i,j)] = WR[(i, j, cnd, cnd)]
        PMs.var(pm, nw, cnd, :wi)[(i,j)] = WI[(i, j, cnd, cnd)]
    end
end


"variables for modeling storage units, includes grid injection and internal variables"
function variable_tp_storage(pm::PMs.GenericPowerModel; kwargs...)
    for c in PMs.conductor_ids(pm)
        PMs.variable_active_storage(pm, cnd=c; kwargs...)
        PMs.variable_reactive_storage(pm, cnd=c; kwargs...)
    end
    PMs.variable_storage_energy(pm; kwargs...)
    PMs.variable_storage_charge(pm; kwargs...)
    PMs.variable_storage_discharge(pm; kwargs...)
end


"generates variables for both `active` and `reactive` slack at each bus"
function variable_bus_power_slack(pm::PMs.GenericPowerModel; kwargs...)
    variable_active_bus_power_slack(pm; kwargs...)
    variable_reactive_bus_power_slack(pm; kwargs...)
end


""
function variable_active_bus_power_slack(pm::PMs.GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    PMs.var(pm, nw, cnd)[:p_slack] = JuMP.@variable(pm.model,
        [i in PMs.ids(pm, nw, :bus)], basename="$(nw)_$(cnd)_p_slack",
        start = PMs.getval(PMs.ref(pm, nw, :bus, i), "p_slack_start", cnd)
    )
end


""
function variable_reactive_bus_power_slack(pm::PMs.GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    PMs.var(pm, nw, cnd)[:q_slack] = JuMP.@variable(pm.model,
        [i in PMs.ids(pm, nw, :bus)], basename="$(nw)_$(cnd)_q_slack",
        start = PMs.getval(PMs.ref(pm, nw, :bus, i), "q_slack_start", cnd)
    )
end

"Creates variables for both `active` and `reactive` power flow at each transformer."
function variable_tp_trans_flow(pm::PMs.GenericPowerModel; kwargs...)
    variable_tp_trans_active_flow(pm; kwargs...)
    variable_tp_trans_reactive_flow(pm; kwargs...)
end


"Create variables for the active power flowing into all transformer windings."
function variable_tp_trans_active_flow(pm::PMs.GenericPowerModel; nw::Int=pm.cnw, bounded=true)
    for cnd in PMs.conductor_ids(pm)
        PMs.var(pm, nw, cnd)[:pt] = JuMP.@variable(pm.model,
            [(l,i,j) in PMs.ref(pm, nw, :arcs_trans)],
            basename="$(nw)_$(cnd)_p_trans",
            start=0
        )
        if bounded
            for arc in PMs.ref(pm, nw, :arcs_trans)
                tr_id = arc[1]
                flow_lb  = -PMs.ref(pm, nw, :trans, tr_id, "rate_a")[cnd]
                flow_ub  =  PMs.ref(pm, nw, :trans, tr_id, "rate_a")[cnd]
                JuMP.setlowerbound(PMs.var(pm, nw, cnd, :pt, arc), flow_lb)
                JuMP.setupperbound(PMs.var(pm, nw, cnd, :pt, arc), flow_ub)
            end
        end
    end
end


"Create variables for the reactive power flowing into all transformer windings."
function variable_tp_trans_reactive_flow(pm::PMs.GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded=true)
    for cnd in PMs.conductor_ids(pm)
        PMs.var(pm, nw, cnd)[:qt] = JuMP.@variable(pm.model,
            [(l,i,j) in PMs.ref(pm, nw, :arcs_trans)],
            basename="$(nw)_$(cnd)_q_trans",
            start=0
        )
        if bounded
            for arc in PMs.ref(pm, nw, :arcs_trans)
                tr_id = arc[1]
                flow_lb  = -PMs.ref(pm, nw, :trans, tr_id, "rate_a")[cnd]
                flow_ub  = PMs.ref(pm, nw, :trans, tr_id, "rate_a")[cnd]
                JuMP.setlowerbound(PMs.var(pm, nw, cnd, :qt, arc), flow_lb)
                JuMP.setupperbound(PMs.var(pm, nw, cnd, :qt, arc), flow_ub)
            end
        end
    end
end


"Create tap variables."
function variable_tp_oltc_tap(pm::PMs.GenericPowerModel; nw::Int=pm.cnw, bounded=true)
    nphases = 3
    oltc_ids = PMs.ids(pm, pm.cnw, :trans)
    for c in 1:nphases
        PMs.var(pm, nw, c)[:tap] = JuMP.@variable(pm.model,
            [i in oltc_ids],
            basename="$(nw)_tm",
            start=PMs.ref(pm, nw, :trans, i, "tm")[c]
        )
        if bounded
            for tr_id in oltc_ids
                JuMP.setlowerbound(PMs.var(pm, nw, c)[:tap][tr_id], PMs.ref(pm, nw, :trans, tr_id, "tm_min")[c])
                JuMP.setupperbound(PMs.var(pm, nw, c)[:tap][tr_id], PMs.ref(pm, nw, :trans, tr_id, "tm_max")[c])
            end
        end
    end
end
