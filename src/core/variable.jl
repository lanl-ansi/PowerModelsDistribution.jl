""
function variable_tp_voltage(pm::GenericPowerModel; kwargs...)
    for c in PMs.conductor_ids(pm)
        PMs.variable_voltage(pm, cnd=c; kwargs...)
    end
end

""
function variable_tp_branch_flow(pm::GenericPowerModel; kwargs...)
    for c in PMs.conductor_ids(pm)
        PMs.variable_branch_flow(pm, cnd=c; kwargs...)
    end
end



""
function variable_tp_voltage(pm::GenericPowerModel{T}; kwargs...) where T <: PMs.AbstractWRForm
    for c in PMs.conductor_ids(pm)
        variable_tp_voltage_magnitude_sqr(pm, cnd=c; kwargs...)
        variable_tp_voltage_product(pm, cnd=c; kwargs...)
    end
end


"variable: `w[i] >= 0` for `i` in `bus`es"
function variable_tp_voltage_magnitude_sqr(pm::GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded=true)
    bus_cnd = [(i, c) for i in ids(pm, nw, :bus) for c in PMs.conductor_ids(pm)]

    if bounded
        W = var(pm, nw)[:w] = @variable(pm.model,
            [i in bus_cnd], basename="$(nw)_w",
            lowerbound = ref(pm, nw, :bus, i[1], "vmin", i[2])^2,
            upperbound = ref(pm, nw, :bus, i[1], "vmax", i[2])^2,
            start = PMs.getval(ref(pm, nw, :bus, i[1]), "w_start", i[2], 1.001)
        )
    else
        W = var(pm, nw)[:w] = @variable(pm.model,
            [i in bus_cnd], basename="$(nw)_w",
            lowerbound = 0,
            start = PMs.getval(ref(pm, nw, :bus, i[1]), "w_start", i[2], 1.001)
        )
    end

    var(pm, nw, cnd)[:w] = Dict{Int,Any}()
    for i in ids(pm, nw, :bus)
        var(pm, nw, cnd, :w)[i] = W[(i, cnd)]
    end
end


""
function variable_tp_voltage_product(pm::GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded=true)
    bp_cndf_cndt = [(i, j, c, d) for (i,j) in keys(ref(pm, nw, :buspairs)) for c in PMs.conductor_ids(pm) for d in PMs.conductor_ids(pm)]
    bus_cnd = [(i, i, c, d) for i in ids(pm, nw, :bus) for c in PMs.conductor_ids(pm) for d in PMs.conductor_ids(pm) if c != d]
    append!(bus_cnd, bp_cndf_cndt)

    WR = var(pm, nw)[:wr] = @variable(pm.model,
        [b in bus_cnd], basename="$(nw)_wr",
        start = PMs.getval(b[1] != b[2] ? ref(pm, nw, :buspairs, b[1:2]) : ref(pm, nw, :bus, b[1]), "wr_start", b[3], 1.0)
    )

    WI = var(pm, nw)[:wi] = @variable(pm.model,
        [b in bus_cnd], basename="$(nw)_wi",
        start = PMs.getval(b[1] != b[2] ? ref(pm, nw, :buspairs, b[1:2]) : ref(pm, nw, :bus, b[1]), "wi_start", b[3])
    )

    if bounded
        # Diagonal bounds
        wr_min, wr_max, wi_min, wi_max = PMs.calc_voltage_product_bounds(ref(pm, nw, :buspairs), cnd)
        for (i, j) in ids(pm, nw, :buspairs)
            PMs.setupperbound(WR[(i, j, cnd, cnd)], wr_max[(i,j)])
            PMs.setupperbound(WI[(i, j, cnd, cnd)], wi_max[(i,j)])

            PMs.setlowerbound(WR[(i, j, cnd, cnd)], wr_min[(i,j)])
            PMs.setlowerbound(WI[(i, j, cnd, cnd)], wi_min[(i,j)])
        end

        # Off-diagonal bounds
        for c in PMs.conductor_ids(pm)
            if c != cnd
                wr_min, wr_max, wi_min, wi_max = calc_tp_voltage_product_bounds(pm, bus_cnd)
                for k in bus_cnd
                    PMs.setupperbound(WR[k], wr_max[k])
                    PMs.setupperbound(WI[k], wi_max[k])

                    PMs.setlowerbound(WR[k], wr_min[k])
                    PMs.setlowerbound(WI[k], wi_min[k])
                end
            end
        end
    end

    var(pm, nw, cnd)[:wr] = Dict{Tuple{Int,Int},Any}()
    var(pm, nw, cnd)[:wi] = Dict{Tuple{Int,Int},Any}()
    for (i, j) in ids(pm, nw, :buspairs)
        var(pm, nw, cnd, :wr)[(i,j)] = WR[(i, j, cnd, cnd)]
        var(pm, nw, cnd, :wi)[(i,j)] = WI[(i, j, cnd, cnd)]
    end
end


"variables for modeling storage units, includes grid injection and internal variables"
function variable_tp_storage(pm::GenericPowerModel; kwargs...)
    for c in PMs.conductor_ids(pm)
        PMs.variable_active_storage(pm, cnd=c; kwargs...)
        PMs.variable_reactive_storage(pm, cnd=c; kwargs...)
    end
    PMs.variable_storage_energy(pm; kwargs...)
    PMs.variable_storage_charge(pm; kwargs...)
    PMs.variable_storage_discharge(pm; kwargs...)
end


"generates variables for both `active` and `reactive` slack at each bus"
function variable_bus_power_slack(pm::GenericPowerModel; kwargs...)
    variable_active_bus_power_slack(pm; kwargs...)
    variable_reactive_bus_power_slack(pm; kwargs...)
end

""
function variable_active_bus_power_slack(pm::GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    var(pm, nw, cnd)[:p_slack] = @variable(pm.model,
        [i in ids(pm, nw, :bus)], basename="$(nw)_$(cnd)_p_slack",
        start = PMs.getval(ref(pm, nw, :bus, i), "p_slack_start", cnd)
    )
end

""
function variable_reactive_bus_power_slack(pm::GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    var(pm, nw, cnd)[:q_slack] = @variable(pm.model,
        [i in ids(pm, nw, :bus)], basename="$(nw)_$(cnd)_q_slack",
        start = PMs.getval(ref(pm, nw, :bus, i), "q_slack_start", cnd)
    )
end

"generates variables for both `active` and `reactive` power flow at each transformer"
function variable_tp_trans_flow(pm::GenericPowerModel; kwargs...)
    variable_tp_trans_active_flow(pm; kwargs...)
    variable_tp_trans_reactive_flow(pm; kwargs...)
end
function variable_tp_trans_active_flow(pm::GenericPowerModel; nw::Int=pm.cnw, bounded = true)
    # TODO add start
    for cnd in PMs.conductor_ids(pm)
        var(pm, nw, cnd)[:p_trans] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs_trans)],
            basename="$(nw)_$(cnd)_p_trans"
        )
    end
    if bounded
        #TODO add bounds
    end
end
function variable_tp_trans_reactive_flow(pm::GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded = true)
    # TODO add start
    for cnd in PMs.conductor_ids(pm)
        var(pm, nw, cnd)[:q_trans] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs_trans)],
            basename="$(nw)_$(cnd)_q_trans"
        )
    end
    if bounded
        #TODO add bounds
    end
end
function variable_tp_trans_tap(pm::GenericPowerModel; nw=pm.cnw, kwargs...)
    tr_ids = [tr_id for tr_id in ids(pm, pm.cnw, :trans)
        if ref(pm, nw, :trans, tr_id, "type")=="tap"
           && !(all(ref(pm, pm.cnw, :trans, tr_id, "tapfix")))
    ]
    if !isempty(tr_ids)
        variable_tp_trans_tap(pm::GenericPowerModel, tr_ids; kwargs...)
    end
end
function variable_tp_trans_tap(pm::GenericPowerModel, tr_ids::Array{Int,1}; nw::Int=pm.cnw, bounded=true)
    tap_tuples = [(tr_id, c) for tr_id in tr_ids for c in 1:3]
    var(pm, nw)[:tap] = @variable(pm.model,
        [(tr_id, c) in tap_tuples],
        basename="$(nw)_tap",
        start=ref(pm, nw, :trans, tr_id, "tapset")[c]
    )
    if bounded
        for (tr_id, c) in tap_tuples
            PMs.setlowerbound(var(pm, nw)[:tap][(tr_id, c)], ref(pm, nw, :trans, tr_id, "tapmin")[c])
            PMs.setupperbound(var(pm, nw)[:tap][(tr_id, c)], ref(pm, nw, :trans, tr_id, "tapmax")[c])
        end
    end
end
