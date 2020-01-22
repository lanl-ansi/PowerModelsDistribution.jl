
function comp_start_value(comp::Dict{String,<:Any}, key::String, conductor=1, default=0.0)
    val = _PMs.comp_start_value(comp, key, default)
    if length(val) > 1
        return val[conductor]
    else
        return val
    end
end

"voltage variables, delegated back to PowerModels"
function variable_mc_voltage(pm::_PMs.AbstractPowerModel; kwargs...)
    variable_mc_voltage_angle(pm; kwargs...)
    variable_mc_voltage_magnitude(pm; kwargs...)
end


""
function variable_mc_voltage_angle(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)

    va = _PMs.var(pm, nw)[:va] = Dict(i => JuMP.@variable(pm.model,
            [c in cnds], base_name="$(nw)_va_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "va_start", c)
        ) for i in _PMs.ids(pm, nw, :bus)
    )

    report && _PMs.sol_component_value(pm, nw, :bus, :va, _PMs.ids(pm, nw, :bus), va)
end

""
function variable_mc_voltage_magnitude(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)

    vm = _PMs.var(pm, nw)[:vm] = Dict(i => JuMP.@variable(pm.model,
            [c in cnds], base_name="$(nw)_vm_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "vm_start", c)
        ) for i in _PMs.ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in _PMs.ref(pm, nw, :bus), c in cnds
            JuMP.set_lower_bound(vm[i][c], bus["vmin"][c])
            JuMP.set_upper_bound(vm[i][c], bus["vmax"][c])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :bus, :vm, _PMs.ids(pm, nw, :bus), vm)
end



"branch flow variables, delegated back to PowerModels"
function variable_mc_branch_flow(pm::_PMs.AbstractPowerModel; kwargs...)
    variable_mc_branch_flow_active(pm; kwargs...)
    variable_mc_branch_flow_reactive(pm; kwargs...)
end

"variable: `p[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_flow_active(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)

    p = _PMs.var(pm, nw)[:p] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in cnds], base_name="$(nw)_p_$((l,i,j))",
            start = comp_start_value(_PMs.ref(pm, nw, :branch, l), "p_start", c)
        ) for (l,i,j) in _PMs.ref(pm, nw, :arcs)
    )

    if bounded
        for c in cnds
            flow_lb, flow_ub = _PMs.ref_calc_branch_flow_bounds(_PMs.ref(pm, nw, :branch), _PMs.ref(pm, nw, :bus), c)

            for arc in _PMs.ref(pm, nw, :arcs)
                l,i,j = arc
                if !isinf(flow_lb[l])
                    JuMP.set_lower_bound(p[arc][c], flow_lb[l])
                end
                if !isinf(flow_ub[l])
                    JuMP.set_upper_bound(p[arc][c], flow_ub[l])
                end
            end
        end
    end

    for (l,branch) in _PMs.ref(pm, nw, :branch)
        if haskey(branch, "pf_start")
            f_idx = (l, branch["f_bus"], branch["t_bus"])
            JuMP.set_start_value(p[f_idx], branch["pf_start"])
        end
        if haskey(branch, "pt_start")
            t_idx = (l, branch["t_bus"], branch["f_bus"])
            JuMP.set_start_value(p[t_idx], branch["pt_start"])
        end
    end

    report && _PMs.sol_component_value_edge(pm, nw, :branch, :pf, :pt, _PMs.ref(pm, nw, :arcs_from), _PMs.ref(pm, nw, :arcs_to), p)
end

"variable: `q[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_flow_reactive(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)

    q = _PMs.var(pm, nw)[:q] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in cnds], base_name="$(nw)_q_$((l,i,j))",
            start = comp_start_value(_PMs.ref(pm, nw, :branch, l), "q_start", c)
        ) for (l,i,j) in _PMs.ref(pm, nw, :arcs)
    )

    if bounded
        for c in cnds
            flow_lb, flow_ub = _PMs.ref_calc_branch_flow_bounds(_PMs.ref(pm, nw, :branch), _PMs.ref(pm, nw, :bus), c)

            for arc in _PMs.ref(pm, nw, :arcs)
                l,i,j = arc
                if !isinf(flow_lb[l])
                    JuMP.set_lower_bound(q[arc][c], flow_lb[l])
                end
                if !isinf(flow_ub[l])
                    JuMP.set_upper_bound(q[arc][c], flow_ub[l])
                end
            end
        end
    end

    for (l,branch) in _PMs.ref(pm, nw, :branch)
        if haskey(branch, "qf_start")
            f_idx = (l, branch["f_bus"], branch["t_bus"])
            JuMP.set_start_value(q[f_idx], branch["qf_start"])
        end
        if haskey(branch, "qt_start")
            t_idx = (l, branch["t_bus"], branch["f_bus"])
            JuMP.set_start_value(q[t_idx], branch["qt_start"])
        end
    end

    report && _PMs.sol_component_value_edge(pm, nw, :branch, :qf, :qt, _PMs.ref(pm, nw, :arcs_from), _PMs.ref(pm, nw, :arcs_to), q)
end



"voltage variables, relaxed form"
function variable_mc_voltage(pm::_PMs.AbstractWRModel; nw::Int=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm; nw=nw)
        variable_mc_voltage_magnitude_sqr(pm; cnd=c, nw=nw, kwargs...)
        variable_mc_voltage_product(pm; cnd=c, nw=nw, kwargs...)
    end
end


"variable: `w[i] >= 0` for `i` in `bus`es"
function variable_mc_voltage_magnitude_sqr(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded=true)
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
function variable_mc_voltage_product(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded=true)
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
                wr_min, wr_max, wi_min, wi_max = _calc_mc_voltage_product_bounds(pm, bus_cnd)
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
function variable_mc_storage(pm::_PMs.AbstractPowerModel; kwargs...)
    variable_mc_storage_active(pm; kwargs...)
    variable_mc_storage_reactive(pm; kwargs...)

    _PMs.variable_storage_energy(pm; kwargs...)
    _PMs.variable_storage_charge(pm; kwargs...)
    _PMs.variable_storage_discharge(pm; kwargs...)
end

function variable_mc_storage_active(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)

    ps = _PMs.var(pm, nw)[:ps] = Dict(i => JuMP.@variable(pm.model,
            [c in cnds], base_name="$(nw)_ps_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :storage, i), "ps_start", c)
        ) for i in _PMs.ids(pm, nw, :storage)
    )

    if bounded
        for c in cnds
            flow_lb, flow_ub = _PMs.ref_calc_storage_injection_bounds(_PMs.ref(pm, nw, :storage), _PMs.ref(pm, nw, :bus), c)

            for i in _PMs.ids(pm, nw, :storage)
                if !isinf(flow_lb[i])
                    JuMP.set_lower_bound(ps[i][c], flow_lb[i])
                end
                if !isinf(flow_ub[l])
                    JuMP.set_upper_bound(ps[i][c], flow_ub[i])
                end
            end
        end
    end

    report && _PMs.sol_component_value(pm, nw, :storage, :ps, _PMs.ids(pm, nw, :storage), ps)
end

function variable_mc_storage_reactive(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)

    qs = _PMs.var(pm, nw)[:qs] = Dict(i => JuMP.@variable(pm.model,
            [c in cnds], base_name="$(nw)_qs_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :storage, i), "qs_start", c)
        ) for i in _PMs.ids(pm, nw, :storage)
    )

    if bounded
        for c in cnds
            flow_lb, flow_ub = _PMs.ref_calc_storage_injection_bounds(_PMs.ref(pm, nw, :storage), _PMs.ref(pm, nw, :bus), c)

            for i in _PMs.ids(pm, nw, :storage)
                if !isinf(flow_lb[i])
                    JuMP.set_lower_bound(qs[i][c], flow_lb[i])
                end
                if !isinf(flow_ub[l])
                    JuMP.set_upper_bound(qs[i][c], flow_ub[i])
                end
            end
        end
    end

    report && _PMs.sol_component_value(pm, nw, :storage, :qs, _PMs.ids(pm, nw, :storage), qs)
end



"generates variables for both `active` and `reactive` slack at each bus"
function variable_mc_bus_power_slack(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, kwargs...)
    for cnd in _PMs.conductor_ids(pm; nw=nw)
        variable_mc_active_bus_power_slack(pm; cnd=cnd, nw=nw, kwargs...)
        variable_mc_reactive_bus_power_slack(pm; cnd=cnd, nw=nw, kwargs...)
    end
end


""
function variable_mc_active_bus_power_slack(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    _PMs.var(pm, nw, cnd)[:p_slack] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_$(cnd)_p_slack",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i), "p_slack_start", cnd)
    )
end


""
function variable_mc_reactive_bus_power_slack(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    _PMs.var(pm, nw, cnd)[:q_slack] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_$(cnd)_q_slack",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i), "q_slack_start", cnd)
    )
end


"Creates variables for both `active` and `reactive` power flow at each transformer."
function variable_mc_transformer_flow(pm::_PMs.AbstractPowerModel; kwargs...)
    variable_mc_transformer_active_flow(pm; kwargs...)
    variable_mc_transformer_reactive_flow(pm; kwargs...)
end


"Create variables for the active power flowing into all transformer windings."
function variable_mc_transformer_active_flow(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)

    pt = _PMs.var(pm, nw)[:pt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in cnds], base_name="$(nw)_pt_$((l,i,j))",
        ) for (l,i,j) in _PMs.ref(pm, nw, :arcs_trans)
    )

    if bounded
        for arc in _PMs.ref(pm, nw, :arcs_trans)
            tr_id = arc[1]
            flow_lb  = -_PMs.ref(pm, nw, :transformer, tr_id, "rate_a")
            flow_ub  =  _PMs.ref(pm, nw, :transformer, tr_id, "rate_a")
            for c in cnds
                JuMP.set_lower_bound(pt[arc][c], flow_lb[c])
                JuMP.set_upper_bound(pt[arc][c], flow_ub[c])
            end
        end
    end

    for (l,transformer) in _PMs.ref(pm, nw, :transformer)
        if haskey(transformer, "pf_start")
            f_idx = (l, transformer["f_bus"], transformer["t_bus"])
            JuMP.set_start_value(pt[f_idx], branch["pf_start"])
        end
        if haskey(transformer, "pt_start")
            t_idx = (l, transformer["t_bus"], transformer["f_bus"])
            JuMP.set_start_value(pt[t_idx], transformer["pt_start"])
        end
    end

    report && _PMs.sol_component_value_edge(pm, nw, :transformer, :pf, :pt, _PMs.ref(pm, nw, :arcs_from_trans), _PMs.ref(pm, nw, :arcs_to_trans), pt)
end


"Create variables for the reactive power flowing into all transformer windings."
function variable_mc_transformer_reactive_flow(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)

    qt = _PMs.var(pm, nw)[:qt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in cnds], base_name="$(nw)_qt_$((l,i,j))",
            start = 0.0
        ) for (l,i,j) in _PMs.ref(pm, nw, :arcs_trans)
    )

    if bounded
        for arc in _PMs.ref(pm, nw, :arcs_trans)
            tr_id = arc[1]
            flow_lb  = -_PMs.ref(pm, nw, :transformer, tr_id, "rate_a")
            flow_ub  =  _PMs.ref(pm, nw, :transformer, tr_id, "rate_a")
            for c in cnds
                JuMP.set_lower_bound(qt[arc][c], flow_lb[c])
                JuMP.set_upper_bound(qt[arc][c], flow_ub[c])
            end
        end
    end

    for (l,transformer) in _PMs.ref(pm, nw, :transformer)
        if haskey(transformer, "qf_start")
            f_idx = (l, transformer["f_bus"], transformer["t_bus"])
            JuMP.set_start_value(qt[f_idx], branch["qf_start"])
        end
        if haskey(transformer, "qt_start")
            t_idx = (l, transformer["t_bus"], transformer["f_bus"])
            JuMP.set_start_value(qt[t_idx], transformer["qt_start"])
        end
    end

    report && _PMs.sol_component_value_edge(pm, nw, :transformer, :qf, :qt, _PMs.ref(pm, nw, :arcs_from_trans), _PMs.ref(pm, nw, :arcs_to_trans), qt)
end


"Create tap variables."
function variable_mc_oltc_tap(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded=true)
    # when extending to 4-wire, this should iterate only over the phase conductors
    for p in _PMs.conductor_ids(pm)
        p_oltc_ids = [id for (id,trans) in _PMs.ref(pm, nw, :transformer) if !trans["fixed"][p]]
        _PMs.var(pm, nw, p)[:tap] = JuMP.@variable(pm.model,
            [i in p_oltc_ids],
            base_name="$(nw)_$(p)_tm",
            start=_PMs.ref(pm, nw, :transformer, i, "tm")[p]
        )
        if bounded
            for tr_id in p_oltc_ids
                JuMP.set_lower_bound(_PMs.var(pm, nw, p)[:tap][tr_id], _PMs.ref(pm, nw, :transformer, tr_id, "tm_min")[p])
                JuMP.set_upper_bound(_PMs.var(pm, nw, p)[:tap][tr_id], _PMs.ref(pm, nw, :transformer, tr_id, "tm_max")[p])
            end
        end
    end
end


"""
Create a dictionary with values of type Any for the load.
Depending on the load model, this can be a parameter or a NLexpression.
These will be inserted into KCL.
"""
function variable_mc_load(pm::_PMs.AbstractPowerModel; nw=pm.cnw, bounded=true)
    for cnd in _PMs.conductor_ids(pm; nw=nw)
        _PMs.var(pm, nw, cnd)[:pd] = Dict{Int, Any}()
        _PMs.var(pm, nw, cnd)[:qd] = Dict{Int, Any}()
    end
end


"Create variables for demand status"
function variable_mc_indicator_demand(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, relax=false)
    if relax
        _PMs.var(pm, nw)[:z_demand] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :load)], base_name="$(nw)_z_demand",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :load, i), "z_demand_on_start", cnd, 1.0)
        )
    else
        _PMs.var(pm, nw)[:z_demand] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :load)], base_name="$(nw)_z_demand",
            binary = true,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :load, i), "z_demand_on_start", cnd, 1.0)
        )
    end
end


"Create variables for shunt status"
function variable_mc_indicator_shunt(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, relax=false)
    if relax
        _PMs.var(pm, nw)[:z_shunt] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :shunt)], base_name="$(nw)_z_shunt",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :shunt, i), "z_shunt_on_start", cnd, 1.0)
        )
    else
        _PMs.var(pm, nw)[:z_shunt] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :shunt)], base_name="$(nw)_z_shunt",
            binary=true,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :shunt, i), "z_shunt_on_start", cnd, 1.0)
        )
    end
end


"Create variables for bus status"
function variable_mc_indicator_bus_voltage(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, relax=false)
    if !relax
        _PMs.var(pm, nw)[:z_voltage] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_z_voltage",
            binary = true,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i), "z_voltage_start", 1, 1.0)
        )
    else
        _PMs.var(pm, nw)[:z_voltage] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_z_voltage",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i), "z_voltage_start", 1, 1.0)
        )
    end
end


"Create variables for generator status"
function variable_mc_indicator_generation(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, relax=false)
    _PMs.variable_generation_indicator(pm; nw=nw, relax=relax)
end


"Create variables for storage status"
function variable_mc_indicator_storage(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, relax=false)
    if !relax
        _PMs.var(pm, nw)[:z_storage] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :storage)], base_name="$(nw)-z_storage",
            binary = true,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :storage, i), "z_storage_start", 1, 1.0)
        )
    else
        _PMs.var(pm, nw)[:z_storage] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :storage)], base_name="$(nw)_z_storage",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :storage, i), "z_storage_start", 1, 1.0)
        )
    end
end


"Create variables for `active` and `reactive` storage injection"
function variable_mc_on_off_storage(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, kwargs...)
    for cnd in _PMs.conductor_ids(pm; nw=nw)
        variabe_mc_on_off_storage_active(pm; cnd=cnd, nw=nw, kwargs...)
        variable_mc_on_off_storage_reactive(pm; cnd=cnd, nw=nw, kwargs...)
    end
end


"Create variables for `active` storage injection"
function variabe_mc_on_off_storage_active(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    inj_lb, inj_ub = _PMs.ref_calc_storage_injection_bounds(_PMs.ref(pm, nw, :storage), _PMs.ref(pm, nw, :bus), cnd)

    _PMs.var(pm, nw, cnd)[:ps] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :storage)], base_name="$(nw)_$(cnd)_ps",
        lower_bound = inj_lb[i],
        upper_bound = inj_ub[i],
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :storage, i), "ps_start", cnd)
    )
end


"Create variables for `reactive` storage injection"
function variable_mc_on_off_storage_reactive(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    _PMs.var(pm, nw, cnd)[:qs] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :storage)], base_name="$(nw)_$(cnd)_qs",
        lower_bound = min(0, _PMs.ref(pm, nw, :storage, i, "qmin", cnd)),
        upper_bound = max(0, _PMs.ref(pm, nw, :storage, i, "qmax", cnd)),
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :storage, i), "qs_start", cnd)
    )
end


"voltage variable magnitude squared (relaxed form)"
function variable_mc_voltage_magnitude_sqr_on_off(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    _PMs.var(pm, nw, cnd)[:w] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_$(cnd)_w",
        lower_bound = 0,
        upper_bound = _PMs.ref(pm, nw, :bus, i, "vmax", cnd)^2,
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i), "w_start", cnd, 1.001)
    )
end


"on/off voltage magnitude variable"
function variable_mc_voltage_magnitude_on_off(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    for cnd in _PMs.conductor_ids(pm; nw=nw)
        _PMs.var(pm, nw, cnd)[:vm] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_$(cnd)_vm",
            lower_bound = 0,
            upper_bound = _PMs.ref(pm, nw, :bus, i, "vmax", cnd),
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :bus, i), "vm_start", cnd, 1.0)
        )
    end
end


"create variables for generators, delegate to PowerModels"
function variable_mc_generation(pm::_PMs.AbstractPowerModel; kwargs...)
    variable_mc_generation_active(pm; kwargs...)
    variable_mc_generation_reactive(pm; kwargs...)
end

function variable_mc_generation_active(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)

    pg = _PMs.var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
            [c in cnds], base_name="$(nw)_pg_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :gen, i), "pg_start", c)
        ) for i in _PMs.ids(pm, nw, :gen)
    )

    if bounded
        for (i,gen) in _PMs.ref(pm, nw, :gen), c in cnds
            JuMP.set_lower_bound(pg[i][c], gen["pmin"][c])
            JuMP.set_upper_bound(pg[i][c], gen["pmax"][c])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :gen, :pg, _PMs.ids(pm, nw, :gen), pg)
end

function variable_mc_generation_reactive(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)

    qg = _PMs.var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
            [c in cnds], base_name="$(nw)_qg_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :gen, i), "qg_start", c)
        ) for i in _PMs.ids(pm, nw, :gen)
    )

    if bounded
        for (i,gen) in _PMs.ref(pm, nw, :gen), c in cnds
            JuMP.set_lower_bound(qg[i][c], gen["qmin"][c])
            JuMP.set_upper_bound(qg[i][c], gen["qmax"][c])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :gen, :qg, _PMs.ids(pm, nw, :gen), qg)
end




"create on/off variables for generators, delegate to PowerModels"
function variable_mc_generation_on_off(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm; nw=nw)
        _PMs.variable_generation_on_off(pm; cnd=c, nw=nw, kwargs...)
    end
end
