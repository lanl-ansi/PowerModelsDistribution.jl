
function comp_start_value(comp::Dict{String,<:Any}, key::String, conductor::Int, default)
    if haskey(comp, key)
        return comp[key][conductor]
    else
        return default
    end
end


function comp_start_value(comp::Dict{String,<:Any}, key::String, default)
    return _PMs.comp_start_value(comp, key, default)
end


""
function variable_mc_voltage_angle(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    va = _PMs.var(pm, nw)[:va] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_va_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "va_start", 0.0)
        ) for i in _PMs.ids(pm, nw, :bus)
    )

    report && _PMs.sol_component_value(pm, nw, :bus, :va, _PMs.ids(pm, nw, :bus), va)
end

""
function variable_mc_voltage_magnitude(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    vm = _PMs.var(pm, nw)[:vm] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_vm_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "vm_start", c, 1.0)
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

""
function variable_mc_voltage_real(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    vr = _PMs.var(pm, nw)[:vr] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_vr_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "vr_start", c, 0.0)
        ) for i in _PMs.ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in _PMs.ref(pm, nw, :bus), c in cnds
            JuMP.set_lower_bound(vr[i][c], -bus["vmax"][c])
            JuMP.set_upper_bound(vr[i][c],  bus["vmax"][c])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :bus, :vr, _PMs.ids(pm, nw, :bus), vr)
end

""
function variable_mc_voltage_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    vi = _PMs.var(pm, nw)[:vi] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_vi_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "vi_start", c, 0.0)
        ) for i in _PMs.ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in _PMs.ref(pm, nw, :bus), c in cnds
            JuMP.set_lower_bound(vi[i][c], -bus["vmax"][c])
            JuMP.set_upper_bound(vi[i][c],  bus["vmax"][c])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :bus, :vi, _PMs.ids(pm, nw, :bus), vi)
end


"branch flow variables, delegated back to PowerModels"
function variable_mc_branch_flow(pm::_PMs.AbstractPowerModel; kwargs...)
    variable_mc_branch_flow_active(pm; kwargs...)
    variable_mc_branch_flow_reactive(pm; kwargs...)
end

"variable: `p[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_flow_active(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    p = _PMs.var(pm, nw)[:p] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_p_$((l,i,j))",
            start = comp_start_value(_PMs.ref(pm, nw, :branch, l), "p_start", c, 0.0)
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
    ncnds = length(cnds)

    q = _PMs.var(pm, nw)[:q] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_q_$((l,i,j))",
            start = comp_start_value(_PMs.ref(pm, nw, :branch, l), "q_start", c, 0.0)
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


"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_current_real(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    branch = _PMs.ref(pm, nw, :branch)
    bus = _PMs.ref(pm, nw, :bus)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    cr = _PMs.var(pm, nw)[:cr] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_cr_$((l,i,j))",
            start = comp_start_value(_PMs.ref(pm, nw, :branch, l), "cr_start", c, 0.0)
        ) for (l,i,j) in _PMs.ref(pm, nw, :arcs)
    )

    if bounded
        ub = Dict()
        for (l,i,j) in _PMs.ref(pm, nw, :arcs_from)
            b = branch[l]
            # ub[l] = Inf
            if haskey(b, "rate_a")
                rate_fr = b["rate_a"].*b["tap"]
                rate_to = b["rate_a"]
                ub[l]  = max.(rate_fr./bus[i]["vmin"], rate_to./bus[j]["vmin"])
            end
            if haskey(b, "c_rating_a")
                ub[l] = b["c_rating_a"]
            end
        end

        for (l,i,j) in _PMs.ref(pm, nw, :arcs)
            for c in _PMs.conductor_ids(pm; nw=nw)
                if !isinf(ub[l][c])
                    JuMP.set_lower_bound(cr[(l,i,j)][c], -ub[l][c])
                    JuMP.set_upper_bound(cr[(l,i,j)][c],  ub[l][c])
                end
            end
        end
    end

    report && _PMs.sol_component_value_edge(pm, nw, :branch, :cr_fr, :cr_to, _PMs.ref(pm, nw, :arcs_from), _PMs.ref(pm, nw, :arcs_to), cr)
end


"variable: `ci[l,i,j] ` for `(l,i,j)` in `arcs`"
function variable_mc_branch_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    branch = _PMs.ref(pm, nw, :branch)
    bus = _PMs.ref(pm, nw, :bus)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    ci = _PMs.var(pm, nw)[:ci] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_ci_$((l,i,j))",
            start = comp_start_value(_PMs.ref(pm, nw, :branch, l), "ci_start", c, 0.0)
        ) for (l,i,j) in _PMs.ref(pm, nw, :arcs)
    )

    if bounded
        ub = Dict()
        for (l,i,j) in _PMs.ref(pm, nw, :arcs_from)
            b = branch[l]
            # ub[l] = Inf
            if haskey(b, "rate_a")
                rate_fr = b["rate_a"].*b["tap"]
                rate_to = b["rate_a"]
                ub[l]  = max.(rate_fr./bus[i]["vmin"], rate_to./bus[j]["vmin"])
            end
            if haskey(b, "c_rating_a")
                ub[l] = b["c_rating_a"]
            end
        end

        for (l,i,j) in _PMs.ref(pm, nw, :arcs)
            for c in _PMs.conductor_ids(pm; nw=nw)
                if !isinf(ub[l][c])
                    JuMP.set_lower_bound(ci[(l,i,j)][c], -ub[l][c])
                    JuMP.set_upper_bound(ci[(l,i,j)][c],  ub[l][c])
                end
            end
        end
    end

    report && _PMs.sol_component_value_edge(pm, nw, :branch, :ci_fr, :ci_to, _PMs.ref(pm, nw, :arcs_from), _PMs.ref(pm, nw, :arcs_to), ci)
end


"variable: `csr[l]` for `l` in `branch`"
function variable_mc_branch_series_current_real(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    branch = _PMs.ref(pm, nw, :branch)
    bus = _PMs.ref(pm, nw, :bus)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    csr = _PMs.var(pm, nw)[:csr] = Dict(l => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_csr_$(l)",
            start = comp_start_value(_PMs.ref(pm, nw, :branch, l), "csr_start", c, 0.0)
        ) for l in _PMs.ids(pm, nw, :branch)
    )

    if bounded
        ub = Dict()
        for (l,i,j) in _PMs.ref(pm, nw, :arcs_from)
            b = branch[l]
            ub[l] = Inf
            if haskey(b, "rate_a")
                rate = b["rate_a"].*b["tap"]
                y_fr = abs.(b["g_fr"] + im*b["b_fr"])
                y_to = abs.(b["g_to"] + im*b["b_to"])
                shuntcurrent = max.(y_fr*bus[i]["vmax"].^2, y_to*bus[j]["vmax"].^2)
                seriescurrent = max.(rate./bus[i]["vmin"], rate./bus[j]["vmin"])
                ub[l] = seriescurrent + shuntcurrent
            end
            if haskey(b, "c_rating_a")
                totalcurrent = b["c_rating_a"]
                y_fr = abs.(b["g_fr"] + im*b["b_fr"])
                y_to = abs.(b["g_to"] + im*b["b_to"])
                shuntcurrent = max.(y_fr*bus[i]["vmax"].^2, y_to*bus[j]["vmax"].^2)
                ub[l] = totalcurrent + shuntcurrent
            end
        end

        for l in _PMs.ids(pm, nw, :branch)
            for c in _PMs.conductor_ids(pm; nw=nw)
                if !isinf(ub[l][c])
                    JuMP.set_lower_bound(csr[l][c], -ub[l][c])
                    JuMP.set_upper_bound(csr[l][c],  ub[l][c])
                end
            end
        end
    end
    report && _PMs.sol_component_value(pm, nw, :branch, :csr_fr, _PMs.ids(pm, nw, :branch), csr)
end


"variable: `csi[l]` for `l` in `branch`"
function variable_mc_branch_series_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    branch = _PMs.ref(pm, nw, :branch)
    bus = _PMs.ref(pm, nw, :bus)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    csi = _PMs.var(pm, nw)[:csi] = Dict(l => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_csi_$(l)",
            start = comp_start_value(_PMs.ref(pm, nw, :branch, l), "csi_start", c, 0.0)
        ) for l in _PMs.ids(pm, nw, :branch)
    )

    if bounded
        ub = Dict()
        for (l,i,j) in _PMs.ref(pm, nw, :arcs_from)
            b = branch[l]
            ub[l] = Inf
            if haskey(b, "rate_a")
                rate = b["rate_a"].*b["tap"]
                y_fr = abs.(b["g_fr"] + im*b["b_fr"])
                y_to = abs.(b["g_to"] + im*b["b_to"])
                shuntcurrent = max.(y_fr*bus[i]["vmax"].^2, y_to*bus[j]["vmax"].^2)
                seriescurrent = max.(rate./bus[i]["vmin"], rate./bus[j]["vmin"])
                ub[l] = seriescurrent + shuntcurrent
            end
            if haskey(b, "c_rating_a")
                totalcurrent = b["c_rating_a"]
                y_fr = abs.(b["g_fr"] + im*b["b_fr"])
                y_to = abs.(b["g_to"] + im*b["b_to"])
                shuntcurrent = max.(y_fr*bus[i]["vmax"].^2, y_to*bus[j]["vmax"].^2)
                ub[l] = totalcurrent + shuntcurrent
            end
        end

        for l in  _PMs.ids(pm, nw, :branch)
            for c in _PMs.conductor_ids(pm; nw=nw)
                if !isinf(ub[l][c])
                    JuMP.set_lower_bound(csi[l][c], -ub[l][c])
                    JuMP.set_upper_bound(csi[l][c],  ub[l][c])
                end
            end
        end
    end
    report && _PMs.sol_component_value(pm, nw, :branch, :csi_fr, _PMs.ids(pm, nw, :branch), csi)
end


"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_transformer_current_real(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    #trans = _PMs.ref(pm, nw, :transformer)
    #bus = _PMs.ref(pm, nw, :bus)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    cr = _PMs.var(pm, nw)[:crt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_crt_$((l,i,j))",
            start = comp_start_value(_PMs.ref(pm, nw, :transformer, l), "cr_start", c, 0.0)
        ) for (l,i,j) in _PMs.ref(pm, nw, :arcs_trans)
    )

    #TODO add bounds
    # if bounded
    #     ub = Dict()
    #     for (l,i,j) in _PMs.ref(pm, nw, :arcs_from)
    #         b = branch[l]
    #         # ub[l] = Inf
    #         if haskey(b, "rate_a")
    #             rate_fr = b["rate_a"].*b["tap"]
    #             rate_to = b["rate_a"]
    #             ub[l]  = max.(rate_fr./bus[i]["vmin"], rate_to./bus[j]["vmin"])
    #         end
    #         if haskey(b, "c_rating_a")
    #             ub[l] = b["c_rating_a"]
    #         end
    #     end
    #
    #     for (l,i,j) in _PMs.ref(pm, nw, :arcs)
    #         for c in _PMs.conductor_ids(pm; nw=nw)
    #             if !isinf(ub[l][c])
    #                 JuMP.set_lower_bound(cr[(l,i,j)][c], -ub[l][c])
    #                 JuMP.set_upper_bound(cr[(l,i,j)][c],  ub[l][c])
    #             end
    #         end
    #     end
    # end

    report && _PMs.sol_component_value_edge(pm, nw, :transformer, :cr_fr, :cr_to, _PMs.ref(pm, nw, :arcs_from_trans), _PMs.ref(pm, nw, :arcs_to_trans), cr)
end


"variable: `ci[l,i,j] ` for `(l,i,j)` in `arcs`"
function variable_mc_transformer_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    #trans = _PMs.ref(pm, nw, :transformer)
    #bus = _PMs.ref(pm, nw, :bus)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    ci = _PMs.var(pm, nw)[:cit] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_cit_$((l,i,j))",
            start = comp_start_value(_PMs.ref(pm, nw, :transformer, l), "ci_start", c, 0.0)
        ) for (l,i,j) in _PMs.ref(pm, nw, :arcs_trans)
    )

    #TODO add bounds
    # if bounded
    #     ub = Dict()
    #     for (l,i,j) in _PMs.ref(pm, nw, :arcs_from)
    #         b = branch[l]
    #         # ub[l] = Inf
    #         if haskey(b, "rate_a")
    #             rate_fr = b["rate_a"].*b["tap"]
    #             rate_to = b["rate_a"]
    #             ub[l]  = max.(rate_fr./bus[i]["vmin"], rate_to./bus[j]["vmin"])
    #         end
    #         if haskey(b, "c_rating_a")
    #             ub[l] = b["c_rating_a"]
    #         end
    #     end
    #
    #     for (l,i,j) in _PMs.ref(pm, nw, :arcs)
    #         for c in _PMs.conductor_ids(pm; nw=nw)
    #             if !isinf(ub[l][c])
    #                 JuMP.set_lower_bound(ci[(l,i,j)][c], -ub[l][c])
    #                 JuMP.set_upper_bound(ci[(l,i,j)][c],  ub[l][c])
    #             end
    #         end
    #     end
    # end

    report && _PMs.sol_component_value_edge(pm, nw, :transformer, :ci_fr, :ci_to, _PMs.ref(pm, nw, :arcs_from_trans), _PMs.ref(pm, nw, :arcs_to_trans), ci)
end



# "voltage variables, relaxed form"
# function variable_mc_voltage(pm::_PMs.AbstractWRModel; kwargs...)
#     variable_mc_voltage_magnitude_sqr(pm; kwargs...)
#     variable_mc_voltage_product(pm; kwargs...)
# end


"variable: `w[i] >= 0` for `i` in `buses"
function variable_mc_voltage_magnitude_sqr(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    w = _PMs.var(pm, nw)[:w] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_w_$(i)",
            lower_bound = 0.0,
            start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "w_start", 1.0)
        ) for i in _PMs.ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in _PMs.ref(pm, nw, :bus), c in cnds
            JuMP.set_lower_bound(w[i][c], bus["vmin"][c]^2)
            JuMP.set_upper_bound(w[i][c], bus["vmax"][c]^2)
        end
    end

    report && _PMs.sol_component_value(pm, nw, :bus, :w, _PMs.ids(pm, nw, :bus), w)
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
    ncnds = length(cnds)

    ps = _PMs.var(pm, nw)[:ps] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_ps_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :storage, i), "ps_start", c, 0.0)
        ) for i in _PMs.ids(pm, nw, :storage)
    )

    if bounded
        for c in cnds
            flow_lb, flow_ub = _PMs.ref_calc_storage_injection_bounds(_PMs.ref(pm, nw, :storage), _PMs.ref(pm, nw, :bus), c)

            for i in _PMs.ids(pm, nw, :storage)
                if !isinf(flow_lb[i])
                    JuMP.set_lower_bound(ps[i][c], flow_lb[i])
                end
                if !isinf(flow_ub[i])
                    JuMP.set_upper_bound(ps[i][c], flow_ub[i])
                end
            end
        end
    end

    report && _PMs.sol_component_value(pm, nw, :storage, :ps, _PMs.ids(pm, nw, :storage), ps)
end

function variable_mc_storage_reactive(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    qs = _PMs.var(pm, nw)[:qs] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_qs_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :storage, i), "qs_start", c, 0.0)
        ) for i in _PMs.ids(pm, nw, :storage)
    )

    if bounded
        for c in cnds
            flow_lb, flow_ub = _PMs.ref_calc_storage_injection_bounds(_PMs.ref(pm, nw, :storage), _PMs.ref(pm, nw, :bus), c)

            for i in _PMs.ids(pm, nw, :storage)
                if !isinf(flow_lb[i])
                    JuMP.set_lower_bound(qs[i][c], flow_lb[i])
                end
                if !isinf(flow_ub[i])
                    JuMP.set_upper_bound(qs[i][c], flow_ub[i])
                end
            end
        end
    end

    report && _PMs.sol_component_value(pm, nw, :storage, :qs, _PMs.ids(pm, nw, :storage), qs)
end



"generates variables for both `active` and `reactive` slack at each bus"
function variable_mc_bus_power_slack(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, kwargs...)
    variable_mc_active_bus_power_slack(pm; nw=nw, kwargs...)
    variable_mc_reactive_bus_power_slack(pm; nw=nw, kwargs...)
end


""
function variable_mc_active_bus_power_slack(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    p_slack = _PMs.var(pm, nw)[:p_slack] = Dict(i => JuMP.@variable(pm.model,
            [cnd in 1:ncnds], base_name="$(nw)_p_slack_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "p_slack_start", cnd, 0.0)
        ) for i in _PMs.ids(pm, nw, :bus)
    )

    report && _PMs.sol_component_value(pm, nw, :bus, :p_slack, _PMs.ids(pm, nw, :bus), p_slack)
end


""
function variable_mc_reactive_bus_power_slack(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    q_slack = _PMs.var(pm, nw)[:q_slack] = Dict(i => JuMP.@variable(pm.model,
            [cnd in 1:ncnds], base_name="$(nw)_q_slack_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "q_slack_start", cnd, 0.0)
        ) for i in _PMs.ids(pm, nw, :bus)
    )

    report && _PMs.sol_component_value(pm, nw, :bus, :q_slack, _PMs.ids(pm, nw, :bus), q_slack)
end


"Creates variables for both `active` and `reactive` power flow at each transformer."
function variable_mc_transformer_flow(pm::_PMs.AbstractPowerModel; kwargs...)
    variable_mc_transformer_flow_active(pm; kwargs...)
    variable_mc_transformer_flow_reactive(pm; kwargs...)
end


"Create variables for the active power flowing into all transformer windings."
function variable_mc_transformer_flow_active(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    pt = _PMs.var(pm, nw)[:pt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_pt_$((l,i,j))",
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
function variable_mc_transformer_flow_reactive(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    qt = _PMs.var(pm, nw)[:qt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_qt_$((l,i,j))",
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
function variable_mc_oltc_tap(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    # when extending to 4-wire, this should iterate only over the phase conductors
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    nph = 3
    p_oltc_ids = [id for (id,trans) in _PMs.ref(pm, nw, :transformer) if !all(trans["fixed"])]
    tap = _PMs.var(pm, nw)[:tap] = Dict(i => JuMP.@variable(pm.model,
        [p in 1:nph],
        base_name="$(nw)_tm_$(i)",
        start=_PMs.ref(pm, nw, :transformer, i, "tm")[p]
    ) for i in p_oltc_ids)
    if bounded
        for tr_id in p_oltc_ids, p in 1:nph
            JuMP.set_lower_bound(_PMs.var(pm, nw)[:tap][tr_id][p], _PMs.ref(pm, nw, :transformer, tr_id, "tm_min")[p])
            JuMP.set_upper_bound(_PMs.var(pm, nw)[:tap][tr_id][p], _PMs.ref(pm, nw, :transformer, tr_id, "tm_max")[p])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :transformer, :tap, _PMs.ids(pm, nw, :transformer), tap)
end


"""
Create a dictionary with values of type Any for the load.
Depending on the load model, this can be a parameter or a NLexpression.
These will be inserted into KCL.
"""
function variable_mc_load(pm::_PMs.AbstractPowerModel; nw=pm.cnw, bounded::Bool=true, report::Bool=true)
    pd = _PMs.var(pm, nw)[:pd] = Dict{Int, Any}()
    qd = _PMs.var(pm, nw)[:qd] = Dict{Int, Any}()
end


"Create variables for demand status"
function variable_mc_indicator_demand(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, relax::Bool=false, report::Bool=true)
    # this is not indexedon cnd; why used in start value?
    cnd = 1
    if relax
        z_demand = _PMs.var(pm, nw)[:z_demand] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :load)], base_name="$(nw)_z_demand",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(_PMs.ref(pm, nw, :load, i), "z_demand_on_start", 1.0)
        )
    else
        z_demand = _PMs.var(pm, nw)[:z_demand] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :load)], base_name="$(nw)_z_demand",
            binary = true,
            start = comp_start_value(_PMs.ref(pm, nw, :load, i), "z_demand_on_start", 1.0)
        )
    end

    # expressions for pd and qd
    pd = _PMs.var(pm, nw)[:pd] = Dict(i => _PMs.var(pm, nw)[:z_demand][i].*_PMs.ref(pm, nw, :load, i)["pd"].values
     for i in _PMs.ids(pm, nw, :load))
    qd = _PMs.var(pm, nw)[:qd] = Dict(i => _PMs.var(pm, nw)[:z_demand][i].*_PMs.ref(pm, nw, :load, i)["qd"].values
     for i in _PMs.ids(pm, nw, :load))

    report && _PMs.sol_component_value(pm, nw, :load, :status, _PMs.ids(pm, nw, :load), z_demand)
    report && _PMs.sol_component_value(pm, nw, :load, :pd, _PMs.ids(pm, nw, :load), pd)
    report && _PMs.sol_component_value(pm, nw, :load, :qd, _PMs.ids(pm, nw, :load), qd)
end


"Create variables for shunt status"
function variable_mc_indicator_shunt(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, relax=false, report::Bool=true)
    # this is not indexedon cnd; why used in start value?
    cnd = 1
    if relax
        z_shunt = _PMs.var(pm, nw)[:z_shunt] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :shunt)], base_name="$(nw)_z_shunt",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(_PMs.ref(pm, nw, :shunt, i), "z_shunt_on_start", 1.0)
        )
    else
        z_shunt = _PMs.var(pm, nw)[:z_shunt] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :shunt)], base_name="$(nw)_z_shunt",
            binary=true,
            start = comp_start_value(_PMs.ref(pm, nw, :shunt, i), "z_shunt_on_start", 1.0)
        )
    end

    report && _PMs.sol_component_value(pm, nw, :shunt, :status, _PMs.ids(pm, nw, :shunt), z_shunt)
end


"Create variables for bus status"
function variable_mc_indicator_bus_voltage(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_voltage = _PMs.var(pm, nw)[:z_voltage] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_z_voltage",
            binary = true,
            start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "z_voltage_start", 1.0)
        )
    else
        z_voltage =_PMs.var(pm, nw)[:z_voltage] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_z_voltage",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "z_voltage_start", 1.0)
        )
    end

    report && _PMs.sol_component_value(pm, nw, :bus, :status, _PMs.ids(pm, nw, :bus), z_voltage)
end


"Create variables for generator status"
function variable_mc_indicator_generation(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_gen = _PMs.var(pm, nw)[:z_gen] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :gen)], base_name="$(nw)_z_gen",
            binary = true,
            start = comp_start_value(_PMs.ref(pm, nw, :gen, i), "z_gen_start", 1.0)
        )
    else
        z_gen = _PMs.var(pm, nw)[:z_gen] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :gen)], base_name="$(nw)_z_gen",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(_PMs.ref(pm, nw, :gen, i), "z_gen_start", 1.0)
        )
    end

    report && _PMs.sol_component_value(pm, nw, :gen, :gen_status, _PMs.ids(pm, nw, :gen), z_gen)
end


"Create variables for storage status"
function variable_mc_indicator_storage(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_storage = _PMs.var(pm, nw)[:z_storage] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :storage)], base_name="$(nw)-z_storage",
            binary = true,
            start = comp_start_value(_PMs.ref(pm, nw, :storage, i), "z_storage_start", 1.0)
        )
    else
        z_storage = _PMs.var(pm, nw)[:z_storage] = JuMP.@variable(pm.model,
            [i in _PMs.ids(pm, nw, :storage)], base_name="$(nw)_z_storage",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(_PMs.ref(pm, nw, :storage, i), "z_storage_start", 1.0)
        )
    end

    report && _PMs.sol_component_value(pm, nw, :storage, :status, _PMs.ids(pm, nw, :storage), z_storage)
end


"Create variables for `active` and `reactive` storage injection"
function variable_mc_on_off_storage(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, kwargs...)
    variable_mc_on_off_storage_active(pm; nw=nw, kwargs...)
    variable_mc_on_off_storage_reactive(pm; nw=nw, kwargs...)
end


"Create variables for `active` storage injection"
function variable_mc_on_off_storage_active(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    inj_lb = Dict()
    inj_ub = Dict()
    for cnd in 1:ncnds
        inj_lb[cnd], inj_ub[cnd] = _PMs.ref_calc_storage_injection_bounds(_PMs.ref(pm, nw, :storage), _PMs.ref(pm, nw, :bus), cnd)
    end

    ps = _PMs.var(pm, nw)[:ps] = Dict(i => JuMP.@variable(pm.model,
        [cnd in 1:ncnds], base_name="$(nw)_ps_$(i)",
        lower_bound = min(0, inj_lb[cnd][i]),
        upper_bound = max(0, inj_ub[cnd][i]),
        start = comp_start_value(_PMs.ref(pm, nw, :storage, i), "ps_start", cnd, 0.0)
    ) for i in _PMs.ids(pm, nw, :storage))

    report && _PMs.sol_component_value(pm, nw, :storage, :ps, _PMs.ids(pm, nw, :storage), ps)
end


"Create variables for `reactive` storage injection"
function variable_mc_on_off_storage_reactive(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    qs = _PMs.var(pm, nw)[:qs] = Dict(i => JuMP.@variable(pm.model,
        [cnd in 1:ncnds], base_name="$(nw)_qs_$(i)",
        lower_bound = min(0, _PMs.ref(pm, nw, :storage, i, "qmin")[cnd]),
        upper_bound = max(0, _PMs.ref(pm, nw, :storage, i, "qmax")[cnd]),
        start = comp_start_value(_PMs.ref(pm, nw, :storage, i), "qs_start", cnd, 0.0)
    ) for i in _PMs.ids(pm, nw, :storage))

    report && _PMs.sol_component_value(pm, nw, :storage, :qs, _PMs.ids(pm, nw, :storage), qs)
end


"voltage variable magnitude squared (relaxed form)"
function variable_mc_voltage_magnitude_sqr_on_off(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    w = _PMs.var(pm, nw)[:w] = Dict(i => JuMP.@variable(pm.model,
        [c in 1:ncnds], base_name="$(nw)_w_$(i)",
        lower_bound = 0,
        upper_bound = _PMs.ref(pm, nw, :bus, i, "vmax")[c]^2,
        start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "w_start", c, 1.001)
    ) for i in _PMs.ids(pm, nw, :bus))

    report && _PMs.sol_component_value(pm, nw, :bus, :w, _PMs.ids(pm, nw, :bus), w)
end


"on/off voltage magnitude variable"
function variable_mc_voltage_magnitude_on_off(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    vm = _PMs.var(pm, nw)[:vm] = Dict(i => JuMP.@variable(pm.model,
        [c in 1:ncnds], base_name="$(nw)_vm_$(i)",
        lower_bound = 0,
        upper_bound = _PMs.ref(pm, nw, :bus, i, "vmax")[c],
        start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "vm_start", c, 1.0)
    ) for i in _PMs.ids(pm, nw, :bus))

    report && _PMs.sol_component_value(pm, nw, :bus, :vm, _PMs.ids(pm, nw, :bus), vm)

end


"create variables for generators, delegate to PowerModels"
function variable_mc_generation(pm::_PMs.AbstractPowerModel; kwargs...)
    variable_mc_generation_active(pm; kwargs...)
    variable_mc_generation_reactive(pm; kwargs...)
end

function variable_mc_generation_active(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    pg = _PMs.var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_pg_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :gen, i), "pg_start", c, 0.0)
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
    ncnds = length(cnds)

    qg = _PMs.var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_qg_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :gen, i), "qg_start", c, 0.0)
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


"variable: `crg[j]` for `j` in `gen`"
function variable_mc_generation_current_real(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    gen = _PMs.ref(pm, nw, :gen)
    bus = _PMs.ref(pm, nw, :bus)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    crg = _PMs.var(pm, nw)[:crg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_crg_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :gen, i), "crg_start", c, 0.0)
        ) for i in _PMs.ids(pm, nw, :gen)
    )
    if bounded
        ub = Dict()
        for (i, g) in gen
            vmin = bus[g["gen_bus"]]["vmin"]
            s = abs.(max.(abs.(g["pmax"]),abs.(g["pmin"])) + im.*max.(abs.(g["qmax"]), abs.(g["qmin"])))
            ub[i] = s./vmin
        end

        for (i, g) in gen
            for c in cnds
                JuMP.set_lower_bound(crg[i][c], -ub[i][c])
                JuMP.set_upper_bound(crg[i][c], ub[i][c])
            end
        end
    end

    report && _PMs.sol_component_value(pm, nw, :gen, :crg, _PMs.ids(pm, nw, :gen), crg)
end

"variable: `cig[j]` for `j` in `gen`"
function variable_mc_generation_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    gen = _PMs.ref(pm, nw, :gen)
    bus = _PMs.ref(pm, nw, :bus)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    cig = _PMs.var(pm, nw)[:cig] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_cig_$(i)",
            start = comp_start_value(_PMs.ref(pm, nw, :gen, i), "cig_start", c, 0.0)
        ) for i in _PMs.ids(pm, nw, :gen)
    )
    if bounded
        ub = Dict()
        for (i, g) in gen
            vmin = bus[g["gen_bus"]]["vmin"]
            s = abs.(max.(abs.(g["pmax"]),abs.(g["pmin"])) + im.*max.(abs.(g["qmax"]), abs.(g["qmin"])))
            ub[i] = s./vmin
        end

        for (i, g) in gen
            for c in cnds
                JuMP.set_lower_bound(cig[i][c], -ub[i][c])
                JuMP.set_upper_bound(cig[i][c], ub[i][c])
            end
        end
    end

    report && _PMs.sol_component_value(pm, nw, :gen, :cig, _PMs.ids(pm, nw, :gen), cig)
end


function variable_mc_generation_on_off(pm::_PMs.AbstractPowerModel; kwargs...)
    variable_mc_active_generation_on_off(pm; kwargs...)
    variable_mc_reactive_generation_on_off(pm; kwargs...)
end


function variable_mc_active_generation_on_off(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(_PMs.conductor_ids(pm, nw))

    pg = _PMs.var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
        [cnd in 1:ncnds], base_name="$(nw)_pg_$(i)",
        lower_bound = min(0, _PMs.ref(pm, nw, :gen, i, "pmin")[cnd]),
        upper_bound = max(0, _PMs.ref(pm, nw, :gen, i, "pmax")[cnd]),
        start = comp_start_value(_PMs.ref(pm, nw, :gen, i), "pg_start", cnd, 0.0)
    ) for i in _PMs.ids(pm, nw, :gen))

    report && _PMs.sol_component_value(pm, nw, :gen, :pg, _PMs.ids(pm, nw, :gen), pg)
end


function variable_mc_reactive_generation_on_off(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, report::Bool=true)
    cnds = _PMs.conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    qg = _PMs.var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
        [cnd in 1:ncnds], base_name="$(nw)_qg_$(i)",
        lower_bound = min(0, _PMs.ref(pm, nw, :gen, i, "qmin")[cnd]),
        upper_bound = max(0, _PMs.ref(pm, nw, :gen, i, "qmax")[cnd]),
        start = comp_start_value(_PMs.ref(pm, nw, :gen, i), "qg_start", cnd, 0.0)
    ) for i in _PMs.ids(pm, nw, :gen))

    report && _PMs.sol_component_value(pm, nw, :gen, :qg, _PMs.ids(pm, nw, :gen), qg)
end
