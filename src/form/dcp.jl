### simple active power only approximations (e.g. DC Power Flow)


######## AbstractDCPForm Models (has va but assumes vm is 1.0) ########

"nothing to do, these models do not have complex voltage constraints"
function constraint_tp_voltage(pm::PMs.GenericPowerModel{T}, n::Int, c::Int) where T <: PMs.AbstractDCPForm
end




######## Lossless Models ########

"Create variables for the active power flowing into all transformer windings."
function variable_tp_trans_active_flow(pm::PMs.GenericPowerModel{T}; nw::Int=pm.cnw, bounded=true) where T <: PMs.DCPlosslessForm
    for cnd in PMs.conductor_ids(pm)
        pt = PMs.var(pm, nw, cnd)[:pt] = JuMP.@variable(pm.model,
            [(l,i,j) in PMs.ref(pm, nw, :arcs_from_trans)],
            base_name="$(nw)_$(cnd)_p_trans",
            start=0
        )
        if bounded
            for arc in PMs.ref(pm, nw, :arcs_from_trans)
                tr_id = arc[1]
                flow_lb  = -PMs.ref(pm, nw, :trans, tr_id, "rate_a")[cnd]
                flow_ub  =  PMs.ref(pm, nw, :trans, tr_id, "rate_a")[cnd]
                JuMP.set_lower_bound(pt[arc], flow_lb)
                JuMP.set_upper_bound(pt[arc], flow_ub)
            end
        end

        for (l,branch) in PMs.ref(pm, nw, :branch)
            if haskey(branch, "pf_start")
                f_idx = (l, branch["f_bus"], branch["t_bus"])
                JuMP.set_value(p[f_idx], branch["pf_start"])
            end
        end

        # this explicit type erasure is necessary
        p_expr = Dict{Any,Any}( ((l,i,j), pt[(l,i,j)]) for (l,i,j) in PMs.ref(pm, nw, :arcs_from_trans) )
        p_expr = merge(p_expr, Dict( ((l,j,i), -1.0*pt[(l,i,j)]) for (l,i,j) in PMs.ref(pm, nw, :arcs_from_trans)))
        PMs.var(pm, nw, cnd)[:pt] = p_expr
    end
end

"Do nothing, this model is symmetric"
function constraint_ohms_tp_yt_to(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.DCPlosslessForm
end

"Do nothing, this model is symmetric"
function constraint_ohms_tp_yt_to_on_off(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.DCPlosslessForm
end




### DC Power Flow Approximation ###

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] == -b*(t[f_bus] - t[t_bus])
```
"""
function constraint_ohms_tp_yt_from(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractDCPForm
    p_fr  = PMs.var(pm, n, c,  :p, f_idx)
    va_fr = [PMs.var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [PMs.var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]

    JuMP.@constraint(pm.model, p_fr == -sum(b[c,d]*(va_fr[c] - va_to[d]) for d in PMs.conductor_ids(pm)))
    # omit reactive constraint
end


"`-b*(t[f_bus] - t[t_bus] + vad_min*(1-branch_z[i])) <= p[f_idx] <= -b*(t[f_bus] - t[t_bus] + vad_max*(1-branch_z[i]))`"
function constraint_ohms_tp_yt_from_on_off(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractDCPForm
    p_fr  = PMs.var(pm, n, c,  :p, f_idx)
    va_fr = [PMs.var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [PMs.var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]
    z = [PMs.var(pm, n, d, :branch_z, i) for d in PMs.conductor_ids(pm)]

    JuMP.@constraint(pm.model, p_fr <= sum(-b[c,d]*(va_fr[d] - va_to[d] + vad_max[d]*(1-z[d])) for d in PMs.conductor_ids(pm)) )
    JuMP.@constraint(pm.model, p_fr >= sum(-b[c,d]*(va_fr[d] - va_to[d] + vad_min[d]*(1-z[d])) for d in PMs.conductor_ids(pm)) )
end



### Network Flow Approximation ###

"nothing to do, no voltage angle variables"
function constraint_tp_theta_ref(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, d) where T <: PMs.NFAForm
end

"nothing to do, no voltage angle variables"
function constraint_ohms_tp_yt_from(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.NFAForm
end

"nothing to do, this model is symmetric"
function constraint_ohms_tp_yt_to(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.NFAForm
end

"nothing to do, no voltage variables"
function constraint_tp_trans_voltage(pm::PMs.GenericPowerModel{T}, nw::Int, i::Int, f_bus::Int, t_bus::Int, tm::PMs.MultiConductorVector, Tv_fr, Tv_im, Cv_to) where T <: PMs.NFAForm
end

"nothing to do, this model is symmetric"
function constraint_tp_trans_flow(pm::PMs.GenericPowerModel{T}, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, tm::PMs.MultiConductorVector, Ti_fr, Ti_im, Cv_to) where T <: PMs.NFAForm
end

