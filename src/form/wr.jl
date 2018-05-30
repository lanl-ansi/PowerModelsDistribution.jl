# Three-phase specific constraints


""
function variable_tp_voltage(pm::GenericPowerModel{T}; kwargs...) where T <: PMs.AbstractWRForm
    variable_tp_voltage_magnitude_sqr(pm; kwargs...)
    variable_tp_voltage_product(pm; kwargs...)
end


""
function constraint_tp_voltage(pm::GenericPowerModel{T}, n::Int, h::Int) where T <: PMs.AbstractWRForm
    w  = var(pm, n,  :w)
    wr = var(pm, n, :wr)
    wi = var(pm, n, :wi)

    for g in h:length(PMs.phase_ids(pm))
        for (i,j) in ids(pm, n, :buspairs)
            InfrastructureModels.relaxation_complex_product(pm.model, w[(i,g)], w[(j,h)], wr[(i,j,h,g)], wi[(i,j,h,g)])
        end

        if g != h
            for i in ids(pm, n, :bus)
                InfrastructureModels.relaxation_complex_product(pm.model, w[(i,g)], w[(i,h)], wr[(i,i,h,g)], wi[(i,i,h,g)])
            end
        end
    end
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] ==        g/tm*w_fr[i] + (-g*tr+b*ti)/tm*(wr[i]) + (-b*tr-g*ti)/tm*(wi[i])
q[f_idx] == -(b+c/2)/tm*w_fr[i] - (-b*tr-g*ti)/tm*(wr[i]) + (-g*tr+b*ti)/tm*(wi[i])
```
"""
function constraint_ohms_tp_yt_from_on_off(pm::GenericPowerModel{T}, n::Int, h::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractWRForm
    p_fr = var(pm, n, h, :p, f_idx)
    q_fr = var(pm, n, h, :q, f_idx)
    w    = var(pm, n, :w)
    wr   = var(pm, n, :wr)
    wi   = var(pm, n, :wi)

    @constraint(pm.model, p_fr ==  ( g_fr[h]+g[h,h]) * w[(f_bus, h)] +
                                sum( g[h,i] * wr[(f_bus, f_bus, h, i)] +
                                     b[h,i] * wi[(f_bus, f_bus, h, i)] for i in PMs.phase_ids(pm) if i != h) +
                                sum(-g[h,i] * wr[(f_bus, t_bus, h, i)] +
                                    -b[h,i] * wi[(f_bus, t_bus, h, i)] for i in PMs.phase_ids(pm)) )
    @constraint(pm.model, q_fr == -( b_fr[h]+b[h,h]) * w[(f_bus, h)] -
                                sum( b[h,i] * wr[(f_bus, f_bus, h, i)] -
                                     g[h,i] * wi[(f_bus, f_bus, h, i)] for i in PMs.phase_ids(pm) if i != h) -
                                sum(-b[h,i] * wr[(f_bus, t_bus, h, i)] +
                                     g[h,i] * wi[(f_bus, t_bus, h, i)] for i in PMs.phase_ids(pm)) )
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==        g*w_to[i] + (-g*tr-b*ti)/tm*(wr[i]) + (-b*tr+g*ti)/tm*(-wi[i])
q[t_idx] == -(b+c/2)*w_to[i] - (-b*tr+g*ti)/tm*(wr[i]) + (-g*tr-b*ti)/tm*(-wi[i])
```
"""
function constraint_ohms_tp_yt_to_on_off(pm::GenericPowerModel{T}, n::Int, h::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractWRForm
    p_to = var(pm, n, h, :p, t_idx)
    q_to = var(pm, n, h, :q, t_idx)
    w    = var(pm, n, :w)
    wr   = var(pm, n, :wr)
    wi   = var(pm, n, :wi)

    @constraint(pm.model, p_to ==  ( g_to[h]+g[h,h]) * w[(t_bus, h)] +
                                sum( g[h,i] * wr[(t_bus, t_bus, h, i)] +
                                     b[h,i] *-wi[(t_bus, t_bus, h, i)] for i in PMs.phase_ids(pm) if i != h) +
                                sum(-g[h,i] * wr[(f_bus, t_bus, h, i)] +
                                    -b[h,i] *-wi[(f_bus, t_bus, h, i)] for i in PMs.phase_ids(pm)) )
    @constraint(pm.model, q_to == -( b_to[h]+b[h,h]) * w[(t_bus, h)] -
                                sum( b[h,i] * wr[(t_bus, t_bus, h, i)] -
                                     g[h,i] *-wi[(t_bus, t_bus, h, i)] for i in PMs.phase_ids(pm) if i != h) -
                                sum(-b[h,i] * wr[(f_bus, t_bus, h, i)] +
                                     g[h,i] *-wi[(f_bus, t_bus, h, i)] for i in PMs.phase_ids(pm)) )
end