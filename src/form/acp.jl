# Three-phase specific constraints

"do nothing, this model does not have complex voltage constraints"
function constraint_tp_voltage(pm::GenericPowerModel{T}, n::Int, c::Int) where T <: PMs.AbstractACPForm
end

"Set loose voltage bounds; they will only be set if the voltage magnitude was unconstrained before."
function constraint_tp_voltage_mag_unbound(pm::GenericPowerModel{T}, i::Int, vmin::Float64, vmax::Float64; nw::Int=pm.cnw) where T <: PMs.AbstractACPForm
    for cnd in PMs.conductor_ids(pm)
        vmin_old = PMs.getlowerbound(var(pm, nw, cnd, :vm, i))
        vmax_old = PMs.getupperbound(var(pm, nw, cnd, :vm, i))
        # a lower bound of zero is considered to be unbound
        if vmin_old <= 0
            PMs.setlowerbound(var(pm, nw, cnd, :vm, i), vmin)
        end
        if vmax == Inf
            PMs.setupperbound(var(pm, nw, cnd, :vm, i), vmax)
        end
    end
end

""
function constraint_kcl_shunt_slack(pm::GenericPowerModel{T}, n::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: PMs.AbstractACPForm
    vm = var(pm, n, c, :vm, i)
    p_slack = var(pm, n, c, :p_slack, i)
    q_slack = var(pm, n, c, :q_slack, i)
    p = var(pm, n, c, :p)
    q = var(pm, n, c, :q)
    pg = var(pm, n, c, :pg)
    qg = var(pm, n, c, :qg)
    p_dc = var(pm, n, c, :p_dc)
    q_dc = var(pm, n, c, :q_dc)

    con(pm, n, c, :kcl_p)[i] = @constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*vm^2 + p_slack)
    con(pm, n, c, :kcl_q)[i] = @constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*vm^2 + q_slack)
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] ==  (g+g_fr)/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
q[f_idx] == -(b+b_fr)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
```
"""
function constraint_ohms_tp_yt_from(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) where T <: PMs.AbstractACPForm
    p_fr  = var(pm, n, c,  :p, f_idx)
    q_fr  = var(pm, n, c,  :q, f_idx)
    vm_fr = [var(pm, n, d, :vm, f_bus) for d in PMs.conductor_ids(pm)]
    vm_to = [var(pm, n, d, :vm, t_bus) for d in PMs.conductor_ids(pm)]
    va_fr = [var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]

    @NLconstraint(pm.model, p_fr ==  (g_fr[c]+g[c,c]) * vm_fr[c]^2 +
                                    sum( g[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                                         b[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in PMs.conductor_ids(pm) if d != c) +
                                    sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                        -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in PMs.conductor_ids(pm)) )
    @NLconstraint(pm.model, q_fr == -(b_fr[c]+b[c,c]) *vm_fr[c]^2 -
                                    sum( b[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                                         g[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in PMs.conductor_ids(pm) if d != c) -
                                    sum(-b[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                         g[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in PMs.conductor_ids(pm)) )
end


"""
Creates Ohms constraints for zero series impedance branches

```
p[f_idx] - g_fr/tm*v[f_bus]^2 + p[t_idx] - g_to*v[t_bus]^2 == 0
q[f_idx] + b_fr/tm*v[f_bus]^2 + q[t_idx] + b_to*v[t_bus]^2 == 0
vm_fr[c] == vm_to[c]
va_fr[c] == va_to[c]
"""
function constraint_ohms_tp_yt_from_impzero(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g_fr, b_fr, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractACPForm
    p_fr  = var(pm, n, c,  :p, f_idx)
    p_to  = var(pm, n, c,  :p, t_idx)
    q_fr  = var(pm, n, c,  :q, f_idx)
    q_to  = var(pm, n, c,  :q, t_idx)
    vm_fr = [var(pm, n, d, :vm, f_bus) for d in PMs.conductor_ids(pm)]
    vm_to = [var(pm, n, d, :vm, t_bus) for d in PMs.conductor_ids(pm)]
    va_fr = [var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]

    @NLconstraint(pm.model, p_fr - g_fr[c]*vm_fr[c]^2 + p_to - g_to[c]*vm_to[c]^2 == 0)
    @NLconstraint(pm.model, q_fr + b_fr[c]*vm_fr[c]^2 + q_to + b_to[c]*vm_to[c]^2 == 0)
    # link the voltages on both sides
    @constraint(pm.model, vm_fr[c] == vm_to[c])
    @constraint(pm.model, va_fr[c] == va_to[c])
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_ohms_tp_yt_to(pm::GenericPowerModel{T}, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm) where T <: PMs.AbstractACPForm
    p_to  = var(pm, n, c,  :p, t_idx)
    q_to  = var(pm, n, c,  :q, t_idx)
    vm_fr = [var(pm, n, d, :vm, f_bus) for d in PMs.conductor_ids(pm)]
    vm_to = [var(pm, n, d, :vm, t_bus) for d in PMs.conductor_ids(pm)]
    va_fr = [var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]

    @NLconstraint(pm.model, p_to == (g_to[c]+g[c,c])*vm_to[c]^2 +
                                   sum( g[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) +
                                        b[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in PMs.conductor_ids(pm) if d != c) +
                                   sum(-g[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                       -b[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in PMs.conductor_ids(pm)) )
    @NLconstraint(pm.model, q_to == -(b_to[c]+b[c,c])*vm_to[c]^2 -
                                   sum( b[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) -
                                        g[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in PMs.conductor_ids(pm) if d != c) -
                                   sum(-b[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                        g[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in PMs.conductor_ids(pm)) )
end


"""
```
p[f_idx] == z*(g/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
q[f_idx] == z*(-(b+c/2)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))
```
"""
function constraint_ohms_tp_yt_from_on_off(pm::GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractACPForm
    p_fr  = var(pm, n, c,  :p, f_idx)
    q_fr  = var(pm, n, c,  :q, f_idx)
    vm_fr = [var(pm, n, d, :vm, f_bus) for d in PMs.conductor_ids(pm)]
    vm_to = [var(pm, n, d, :vm, t_bus) for d in PMs.conductor_ids(pm)]
    va_fr = [var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]
    z = var(pm, n, c, :branch_z, i)

    @NLconstraint(pm.model, p_fr == z*((g_fr[c]+g[c,c]) * vm_fr[c]^2 +
                                        sum( g[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                                             b[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in PMs.conductor_ids(pm) if d != c) +
                                        sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                            -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in PMs.conductor_ids(pm))) )
    @NLconstraint(pm.model, q_fr == z*(-(b_fr[c]+b[c,c]) *vm_fr[c]^2 -
                                        sum( b[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                                             g[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in PMs.conductor_ids(pm) if d != c) -
                                        sum(-b[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                                             g[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in PMs.conductor_ids(pm))) )
end

"""
```
p[t_idx] == z*(g*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
q[t_idx] == z*(-(b+c/2)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))
```
"""
function constraint_ohms_tp_yt_to_on_off(pm::GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max) where T <: PMs.AbstractACPForm
    p_to  = var(pm, n, c,  :p, t_idx)
    q_to  = var(pm, n, c,  :q, t_idx)
    vm_fr = [var(pm, n, d, :vm, f_bus) for d in PMs.conductor_ids(pm)]
    vm_to = [var(pm, n, d, :vm, t_bus) for d in PMs.conductor_ids(pm)]
    va_fr = [var(pm, n, d, :va, f_bus) for d in PMs.conductor_ids(pm)]
    va_to = [var(pm, n, d, :va, t_bus) for d in PMs.conductor_ids(pm)]
    z = var(pm, n, c, :branch_z, i)

    @NLconstraint(pm.model, p_to == z*((g_to[c]+g[c,c])*vm_to[c]^2 +
                                        sum( g[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) +
                                             b[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in PMs.conductor_ids(pm) if d != c) +
                                        sum(-g[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                            -b[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in PMs.conductor_ids(pm))) )
    @NLconstraint(pm.model, q_to == z*(-(b_to[c]+b[c,c])*vm_to[c]^2 -
                                        sum( b[c,d]*vm_to[c]*vm_to[d]*cos(va_to[c]-va_to[d]) -
                                             g[c,d]*vm_to[c]*vm_to[d]*sin(va_to[c]-va_to[d]) for d in PMs.conductor_ids(pm) if d != c) -
                                        sum(-b[c,d]*vm_to[c]*vm_fr[d]*cos(va_to[c]-va_fr[d]) +
                                             g[c,d]*vm_to[c]*vm_fr[d]*sin(va_to[c]-va_fr[d]) for d in PMs.conductor_ids(pm))) )
end

""
function constraint_tp_trans_voltage(pm::GenericPowerModel, i::Int, f_bus::Int, t_bus::Int, tapset::MultiConductorVector, Tv_fr, Tv_im, Cv_to; nw::Int=pm.cnw)
    ncnd  = 3
    # intermediate bus voltage, for now ignore tap changer
    vm_to = [var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
    vm_im = @NLexpression(pm.model, [c in 1:ncnd], vm_to[c]*tapset[c]*Cv_to)
    va_im = [var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
    for c in 1:3
        if !haskey(var(pm, pm.cnw, c), :vm_trans)
            var(pm, pm.cnw, c)[:vm_trans] = Dict{Int, Any}()
            var(pm, pm.cnw, c)[:va_trans] = Dict{Int, Any}()
        end
        var(pm, pm.cnw, c, :vm_trans)[i] = vm_im[c]
        var(pm, pm.cnw, c, :va_trans)[i] = va_im[c]
    end
    # from side
    vm_fr = [var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    for n in 1:size(Tv_fr)[1]
        @NLconstraint(pm.model,
              sum(Tv_fr[n,c]*vm_fr[c]*cos(va_fr[c]) for c in 1:ncnd)
            ==sum(Tv_im[n,c]*vm_im[c]*cos(va_im[c]) for c in 1:ncnd)
        )
        @NLconstraint(pm.model,
              sum(Tv_fr[n,c]*vm_fr[c]*sin(va_fr[c]) for c in 1:ncnd)
            ==sum(Tv_im[n,c]*vm_im[c]*sin(va_im[c]) for c in 1:ncnd)
        )
    end
end

""
function constraint_tp_trans_voltage_var(pm::GenericPowerModel, i::Int, f_bus::Int, t_bus::Int, Tv_fr, Tv_im, Cv_to; nw::Int=pm.cnw)
    ncnd  = 3
    # intermediate bus voltage, for now ignore tap changer
    vm_to = [var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
    tap = [var(pm, nw, c, :tap)[i] for c in 1:ncnd]
    vm_im = @NLexpression(pm.model, [c in 1:ncnd], vm_to[c]*tap[c]*Cv_to)
    va_im = [var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
    for c in 1:3
        if !haskey(var(pm, pm.cnw, c), :vm_trans)
            var(pm, pm.cnw, c)[:vm_trans] = Dict{Int, Any}()
            var(pm, pm.cnw, c)[:va_trans] = Dict{Int, Any}()
        end
        var(pm, pm.cnw, c, :vm_trans)[i] = vm_im[c]
        var(pm, pm.cnw, c, :va_trans)[i] = va_im[c]
    end
    # from side
    vm_fr = [var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    for n in 1:size(Tv_fr)[1]
        @NLconstraint(pm.model,
              sum(Tv_fr[n,c]*vm_fr[c]*cos(va_fr[c]) for c in 1:ncnd)
            ==sum(Tv_im[n,c]*vm_im[c]*cos(va_im[c]) for c in 1:ncnd)
        )
        @NLconstraint(pm.model,
              sum(Tv_fr[n,c]*vm_fr[c]*sin(va_fr[c]) for c in 1:ncnd)
            ==sum(Tv_im[n,c]*vm_im[c]*sin(va_im[c]) for c in 1:ncnd)
        )
    end
end

""
function constraint_tp_trans_tap_fix(pm::GenericPowerModel, i::Int, tapfix::MultiConductorVector, tapset::MultiConductorVector; nw=pm.cnw)
    for (c,fixed) in enumerate(tapfix)
        if fixed
            @constraint(pm.model, var(pm, nw, c, :tap)[i]==tapset[c])
        end
    end
end

""
function constraint_tp_trans_flow(pm::GenericPowerModel, i::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, Ti_fr, Ti_to; nw::Int=pm.cnw)
    # the intermediate bus voltage is saved as an expression
    ncnd  = 3
    vm_im = [var(pm, nw, c, :vm_trans, i) for c in 1:ncnd]
    va_im = [var(pm, nw, c, :va_trans, i) for c in 1:ncnd]
    # power is unaffected by tap-changer
    p_to = [var(pm, nw, c, :p_trans, t_idx) for c in 1:ncnd]
    q_to = [var(pm, nw, c, :q_trans, t_idx) for c in 1:ncnd]
    # from side variables
    vm_fr = [var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
    va_fr = [var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
    p_fr = [var(pm, nw, c, :p_trans, f_idx) for c in 1:ncnd]
    q_fr = [var(pm, nw, c, :q_trans, f_idx) for c in 1:ncnd]
    for n in 1:size(Ti_fr)[1]
        # i_fr_re[c] = 1/vm_fr[c]*(p_fr[c]*cos(va_fr[c])+q_fr[c]*sin(va_fr[c]))
        # i_fr_im[c] = 1/vm_fr[c]*(p_fr[c]*sin(va_fr[c])-q_fr[c]*cos(va_fr[c]))
        # i_to_re[c] = 1/vm_to[c]*(p_to[c]*cos(va_im[c])+q_to[c]*sin(va_im[c]))
        # i_to_im[c] = 1/vm_to[c]*(p_to[c]*sin(va_im[c])-q_to[c]*cos(va_im[c]))
        @NLconstraint(pm.model,
              sum(Ti_fr[n,c]*
                    1/vm_fr[c]*(p_fr[c]*cos(va_fr[c])+q_fr[c]*sin(va_fr[c])) # i_fr_re[c]
              for c in 1:ncnd)
            + sum(Ti_to[n,c]*
                    1/vm_im[c]*(p_to[c]*cos(va_im[c])+q_to[c]*sin(va_im[c])) # i_to_re[c]
              for c in 1:ncnd)
            == 0
        )
        @NLconstraint(pm.model,
              sum(Ti_fr[n,c]*
                    1/vm_fr[c]*(p_fr[c]*sin(va_fr[c])-q_fr[c]*cos(va_fr[c])) # i_fr_im[c]
              for c in 1:ncnd)
            + sum(Ti_to[n,c]*
                    1/vm_im[c]*(p_to[c]*sin(va_im[c])-q_to[c]*cos(va_im[c])) # i_to_im[c]
              for c in 1:ncnd)
            == 0
        )
    end
end

# ""
# function constraint_tp_trans_vartap_fix(pm::GenericPowerModel, i::Int, tapset, tapfix; nw::Int=pm.cnw)
#     for cnd in 1:3
#         if tapfix[cnd]
#             @constraint(pm.model, var(pm, nw, :tap)[(i,cnd)] == tapset[cnd])
#         end
#     end
# end
#
# ""
# function constraint_tp_trans_voltage_vartap(pm::GenericPowerModel, i::Int, f_bus::Int, t_bus::Int, bkv_fr::Float64, bkv_to::Float64, vnom_kv_fr::Float64, vnom_kv_to::Float64; nw::Int=pm.cnw)
#     ncnd  = 3
#     # vm_fr = append_zero([var(pm, nw, c, :vm, f_bus) for c in 1:ncnd])
#     # va_fr = append_zero([var(pm, nw, c, :va, f_bus) for c in 1:ncnd])
#     # vm_to = append_zero([var(pm, nw, c, :vm, t_bus) for c in 1:ncnd])
#     # va_to = append_zero([var(pm, nw, c, :va, t_bus) for c in 1:ncnd])
#     vm_fr = [var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
#     va_fr = [var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
#     vm_to = [var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
#     va_to = [var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
#     for n in 1:ncnd
#         @NLconstraint(pm.model,
#             bkv_fr/bkv_to*vm_fr[n]
#             # 4 wire: bkv_fr/bkv_to*(vm_fr[n]-vm_fr[end])
#             ==
#             var(pm, nw, :tap)[(i,n)]*vnom_kv_fr/vnom_kv_to*vm_to[n]
#             # 4 wire: var(pm, nw, :tap)[(i,n)]*vnom_kv_fr/vnom_kv_to*(vm_to[n]-vm_to[end])
#         )
#         @constraint(pm.model,
#             va_fr[n]
#             # 4 wire: va_fr[n]-va_fr[end]
#             ==
#             va_to[n]
#             # 4 wire: va_to[n]-va_to[end]
#         )
#     end
# end
#
# ""
# function constraint_tp_trans_power_vartap(pm::GenericPowerModel, i::Int, f_bus::Int, t_bus::Int, f_idx, t_idx, bkv_fr::Float64, bkv_to::Float64, f_vnom_kv::Float64, t_vnom_kv::Float64; nw::Int=pm.cnw)
#     ncnd  = 3
#     vm_fr = [var(pm, nw, c, :vm, f_bus) for c in 1:ncnd]
#     va_fr = [var(pm, nw, c, :va, f_bus) for c in 1:ncnd]
#     vm_to = [var(pm, nw, c, :vm, t_bus) for c in 1:ncnd]
#     va_to = [var(pm, nw, c, :va, t_bus) for c in 1:ncnd]
#     p_fr = [var(pm, nw, c, :p_trans, f_idx) for c in 1:ncnd]
#     q_fr = [var(pm, nw, c, :q_trans, f_idx) for c in 1:ncnd]
#     p_to = [var(pm, nw, c, :p_trans, t_idx) for c in 1:ncnd]
#     q_to = [var(pm, nw, c, :q_trans, t_idx) for c in 1:ncnd]
#
#     for n in 1:ncnd
#         @NLconstraint(pm.model,
#             var(pm, nw, :tap)[(i,n)]*bkv_to/bkv_fr*f_vnom_kv/t_vnom_kv
#             *vm_to[n]*(p_fr[n]*cos(va_to[n])-q_fr[n]*sin(va_to[n]))
#             +vm_fr[n]*(p_to[n]*cos(va_fr[n])-q_to[n]*sin(va_fr[n]))
#             == 0
#         )
#         @NLconstraint(pm.model,
#             var(pm, nw, :tap)[(i,n)]*bkv_to/bkv_fr*f_vnom_kv/t_vnom_kv
#             *vm_to[n]*(p_fr[n]*sin(va_to[n])+q_fr[n]*cos(va_to[n]))
#             +vm_fr[n]*(p_to[n]*sin(va_fr[n])+q_to[n]*cos(va_fr[n]))
#             == 0
#         )
#     end
# end
