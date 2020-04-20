import LinearAlgebra: tr


function variable_mc_branch_flow(pm::_PMs.AbstractWRModels; kwargs...)
    variable_mc_branch_flow_full(pm; kwargs...)
end

""
function variable_mc_voltage(pm::_PMs.AbstractWRModels; kwargs...)
    variable_mc_voltage_prod_hermitian(pm; kwargs...)
    variable_mc_voltage_cross_product(pm; kwargs...)
    constraint_mc_voltage_leaf_buses_psd(pm)
end

""
function variable_mc_voltage_cross_product(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true,  report::Bool=true)
    n_cond = 3

    bus_pairs = [key for (key, value) in _PMs.ref(pm, nw, :buspairs)]

    if bounded
        bound = Dict{eltype(bus_pairs), Array{Real,2}}()
        for (i,j) in bus_pairs
            bus_fr = _PMs.ref(pm, nw, :bus, i)
            bus_to = _PMs.ref(pm, nw, :bus, j)
            bound[(i,j)] = bus_fr["vmax"] * bus_to["vmax"]'
        end
        # create matrix variables
        (Wijr,Wiji) = variable_mx_complex(pm.model, bus_pairs, n_cond, n_cond;
            symm_bound=bound, name=("Wijr", "Wiji"), prefix="$nw")
    else
        (Wijr,Wiji) = variable_mx_complex(pm.model, bus_pairs, n_cond, n_cond;
            name=("Wijr", "Wiji"), prefix="$nw")
    end
    # save reference

    Wjir = Dict((j,i) =>  M' for ((i,j), M) in Wijr)
    Wjii = Dict((j,i) => -M' for ((i,j), M) in Wiji)

    _PMs.var(pm, nw)[:Wijr] = merge(Wijr, Wjir)
    _PMs.var(pm, nw)[:Wiji] = merge(Wiji, Wjii)

    # TODO
    # report && _PMs.sol_component_value_edge(pm, nw, :branch, :Wijr, :Wjir, _PMs.ref(pm, nw, :buspairs), _PMs.ref(pm, nw, :buspairs), Wijr)
    # report && _PMs.sol_component_value_edge(pm, nw, :branch, :Wiji, :Wjii, _PMs.ref(pm, nw, :buspairs), _PMs.ref(pm, nw, :buspairs), Wiji)
end



# """
# For the matrix KCL formulation, the generator needs an explicit current and
# power variable.
# """
# function variable_mc_generation(pm::_PMs.AbstractWRModels; nw=pm.cnw)
#     variable_mc_generation_current(pm; nw=nw)
#     variable_mc_generation_power(pm; nw=nw)
# end


"Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)"
function constraint_mc_ohms_yt_from(pm::_PMs.AbstractWRModels, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    (l,i,j) = f_idx
    Wr = _PMs.var(pm, n, :Wr, i)
    Wi = _PMs.var(pm, n, :Wi, i)

    Wijr = _PMs.var(pm, n, :Wijr)[(i,j)]
    Wiji = _PMs.var(pm, n, :Wiji)[(i,j)]

    P = _PMs.var(pm, n, :P)[f_idx]
    Q = _PMs.var(pm, n, :Q)[f_idx]

    JuMP.@constraint(pm.model, P .== Wr*(g_fr + g)' + Wi*(b_fr+b)' - Wijr*g' - Wiji*b')
    JuMP.@constraint(pm.model, Q .== Wi*(g_fr + g)' - Wr*(b_fr+b)' - Wiji*g' + Wijr*b')
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
#TODO   not assume it is symmetric?
function constraint_mc_ohms_yt_to(pm::_PMs.AbstractWRModels, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    constraint_mc_ohms_yt_from(pm, n, t_bus, f_bus, t_idx, f_idx, g, b, g_to, b_to, tr, ti, tm)
end

function constraint_mc_voltage_angle_difference(pm::_PMs.AbstractWRModels, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx

    ncnds = length(_PMs.conductor_ids(pm, n))

    w_fr = _PMs.var(pm, n, :w, f_bus)
    w_to = _PMs.var(pm, n, :w, t_bus)

    wr   = [_PMs.var(pm, n, :Wijr)[(f_bus, t_bus)][c,c] for c in 1:ncnds]
    wi   = [_PMs.var(pm, n, :Wiji)[(f_bus, t_bus)][c,c] for c in 1:ncnds]

    JuMP.@constraint(pm.model, wi .<= tan.(angmax).*wr)
    JuMP.@constraint(pm.model, wi .>= tan.(angmin).*wr)

    for c in 1:ncnds
        _PMs.cut_complex_product_and_angle_difference(pm.model, w_fr[c], w_to[c], wr[c], wi[c], angmin[c], angmax[c])
    end
end


""
function constraint_mc_theta_ref(pm::_PMs.AbstractWRModels, n::Int, i::Int, va_ref)
    nconductors = length(_PMs.conductor_ids(pm))

    Wr = _PMs.var(pm, n, :Wr)[i]
    Wi = _PMs.var(pm, n, :Wi)[i]

    beta = exp.(im.*va_ref)
    gamma = beta*beta'

    Wr_ref = real(gamma).*Wr[1,1]
    Wi_ref = imag(gamma).*Wr[1,1]
    JuMP.@constraint(pm.model, diag(Wr)[2:nconductors]        .== diag(Wr_ref)[2:nconductors]) # first equality is implied
    JuMP.@constraint(pm.model, _mat2utrivec!(Wr) .== _mat2utrivec!(Wr_ref))
    JuMP.@constraint(pm.model, _mat2utrivec!(Wi) .== _mat2utrivec!(Wi_ref))
end
