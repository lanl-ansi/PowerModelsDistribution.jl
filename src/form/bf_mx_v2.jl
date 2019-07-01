#TODO
#- Should unbounded be supported as well by variable_mx? If so, redesign signatures
#- Fix matrix variable dict keys once merging
#- Add formulation types
#- Add matrix KCL
#- Allow to pass vectors for upperbound/lowerbound for Hermitian
#- Full matrix fixed diagonal
import LinearAlgebra: diag, diagm

function variable_tp_branch_current_v2(pm::PMs.GenericPowerModel{T}; kwargs...) where T <: AbstractUBFForm
    variable_tp_branch_series_current_prod_hermitian_v2(pm; kwargs...)
end


function variable_tp_voltage_v2(pm::PMs.GenericPowerModel{T}; kwargs...) where T <: AbstractUBFForm
    variable_tp_voltage_prod_hermitian_v2(pm; kwargs...)
end


function variable_tp_voltage_prod_hermitian_v2(pm::PMs.GenericPowerModel{T}; n_cond::Int=3, nw::Int=pm.cnw, bounded = true) where T <: AbstractUBFForm
    bus_ids = collect(PMs.ids(pm, nw, :bus))

    if bounded
        # get bounds
        vmax = Dict([(id, PMs.ref(pm, nw, :bus, id, "vmax").values) for id in bus_ids])
        vmin = Dict([(id, PMs.ref(pm, nw, :bus, id, "vmin").values) for id in bus_ids])
        # create bounded Hermitian matrix variables
        (Wre,Wim) = variable_mx_hermitian_sqrt_bounds(pm, bus_ids, n_cond,
            vmax, vmin; name="W", prefix="$nw")
    else
        # create unbounded Hermitian matrix variables
        (Wre,Wim) = variable_mx_hermitian(pm, bus_ids, n_cond; name="W", prefix="$nw")
    end

    # save references in dict
    PMs.var(pm, nw)[:W_re] = Wre
    PMs.var(pm, nw)[:W_im] = Wim
    for c in 1:n_cond
        PMs.var(pm, nw, c)[:w] = Dict{Int, Any}([(id, Wre[id][c,c]) for id in bus_ids])
    end
end


function variable_tp_branch_series_current_prod_hermitian_v2(pm::PMs.GenericPowerModel{T}; n_cond::Int=3, nw::Int=pm.cnw, bounded = true) where T <: AbstractUBFForm
    branches = PMs.ref(pm, nw, :branch)
    buses = PMs.ref(pm, nw, :bus)

    branch_ids = collect(keys(branches))

    if bounded
        # calculate max series current for each branch
        cmax = Dict{eltype(branch_ids), Array{Real,1}}()
        for (key, branch) in branches
            bus_fr = buses[branch["f_bus"]]
            bus_to = buses[branch["t_bus"]]

            vmin_fr = bus_fr["vmin"].values
            vmin_to = bus_to["vmin"].values

            vmax_fr = bus_fr["vmax"].values
            vmax_to = bus_to["vmax"].values

            # assumed to be matrices already
            # temportary fix by shunts_diag2mat!
            y_fr = branch["g_fr"].values + im* branch["b_fr"].values
            y_to = branch["g_to"].values + im* branch["b_to"].values

            smax = branch["rate_a"].values
            cmaxfr = smax./vmin_fr + abs.(y_fr)*vmax_fr
            cmaxto = smax./vmin_to + abs.(y_to)*vmax_to

            cmax[key] = max.(cmaxfr, cmaxto)
        end
        # create matrix variables
        (Lre,Lim) = variable_mx_hermitian_sqrt_bounds(pm, branch_ids, n_cond, cmax; name="CC", prefix="$nw")
    else
        (Lre,Lim) = variable_mx_hermitian(pm, branch_ids, n_cond; name="CC", prefix="$nw")
    end

    # save reference
    PMs.var(pm, nw)[:CC_re] = Lre
    PMs.var(pm, nw)[:CC_im] = Lim
    for c in 1:n_cond
        PMs.var(pm, nw, c)[:cm] = Dict([(id, Lre[id][c,c]) for id in branch_ids])
    end
end


""
function variable_tp_branch_flow_v2(pm::PMs.GenericPowerModel{T}; n_cond::Int=3, nw::Int=pm.cnw, bounded = true) where T <: AbstractUBFForm
    # calculate S bound
    branch_arcs = PMs.ref(pm, nw, :arcs)
    bound = Dict{eltype(branch_arcs), Array{Real,2}}()
    for (l,i,j) in branch_arcs
        vmin = PMs.ref(pm, nw, :bus, i)["vmin"].values
        vmax = PMs.ref(pm, nw, :bus, i)["vmax"].values
        smax = PMs.ref(pm, nw, :branch, l)["rate_a"].values
        cmax = smax./vmin
        bound[(l,i,j)] = vmax*cmax'
    end
    # create matrix variables
    (P,Q) = variable_mx_complex(pm, branch_arcs, n_cond, n_cond, bound; name=("P", "Q"), prefix="$nw")
    # save reference
    PMs.var(pm, nw)[:P_mx] = P
    PMs.var(pm, nw)[:Q_mx] = Q
    for c in 1:n_cond
        PMs.var(pm, nw, c)[:p] = Dict([(id,P[id][c,c]) for id in branch_arcs])
        PMs.var(pm, nw, c)[:q] = Dict([(id,Q[id][c,c]) for id in branch_arcs])
    end
end


"Converts all shunt values from vector to matrix."
function shunts_diag2mat!(pm)
    for (nw, ref_nw) in  PMs.nws(pm)
        for (id, br) in ref_nw[:branch]
            for key in ["g_fr", "g_to", "b_fr", "b_to"]
                if !isa(br[key], PMs.MultiConductorMatrix)
                    br[key] = LinearAlgebra.diagm(0=>br[key])
                end
            end
        end
    end
end


function constraint_tp_voltage_magnitude_difference_v2(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    branch = PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    r = branch["br_r"].values
    x = branch["br_x"].values
    g_sh_fr = branch["g_fr"].values
    b_sh_fr = branch["b_fr"].values
    tm = branch["tap"].values

    constraint_tp_voltage_magnitude_difference(pm, nw, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm)
end


function constraint_tp_branch_current_v2(pm::PMs.GenericPowerModel{T}, i::Int; nw::Int=pm.cnw) where T <: AbstractUBFForm
    branch = PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    g_sh_fr = branch["g_fr"].values
    b_sh_fr = branch["b_fr"].values

    constraint_tp_branch_current(pm, nw, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
end


function constraint_tp_flow_losses_v2(pm::PMs.GenericPowerModel, i::Int; nw::Int=pm.cnw)
    branch = PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    r = branch["br_r"].values
    x = branch["br_x"].values
    g_sh_fr = branch["g_fr"].values
    g_sh_to = branch["g_to"].values
    b_sh_fr = branch["b_fr"].values
    b_sh_to = branch["b_to"].values

    constraint_tp_flow_losses(pm::PMs.GenericPowerModel, nw, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to)
end


# MATRIX POWER BALANCE
function variable_tp_generation_power_mx(pm::PMs.GenericPowerModel; nw=pm.cnw)
    gen_ids = collect(PMs.ids(pm, nw, :gen))
    ncnds = length(PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(gen_ids), Array{Real,2}}()
    for (id, gen) in PMs.ref(pm, nw, :gen)
        bus_id = gen["gen_bus"]
        vmax = PMs.ref(pm, nw, :bus, bus_id)["vmax"].values
        vmin = PMs.ref(pm, nw, :bus, bus_id)["vmin"].values
        pmax = max(abs.(gen["pmax"].values), abs.(gen["pmin"].values))
        qmax = max(abs.(gen["qmax"].values), abs.(gen["qmin"].values))
        smax = sqrt.(pmax.^2 + qmax.^2)
        cmax = smax./vmin
        bound[id] = vmax*cmax'
    end
    # create matrix variables
    (Pg,Qg) = variable_mx_complex(pm, gen_ids, ncnds, bound; name=("Pg", "Qg"), prefix="$nw")
    # save references
    PMs.var(pm, nw)[:Pg] = Pg
    PMs.var(pm, nw)[:Qg] = Qg
    for c in 1:ncnds
        PMs.var(pm, nw, c)[:pg] =Dict([(id, Pg[id][c,c]) for id in gen_ids])
        PMs.var(pm, nw, c)[:qg] =Dict([(id, Qg[id][c,c]) for id in gen_ids])
    end
end


# MATRIX POWER BALANCE
function variable_tp_generation_current_mx(pm::PMs.GenericPowerModel; nw=pm.cnw)
    gen_ids = collect(PMs.ids(pm, nw, :gen))
    ncnds = length(PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(gen_ids), Array{Real,2}}()
    for (id, gen) in PMs.ref(pm, nw, :gen)
        bus_id = gen["gen_bus"]
        vmax = PMs.ref(pm, nw, :bus, bus_id)["vmax"].values
        vmin = PMs.ref(pm, nw, :bus, bus_id)["vmin"].values
        pmax = max(abs.(gen["pmax"].values), abs.(gen["pmin"].values))
        qmax = max(abs.(gen["qmax"].values), abs.(gen["qmin"].values))
        smax = sqrt.(pmax.^2 + qmax.^2)
        bound[id] = cmax*cmax'
    end
    # create matrix variables
    (Lre,Lim) = variable_mx_hermitian(pm, gen_ids, ncnds, bound; name="Ld", prefix="$nw")
    # save references
    PMs.var(pm, nw)[:Lgre] = Lre
    PMs.var(pm, nw)[:Lgre] = Lim
end


function variable_tp_load_power_mx(pm::PMs.GenericPowerModel; nw=pm.cnw)
    load_ids = collect(PMs.ids(pm, nw, :load))
    ncnds = length(PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(load_ids), Array{Real,2}}()
    for (id, load) in PMs.ref(pm, nw, :load)
        bus_id = load["load_bus"]
        vmax = PMs.ref(pm, nw, :bus, bus_id)["vmax"].values
        vmin = PMs.ref(pm, nw, :bus, bus_id)["vmin"].values
        # this presumes constant power, wye loads!
        @assert(load["model"]=="constant_power")
        #TODO extend to other load models
        pmax = abs.(load["pd"].values)
        qmax = abs.(load["qd"].values)
        smax = sqrt.(pmax.^2 + qmax.^2)
        cmax = smax./vmin
        bound[id] = vmax*cmax'
    end
    # create matrix variables
    (Pd,Qd) = variable_mx_complex(pm, load_ids, ncnds, bound; name=("Pd", "Qd"), prefix="$nw")
    # save references
    PMs.var(pm, nw)[:Pd] = Pd
    PMs.var(pm, nw)[:Qd] = Qd
    for c in 1:ncnds
        PMs.var(pm, nw, c)[:pd] =Dict([(id, Pd[id][c,c]) for id in gen_ids])
        PMs.var(pm, nw, c)[:qd] =Dict([(id, Qd[id][c,c]) for id in gen_ids])
    end
end


function variable_tp_load_current_mx(pm::PMs.GenericPowerModel; nw=pm.cnw)
    load_ids = collect(PMs.ids(pm, nw, :load))
    ncnds = length(PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(load_ids), Array{Real,2}}()
    for (id, load) in PMs.ref(pm, nw, :load)
        bus_id = load["load_bus"]
        vmax = PMs.ref(pm, nw, :bus, bus_id)["vmax"].values
        vmin = PMs.ref(pm, nw, :bus, bus_id)["vmin"].values
        # this presumes constant power, wye loads!
        @assert(load["model"]=="constant_power")
        #TODO extend to other load models
        pmax = abs.(load["pd"].values)
        qmax = abs.(load["qd"].values)
        smax = sqrt.(pmax.^2 + qmax.^2)
        cmax = smax./vmin
        bound[id] = cmax*cmax'
    end
    # create matrix variables
    (Ldre, Ldim) = variable_mx_hermitian(pm, load_ids, ncnds, bound; name="Ld", prefix="$nw")
    # save references
    PMs.var(pm, nw)[:Ldre] = Ldre
    PMs.var(pm, nw)[:Ldim] = Ldim
end


function constraint_tp_generation_current(pm::PMs.GenericPowerModel{T}, gen_id::Int; nw::Int=pm.cnw) where T <: AbstractUBFForm
    Pg = PMs.var(pm, nw, :Pg, gen_id)
    Qg = PMs.var(pm, nw, :Qg, gen_id)
    bus_id = PMs.ref(pm, nw, :gen)["gen_bus"]
    W_re = PMs.var(pm, nw, :W_re, bus_id)
    W_im = PMs.var(pm, nw, :W_im, bus_id)
    Ldre = PMs.var(pm, nw, :Lgre, load_id)
    Ldim = PMs.var(pm, nw, :Lgim, load_id)
    constraint_SWL_psd(Pg, Qg, W_re, W_im, Lgre, Lgim)
end


function constraint_tp_load_current(pm::PMs.GenericPowerModel{T}, load_id::Int; nw::Int=pm.cnw) where T <: AbstractUBFForm
    Pd = PMs.var(pm, nw, :Pd, load_id)
    Qd = PMs.var(pm, nw, :Qd, load_id)
    bus_id = PMs.ref(pm, nw, :load)["load_bus"]
    W_re = PMs.var(pm, nw, :W_re, bus_id)
    W_im = PMs.var(pm, nw, :W_im, bus_id)
    Ldre = PMs.var(pm, nw, :Ldre, load_id)
    Ldim = PMs.var(pm, nw, :Ldim, load_id)
    constraint_SWL_psd(pm, Pd, Qd, W_re, W_im, Ldre, Ldim)
end


function constraint_tp_voltage_psd(pm::PMs.GenericPowerModel; nw=pm.cnw)
    buses_covered = [i for (l,i,j) in PMs.ref(pm, nw, :arcs)]
    buses_psd = [i for i in PMs.ids(pm, nw, :bus) if !(i in buses_covered)]
    for bus_id in buses_psd
        W_re = PMs.var(pm, nw, :W_re, bus_id)
        W_im = PMs.var(pm, nw, :W_im, bus_id)
        constraint_M_psd(W_re, W_im)
    end
end


function constraint_SWL_psd(pm::PMs.GenericPowerModel, P, Q, W_re, W_im, L_re, L_im)
    M_re = [W_re P; P' L_re]
    M_im = [W_im Q; -Q' L_im]
    constraint_M_psd(pm, M_re, M_im)
end


function constraint_M_psd(pm, M_re, M_im)
    JuMP.@constraint(pm.model, [M_re -M_im; M_im M_re] in JuMP.PSDCone())
end


# UTILITY MATRIX VARIABLE FUNCTIONS
"""
This function creates a set of real matrix variables of size NxM,
indexed over the elements of the indices argument.
The upper and lower bounds have to be specified,
and are dictionaries with the indices as keys and the matrix bounds as values.
The name and prefix arguments will be combined into the base_name argument for
JuMP; the prefix will typically be the network number nw.
Instead of sequentially creating the matrix variables, the elements of the
matrices are created sequentially for all matrices at once. I.e., we loop
over the elements, and not over the indices. This is needed so that the
variable names printed by JuMP are in line with the current design.

Returns a dictionary of (index, matrix  variable) pairs
"""
function variable_mx_real(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int, M::Int,
        lower_bound::Dict{T,Array{LB,2}}, upper_bound::Dict{T,Array{UB,2}};
        name="", prefix="") where {T, LB<:Real, UB<:Real}
    # the output is a dictionary of (index, matrix) pairs
    dict_mat_vars_temp = Dict{T, Array{Any, 2}}([(i,fill(missing,N,M)) for i in indices])
    for n in 1:N
        for m in 1:M
            varname = isempty(prefix) ? "$(n)$(m)_$name" : "$(prefix)_$(n)$(m)_$name"
            # create the element (n,m) for all indices
            mat_nm = JuMP.@variable(pm.model, [index in indices],
                base_name=varname,
                lower_bound=lower_bound[index][n,m],
                upper_bound=upper_bound[index][n,m]
            )
            # unpack element (n,m) to the correct place in the ouput dict
            for index in indices
                dict_mat_vars_temp[index][n,m] = mat_nm[index]
            end
        end
    end
    dict_mat_vars = Dict{T, Array{JuMP.VariableRef, 2}}([(id, arr) for (id, arr) in dict_mat_vars_temp])

    return dict_mat_vars
end


"""
Same as variable_mx_real, but has to be square and the diagonal of the matrix
variables consists of the constants passed as the diag argument.
The diag argument is a dictionary of (index, 1d-array) pairs.
Useful for power matrices with specified diagonals (constant power wye loads).
If not specified, the diagonal elements are set to zero.
"""
function variable_mx_real_const_diag(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int,
        lower_bound::Dict{T,Array{LB,2}}, upper_bound::Dict{T,Array{UB,2}};
        diag::Dict{T,Array{D,1}}=Dict([(i, zeros(N)) for i in indices]),
        name="", prefix="") where {T, LB<:Real, UB<:Real, D<:Real}
    # the output is a dictionary of (index, matrix) pairs
    dict_mat_vars = Dict{T, Array{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, 2}}([(index, zeros(N,N)) for index in indices])
    for n in 1:N
        for m in 1:N
            if n==m
                for index in indices
                    dict_mat_vars[index][n,n] = diag[index][n]
                end
            else
                varname = isempty(prefix) ? "$(n)$(m)_$name" : "$(prefix)_$(n)$(m)_$name"
                # create the element (n,m) for all indices
                mat_nm = JuMP.@variable(pm.model, [index in indices],
                    base_name=varname,
                    lower_bound=lower_bound[index][n,m],
                    upper_bound=upper_bound[index][n,m]
                )
                # unpack element (n,m) to the correct place in the ouput dict
                for index in indices
                    dict_mat_vars[index][n,m] = mat_nm[index]
                end
            end
        end
    end
    return dict_mat_vars
end


"""
Shorthand to create two real matrix variables, where the first is the real part
and the second the imaginary part.
If the name argument is a String, it will be suffixed with 're' and  'im'.
It is possible to  specify the names of the real and imaginary part directly
as a Tuple as well (to achieve P and Q instead of Sre and Sim for example).
"""
function variable_mx_complex(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int, M::Int,
        lower_bound::Dict{T,Array{LB,2}}, upper_bound::Dict{T,Array{UB,2}};
        name::Union{String, Tuple{String,String}}="", prefix="") where {T, LB<:Real,UB<:Real}
    name_real = isa(name, Tuple) ? name[1] : "$(name)re"
    name_imag = isa(name, Tuple) ? name[2] : "$(name)im"
    Mre = variable_mx_real(pm, indices, N, M, lower_bound, upper_bound, prefix=prefix,
        name=name_real
    )
    Mim = variable_mx_real(pm, indices, N, M, lower_bound, upper_bound;
        prefix=prefix, name=name_imag)
    return (Mre,Mim)
end


"""
Same as variable_mx_complex, but square and with constant diagonals.
If not specified, the diagonal elements are set to zero.
"""
function variable_mx_complex_const_diag(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int,
        lower_bound::Dict{T,Array{LB,2}}, upper_bound::Dict{T,Array{UB,2}};
        diag_re::Dict{T,Array{DR,1}}=Dict([(i, zeros(N)) for i in indices]),
        diag_im::Dict{T,Array{DI,1}}=Dict([(i, zeros(N)) for i in indices]),
        name::Union{String, Tuple{String,String}}="", prefix="") where {T, LB<:Real,UB<:Real, DR<:Real, DI<:Real}
    name_real = isa(name, Tuple) ? name[1] : "$(name)re"
    name_imag = isa(name, Tuple) ? name[2] : "$(name)im"
    Mre = variable_mx_real_const_diag(pm, indices, N, lower_bound, upper_bound, diag_re;
        prefix=prefix, name=name_real
    )
    Mim = variable_mx_real_const_diag(pm, indices, N, lower_bound, upper_bound, diag_im;
        prefix=prefix, name=name_imag)
    return (Mre,Mim)
end


"""
Same as variable_mx_complex, but the diagonal of the matrix variables consists
of the constants passed as the diag_re and diag_im argument. The diag argument
is a dictionary of (index, 1d-array) pairs.
Useful for power matrices with specified diagonals (constant power wye loads).
"""
function variable_mx_complex(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int, M::Int,
        bound::Dict{T,Array{B,2}};
        name::Union{String, Tuple{String,String}}="", prefix="") where {T, B<:Real}
    upper_bound = bound
    lower_bound = Dict([(k,-v) for (k,v) in bound])
    return variable_mx_complex(pm, indices, N, M, lower_bound, upper_bound;
        name=name, prefix=prefix)
end


"""
Same as variable_mx_complex_const_diag, but with symmetric bounds.
"""
function variable_mx_complex_const_diag(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int,
        bound::Dict{T,Array{B,2}}, diag_re::Dict{T,Array{DR,1}}, diag_im::Dict{T,Array{DI,1}};
        name::Union{String, Tuple{String,String}}="", prefix="") where {T, B<:Real, DR<:Real, DI<:Real}
    upper_bound = bound
    lower_bound = Dict([(k,-v) for (k,v) in bound])
    return variable_mx_complex_const_diag(pm, indices, N, lower_bound, upper_bound, diag_re, diag_im;
        name=name, prefix=prefix)
end


"""
Same as variable_mx_real, but adds symmetry structure
"""
function variable_mx_symmetric(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int,
        lower_bound::Dict{T,Array{LB,2}}, upper_bound::Dict{T,Array{UB,2}};
        name="", prefix="") where {T, LB<:Real, UB<:Real}
    # the output is a dictionary of (index, matrix) pairs
    # the data type of the matrix has to be JuMP.GenericAffExpr, because it
    # will also contain inverted variables (upper triangle) and potentially
    # constants (zero) on the diagonal
    dict_mat_vars = Dict{T, Array{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, 2}}([(index, zeros(N,N)) for index in indices])
    for n in 1:N
        for m in 1:n
            varname = isempty(prefix) ? "$(n)$(m)_$name" : "$(prefix)_$(n)$(m)_$name"
            # create the lower triangle element (n,m) for all indices
            mat_nm = JuMP.@variable(pm.model, [index in indices],
                base_name=varname,
                lower_bound=lower_bound[index][n,m],
                upper_bound=upper_bound[index][n,m]
            )
            # unpack element (n,m) to the correct place in the ouput dict
            for index in indices
                # a diagonal element only appears on the diagonal
                if n==m
                    dict_mat_vars[index][n,n] = mat_nm[index]
                # a lower triangle element apears in the lower triangle
                # and also in the upper triangle
                else
                    dict_mat_vars[index][n,m] = mat_nm[index]
                    dict_mat_vars[index][m,n] = mat_nm[index]
                end
            end
        end
    end
    return dict_mat_vars
end


"""
Alternative method of variable_mx_symmetric,
that does not take any bounds, and creates unbounded matrix variables.
"""
function variable_mx_symmetric(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int;
        name="", prefix="") where {T, LB<:Real, UB<:Real}
    # the output is a dictionary of (index, matrix) pairs
    # the data type of the matrix has to be JuMP.GenericAffExpr, because it
    # will also contain inverted variables (upper triangle) and potentially
    # constants (zero) on the diagonal
    dict_mat_vars = Dict{T, Array{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, 2}}([(index, zeros(N,N)) for index in indices])
    for n in 1:N
        for m in 1:n
            varname = isempty(prefix) ? "$(n)$(m)_$name" : "$(prefix)_$(n)$(m)_$name"
            # create the lower triangle element (n,m) for all indices
            mat_nm = JuMP.@variable(pm.model, [index in indices],
                base_name=varname
            )
            # unpack element (n,m) to the correct place in the ouput dict
            for index in indices
                # a diagonal element only appears on the diagonal
                if n==m
                    dict_mat_vars[index][n,n] = mat_nm[index]
                # a lower triangle element apears in the lower triangle
                # and also in the upper triangle
                else
                    dict_mat_vars[index][n,m] = mat_nm[index]
                    dict_mat_vars[index][m,n] = mat_nm[index]
                end
            end
        end
    end
    return dict_mat_vars
end


"""
Same as variable_mx_real, but adds skew-symmetry structure.
"""
function variable_mx_skewsymmetric(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int,
        lower_bound::Dict{T,Array{LB,2}}, upper_bound::Dict{T,Array{UB,2}};
        name="", prefix="") where {T, LB<:Real, UB<:Real}
    # the output is a dictionary of (index, matrix) pairs
    # the data type of the matrix has to be JuMP.GenericAffExpr, because it
    # will also contain inverted variables (upper triangle) and potentially
    # constants (zero) on the diagonal
    dict_mat_vars = Dict{T, Array{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, 2}}([(index, zeros(N,N)) for index in indices])
    for n in 1:N
        for m in 1:n
            varname = isempty(prefix) ? "$(n)$(m)_$name" : "$(prefix)_$(n)$(m)_$name"
            # create the lower triangle element (n,m) for all indices
            # if diagonal element (n,m) is zero
            if m==n
                mat_nm = Dict{T, Float64}([(index, 0.0) for index in indices])
            # otherwise, create JuMP variables
            else
                mat_nm = JuMP.@variable(pm.model, [index in indices],
                    base_name=varname,
                    lower_bound=lower_bound[index][n,m],
                    upper_bound=upper_bound[index][n,m]
                )
            end
            # unpack element (n,m) to the correct place in the ouput dict
            for index in indices
                # a diagonal element only appears on the diagonal
                if n==m
                    dict_mat_vars[index][n,n] = mat_nm[index]
                # a lower triangle element apears in the lower triangle
                # and also in the upper triangle (but inverted)
                else
                    dict_mat_vars[index][n,m] =  mat_nm[index]
                    dict_mat_vars[index][m,n] = -mat_nm[index]
                end
            end
        end
    end
    return dict_mat_vars
end


"""
Alternative method of variable_mx_skewsymmetric,
that does not take any bounds, and creates unbounded matrix variables.
"""
function variable_mx_skewsymmetric(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int;
        name="", prefix="") where {T, LB<:Real, UB<:Real}
    # the output is a dictionary of (index, matrix) pairs
    # the data type of the matrix has to be JuMP.GenericAffExpr, because it
    # will also contain inverted variables (upper triangle) and potentially
    # constants (zero) on the diagonal
    dict_mat_vars = Dict{T, Array{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, 2}}([(index, zeros(N,N)) for index in indices])
    for n in 1:N
        for m in 1:n
            varname = isempty(prefix) ? "$(n)$(m)_$name" : "$(prefix)_$(n)$(m)_$name"
            # create the lower triangle element (n,m) for all indices
            # if diagonal element (n,m) is zero
            if m==n
                mat_nm = Dict{T, Float64}([(index, 0.0) for index in indices])
            # otherwise, create JuMP variables
            else
                mat_nm = JuMP.@variable(pm.model, [index in indices],
                    base_name=varname
                )
            end
            # unpack element (n,m) to the correct place in the ouput dict
            for index in indices
                # a diagonal element only appears on the diagonal
                if n==m
                    dict_mat_vars[index][n,n] = mat_nm[index]
                # a lower triangle element apears in the lower triangle
                # and also in the upper triangle (but inverted)
                else
                    dict_mat_vars[index][n,m] =  mat_nm[index]
                    dict_mat_vars[index][m,n] = -mat_nm[index]
                end
            end
        end
    end
    return dict_mat_vars
end


"""
Returns a pair of symmetric and skew-symmetric matrix variables.
"""
function variable_mx_hermitian(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int,
        lower_bound::Dict{T,Array{LB,2}}, upper_bound::Dict{T,Array{UB,2}};
        name::Union{String, Tuple{String,String}}="", prefix="") where {T, LB<:Real, UB<:Real}
    name_real = isa(name, Tuple) ? name[1] : "$(name)re"
    name_imag = isa(name, Tuple) ? name[2] : "$(name)im"
    Mre = variable_mx_symmetric(pm, indices, N, lower_bound, upper_bound;
        prefix=prefix, name=name_real)
    Mim = variable_mx_skewsymmetric(pm, indices, N, lower_bound, upper_bound;
        prefix=prefix, name=name_imag)
    return (Mre,Mim)
end

"""
Alternative method of variable_mx_hermitian,
that does not take any bounds, and creates unbounded matrix variables.
"""
function variable_mx_hermitian(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int;
        name::Union{String, Tuple{String,String}}="", prefix="") where {T, LB<:Real, UB<:Real}
    name_real = isa(name, Tuple) ? name[1] : "$(name)re"
    name_imag = isa(name, Tuple) ? name[2] : "$(name)im"
    Mre = variable_mx_symmetric(pm, indices, N;
        prefix=prefix, name=name_real)
    Mim = variable_mx_skewsymmetric(pm, indices, N;
        prefix=prefix, name=name_imag)
    return (Mre,Mim)
end


"""
Often, lower_bound=-upper_bound, so this method takes only a magnitude bound,
which is then applied symmetrically to the rectangular coordinates.
"""
function variable_mx_hermitian(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int,
        bound::Dict{T,Array{B,2}};
        name::Union{String, Tuple{String,String}}="", prefix="") where {T, B<:Real}
    upper_bound = bound
    lower_bound = Dict([(k,-v) for (k,v) in bound])
    return variable_mx_hermitian(pm, indices, N,
        lower_bound, upper_bound,
        name=name, prefix=prefix)
end


"""
Take a Hermitian matrix of the form W = U*U'
Often, bounds on W are specified through magnitude upper bounds on Umax
√(|W[c,c]|) <= Umax, hence Umax is called bound_diag_mag_sqrt_max
In that case, Wmax = Umax*Umax' and Wmin = -Umax*Umax'
Except for the diagonal which is bounded below by 0, diag(Wmin) = 0
"""
function variable_mx_hermitian_sqrt_bounds(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int,
        bound_diag_mag_sqrt_max::Dict{T,Array{B,1}};
        name::Union{String, Tuple{String,String}}="", prefix="") where {T, B<:Real}
    upper_bound = Dict([(k,v*v') for (k,v) in bound_diag_mag_sqrt_max])
    lower_bound = Dict([(k,-w) for (k,w) in upper_bound])
    for c in 1:N
        for (_,w) in lower_bound
            w[c,c] = 0
        end
    end
    return variable_mx_hermitian(pm, indices, N,
        lower_bound, upper_bound,
        name=name, prefix=prefix)
end


"""
Same as variable_mx_hermitian, but has an additional lower bound on the
magnitude of U, Umin>=0
√(|W[c,c]|) >= Umin, hence Umax is called bound_diag_mag_sqrt_min
"""
function variable_mx_hermitian_sqrt_bounds(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int,
        bound_diag_mag_sqrt_max::Dict{T,Array{B,1}},
        bound_diag_mag_sqrt_min::Dict{T,Array{B,1}};
        name::Union{String, Tuple{String,String}}="", prefix="") where {T, B<:Real}
    upper_bound = Dict([(k,v*v') for (k,v) in bound_diag_mag_sqrt_max])
    lower_bound = Dict([(k,-w) for (k,w) in upper_bound])
    for c in 1:N
        for (id, w) in lower_bound
            w[c,c] = bound_diag_mag_sqrt_min[id][c]^2
            @assert(w[c,c]>=0)
        end
    end
    return variable_mx_hermitian(pm, indices, N,
        lower_bound, upper_bound,
        name=name, prefix=prefix)
end
