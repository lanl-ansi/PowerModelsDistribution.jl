"""
TODO
- Should unbounded be supported as well by variable_mx? If so, redesign signatures
"""
import LinearAlgebra: diag, diagm


function variable_tp_branch_current_v2(pm::PMs.GenericPowerModel{T}; kwargs...) where T <: AbstractUBFForm
    variable_tp_branch_series_current_prod_hermitian_v2(pm; kwargs...)
end


function variable_tp_voltage_v2(pm::PMs.GenericPowerModel{T}; kwargs...) where T <: AbstractUBFForm
    variable_tp_voltage_prod_hermitian_v2(pm; kwargs...)
end


function variable_tp_voltage_prod_hermitian_v2(pm::PMs.GenericPowerModel{T}; n_cond::Int=3, nw::Int=pm.cnw, bounded = true) where T <: AbstractUBFForm
    bus_ids = collect(PMs.ids(pm, nw, :bus))
    # calculate bounds
    lower_bound = Dict{eltype(bus_ids), Array{Real,2}}()
    upper_bound = Dict{eltype(bus_ids), Array{Real,2}}()
    for id in bus_ids
        vmax = PMs.ref(pm, nw, :bus, id)["vmax"].values
        vmin = PMs.ref(pm, nw, :bus, id)["vmin"].values
        lower_bound[id] = -vmax*vmax'
        for c in 1:n_cond
            lower_bound[id][c,c] = vmin[c]^2
        end
        upper_bound[id] =  vmax*vmax'
    end
    # create  matrix variables
    # TODO is it okay to always be bounded?
    (Wre,Wim) = variable_mx_hermitian(pm, bus_ids, n_cond, lower_bound, upper_bound; name="W")
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

    # calculate max series current for each branch
    cmax = Dict{Int, Array{Real,1}}()
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
    # create  L bound
    branch_ids = collect(keys(branches))
    bound = Dict{eltype(branch_ids), Array{Real,2}}([(id, cmax[id]*cmax[id]') for id in branch_ids])
    # create matrix variables
    (Lre,Lim) = variable_mx_hermitian(pm, branch_ids, n_cond, bound; name="L")
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
    (P,Q) = variable_mx_complex(pm, branch_arcs, n_cond, bound; name=("P", "Q"))
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


# UTILITY MATRIX VARIABLE FUNCTIONS
"""
This function creates a set of real matrix variables of size NxN,
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
function variable_mx_real(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int,
        lower_bound::Dict{T,Array{LB,2}}, upper_bound::Dict{T,Array{UB,2}};
        name="", prefix="") where {T, LB<:Real, UB<:Real}
    # the output is a dictionary of (index, matrix) pairs
    dict_mat_vars = Dict{T, Array{Union{JuMP.VariableRef,Nothing}, 2}}([(i,fill(nothing,N,N)) for i in indices])
    for n in 1:N
        for m in 1:N
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
    return dict_mat_vars
end


"""
Shorthand to create two real matrix variables, where the first is the real part
and the second the imaginary part.
If the name argument is a String, it will be suffixed with 're' and  'im'.
It is possible to  specify the names of the real and imaginary part directly
as a Tuple as well (to achieve P and Q instead of Sre and Sim for example).
"""
function variable_mx_complex(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int,
        lower_bound::Dict{T,Array{LB,2}}, upper_bound::Dict{T,Array{UB,2}};
        name::Union{String, Tuple{String,String}}="", prefix="") where {T, LB<:Real,UB<:Real}
    name_real = isa(name, Tuple) ? name[1] : "$(name)re"
    name_imag = isa(name, Tuple) ? name[2] : "$(name)im"
    Mre = variable_mx_real(pm, indices, N,
        lower_bound, upper_bound, prefix=prefix,
        name=name_real
    )
    Mim = variable_mx_real(pm, indices, N, lower_bound, upper_bound;
        prefix=prefix, name=name_imag)
    return (Mre,Mim)
end


"""
Often, lower_bound=-upper_bound, so this method takes only a magnitude bound,
which is then applied symmetrically to the rectangular coordinates.
"""
function variable_mx_complex(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int,
        bound::Dict{T,Array{B,2}};
        name::Union{String, Tuple{String,String}}="", prefix="") where {T, B<:Real}
    upper_bound = bound
    lower_bound = Dict([(k,-v) for (k,v) in bound])
    return variable_mx_complex(pm, indices, N,
        lower_bound, upper_bound,
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
Same as variable_mx_real, but adds skew-symmetry structure.
Also, the diagonal can be fixed to zero with the flag diag_is_zero. This is
useful for the imaginary part of Hermitian matrices.
"""
function variable_mx_skewsymmetric(pm::PMs.GenericPowerModel, indices::Array{T,1}, N::Int,
        lower_bound::Dict{T,Array{LB,2}}, upper_bound::Dict{T,Array{UB,2}};
        diag_is_zero=false, name="", prefix="") where {T, LB<:Real, UB<:Real}
    # the output is a dictionary of (index, matrix) pairs
    # the data type of the matrix has to be JuMP.GenericAffExpr, because it
    # will also contain inverted variables (upper triangle) and potentially
    # constants (zero) on the diagonal
    dict_mat_vars = Dict{T, Array{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, 2}}([(index, zeros(N,N)) for index in indices])
    for n in 1:N
        for m in 1:n
            varname = isempty(prefix) ? "$(n)$(m)_$name" : "$(prefix)_$(n)$(m)_$name"
            # create the lower triangle element (n,m) for all indices
            # if diagonal element and diag_is_zero is flagged, (n,m) is zero
            if m==n && diag_is_zero
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
        diag_is_zero=true, prefix=prefix, name=name_imag)
    return (Mre,Mim)
end
