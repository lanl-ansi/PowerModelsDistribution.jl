# UTILITY MATRIX VARIABLE FUNCTIONS

"""
    _make_matrix_variable_element(
        model::JuMP.Model,
        indices::Array{T,1},
        n::Int,
        m::Int;
        upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        varname::String=""
    ) where T

Sometimes we want to bound only a subset of the elements of a matrix variable.
For example, an unbounded Hermitian variable usually still has a lower bound
of zero on the real diagonal elements. When there is a mix of bounded and unbounded elements,
the unboundedness is encoded as 'Inf' and '-Inf' in the bound parameters.
This cannot be passed directlty to JuMP, because it would lead to an error in Mosek for example.
Instead, this method checks whether all bounds for an element (n,m) are Inf, and if so,
does not pass a bound to JuMP.
"""
function _make_matrix_variable_element(model::JuMP.Model, indices::Array{T,1}, n::Int, m::Int;
    upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
    lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
    varname::String="") where T
    nm_has_lb = !(ismissing(lower_bound) || all([lower_bound[index][n,m]==-Inf || isnan(lower_bound[index][n,m]) for index in indices]))
    nm_has_ub = !(ismissing(upper_bound) || all([upper_bound[index][n,m]== Inf || isnan(upper_bound[index][n,m]) for index in indices]))

    if !nm_has_ub && !nm_has_lb
        mat_nm = JuMP.@variable(model, [index in indices],
            base_name=varname
        )
    elseif nm_has_ub && !nm_has_lb
        mat_nm = JuMP.@variable(model, [index in indices],
            base_name=varname,
            upper_bound=upper_bound[index][n,m]
        )
    elseif !nm_has_ub && nm_has_lb
        mat_nm = JuMP.@variable(model, [index in indices],
            base_name=varname,
            lower_bound=lower_bound[index][n,m]
        )
    else #if nm_has_ub && nm_has_lb
        mat_nm = JuMP.@variable(model, [index in indices],
            base_name=varname,
            lower_bound=lower_bound[index][n,m],
            upper_bound=upper_bound[index][n,m]
        )
    end

    return mat_nm
end

# GENERAL MATRICES
##################

"""
    variable_mx_real(
        model::JuMP.Model,
        indices::Array{T,1},
        N::Dict{T,Vector{Int}},
        M::Dict{T,Vector{Int}};
        upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        name::String="",
        prefix::String=""
    ) where T

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
function variable_mx_real(model::JuMP.Model, indices::Array{T,1}, N::Dict{T,Vector{Int}}, M::Dict{T,Vector{Int}};
    upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
    lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
    name::String="",
    prefix::String="") where T
    # the output is a dictionary of (index, matrix) pairs
    dict_mat_vars_temp = Dict{T, Array{Any, 2}}([(i,fill(missing,length(N[i]),length(M[i]))) for i in indices])
    for index in indices
        for (n,t) in enumerate(N[index])
            for (m,u) in enumerate(M[index])
                varname = isempty(prefix) ? "$(name)_$(t)$(u)" : "$(prefix)_$(name)_$(t)$(u)"
                # create the element (n,m) for all indices
                mat_nm = _make_matrix_variable_element(model, [index], n, m;
                    upper_bound=upper_bound, lower_bound=lower_bound, varname=varname)
                # unpack element (n,m) to the correct place in the ouput dict
                dict_mat_vars_temp[index][n,m] = mat_nm[index]
            end
        end
    end
    dict_mat_vars = Dict{T, Array{JuMP.VariableRef, 2}}([(id, arr) for (id, arr) in dict_mat_vars_temp])

    return dict_mat_vars
end


"""
    variable_mx_real_with_diag(
        model::JuMP.Model,
        indices::Array{T,1},
        N::Dict{T,Vector{Int}};
        upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        diag::Dict{T,<:Vector{<:Any}}=Dict([(i, fill(0, length(N[i]))) for i in indices]),
        name::String="",
        prefix::String=""
    ) where T

Same as [`variable_mx_real`](@ref variable_mx_real), but has to be square and the diagonal of the matrix
variables consists of the elements passed as the diag argument.
The diag argument is a dictionary of (index, 1d-array) pairs.
Useful for power matrices with specified diagonals (constant power wye loads).
If not specified, the diagonal elements are set to zero.
"""
function variable_mx_real_with_diag(model::JuMP.Model, indices::Array{T,1}, N::Dict{T,Vector{Int}};
    upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
    lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
    diag::Dict{T,<:Vector{<:Any}}=Dict([(i, fill(0, length(N[i]))) for i in indices]),
    name::String="",
    prefix::String="") where T
    # the output is a dictionary of (index, matrix) pairs
    dict_mat_vars = Dict{T,Array{JuMP.GenericAffExpr{Float64,JuMP.VariableRef},2}}([(index, zeros(length(N[index]),length(N[index]))) for index in indices])
    for index in indices
        for (n,t) in enumerate(N[index])
            for (m,u) in enumerate(N[index])
                if n==m
                    dict_mat_vars[index][n,n] = diag[index][n]
                else
                    varname = isempty(prefix) ? "$(name)_$(t)$(u)" : "$(prefix)_$(name)_$(t)$(u)"
                    # create the element (n,m) for all indices
                    mat_nm = _make_matrix_variable_element(model, [index], n, m;
                        upper_bound=upper_bound, lower_bound=lower_bound, varname=varname)
                    # unpack element (n,m) to the correct place in the ouput dict

                    dict_mat_vars[index][n,m] = mat_nm[index]
                end
            end
        end
    end
    return dict_mat_vars
end


"""
    variable_mx_complex(
        model::JuMP.Model,
        indices::Array{T,1},
        N::Dict{T,Vector{Int}},
        M::Dict{T,Vector{Int}};
        upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        symm_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        name::Union{String, Tuple{String,String}}="",
        prefix::String=""
    )::Tuple where T

Shorthand to create two real matrix variables, where the first is the real part
and the second the imaginary part.

If the name argument is a String, it will be suffixed with 're' and  'im'.
It is possible to  specify the names of the real and imaginary part directly
as a Tuple as well (to achieve P and Q instead of Sre and Sim for example).
"""
function variable_mx_complex(model::JuMP.Model, indices::Array{T,1}, N::Dict{T,Vector{Int}}, M::Dict{T,Vector{Int}};
    upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
    lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
    symm_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
    name::Union{String, Tuple{String,String}}="",
    prefix::String="")::Tuple where T
    if !ismissing(symm_bound)
        @assert(ismissing(upper_bound) && ismissing(lower_bound), "When a symmetric bound is specified, no lower or upper bound can be specified.")
        upper_bound = symm_bound
        lower_bound = Dict{T, Array{Real,2}}([(k,-v) for (k,v) in symm_bound])
    end

    name_real = isa(name, Tuple) ? name[1] : "$(name)r"
    name_imag = isa(name, Tuple) ? name[2] : "$(name)i"
    Mre = variable_mx_real(model, indices, N, M; upper_bound=upper_bound, lower_bound=lower_bound, prefix=prefix, name=name_real)
    Mim = variable_mx_real(model, indices, N, M; upper_bound=upper_bound, lower_bound=lower_bound, prefix=prefix, name=name_imag)

    return (Mre,Mim)
end


"""
    variable_mx_complex_with_diag(
        model::JuMP.Model,
        indices::Array{T,1},
        N::Dict{T,Vector{Int}};
        upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        symm_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        diag_re::Dict{T,<:Vector{<:Any}}=Dict([(i, zeros(length(N[i]))) for i in indices]),
        diag_im::Dict{T,<:Vector{<:Any}}=Dict([(i, zeros(length(N[i]))) for i in indices]),
        name::Union{String, Tuple{String,String}}="",
        prefix::String=""
    )::Tuple where T

Same as [`variable_mx_complex`](@ref variable_mx_complex), but square and the diagonal of the matrix variables
consists of the constants passed as the diag_re and diag_im argument. The diag
argument is a dictionary of (index, 1d-array) pairs.

Useful for power matrices with specified diagonals (constant power wye loads).
"""
function variable_mx_complex_with_diag(model::JuMP.Model, indices::Array{T,1}, N::Dict{T,Vector{Int}};
        upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing, lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        symm_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        diag_re::Dict{T,<:Vector{<:Any}}=Dict([(i, zeros(length(N[i]))) for i in indices]),
        diag_im::Dict{T,<:Vector{<:Any}}=Dict([(i, zeros(length(N[i]))) for i in indices]),
        name::Union{String, Tuple{String,String}}="",
        prefix::String="")::Tuple where T

    if !ismissing(symm_bound)
        @assert(ismissing(upper_bound) && ismissing(lower_bound), "When a symmetric bound is specified, no lower or upper bound can be specified.")
        upper_bound = symm_bound
        lower_bound = Dict{T, Array{Real,2}}([(k,-v) for (k,v) in symm_bound])
    end

    name_real = isa(name, Tuple) ? name[1] : "$(name)r"
    name_imag = isa(name, Tuple) ? name[2] : "$(name)i"
    Mre = variable_mx_real_with_diag(model, indices, N; upper_bound=upper_bound, lower_bound=lower_bound, diag=diag_re, prefix=prefix, name=name_real)
    Mim = variable_mx_real_with_diag(model, indices, N; upper_bound=upper_bound, lower_bound=lower_bound, diag=diag_im, prefix=prefix, name=name_imag)

    return (Mre,Mim)
end

# HERMITIAN MATRICES
####################

"""
    variable_mx_real_symmetric(
        model::JuMP.Model,
        indices::Array{T,1},
        N::Dict{T,Vector{Int}};
        upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        name::String="",
        prefix::String=""
    )::Dict where T

Same as [`variable_mx_real`](@ref variable_mx_real), but adds symmetry structure
"""
function variable_mx_real_symmetric(model::JuMP.Model, indices::Array{T,1}, N::Dict{T,Vector{Int}};
        upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        name::String="",
        prefix::String="")::Dict where T
    # the output is a dictionary of (index, matrix) pairs
    # the data type of the matrix has to be JuMP.GenericAffExpr, because it
    # will also contain inverted variables (upper triangle) and potentially
    # constants (zero) on the diagonal
    dict_mat_vars = Dict{T, Array{Any, 2}}([(index, zeros(length(N[index]),length(N[index]))) for index in indices])
    for index in indices
        for (n,t) in enumerate(N[index])
            for (m,u) in enumerate(N[index][1:n])
                varname = isempty(prefix) ? "$(name)_$(t)$(u)" : "$(prefix)_$(name)_$(t)$(u)"
                # create the lower triangle element (n,m) for all indices
                mat_nm = _make_matrix_variable_element(model, [index], n, m; upper_bound=upper_bound, lower_bound=lower_bound, varname=varname)
                # unpack element (n,m) to the correct place in the ouput dict

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
    # create explicit Symmetric matrices, supported by JuMP
    # cannot do it directly through JuMP because of name generation
    # out = Dict{T, Array{JuMP.VariableRef, 2}}([x for x in dict_mat_vars])
    out = convert(Dict{T,Array{JuMP.VariableRef,2}}, dict_mat_vars)
    return Dict([(k, LinearAlgebra.Symmetric(v)) for (k, v) in out])
end


"""
    variable_mx_real_skewsymmetric(
        model::JuMP.Model,
        indices::Array{T,1},
        N::Dict{T,Vector{Int}};
        upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        set_diag_to_zero::Bool=true,
        name::String="",
        prefix::String=""
    )::Dict where T

Same as [`variable_mx_real`](@ref variable_mx_real), but adds skew-symmetry structure.
"""
function variable_mx_real_skewsymmetric(model::JuMP.Model, indices::Array{T,1}, N::Dict{T,Vector{Int}};
        upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        set_diag_to_zero::Bool=true,
        name::String="",
        prefix::String="")::Dict where T
    # the output is a dictionary of (index, matrix) pairs
    # the data type of the matrix has to be JuMP.GenericAffExpr, because it
    # will also contain inverted variables (upper triangle) and potentially
    # constants (zero) on the diagonal
    dict_mat_vars = Dict{T, Array{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, 2}}([(index, zeros(length(N[index]),length(N[index]))) for index in indices])
    for index in indices
        for (n,t) in enumerate(N[index])
            for (m,u) in enumerate(N[index][1:n])
                varname = isempty(prefix) ? "$(name)_$(t)$(u)" : "$(prefix)_$(name)_$(t)$(u)"
                # create the lower triangle element (n,m) for all indices
                # if diagonal element (n,m) is zero
                if m==n && set_diag_to_zero
                    mat_nm = Dict{T, Float64}(index => 0.0)
                # otherwise, create JuMP variables
                else
                    mat_nm = _make_matrix_variable_element(model, [index], n, m; upper_bound=upper_bound, lower_bound=lower_bound, varname=varname)
                end
                # unpack element (n,m) to the correct place in the ouput dict
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
    variable_mx_hermitian(
        model::JuMP.Model,
        indices::Array{T,1},
        N::Dict{T,Vector{Int}};
        upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        symm_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        sqrt_upper_bound::Union{Missing, Dict{T,<:Vector{<:Real}}}=missing,
        sqrt_lower_bound::Union{Missing, Dict{T,<:Vector{<:Real}}}=missing,
        set_lower_bound_diag_to_zero::Bool=false,
        imag_set_diag_to_zero::Bool=true,
        name::Union{String,Tuple{String,String}}="",
        prefix::String=""
    )::Tuple where T

Returns a pair of symmetric and skew-symmetric matrix variables.
"""
function variable_mx_hermitian(model::JuMP.Model, indices::Array{T,1}, N::Dict{T,Vector{Int}};
        upper_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        lower_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        symm_bound::Union{Missing, Dict{T,<:Matrix{<:Real}}}=missing,
        sqrt_upper_bound::Union{Missing, Dict{T,<:Vector{<:Real}}}=missing,
        sqrt_lower_bound::Union{Missing, Dict{T,<:Vector{<:Real}}}=missing,
        set_lower_bound_diag_to_zero::Bool=false,
        imag_set_diag_to_zero::Bool=true,
        name::Union{String,Tuple{String,String}}="",
        prefix::String="")::Tuple where T

    # translate special bound options and ensure not in conflict
    if !ismissing(symm_bound)
        @assert(ismissing(upper_bound) && ismissing(lower_bound), "When a symmetric bound is specified, no lower or upper bound can be specified.")
        upper_bound = symm_bound
        lower_bound = Dict([(k,-v) for (k,v) in symm_bound])
    end

    if !ismissing(sqrt_upper_bound)
        @assert(ismissing(upper_bound) && ismissing(lower_bound) && ismissing(symm_bound), "When a square root bound is specified, no lower, upper or symmetric bound can be specified.")
        upper_bound = Dict([(k,v*v') for (k,v) in sqrt_upper_bound])
        lower_bound = Dict([(k,-w) for (k,w) in upper_bound])

        if !ismissing(sqrt_lower_bound)
            for (id, w) in lower_bound
                for (idx,t) in enumerate(N[id])
                    w[idx,idx] = sqrt_lower_bound[id][idx]^2
                    #@assert(w[c,c]>=0)
                end
            end
        end
    end

    if set_lower_bound_diag_to_zero
        if !ismissing(lower_bound)
            for (id, M) in lower_bound
                for (idx,t) in enumerate(N[id])
                    M[idx, idx] = 0.0
                end
            end
        end
    end

    name_real = isa(name, Tuple) ? name[1] : "$(name)r"
    name_imag = isa(name, Tuple) ? name[2] : "$(name)i"
    Mre = variable_mx_real_symmetric(model, indices, N; upper_bound=upper_bound, lower_bound=lower_bound, prefix=prefix, name=name_real)
    Mim = variable_mx_real_skewsymmetric(model, indices, N; upper_bound=upper_bound, lower_bound=lower_bound, set_diag_to_zero=imag_set_diag_to_zero, prefix=prefix, name=name_imag)

    return (Mre,Mim)
end
