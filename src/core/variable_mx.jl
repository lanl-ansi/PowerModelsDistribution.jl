# UTILITY MATRIX VARIABLE FUNCTIONS

# Each function has several methods, allowing to
# + pass upper and lower bounds,
# + a single symmetric bound
# + or no bounds at all.

# General matrix variables of size (N,M) are created by
# + variable_mx_real
# + variable_mx_complex
# An alternative is provided which allows to pass the diagonal elements as
# a keyword argument; this is a use case that appears often for P,Q matrices
# + variable_mx_real_with_diag
# + variable_mx_complex_with_diag

# Symmetric, skewsymmetric and Hermitian matrix variables are always square
# variable_mx_symmetric
# variable_mx_skewsymmetric
# variable_mx_hermitian

# A special function is available for the Hermitian matrix of the form W=UU',
# allowing to pass magnitude bounds on U instead
# variable_mx_hermitian_sqrt_bounds

# Most of the method arguments are strongly typed, to make it clear what is expected.
# These functions only depend on JuMP models, not PowerModels.


# GENERAL MATRICES
#%##########################


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
function variable_mx_real(model::JuMP.Model, indices::Array{T,1}, N::Int, M::Int,
        upper_bound::Dict{T,Array{UB,2}}, lower_bound::Dict{T,Array{LB,2}};
        name="", prefix="") where {T, LB<:Real, UB<:Real}
    # the output is a dictionary of (index, matrix) pairs
    dict_mat_vars_temp = Dict{T, Array{Any, 2}}([(i,fill(missing,N,M)) for i in indices])
    for n in 1:N
        for m in 1:M
            varname = isempty(prefix) ? "$(n)$(m)_$name" : "$(prefix)_$(n)$(m)_$name"
            # create the element (n,m) for all indices
            mat_nm = JuMP.@variable(model, [index in indices],
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
Alternative method of variable_mx_real,
that does not take any bounds, and creates unbounded matrix variables.
"""
function variable_mx_real(model::JuMP.Model, indices::Array{T,1}, N::Int, M::Int;
        name="", prefix="") where T
    # the output is a dictionary of (index, matrix) pairs
    dict_mat_vars_temp = Dict{T, Array{Any, 2}}([(i,fill(missing,N,M)) for i in indices])
    for n in 1:N
        for m in 1:M
            varname = isempty(prefix) ? "$(n)$(m)_$name" : "$(prefix)_$(n)$(m)_$name"
            # create the element (n,m) for all indices
            mat_nm = JuMP.@variable(model, [index in indices],
                base_name=varname
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
variables consists of the elements passed as the diag argument.
The diag argument is a dictionary of (index, 1d-array) pairs.
Useful for power matrices with specified diagonals (constant power wye loads).
If not specified, the diagonal elements are set to missing.
"""
function variable_mx_real_with_diag(model::JuMP.Model, indices::Array{T,1}, N::Int,
        upper_bound::Dict{T,Array{UB,2}}, lower_bound::Dict{T,Array{LB,2}};
        diag::Dict{T,Array{Any,1}}=Dict([(i, fill(missing, N)) for i in indices]),
        name="", prefix="") where {T, LB<:Real, UB<:Real}
    # the output is a dictionary of (index, matrix) pairs
    dict_mat_vars = Dict{T, Array{Any, 2}}([(index, zeros(N,N)) for index in indices])
    for n in 1:N
        for m in 1:N
            if n==m
                for index in indices
                    dict_mat_vars[index][n,n] = diag[index][n]
                end
            else
                varname = isempty(prefix) ? "$(n)$(m)_$name" : "$(prefix)_$(n)$(m)_$name"
                # create the element (n,m) for all indices
                mat_nm = JuMP.@variable(model, [index in indices],
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
Same as variable_mx_real_with diag, but without bounds.
"""
function variable_mx_real_with_diag(model::JuMP.Model, indices::Array{T,1}, N::Int;
        diag::Dict{T,Array{Any,1}}=Dict([(i, fill(missing, N)) for i in indices]),
        name="", prefix="") where T
    # the output is a dictionary of (index, matrix) pairs
    dict_mat_vars = Dict{T, Array{Any, 2}}([(index, zeros(N,N)) for index in indices])
    for n in 1:N
        for m in 1:N
            if n==m
                for index in indices
                    dict_mat_vars[index][n,n] = diag[index][n]
                end
            else
                varname = isempty(prefix) ? "$(n)$(m)_$name" : "$(prefix)_$(n)$(m)_$name"
                # create the element (n,m) for all indices
                mat_nm = JuMP.@variable(model, [index in indices],
                    base_name=varname
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
function variable_mx_complex(model::JuMP.Model, indices::Array{T,1}, N::Int, M::Int,
        upper_bound::Dict{T,Array{UB,2}}, lower_bound::Dict{T,Array{LB,2}};
        name::Union{String, Tuple{String,String}}="", prefix="") where {T, LB<:Real,UB<:Real}
    name_real = isa(name, Tuple) ? name[1] : "$(name)re"
    name_imag = isa(name, Tuple) ? name[2] : "$(name)im"
    Mre = variable_mx_real(model, indices, N, M, lower_bound, upper_bound;
        prefix=prefix, name=name_real)
    Mim = variable_mx_real(model, indices, N, M, lower_bound, upper_bound;
        prefix=prefix, name=name_imag)
    return (Mre,Mim)
end


"""
Same as variable_mx_complex_with_diag, but with symmetric bounds.
"""
function variable_mx_complex(model::JuMP.Model, indices::Array{T,1}, N::Int, M::Int,
        bound::Dict{T,Array{B,2}}; kwargs...) where {T, B<:Real}
    upper_bound = bound
    lower_bound = Dict([(k,-v) for (k,v) in bound])
    return variable_mx_complex(model, indices, N, M, lower_bound, upper_bound; kwargs...)
end


"""
Alternative method of variable_mx_complex,
that does not take any bounds, and creates unbounded matrix variables.
"""
function variable_mx_complex(model::JuMP.Model, indices::Array{T,1}, N::Int, M::Int;
        name::Union{String, Tuple{String,String}}="", prefix="") where T
    name_real = isa(name, Tuple) ? name[1] : "$(name)re"
    name_imag = isa(name, Tuple) ? name[2] : "$(name)im"
    Mre = variable_mx_real(model, indices, N, M;
        prefix=prefix, name=name_real)
    Mim = variable_mx_real(model, indices, N, M;
        prefix=prefix, name=name_imag)
    return (Mre,Mim)
end


"""
Same as variable_mx_complex, but square and the diagonal of the matrix variables
consists of the constants passed as the diag_re and diag_im argument. The diag
argument is a dictionary of (index, 1d-array) pairs.
Useful for power matrices with specified diagonals (constant power wye loads).
"""
function variable_mx_complex_with_diag(model::JuMP.Model, indices::Array{T,1}, N::Int,
        upper_bound::Dict{T,Array{UB,2}}, lower_bound::Dict{T,Array{LB,2}};
        diag_re::Dict{T,Array{Any,1}}=Dict([(i, zeros(N)) for i in indices]),
        diag_im::Dict{T,Array{Any,1}}=Dict([(i, zeros(N)) for i in indices]),
        name::Union{String, Tuple{String,String}}="", prefix="") where {T, LB<:Real,UB<:Real}
    name_real = isa(name, Tuple) ? name[1] : "$(name)re"
    name_imag = isa(name, Tuple) ? name[2] : "$(name)im"
    Mre = variable_mx_real_with_diag(model, indices, N, lower_bound, upper_bound, diag_re;
        prefix=prefix, name=name_real
    )
    Mim = variable_mx_real_with_diag(model, indices, N, lower_bound, upper_bound, diag_im;
        prefix=prefix, name=name_imag)
    return (Mre,Mim)
end


"""
Same as variable_mx_complex_with_diag, but symmetric bounds.
"""
function variable_mx_complex_with_diag(model::JuMP.Model, indices::Array{T,1}, N::Int,
        bound::Dict{T,Array{B,2}}; kwargs...) where {T,B<:Real}
    upper_bound = bound
    lower_bound = Dict([(k,-v) for (k,v) in bound])
    return variable_mx_complex_with_diag(model, indices, N, lower_bound, upper_bound; kwargs...)
end


"""
Alternative method of variable_mx_complex_with_diag,
that does not take any bounds, and creates unbounded matrix variables.
"""
function variable_mx_complex_with_diag(model::JuMP.Model, indices::Array{T,1}, N::Int;
        diag_re::Dict{T,Array{Any,1}}=Dict([(i, zeros(N)) for i in indices]),
        diag_im::Dict{T,Array{Any,1}}=Dict([(i, zeros(N)) for i in indices]),
        name::Union{String, Tuple{String,String}}="", prefix="") where T
    name_real = isa(name, Tuple) ? name[1] : "$(name)re"
    name_imag = isa(name, Tuple) ? name[2] : "$(name)im"
    Mre = variable_mx_real_with_diag(model, indices, N; diag=diag_re,
        prefix=prefix, name=name_real
    )
    Mim = variable_mx_real_with_diag(model, indices, N, diag=diag_im,
        prefix=prefix, name=name_imag)
    return (Mre,Mim)
end


# HERMITIAN MATRICES
#%##########################


"""
Same as variable_mx_real, but adds symmetry structure
"""
function variable_mx_symmetric(model::JuMP.Model, indices::Array{T,1}, N::Int,
        upper_bound::Dict{T,Array{UB,2}}, lower_bound::Dict{T,Array{LB,2}};
        name="", prefix="") where {T, LB<:Real, UB<:Real}
    # the output is a dictionary of (index, matrix) pairs
    # the data type of the matrix has to be JuMP.GenericAffExpr, because it
    # will also contain inverted variables (upper triangle) and potentially
    # constants (zero) on the diagonal
    dict_mat_vars = Dict{T, Array{Any, 2}}([(index, zeros(N,N)) for index in indices])
    for n in 1:N
        for m in 1:n
            varname = isempty(prefix) ? "$(n)$(m)_$name" : "$(prefix)_$(n)$(m)_$name"
            # create the lower triangle element (n,m) for all indices
            mat_nm = JuMP.@variable(model, [index in indices],
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
    # recast as Array of VariableRef, useful for JuMP inspection
    return Dict{T, Array{JuMP.VariableRef, 2}}([x for x in dict_mat_vars])
end


"""
Alternative method of variable_mx_symmetric,
that does not take any bounds, and creates unbounded matrix variables.
"""
function variable_mx_symmetric(model::JuMP.Model, indices::Array{T,1}, N::Int;
        name="", prefix="", lb_diag_zero=false) where T
    # the output is a dictionary of (index, matrix) pairs
    # the data type of the matrix has to be JuMP.GenericAffExpr, because it
    # will also contain inverted variables (upper triangle) and potentially
    # constants (zero) on the diagonal
    dict_mat_vars = Dict{T, Array{Any, 2}}([(index, fill(missing, N, N)) for index in indices])
    for n in 1:N
        for m in 1:n
            varname = isempty(prefix) ? "$(n)$(m)_$name" : "$(prefix)_$(n)$(m)_$name"
            # create the lower triangle element (n,m) for all indices
            # even when not bounded, still bound diagonal below at zero if flagged
            if m==n && lb_diag_zero
                mat_nm = JuMP.@variable(model, [index in indices],
                    base_name=varname,
                    lower_bound=0
                )
            else
                mat_nm = JuMP.@variable(model, [index in indices],
                    base_name=varname
                )
            end
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
    # recast as Array of VariableRef, useful for JuMP inspection
    return Dict{T, Array{JuMP.VariableRef, 2}}([x for x in dict_mat_vars])
end


"""
Same as variable_mx_real, but adds skew-symmetry structure.
"""
function variable_mx_skewsymmetric(model::JuMP.Model, indices::Array{T,1}, N::Int,
        upper_bound::Dict{T,Array{UB,2}}, lower_bound::Dict{T,Array{LB,2}};
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
                mat_nm = JuMP.@variable(model, [index in indices],
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
function variable_mx_skewsymmetric(model::JuMP.Model, indices::Array{T,1}, N::Int;
        name="", prefix="") where T
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
                mat_nm = JuMP.@variable(model, [index in indices],
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
function variable_mx_hermitian(model::JuMP.Model, indices::Array{T,1}, N::Int,
        lower_bound::Dict{T,Array{LB,2}}, upper_bound::Dict{T,Array{UB,2}};
        name::Union{String, Tuple{String,String}}="", prefix="") where {T, LB<:Real, UB<:Real}
    name_real = isa(name, Tuple) ? name[1] : "$(name)re"
    name_imag = isa(name, Tuple) ? name[2] : "$(name)im"
    Mre = variable_mx_symmetric(model, indices, N, lower_bound, upper_bound;
        prefix=prefix, name=name_real)
    Mim = variable_mx_skewsymmetric(model, indices, N, lower_bound, upper_bound;
        prefix=prefix, name=name_imag)
    return (Mre,Mim)
end


"""
Alternative method of variable_mx_hermitian,
that does not take any bounds, and creates unbounded matrix variables.
Note that diagonal is still bounded below at zero if flagged.
"""
function variable_mx_hermitian(model::JuMP.Model, indices::Array{T,1}, N::Int;
        name::Union{String, Tuple{String,String}}="", prefix="", lb_diag_zero=false) where {T, LB<:Real, UB<:Real}
    name_real = isa(name, Tuple) ? name[1] : "$(name)re"
    name_imag = isa(name, Tuple) ? name[2] : "$(name)im"
    Mre = variable_mx_symmetric(model, indices, N;
        prefix=prefix, name=name_real, lb_diag_zero=lb_diag_zero)
    Mim = variable_mx_skewsymmetric(model, indices, N;
        prefix=prefix, name=name_imag)
    return (Mre,Mim)
end


"""
Often, lower_bound=-upper_bound, so this method takes only a magnitude bound,
which is then applied symmetrically to the rectangular coordinates.
"""
function variable_mx_hermitian(model::JuMP.Model, indices::Array{T,1}, N::Int,
        bound::Dict{T,Array{B,2}}; kwargs...) where {T, B<:Real}
    upper_bound = bound
    lower_bound = Dict([(k,-v) for (k,v) in bound])
    return variable_mx_hermitian(model, indices, N, upper_bound, lower_bound; kwargs...)
end


"""
Take a Hermitian matrix of the form W = U*U'
Often, bounds on W are specified through magnitude upper bounds on Umax
√(|W[c,c]|) <= Umax, hence Umax is called bound_diag_mag_sqrt_max
In that case, Wmax = Umax*Umax' and Wmin = -Umax*Umax'
Except for the diagonal which is bounded below by 0, diag(Wmin) = 0
"""
function variable_mx_hermitian_sqrt_bounds(model::JuMP.Model, indices::Array{T,1}, N::Int,
        bound_diag_mag_sqrt_max::Dict{T,Array{B,1}}; kwargs...) where {T, B<:Real}
    upper_bound = Dict([(k,v*v') for (k,v) in bound_diag_mag_sqrt_max])
    lower_bound = Dict([(k,-w) for (k,w) in upper_bound])
    for c in 1:N
        for (_,w) in lower_bound
            w[c,c] = 0
        end
    end
    return variable_mx_hermitian(model, indices, N, upper_bound, lower_bound; kwargs...)
end


"""
Same as variable_mx_hermitian, but has an additional lower bound on the
magnitude of U, Umin>=0
√(|W[c,c]|) >= Umin, hence Umax is called bound_diag_mag_sqrt_min
"""
function variable_mx_hermitian_sqrt_bounds(model::JuMP.Model, indices::Array{T,1}, N::Int,
        bound_diag_mag_sqrt_max::Dict{T,Array{B,1}},
        bound_diag_mag_sqrt_min::Dict{T,Array{B,1}}; kwargs...) where {T, B<:Real}
    upper_bound = Dict([(k,v*v') for (k,v) in bound_diag_mag_sqrt_max])
    lower_bound = Dict([(k,-w) for (k,w) in upper_bound])
    for c in 1:N
        for (id, w) in lower_bound
            w[c,c] = bound_diag_mag_sqrt_min[id][c]^2
            @assert(w[c,c]>=0)
        end
    end
    return variable_mx_hermitian(model, indices, N, upper_bound, lower_bound; kwargs...)
end
