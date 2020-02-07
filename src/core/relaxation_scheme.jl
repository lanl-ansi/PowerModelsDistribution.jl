function partition_matrix_with_scalar(Mre, Mim)
    n = size(Mre, 1)
    alphare = Mre[1,1]
    Are = Mre[2:n, 2:n]
    Aim = Mim[2:n, 2:n]
    are = Mre[2:n, 1]
    aim = Mim[2:n, 1]
    return Are, Aim, are, aim, alphare
end

function relaxation_psd_to_soc_complex_kim_kojima_3x3(model, Mre, Mim)
    @assert size(Mre) == size(Mim) == (3,3)
    cre = [1;1]
    cim = [0;0]

    pm1 = [1; 2; 3]
    Are, Aim, are, aim, alphare = partition_matrix_with_scalar(Mre[pm1,pm1], Mim[pm1,pm1])
    relaxation_psd_to_soc_complex_kim_kojima(model, Are, Aim, are, aim, alphare, cre, cim)

    pm2 = [2; 3; 1]
    Are, Aim, are, aim, alphare = partition_matrix_with_scalar(Mre[pm2,pm2], Mim[pm2,pm2])
    relaxation_psd_to_soc_complex_kim_kojima(model, Are, Aim, are, aim, alphare, cre, cim)

    pm3 = [3; 1; 2]
    Are, Aim, are, aim, alphare = partition_matrix_with_scalar(Mre[pm3,pm3], Mim[pm3,pm3])
    relaxation_psd_to_soc_complex_kim_kojima(model, Are, Aim, are, aim, alphare, cre, cim)

    cre = [1;0]
    cim = [0;1]
    pm1 = [1; 2; 3]
    Are, Aim, are, aim, alphare = partition_matrix_with_scalar(Mre[pm1,pm1], Mim[pm1,pm1])
    relaxation_psd_to_soc_complex_kim_kojima(model, Are, Aim, are, aim, alphare, cre, cim)

    pm2 = [2; 3; 1]
    Are, Aim, are, aim, alphare = partition_matrix_with_scalar(Mre[pm2,pm2], Mim[pm2,pm2])
    relaxation_psd_to_soc_complex_kim_kojima(model, Are, Aim, are, aim, alphare, cre, cim)

    pm3 = [3; 1; 2]
    Are, Aim, are, aim, alphare = partition_matrix_with_scalar(Mre[pm3,pm3], Mim[pm3,pm3])
    relaxation_psd_to_soc_complex_kim_kojima(model, Are, Aim, are, aim, alphare, cre, cim)


end


function relaxation_psd_to_soc_complex_kim_kojima_3x3_conic(model, Mre, Mim)
    @assert size(Mre) == size(Mim) == (3,3)
    cre = [1;1]
    cim = [0;0]

    pm1 = [1; 2; 3]
    Are, Aim, are, aim, alphare = partition_matrix_with_scalar(Mre[pm1,pm1], Mim[pm1,pm1])
    relaxation_psd_to_soc_complex_kim_kojima_conic(model, Are, Aim, are, aim, alphare, cre, cim)

    pm2 = [2; 3; 1]
    Are, Aim, are, aim, alphare = partition_matrix_with_scalar(Mre[pm2,pm2], Mim[pm2,pm2])
    relaxation_psd_to_soc_complex_kim_kojima_conic(model, Are, Aim, are, aim, alphare, cre, cim)

    pm3 = [3; 1; 2]
    Are, Aim, are, aim, alphare = partition_matrix_with_scalar(Mre[pm3,pm3], Mim[pm3,pm3])
    relaxation_psd_to_soc_complex_kim_kojima_conic(model, Are, Aim, are, aim, alphare, cre, cim)

    cre = [1;0]
    cim = [0;1]

    pm1 = [1; 2; 3]
    Are, Aim, are, aim, alphare = partition_matrix_with_scalar(Mre[pm1,pm1], Mim[pm1,pm1])
    relaxation_psd_to_soc_complex_kim_kojima_conic(model, Are, Aim, are, aim, alphare, cre, cim)

    pm2 = [2; 3; 1]
    Are, Aim, are, aim, alphare = partition_matrix_with_scalar(Mre[pm2,pm2], Mim[pm2,pm2])
    relaxation_psd_to_soc_complex_kim_kojima_conic(model, Are, Aim, are, aim, alphare, cre, cim)

    pm3 = [3; 1; 2]
    Are, Aim, are, aim, alphare = partition_matrix_with_scalar(Mre[pm3,pm3], Mim[pm3,pm3])
    relaxation_psd_to_soc_complex_kim_kojima_conic(model, Are, Aim, are, aim, alphare, cre, cim)

end



"""
SDP to SOC relaxation of type 1, applied to complex-value matrix, extended from:
```
@article{Kim2003,
author = {Kim, S and Kojima, M and Yamashita, M},
title = {{Second order cone programming relaxation of a positive semidefinite constraint}},
doi = {10.1080/1055678031000148696},
journal = {Optimization Methods and Software},
number = {5},
pages = {535--541},
volume = {18},
year = {2003}
}
```
"""
function relaxation_psd_to_soc_complex_kim_kojima(model, Are, Aim, are, aim, alphare, cre, cim; tol=1e-8)
    @assert size(Are) == size(Aim)
    @assert size(are) == size(aim)

    @assert size(Are)[1] == size(are)[1]

    cre[abs.(cre).<=tol] .=0
    cim[abs.(cim).<=tol] .=0

    c = cre+im*cim
    C = c*c'
    Cre = real(C)
    Cim = imag(C)
    Cre[abs.(Cre).<=tol] .=0
    Cim[abs.(Cim).<=tol] .=0

    @assert size(Cre) == size(Are)

    rhs_1 = alphare
    rhs_2 = sum(Cre.*Are) + sum(Cim.*Aim)

    lhs_re = cre'* are + cim'* aim
    lhs_im = cre'* aim - cim'* are


    JuMP.@constraint(model, rhs_2 >= 0)
    JuMP.@constraint(model, lhs_re'*lhs_re + lhs_im'*lhs_im <= rhs_1*rhs_2)

end


"""
SDP to SOC relaxation of type 1, applied to complex-value matrix, extended from:
```
@article{Kim2003,
author = {Kim, S and Kojima, M and Yamashita, M},
title = {{Second order cone programming relaxation of a positive semidefinite constraint}},
doi = {10.1080/1055678031000148696},
journal = {Optimization Methods and Software},
number = {5},
pages = {535--541},
volume = {18},
year = {2003}
}
```
"""
function relaxation_psd_to_soc_complex_kim_kojima_conic(model, Are, Aim, are, aim, alphare, cre, cim; tol=1e-8)
    @assert size(Are) == size(Aim)
    @assert size(are) == size(aim)

    @assert size(Are)[1] == size(are)[1]

    cre[abs.(cre).<=tol] .=0
    cim[abs.(cim).<=tol] .=0

    c = cre+im*cim
    C = c*c'
    Cre = real(C)
    Cim = imag(C)
    Cre[abs.(Cre).<=tol] .=0
    Cim[abs.(Cim).<=tol] .=0

    @assert size(Cre) == size(Are)

    rhs_1 = alphare
    rhs_2 = sum(Cre.*Are) + sum(Cim.*Aim)

    lhs_re = cre'* are + cim'* aim
    lhs_im = cre'* aim - cim'* are

    JuMP.@constraint(model, rhs_2 >= 0)
    JuMP.@constraint(model,     [rhs_1+rhs_2;
                                 rhs_1-rhs_2;
                                 2*lhs_re;
                                 2*lhs_im] in JuMP.SecondOrderCone())

end

"""
SDP to SOC relaxation of type 2, applied to real-value matrix,  as described in:
```
@article{Kim2003,
author = {Kim, S and Kojima, M and Yamashita, M},
title = {{Second order cone programming relaxation of a positive semidefinite constraint}},
doi = {10.1080/1055678031000148696},
journal = {Optimization Methods and Software},
number = {5},
pages = {535--541},
volume = {18},
year = {2003}
}
```
"""
function relaxation_psd_to_soc_real(m, mx)
    @assert size(mx,1) == size(mx,2)
    n_elements = size(mx,1)
    for i in 1:n_elements-1
        for j in i+1:n_elements
            JuMP.@constraint(m, mx[i,j]^2 <= mx[i,i]*mx[j,j])
        end
    end
end


"""
SDP to SOC relaxation of type 2, applied to complex-value matrix,  as described in:
```
@article{Kim2003,
author = {Kim, S and Kojima, M and Yamashita, M},
title = {{Second order cone programming relaxation of a positive semidefinite constraint}},
doi = {10.1080/1055678031000148696},
journal = {Optimization Methods and Software},
number = {5},
pages = {535--541},
volume = {18},
year = {2003}
}
```
"""
function relaxation_psd_to_soc_complex(m, mxreal, mximag)
    @assert size(mxreal) == size(mximag)
    n_elements = size(mxreal,1)
    for i in 1:n_elements-1
        for j in i+1:n_elements
            JuMP.@constraint(m, mxreal[i,j]^2 + mximag[i,j]^2 <= mxreal[i,i]*mxreal[j,j])
        end
    end
end


"""
SDP to SOC relaxation of type 2, applied to real-value matrix,  as described in:
```
@article{Kim2003,
author = {Kim, S and Kojima, M and Yamashita, M},
title = {{Second order cone programming relaxation of a positive semidefinite constraint}},
doi = {10.1080/1055678031000148696},
journal = {Optimization Methods and Software},
number = {5},
pages = {535--541},
volume = {18},
year = {2003}
}
```
"""
function relaxation_psd_to_soc_real_conic(m, mx)
    @assert size(mx,1) == size(mx,2)
    n_elements = size(mx,1)
    for i in 1:n_elements-1
        for j in i+1:n_elements
            JuMP.@constraint(m, [mx[i,i]+mx[j,j], 2*mx[i,j], mx[i,i]-mx[j,j]] in JuMP.SecondOrderCone())
        end
    end
end


"""
SDP to SOC relaxation of type 2, applied to complex-value matrix,  as described in:
```
@article{Kim2003,
author = {Kim, S and Kojima, M and Yamashita, M},
title = {{Second order cone programming relaxation of a positive semidefinite constraint}},
doi = {10.1080/1055678031000148696},
journal = {Optimization Methods and Software},
number = {5},
pages = {535--541},
volume = {18},
year = {2003}
}
```
"""
function relaxation_psd_to_soc_complex_conic(m, mxreal, mximag)
    @assert size(mxreal) == size(mximag)
    n_elements = size(mxreal,1)
    for i in 1:n_elements-1
        for j in i+1:n_elements
            JuMP.@constraint(m, [mxreal[i,i]+mxreal[j,j], 2*mxreal[i,j], 2*mximag[i,j], mxreal[i,i]-mxreal[j,j]] in JuMP.SecondOrderCone())
        end
    end
end


"""
See section 4.3 for complex to real PSD constraint transformation:
@article{Fazel2001,
author = {Fazel, M. and Hindi, H. and Boyd, S.P.},
title = {{A rank minimization heuristic with application to minimum order system approximation}},
doi = {10.1109/ACC.2001.945730},
journal = {Proc. American Control Conf.},
number = {2},
pages = {4734--4739},
url = {http://ieeexplore.ieee.org/lpdocs/epic03/wrapper.htm?arnumber=945730},
volume = {6},
year = {2001}
}
"""
function relaxation_psd_to_soc(m, mxreal, mximag; complex=true)
    if complex==false
        @assert size(mxreal) == size(mximag)
        mx =
            [
            mxreal -mximag;
            mximag  mxreal
            ]

        relaxation_psd_to_soc_real(m, mx)
    else
        relaxation_psd_to_soc_complex(m, mxreal, mximag)
    end
end


"""
See section 4.3 for complex to real PSD constraint transformation:
@article{Fazel2001,
author = {Fazel, M. and Hindi, H. and Boyd, S.P.},
title = {{A rank minimization heuristic with application to minimum order system approximation}},
doi = {10.1109/ACC.2001.945730},
journal = {Proc. American Control Conf.},
number = {2},
pages = {4734--4739},
url = {http://ieeexplore.ieee.org/lpdocs/epic03/wrapper.htm?arnumber=945730},
volume = {6},
year = {2001}
}
"""
function relaxation_psd_to_soc_conic(m, mxreal, mximag; complex=true)
    if complex==false
        @assert size(mxreal) == size(mximag)
        mx =
            [
            mxreal -mximag;
            mximag  mxreal
            ]

        relaxation_psd_to_soc_real_conic(m, mx)
    else
        relaxation_psd_to_soc_complex_conic(m, mxreal, mximag)
    end
end


"""
For debugging / exploration: real-valued SDP to SDP relaxation based on PSDness of principal minors, default is 3x3 SDP relaxation
"""
function relaxation_psd_to_psd_real(m, mxreal, mximag; ndim=3)
    @assert size(mxreal) == size(mximag)
    @assert size(mxreal,1) >= ndim
    n_elements = size(mxreal,1)
    for i in 1:n_elements-(ndim-1)
        j = i+(ndim-1)
        mr = mxreal[i:j, i:j]
        mi = mximag[i:j, i:j]
        JuMP.@constraint(m, [mr -mi; mi mr] in JuMP.PSDCone())
    end
end


"""
This constraint models
M = [ Ar+im*Ai    Cr+im*Ci
     (Cr+im*Ci)'  Br+im*Bi]
M in PSDCone, rank(M) = 1
as a set of polynomial matrix equations
(Cr+im*Ci) (Cr+im*Ci)' = tr(Br)(Ar+im*Ai)
(Cr+im*Ci)'(Cr+im*Ci)  = tr(Ar)(Br+im*Bi)
"""
function reformulation_psd_rank1_to_quadratic(model, Ar, Ai, Br, Bi, Cr, Ci)
    @assert size(Ar) == size(Ai)
    @assert size(Br) == size(Bi)
    @assert size(Cr) == size(Ci)
    @assert size(Cr,1) == size(Ar,1)
    @assert size(Cr,2) == size(Br,2)

    # (Cr+im*Ci) (Cr+im*Ci)' = tr(Br)(Ar+im*Ai)
    matrix_product_real(model, Ar, Br, Cr, Ci)
    matrix_product_imag(model, Ai, Br, Cr, Ci)

    # equations are redundant when size(Cr) == (1,1)
    if size(Cr) != (1, 1)
        # (Cr+im*Ci)'(Cr+im*Ci)  = tr(Ar)(Br+im*Bi)
        # swap A and B, take complex conjugate transpose of C
        matrix_product_real(model, Br, Ar, Cr', -Ci')
        matrix_product_imag(model, Bi, Ar, Cr', -Ci')
    end
end

"""
This constraints models
M = Mr + im*Mi in PSDCone, rank(M) = 1
as a set of quadratic matrix equations by partitioning M in four equal quadrants
Ar+im*Ai the upper diagonal block
Br+im*Bi the lower diagonal block
Cr+im*Ci the upper right block

and then calls reformulation_psd_rank1_to_quadratic(model, Ar, Ai, Br, Bi, Cr, Ci)
"""
function reformulation_psd_rank1_to_quadratic(model, Mr, Mi)
    @assert size(Mr) == size(Mi)
    @assert size(Mr,1) == size(Mr,2)
    @assert iseven(size(Mr,1))

    n = Int(size(Mr,1)/2)

    Ar = Mr[1:n,    1:n]
    Ai = Mi[1:n,    1:n]
    Br = Mr[n+1:2n, n+1:2n]
    Bi = Mi[n+1:2n, n+1:2n]
    Cr = Mr[1:n,    n+1:2n]
    Ci = Mi[1:n,    n+1:2n]

    reformulation_psd_rank1_to_quadratic(model, Ar, Ai, Br, Bi, Cr, Ci)
end

"""
This constraints models the matrix equation
(Cr)(Cr)' + (Ci)(Ci)' = tr(Br)(Ar)
as a set of scalar equations.
The lower triangular elements generate redundant constraints due to symmetry, and are therefore not constructed
"""
function matrix_product_real(model, Ar, Br, Cr, Ci)
    n = size(Cr,1)
    m = size(Cr,2)
    utridiag   = [(i,j) for i=1:n, j=1:n if i<=j] # upper triangle + diagonal elements of A

    for (a,b) in utridiag
            JuMP.@constraint(model, sum(Cr[a,j]*Cr[b,j] + Ci[a,j]*Ci[b,j] for j in 1:m) == Ar[a,b] * sum(Br[j,j] for j in 1:m))
    end
end

"""
This constraints models the matrix equation
(Ci)(Cr)' - (Cr)(Ci)' = tr(Br)(Ai)
as a set of scalar equations.
The lower triangular and diagonal elements generate redundant constraints due to symmetry, and are therefore not constructed
"""
function matrix_product_imag(model, Ai, Br, Cr, Ci)
    n = size(Cr,1)
    m = size(Cr,2)
    utri       = [(i,j) for i=1:n, j=1:n if i< j] # upper triangle elements of A

    for (a,b) in utri
        JuMP.@constraint(model, sum(Ci[a,j]*Cr[b,j] - Cr[a,j]*Ci[b,j] for j in 1:m) == Ai[a,b] * sum(Br[j,j] for j in 1:m))
    end
end
