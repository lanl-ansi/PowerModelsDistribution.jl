"""
SDP to SOC relaxation of type 2 as described in:
Kim, S., Kojima, M., & Yamashita, M. (2003). Second order cone programming relaxation of a positive semidefinite constraint. Optimization Methods and Software, 18(5), 535–541. https://doi.org/10.1080/1055678031000148696
"""

function psd_to_soc(pm::GenericPowerModel{T}, mat) where T<:AbstractConicUBFForm
    assert(size(mat,1) == size(mat,2))
    n_cond = Int(size(mat,1)/4)
    for i in 1:n_cond-1
        for j in i+1:n_cond
            # if j!=i+2*n_cond # trivial constraints
                @constraint(pm.model, norm([2*mat[i,j], mat[i,i]-mat[j,j]]) <= mat[i,i]+mat[j,j])
            # end
        end
    end
end

"""
SDP to SOC relaxation of type 2 as described in:
Kim, S., Kojima, M., & Yamashita, M. (2003). Second order cone programming relaxation of a positive semidefinite constraint. Optimization Methods and Software, 18(5), 535–541. https://doi.org/10.1080/1055678031000148696
"""
function psd_to_soc(pm::GenericPowerModel{T}, mat) where T<:AbstractNLPUBFForm
    assert(size(mat,1) == size(mat,2))
    n_cond = Int(size(mat,1)/4)
    for i in 1:n_cond-1
        for j in i+1:n_cond
            # if j!=i+2*n_cond # trivial constraints
                @constraint(pm.model, mat[i,j]^2 <= mat[i,i]*mat[j,j])
            # end
        end
    end
end

"""
SOC relaxation as described in:
Andersen, M. S., Hansson, A., & Vandenberghe, L. (2014). Reduced-complexity semidefinite relaxations of optimal power flow problems. IEEE Trans. Power Syst., 29(4), 1855–1863.
"""
function psd_to_soc_complex(pm::GenericPowerModel{T}, matreal, matimag) where T<:AbstractConicUBFForm
    assert(size(matreal) == size(matimag))
    n_cond = Int(size(matreal,1)/2)
    for i in 1:n_cond-1
        for j in i+1:n_cond
            # if j!=i+n_cond # trivial constraints?
                @constraint(pm.model, norm([2*matreal[i,j], 2*matimag[i,j], matreal[i,i]-matreal[j,j]]) <= matreal[i,i]+matreal[j,j])
            # end
        end
    end
end

"""
SOC relaxation as described in:
Andersen, M. S., Hansson, A., & Vandenberghe, L. (2014). Reduced-complexity semidefinite relaxations of optimal power flow problems. IEEE Trans. Power Syst., 29(4), 1855–1863.
"""
function psd_to_soc_complex(pm::GenericPowerModel{T}, matreal, matimag) where T<:AbstractNLPUBFForm
    assert(size(matreal) == size(matimag))
    n_cond = Int(size(matreal,1)/2)
    for i in 1:n_cond-1
        for j in i+1:n_cond
            # if j!=i+n_cond # trivial constraints?
                @constraint(pm.model, matreal[i,j]^2 + matimag[i,j]^2 <= matreal[i,i]*matreal[j,j])
            # end
        end
    end
end


"""
See section 4.3 in:
Fazel, M., Hindi, H., & Boyd, S. P. (2001). A rank minimization heuristic with application to minimum order system approximation. Proc. American Control Conf., 6(2), 4734–4739. https://doi.org/10.1109/ACC.2001.945730
"""
function psd_to_soc(pm::GenericPowerModel, matrixreal, matriximag; complex=true)
    if complex==false
        assert(size(matrixreal) == size(matriximag))
        matrix =
            [
            matrixreal -matriximag;
            matriximag  matrixreal
            ]

        psd_to_soc(pm, matrix)
    else
        psd_to_soc_complex(pm, matrixreal, matriximag)
    end
end
