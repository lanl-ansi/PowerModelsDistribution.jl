# SOC to SDP relaxation of type 2 as described in:
# Kim, S., Kojima, M., & Yamashita, M. (2003). Second order cone programming relaxation of a positive semidefinite constraint. Optimization Methods and Software, 18(5), 535–541. https://doi.org/10.1080/1055678031000148696
function psd_to_soc(pm::GenericPowerModel{T}, mat) where T<:AbstractConicUBFForm
    assert(size(mat,1) == size(mat,2))
    n_cond = 3
    for i in 1:length(diag(mat))-1
        for j in i+1:length(diag(mat))
            if j!=i+2*n_cond # trivial constraints
                @constraint(pm.model, norm([2*mat[i,j], mat[i,i]-mat[j,j]]) <= mat[i,i]+mat[j,j])
            end
        end
    end
end

# SOC to SDP relaxation of type 2 as described in:
# Kim, S., Kojima, M., & Yamashita, M. (2003). Second order cone programming relaxation of a positive semidefinite constraint. Optimization Methods and Software, 18(5), 535–541. https://doi.org/10.1080/1055678031000148696
function psd_to_soc(pm::GenericPowerModel{T}, mat) where T<:AbstractNLPUBFForm
    assert(size(mat,1) == size(mat,2))
    n_cond = 3
    for i in 1:length(diag(mat))-1
        for j in i+1:length(diag(mat))
            if j!=i+2*n_cond # trivial constraints
                @constraint(pm.model, mat[i,j]^2 <= mat[i,i]*mat[j,j])
            end
        end
    end
end

# See section 4.3 in:
#Fazel, M., Hindi, H., & Boyd, S. P. (2001). A rank minimization heuristic with application to minimum order system approximation. Proc. American Control Conf., 6(2), 4734–4739. https://doi.org/10.1109/ACC.2001.945730
function psd_to_soc(mp, varreal, varimag)
    assert(size(varreal) == size(varimag))
    mat =
        [
        varreal -varimag;
        varimag varreal
        ]

    psd_to_soc(mp, mat)
end

# convexification of complex power definition in indivual phases
function psd_to_soc_diag(pm::GenericPowerModel{T}, ps, qs, wi, ccms) where T<:AbstractConicUBFForm
    for cond in 1:length(diag(ps))
        @constraint(pm.model, norm([2*ps[cond,cond],  2*qs[cond,cond],  wi[cond,cond]-ccms[cond,cond]]) <=  wi[cond,cond]+ccms[cond,cond])
    end
end

# convexification of complex power definition in indivual phases
function psd_to_soc_diag(pm::GenericPowerModel{T}, ps, qs, wi, ccms) where T<:AbstractNLPUBFForm
    for cond in 1:length(diag(ps))
        @constraint(pm.model, ps[cond,cond]^2 + qs[cond,cond]^2 <=  wi[cond,cond]*ccms[cond,cond])
    end
end
