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
            @constraint(m, mx[i,j]^2 <= mx[i,i]*mx[j,j])
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
            @constraint(m, mxreal[i,j]^2 + mximag[i,j]^2 <= mxreal[i,i]*mxreal[j,j])
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
            @constraint(m, norm([2*mx[i,j], mx[i,i]-mx[j,j]]) <= mx[i,i]+mx[j,j])
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
            @constraint(m, norm([2*mxreal[i,j], 2*mximag[i,j], mxreal[i,i]-mxreal[j,j]]) <= mxreal[i,i]+mxreal[j,j])
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


#for debugging / exploration:
"""
real-valued SDP to SDP relaxation based on PSDness of principal minors, default is 3x3 SDP relaxation
"""
function relaxation_psd_to_psd_real(m, mxreal, mximag; ndim=3)
    @assert size(mxreal) == size(mximag)
    @assert size(mxreal,1) >= ndim
    n_elements = size(mxreal,1)
    for i in 1:n_elements-(ndim-1)
        j = i+(ndim-1)
        mr = mxreal[i:j, i:j]
        mi = mximag[i:j, i:j]
        @SDconstraint(m, [mr -mi; mi mr] >=0)
    end
end
