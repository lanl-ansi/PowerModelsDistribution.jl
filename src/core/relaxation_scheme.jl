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
function relaxation_psd_to_soc(m::JuMP.Model, mxreal, mximag; complex::Bool=true)
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
A valid inequality for the product of two complex variables with magnitude and
angle difference bounds.
In the literature this constraints are called the Lifted Nonlinear Cuts (LNCs).
@misc{1512.04644,
    Author = {Carleton Coffrin and Hassan Hijazi and Pascal Van Hentenryck},
    Title = {Strengthening the SDP Relaxation of AC Power Flows with Convex
        Envelopes, Bound Tightening, and Lifted Nonlinear Cuts},
    Year = {2015},
    Eprint = {arXiv:1512.04644},
}
"""
function cut_complex_product_and_angle_difference(m, wf, wt, wr, wi, angmin, angmax)
    @assert angmin >= -pi/2 && angmin <= pi/2
    @assert angmax >= -pi/2 && angmax <= pi/2
    @assert angmin < angmax

    wf_lb, wf_ub = _IM.variable_domain(wf)
    wt_lb, wt_ub = _IM.variable_domain(wt)

    vf_lb, vf_ub = sqrt(wf_lb), sqrt(wf_ub)
    vt_lb, vt_ub = sqrt(wt_lb), sqrt(wt_ub)
    td_ub = angmax
    td_lb = angmin

    phi = (td_ub + td_lb)/2
    d   = (td_ub - td_lb)/2

    sf = vf_lb + vf_ub
    st = vt_lb + vt_ub

    JuMP.@constraint(m, sf*st*(cos(phi)*wr + sin(phi)*wi) - vt_ub*cos(d)*st*wf - vf_ub*cos(d)*sf*wt >=  vf_ub*vt_ub*cos(d)*(vf_lb*vt_lb - vf_ub*vt_ub))
    JuMP.@constraint(m, sf*st*(cos(phi)*wr + sin(phi)*wi) - vt_lb*cos(d)*st*wf - vf_lb*cos(d)*sf*wt >= -vf_lb*vt_lb*cos(d)*(vf_lb*vt_lb - vf_ub*vt_ub))
end
