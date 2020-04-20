
"SDP relaxation of BIM"
function constraint_mc_model_voltage(pm::_PMs.AbstractWRMModel, n::Int)
    for (l,i,j) in _PMs.ref(pm, n, :arcs_from)
        Wr_fr = _PMs.var(pm, n, :Wr)[i]
        Wi_fr = _PMs.var(pm, n, :Wi)[i]

        Wr_to = _PMs.var(pm, n, :Wr)[j]
        Wi_to = _PMs.var(pm, n, :Wi)[j]

        Wijr = _PMs.var(pm, n, :Wijr)[(i,j)]
        Wiji = _PMs.var(pm, n, :Wiji)[(i,j)]

        mat_real = [
        Wr_fr     Wijr  ;
        Wijr'    Wr_to  ;
        ]

        mat_imag = [
        Wi_fr     Wiji  ;
        -Wiji'    Wi_to  ;
        ]

        JuMP.@constraint(pm.model,
        [
        mat_real  mat_imag;
        -mat_imag   mat_real
        ] in JuMP.PSDCone())

        P_fr = _PMs.var(pm, n, :P)[(l,i,j)]
        P_to = _PMs.var(pm, n, :P)[(l,j,i)]
        #cuts improves convergence
        JuMP.@constraint(pm.model, tr(P_fr) + tr(P_to) >= 0 )
    end
end


"Add explicit PSD-ness of W for nodes where it is not implied"
function constraint_mc_voltage_psd(pm::_PMs.AbstractWRMModel, n::Int, i)
    Wr = _PMs.var(pm, n, :Wr)[i]
    Wi = _PMs.var(pm, n, :Wi)[i]

    JuMP.@constraint(pm.model,
    [
    Wr  Wi;
    -Wi   Wr
    ] in JuMP.PSDCone())
end
