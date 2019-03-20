using ThreePhasePowerModels
TPPMs = ThreePhasePowerModels
using PowerModels
PMs = PowerModels
using Ipopt
using JuMP
using Compat.Test

@testset "validate loss model" begin
    file = "data/ut_trans_1.dss"
    data_tppm = TPPMs.parse_file(file)
    pm = PMs.build_generic_model(data_tppm, PMs.ACPPowerModel, TPPMs.post_tp_pf, multiconductor=true)
    # lower bound is needed for convergence
    for i in ids(pm,  :bus)
        for c in 1:3
            setlowerbound(var(pm, pm.cnw, c, :vm)[i], 0.5)
        end
    end
    sol = PMs.solve_generic_model(pm, Ipopt.IpoptSolver())
    sol
    sol["solution"]["gen"]["1"]["pg"]
    ##
    vm4 = sol["solution"]["bus"]["4"]["vm"]
    @test isapprox(vm4, [0.9147, 0.8723, 0.8701], atol=1E-4)
    va4 = sol["solution"]["bus"]["4"]["va"]
    @test isapprox(round.(va4*180/pi; digits=1), [30.1, -90.7, 151.2], atol=1E-2)
    vm3 = sol["solution"]["bus"]["3"]["vm"]
    @test isapprox(vm3, [0.9134, 0.8718, 0.8701], atol=1E-4)
    va3 = sol["solution"]["bus"]["3"]["va"]
    @test isapprox(round.(va3*180/pi; digits=1), [30.2, -90.7, 151.2], atol=1E-2)
end
