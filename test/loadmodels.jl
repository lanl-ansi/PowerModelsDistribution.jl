@info "running load models tests"

@testset "test loadmodels pf" begin
    @testset "loadmodels connection variations" begin
        pmd = parse_file("../test/data/opendss/case3_lm_1230.dss")
        sol = run_mc_pf(pmd, ACPPowerModel, ipopt_solver; make_si=false)
        # voltage magnitude at load bus
        @test isapprox(vm(sol, "loadbus"), [0.999993, 0.999992, 0.999993], atol=1E-5)
        # single-phase delta loads
        @test isapprox(pd(sol, "d1ph23"), [0.0, 0.2866, 0.1134], atol=1E-4)
        @test isapprox(qd(sol, "d1ph23"), [0.0, 0.0345, 0.2655], atol=1E-4)
        # single-phase wye loads
        @test isapprox(pd(sol, "y1ph2"), [0.4000], atol=1E-4)
        @test isapprox(qd(sol, "y1ph2"), [0.3000], atol=1E-4)
        # three-phase loads
        @test isapprox(pd(sol, "d3ph"), [0.1333, 0.1333, 0.1333], atol=1E-4)
        @test isapprox(qd(sol, "d3ph"), [0.100, 0.100, 0.100], atol=1E-4)
        @test isapprox(pd(sol, "d3ph123"), [0.1333, 0.1333, 0.1333], atol=1E-4)
        @test isapprox(qd(sol, "d3ph123"), [0.100, 0.100, 0.100], atol=1E-4)
        @test isapprox(pd(sol, "d3ph213"), [0.1333, 0.1333, 0.1333], atol=1E-4)
        @test isapprox(qd(sol, "d3ph213"), [0.100, 0.100, 0.100], atol=1E-4)
    end
    @testset "loadmodels 1/2/5 in acp pf" begin
        pmd = parse_file("../test/data/opendss/case3_lm_models.dss")
        sol = run_mc_pf(pmd, ACPPowerModel, ipopt_solver; make_si=false)
        # voltage magnitude at load bus
        @test isapprox(vm(sol, "loadbus"), [0.83072, 0.99653, 1.0059], atol=1.5E-4)
        # delta and wye single-phase load models
        @test isapprox(pd(sol, "y1phm1"), [0.4000], atol=1E-4)
        @test isapprox(qd(sol, "y1phm1"), [0.3000], atol=1E-4)
        @test isapprox(pd(sol, "y1phm2"), [0.2783], atol=1E-4)
        @test isapprox(qd(sol, "y1phm2"), [0.2087], atol=1E-4)
        @test isapprox(pd(sol, "y1phm5"), [0.3336], atol=1E-4)
        @test isapprox(qd(sol, "y1phm5"), [0.2502], atol=1E-4)
        # delta three-phase loads
        @test isapprox(pd(sol, "d3phm1"), [0.1160, 0.1465, 0.1375], atol=1E-4)
        @test isapprox(qd(sol, "d3phm1"), [0.0896, 0.0977, 0.1127], atol=1E-4)
        @test isapprox(pd(sol, "d3phm2"), [0.1005, 0.1348, 0.1212], atol=1E-4)
        @test isapprox(qd(sol, "d3phm2"), [0.0771, 0.0854, 0.1048], atol=1E-4)
        @test isapprox(pd(sol, "d3phm5"), [0.1080, 0.1405, 0.1291], atol=1E-4)
        @test isapprox(qd(sol, "d3phm5"), [0.0831, 0.0914, 0.1087], atol=1E-4)
        # wye three-phase loads
        @test isapprox(pd(sol, "y3phm1"), [0.1333, 0.1333, 0.1333], atol=1E-4)
        @test isapprox(qd(sol, "y3phm1"), [0.1000, 0.1000, 0.1000], atol=1E-4)
        @test isapprox(pd(sol, "y3phm2"), [0.0920, 0.1324, 0.1349], atol=1E-4)
        @test isapprox(qd(sol, "y3phm2"), [0.0690, 0.0993, 0.1012], atol=1E-4)
        @test isapprox(pd(sol, "y3phm5"), [0.1108, 0.1329, 0.1341], atol=1E-4)
        @test isapprox(qd(sol, "y3phm5"), [0.0831, 0.0997, 0.1006], atol=1E-4)
    end
    @testset "loadmodels 1/2/5 in acr pf" begin
        pmd = parse_file("../test/data/opendss/case3_lm_models.dss")
        sol = run_mc_pf(pmd, ACRPowerModel, ipopt_solver; make_si=false)
        # voltage magnitude at load bus
        @test isapprox(calc_vm_acr(sol, "loadbus"), [0.83072, 0.99653, 1.0059], atol=1.5E-4)
        # delta and wye single-phase load models
        @test isapprox(pd(sol, "y1phm1"), [0.4000], atol=1E-4)
        @test isapprox(qd(sol, "y1phm1"), [0.3000], atol=1E-4)
        @test isapprox(pd(sol, "y1phm2"), [0.2783], atol=1E-4)
        @test isapprox(qd(sol, "y1phm2"), [0.2087], atol=1E-4)
        @test isapprox(pd(sol, "y1phm5"), [0.3336], atol=1E-4)
        @test isapprox(qd(sol, "y1phm5"), [0.2502], atol=1E-4)
        # delta three-phase loads
        @test isapprox(pd(sol, "d3phm1"), [0.1160, 0.1465, 0.1375], atol=1E-4)
        @test isapprox(qd(sol, "d3phm1"), [0.0896, 0.0977, 0.1127], atol=1E-4)
        @test isapprox(pd(sol, "d3phm2"), [0.1005, 0.1348, 0.1212], atol=1E-4)
        @test isapprox(qd(sol, "d3phm2"), [0.0771, 0.0854, 0.1048], atol=1E-4)
        @test isapprox(pd(sol, "d3phm5"), [0.1080, 0.1405, 0.1291], atol=1E-4)
        @test isapprox(qd(sol, "d3phm5"), [0.0831, 0.0914, 0.1087], atol=1E-4)
        # wye three-phase loads
        @test isapprox(pd(sol, "y3phm1"), [0.1333, 0.1333, 0.1333], atol=1E-4)
        @test isapprox(qd(sol, "y3phm1"), [0.1000, 0.1000, 0.1000], atol=1E-4)
        @test isapprox(pd(sol, "y3phm2"), [0.0920, 0.1324, 0.1349], atol=1E-4)
        @test isapprox(qd(sol, "y3phm2"), [0.0690, 0.0993, 0.1012], atol=1E-4)
        @test isapprox(pd(sol, "y3phm5"), [0.1108, 0.1329, 0.1341], atol=1E-4)
        @test isapprox(qd(sol, "y3phm5"), [0.0831, 0.0997, 0.1006], atol=1E-4)
    end
    @testset "loadmodels 1/2/5 in ivr pf" begin
        pmd = parse_file("../test/data/opendss/case3_lm_models.dss")
        sol = run_mc_pf(pmd, IVRPowerModel, ipopt_solver; make_si=false)
        # voltage magnitude at load bus
        @test isapprox(calc_vm_acr(sol, "loadbus"), [0.83072, 0.99653, 1.0059], atol=1.5E-4)
        # delta and wye single-phase load models
        @test isapprox(pd(sol, "y1phm1"), [0.4000], atol=1E-4)
        @test isapprox(qd(sol, "y1phm1"), [0.3000], atol=1E-4)
        @test isapprox(pd(sol, "y1phm2"), [0.2783], atol=1E-4)
        @test isapprox(qd(sol, "y1phm2"), [0.2087], atol=1E-4)
        @test isapprox(pd(sol, "y1phm5"), [0.3336], atol=1E-4)
        @test isapprox(qd(sol, "y1phm5"), [0.2502], atol=1E-4)
        # delta three-phase loads
        @test isapprox(pd(sol, "d3phm1"), [0.1160, 0.1465, 0.1375], atol=1E-4)
        @test isapprox(qd(sol, "d3phm1"), [0.0896, 0.0977, 0.1127], atol=1E-4)
        @test isapprox(pd(sol, "d3phm2"), [0.1005, 0.1348, 0.1212], atol=1E-4)
        @test isapprox(qd(sol, "d3phm2"), [0.0771, 0.0854, 0.1048], atol=1E-4)
        @test isapprox(pd(sol, "d3phm5"), [0.1080, 0.1405, 0.1291], atol=1E-4)
        @test isapprox(qd(sol, "d3phm5"), [0.0831, 0.0914, 0.1087], atol=1E-4)
        # wye three-phase loads
        @test isapprox(pd(sol, "y3phm1"), [0.1333, 0.1333, 0.1333], atol=1E-4)
        @test isapprox(qd(sol, "y3phm1"), [0.1000, 0.1000, 0.1000], atol=1E-4)
        @test isapprox(pd(sol, "y3phm2"), [0.0920, 0.1324, 0.1349], atol=1E-4)
        @test isapprox(qd(sol, "y3phm2"), [0.0690, 0.0993, 0.1012], atol=1E-4)
        @test isapprox(pd(sol, "y3phm5"), [0.1108, 0.1329, 0.1341], atol=1E-4)
        @test isapprox(qd(sol, "y3phm5"), [0.0831, 0.0997, 0.1006], atol=1E-4)
    end
end
