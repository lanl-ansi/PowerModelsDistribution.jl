@info "running load models tests"

@testset "test loadmodels pf" begin
    @testset "loadmodels connection variations" begin
        result = solve_mc_pf(case3_lm_1230, ACPUPowerModel, ipopt_solver)
        # voltage magnitude at load bus
        vbase = case3_lm_1230["settings"]["vbases_default"]["sourcebus"]
        @test isapprox(vm(result, "loadbus") ./ vbase, [0.999993, 0.999992, 0.999993]; atol=1E-5)
        # single-phase delta loads
        @test isapprox(pd(result, "d1ph23"), [0.0, 286.6, 113.4], atol=1E-1)
        @test isapprox(qd(result, "d1ph23"), [0.0, 34.5, 265.5], atol=1E-1)
        # single-phase wye loads
        @test isapprox(pd(result, "y1ph2"), [400.0], atol=1E-1)
        @test isapprox(qd(result, "y1ph2"), [300.0], atol=1E-1)
        # three-phase loads
        @test isapprox(pd(result, "d3ph"), [133.3, 133.3, 133.3], atol=1E-1)
        @test isapprox(qd(result, "d3ph"), [100.0, 100.0, 100.0], atol=1E-1)
        @test isapprox(pd(result, "d3ph123"), [133.3, 133.3, 133.3], atol=1E-1)
        @test isapprox(qd(result, "d3ph123"), [100.0, 100.0, 100.0], atol=1E-1)
        @test isapprox(pd(result, "d3ph213"), [133.3, 133.3, 133.3], atol=1E-1)
        @test isapprox(qd(result, "d3ph213"), [100.0, 100.0, 100.0], atol=1E-1)
    end

    @testset "loadmodels 1/2/5 in acp pf" begin
        result = solve_mc_pf(case3_lm_models, ACPUPowerModel, ipopt_solver)
        # voltage magnitude at load bus
        vbase = case3_lm_models["settings"]["vbases_default"]["sourcebus"]
        @test isapprox(vm(result, "loadbus") ./ vbase, [0.83072, 0.99653, 1.0059], atol=1.5E-4)
        # delta and wye single-phase load models
        @test isapprox(pd(result, "y1phm1"), [400.0], atol=1E-1)
        @test isapprox(qd(result, "y1phm1"), [300.0], atol=1E-1)
        @test isapprox(pd(result, "y1phm2"), [278.3], atol=1E-1)
        @test isapprox(qd(result, "y1phm2"), [208.7], atol=1E-1)
        @test isapprox(pd(result, "y1phm5"), [333.6], atol=1E-1)
        @test isapprox(qd(result, "y1phm5"), [250.2], atol=1E-1)
        # delta three-phase loads
        @test isapprox(pd(result, "d3phm1"), [116.0, 146.5, 137.5], atol=1e-1)
        @test isapprox(qd(result, "d3phm1"), [ 89.6,  97.7, 112.7], atol=1e-1)
        @test isapprox(pd(result, "d3phm2"), [100.5, 134.8, 121.2], atol=1e-1)
        @test isapprox(qd(result, "d3phm2"), [ 77.1,  85.4, 104.8], atol=1e-1)
        @test isapprox(pd(result, "d3phm5"), [108.0, 140.5, 129.1], atol=1e-1)
        @test isapprox(qd(result, "d3phm5"), [ 83.1,  91.4, 108.7], atol=1e-1)
        # wye three-phase loads
        @test isapprox(pd(result, "y3phm1"), [133.3, 133.3, 133.3], atol=1e-1)
        @test isapprox(qd(result, "y3phm1"), [100.0, 100.0, 100.0], atol=1e-1)
        @test isapprox(pd(result, "y3phm2"), [ 92.0, 132.4, 134.9], atol=1e-1)
        @test isapprox(qd(result, "y3phm2"), [ 69.0,  99.3, 101.2], atol=1e-1)
        @test isapprox(pd(result, "y3phm5"), [110.8, 132.9, 134.1], atol=1e-1)
        @test isapprox(qd(result, "y3phm5"), [ 83.1,  99.7, 100.6], atol=1e-1)
    end

    @testset "loadmodels 1/2/5 in acr pf" begin
        result = solve_mc_pf(case3_lm_models, ACRUPowerModel, ipopt_solver)
        # voltage magnitude at load bus
        vbase = case3_lm_models["settings"]["vbases_default"]["sourcebus"]
        @test isapprox(calc_vm_acr(result, "loadbus") ./ vbase, [0.83072, 0.99653, 1.0059], atol=1.5E-4)
        # delta and wye single-phase load models
        @test isapprox(pd(result, "y1phm1"), [400.0], atol=1e-1)
        @test isapprox(qd(result, "y1phm1"), [300.0], atol=1e-1)
        @test isapprox(pd(result, "y1phm2"), [278.3], atol=1e-1)
        @test isapprox(qd(result, "y1phm2"), [208.7], atol=1e-1)
        @test isapprox(pd(result, "y1phm5"), [333.6], atol=1e-1)
        @test isapprox(qd(result, "y1phm5"), [250.2], atol=1e-1)
        # delta three-phase loads
        @test isapprox(pd(result, "d3phm1"), [116.0, 146.5, 137.5], atol=1e-1)
        @test isapprox(qd(result, "d3phm1"), [ 89.6,  97.7, 112.7], atol=1e-1)
        @test isapprox(pd(result, "d3phm2"), [100.5, 134.8, 121.2], atol=1e-1)
        @test isapprox(qd(result, "d3phm2"), [ 77.1,  85.4, 104.8], atol=1e-1)
        @test isapprox(pd(result, "d3phm5"), [108.0, 140.5, 129.1], atol=1e-1)
        @test isapprox(qd(result, "d3phm5"), [ 83.1,  91.4, 108.7], atol=1e-1)
        # wye three-phase loads
        @test isapprox(pd(result, "y3phm1"), [133.3, 133.3, 133.3], atol=1e-1)
        @test isapprox(qd(result, "y3phm1"), [100.0, 100.0, 100.0], atol=1e-1)
        @test isapprox(pd(result, "y3phm2"), [092.0, 132.4, 134.9], atol=1e-1)
        @test isapprox(qd(result, "y3phm2"), [069.0,  99.3, 101.2], atol=1e-1)
        @test isapprox(pd(result, "y3phm5"), [110.8, 132.9, 134.1], atol=1e-1)
        @test isapprox(qd(result, "y3phm5"), [083.1,  99.7, 100.6], atol=1e-1)
    end

    @testset "loadmodels 1/2/5 in ivr pf" begin
        result = solve_mc_pf(case3_lm_models, IVRUPowerModel, ipopt_solver)
        # voltage magnitude at load bus
        vbase = case3_lm_models["settings"]["vbases_default"]["sourcebus"]
        @test isapprox(calc_vm_acr(result, "loadbus") ./ vbase, [0.83072, 0.99653, 1.0059], atol=1.5E-4)
        # delta and wye single-phase load models
        @test isapprox(pd(result, "y1phm1"), [400.0], atol=1e-1)
        @test isapprox(qd(result, "y1phm1"), [300.0], atol=1e-1)
        @test isapprox(pd(result, "y1phm2"), [278.3], atol=1e-1)
        @test isapprox(qd(result, "y1phm2"), [208.7], atol=1e-1)
        @test isapprox(pd(result, "y1phm5"), [333.6], atol=1e-1)
        @test isapprox(qd(result, "y1phm5"), [250.2], atol=1e-1)
        # delta three-phase loads
        @test isapprox(pd(result, "d3phm1"), [116.0, 146.5, 137.5], atol=1e-1)
        @test isapprox(qd(result, "d3phm1"), [ 89.6,  97.7, 112.7], atol=1e-1)
        @test isapprox(pd(result, "d3phm2"), [100.5, 134.8, 121.2], atol=1e-1)
        @test isapprox(qd(result, "d3phm2"), [ 77.1,  85.4, 104.8], atol=1e-1)
        @test isapprox(pd(result, "d3phm5"), [108.0, 140.5, 129.1], atol=1e-1)
        @test isapprox(qd(result, "d3phm5"), [ 83.1,  91.4, 108.7], atol=1e-1)
        # wye three-phase loads
        @test isapprox(pd(result, "y3phm1"), [133.3, 133.3, 133.3], atol=1e-1)
        @test isapprox(qd(result, "y3phm1"), [100.0, 100.0, 100.0], atol=1e-1)
        @test isapprox(pd(result, "y3phm2"), [ 92.0, 132.4, 134.9], atol=1e-1)
        @test isapprox(qd(result, "y3phm2"), [ 69.0,  99.3, 101.2], atol=1e-1)
        @test isapprox(pd(result, "y3phm5"), [110.8, 132.9, 134.1], atol=1e-1)
        @test isapprox(qd(result, "y3phm5"), [ 83.1,  99.7, 100.6], atol=1e-1)
    end

    @testset "ZIP loadmodels 8" begin
        result = solve_mc_pf(case3_unbalanced_ZIPloads, ACPUPowerModel, ipopt_solver)
        # voltage magnitude at load bus
        vbase = case3_unbalanced_ZIPloads["settings"]["vbases_default"]["sourcebus"]
        @test isapprox(vm(result, "loadbus") ./ vbase, [0.96062, 0.96074, 0.96169]; atol=1E-5)
        # delta loads
        @test isapprox(pd(result, "l1"), [2.9, 2.9, 2.9], atol=1E-1)
        @test isapprox(qd(result, "l1"), [1.0, 1.0, 1.0], atol=1E-1)
        # wye loads
        @test isapprox(pd(result, "l2"), [8.7], atol=1E-1)
        @test isapprox(qd(result, "l2"), [2.9], atol=1E-1)
        @test isapprox(pd(result, "l3"), [8.7], atol=1E-1)
        @test isapprox(qd(result, "l3"), [2.9], atol=1E-1)
        @test isapprox(pd(result, "l4"), [8.5], atol=1E-1)
        @test isapprox(qd(result, "l4"), [2.8], atol=1E-1) 
    end
end
