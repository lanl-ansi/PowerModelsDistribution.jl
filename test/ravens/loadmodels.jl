@info "ravens - running load models tests"

@testset "ravens - test loadmodels pf" begin
    @testset "loadmodels connection variations" begin
        math = transform_data_model(ravens_case3_lm_1230)
        result = solve_mc_pf(math, ACPUPowerModel, ipopt_solver)
        # voltage magnitude at load bus
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]
        @test isapprox(vm(result, string(sourcebus_id)) ./ vbase, [4.330113905476207, 4.330103287702802, 4.330113453504967]; atol=1E-5)
        # single-phase delta loads
        @test isapprox(pd(result, "1")*1000, [400], atol=1E-1)
        @test isapprox(qd(result, "1")*1000, [300], atol=1E-1)
        # single-phase wye loads
        @test isapprox(pd(result, "10")*1000, [133.3, 133.3, 133.3], atol=1E-1)
        @test isapprox(qd(result, "10")*1000, [100.0, 100, 100], atol=1E-1)
        # three-phase loads
        @test isapprox(qd(result, "2")*1000, [100.0, 100.0, 100.0], atol=1E-1)
        @test isapprox(pd(result, "2")*1000, [133.3, 133.3, 133.3], atol=1E-1)
        @test isapprox(pd(result, "6")*1000, [133.3, 133.3, 133.3], atol=1E-1)
        @test isapprox(qd(result, "6")*1000, [100.0, 100.0, 100.0], atol=1E-1)
        @test isapprox(pd(result, "8")*1000, [133.3, 133.3, 133.3], atol=1E-1)
        @test isapprox(qd(result, "8")*1000, [100.0, 100.0, 100.0], atol=1E-1)
    end

    @testset "loadmodels 1/2/5 in acp pf" begin
        math = transform_data_model(ravens_case3_lm_models)
        result = solve_mc_pf(math, ACPUPowerModel, ipopt_solver)
        # voltage magnitude at load bus
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]
        @test isapprox(vm(result, string(sourcebus_id)) ./ vbase, [4.330086349110862, 4.330116429951092, 4.330118057511072], atol=1.5E-4)
        # loads
        @test isapprox(pd(result, "1"), (1/1000)*[400], atol=1e-1)
        @test isapprox(qd(result, "1"), (1/1000)*[300], atol=1e-1)
        @test_broken isapprox(pd(result, "2"), (1/1000)*[332], atol=1e-1)
        @test isapprox(qd(result, "2"), (1/1000)*[249], atol=1e-1)
        @test isapprox(pd(result, "5"), (1/1000)*[110.8, 132.9, 134.1], atol=1e-1)
        @test isapprox(qd(result, "5"), (1/1000)*[ 83.1,  99.7, 100.6], atol=1e-1)
    end

    @testset "loadmodels 1/2/5 in acr pf" begin
        math = transform_data_model(ravens_case3_lm_models)
        result = solve_mc_pf(math, ACRUPowerModel, ipopt_solver)
        # voltage magnitude at load bus
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]
        @test isapprox(calc_vm_acr(result, string(sourcebus_id)) ./ vbase, [4.330086349110862, 4.330116429951092, 4.330118057511072], atol=1.5E-4)
        # loads
        @test isapprox(pd(result, "1"), (1/1000)*[400], atol=1e-1)
        @test isapprox(qd(result, "1"), (1/1000)*[300], atol=1e-1)
        @test_broken isapprox(pd(result, "2"), (1/1000)*[332], atol=1e-1)
        @test isapprox(qd(result, "2"), (1/1000)*[249], atol=1e-1)
        @test isapprox(pd(result, "5"), (1/1000)*[110.8, 132.9, 134.1], atol=1e-1)
        @test isapprox(qd(result, "5"), (1/1000)*[083.1,  99.7, 100.6], atol=1e-1)
    end

    @testset "loadmodels 1/2/5 in ivr pf" begin
        math = transform_data_model(ravens_case3_lm_models)
        result = solve_mc_pf(math, IVRUPowerModel, ipopt_solver)
        # voltage magnitude at load bus
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]
        @test isapprox(calc_vm_acr(result, string(sourcebus_id)) ./ vbase, [4.330086349110861, 4.330116429951092, 4.330118057511072], atol=1.5E-4)
        # loads
        @test isapprox(pd(result, "1"), (1/1000)*[400], atol=1e-1)
        @test isapprox(qd(result, "1"), (1/1000)*[300], atol=1e-1)
        @test_broken isapprox(pd(result, "2"), (1/1000)*[332], atol=1e-1)
        @test isapprox(qd(result, "2"), (1/1000)*[240], atol=1e-1)
        @test isapprox(pd(result, "5"), (1/1000)*[110.8, 132.9, 134.1], atol=1e-1)
        @test isapprox(qd(result, "5"), (1/1000)*[ 83.1,  99.7, 100.6], atol=1e-1)
    end

    @testset "ZIP loadmodels 8" begin
        math = transform_data_model(ravens_case3_unbalanced_ZIPloads)
        result = solve_mc_pf(math, ACPUPowerModel, ipopt_solver)
        # voltage magnitude at load bus
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]
        # delta loads
        @test isapprox(pd(result, "1"), [0.01], atol=1E-3)
        @test isapprox(qd(result, "1"), [0.003], atol=1E-3)
        # wye loads
        @test isapprox(pd(result, "2"), [0.01], atol=1E-3)
        @test isapprox(qd(result, "2"), [0.003], atol=1E-3)
        @test isapprox(pd(result, "3"), [0.01], atol=1E-3)
        @test isapprox(qd(result, "3"), [0.003], atol=1E-3)
        @test isapprox(pd(result, "4"), [0.003, 0.003, 0.003], atol=1E-3)

        @test isapprox(qd(result, "4"), [0.001, 0.001, 0.001], atol=1E-3)  
    end
end
