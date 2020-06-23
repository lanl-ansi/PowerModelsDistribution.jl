@info "running minimum load delta (mld) tests"

@testset "test mld" begin
    case5 = PM.parse_file("$(pms_path)/test/data/matpower/case5.m"); make_multiconductor!(case5, 3)
    case5_strg = PM.parse_file("$(pms_path)/test/data/matpower/case5_strg.m"); make_multiconductor!(case5_strg, 3)
    case3_ml = PM.parse_file("../test/data/matpower/case3_ml.m"); make_multiconductor!(case3_ml, 3)
    ut_trans_2w_yy = PowerModelsDistribution.parse_file("../test/data/opendss/ut_trans_2w_yy.dss")

    @testset "5-bus acp mld" begin
        result = run_mc_mld(case5, ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.3377; atol = 1e-4)
        @test all(isapprox(result["solution"]["load"]["1"]["pd"], [3.0, 3.0, 3.0]; atol=1e-4))
        @test all(isapprox(result["solution"]["load"]["1"]["qd"], [0.9861, 0.9861, 0.9861]; atol=1e-4))

        @test isapprox(result["solution"]["load"]["1"]["status"], 1.000; atol = 1e-3)
    end

    @testset "5-bus storage acp mld" begin
        result = run_mc_mld(case5_strg, ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.3432; atol=1e-2)

        @test isapprox(result["solution"]["load"]["1"]["status"], 1.000; atol=1e-3)
    end

    @testset "5-bus nfa mld" begin
        result = run_mc_mld(case5, NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.1557; atol = 1e-4)

        @test isapprox(result["solution"]["load"]["1"]["status"], 1.000; atol = 1e-3)
    end

    @testset "5-bus storage nfa mld" begin
        result = run_mc_mld(case5_strg, NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.1557; atol=1e-4)

        @test isapprox(result["solution"]["load"]["1"]["status"], 1.000; atol=1e-3)
    end

    @testset "5-bus mn acp mld" begin
        case5_mn = InfrastructureModels.replicate(case5, 3, Set(["per_unit"]))
        result = run_mn_mc_mld(case5_mn, ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 3*0.3377; atol=2.0e-4)
        @test all(isapprox(result["solution"]["nw"]["1"]["load"]["1"]["pd"], [3.0, 3.0, 3.0]; atol=1.0e-4))
        @test all(isapprox(result["solution"]["nw"]["2"]["load"]["1"]["qd"], [0.9861, 0.9861, 0.9861]; atol=1.0e-4))
        @test isapprox(result["solution"]["nw"]["3"]["load"]["1"]["status"], 1.0; atol=1.0e-3)
    end

    @testset "5-bus mn storage acp mld" begin
        case5_strg_mn = InfrastructureModels.replicate(case5_strg, 3, Set(["per_unit"]))
        result = run_mn_mc_mld(case5_strg_mn, ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test result["objective"] >= 3.0*0.3432
        @test isapprox(result["solution"]["nw"]["1"]["load"]["1"]["status"], 1.0; atol=1.0e-3)
    end

    @testset "5-bus mn nfa mld" begin
        case5_mn = InfrastructureModels.replicate(case5, 3, Set(["per_unit"]))
        result = run_mn_mc_mld(case5_mn, NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 3.0*0.1557; atol=2.0e-4)
        @test isapprox(result["solution"]["nw"]["1"]["load"]["1"]["status"], 1.0; atol=1.0e-3)
    end

    @testset "5-bus mn storage nfa mld" begin
        case5_strg_mn = InfrastructureModels.replicate(case5_strg, 3, Set(["per_unit"]))
        result = run_mn_mc_mld(case5_strg_mn, NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test result["objective"] >= 3.0*0.1557
        @test isapprox(result["solution"]["nw"]["1"]["load"]["1"]["status"], 1.0; atol=1.0e-3)
    end

    @testset "5-bus mn lpubfdiag mld" begin
        case5_mn = InfrastructureModels.replicate(case5, 3, Set(["per_unit"]))
        result = run_mn_mc_mld(case5_mn, LPUBFDiagPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 3.0*0.1557; atol=2.0e-4)
        @test isapprox(result["solution"]["nw"]["1"]["load"]["1"]["status"], 1.0; atol=1.0e-3)
    end

    @testset "3-bus nfa mld" begin
        mp_data = PM.parse_file("../test/data/matpower/case3_ml.m"); make_multiconductor!(mp_data, 3)
        result = run_mc_mld(mp_data, NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 24.9; atol = 1e-1)

        @test isapprox(result["solution"]["load"]["1"]["status"], 0.689; atol = 1e-3)
    end

    @testset "3-bus shunt nfa mld" begin
        mp_data = PM.parse_file("../test/data/matpower/case3_ml_s.m"); make_multiconductor!(mp_data, 3)
        result = run_mc_mld(mp_data, NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 22.2; atol = 1e-1)

        @test isapprox(result["solution"]["load"]["1"]["status"], 0.6; atol = 1e-3)
    end

    @testset "3-bus line charge nfa mld" begin
        mp_data = PM.parse_file("../test/data/matpower/case3_ml_lc.m"); make_multiconductor!(mp_data, 3)
        result = run_mc_mld(mp_data, NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 43.8; atol = 1e-1)

        @test isapprox(result["solution"]["load"]["1"]["status"], 0.544; atol = 1e-3)
    end

    @testset "transformer nfa mld" begin
        result = run_mc_mld(ut_trans_2w_yy, NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.411, atol = 1e-3)
        @test isapprox(result["solution"]["load"]["load1"]["status"], 1.0, atol = 1e-3)
    end

    @testset "5-bus lpubfdiag mld" begin
        result = run_mc_mld(case5, LPUBFDiagPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.1557; atol = 1e-4)

        @test isapprox(result["solution"]["load"]["1"]["status"], 1.000; atol = 1e-3)

        @test_throws(TESTLOG, MethodError, run_mc_mld_bf(case5, NFAPowerModel, ipopt_solver))
    end

    @testset "3-bus lpubfdiag mld" begin

        result = run_mc_mld(case3_ml, LPUBFDiagPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 24.98; atol = 1e-1)

        @test isapprox(result["solution"]["load"]["1"]["status"], 0.313; atol = 1e-3)
    end

    @testset "transformer case" begin
        result = run_mc_mld(ut_trans_2w_yy, LPUBFDiagPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-3)
        @test isapprox(result["solution"]["load"]["load1"]["status"], 1.0; atol=1e-3)
    end

    @testset "5-bus acp mld_uc" begin
        result = run_mc_mld_uc(case5, ACPPowerModel, juniper_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.338; atol=1e-3)
        @test all_gens_on(result)
        @test all_voltages_on(result)
    end

    @testset "5-bus nfa mld_uc" begin
        result = run_mc_mld_uc(case5, NFAPowerModel, juniper_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        for (i, gen) in result["solution"]["gen"]
            if i != "4"
                @test isapprox(gen["gen_status"], 1.0; atol=1e-5)
            end
        end
        @test all_voltages_on(result)
    end
end
