@testset "test mld" begin
    @testset "test ACP" begin
        @testset "5-bus case" begin
            mp_data = PMs.parse_file("$(pms_path)/test/data/matpower/case5.m"); PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_mld(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.3295; atol = 1e-4)
            @test all(isapprox(result["solution"]["load"]["1"]["pd"], [3.0, 3.0, 3.0]; atol=1e-4))
            @test all(isapprox(result["solution"]["load"]["1"]["qd"], [0.9861, 0.9861, 0.9861]; atol=1e-4))

            @test isapprox(result["solution"]["load"]["1"]["status"], 1.000; atol = 1e-3)
        end
        @testset "5-bus storage case" begin
            mp_data = PMs.parse_file("$(pms_path)/test/data/matpower/case5_strg.m"); PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_mld_strg(mp_data, PMs.ACPPowerModel, ipopt_solver)  # why higher tol required?

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.3553; atol=1e-4)

            @test isapprox(result["solution"]["load"]["1"]["status"], 1.000; atol=1e-3)
        end
    end
    @testset "test NFA" begin
        @testset "5-bus case" begin
            mp_data = PMs.parse_file("$(pms_path)/test/data/matpower/case5.m"); PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_mld(mp_data, PMs.NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.1557; atol = 1e-4)

            @test isapprox(result["solution"]["load"]["1"]["status"], 1.000; atol = 1e-3)
        end
        @testset "5-bus storage case" begin
            mp_data = PMs.parse_file("$(pms_path)/test/data/matpower/case5_strg.m"); PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_mld_strg(mp_data, PMs.NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.1557; atol=1e-4)

            @test isapprox(result["solution"]["load"]["1"]["status"], 1.000; atol=1e-3)
        end
        @testset "3-bus case" begin
            mp_data = PMs.parse_file("../test/data/matpower/case3_ml.m"); PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_mld(mp_data, PMs.NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 24.9; atol = 1e-1)

            @test isapprox(result["solution"]["load"]["1"]["status"], 0.612; atol = 1e-3)
        end
        @testset "3-bus shunt case" begin
            mp_data = PMs.parse_file("../test/data/matpower/case3_ml_s.m"); PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_mld(mp_data, PMs.NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 22.2; atol = 1e-1)

            @test isapprox(result["solution"]["load"]["1"]["status"], 0.6; atol = 1e-3)
        end
        @testset "3-bus line charge case" begin
            mp_data = PMs.parse_file("../test/data/matpower/case3_ml_lc.m"); PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_mld(mp_data, PMs.NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 43.8; atol = 1e-1)

            @test isapprox(result["solution"]["load"]["1"]["status"], 0.496; atol = 1e-3)
        end
        @testset "transformer" begin
            mp_data = PowerModelsDistribution.parse_file("../test/data/opendss/ut_trans_2w_yy.dss")
            result = run_mc_mld(mp_data, PMs.NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.411, atol = 1e-3)
            @test isapprox(result["solution"]["load"]["1"]["status"], 1.0, atol = 1e-3)
        end
    end
    @testset "test LPLinUBFPowerModel" begin
        @testset "5-bus case" begin
            mp_data = PMs.parse_file("$(pms_path)/test/data/matpower/case5.m"); PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_mld_bf(mp_data, LPLinUBFPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.1557; atol = 1e-4)

            @test isapprox(result["solution"]["load"]["1"]["status"], 1.000; atol = 1e-3)

            @test_throws(TESTLOG, ErrorException, run_mc_mld_bf(mp_data, PMs.NFAPowerModel, ipopt_solver))
        end
        @testset "3-bus case" begin
            mp_data = PMs.parse_file("../test/data/matpower/case3_ml.m"); PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_mld_bf(mp_data, LPLinUBFPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 24.98; atol = 1e-1)

            @test isapprox(result["solution"]["load"]["1"]["status"], 0.313; atol = 1e-3)
        end
        @testset "transformer case" begin
            dss = PowerModelsDistribution.parse_file("../test/data/opendss/ut_trans_2w_yy.dss")
            result = run_mc_mld_bf(dss, LPLinUBFPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0; atol=1e-3)
            @test isapprox(result["solution"]["load"]["1"]["status"], 1.0; atol=1e-3)
        end
    end
    @testset "test mld_uc" begin
        @testset "test ACPPowerModel" begin
            @testset "test 5-bus case" begin
                mp_data = PMs.parse_file("$(pms_path)/test/data/matpower/case5.m"); PMD.make_multiconductor!(mp_data, 3)
                result = run_mc_mld_uc(mp_data, PMs.ACPPowerModel, juniper_solver)

                @test result["termination_status"] == PMs.LOCALLY_SOLVED
                @test isapprox(result["objective"], 0.330; atol=1e-3)
                @test all_gens_on(result)
                @test all_voltages_on(result)
            end
        end
        @testset "test NFAPowerModel" begin
            @testset "test 5-bus case" begin
                mp_data = PMs.parse_file("$(pms_path)/test/data/matpower/case5.m"); PMD.make_multiconductor!(mp_data, 3)
                result = run_mc_mld_uc(mp_data, PMs.NFAPowerModel, juniper_solver)

                @test result["termination_status"] == PMs.LOCALLY_SOLVED
                for (i, gen) in result["solution"]["gen"]
                    if i != "4"
                        @test isapprox(gen["gen_status"], 1.0; atol=1e-5)
                    end
                end
                @test all_voltages_on(result)
            end
        end
    end
end
