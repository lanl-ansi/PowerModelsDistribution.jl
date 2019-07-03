@testset "test nfa opf" begin
    @testset "5-bus case" begin
        mp_data = PMs.parse_file("$(pms_path)/test/data/matpower/case5.m"); PMs.make_multiconductor!(mp_data, 3)
        result = run_tp_mld(mp_data, PMs.NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol = 1e-4)

        @test isapprox(result["solution"]["load"]["1"]["status"], 1.000; atol = 1e-3)
    end
    @testset "3-bus case" begin
        mp_data = PMs.parse_file("../test/data/matpower/case3_ml.m"); PMs.make_multiconductor!(mp_data, 3)
        result = run_tp_mld(mp_data, PMs.NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 2.00; atol = 1e-1)

        @test isapprox(result["solution"]["load"]["1"]["status"], 0.607; atol = 1e-3)
    end
    @testset "3-bus shunt case" begin
        mp_data = PMs.parse_file("../test/data/matpower/case3_ml_s.m"); PMs.make_multiconductor!(mp_data, 3)
        result = run_tp_mld(mp_data, PMs.NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 1.83; atol = 1e-1)

        @test isapprox(result["solution"]["load"]["1"]["status"], 0.589; atol = 1e-3)
    end
    @testset "3-bus line charge case" begin
        mp_data = PMs.parse_file("../test/data/matpower/case3_ml_lc.m"); PMs.make_multiconductor!(mp_data, 3)
        result = run_tp_mld(mp_data, PMs.NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 4.19; atol = 1e-1)

        @test isapprox(result["solution"]["load"]["1"]["status"], 0.486; atol = 1e-3)
    end
    @testset "transformer" begin
        mp_data = PowerModelsDistribution.parse_file("$(pmd_path)/test/data/opendss/ut_trans_2w_yy.dss")
        result = run_tp_mld(mp_data, PMs.NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0, atol = 1e-3)
        @test isapprox(result["solution"]["load"]["1"]["status"], 1.0, atol = 1e-3)
    end
end