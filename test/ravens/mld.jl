@info "ravens - running minimum load delta (mld) tests"

@testset "test mld" begin
    @testset "transformer nfa mld" begin
        math = transform_data_model(ravens_ut_trans_2w_yy)
        result = solve_mc_mld(math, NFAUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        @test isapprox(result["objective"], 0, atol=1)
        @test isapprox(result["solution"]["load"]["1"]["status"], 1.0, atol = 1e-3)
    end

    @testset "transformer lpubfdiag mld" begin
        math = transform_data_model(ravens_ut_trans_2w_yy)
        result = solve_mc_mld(math, LPUBFDiagPowerModel, ipopt_solver)

        @test result["termination_status"] == NORM_LIMIT
        #ACTUAL: @test result["termination_status"] == OTHER_ERROR

        @test isapprox(result["objective"], 36; atol=1) #evaluates 0.43==36
        #ACTUAL: @test isapprox(result["objective"], 46; atol=1)
        @test isapprox(result["solution"]["load"]["1"]["status"], .75; atol=1e-3) #0.75 == 1
        #ACTUAL: @test isapprox(result["solution"]["load"]["1"]["status"], 0.623; atol=1e-3)
    end
end
