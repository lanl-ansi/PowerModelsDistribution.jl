@info "running minimum load delta (mld) tests"

@testset "test mld" begin
    @testset "transformer nfa mld" begin
        result = solve_mc_mld(ut_trans_2w_yy, NFAUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 411, atol=1)
        @test isapprox(result["solution"]["load"]["load1"]["status"], 1.0, atol = 1e-3)
    end

    @testset "transformer lpubfdiag mld" begin
        result = solve_mc_mld(ut_trans_2w_yy, LPUBFDiagPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 421; atol=1)
        @test isapprox(result["solution"]["load"]["load1"]["status"], 1.0; atol=1e-3)
    end
end
