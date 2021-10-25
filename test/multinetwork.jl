@info "running multinetwork tests"

@testset "test multinetwork" begin
    @testset "3-bus balanced multinetwork nfa opb" begin
        eng_ts = make_multinetwork(case3_balanced)
        result_mn = PowerModelsDistribution._solve_mn_mc_opb(eng_ts, NFAUPowerModel, ipopt_solver)

        @test result_mn["termination_status"] == LOCALLY_SOLVED
    end
end
