@info "running multinetwork tests"

@testset "test multinetwork" begin
    @testset "3-bus balanced multinetwork nfa opb" begin
        eng_ts = parse_file("../test/data/opendss/case3_balanced.dss"; time_series="daily")
        result_mn = PowerModelsDistribution._solve_mn_mc_opb(eng_ts, NFAPowerModel, ipopt_solver)

        @test result_mn["termination_status"] == LOCALLY_SOLVED
    end
end
