@info "running debug (pdf) tests"

@testset "test pbs" begin
    @testset "case 3 unbalanced - ac pf pbs" begin
        result = solve_mc_pf_pbs(case3_unbalanced, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-4)
    end
end
