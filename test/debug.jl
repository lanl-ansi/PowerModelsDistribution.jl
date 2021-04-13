@info "running debug (pdf) tests"

@testset "test pbs" begin
    @testset "case 3 unbalanced - ac pf pbs" begin
        mp_data = parse_file("../test/data/opendss/case3_unbalanced.dss")
        result = solve_mc_pf_pbs(mp_data, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-4)
    end
end
