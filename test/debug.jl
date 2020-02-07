@info "running debug (pdf) tests"

@testset "test pbs" begin
    @testset "case 3 unbalanced - ac pf pbs" begin
        mp_data = PMD.parse_file("../test/data/opendss/case3_unbalanced.dss")
        result = PMD.run_mc_pf_pbs(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-4)
    end
end
