@info "running debug (pdf) tests"

@testset "test pbs" begin
    @testset "5-bus coupled meshed infeasible case - acp opf pbs" begin
        mp_data = PMD.parse_file("../test/data/matlab/case5_c_m_b.m")
        result = PMD.run_mc_opf_pbs(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 26.4471; atol = 1e-2)
    end

    @testset "5-bus coupled meshed infeasible case - soc opf pbs" begin
        mp_data = PMD.parse_file("../test/data/matlab/case5_c_m_b.m")
        result = PMD.run_mc_pf_pbs(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol = 1e-4)
    end

    @testset "5-bus coupled meshed network (b) - ac pf pbs" begin
        mp_data = PMD.parse_file("../test/data/matlab/case5_c_m_b.m")
        result = PMD.run_mc_pf_pbs(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-4)
    end

    @testset "5-bus coupled meshed network (b) - soc pf pbs" begin
        mp_data = PMD.parse_file("../test/data/matlab/case5_c_m_b.m")
        result = PMD.run_mc_pf_pbs(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-4)
    end
end