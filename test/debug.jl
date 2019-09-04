@info "running debug.jl tests"
@testset "test ac opf pbs" begin

    @testset "5-bus coupled meshed infeasible case" begin
        @testset "ac case" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_c_m_b.m")
            result = PMD.run_mc_opf_pbs(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 26.4471; atol = 1e-2)
        end
    end

end

@testset "test soc opf pbs" begin

    @testset "5-bus coupled meshed infeasible case" begin
        mp_data = PMD.parse_file("../test/data/matlab/case5_c_m_b.m")
        result = PMD.run_mc_pf_pbs(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol = 1e-4)
    end

end


@testset "test soc pf pbs" begin

    @testset "5-bus coupled meshed network (b)" begin
        mp_data = PMD.parse_file("../test/data/matlab/case5_c_m_b.m")
        result = PMD.run_mc_pf_pbs(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-4)
    end

end

@testset "test soc pf pbs" begin

    @testset "5-bus coupled meshed network (b)" begin
        mp_data = PMD.parse_file("../test/data/matlab/case5_c_m_b.m")
        result = PMD.run_mc_pf_pbs(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-4)
    end

end