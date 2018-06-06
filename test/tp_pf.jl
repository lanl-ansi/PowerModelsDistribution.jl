@testset "test ac polar pf" begin
    @testset "5-bus independent meshed network" begin
        result = run_ac_tp_pf("../test/data/matlab/case5_i_m_b.m", ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "5-bus coupled meshed network" begin
        result = run_ac_tp_pf("../test/data/matlab/case5_c_m_a.m", ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "5-bus independent radial w/ shunts" begin
        result = run_ac_tp_pf("../test/data/matlab/case5_i_r_b.m", ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end
end