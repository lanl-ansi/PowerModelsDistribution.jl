
@testset "test ac opf pbs" begin

    @testset "5-bus coupled meshed infeasible case" begin
        @testset "ac case" begin
            mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_m_b.m")
            result = ThreePhasePowerModels.run_tp_opf_pbs(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["status"] == :LocalOptimal
            @test isapprox(result["objective"], 26.4471; atol = 1e-2)
        end
    end

end

@testset "test soc opf pbs" begin

    @testset "5-bus coupled meshed infeasible case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_m_b.m")
        result = ThreePhasePowerModels.run_tp_pf_pbs(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0; atol = 1e-4)
    end

end


@testset "test soc pf pbs" begin

    @testset "5-bus coupled meshed network (b)" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_m_b.m")
        result = ThreePhasePowerModels.run_tp_pf_pbs(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0; atol=1e-4)
    end

end

@testset "test soc pf pbs" begin

    @testset "5-bus coupled meshed network (b)" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_m_b.m")
        result = ThreePhasePowerModels.run_tp_pf_pbs(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0; atol=1e-4)
    end

end