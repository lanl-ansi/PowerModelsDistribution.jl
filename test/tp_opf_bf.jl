@testset "test soc distflow opf_bf" begin
    @testset "3-bus case" begin
        mp_data = PowerModels.parse_file("../../PowerModels/test/data/matpower/case3.m")
        PowerModels.make_multiphase(mp_data, 3)
        result = run_tp_opf_bf(mp_data, PMs.SOCBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45829.6 ; atol = 1e0)
    end
    @testset "5-bus case" begin
        mp_data = PowerModels.parse_file("../test/data/matlab/case5.m")
        PowerModels.make_multiphase(mp_data, 3)
        result = run_tp_opf_bf(mp_data, PMs.SOCBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45040; atol = 1e0)
    end
    @testset "5-bus independent radial identical case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_a.m")
        result = run_tp_opf_bf(mp_data, PMs.SOCBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 55063.7; atol = 1e-1)
    end
    @testset "5-bus independent radial different case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_b.m")
        result = run_tp_opf_bf(mp_data, PMs.SOCBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 55572.1; atol = 1e-1)
    end
end